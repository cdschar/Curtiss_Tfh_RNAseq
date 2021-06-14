library("DESeq2");

experimentName = "Lund.TFH.Balbc"

homeDir = "~/Lund_Curtis_TFH/deseq2/"
setwd(homeDir);

#set significance thresholds
sigFDR = 0.05
sigFC = 1

#read in manifest file
fqDir = "~/Lund_Curtis_TFH/pipeline/";
fqFiles = "RNAseq.sample.manifest.txt";
files = read.table(paste(fqDir, fqFiles, sep = ""), sep = "\t", header = TRUE, as.is = TRUE);
#files = files[files$include,]

#read in raw and rpm normalized files
covDir = "~/Lund_Curtis_TFH/coverage/";
cts = read.table(paste0(covDir, experimentName, ".geneCts.detected.RPM.3.exon.csv"), sep = ",", header = T, check.names=F)
rpkm = read.table(paste0(covDir, experimentName, ".geneRpkm.detected.RPM.3.exon.csv"), sep = ",", header = T, check.names=F);

#establish groups for comparisons
grps = unique(files$group)
grpPairs = combn(grps, 2)

##########################
# Establish DESeq object and filter genes

coldata = files[, c("sample", "group")]
row.names(coldata) = coldata$sample

#establish DEseq DGE object
dds = DESeqDataSetFromMatrix(countData=as.matrix(cts[, files$sample]),
                             colData=coldata,
                             design = ~ group)

#run deseq fitting and get results
dds <- DESeq(dds)

#######################################################################
## GLM Function
#loop through each comparison and append stats to diff matrix

for (i in 1:dim(grpPairs)[2]) {

    #set comparisons
    comp = paste0(grpPairs[2, i], ".v.", grpPairs[1, i]);
    print(paste(comp,  Sys.time()))

    #get the results table by this pairwise comparison
    res = results(dds, contrast=c("group",grpPairs[2, i],grpPairs[1, i]), pAdjustMethod="fdr")
    res = as.data.frame(res)

    res = res[, c(1,2,5,6)]
    colnames(res) = c("baseMean", "logFC", "pvalue", "fdr")

    #mark genes that are significant by FDR criteria or by FDR and logFC
    #res$sig = (res$fdr < sigFDR & abs(res$logFC) >= sigFC)
    res$sig = res$fdr < sigFDR
    
    #add header to the lrt$table slot that describes the comparison performed
    colnames(res) = paste0(comp, ".", names(res))

    #cbind diff test onto the diff matrix
    rpkm[names(res)] = NA
    rpkm[, names(res)] = res
}

#add column that tells if a gene is significant in any comparison
if( sum(grepl(".sig", colnames(rpkm)))>1 ){
	rpkm$sigAny <- apply(rpkm[,(grep(".sig"== "TRUE", rpkm))], 1, any)
} else{
	rpkm$sigAny <- rpkm[,grep(".sig"== "TRUE", rpkm)]
}

#write the whole of table of detected genes to file
write.table(rpkm, file = paste0(homeDir, "diff.glm.", experimentName, ".txt"), sep = "\t", quote = F, row.names = F)

#filter only the significant genes
sig.rpkm = rpkm[which(rpkm$sigAny == "TRUE"), ]
write.table(sig.rpkm, file = paste0("diff.significant.glm.", experimentName, ".txt"), sep = "\t", quote = F, row.names = F)

#################################################################
#write table that tells how many genes are significant for each comparison
cols = colnames(rpkm)[grepl(".sig", colnames(rpkm))]
stats = data.frame(Comp = gsub(".sig", "", cols))

for(i in 1:dim(stats)[1]){

  print(stats$Comp[i])
  print(table(rpkm[, cols[i]]))

  stats$notSig[i] = as.numeric(table(rpkm[, cols[i]])[1])
  stats$Sig[i] = as.numeric(table(rpkm[, cols[i]])[2])

}

write.table(stats, file = paste0("Sig.stats.glm.", experimentName, ".txt"), sep = "\t", quote = F, row.names = F);

#################################################################
#write rnk files for GSEA for all comps
#specify output folder
rnkDir = "rnk_files/"
if (!file.exists(rnkDir)) dir.create(rnkDir);

#identify data columns
fc_cols = grep(".logFC", names(rpkm))
p_cols = grep(".pvalue", names(rpkm))

for(i in 1:length(fc_cols)){

	fc_col = fc_cols[i]
	p_col = p_cols[i]
  valid = !is.na(rpkm[,p_col])
	comp_name = gsub(".logFC", "", names(rpkm)[fc_col])

	rank = sign(rpkm[,fc_col]) * -log10(rpkm[,p_col]+1e-200)
	out_df = data.frame( GENE=toupper(rpkm$SYMBOL), rnk=rank )
  out_df = out_df[valid, ]
	write.table(out_df, paste0(rnkDir, comp_name, ".rnk"), row.names=F, col.names=F, sep="\t", quote=F)
}

###########print software versions #######
#output R session
capture.output(sessionInfo(), file = paste0("Rsession.Info.", gsub("\\D", "", Sys.time()), ".txt"))

