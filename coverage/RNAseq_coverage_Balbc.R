##############
#
#Requires the manifest file to have a group designations for each sample for the detection filtering

library("Rsamtools")
library("rtracklayer")
library("GenomicRanges")
library("GenomicFeatures")
library("plotly")
library("htmlwidgets")

experimentName = "Lund.TFH.Balbc"

#set working directory
homeDir = "~/Lund_Curtis_TFH/coverage/";
setwd(homeDir);

#read in manifest file
fqDir = "~/Lund_Curtis_TFH/pipeline/";
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);

#make custom TxDb object from Balbc gff file

#gff3 annotation file
#rename chromosomes to have chr prefix, only do this once
#gffinputFile = "Mus_musculus_balbcj.BALB_cJ_v1.101.gff3"
#gffFile = "Mus_musculus_balbcj.BALB_cJ_v1.101.v2.gff3"

#gff = readGFF(gffinputFile)
#gff$seqid = paste0("chr", gff$seqid)
#export.gff3(gff, gffFile, format = "GFF3")

#chromosome names/lengths
si = seqinfo(BamFile(paste0(files$dir[1], files$bamFile[1])))
chrominfo = data.frame(si@seqnames, si@seqlengths, si@is_circular)
colnames(chrominfo) = c("chrom", "length", "is_circular")
#chrominfo$chrom = gsub("chr", "", chrominfo$chrom)

metadata <- data.frame(name="Resource URL",
                       value=paste0("ftp://ftp.ensembl.org/pub/release-101/gff3/mus_musculus_balbcj"))

balbc.txdb = makeTxDbFromGFF(file=gffFile,
                         chrominfo=chrominfo,
                         dataSource="ensemblgenomes",
                         organism="Mus musculus",
                         metadata=metadata)

#set up exon and gene tables for annotation of reads
genes = genes(balbc.txdb)
exons = exonsBy(balbc.txdb, by = "gene")

geneCts = matrix(0, ncol = dim(files)[1], nrow = length(exons), dimnames = list(names(exons), files$sample))

for (i in 1:dim(files)[1]) {

	print(paste(files$sample[i], Sys.time()))

	si = seqinfo(BamFile(paste0(files$dir[i], files$bamFile[i])))

	sbp = ScanBamParam(flag = scanBamFlag(isDuplicate = F))

	#read in paired alignment
	reads = readGAlignmentPairs(paste0(files$dir[i], files$bamFile[i]), param = sbp)

	#Get counts
	so = summarizeOverlaps(exons, reads, mode = "IntersectionNotEmpty", singleEnd = F, ignore.strand = T)
	cts = assays(so)$counts

	geneCts[rownames(cts), files$sample[i]] = cts
	
	rm(reads)
	gc()
}

#Get gene symbol from entrez id row names
#geneSym = select(org, columns = "SYMBOL", keytype = "ENTREZID", keys = row.names(geneCts))
geneSym = data.frame(SYMBOL = row.names(geneCts))
geneCts = cbind(geneSym, geneCts)

#add in gene length info
width = sum(width(exons))
geneCts$length = width

#save Cts for all genes
write.table(geneCts, file = "geneCts.exon.csv", sep = ",", quote = F, row.names = F)

#make RPM normalized gene counts file
geneRpm = geneCts
geneRpm = geneRpm[, grepl(paste(files$sample, collapse = "|"), colnames(geneRpm))]
colnames(geneRpm) = paste0(files$sample, ".rpm")
geneRpm = t(t(geneRpm) * 1e6 / colSums(geneRpm))
geneRpm = round(geneRpm, 3)
geneRpm = cbind(geneCts[, !grepl(paste(files$sample, collapse = "|"), colnames(geneCts))], geneRpm)

#make RPKM normalized gene counts file
geneRpkm = geneRpm
geneRpkm = geneRpkm[, grepl(paste(files$sample, collapse = "|"), colnames(geneRpkm))]
colnames(geneRpkm) = gsub(".rpm", ".rpkm", names(geneRpkm))
geneRpkm = (geneRpkm * 1000)/ geneCts$length
geneRpkm = round(geneRpkm, 3)
geneRpkm = cbind(geneCts[, !grepl(paste(files$sample, collapse = "|"), colnames(geneCts))], geneRpkm)

#plot distribution of RPM by pdf and interactive html
plot = geneRpm[,grepl(paste(files$sample, collapse = "|"), colnames(geneRpm))]

#plotly interactive html
#axis options
#x-axis
ax <- list(
  zeroline = FALSE,
  showline = TRUE,
  showgrid = FALSE,
  linecolor = toRGB("black"),
  linewidth = 1, 
  range = c(-10,15)
)
#y-axis
ay <- list(
  zeroline = FALSE,
  showticklabels = TRUE,
  showline = TRUE,
  showgrid = FALSE,
  linecolor = toRGB("black"),
  linewidth = 1, 
  range = c(0, 0.2)
)

#establish blank plot
fig = plotly_empty() %>%
  layout(xaxis = ax, yaxis = ay)

#loop through and add each sample to the plot
for(i in 1:dim(plot)[2]){
  d = density(log2(plot[,i]+0.01))
  fig = fig %>% add_trace(x = d$x, y = d$y, type = 'scatter', mode = 'lines', name = as.character(names(plot[i])))
}

#save to html
saveWidget(fig, paste0(experimentName, ".Transcript.Density.plot.html"), selfcontained = F, libdir = "lib")

#static pdf
#cutoffs
rpm10 = log2(10+0.01)
rpm5 = log2(5+0.01)
rpm3 = log2(3+0.01)

cairo_pdf(file = paste0(experimentName, ".Transcript.Density.plot.pdf"), height = 5, width = 5);
	par(mai = c(1.5, 0.75, 0.75, 0.25), mgp = c(2, 0.5, 0), family = "Arial");
	plot(NA, ylim = c(0, 0.2), xlim = c(-10,15), xlab = "log2 RPM+0.01", ylab = "density")
	#annotate where 3 cut offs are
	lines(c(rpm10, rpm10), c(0, 0.2), lty = 2, col = "blue", lwd = 1)
	lines(c(rpm5, rpm5), c(0, 0.2), lty = 2, col = "blue", lwd = 1)
	lines(c(rpm3, rpm3), c(0, 0.2), lty = 2, col = "blue", lwd = 1)
	text(rpm10+0.8, 0.15, labels = "10", col = "blue", cex = 0.8)
	text(rpm5+0.5, 0.12, labels = "5", col = "blue", cex = 0.8)
	text(rpm3-0.5, 0.11, labels = "3", col = "blue", cex = 0.8)
	#plot data
	for(i in 1:dim(plot)[2]){
		lines(density(log2(plot[,i]+0.01)), lwd = 0.2)
		}
dev.off();

#use cutoff of 3, filter for detected gene based on groups
cutoff = 3
detected = apply(geneRpm[,grepl(paste(files$sample, collapse = "|"), colnames(geneRpm))], 1, function(x) filter_detected(x, files$group, cutoff) )
geneCts_detected = geneCts[detected,]
write.table(geneCts_detected, file = paste0(experimentName, ".geneCts.detected.RPM.", cutoff, ".exon.csv"), sep = ",", quote = F, row.names = F)

geneRpm_detected = geneRpm[detected,]
write.table(geneRpm_detected, file = paste0(experimentName, ".geneRpm.detected.RPM.", cutoff, ".exon.csv"), sep = ",", quote = F, row.names = F);

geneRpkm_detected = geneRpkm[detected,]
write.table(geneRpkm_detected, file = paste0(experimentName, ".geneRpkm.detected.RPM.", cutoff, ".exon.csv"), sep = ",", quote = F, row.names = F);

###########print software versions #######
#output R session
capture.output(sessionInfo(), file = paste0("Rsession.Info.", gsub("\\D", "", Sys.time()), ".txt"))
