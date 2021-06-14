#RNA-seq data QC/Mapping/post-processing pipeline

#Balbc reference genome info
#https://uswest.ensembl.org/info/data/ftp/index.html

#set working directory
homeDir = "~/Lund_Curtis_TFH/pipeline/";
setwd(homeDir);

#read in bistools
source("bisTools.R");

#set output directory
baseDir = "~/Lund_Curtis_TFH/data/"

#set STAR genome specific options, choose one from available genomes: hg19, hg38, mm9, mm10, rn6, MacaM
starGenome = select_ref_STAR( "BALBc" )

#manifest file
fqDir = homeDir;
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);

#set threads amount
threads = 8
sortMem = "3G"

#STAR call
starCmd = paste("STAR --runThreadN", threads, "--outSAMtype BAM Unsorted --outSAMunmapped Within KeepPairs --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --alignIntronMax 1000000 --outSAMmultNmax 1 --genomeDir", starGenome,"--sjdbGTFfile /Volumes/GRAID/seqTools/genomes/BALBc/Mus_musculus_balbcj.BALB_cJ_v1.101.gtf --sjdbGTFchrPrefix chr --outFileNamePrefix")

#picard call
picardCmd = "java -Xmx16G -jar picard.jar MarkDuplicates "

#file extensions
bamSortExt = ".sort"
starBamFile = "Aligned.out.bam"

#adapter trimming options, uncomment one option
#seq_platform = "nextera"  #Nextera Tn5 adapters
seq_platform = "illumina" #TRUSEQ adapters

#set read counts data frame for each chromosome
chrCounts = NA

for (i in 1:dim(files)[1]) {
	if (files$include[i]) {

		print(paste(files$sample[i], Sys.time()))
		start = Sys.time();

		#files$dir[i] = mvFiles(files$fqFile[i], files$dir[i], paste0(baseDir, files$sample[i], "/"))
		mvFiles(files$fqMate1[i], files$dir[i], paste0(baseDir, files$sample[i], "/"))
		files$dir[i] = mvFiles(files$fqMate2[i], files$dir[i], paste0(baseDir, files$sample[i], "/"))

		#format and concatenate fastq files
		#files$fqFile_input[i] = files$fqFile[i]
		#files$fqFile[i] = formatFastq(files$fqFile[i], files$dir[i], files$sample[i], threads=threads)
		files$fqMate1_input[i] = files$fqMate1[i]
		files$fqMate2_input[i] = files$fqMate2[i]
		files$fqMate1[i] = formatFastq(files$fqMate1[i], files$dir[i], paste0(files$sample[i], "_1"), threads=threads)
		files$fqMate2[i] = formatFastq(files$fqMate2[i], files$dir[i], paste0(files$sample[i], "_2"), threads=threads)

		#system(paste0("fastqc ", files$dir[i],files$fqFile[i], " -o ", files$dir[i]))
		system(paste0("fastqc ", files$dir[i],files$fqMate1[i], " -o ", files$dir[i]))
		system(paste0("fastqc ", files$dir[i],files$fqMate2[i], " -o ", files$dir[i]))
		
		#cut 3 prime adapter sequences
		#files$fqFile[i] = fqCutadapt(files$fqFile[i], files$dir[i], seq_platform, threads=threads)
		pelist = fqPECutadapt(files$fqMate1[i], files$fqMate2[i], files$dir[i], seq_platform, threads=threads)
		files$fqMate1[i] = pelist[[1]]
		files$fqMate2[i] = pelist[[2]]

		#single-end or paired-end STAR call
		if (!is.null(files$fqFile[i]) && !is.na(files$fqFile[i]) && files$fqFile[i] != "") {
			starCall = paste(starCmd, files$dir[i], "--readFilesIn", paste0(files$dir[i], files$fqFile[i]))
			system(starCall)
			files$fqFile[i] = gzipFastq(paste0(files$dir[i], files$fqFile[i]), threads=threads);
		} else if (!is.null(files$fqMate1[i]) && !is.na(files$fqMate1[i]) && files$fqMate1[i] != "" && !is.null(files$fqMate2[i]) && !is.na(files$fqMate2[i]) && files$fqMate2[i] != "") {
			starCall = paste(starCmd, files$dir[i], "--readFilesIn", paste0(files$dir[i], files$fqMate1[i], " ", files$dir[i], files$fqMate2[i]))
			system(starCall)
			files$fqMate1[i] = gzipFastq(paste0(files$dir[i], files$fqMate1[i]), threads=threads);
			files$fqMate2[i] = gzipFastq(paste0(files$dir[i], files$fqMate2[i]), threads=threads);
		}

		#sort bam
		files$bamFile[i] = sortBam(paste0(files$dir[i], starBamFile), bamSortFile = paste0(files$dir[i], files$sample[i], bamSortExt), delBam = T, threads=threads, mem = sortMem)

		#mark duplicates
		files$bamFile[i] = markDups(paste0(files$dir[i], files$bamFile[i]), picardCmd, delBam = T)

		#make BigWig
		#bamToBigWig(paste0(files$dir[i], files$bamFile[i]), isDup = FALSE)

		#get counts
		cts = getBamCts(paste0(files$dir[i], files$bamFile[i]) )

		files$unmapped.reads[i] = cts[1];
		files$mapped.reads[i] = cts[2];
		files$unique.reads[i] = cts[3];
		files$paired.reads[i] = cts[4];

		#calculate time to run loop
		end = Sys.time();
		files$timeTaken[i] = difftime(end, start, units= "hours");

		#update files manifest
		write.table(files, file = paste0(fqDir, fqFile), sep = "\t", row.names = F, quote = F);
	}
}

###################
#plot mapping stats

#set up data frame
stats = t(data.frame(unmapped = files$unmapped.reads, dupe = files$mapped.reads - files$unique.reads, unique = files$unique.reads))
row.names(stats) = c("Unmapped", "Duplicate", "Unique")
#stats = stats[, order(unique, decreasing=T)]

#plot options
mai = c(1.4, 0.8, 0.3, 0.3)
mgp = c(2, 0.5, 0.5)
cols = c(rgb(99,99,99, maxColorValue=255), rgb(31,120,180, maxColorValue=255), rgb(77,175,74, maxColorValue=255))
options(scipen=9)

#make the barplot
cairo_pdf(file = "Mapping_Stats.pdf", height = 5, width = 5)
par(mai = mai, mgp = mgp)
barplot(stats/1e6, ylab = "Reads (millions)", col = cols, legend = rownames(stats), cex.lab=0.8, cex.axis = 0.8, 
	cex.names=0.8, names.arg = gsub("[[:punct:]]", " ", files$sample), las = 3)
dev.off()

###########print software versions #######
#output R session
capture.output(sessionInfo(), file = paste0("Rsession.Info.", gsub("\\D", "", Sys.time()), ".txt"))

#capture versions of other tools
tools = c("skewer",
	    "fastqc",
	    "star",
	    "samtools",
	    "java")

args = c("-version", "--version", "--version", "--version",
         "-jar picard.jar MarkDuplicates --version")

names = c("skewer", "fastqc", "STAR", "samtools", "MarkDuplicates")

version = list()

for(i in 1:length(tools)) {
	version[[names[i]]] = system2(tools[i], args[i], stdout=TRUE, stderr=TRUE)
}

capture.output(version, file = paste0("Software.versions.", gsub("\\D", "", Sys.time()), ".txt"))
