
# /cm/shared/apps/R/3.2.4/bin/R --save --args IBS_family.txt IBSnames.csv /uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/IBS_family < /home/jding0/siCHILD_toGC/specialIBSfamily.R


#rm(list=ls())
unlink(".RData")



args <- commandArgs(TRUE)
Gtypes_upload <- args[1]
compareCol <- args[2]
inputPath <- args[3]



print("Setting up...")
chromOpt <- "GenomeByChromosome"
genomebuild <- "vGRCh37"
siCHILD_DIR <- "/home/jding0/siCHILD_toGC/"
compiled = siCHILD_DIR
Gtypes_upload <- paste(inputPath,Gtypes_upload,sep="/")
compareCol <- paste(inputPath,compareCol,sep="/")
outPath <- paste0(inputPath,"_out/")

setwd <- outPath

print("Loading library...")
library(limma)
library(signal)
library(plotrix)
library(MASS)
library(signal)
library(plyr)
library(GenomicRanges)
library(kimisc)

print("Loading sources...")
source(paste(siCHILD_DIR,"SNPduoFunctions.R",sep=""))
source(paste(siCHILD_DIR,"snpduo_for_IBS2.R",sep=""))


print("Loading input file ...")
compareName <- read.table(compareCol,header=T,stringsAsFactors=F)
print(compareName)
dataraw <- read.table(Gtypes_upload,sep="\t",header=T)



for (n in 1:nrow(compareName)){
		col1 <- grep(paste0(compareName[n,1],".GType"),colnames(dataraw))
		col2 <- grep(paste0(compareName[n,2],".GType"),colnames(dataraw))
		snpduo_for_IBS2( Gtypes_upload, chromOpt, col1, col2, genomebuild )
		Family <- paste(compareName[n,1],compareName[n,2],sep="_")

		IBS0_file <- paste(inputPath,list.files(path = inputPath, pattern = "IBS0"),sep="/")
		IBS0 <- do.call("rbind",lapply(IBS0_file,read.csv, header=TRUE, sep="\t", stringsAsFactors=F))
		IBS0 <- IBS0[order(IBS0$Chr, IBS0$StartPosition),]
		IBS0 <- IBS0[,-4]
		colnames(IBS0) <- c("chrom","start","end")
		write.table(IBS0,paste(outPath,Family,"_IBS0.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

		IBS1_file <- paste(inputPath,list.files(path = inputPath, pattern = "IBS1"),sep="/")
		IBS1 <- do.call("rbind",lapply(IBS1_file,read.csv, header=TRUE, sep="\t", stringsAsFactors=F))
		IBS1 <- IBS1[order(IBS1$Chr, IBS1$StartPosition),]
		IBS1 <- IBS1[,-4]
		colnames(IBS1) <- c("chrom","start","end")
		write.table(IBS1,paste(outPath,Family,"_IBS1.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

		IBS2_file <- paste(inputPath,list.files(path = inputPath, pattern = "IBS2"),sep="/")
		IBS2 <- do.call("rbind",lapply(IBS2_file,read.csv, header=TRUE, sep="\t", stringsAsFactors=F))
		IBS2 <- IBS2[order(IBS2$Chr, IBS2$StartPosition),]
		IBS2 <- IBS2[,-4]
		colnames(IBS2) <- c("chrom","start","end")
		write.table(IBS2,paste(outPath,Family,"_IBS2.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

print(paste(outPath,Family,"_IBS2.txt",sep=""))
		system_cmd <- paste0("mv ", paste0(inputPath,"/IBS_family.txtchrX_1.png ", outPath, Family, "ChrX.png"))
		print(system_cmd)
		system(system_cmd)
		system_cmd <- paste0("mv ", paste0(inputPath,"/IBS_family.txtchrX_1.bedsummaryIBS0.txt ", outPath, Family, "ChrXIBS0.txt"))
		print(system_cmd)
		system(system_cmd)
		system_cmd <- paste0("mv ", paste0(inputPath,"/IBS_family.txtchrX_1.bedsummaryIBS1.txt ", outPath, Family, "ChrXIBS1.txt"))
		print(system_cmd)
		system(system_cmd)
		system_cmd <- paste0("mv ", paste0(inputPath,"/IBS_family.txtchrX_1.bedsummaryIBS2.txt ", outPath, Family, "ChrXIBS2.txt"))
		print(system_cmd)
		system(system_cmd)


#################################################
########## overlap IBS2 with raw data	########
#################################################

#		df <- dataraw[,c(2,3,3)]
#		colnames(df) <- c("chrom","start","end")
#		gr1 = as(df, "GRanges")
#		gr2 = as(IBS2, "GRanges")
#		gr3 = as(IBS1, "GRanges")
#		gr4 = as(IBS0, "GRanges")

#		OL2<-subsetByOverlaps(gr1, gr2)
#		OL2 <- as.data.frame(OL2)
#		OL1<-subsetByOverlaps(gr1, gr3)
#		OL1 <- as.data.frame(OL1)
#		OL0<-subsetByOverlaps(gr1, gr4)
#		OL0 <- as.data.frame(OL0)

#		olaps2 <- rownames(dataraw[rownames(dataraw) %in% (rownames(as.data.frame(OL2))),])
#		olaps1 <- rownames(dataraw[rownames(dataraw) %in% (rownames(as.data.frame(OL1))),])
#		olaps0 <- rownames(dataraw[rownames(dataraw) %in% (rownames(as.data.frame(OL0))),])

#		#dataraw <- olapsdataraw
#		write.table(olaps2,paste(outPath,Family,"_olaps2.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
#		write.table(olaps1,paste(outPath,Family,"_olaps1.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
#		write.table(olaps0,paste(outPath,Family,"_olaps0.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

}

