#rm(list=ls())
unlink(".RData")



args <- commandArgs(TRUE)
Family <- args[1]
SibPattern = c("E*_Bl*|E*_TE*")
PGD_dir <- args[2]
outPath <- args[3]
pick_seed <- args[4]    # "NA"	# as.character(args[3])	# default NA, otherwise "E11_Bl001"

#check_IBS_only <- args[7]	# default "", otherwise "yes"


print("Loading parameters and intervals...")
parametersFile <- paste0(PGD_dir, "/", Family, "_Parameters.adj")
Params <- read.table(parametersFile, sep = "\t", header = T, stringsAsFactors = F)
Parent = Params[Params$Param == "Parent", "Value"]

Parent1 = paste(Parent, "_", Family, sep = "")
GC_File = Params[Params$Param == "GC_File", "Value"]
Condition = Params[Params$Param == "Condition", "Value"]    # only for uncle/ stepparents situations
siCHILD_DIR = Params[Params$Param == "siCHILD_DIR", "Value"]


dataFile <- paste0(PGD_dir, "/", Family, ".adj")



## READ PARAMTERS ##
for (p in c("GC_File")) {Params <- Params[Params$Param != p,]}

for (i in 1 : nrow(Params)) {

    if (Params[i, "Param"] == "siCHILD_DIR" |
        Params[i, "Param"] == "Parent" |
        Params[i, "Param"] == "Seed" |
        Params[i, "Param"] == "GC_File" |
        Params[i, "Param"] == "Condition") {eval(parse(text = paste(Params[i, "Param"], "='", Params[i, "Value"], "'", sep = "")))
    } else {eval(parse(text = paste(Params[i, "Param"], "=", Params[i, "Value"])))}
}

PGD_EXPORTED_DIR <- paste(PGD_dir, "/", Family, sep = "")
Interval_file = paste(PGD_EXPORTED_DIR, "_Intervals.txt", sep = "")

if (ExcInt == 1) {Int <- read.table(Interval_file, sep = "\t", header = T, stringsAsFactors = F)}


#outPath <- paste0(PGD_dir,Family,"_out_",Seed,"_",version,"/")

if (! grepl('/$', outPath)) {outPath <- paste(outPath, "/", sep = "")}
if (! file.exists(outPath)) {
    print(paste0("create ", outPath, " directory!"))
    dir.create(outPath, showWarnings = TRUE)
}

siCHILD_DIR = paste0(siCHILD_DIR, "/")


library(limma)
library(signal)
library(plotrix)
library(MASS)
library(signal)
library(plyr)
library(GenomicRanges)
library(kimisc)

print("Loading sources...")
load(paste(siCHILD_DIR, "Ideogram_hg19.rda", sep = ""))
load(paste(siCHILD_DIR, "REF_24h_QC_illuminaCytoSNP12.rda", sep = ""))
source(paste(siCHILD_DIR, "PlotCyto.R", sep = ""))
source(paste(siCHILD_DIR, "fastPCF.R", sep = ""))
source(paste(siCHILD_DIR, "po2.R", sep = ""))
source(paste(siCHILD_DIR, "qcgtype.R", sep = ""))
source(paste(siCHILD_DIR, "qcbyparents.R", sep = ""))
source(paste(siCHILD_DIR, "avgint.R", sep = ""))
source(paste(siCHILD_DIR, "callpcf.R", sep = ""))
source(paste(siCHILD_DIR, "callpcfBAF.R", sep = ""))
source(paste(siCHILD_DIR, "callpcfBAFtry.R", sep = ""))
source(paste(siCHILD_DIR, "chrxhtypingOpt1.R", sep = ""))
source(paste(siCHILD_DIR, "chrxhtypingOpt2.R", sep = ""))
source(paste(siCHILD_DIR, "chrxhtypingOpt3.R", sep = ""))
source(paste(siCHILD_DIR, "chrxhtypingOpt4.R", sep = ""))
source(paste(siCHILD_DIR, "chrxhtypingOpt5.R", sep = ""))
source(paste(siCHILD_DIR, "chrxhtypingOpt6.R", sep = ""))
source(paste(siCHILD_DIR, "testmedfilt.R", sep = ""))
source(paste(siCHILD_DIR, "inthapnew1.R", sep = ""))
source(paste(siCHILD_DIR, "patscore2.R", sep = ""))
source(paste(siCHILD_DIR, "artscsnpspar.R", sep = ""))
source(paste(siCHILD_DIR, "htypingOpt2Aut.R", sep = ""))
source(paste(siCHILD_DIR, "meanwindow.R", sep = ""))
source(paste(siCHILD_DIR, "htypingOpt1.R", sep = ""))
source(paste(siCHILD_DIR, "htypingOpt2.R", sep = ""))
source(paste(siCHILD_DIR, "htypingOpt3.R", sep = ""))
source(paste(siCHILD_DIR, "htypingOpt4.R", sep = ""))
source(paste(siCHILD_DIR, "htypingOpt5.R", sep = ""))
source(paste(siCHILD_DIR, "htypingOpt6.R", sep = ""))
source(paste(siCHILD_DIR, "intphappropser.R", sep = ""))
source(paste(siCHILD_DIR, "intphappropserOpt3.R", sep = ""))
source(paste(siCHILD_DIR, "phasebafOpt1.R", sep = ""))
source(paste(siCHILD_DIR, "phasebafOpt2.R", sep = ""))
source(paste(siCHILD_DIR, "phasebafOpt3.R", sep = ""))
source(paste(siCHILD_DIR, "pickBestSingleCellbyChr.R", sep = ""))
source(paste(siCHILD_DIR, "phasebafOpt4.R", sep = ""))
source(paste(siCHILD_DIR, "phasebafOpt5.R", sep = ""))
source(paste(siCHILD_DIR, "phasebafOpt6.R", sep = ""))
source(paste(siCHILD_DIR, "distHap.R", sep = ""))
source(paste(siCHILD_DIR, "SNPduoFunctions.R", sep = ""))
source(paste(siCHILD_DIR, "snpduo_for_IBS2.R", sep = ""))
source(paste(siCHILD_DIR, "pickseed.R", sep = ""))
source(paste(siCHILD_DIR, "file_rename.R", sep = ""))
source(paste(siCHILD_DIR, "callrate.R", sep = ""))
source(paste(siCHILD_DIR, "mendinc.R", sep = ""))
source(paste(siCHILD_DIR, "wavecorrtmean2.R", sep = ""))
source(paste(siCHILD_DIR, "testmedfilt.R", sep = ""))
source(paste(siCHILD_DIR, "plotCytoGenome.R", sep = ""))
source(paste(siCHILD_DIR, "chrspecbafplot.R", sep = ""))
source(paste(siCHILD_DIR, "genomebafplotsc.R", sep = ""))


P1 <- NULL
P1Seg <- NULL
P2Seg <- NULL
M1 <- NULL
M1Seg <- NULL
M2 <- NULL
M2Seg <- NULL
logRsSeg <- NULL
logRs <- NULL

`%ni%` = Negate(`%in%`)

#Default parameter settings
options(scipen = 999)
Func = "mean"


is.integer0 <- function(x) {is.integer(x) && length(x) == 0L}

Chroms <- c(1 : 22, "X", "Y")

centro <- ideogram[ideogram[, 5] %in% c("acen", "gvar"),]
centro[, 1] = as.character(gsub("chr", "", centro[, 1]))
write.table(centro, paste(outPath, Family, "_centro.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")


system_cmd <- paste0("cp ", paste0(siCHILD_DIR, "human_genetic.map"), " " , paste0(outPath, Family, ".map"))
system(system_cmd)
genetic_map <- read.table(paste0(siCHILD_DIR, "human_genetic.map"), sep = " ", stringsAsFactors = F)
colnames(genetic_map) <- c("Chr", "Name", "Dist", "Position")
centro_snp <- NULL
for (chr in Chroms) {
    cn <- centro[centro$Chr %in% chr,]
    cname <- genetic_map[genetic_map$Chr %in% chr &
        genetic_map$Position >= cn[1, 2] &
        genetic_map$Position <= cn[2, 3],]$Name
    centro_snp <- c(centro_snp, cname)
}




# E02_Bl001.GType : Embryo2 blastomere1 (var nr of cols, x per blastomere)
print(paste("Reading data file:", dataFile))
data <- read.table(dataFile, sep = "\t", header = T, stringsAsFactors = F)
if (any(grepl("Affected", colnames(data)))) {
    colnames(data) <- gsub("Affected", "Sibling", colnames(data))
}
#print(str(data))
if (Parent1 == paste("Mother_", Family, sep = "")) {Parent2 = paste("Father_", Family, sep = "")}else {Parent2 = paste("Mother_", Family, sep = "")}


data <- data[, colSums(is.na(data)) != nrow(data)]    # remove columns with NA values

#for(c in c(1,2,grep("GType",colnames(data)))){data[,c]<-as.character(data[,c])}
for (c in c(grep("Log.R.Ratio", colnames(data)))) {data[, c] <- as.numeric(as.character(data[, c]))}
for (c in c(grep("B.Allele.Freq", colnames(data)))) {data[, c] <- as.numeric(as.character(data[, c]))}
dataraw <- na.omit(data)

dataraw <- dataraw[grep("cnvi", as.character(dataraw$Name), invert = TRUE),]
dataraw <- dataraw[dataraw[, "Position"] != 0,]
#dataraw <- dataraw[dataraw[,"Chr"]=="Y",]
#dataraw <- dataraw[dataraw[,"Chr"]!="XY",]



############################################################
##### QA calculation for gene and its flanking regions #####
############################################################
Parent_1 <- Parent
if (Parent_1 == "Mother") {Parent_2 <- "Father"} else {Parent_2 <- "Mother"}
genQA <- as.data.frame(matrix(nrow = 5, ncol = 2))
colnames(genQA) <- c("SNPs", "iSNPs")
rownames(genQA) <- c("up2MB","up1MB", "gene", "down1MB", "down2MB")
IntRow <- c()
for (i in 1 : nrow(Int)) {
    if (Int[i, 1] != 0) {
        IntRow <- i
    }
}

genQA_cutoff <- as.data.frame(matrix(c(-2000000,-1000000,0, 0, 0, 0, 0,0, 1000000, 2000000), nrow = 5, ncol = 2))
for (i in 1 : nrow(genQA)) {    # take upstream & downstream 2M of flanking disease gene
    genQA[i, 1] <- nrow(dataraw[(dataraw$Chr == Int[IntRow, 1]) &
        dataraw$Position > (Int[IntRow, 2] + genQA_cutoff[i, 1]) &
        dataraw$Position < (Int[IntRow, 3] + genQA_cutoff[i, 2]),])
    genQA[i, 2] <- nrow(dataraw[(dataraw$Chr == Int[IntRow, 1]) &
        dataraw$Position > (Int[IntRow, 2] + genQA_cutoff[i, 1]) &
        dataraw$Position < (Int[IntRow, 3] + genQA_cutoff[i, 2]) &
        (dataraw[grep(Parent_1, colnames(dataraw))[1]] == "AB") &
        ((dataraw[grep(Parent_2, colnames(dataraw))[1]] == "AA") | (dataraw[grep(Parent_2, colnames(dataraw))[1]] == "BB")),])
}

write.table(genQA, paste(outPath, Family, "_genQA.txt", sep = ""), col.names = T, row.names = T, quote = F, sep = "\t")





GC <- read.table(GC_File, header = F, sep = "\t")
rownames(dataraw) <- as.character(dataraw[, 1])
rownames(GC) <- as.character(GC[, 4])
rowsTot <- intersect(rownames(dataraw), rownames(GC))
GC <- GC[rowsTot,]

dataraw <- dataraw[rowsTot,]
# dataraw <- rbind(dataraw,data[data$Chr == "Y",][grep("rs",(data[data$Chr == "Y",1])),])


write.table(dataraw, paste(outPath, Family, "_dataraw.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")


if (! is.na((grep("Uncle", Seed) == 1) || (grep("Aunt", Seed) == 1))) { # || (grep("Step", Seed) == 1) )) {

    if (! (file.exists(paste(outPath, Family, "_olaps1.txt", sep = "")))) {

        ########################################
        ######### calculate IBS2	########
        ########################################
        library(GenomicRanges)
        library(plyr)
        chromOpt <- "GenomeByChromosome"

        #Gtypes_upload <- paste(outPath,Family,"_0.75.gtp",sep="")
        #col1 <-  which(colnames(Gtypes) %in% paste(Parent,"_",Family,sep=""))
        #col2 <-  which(colnames(Gtypes) %in% paste(Seed,"_",Family,sep=""))

        if (length(grep("Step", Seed) == 1) == 1) {
            Seed1 <- (unlist(strsplit(Seed, "_"))[2])
            Seed2 <- (unlist(strsplit(Seed, "_"))[3])
            c1 <- Seed1
            c2 <- Seed2
        } else if (! is.na((grep("Uncle", Seed) == 1) || (grep("Aunt", Seed) == 1))) {
            c1 <- Parent
            c2 <- Seed
        } else {
            print("check!")
        }


        col1 <- which(colnames(dataraw) %in% paste(c1, ".GType", sep = ""))
        col2 <- which(colnames(dataraw) %in% paste(c2, ".GType", sep = ""))


        #dataclean <- dataraw[dataraw[,col1] != "NC" | dataraw[,col2] != "NC",]
        #write.table(dataclean,paste(outPath,Family,"_dataclean.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
        #Gtypes_upload <- paste(outPath,Family,"_dataclean.txt",sep="")
        ##Gtypes_upload <- "/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD065/PGD065_short.adj"


        Gtypes_upload <- paste(outPath, Family, "_dataraw.txt", sep = "")


        #print(paste(Gtypes_upload, chrom, col1, col2, compile))
        #		compiled =  "/home/jding0/programs/snpduoweb-master/cgi-bin/"
        compiled = siCHILD_DIR
        genomebuild = "vGRCh37"
        snpduo_for_IBS2(Gtypes_upload, chromOpt, col1, col2, genomebuild)



        if (file.exists(paste0(outPath, Family, "_IBS0.txt"))) {
            print("----->IBS0 file exists, must first delete it!!!")
            system_cmd <- paste0("rm ", paste0(outPath, Family, "_IBS0.txt"))
            print(system_cmd)
            system(system_cmd)
        }

        IBS0_file <- paste(outPath, list.files(path = outPath, pattern = "IBS0"), sep = "/")
        IBS0 <- do.call("rbind", lapply(IBS0_file, read.csv, header = TRUE, sep = "\t", stringsAsFactors = F))
        IBS0 <- IBS0[order(IBS0$Chr, IBS0$StartPosition),]
        IBS0 <- IBS0[, - 4]
        colnames(IBS0) <- c("chrom", "start", "end")
        write.table(IBS0, paste(outPath, Family, "_IBS0.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")


        if (file.exists(paste0(outPath, Family, "_IBS1.txt"))) {
            print("----->IBS1 file exists, must first delete it!!!")
            system_cmd <- paste0("rm ", paste0(outPath, Family, "_IBS1.txt"))
            print(system_cmd)
            system(system_cmd)
        }

        IBS1_file <- paste(outPath, list.files(path = outPath, pattern = "IBS1"), sep = "/")
        IBS1 <- do.call("rbind", lapply(IBS1_file, read.csv, header = TRUE, sep = "\t", stringsAsFactors = F))
        IBS1 <- IBS1[order(IBS1$Chr, IBS1$StartPosition),]
        IBS1 <- IBS1[, - 4]
        colnames(IBS1) <- c("chrom", "start", "end")
        write.table(IBS1, paste(outPath, Family, "_IBS1.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")


        if (file.exists(paste0(outPath, Family, "_IBS2.txt"))) {
            print("----->IBS2 file exists, must first delete it!!!")
            system_cmd <- paste0("rm ", paste0(outPath, Family, "_IBS2.txt"))
            print(system_cmd)
            system(system_cmd)
        }

        IBS2_file <- paste(outPath, list.files(path = outPath, pattern = "IBS2"), sep = "/")
        IBS2 <- do.call("rbind", lapply(IBS2_file, read.csv, header = TRUE, sep = "\t", stringsAsFactors = F))
        IBS2 <- IBS2[order(IBS2$Chr, IBS2$StartPosition),]
        IBS2 <- IBS2[, - 4]
        colnames(IBS2) <- c("chrom", "start", "end")
        write.table(IBS2, paste(outPath, Family, "_IBS2.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")


        ################################################
        ######### overlap IBS2 with raw data	########
        ################################################

        df <- dataraw[, c(2, 3, 3)]
        colnames(df) <- c("chrom", "start", "end")
        gr1 = as(df, "GRanges")
        gr2 = as(IBS2, "GRanges")
        gr3 = as(IBS1, "GRanges")
        gr4 = as(IBS0, "GRanges")

        OL2 <- subsetByOverlaps(gr1, gr2)
        OL2 <- as.data.frame(OL2)
        OL1 <- subsetByOverlaps(gr1, gr3)
        OL1 <- as.data.frame(OL1)
        OL0 <- subsetByOverlaps(gr1, gr4)
        OL0 <- as.data.frame(OL0)

        olaps2 <- rownames(dataraw[rownames(dataraw) %in% (rownames(as.data.frame(OL2))),])
        olaps1 <- rownames(dataraw[rownames(dataraw) %in% (rownames(as.data.frame(OL1))),])
        olaps0 <- rownames(dataraw[rownames(dataraw) %in% (rownames(as.data.frame(OL0))),])

        #dataraw <- olapsdataraw
        write.table(olaps2, paste(outPath, Family, "_olaps2.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
        write.table(olaps1, paste(outPath, Family, "_olaps1.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
        write.table(olaps0, paste(outPath, Family, "_olaps0.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
    } else {
        olaps2 <- read.table(paste(outPath, Family, "_olaps2.txt", sep = ""), header = T, stringsAsFactors = F, sep = "\t")[, 1]
        olaps1 <- read.table(paste(outPath, Family, "_olaps1.txt", sep = ""), header = T, stringsAsFactors = F, sep = "\t")[, 1]
        olaps0 <- read.table(paste(outPath, Family, "_olaps0.txt", sep = ""), header = T, stringsAsFactors = F, sep = "\t")[, 1]
        IBS1 <- read.table(paste(outPath, Family, "_IBS1.txt", sep = ""), header = T, stringsAsFactors = F, sep = "\t")
        IBS0 <- read.table(paste(outPath, Family, "_IBS0.txt", sep = ""), header = T, stringsAsFactors = F, sep = "\t")
        IBS2 <- read.table(paste(outPath, Family, "_IBS2.txt", sep = ""), header = T, stringsAsFactors = F, sep = "\t")
    }
}






####################################################
### not only checking IBS, but the whole process ###
####################################################
#if (check_IBS_only == "no") {	### not only checking IBS, but the whole process ###



logRsRaw <- cbind(dataraw[, c("Name", "Chr", "Position")], dataraw[, grep(".Log.R.Ratio", colnames(dataraw))])
colnames(logRsRaw)[- c(1 : 3)] <- paste(gsub(".Log.R.Ratio", "", colnames(logRsRaw)[- c(1 : 3)]), "_", Family, sep = "")

BAFs <- cbind(dataraw[, c("Name", "Chr", "Position")], dataraw[, grep(".B.Allele.Freq", colnames(dataraw))])
colnames(BAFs)[- c(1 : 3)] <- paste(gsub(".B.Allele.Freq", "", colnames(BAFs)[- c(1 : 3)]), "_", Family, sep = "")

Gtypes <- cbind(dataraw[, c("Name", "Chr", "Position")], dataraw[, grep(".GType", colnames(dataraw))])
colnames(Gtypes)[- c(1 : 3)] <- paste(gsub(".GType", "", colnames(Gtypes)[- c(1 : 3)]), "_", Family, sep = "")



ChrPos <- Gtypes[, c("Chr", "Position")]
QC <- qcgtype(Gtypes, ChrPos, Family, SibPattern, outPath)
save(QC, file = paste(outPath, "QC_", Family, ".rda", sep = ""))



QCbyParents <- qcbyparents(Gtypes, SibPattern)
SnpSpecArts <- artscsnpspar(Gtypes, SibPattern, outPath)
SnpSpecArts[, "Position"] <- as.numeric(SnpSpecArts[, "Position"])
dataPo <- po2(Gtypes, Family, SibPattern, outPath)

write.table(dataPo, paste(outPath, Family, "_debugParScoreGtypes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
ParScore <- patscore2(dataPo, QC, Chroms, Gtypes, Family)


write.table(logRsRaw, paste(outPath, Family, "_logRsRaw.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(BAFs, paste(outPath, Family, "_BAFs.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(Gtypes, paste(outPath, Family, "_0.75.gtp", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(dataPo, paste(outPath, Family, "_dataPo.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")



save(ParScore, file = paste(outPath, "ParScore_", Family, ".rda", sep = ""))


ParScore_exc <- ParScore[c(names(ParScore)[grep("^E", names(ParScore))])]    # only embryos

exc <- NULL
for (n in 1 : length(ParScore_exc)) {
    if (length(unique(as.numeric(ParScore_exc[[n]][1 : 22, 6]))) > 1) {
        exc <- names(ParScore_exc)[n]
    }
}





############################################################
##### sort single cell by the mean & sd of callrate #####
############################################################
order_embryo <- cbind(t(as.data.frame(lapply(ParScore, function(x) mean(x[1 : 22, 5])))), t(as.data.frame(lapply(ParScore, function(x) sd(x[1 : 22, 5])))))
colnames(order_embryo) <- c("mean", "sd")
write.table(order_embryo, paste(outPath, Family, "_orderEmbryoByCallrate.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")


library(limma)
window = Win / 2

if (file.exists(paste(outPath, Family, "_logRsAvgWindow.txt", sep = ""))) {
    print("logRsAvgWindow file exists!")
    logRs <- read.table(paste(outPath, Family, "_logRsAvgWindow.txt", sep = ""), sep = "\t", header = T, stringsAsFactors = F)
} else {
    logRs <- meanwindow(logRsRaw[, c(1 : 3, grep("^Sibling", colnames(logRsRaw)), grep("^E", colnames(logRsRaw)), grep("Father", colnames(logRsRaw)), grep("Mother", colnames(logRsRaw)))], GC, window, Func, Family, ParScore, outPath)
    #	if (!(Seed %in% c("Grandparents","Sibling"))) {
    #		logRs <- logRs[logRs$Name %in% olaps1,]
    #	}
    write.table(logRs, paste(outPath, Family, "_logRsAvgWindow.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}
data <- dataraw
gamma = gammaBAF





##########################################
##### start haplotyping and call BAF #####
##########################################

Haps <- NULL
NewSeed <- c()
dataPath <- c()
datasGtypes <- NULL

print(Parent1)
if (Seed == "Grandparents") {
   Haps <- htypingOpt1(Gtypes, dataPo, ParScore, Parent1, SibPattern, Family, outPath)
   print("-------------------------------------")
   print("Option 1 haplotyping was applied...")
   print(Parent1)
   print("-------------------------------------")

   write.table(Gtypes, paste(outPath, Family, "debugGtypes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
   write.table(dataraw, paste(outPath, Family, "debugdataraw.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
   Parents <- Haps[["Parents"]]

   PhBAF <- phasebafOpt1(Gtypes, dataraw, Family, Parent1, gamma, outPath, Parents)
} else if (Seed == "Sibling") {

   if (any(grepl("Affected", colnames(Gtypes)))) {
       colnames(Gtypes) <- gsub("Affected", "Sibling", colnames(Gtypes))
       colnames(dataraw) <- gsub("Affected", "Sibling", colnames(dataraw))
   }

   if (! any(grepl("Sibling", colnames(Gtypes)))) {
       pick_seed <- "new"
       datasGtypes <- pickseed(pick_seed, ParScore, Gtypes, dataraw, chr)
       Gtypes <- datasGtypes[["Gtypes"]]
       dataraw <- datasGtypes[["dataraw"]]
   }

   dataPo <- po2(Gtypes, Family, SibPattern, dataPath)    #### this shall not be used for the Chr plot !!!
   QC <- qcgtype(Gtypes, ChrPos, Family, SibPattern, outPath)
   ParScore <- patscore2(dataPo, QC, Chroms, Gtypes, Family)


   Haps <- htypingOpt2(Gtypes, Window, Int, dataPo, ParScore, SibPattern, Family, outPath)
   print("-------------------------------------")
   print("Option 2 haplotyping was applied...")
   print(Parent1)
   print("-------------------------------------")
   PhBAF <- phasebafOpt2(Gtypes, dataraw, Family, gamma, ParScore, outPath)
} else if (! is.na((grep("Uncle", Seed) == 1) || (grep("Aunt", Seed) == 1))) {

   datarawOrg <- dataraw
   GtypesOrg <- Gtypes
   ParScoreOrg <- ParScore

   print("-------------------------------------")
   print("Option 3.2 haplotyping (uncle/aunt) was applied...")
   print("-------------------------------------")
   print(Parent1)
   Seed <- paste(Seed, Family, sep = "_")
   Haps <- htypingOpt3(Gtypes, dataPo, ParScore, Parent1, SibPattern, Family, outPath, olaps0, olaps1, olaps2, IBS1, IBS0)    #? ParScore => ParScoreOrg?
   if (Parent1 == paste("Father_", Family, sep = "")) { par1 = "_Pat"; par2 = "_Mat"} else { par1 = "_Mat"; par2 = "_Pat"}

   dataHap <- Haps[["dataHap"]]
   dataHapRaw <- Haps[["dataHapRaw"]]
   Htype <- Haps[["Htype"]]


   library(plyr)
   pp <- rbind(dataHapRaw, dataHap)
   HapCompare <- ddply(pp, .(Name), function(x) colSums(x[, - c(1 : 3)]))
   m <- merge(HapCompare, dataHapRaw, by = "Name")
   HapCompare <- m[, c(1, (ncol(HapCompare) + 1) : (ncol(HapCompare) + 2), 2 : ncol(HapCompare))]
   colnames(HapCompare) <- gsub(".x", "", colnames(HapCompare))

   HapCompare <- HapCompare[order(HapCompare$Chr, HapCompare$Position),]
   rownames(HapCompare) <- NULL
   HapCompare <- HapCompare[, c(1 : 3, grep(par1, colnames(HapCompare)))]


   NewSeed <- "IBSseed"
   dataPath <- paste0(outPath, "/", NewSeed, "/")
   if (! file.exists(dataPath)) {
       dir.create(dataPath)
   }

   write.table(dataHap, paste(dataPath, Family, "_Itp.hap", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
   write.table(dataHapRaw, paste(dataPath, Family, "_Raw.hap", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
   write.table(GtypesOrg, paste(dataPath, Family, "_0.75.gtp", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
   write.table(Htype, paste(dataPath, Family, "_0.75.htp", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
   write.table(HapCompare, paste(dataPath, Family, "_HapCompare_before.hap", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")

   Intp <- intphappropserOpt3(dataPath, Family, ParScore, SibPattern)

   PhBAF <- phasebafOpt3(GtypesOrg, Htype, dataraw, Family, Parent1, gamma, outPath, olaps0, olaps1, olaps2, IBS1, IBS0, NewSeed, HapCompare)



   ####################################
   ##### select best single cell ######
   ####################################


   Gtypes <- GtypesOrg
   dataraw <- datarawOrg


   # <-------------------------

   if (is.na(pick_seed)) {    # select the best single cell per chr, build artificial single cell

       Gtypes[, ncol(Gtypes) + 1] <- NA
       colnames(Gtypes)[ncol(Gtypes)] <- "Sibling"
       dataraw[, ncol(dataraw) + 1] <- NA
       colnames(dataraw)[ncol(dataraw)] <- "Sibling.GType"
       dataraw[, ncol(dataraw) + 1] <- NA
       colnames(dataraw)[ncol(dataraw)] <- "Sibling.B.Allele.Freq"
       dataraw[, ncol(dataraw) + 1] <- NA
       colnames(dataraw)[ncol(dataraw)] <- "Sibling.Log.R.Ratio"

       empty_chr <- NULL
       for (chr in c(1 : 22, "X", "Y")) {
           if (! is.na(PhBAF[["pick"]][chr, 2])) {
               Gtypes[Gtypes[, 2] == chr, ncol(Gtypes)] <- Gtypes[Gtypes[, 2] == chr, PhBAF[["pick"]][chr, 2]]
               NewSiblingCol <- paste(strsplit(PhBAF[["pick"]][chr, 2], "_")[[1]][1], strsplit(PhBAF[["pick"]][chr, 2], "_")[[1]][2], sep = "_")
               dataraw[dataraw[, 2] == chr, (ncol(dataraw) - 2) : ncol(dataraw)] <- dataraw[dataraw[, 2] == chr, grep(NewSiblingCol, colnames(dataraw))]
           } else {
               empty_chr <- chr
           }
       }

       Chroms <- c(1 : 22, "X", "Y")


       #		#############################################################################################
       #		##### ChrX has problem, choose the best single cell by call rate, then choose its ChrX ######
       #		#############################################################################################
       #		#unlist(lapply(ParScore, function(x) max(x["X",5])))
       #		NewSeed_select <- cbind(t(as.data.frame(lapply(ParScore, function(x) mean(x[1:22,5])))), t(as.data.frame(lapply(ParScore, function(x) sd(x[1:22,5])))))
       #		colnames(NewSeed_select) <- c("mean","sd")
       #		NewSeedX <- rownames(NewSeed_select[order(-NewSeed_select[,1],-NewSeed_select[,2]),])[1]
       #
       #		NewAffectedCol <- paste(strsplit(NewSeedX,"_")[[1]][1], strsplit(NewSeedX,"_")[[1]][2], sep="_")

       #		Gtypes[Gtypes[,2] %in% c("X","Y",empty_chr),ncol(Gtypes)] <- Gtypes[Gtypes[,2] %in% c("X","Y",empty_chr),NewSeedX]
       #		dataraw[dataraw[,2] %in% c("X","Y",empty_chr),(ncol(dataraw)-2):ncol(dataraw)] <- dataraw[dataraw[,2] %in% c("X","Y",empty_chr),grep(NewAffectedCol,colnames(dataraw))]
   } else {    # know which single sell to use in advance,e.g. E01
       Gtypes[, ncol(Gtypes) + 1] <- Gtypes[, grep(pick_seed, colnames(Gtypes))]
       colnames(Gtypes)[ncol(Gtypes)] <- "Sibling"
       nc <- ncol(dataraw)
       dataraw[, (nc + 1) : (nc + 3)] <- dataraw[, grep(pick_seed, colnames(dataraw))]
       colnames(dataraw)[(nc + 1) : (nc + 3)] <- c("Sibling.GType", "Sibling.B.Allele.Freq", "Sibling.Log.R.Ratio")
   }


   dataPo <- po2(Gtypes, Family, SibPattern, dataPath)    #### this shall not be used for the Chr plot !!!
   QC <- qcgtype(Gtypes, ChrPos, Family, SibPattern, outPath)
   ParScore <- patscore2(dataPo, QC, Chroms, Gtypes, Family)

   Haps <- htypingOpt2(Gtypes, Window, Int, dataPo, ParScore, SibPattern, Family, outPath)
   PhBAF <- phasebafOpt2(Gtypes, dataraw, Family, gamma, ParScore, outPath)

   write.table(Gtypes, paste(dataPath, Family, "_0.75Sibling.gtp", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
   write.table(dataraw, paste(dataPath, Family, "_datarawSibling.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
} else if (length(grep("Step", Seed) == 1) == 1) {    #  no longer in use! "Stepparents_Sibling"

   Seed <- paste((unlist(strsplit(Seed, "_"))[2]), Family, sep = "_")

   Haps <- htypingOpt4(Gtypes, dataPo, ParScore, Parent1, SibPattern, Family, outPath, Condition)
   print("-------------------------------------")
   print("Option 4.2 haplotyping was applied...")
   print(Parent1)
   print("-------------------------------------")
   Parents <- Haps[["Parents"]]

   PhBAF <- phasebafOpt4(Gtypes, dataraw, Family, Parent1, gamma, outPath, Parents)
} else if (! is.na((grep("Grandmother", Seed) == 1) || (grep("Grandfather", Seed) == 1))) {
   Seed <- paste(Seed, Family, sep = "_")
   #	Seed <- paste((unlist(strsplit(Seed,"_"))[2]),Family,sep="_")
   print(paste0("====== ", Seed, "======"))
   Haps <- htypingOpt5(Gtypes, dataPo, ParScore, Parent1, SibPattern, Family, outPath)
   print("-------------------------------------")
   print("Option 5.2 haplotyping was applied...")
   print(Parent1)
   print("-------------------------------------")
   Parents <- Haps[["Parents"]]

   PhBAF <- phasebafOpt5(Gtypes, dataraw, Family, Parent1, gamma, outPath, Parents)
} else if (length(grep("Halfsib", Seed) == 1) == 1) {    # "divorced family"

   Seed <- paste("Sibling", Family, sep = "_")
   print(Seed)
   #	dataPo <- dataPo[ , !(names(dataPo) %in% Seed)]
   #	ParScore <- ParScore[ !(names(ParScore) %in% Seed)]
   Haps <- htypingOpt6(Gtypes, dataPo, ParScore, Parent1, SibPattern, Family, outPath)
   print("-------------------------------------")
   print("Option 6.2 haplotyping was applied...")
   print(Parent1)
   print("-------------------------------------")
   Parents <- Haps[["Parents"]]

   PhBAF <- phasebafOpt6(Gtypes, dataraw, Family, Parent1, gamma, outPath, Parents)
   #	NewSeed <- "HalfSib"
}




dataHap <- Haps[["dataHap"]]
dataHapRaw <- Haps[["dataHapRaw"]]

write.table(dataHap, paste(outPath, Family, "_Itp.hap", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(dataHapRaw, paste(outPath, Family, "_Raw.hap", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")


library(MASS)
if (! is.na((grep("Uncle", Seed) == 1) || (grep("Aunt", Seed) == 1))) {
   ParScore <- ParScoreOrg
   #	logRs <- logRs[logRs$Name %in% olaps1,]
   Intp <- intphappropserOpt3(outPath, Family, ParScore, SibPattern)
} else { # if ((Seed %in% c("Grandparents","Sibling")) || (grep("Step", Seed) == 1) || (grep("SingleGP", Seed) == 1) ) {
   Intp <- intphappropser(outPath, Family, ParScore, SibPattern)
}

#Intp <- intphappropser(outPath,Family,ParScore,SibPattern)
print("Finish Intp ........ ")




AvgLogRs <- avgint(logRs, Int, Family, outPath)
SegLogRs <- callpcf(logRs, gammaSC, gammaMC, plateau, Family, outPath)    # logR segment is from bestsinglecell.
print("Finish segment LogRs .... ")

if (is.null(dataPath)) {
   dataPath <- outPath
}

disthap(Family, Int, dataPath)

print("-----------------------------")
print("....... Finished !!!! .......")
print("-----------------------------")

# PpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPp
print("Plotting ......")
print(Int)
print("---------------")

P1 <- PhBAF[["P1"]]
P2 <- PhBAF[["P2"]]
M1 <- PhBAF[["M1"]]
M2 <- PhBAF[["M2"]]
P1Seg <- PhBAF[["P1Seg"]]
P2Seg <- PhBAF[["P2Seg"]]
M1Seg <- PhBAF[["M1Seg"]]
M2Seg <- PhBAF[["M2Seg"]]
logRsSeg <- read.table(paste(outPath, list.files(path = outPath, pattern = "gamma", all.files = FALSE), sep = ""), header = T, sep = "\t", stringsAsFactors = F)


##############################################################################################
#### for uncle/aunt case, use file in IBSseed; otherwise, dataPath is the same as outPath ####
##############################################################################################
dataHap <- read.table(paste(dataPath, Family, "_Itp2.hap", sep = ""), header = T, sep = "\t", stringsAsFactors = F)
dataHapRaw <- read.table(paste(dataPath, Family, "_Raw.hap", sep = ""), header = T, sep = "\t", stringsAsFactors = F)
dataPo <- read.table(paste(dataPath, Family, ".poo", sep = ""), header = T, sep = "\t", stringsAsFactors = F)
#HapCompare <- read.table(paste(dataPath,Family,"_HapCompare_after.hap",sep=""),header=T,sep="\t",stringsAsFactors=F)

# for debugging #
# dataraw <- read.table(paste(outPath,Family,"_dataraw.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)



if (! is.na((grep("Uncle", Seed) == 1) || (grep("Aunt", Seed) == 1))) {} else {
   IBS1 <- 99999
   IBS0 <- 99999
   IBS2 <- 99999
   olaps1 <- 99999
   olaps0 <- 99999
}

print("=========== Plotting ==============")


print("-------- save objects -------------")
save(dataHap, dataHapRaw, dataPo, BAFs, logRs, logRsSeg, P1, P1Seg, P2, P2Seg, M1, M1Seg, M2, M2Seg, ChrsLength, ideogram, Family, outPath, Chroms, Int, Seed, IBS2, IBS1, IBS0, olaps1, olaps0, centro, NewSeed, siCHILD_DIR, file = paste0(dataPath, "/chrspecbafplot.RData"))
chrspecbafplot(dataHap, dataHapRaw, dataPo, BAFs, logRs, logRsSeg, P1, P1Seg, P2, P2Seg, M1, M1Seg, M2, M2Seg, ChrsLength, ideogram, Family, outPath, Chroms, Int, Seed, IBS2, IBS1, IBS0, olaps1, olaps0, centro, NewSeed)

save(BAFs, logRs, logRsSeg, P1, P1Seg, P2, P2Seg, M1, M1Seg, M2, M2Seg, ChrsLength, ideogram, Family, outPath, Seed, dataHap, siCHILD_DIR, file = paste0(dataPath, "/genomebafplot.RData"))
genomebafplot(BAFs, logRs, logRsSeg, P1, P1Seg, P2, P2Seg, M1, M1Seg, M2, M2Seg, ChrsLength, ideogram, Family, outPath, Seed, dataHap)
# PpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPpPp






