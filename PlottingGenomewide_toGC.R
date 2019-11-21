
args <- commandArgs(TRUE)
Family <- args[1]
Seed <- args[2]
version <- args[3]
SibPattern <- args[4]	# "GM*"
PGD_dir <- args[5]



parametersFile <- paste0(PGD_dir,Family,"/",Family,"_Parameters_",Seed,".txt")
outPath <- paste0(PGD_dir,Family,"_out_",Seed,"_",version,"/")


#PGD_EXPORTED_DIR <- paste("/uz/data/Admin/cme_genome_raw/PGD_Exported/",Family,"/",sep="")
PGD_EXPORTED_DIR <- paste(PGD_dir,"/",Family,"/",sep="")
if (!grepl('/$',outPath)) {outPath <- paste(outPath,"/",sep="")}

siCHILD_DIR <- "/home/jding0/siCHILD_toGC/"
source(paste(siCHILD_DIR,"genomebafplotsc.R",sep=""))
load(paste(siCHILD_DIR,"Ideogram_hg19.rda",sep=""))
Chroms = c(1:22,"X","Y")

Params <- read.table(parametersFile,sep="\t",header=T,stringsAsFactors=F)

for(i in 1:nrow(Params)) { if(Params[i,"Param"]=="Parent" | Params[i,"Param"]=="Seed"  | Params[i,"Param"]=="GC_File" | Params[i,"Param"]=="Condition") {eval(parse(text=paste(Params[i,"Param"],"='",Params[i,"Value"],"'",sep="")))} 
			   else {eval(parse(text=paste(Params[i,"Param"],"=",Params[i,"Value"])))}
			 }

if(ExcInt==1){Int <- read.table(paste(PGD_EXPORTED_DIR,Family,"_Intervals.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)}


outPath <- paste(outPath,"GenomewidePlots/",sep="")

if (!file.exists(outPath)){
  dir.create(outPath)
}



#if (is.null(dataPath)) {
	dataPath <- outPath
#} else {
#	NewSeed <- "IBSseed"
#	dataPath <- paste0(outPath,"/",NewSeed,"/")
#}


P1 <- read.table(paste(dataPath,"P1_",Family,".txt",sep=""),header=T,sep="\t",,stringsAsFactors=F)
P2 <- read.table(paste(dataPath,"P2_",Family,".txt",sep=""),header=T,sep="\t",,stringsAsFactors=F)
M1 <- read.table(paste(dataPath,"M1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M2 <- read.table(paste(dataPath,"M2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P1Seg <- read.table(paste(dataPath,"P1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P2Seg <- read.table(paste(dataPath,"P2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M1Seg <- read.table(paste(dataPath,"M1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M2Seg <- read.table(paste(dataPath,"M2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
BAFs <- read.table(paste(dataPath,Family,"_BAFs.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
logRs <- read.table(paste(dataPath,Family,"_logRsAvgWindow.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
logRsSeg <- read.table(paste(dataPath,Family,"_gammaSc300_gammaMc50_logRseg.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
#dataHap <- read.table(paste(dataPath,Family,"_Itp.hap",sep=""),header=T,sep="\t",stringsAsFactors=F)

genomebafplot(BAFs,logRs,logRsSeg,P1,P1Seg,P2,P2Seg,M1,M1Seg,M2,M2Seg,ChrsLength,ideogram,Family,outPath,Seed,dataHap)


