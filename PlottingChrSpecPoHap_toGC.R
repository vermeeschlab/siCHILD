
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
#source(paste(siCHILD_DIR,"ChrSpecbafplotwithpohap.R",sep=""))

source(paste(siCHILD_DIR,"chrspecbafplot.R",sep=""))

load(paste(siCHILD_DIR,"Ideogram_hg19.rda",sep=""))
Chroms = c(1:22,"X","Y")
#Chroms = 12

Params <- read.table(parametersFile,sep="\t",header=T,stringsAsFactors=F)

for(i in 1:nrow(Params)) { if(Params[i,"Param"]=="Parent" | Params[i,"Param"]=="Seed" | Params[i,"Param"]=="GC_File" | Params[i,"Param"]=="Condition") {eval(parse(text=paste(Params[i,"Param"],"='",Params[i,"Value"],"'",sep="")))} 
			   else {eval(parse(text=paste(Params[i,"Param"],"=",Params[i,"Value"])))}
			 }

if(ExcInt==1){Int <- read.table(paste(PGD_EXPORTED_DIR,Family,"_Intervals.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)}

dataPath <- outPath
#outPath <- paste(outPath,"ChrSpecPlots/",sep="")


if (!file.exists(outPath)){
	dir.create(outPath)
}

P1 <- read.table(paste(dataPath,"P1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P2 <- read.table(paste(dataPath,"P2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M1 <- read.table(paste(dataPath,"M1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M2 <- read.table(paste(dataPath,"M2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P1Seg <- read.table(paste(dataPath,"P1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P2Seg <- read.table(paste(dataPath,"P2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M1Seg <- read.table(paste(dataPath,"M1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M2Seg <- read.table(paste(dataPath,"M2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
BAFs <- read.table(paste(dataPath,Family,"_BAFs.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
logRs <- read.table(paste(dataPath,Family,"_logRsAvgWindow.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
logRsSeg <- read.table(paste(dataPath,Family,"_gammaSc300_gammaMc50_logRseg.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
dataPo <- read.table(paste(dataPath,Family,".poo",sep=""),header=T,sep="\t",stringsAsFactors=F)


centro <- read.table(paste(dataPath,Family,"_centro.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)

if (!is.na( (grep("Uncle", Seed) == 1) || (grep("Aunt", Seed) == 1) )) { 
	IBS1 <- read.table(paste(dataPath,Family,"_IBS1.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	IBS0 <- read.table(paste(dataPath,Family,"_IBS0.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)	
	IBS2 <- read.table(paste(dataPath,Family,"_IBS2.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	IBS1 <- IBS1[order(IBS1[,1],IBS1[,2]),] 
	IBS0 <- IBS0[order(IBS0[,1],IBS0[,2]),] 
	IBS2 <- IBS2[order(IBS2[,1],IBS2[,2]),] 
	
	olaps2 <- read.table(paste(dataPath,Family,"_olaps2.txt",sep=""),header=T,stringsAsFactors=F,sep="\t")[,1]
	olaps1 <- read.table(paste(dataPath,Family,"_olaps1.txt",sep=""),header=T,stringsAsFactors=F,sep="\t")[,1]
	olaps0 <- read.table(paste(dataPath,Family,"_olaps0.txt",sep=""),header=T,stringsAsFactors=F,sep="\t")[,1]
	NewSeed <- "IBSseed"
} else {
	IBS1 <- 99999
	IBS0 <- 99999
	olaps1 <- 99999
	olaps0 <- 99999
	NewSeed <- c()
} 


#ADIO <- read.table(paste(dataPath,"ADI_ADO_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)[,1]
dataPath <- paste0(dataPath,"/",NewSeed,"/")

dataHap <- read.table(paste(dataPath,Family,"_Itp2.hap",sep=""),header=T,sep="\t",stringsAsFactors=F)
dataHapRaw <- read.table(paste(dataPath,Family,"_Raw.hap",sep=""),header=T,sep="\t",stringsAsFactors=F)
#print(paste(dataPath,Family,"_Itp2.hap",sep=""))



chrspecbafplot(dataHap,dataHapRaw,dataPo,BAFs,logRs,logRsSeg,P1,P1Seg,P2,P2Seg,M1,M1Seg,M2,M2Seg,ChrsLength,ideogram,Family,outPath,Chroms,Int,Seed,IBS2,IBS1,IBS0,olaps1,olaps0,centro,NewSeed)




