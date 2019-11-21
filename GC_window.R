load("/.../Ideogram_hg19.rda")

Name <- ""
Window <- 10000
outPath <- ""

SnpAnotFile <- as.data.frame(Gtype[,c("Name","Chr","Position")],stringsAsFactors=F)

SnpAnotFile$Position <- as.numeric(as.character(SnpAnotFile$Position))
SnpAnotFile$Chr <- as.character(SnpAnotFile$Chr)

Interval <- cbind(SnpAnotFile$Position-round(Window/2),SnpAnotFile$Position+round(Window/2))

#Bed <- cbind(paste("chr",SnpAnotFile$Chr,sep=""),Interval)
#Bed <- cbind(SnpAnotFile$Chr,Interval)
Bed <- cbind(SnpAnotFile$Chr,Interval,as.character(SnpAnotFile$Name))
Bed <- data.frame(Bed)
Bed[,1] <- as.character(Bed[,1])
Bed[,4] <- as.character(Bed[,4])
Bed[,2] <- as.numeric(as.character(Bed[,2]))
Bed[,3] <- as.numeric(as.character(Bed[,3]))
ChrsLength <- data.frame(ChrsLength)
ChrsLength[,1] <- as.character(ChrsLength[,1]) 
ChrsLength[,2] <- as.numeric(as.character(ChrsLength[,2])) 
#Bed <- gsub(24,"Y")
Chroms <- c(1:22,"X","Y")
for(chr in Chroms){

BedChr <- Bed[Bed[,1]==chr,] 
ChrLength <- ChrsLength[ChrsLength[,1]==paste("chr",chr,sep=""),2]  

if(sum(BedChr[,2]<0)>=1){BedChr[BedChr[,2]<0,2]<-0}
if(BedChr[nrow(BedChr),3]>ChrLength){BedChr[nrow(BedChr),3]<-ChrLength}
if(BedChr[nrow(BedChr),2]>ChrLength){BedChr[nrow(BedChr),2]<-ChrLength}

if(chr==Chroms[1]){BedGenome<-BedChr}else{BedGenome<-rbind(BedGenome,BedChr)}

}#end chr loop

BedGenome[,1] <- paste0("chr",BedGenome[,1],sep="") 

write.table(BedGenome,paste(outPath,Name,"_Window",Window,".bed",sep=""),sep="\t",quote=F,col.names=F,row.names=F)


nucBed -fi /.../ucsc.hg19.fasta -bed /…/Family_Window10000.bed > /../Seqstat_Family_Window10000.txt
