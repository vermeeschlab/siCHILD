ind <- "E01_Bl001_PGD090"
require("plotrix")
Lab <- 2.8
Main <- 6
Ax <- 2.8
siCHILD_DIR <- "/home/jding0/siCHILD_toGC/"
load(paste(siCHILD_DIR,"Ideogram_hg19.rda",sep=""))
print(ind)

Chr=6
C=2


##################################
Family <- "PGD090"
dataPath <- "/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD090_out_Grandparents_v1/"
P1 <- read.table(paste(dataPath,"P1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P2 <- read.table(paste(dataPath,"P2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M1 <- read.table(paste(dataPath,"M1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M2 <- read.table(paste(dataPath,"M2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P1Seg <-  read.table(paste(dataPath,"P1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P2Seg <- read.table(paste(dataPath,"P2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M1Seg <- read.table(paste(dataPath,"M1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M2Seg <- read.table(paste(dataPath,"M2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
BAFs <- read.table(paste(dataPath,Family,"_BAFs.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
logRs <- read.table(paste(dataPath,Family,"_logRsAvgWindow.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
logRsSeg <- read.table(paste(dataPath,Family,"_gammaSc300_gammaMc50_logRseg.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
dataPo <- read.table(paste(dataPath,Family,".poo",sep=""),header=T,sep="\t",stringsAsFactors=F)
dataHap <- read.table(paste(dataPath,Family,"_Itp2.hap",sep=""),header=T,sep="\t",stringsAsFactors=F)
dataHapRaw <- read.table(paste(dataPath,Family,"_Raw.hap",sep=""),header=T,sep="\t",stringsAsFactors=F)

centro <- read.table(paste(dataPath,Family,"_centro.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
##################################




Time <- gsub(" ","_",gsub(":","",gsub("-","",Sys.time())))
jpeg_fn <- "~/tmp2.jpg"
jpeg(jpeg_fn ,width=20000,height=18000,res=600)

print(Chr)
#layout(rbind(matrix(1,2,10),matrix(2,3.5,10),matrix(3,1,10),matrix(4,2,10),matrix(5,4,10),matrix(6,4,10),matrix(7,4,10),matrix(8,5,10),matrix(9,4,10),matrix(10,5,10),matrix(11,4,10),matrix(12,5,10),matrix(13,2,10)))

layout(rbind(matrix(1,2,10),matrix(2,3.5,10),matrix(3,1,10),matrix(4,2,10),matrix(5,4,10),matrix(14,1,10),matrix(15,1.5,10),matrix(6,4,10),matrix(16,1,10),matrix(17,1.5,10),matrix(7,4,10),matrix(8,5,10),matrix(9,4,10),matrix(10,5,10),matrix(11,4,10),matrix(12,5,10),matrix(13,2,10)))

ChrLength <- as.numeric(ChrsLength[ChrsLength[,1]==paste("chr",Chr,sep=""),2])




#Seed <- strsplit(Seed,"_")[[1]][1]
#if ((Seed %in% c("Grandparents","Siblings"))) {
#	colBar <- c("blue","cornflowerblue","red","pink")
#} else {
#	colBar <- c("dark green","light green","orange","yellow")
#}


## ploting color bar base on Condition of seeds. (stepparents & uncle/aunt options)
#if (Condition == "Sick") {
#	hapCond <- c(1,2)
#} else if (Condition == "Healthy") {
#	hapCond <- c(2,1)	
#} else {
#	hapCond <- c(1,2)
#}


#		print(Chr)
#		chr = Chr
#		Time <- gsub(" ","_",gsub(":","",gsub("-","",Sys.time())))

#jpeg_fn <- "~/tmp1.jpg"


#		jpeg(jpeg_fn,width=20000,height=18000,res=600)

#		#layout(rbind(matrix(1,2,10),matrix(2,3.5,10),matrix(3,1,10),matrix(4,2,10),matrix(5,4,10),matrix(6,4,10),matrix(7,4,10),matrix(8,5,10),matrix(9,4,10),matrix(10,5,10),matrix(11,4,10),matrix(12,5,10),matrix(13,2,10)))

#		layout(rbind(matrix(1,2,10),matrix(2,3.5,10),matrix(3,1,10),matrix(4,2,10),matrix(5,4,10),matrix(14,1,10),matrix(15,1.5,10),matrix(6,4,10),matrix(16,1,10),matrix(17,1.5,10),matrix(7,4,10),matrix(8,5,10),matrix(9,4,10),matrix(10,5,10),matrix(11,4,10),matrix(12,5,10),matrix(13,2,10)))

#		ChrLength <- as.numeric(ChrsLength[ChrsLength[,1]==paste("chr",Chr,sep=""),2])




print("Pat_BAF")

		 par(mar=c(0,6,0,1))
    	plot(0,ylab="Pat-BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    			abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1),col="grey")
    	points(P1[P1[,"Chr"]==Chr,"Position"],P1[P1[,"Chr"]==Chr,ind],pch=20,col="#0000ff10",cex=C)
		points(P2[P2[,"Chr"]==Chr,"Position"],P2[P2[,"Chr"]==Chr,ind],pch=20,col="#ff000010",cex=C)
		points(P1Seg[P1Seg[,"Chr"]==Chr,"Position"],P1Seg[P1Seg[,"Chr"]==Chr,ind],pch=15,col="#0000ff80",cex=C)
		points(P2Seg[P2Seg[,"Chr"]==Chr,"Position"],P2Seg[P2Seg[,"Chr"]==Chr,ind],pch=15,col="#ff000080",cex=C)
# if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}


# Pat-BAF
		par(mar=c(0,6,0,1))
    		plot(0,ylab="Pat-BAF?????",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1),col="grey")
		if (length(dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Pat",sep="")]==1]) > 0) {
    		points(P1[P1[,"Chr"]==Chr,"Position"],P1[P1[,"Chr"]==Chr,ind],pch=20,col="#0000ff10",cex=C)	# blue
		points(P2[P2[,"Chr"]==Chr,"Position"],P2[P2[,"Chr"]==Chr,ind],pch=20,col="#ff000010",cex=C)	# red
		points(P1Seg[P1Seg[,"Chr"]==Chr,"Position"],P1Seg[P1Seg[,"Chr"]==Chr,ind],pch=15,col="#0000ff80",cex=C)
		points(P2Seg[P2Seg[,"Chr"]==Chr,"Position"],P2Seg[P2Seg[,"Chr"]==Chr,ind],pch=15,col="#ff000080",cex=C)
#		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}
		}	

 
# Mat-BAF #
		par(mar=c(0,6,0,1))
    		plot(0,ylab="Mat-BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1),col="grey")
		if (length(dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Mat",sep="")]==1]) > 0) {		
    		points(M1[M1[,"Chr"]==Chr,"Position"],M1[M1[,"Chr"]==Chr,ind],pch=20,col="#0000ff10",cex=C)
		points(M2[M2[,"Chr"]==Chr,"Position"],M2[M2[,"Chr"]==Chr,ind],pch=20,col="#ff000010",cex=C)
		points(M1Seg[M1Seg[,"Chr"]==Chr,"Position"],M1Seg[M1Seg[,"Chr"]==Chr,ind],pch=15,col="#0000ff80",cex=C)
		points(M2Seg[M2Seg[,"Chr"]==Chr,"Position"],M2Seg[M2Seg[,"Chr"]==Chr,ind],pch=15,col="#ff000080",cex=C)
#		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}
		}	
	

		dev.off()


# resize jpg file
jpeg_fn_resize <- paste(jpeg_fn,".resize.jpg",sep="")
print(jpeg_fn_resize)
system(paste0("/usr/bin/convert -resize 20% ",jpeg_fn, " ", jpeg_fn_resize))


	
