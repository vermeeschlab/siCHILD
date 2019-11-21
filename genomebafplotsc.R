genomebafplot <- function(BAFs,logRs,logRsSeg,P1,P1Seg,P2,P2Seg,M1,M1Seg,M2,M2Seg,CumLengths,ideogram,Family,outPath,Seed,dataHap){

source(paste(siCHILD_DIR,"plotCytoGenome.R",sep=""))

require("plotrix")
Lab <- 1.2
Main <- 2.2
Ax <- 1.2

CumLengths <- ChrsLength
CumLengths[,"Length"] <- cumsum(as.numeric(ChrsLength[,2]))

GenomeLength <- as.numeric(CumLengths[CumLengths[,"Chromosome"]=="chrX","Length"])

P1 <- na.omit(P1)
P1Seg <- na.omit(P1Seg)

P2 <- na.omit(P2)
P2Seg <- na.omit(P2Seg)

M1 <- na.omit(M1)
M1Seg <- na.omit(M1Seg)

M2 <- na.omit(M2)
M2Seg <- na.omit(M2Seg)

for(chr in CumLengths[,"Chromosome"][2:nrow(CumLengths)]){


	ToAdd <- as.numeric(CumLengths[grep(paste(chr,"$",sep=""),CumLengths[,"Chromosome"])-1,"Length"])
	
	BAFs[BAFs$Chr==gsub("chr","",chr),"Position"] <- BAFs[BAFs$Chr==gsub("chr","",chr),"Position"] +ToAdd
	logRs[logRs$Chr==gsub("chr","",chr),"Position"] <- logRs[logRs$Chr==gsub("chr","",chr),"Position"] +ToAdd
	logRsSeg[logRsSeg$Chr==gsub("chr","",chr),"Position"] <- logRsSeg[logRsSeg$Chr==gsub("chr","",chr),"Position"] +ToAdd

	P1[P1$Chr==gsub("chr","",chr),"Position"] <- P1[P1$Chr==gsub("chr","",chr),"Position"] +ToAdd
	P1Seg[P1Seg$Chr==gsub("chr","",chr),"Position"] <- P1Seg[P1Seg$Chr==gsub("chr","",chr),"Position"] +ToAdd
	P2[P2$Chr==gsub("chr","",chr),"Position"] <- P2[P2$Chr==gsub("chr","",chr),"Position"] +ToAdd
	P2Seg[P2Seg$Chr==gsub("chr","",chr),"Position"] <- P2Seg[P2Seg$Chr==gsub("chr","",chr),"Position"] +ToAdd

	M1[M1$Chr==gsub("chr","",chr),"Position"] <- M1[M1$Chr==gsub("chr","",chr),"Position"] +ToAdd
	M1Seg[M1Seg$Chr==gsub("chr","",chr),"Position"] <- M1Seg[M1Seg$Chr==gsub("chr","",chr),"Position"] +ToAdd
	M2[M2$Chr==gsub("chr","",chr),"Position"] <- M2[M2$Chr==gsub("chr","",chr),"Position"] +ToAdd
	M2Seg[M2Seg$Chr==gsub("chr","",chr),"Position"] <- M2Seg[M2Seg$Chr==gsub("chr","",chr),"Position"] +ToAdd
#print(chr)
}#end chr loop


Poses <- rep(0,24)
Poses[1] <- (as.numeric(ChrsLength[1,"Length"])/2)

for(i in 2:24){
	
	Poses[i] <- as.numeric(CumLengths[i-1,"Length"])+(as.numeric(ChrsLength[i,"Length"])/2)
}


#for(ind in colnames(BAFs)[-c(1:3,grep("ther",colnames(BAFs)))]){
for(ind in colnames(BAFs)[c(grep("^E",colnames(BAFs)))]){

Time <- gsub(" ","_",gsub(":","",gsub("-","",Sys.time())))

dir.create(paste(outPath,"/GenomePlots/",sep=""), showWarnings = FALSE)
fn <- paste(outPath,"/GenomePlots/",Family,"_",ind,"_",Time,"_GenomeMultiProfile",sep="")
pdf(fn,width=10000/600,height=6500/600)
#bitmap(fn,type="jpeg",width=10000,height=6500,units = "px")
#print(fn)

layout(rbind(matrix(1,2,10),matrix(2,3.5,10),matrix(3,1,10),matrix(4,4,10),matrix(5,4,10),matrix(6,4,10),matrix(7,5,10),matrix(8,4,10),matrix(9,5,10),matrix(10,4,10),matrix(11,5,10),matrix(12,2,10)))

ChrsLength<- ChrsLength[ChrsLength[,"Chromosome"]!="chrY",]
ideogram<- ideogram[ideogram[,"Chromosome"]!="chrY",]

par(mar=c(0,6,2,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,GenomeLength))
text((GenomeLength/2),0, paste("Whole-genome profile (",ind,Seed,")"),cex=Main,font=2)



source(paste(siCHILD_DIR,"plotCytoGenome.R",sep=""))
par(mar=c(2.5,6,1.5,1))
plotCytoGenome(BandName=T,ChrsLength,ideogram)



par(mar=c(0,6,0,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,GenomeLength))

for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){
        rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"]),2,col="#84848430",border="#84848430")
        }# end chr loop

text(Poses,0, c(1:22,"X","Y"),cex=Lab,font=2)

C=0.4
		 par(mar=c(0,6,0,1))
    	plot(0,ylab="BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
	for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){		    	 
	rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],2,col="#84848430",border="#84848430") 
	}# end chr loop
		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

    	points(BAFs[,"Position"],BAFs[,ind],pch=20,col="#00000050",cex=C)	# grey

		 par(mar=c(0,6,0,1))
    	plot(0,ylab="Pat-BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)

	for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){		    	 
	rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],2,col="#84848430",border="#84848430") 
	}# end chr loop
    			abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

#	if (length(dataHap$Position[as.character(dataHap$Chr)==gsub("chr","",chr) & dataHap[,paste(ind,"_Pat",sep="")]==1]) > 0) {
    		points(P1[,"Position"],P1[,ind],pch=20,col="#0000ff10",cex=C)	# very light grey
		points(P2[,"Position"],P2[,ind],pch=20,col="#ff000010",cex=C)
		points(P1Seg[,"Position"],P1Seg[,ind],pch=15,col="#0000ff80",cex=C*1.5)	# blue
		points(P2Seg[,"Position"],P2Seg[,ind],pch=15,col="#ff000080",cex=C*1.5)	# red

		par(mar=c(0,6,0,1))
    		plot(0,ylab="Mat-BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
	
		for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){		    	 
		rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],2,col="#84848430",border="#84848430") 
		}# end chr loop
    		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))
#	}


#	if (length(dataHap$Position[as.character(dataHap$Chr)==gsub("chr","",chr) & dataHap[,paste(ind,"_Mat",sep="")]==1]) > 0) {
    		points(M1[,"Position"],M1[,ind],pch=20,col="#0000ff10",cex=C)
		points(M2[,"Position"],M2[,ind],pch=20,col="#ff000010",cex=C)
		points(M1Seg[,"Position"],M1Seg[,ind],pch=15,col="#0000ff80",cex=C*1.5)
		points(M2Seg[,"Position"],M2Seg[,ind],pch=15,col="#ff000080",cex=C*1.5)

 		par(mar=c(0,6,0,1))
    		plot(0,ylab="logR",col="white",main="",xlab="",frame=F,ylim=c(-3,3),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
        	axis(side=2,at=c(-3,-2,-1,0,1,2,3),labels=c("","-2","","0","","2",""),cex.axis=Ax)
	
		for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){		    	 
		rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-5,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],5,col="#84848430",border="#84848430") 
		}# end chr loop
    		abline(h=seq(-2,2,1),lty=2,ylim=c(0,1))
#	}



    	points(logRs[,"Position"],logRs[,ind],pch=20,col="#00000030",cex=C)	# grey
    	points(logRsSeg[,"Position"],logRsSeg[,ind],pch=20,col="orange",cex=C*1.5)



		par(mar=c(0,6,0,1))
    		plot(0,ylab="Father-BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
	
	for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){		    	 
		rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],2,col="#84848430",border="#84848430") 
		}# end chr loop
    		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

    	points(BAFs[,"Position"],BAFs[,paste("Father_",Family,sep="")],pch=20,col="#00000030",cex=C)
	
		par(mar=c(0,6,0,1))
    	plot(0,ylab="Father-logR",col="white",main="",xlab="",frame=F,ylim=c(-3,3),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
        axis(side=2,at=c(-3,-2,-1,0,1,2,3),labels=c("","-2","","0","","2",""),cex.axis=Ax)

	for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){		    	 
	rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-5,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],5,col="#84848430",border="#84848430") 
	}# end chr loop
    		abline(h=seq(-2,2,1),lty=2,ylim=c(0,1))

    	points(logRs[,"Position"],logRs[,grep("Father_",colnames(logRs))],pch=20,col="#00000030",cex=C)
    	points(logRsSeg[,"Position"],logRsSeg[,grep("Father_",colnames(logRs))],pch=20,col="orange",cex=C*1.5)


		 par(mar=c(0,6,0,1))
    	plot(0,ylab="Mother-BAF",col="white",main="",xlab="Position (Gb)",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,yaxt="n",xaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
		#axis(side=1,at=seq(0,3000000000,500000000),labels=c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),cex.axis=Ax,xlab="Position (Gb)")	
		title(xlab="Position (Gb)")	
	for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){		    	 
	rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-5,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],5,col="#84848430",border="#84848430") 
	}# end chr loop
    		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

    	points(BAFs[,"Position"],BAFs[,paste("Mother_",Family,sep="")],pch=20,col="#00000030",cex=C)


		par(mar=c(0,6,0,1))
    	plot(0,ylab="Mother-logR",col="white",main="",xlab="",frame=F,ylim=c(-3,3),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
	axis(side=2,at=c(-3,-2,-1,0,1,2,3),labels=c("","-2","","0","","2",""),cex.axis=Ax)

	for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){		    	 
	rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-5,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],5,col="#84848430",border="#84848430") 
	}# end chr loop
    		abline(h=seq(-2,2,1),lty=2,ylim=c(0,1))

    	points(logRs[,"Position"],logRs[,grep("Mother_",colnames(logRs))],pch=20,col="#00000030",cex=C)
    	points(logRsSeg[,"Position"],logRsSeg[,grep("Mother_",colnames(logRs))],pch=20,col="orange",cex=C*1.5)

        axis(side=1,at=seq(0,3000000000,500000000),labels=c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),cex.axis=Ax,xlab="Position (Gb)")

par(mar=c(0.5,6,2,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,GenomeLength))
text((GenomeLength/2),0, paste("Position (Gb)"),cex=Ax,font=1)

dev.off()


#jpeg_fn <- paste(fn,".jpg",sep="")
#system(paste0("gs -sDEVICE=jpeg -o ",jpeg_fn," -sPAPERSIZE=a4 ",fn))
#system(paste0("rm ", fn))
# print(jpeg_fn)

}#end ind loop

print(paste("The plots were saved at",outPath))

}#end function
