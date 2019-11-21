chrspecbafplot <- function(dataHap,dataHapRaw,dataPo,BAFs,logRs,logRsSeg,P1,P1Seg,P2,P2Seg,M1,M1Seg,M2,M2Seg,CumLengths,ideogram,Family,outPath,Chroms,Int,Seed,IBS2, IBS1,IBS0,olaps1,olaps0,centro, NewSeed){


source(paste(siCHILD_DIR,"PlotCyto.R",sep=""))

require("plotrix")
Lab <- 2.8
Main <- 6
Ax <- 2.8
C <- 2

#if ((NewSeed %in% c("HalfSib"))) {
	inds <- colnames(BAFs)[c(grep("^Sibling",colnames(BAFs)),grep("^E",colnames(BAFs)))]
#	NewSeed <- NULL
#} else {
#	inds <- colnames(BAFs)[c(grep("E",colnames(BAFs)))]
#}

print(inds)

Seed <- strsplit(Seed,"_")[[1]][1]
if ((Seed %in% c("Grandparents","Sibling"))) {
	colBar <- c("blue","cornflowerblue","red","pink")
} else {
	colBar <- c("dark green","light green","orange","yellow")
}


## ploting color bar base on Condition of seeds. (stepparents & uncle/aunt options)
#if (Condition == "Sick") {
#	hapCond <- c(1,2)
#} else if (Condition == "Healthy") {
#	hapCond <- c(2,1)	
#} else {
#	hapCond <- c(1,2)
#}

hapCond <- c(1,2)



olaps <- c(olaps1)


#if (Condition == "Sick") {
#	olaps <- olaps1
#} else if (Condition == "Healthy") {
#	olaps <- c(olaps0,olaps1)
#}

##########################################
## plot IBS0 -> grey thin line
## plot IBS2 -> black thin line
##########################################


##########################################
##### for Uncle/Aunt cases, add both #####
##########################################
dataPath <- paste0(outPath,"/",NewSeed,"/")
if (length(NewSeed) > 0) {	
	P1SegNewSeed <- read.table(paste(dataPath,"P1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	P2SegNewSeed <- read.table(paste(dataPath,"P2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	M1SegNewSeed <- read.table(paste(dataPath,"M1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	M2SegNewSeed <- read.table(paste(dataPath,"M2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
#	HapCompare <- read.table(paste(outPath,Family,"_HapCompare_after.hap",sep=""),header=T,sep="\t",stringsAsFactors=F)
}




for(ind in inds){
	print(ind)

#	print(paste0("====",nrow(dataHap)))

	for(Chr in Chroms){
		print(Chr)
		chr = Chr
		Time <- gsub(" ","_",gsub(":","",gsub("-","",Sys.time())))
	dir.create(paste(outPath,"/ChrPlots/",sep=""), showWarnings = FALSE)
	fn <- paste(outPath,"/ChrPlots/",Family,"_",ind,"_",Time,"_Chr",Chr,sep="")
#	bitmap(fn,type="jpeg",width=20000,height=18000,units="px")
	#par(omi=c(0,0,0,0), mgp=c(0,0,0))
	pdf(fn,width=20000/600,height=18000/600)
	#	jpeg(jpeg_fn,width=20000,height=18000,res=600)

		#layout(rbind(matrix(1,2,10),matrix(2,3.5,10),matrix(3,1,10),matrix(4,2,10),matrix(5,4,10),matrix(6,4,10),matrix(7,4,10),matrix(8,5,10),matrix(9,4,10),matrix(10,5,10),matrix(11,4,10),matrix(12,5,10),matrix(13,2,10)))

		layout(rbind(matrix(1,2,10),matrix(2,3.5,10),matrix(3,1,10),matrix(4,2,10),matrix(5,4,10),matrix(14,1,10),matrix(15,1.5,10),matrix(6,4,10),matrix(16,1,10),matrix(17,1.5,10),matrix(7,4,10),matrix(8,5,10),matrix(9,4,10),matrix(10,5,10),matrix(11,4,10),matrix(12,5,10),matrix(13,2,10)))

		ChrLength <- as.numeric(ChrsLength[ChrsLength[,1]==paste("chr",Chr,sep=""),2])

		par(mar=c(0,6,2,1))

		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))

		text((ChrLength/2),0, paste("Chromosome",Chr," (",ind,Seed,")"),cex=Main,font=2)

		par(mar=c(2.5,6,1.5,1))

#par(mar=c(0,0,0,0))
		plotCyto(Chr,BandName=T,ChrsLength,ideogram)
		#plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))

		par(mar=c(0,6,0,1))

		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))

# dataPo #
		par(mar=c(0.7,6,0,1))

		plot(dataPo$Position[dataPo$Chr == Chr & dataPo[,ind]>0], dataPo[dataPo$Chr == Chr& dataPo[,ind]>0,ind], "h",xaxt="n",yaxt="n",xlab="", xlim = c(0,ChrLength),ylab="PO", col = "pink",frame=F, ylim = c(-1,1),cex.lab=Lab,cex.main=Main,cex.axis=Ax)                
		points(dataPo$Position[dataPo$Chr == Chr & dataPo[,ind]<0], dataPo[dataPo$Chr == Chr& dataPo[,ind]<0,ind], "h", col = "blue")
		points(dataPo$Position[dataPo$Chr == Chr],rep(0,sum(dataPo[,2] == Chr)), "p", col = "grey", pch = 8,cex=0.1)
		abline(h=seq(-1,1,by=0.5),lty=2,col="grey")
		axis(side=2,at=c(-1,0,1),labels=c("-1","0","1"),cex.axis=Ax)
		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}

# BAF #
		C=2
		par(mar=c(0,6,0,1))

		plot(0,ylab="BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1),col="grey")
		points(BAFs[BAFs[,"Chr"]==Chr,"Position"],BAFs[BAFs[,"Chr"]==Chr,ind],pch=20,col="#00000050",cex=C)	# black
		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}	# orange

# Pat-BAF
		par(mar=c(0,6,0,1))

    		plot(0,ylab="Pat-BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1),col="grey")
	if (length(dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Pat",sep="")] %in% c(1,2)]) > 0) {
    		points(P1[P1[,"Chr"]==Chr,"Position"],P1[P1[,"Chr"]==Chr,ind],pch=20,col="#0000ff10",cex=C)	# blue
		points(P2[P2[,"Chr"]==Chr,"Position"],P2[P2[,"Chr"]==Chr,ind],pch=20,col="#ff000010",cex=C)	# red
		points(P1Seg[P1Seg[,"Chr"]==Chr,"Position"],P1Seg[P1Seg[,"Chr"]==Chr,ind],pch=15,col="#0000ff80",cex=C)
		points(P2Seg[P2Seg[,"Chr"]==Chr,"Position"],P2Seg[P2Seg[,"Chr"]==Chr,ind],pch=15,col="#ff000080",cex=C)

		if (!is.na( (grep("Uncle", Seed) == 1) || (grep("Aunt", Seed) == 1) )) {
		points(P1SegNewSeed[P1SegNewSeed[,"Chr"]==Chr & P1SegNewSeed$Name %in% olaps1,"Position"],P1SegNewSeed[P1SegNewSeed[,"Chr"]==Chr & P1SegNewSeed$Name %in% olaps1 ,ind],pch=4,col="orange",cex=C)
		points(P2SegNewSeed[P2SegNewSeed[,"Chr"]==Chr & P2SegNewSeed$Name %in% olaps1,"Position"],P2SegNewSeed[P2SegNewSeed[,"Chr"]==Chr & P2SegNewSeed$Name %in% olaps1,ind],pch=4,col="light green",cex=C)
		}

		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}
	}	

 
# Mat-BAF #
		par(mar=c(0,6,0,1))

    		plot(0,ylab="Mat-BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1),col="grey")
	if (length(dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Mat",sep="")] %in% c(1,2)]) > 0) {		
    		points(M1[M1[,"Chr"]==Chr,"Position"],M1[M1[,"Chr"]==Chr,ind],pch=20,col="#0000ff10",cex=C)
		points(M2[M2[,"Chr"]==Chr,"Position"],M2[M2[,"Chr"]==Chr,ind],pch=20,col="#ff000010",cex=C)
		points(M1Seg[M1Seg[,"Chr"]==Chr,"Position"],M1Seg[M1Seg[,"Chr"]==Chr,ind],pch=15,col="#0000ff80",cex=C)
		points(M2Seg[M2Seg[,"Chr"]==Chr,"Position"],M2Seg[M2Seg[,"Chr"]==Chr,ind],pch=15,col="#ff000080",cex=C)

		if (!is.na( (grep("Uncle", Seed) == 1) || (grep("Aunt", Seed) == 1) )) {
		points(M1SegNewSeed[M1SegNewSeed[,"Chr"]==Chr & M1SegNewSeed$Name %in% olaps1,"Position"],M1SegNewSeed[M1SegNewSeed[,"Chr"]==Chr & M1SegNewSeed$Name %in% olaps1,ind],pch=4,col="orange",cex=C)
		points(M2SegNewSeed[M2SegNewSeed[,"Chr"]==Chr & M2SegNewSeed$Name %in% olaps1,"Position"],M2SegNewSeed[M2SegNewSeed[,"Chr"]==Chr & M2SegNewSeed$Name %in% olaps1,ind],pch=4,col="light green",cex=C)
		}

		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}
	}	
	
# logR #
		par(mar=c(0,6,0,1))

    		plot(0,ylab="logR",col="white",main="",xlab="",frame=F,ylim=c(-1.5,1.5),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5),labels=c("","-2","","-1","","0","","1","","2",""),cex.axis=Ax)
    		abline(h=seq(-2,2,0.5),lty=2,ylim=c(0,1),col="grey")

    		points(logRs[logRs[,"Chr"]==Chr,"Position"],logRs[logRs[,"Chr"]==Chr,ind],pch=20,col="#00000030",cex=C)
    		points(logRsSeg[logRsSeg[,"Chr"]==Chr,"Position"],logRsSeg[logRsSeg[,"Chr"]==Chr,ind],pch=20,col="#ff000080",cex=C)
		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}

		par(mar=c(0,6,0,1))

    		plot(0,ylab="Father-BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1),col="grey")

    		points(BAFs[BAFs[,"Chr"]==Chr,"Position"],BAFs[BAFs[,"Chr"]==Chr,paste("Father_",Family,sep="")],pch=20,col="#00000030",cex=C)
		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}	
		par(mar=c(0,6,0,1))

    		plot(0,ylab="Father-logR",col="white",main="",xlab="",frame=F,ylim=c(-1.5,1.5),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5),labels=c("","-2","","-1","","0","","1","","2",""),cex.axis=Ax)
    		abline(h=seq(-2,2,0.5),lty=2,ylim=c(0,1),col="grey")

    		points(logRs[logRs[,"Chr"]==Chr,"Position"],logRs[logRs[,"Chr"]==Chr,grep("Father_",colnames(logRs))],pch=20,col="#00000030",cex=C)
    		points(logRsSeg[logRsSeg[,"Chr"]==Chr,"Position"],logRsSeg[logRsSeg[,"Chr"]==Chr,grep("Father_",colnames(logRs))],pch=20,col="#ff000080",cex=C)
		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}
		par(mar=c(0,6,0,1))
    		plot(0,ylab="Mother-BAF",col="white",main="",xlab="Position (Gb)",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,yaxt="n",xaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
		title(xlab="Position (Gb)")	
    		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

    		points(BAFs[BAFs[,"Chr"]==Chr,"Position"],BAFs[BAFs[,"Chr"]==Chr,paste("Mother_",Family,sep="")],pch=20,col="#00000030",cex=C)
		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}
		par(mar=c(0,6,0,1))
    		plot(0,ylab="Mother-logR",col="white",main="",xlab="",frame=F,ylim=c(-1.5,1.5),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(-1.5,-1,-0.5,0,0.5,1,1.5),labels=c("","-1","","0","","1",""),cex.axis=Ax)
    		abline(h=seq(-1,1,0.5),lty=2,ylim=c(0,1))

    		points(logRs[logRs[,"Chr"]==Chr,"Position"],logRs[logRs[,"Chr"]==Chr,grep("Mother_",colnames(logRs))],pch=20,col="#00000030",cex=C)
    		points(logRsSeg[logRsSeg[,"Chr"]==Chr,"Position"],logRsSeg[logRsSeg[,"Chr"]==Chr,grep("Mother_",colnames(logRs))],pch=20,col="#ff000080",cex=C)

  		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
		XaxPos = round(seq(0,ChrLength,round(ChrLength/4))/1000000,digits=2)
       		axis(side=1,at=XaxPos*1000000,labels=as.character(XaxPos),cex.axis=Ax)
		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}

		par(mar=c(0.5,6,2,1))
		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
		text((ChrLength/2),0, paste("Position (Mb)"),cex=Ax,font=1)



# hyplotype block - father #

	if (!is.na( (grep("Uncle", Seed) == 1) || (grep("Aunt", Seed) == 1) )) {
		par(mar=c(1,6,1,1))
		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
		abline(v=dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Pat",sep="")]==hapCond[1] & dataHapRaw$Name %in% olaps], col = colBar[1])	# dark
		abline(v=dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Pat",sep="")]==hapCond[2] & dataHapRaw$Name %in% olaps], col = colBar[2])	# light
		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}

	
#		points(dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Pat",sep="")]==hapCond[1] & dataHapRaw$Name %in% olaps & dataHapRaw$Name %in% ADIO], col = "red")	# dark
#		points(dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Pat",sep="")]==hapCond[2] & dataHapRaw$Name %in% olaps & dataHapRaw$Name %in% ADIO], col = "green")	# light




		par(mar=c(0.3,6,0.3,1))
		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
		abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Pat",sep="")]==hapCond[1] & dataHap$Name %in% olaps], col = colBar[1])
		abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Pat",sep="")]==hapCond[2] & dataHap$Name %in% olaps], col = colBar[2])

	} else {

		par(mar=c(1,6,1,1))
		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
		abline(v=dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Pat",sep="")]==hapCond[1] ], col = colBar[1])	# dark
		abline(v=dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Pat",sep="")]==hapCond[2] ], col = colBar[2])	# light
		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}

		par(mar=c(0.3,6,0.3,1))
		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
		abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Pat",sep="")]==hapCond[1] ], col = colBar[1])
		abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Pat",sep="")]==hapCond[2] ], col = colBar[2])
	}


		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}

		if (IBS1 == 99999 || IBS0 == 99999) {	# grandparents & sibling option
#			print(Seed)
		} else {
		    if (length(dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Pat",sep="")] %in% c(1,2)]) > 0) {

			if (nrow(IBS2[as.character(IBS2$chrom) == chr, ]) > 0) {
#print(IBS1[as.character(IBS1$chrom) == chr, ])
				for (i in 1:nrow(IBS2[as.character(IBS2$chrom) == chr, ])) {
					rect(IBS2[as.character(IBS2$chrom)==chr,][i,]$start, -1.09,IBS2[as.character(IBS2$chrom)==chr,][i,]$end, -1, col = "black", border = "black")
#print(IBS1[as.character(IBS1$chrom)==chr,][i,]$start)
#print(IBS1[as.character(IBS1$chrom)==chr,][i,]$end)
				}
			}

			if (nrow(IBS0[as.character(IBS0$chrom) == chr, ]) > 0) {
				for (i in 1:nrow(IBS0[as.character(IBS0$chrom) == chr, ])) {
					rect(IBS0[as.character(IBS0$chrom)==chr,][i,]$start, -1.09,IBS0[as.character(IBS0$chrom)==chr,][i,]$end, -1, col = "grey", border = "grey")
				}
			}
				
			for (i in 1:nrow(centro[as.character(centro[,1])==chr,])) {
				rect(centro[as.character(centro[,1])==chr,][i,]$Start, -1.09,centro[as.character(centro[,1])==chr,][i,]$Stop, -1, col = "white", border = "white")
			}

#			## plot dataHapRaw & dataHap conflict cases!
#			if (nrow(HapCompare[as.character(HapCompare$Chr)==chr & HapCompare[,ind]>0,]) > 0 ) {
#				rect(HapCompare[as.character(HapCompare$Chr)==chr & HapCompare[,ind]>0,]$Position, 0.95,HapCompare[as.character(HapCompare$Chr)==chr & HapCompare[,ind]>0,]$Position+1, 1.01, col = "red", border = "red")
#			}

		    }
		}



# hyplotype block - mother #

	## plot block
	if (!is.na( (grep("Uncle", Seed) == 1) || (grep("Aunt", Seed) == 1) )) {
		par(mar=c(1,6,1,1))
		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
		abline(v=dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Mat",sep="")]==hapCond[1] & dataHapRaw$Name %in% olaps], col = colBar[3])
		abline(v=dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Mat",sep="")]==hapCond[2] & dataHapRaw$Name %in% olaps], col = colBar[4])
		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}

#		points(dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Mat",sep="")]==hapCond[1] & dataHapRaw$Name %in% olaps & dataHapRaw$Name %in% ADIO], col = "red")
#		points(dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Mat",sep="")]==hapCond[2] & dataHapRaw$Name %in% olaps & dataHapRaw$Name %in% ADIO], col = "green")

		par(mar=c(0.3,6,0.3,1))
		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
		abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Mat",sep="")]==hapCond[1] & dataHap$Name %in% olaps], col = colBar[3])
		abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Mat",sep="")]==hapCond[2] & dataHap$Name %in% olaps], col = colBar[4])
	} else {
		par(mar=c(1,6,1,1))
		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
		abline(v=dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Mat",sep="")]==hapCond[1]], col = colBar[3])
		abline(v=dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Mat",sep="")]==hapCond[2]], col = colBar[4])
		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}

		par(mar=c(0.3,6,0.3,1))
		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
		abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Mat",sep="")]==hapCond[1]], col = colBar[3])
		abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Mat",sep="")]==hapCond[2]], col = colBar[4])
	}

	
	## plot IBS
		if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=5,lwd=4,col="#ff800040")}


		if (IBS1 == 99999 || IBS0 == 99999) {	# grandparents & sibling option
#			print(Seed)
		} else {
		    if (length(dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Mat",sep="")] %in% c(1,2)]) > 0) {
#print(IBS1[as.character(IBS1$chrom) == chr, ])
			if (nrow(IBS2[as.character(IBS2$chrom) == chr, ]) > 0) {
				for (i in 1:nrow(IBS2[as.character(IBS2$chrom) == chr, ])) {
					rect(IBS2[as.character(IBS2$chrom)==chr,][i,]$start, -1.09,IBS2[as.character(IBS2$chrom)==chr,][i,]$end, -1, col = "black", border = "black")
#print(IBS1[as.character(IBS1$chrom)==chr,][i,]$start)
#print(IBS1[as.character(IBS1$chrom)==chr,][i,]$end)
				}
			}

			if (nrow(IBS0[as.character(IBS0$chrom) == chr, ]) > 0) {
				for (i in 1:nrow(IBS0[as.character(IBS0$chrom) == chr, ])) {
					rect(IBS0[as.character(IBS0$chrom)==chr,][i,]$start, -1.09,IBS0[as.character(IBS0$chrom)==chr,][i,]$end, -1, col = "grey", border = "grey")
				}
			}
				
			for (i in 1:nrow(centro[as.character(centro[,1])==chr,])) {
				rect(centro[as.character(centro[,1])==chr,][i,]$Start, -1.09,centro[as.character(centro[,1])==chr,][i,]$Stop, -1, col = "white", border = "white")
			}

#			## plot dataHapRaw & dataHap conflict cases!
#			if (nrow(HapCompare[as.character(HapCompare$Chr)==chr & HapCompare[,ind]>0,]) > 0 ) {
#				rect(HapCompare[as.character(HapCompare$Chr)==chr & HapCompare[,ind]>0,]$Position, 0.95,HapCompare[as.character(HapCompare$Chr)==chr & HapCompare[,ind]>0,]$Position+1, 1.01, col = "red", border = "red")
#			}

		    }
		}


#mtext(" 44444 informative markers!", side =4)

		dev.off()

jpeg_fn <- paste(fn,".jpg",sep="")
system(paste0("gs -sDEVICE=jpeg -o ",jpeg_fn," -sPAPERSIZE=a4 ",fn))
system(paste0("rm ", fn))
#system(paste0("/uz/data/hydra/shared_app/apps/imagemagick/7.0.7-27/bin/convert -flatten -density 250 -quality 100 -resize 20% ",fn, " ", jpeg_fn))
print(jpeg_fn)

	}#end chr loop


}#end ind loop

	print(paste("The plots were saved at",outPath))

}#end function
