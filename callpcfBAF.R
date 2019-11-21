callpcfBAF <- function(BAFph,gammaSC,gammaMC,plateau){

	source(paste(siCHILD_DIR,"fastPCF.R",sep=""))
	print("PCF segmentation is applying...")

	informative <- data.frame(Embryo="", Chr="", Nr="", stringsAsFactors=FALSE)
	j=0

	SegGenomes <- NULL
	SegGenomeIntervals <- NULL
	SegGenomeLens <- NULL


 	for(i in 4:ncol(BAFph)){

		if(sum(grep("^E",colnames(BAFph)[i]))==1 | sum(grep("_MC",colnames(BAFph)[i]))==1){gamma=gammaSC}else{gamma=gammaMC}
		SegGenome <- NULL
		SegGenomeInterval <- NULL
		SegGenomeLen <- NULL

		for(chr in unique(BAFph[,"Chr"])){
#print(chr)
   			j=j+1
			informative[j,] <- c(colnames(BAFph)[i], chr, -999)
			logRChr <- BAFph[as.character(BAFph$Chr)==chr,i]
	
			#if(sum(is.na(logRChr))>=1){warning(paste("There are",sum(is.na(logRChr),"missing values in logR-values of Chr.",chr))}

			while(sum(is.na(logRChr))>=1 && length(logRChr) >10){logRChr[which(is.na(logRChr))]<-logRChr[which(is.na(logRChr))-1]}
			informative[j,] <- c(colnames(BAFph)[i], chr, length(logRChr))

			sdev <- getMad(logRChr,k=plateau)
			res <- try(selectFastPcf(logRChr,3,gamma*sdev,T),silent=T)

#			SegLen <- res$Lengde

			if(class(res) == "try-error") {
				warning(paste0(length(logRChr), " Too little informative markers on Chr", chr, ", Gamma=",gamma*sdev, ", ",colnames(BAFph)[i] ))
				SegChr <- logRChr
				SegInterval <- NA
			} else {
				SegChr <- res$yhat
				SegInterval <- res$nIntervals
			}

#			if(chr==1){  
#				SegGenome <- SegChr; 
#				SegGenomeInterval <- SegInterval; 
##				SegGenomeLen <- SegLen;  
#			}else{
				SegGenome <- c(SegGenome,SegChr); 
				SegGenomeInterval <- c(SegGenomeInterval,SegInterval); 
#				SegGenomeLen <- c(SegGenomeLen,SegLen);
#			}
		}#end chr loop
    

#		if(i==4){
#			SegGenomes <- SegGenome;  
#			SegGenomeIntervals <- SegGenomeInterval;  
##			SegGenomeLens <- SegGenomeLen;
#		}else{
			SegGenomes<-cbind(SegGenomes,SegGenome);  
			SegGenomeIntervals <- cbind(SegGenomeIntervals, SegGenomeInterval);  
#			SegGenomeLens <- cbind(SegGenomeLens,SegGenomeLen);
#		}
		print(paste(colnames(BAFph)[i],"==> gamma",gamma, "is applied"))

	}#end file loop




	write.table(informative,paste(outPath,Family,"_informative.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

#	colnames(SegGenomeIntervals) <- colnames(BAFph)[4:8]
	colnames(SegGenomeIntervals) <- colnames(BAFph)[4:ncol(BAFph)]
	SegLogRs <- cbind(BAFph[,c("Name","Chr","Position")],SegGenomes)
	colnames(SegLogRs)<- colnames(BAFph)
#	print("PCF segmentation is finished ////")

	SegLogRs

#	Seg <- vector("list",2)
#	names(Seg) <- c("SegBAF","SegInterval")

#	Seg
}#end function
