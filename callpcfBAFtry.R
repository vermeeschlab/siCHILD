callpcfBAFtry <- function(BAFph,gammaSC,gammaMC,plateau){

	source(paste(siCHILD_DIR,"fastPCF.R",sep=""))
	print("PCF segmentation is applying...")

	informative <- data.frame(Embryo="", Chr="", Nr="", stringsAsFactors=FALSE)
	j=0

	SegGenomes <- NULL
	for(i in 4:ncol(BAFph)){
	#print(colnames(BAFph)[i])
		if(sum(grep("^E",colnames(BAFph)[i]))==1 | sum(grep("_MC",colnames(BAFph)[i]))==1){gamma=gammaSC}else{gamma=gammaMC}
		SegGenome <- NULL
		for(chr in (unique(BAFph[,"Chr"]))){
			j=j+1
			informative[j,] <- c(colnames(BAFph)[i], chr, -999)
			logRChr <- BAFph[as.character(BAFph$Chr)==chr,i]

			#if(sum(is.na(logRChr))>=1){warning(paste("There are",sum(is.na(logRChr),"missing values in logR-values of Chr.",chr))}

			while(sum(is.na(logRChr))>=1){logRChr[which(is.na(logRChr))]<-logRChr[which(is.na(logRChr))-1]}
			informative[j,] <- c(colnames(BAFph)[i], chr, length(logRChr))
			sdev <- getMad(logRChr,k=plateau)
			res <- try(selectFastPcf(logRChr,3,gamma*sdev,T),silent=T)
			if(class(res) == "try-error") {
				warning(paste0(length(logRChr), " Too little informative markers on Chr", chr, ", Gamma=",gamma*sdev, ", ",colnames(BAFph)[i] ))
				SegChr <- logRChr
			} else {
				SegChr <- res$yhat
			}

			SegGenome <- c(SegGenome,SegChr)

		}#end chr loop
    
		SegGenomes<-cbind(SegGenomes,SegGenome)
		print(paste(colnames(BAFph)[i],"==> gamma",gamma, "is applied"))

	}#end file loop


	write.table(informative,paste(outPath,Family,Seed,"_informative.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
	SegLogRs <- cbind(BAFph[,c("Name","Chr","Position")],SegGenomes)
	colnames(SegLogRs)<- colnames(BAFph)
	print("PCF segmentation is finished ////")
	SegLogRs
}#end function
