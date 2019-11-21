pickBestSingleCellbyChr <- function(Gtypes,Seg1,Seg2,Parent1,exc){
	source(paste(siCHILD_DIR,"distcompmedseg2.R",sep=""))

	#	Seg1 <- read.table(paste(dataPath,"P1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	#	Seg2 <- read.table(paste(dataPath,"P2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	#	Seg1 <- read.table(paste(dataPath,"M1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	#	Seg2 <- read.table(paste(dataPath,"M2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)



	s1 <- Seg1[,c(1:3,grep("^E",colnames(Seg1)))]
	s2 <- Seg2[,c(1:3,grep("^E",colnames(Seg2)))]
	#print(s1[1,])
	#print(s2[1,])

	if (!is.null(exc)) {
		s1 <- s1[,-(which(colnames(s1) == exc))]
		s2 <- s2[,-(which(colnames(s2) == exc))]
	}



	#Chroms <- unique(sort(s1[,2]))
	Chroms <- 1:22

	distPat_mean <- s1[0,c(2,4:ncol(s1))]
	pick_singleCell <- data.frame(Chr="",Embryo="", stringsAsFactors=FALSE)

	for(chr in Chroms){
		if ((length(unique(s1[s1[,2] == chr,4])) >2) & (length(unique(s2[s2[,2] == chr,4])) >2)) {
			print(chr)
			distPat <- distcompseg(s1[s1[,2] == chr,],s2[s2[,2] == chr,])
			distPat <- distPat[complete.cases(distPat),]
			if (nrow(distPat) > 1) {
				distPat_mean[chr,1] <- chr
				for (d in 4:ncol(distPat)) {
					distPat_mean[chr,(d-2)] <- abs(mean(density(distPat[,d])$x)-0.5)
				}
				pick_singleCell[chr,1] <- chr
				pick_singleCell[chr,2] <- colnames(distPat_mean)[which(distPat_mean[chr,] == min(as.numeric(distPat_mean[chr,])))]
				# Gtypes[Gtypes[,2] == chr,ncol(Gtypes)] <- Gtypes[Gtypes[,2] == chr,pick_singleCell[chr,2]]
			}

		} else {
			pick_singleCell[chr,1] <- chr
			perChr <- (t(as.data.frame(lapply(ParScore, function(x) mean(x[3, 5])))))
			pick_singleCell[chr,2] <- names(sort(perChr[,1], decreasing = TRUE)[1])
		}
	}

	ParScore <- ParScore_exc
	print(names(ParScore))
	# for chrX, take the embryo with best call rate
	pick_singleCell["X",1] <- "X"
	pick_singleCell["X",2] <- names(which.max(lapply(ParScore, function(x) (x["X","CallRate"]))))
	pick_singleCell["Y",1] <- "Y"
	pick_singleCell["Y",2] <- names(which.max(lapply(ParScore, function(x) (x["X","CallRate"]))))

	fn <- paste(outPath,"PickSingleCell_",Parent1,".txt",sep="")
	write.table(pick_singleCell,fn,quote=F,sep="\t",col.names=T,row.names=F)

	pick_singleCell

}#end Family


