pickBestSingleCellbyChr <- function(dataPath,P1Seg,P2Seg,M1Seg,M2Seg){ 
	source(paste(siCHILD_DIR,"distcompmedseg2.R",sep=""))

#	P1Seg <- read.table(paste(dataPath,"P1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
#	P2Seg <- read.table(paste(dataPath,"P2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
#	M1Seg <- read.table(paste(dataPath,"M1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
#	M2Seg <- read.table(paste(dataPath,"M2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)

	p1 <- P1Seg[,c(1:3,grep("^E",colnames(P1Seg)))]
	p2 <- P2Seg[,c(1:3,grep("^E",colnames(P2Seg)))]
	m1 <- M1Seg[,c(1:3,grep("^E",colnames(M1Seg)))]
	m2 <- M2Seg[,c(1:3,grep("^E",colnames(M2Seg)))]


	Chroms <- c(1:22)

	### Pat ###
	distPat_mean <- p1[0,c(2,4:ncol(p1))]
	pick_singleCell_Pat <- data.frame(Chr="",Embryo="", stringsAsFactors=FALSE)
	if (length(unique(p1[,4])) >2) {
		for(chr in Chroms){
			distPat <- distcompseg(p1[p1[,2] == chr,],p2[p2[,2] == chr,])
			distPat_mean[chr,1] <- chr
			for (d in 4:ncol(distPat)) {
				distPat_mean[chr,(d-2)] <- abs(mean(density(distPat[,d])$x)-0.5)
			}
			pick_singleCell_Pat[chr,1] <- chr
			pick_singleCell_Pat[chr,2] <- colnames(distPat_mean)[which(distPat_mean[chr,] == min(as.numeric(distPat_mean[chr,])))]
		}
	}	


	### Mat ###
	distMat_mean <- p1[0,c(2,4:ncol(p1))]
	pick_singleCell_Mat <- data.frame(Chr="",Embryo="", stringsAsFactors=FALSE)
	if (length(unique(p1[,4])) >2) {
		for(chr in Chroms){
			distMat <- distcompseg(p1[p1[,2] == chr,],p2[p2[,2] == chr,])
			distMat_mean[chr,1] <- chr
			for (d in 4:ncol(distMat)) {
				distMat_mean[chr,(d-2)] <- abs(mean(density(distMat[,d])$x)-0.5)
			}
			pick_singleCell_Mat[chr,1] <- chr
			pick_singleCell_Mat[chr,2] <- colnames(distMat_mean)[which(distMat_mean[chr,] == min(as.numeric(distMat_mean[chr,])))]
		}
	}	


}#end Family


