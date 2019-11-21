#outPath <- "/uz/data/Admin/cme_genome_calculated/"
#Family = "PGD047_C1"

#Int <- rbind(c(0,0,0,"Pat"),c(7,117105838,117356025,"Mat"))

disthap<-function(Family,Int,outPath){

Haps <- read.table(paste(outPath,"/",Family,"_Itp2.hap",sep=""),header=T,sep="\t",stringsAsFactors =F)



Chroms <- Int[Int[,1]!="0",1]
Scs <- colnames(Haps)[grep("^E",colnames(Haps))]

DistFam <- NULL
for(ind in Scs){
	print(ind)
	Dists <- NULL
	for(chr in Chroms){
		
		IntChr <- Int[Int[,1]==chr,]	
		#ChrLength <- as.numeric(ChrsLength[ChrsLength[,1]==paste("chr",chr,sep=""),2])
		HapChr <- Haps[Haps[,"Chr"]==chr,c("Position",ind)]	
		UP <- HapChr[HapChr$Position<=as.numeric(IntChr[2]),]
		if (nrow(UP) <1) {
			UP[1,1] <- NA; 
			UP[1,2] <- NA; 		
		}
		Down <- HapChr[HapChr$Position>=as.numeric(IntChr[3]),]
		if (nrow(Down) <1) {
			Down[1,1] <- NA; 
			Down[1,2] <- NA; 		
		}
		UPRle <- rle(UP[,ind])
		DownRle <- rle(Down[,ind])

		#if(UPRle$values[length(UPRle$values)]==DownRle$values[1]){

		revUP <- UP[seq(dim(UP)[1],1),]
		bpUP <- sort( revUP[c(1,UPRle$lengths[length(UPRle$values)]),1])
		bpDown <- Down[c(1,DownRle$lengths[1]),1]
		Dist <- cbind(ind,chr,UPRle$lengths[length(UPRle$values)],DownRle$lengths[1],(bpUP[2]-bpUP[1]), (bpDown[2]-bpDown[1]), UPRle$values[length(UPRle$values)], DownRle$values[1])#}	

		Dists<- rbind(Dists,Dist)
	
	}#end chr
	
	DistFam <-rbind(DistFam,Dists)
	print(ind)
		
}#end ind loop

colnames(DistFam) <- c("Haplotype","Chr","aSNPup","aSNPdown","LengthUp","LengthDown","HapUp","HapDown")
outfile <- paste(outPath,"DistanceToHR_",Family,".txt",sep="")
print(outfile)
write.table(DistFam,outfile,quote=F,sep="\t",col.names=T,row.names=F)
}
