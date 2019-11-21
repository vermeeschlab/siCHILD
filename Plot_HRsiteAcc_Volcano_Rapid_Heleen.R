#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%       				Volcano Like Plot (HR)          			   	  %%%%%%%%%%%%%#
#	
#Author: Jia Ding (modified from MZE)
#
#cDate: 13.March.2017, original version 01/01/2012
#mDtae: 
#
#(->):  Read in reconstructed haplotypes 
#
#(<-): Volcano like plots 
#
#Description: This shows accuracy of determined HR-stites for WGA-Haplotypes before and after interpretation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rm(list=ls(all=T))

#args <- commandArgs(TRUE)
#File1 <- args[1]	# F1: "PGD_Hapmap_out_Grandparents_f0/PGD_Hapmap_Itp.hap"
#File2 <- args[2]	# F2: "PGD_Hapmap_out_Grandparents_f0/PGD_Hapmap_Raw.hap"

dataPath <- "/uz/data/hydra/genomicscore/gcpi_r_kul_joris_vermeesch/jding/PGD/PGD_Hapmap/MC_SC_compare/"
outPath <-  dataPath

File1 <- "HapMap_RepliG_3dayprotocol_GM12887_Itp2.hap"
File2 <- "HapMap_RepliG_3dayprotocol_GM12887_Raw.hap"


#Interval for plotting
numSC <- 12
AccInt <- 500 
			
Chroms <- as.character(1:22)
Pars=c("Pat","Mat")



for(par in Pars){
	print(par)
	dataItp <- read.table(paste(dataPath,File1,sep="") ,sep="\t",header=T,stringsAsFactors=F)
	dataRaw <- read.table(paste(dataPath,File2,sep="") ,sep="\t",header=T,stringsAsFactors=F)

	SNPChrPos <- dataItp[,1:3]
	rownames(dataItp) <- SNPChrPos[,1]
	rownames(dataRaw) <- SNPChrPos[,1]
	#Chroms <- c(1:22,"X")

	dataItp <- dataItp[,c(1:3,grep(par,colnames(dataItp)))]
	dataRaw <- dataRaw[,c(1:3,grep(par,colnames(dataRaw)))]

	SCItp <- as.data.frame(dataItp[,grep("SC",colnames(dataItp))])
	SCRaw <- as.data.frame(dataRaw[,grep("SC",colnames(dataRaw))])
	MCRef <- as.data.frame(dataItp[,grep("MC",colnames(dataItp))])
	colnames(MCRef) <- colnames(dataItp)[grep("MC",colnames(dataItp))]
	dataItp$Chr <- as.character(dataItp$Chr)

	#Determining HR-sites in the MC reference

	mc <- ncol(MCRef)

	for(chr in Chroms){	
		dataChr <- as.data.frame(MCRef[dataItp$Chr==chr,])
		SNPChrPosChr <- dataItp[dataItp$Chr==chr,1:3]
		Blk <- rle(dataChr[,mc])
		Blk1 <- as.data.frame(cbind(Blk$values,Blk$lengths,cumsum(Blk$lengths),matrix(0,length(Blk$values),1)))
		Blk1 <- Blk1[Blk1[,1]!=0,]
		Blk2 <- Blk1[Blk1[,2]>=5,]
		if(is.vector(Blk2)){print(paste("No HRsite for Chromosome",chr))
		} else if (nrow(Blk2) == 1) {print(paste("One haplotype, No HRsite for Chromosome",chr))
		} else {			
			for(i in 1:(nrow(Blk2)-1)){if( Blk2[i,2]>=AccInt & Blk2[i+1,2]>=AccInt ){ Blk2[i,4] <- 1}}#end i loop
			HRsite <- cbind(SNPChrPosChr[Blk2[Blk2[,4]==1,3],2:3],Blk2[Blk2[,4]==1,3])
			rownames(HRsite) <- SNPChrPosChr[Blk2[Blk2[,4]==1,3],1]
			if (chr==1){ HRsites <- HRsite }else{ HRsites <- rbind(HRsites,HRsite) }	
		}
	}#end chr loop

	rownames(dataItp) <- as.character(dataItp$Name)
	HRPos <- dataItp[rownames(HRsites),]  
        
	for(i in 1:nrow(HRsites)){
		HRPos <- which(rownames(dataItp)==rownames(HRsites)[i])
		HRsiteAcc_Ref <- MCRef[c((HRPos-AccInt):(HRPos-1),(HRPos+1):(HRPos+AccInt)),mc]
		HRsiteItp_SC <- SCItp[c((HRPos-AccInt):(HRPos-1),(HRPos+1):(HRPos+AccInt)),grep("SC",colnames(SCItp))]
		HRsiteRaw_SC <- SCRaw[c((HRPos-AccInt):(HRPos-1),(HRPos+1):(HRPos+AccInt)),grep("SC",colnames(SCItp))]
	
		AccItp_SCs <- matrix(NA,ncol(HRsiteItp_SC),AccInt*2)
		AccRaw_SCs <- matrix(NA,ncol(HRsiteItp_SC),AccInt*2)

		for(sc in 1:ncol(HRsiteItp_SC)){		
			AccItp_SCs[sc,] <- HRsiteAcc_Ref == HRsiteItp_SC[,sc] & HRsiteAcc_Ref!=0
			AccRaw_SCs[sc,] <- HRsiteAcc_Ref == HRsiteRaw_SC[,sc] & HRsiteAcc_Ref!=0	
		}#end sc loop
	
		if(i==1){
			HRsiteAccItp_SCs <- AccItp_SCs
			HRsiteAccRaw_SCs <- AccRaw_SCs
		}else{
			HRsiteAccItp_SCs <- rbind(HRsiteAccItp_SCs,AccItp_SCs)
			HRsiteAccRaw_SCs <- rbind(HRsiteAccRaw_SCs,AccRaw_SCs)
		}			
	}#end i loop


	AccItpSC <- (apply(HRsiteAccItp_SCs,2,sum)/nrow(HRsiteAccItp_SCs))*100
	AccRawSC <- (apply(HRsiteAccRaw_SCs,2,sum)/nrow(HRsiteAccItp_SCs))*100

	outfile <- paste(outPath,"HRsiteDetectionAcc_illumina_Rapid_",par,colnames(MCRef)[mc],".pdf",sep="")
	print(outfile)
	pdf(outfile,w=10,h=5,bg="transparent",family="Helvetica",title ="RB plot")
	#par(mfrow=c(1,2))

	plot(1:(AccInt*2),matrix(0,1,AccInt*2),xlab = "SNP Position",ylab="HR-detection accuracy (%)",pch="*",ylim=c(0,100),frame=FALSE,main= "illumina CytoSNP12-v2.1 + SC")
	points(1:1000, AccItpSC,"h", col = "#ff000020")
	points(1:1000, AccRawSC,"h", col = "#ff000080")

	#abline(v=500,lty=2,col="grey",lwd=2)

	dev.off()

    

}#end par loop




