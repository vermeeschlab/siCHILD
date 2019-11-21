htypingOpt3 <- function(Gtypes,dataPo,ParScore,Parent1,SibPattern,Family,outPath,olaps0,olaps1,olaps2,IBS1,IBS0){

	print("******************************************")
	print("****     Option3Htype analysis...     ****")

	ScSexes <- matrix(NA,length(ParScore),2)
	rownames(ScSexes)<-names(ParScore)


	Inds <-names(ParScore)
	print("******************************************")

#rownames(Gtypes) <- Gtypes$Name


	for(sc in names(ParScore)){
		ScSexes[sc,1] <- sc
		if(is.na(ParScore[[sc]]["X","Par"])){print(paste("Sex of ",sc ,"could not be determined!!"))
		}else if(ParScore[[sc]]["X","Par"]==12){ScSexes[sc,2] <- "female";	print(paste(sc,"is female"))
		}else if(ParScore[[sc]]["X","Par"]==22){ScSexes[sc,2] <- "male";	print(paste(sc,"is male"))
		}else if(ParScore[[sc]]["X","Par"]==11){ScSexes[sc,2] <- "female";	print(paste(sc,"is a female but without maternal X chromosome"))
		}else if(ParScore[[sc]]["X","Par"]==0){ScSexes[sc,2]  <- "ND";		print(paste("Chr.X nullisomy of ",sc))
		}else{print(paste("Sex of ",sc ,"could not be determined!!"))}
	}

	Father <- Gtypes[Gtypes[,"Chr"]=="X",paste("Father_",Family,sep="")]

	if(sum(Father=="AB")>=1){
		print("Checking mother/father sample swop!")
		print(paste("On the chr. X  of the father",sum(Father=="AB"),"(",sum(Father=="AB")/sum(Gtypes$Chr=="X"), "%) of the SNP-calls are heterozygous!!"))
		print("These will be treated as NoCalls")
		Father[Father=="AB"] <- "NC"#There could not be htz SNPs on chromosome X of the father
		Gtypes[Gtypes[,"Chr"]=="X",paste("Father_",Family,sep="")] <- Father
	}




	#------------	Phasing chr X	-------------
	HapXMat <- chrxhtypingOpt3(Gtypes,Parent1,ScSexes,Family,olaps0,olaps1,olaps2)


	if(Parent1==paste("Mother_",Family,sep="")){Parent2 = paste("Father_",Family,sep="")}else{Parent2 = paste("Mother_",Family,sep="")}



##############################
###### haplotyping Parent1 ###
##############################

	PhasedPar <- Gtypes[,Parent1]
	seedolap0RN <- c(rownames(Gtypes[((Gtypes[,Seed]=="AB" & Gtypes[,Parent1] == "BB") | (Gtypes[,Seed]=="AB" & Gtypes[,Parent1] == "AA")) & Gtypes$Name %in% olaps0,]))
	if (!is.integer0(grep(Seed,colnames(Gtypes)))) {	
		PhasedPar_RN <- c(rownames(Gtypes[Gtypes[,Seed]=="BB" & Gtypes[,Parent1]=="AB" & Gtypes$Name %in% c(olaps1),]),rownames(Gtypes[Gtypes[,Seed]=="AA" & Gtypes[,Parent1]=="AB" & Gtypes$Name %in% c(olaps1),]), rownames(Gtypes[((Gtypes[,Seed]=="BB" & Gtypes[,Parent1] == "AB")) & Gtypes$Name %in% olaps0,]),rownames(Gtypes[((Gtypes[,Seed]=="AA" & Gtypes[,Parent1] == "AB")) & Gtypes$Name %in% olaps0,]))
		PhasedPar[Gtypes[,Seed]=="BB" & Gtypes[,Parent1]=="AB" & Gtypes$Name %in% olaps1] <- "BA"
		PhasedPar[Gtypes[,Seed]=="AA" & Gtypes[,Parent1]=="AB" & Gtypes$Name %in% olaps1] <- "AB"	

		PhasedPar[((Gtypes[,Seed]=="BB" & Gtypes[,Parent1] == "AB")) & Gtypes$Name %in% olaps0] <- "AB"
		PhasedPar[((Gtypes[,Seed]=="AA" & Gtypes[,Parent1] == "AB")) & Gtypes$Name %in% olaps0] <- "BA"


	}

	Gtypes[,Parent1] <- PhasedPar
	Gtypes[!(rownames(Gtypes) %in% PhasedPar_RN),][,Parent1] <- "NC"
write.table(seedolap0RN,paste(outPath,"seedolap0RN_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)


	Htype <- Gtypes
	Htype$IBS <- NA
	Htype[Htype$Name %in% olaps0,]$IBS <- 0
	Htype[Htype$Name %in% olaps1,]$IBS <- 1
	Htype[Htype$Name %in% olaps2,]$IBS <- 2




	Sibs <- colnames(Gtypes)[grep(SibPattern,colnames(Gtypes))]
	ParHaps <- NULL
	ADIO <- NULL
	for(sib in Sibs){
		Hap1 <- rep(0,nrow(Gtypes))
		#print(sib)

		GtypeChild <- Gtypes[,sib]
	

		## only allow ADI ADO at IBS1 region ##
		#-------------------------------------------------------------------------------------------------------
		## 1st allele is common allele; single cell didn't inherited common allele'. light color; 
		Hap1[(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AA") |	
			(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AB") | 
			(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AA") |	# & Gtypes$Name %in% olaps1) | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="BB") |
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AB") | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="BB")] <- 2


		## 1st allele is common allele; single cell inherited the common allele; dark color; 
		## if common=sick, H1=dark=sick;H2=light=healthy. 
		## if common=healthy, H1=dark=healthy, H2=light=sick. 
		Hap1[(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AB") |	
			(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="BB") | 
			(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="BB") | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AB") | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AA") | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AA")] <- 1	
		#----------------------------------------------------------------------------------------------------------

ADIO <- c(ADIO,Gtypes[(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AA") | (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="BB") |(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="BB") | (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AA"),]$Name)

		ParHaps <-cbind(ParHaps,Hap1)	 
			 
	}#end sib loop			 

write.table(ADIO,paste(outPath,"ADI_ADO_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)

	if (length(Sibs) ==1) {
		ParHaps <- as.data.frame(ParHaps)
		colnames(ParHaps) <- sib
	} else {
		colnames(ParHaps) <- Sibs
	}


	#------------------------------------------------------------
	#------------------------------------------------------------
	ParHaps[Gtypes$Chr=="X",] <- HapXMat
	library(signal)

	Haps <- vector("list",3)
	names(Haps) <- c("dataHapRaw","dataHap","Htype")

#	Gtypes <- Gtypes[Gtypes$Name %in% olaps1,]
	

	ParHapsAll <- cbind(Gtypes[,c("Name","Chr","Position")],ParHaps)
	dataHapRaw <- ParHapsAll[ParHapsAll$Name %in% c(olaps1,olaps0),]	## plot IBS1 & IBS0 region.
	rownames(dataHapRaw) <- NULL
	ParHapsAll <- dataHapRaw


	dataHap1 <- NULL
	for(chr in Chroms){
#		print(chr)
		IBSblocks <- list(IBS1block = IBS1[IBS1[,1] %in% chr,], IBS0block = IBS0[IBS0[,1] %in% chr,])
		for (b in 1:length(IBSblocks)) {
			IBSblock <- IBSblocks[[b]]
#			print(IBSblock)

			if ((is.data.frame(IBSblock) && nrow(IBSblock)==0)) {} else {
				for (bk in 1:nrow(IBSblock)) {
					dataHapRawSeg <- dataHapRaw[dataHapRaw$Chr %in% chr & dataHapRaw$Position >= IBSblock[bk,2] & dataHapRaw$Position <= IBSblock[bk,3],]
					if (nrow(dataHapRawSeg) > 0) {
						dataHap1Seg<-inthapnew1(dataHapRawSeg,Window,Int)
						dataHap1 <- rbind(dataHap1, dataHap1Seg)	
					}
				}
			}
		}
	}



	dataHap1 <- inthapnew1(dataHap1,Window,Int)


	if(Parent1==paste("Father_",Family,sep="")){ par1="_Pat"; par2="_Mat" }else{ par1="_Mat"; par2="_Pat" }

	Haps[["dataHapRaw"]] <- cbind(ParHapsAll,matrix(0,nrow(ParHapsAll),(ncol(ParHapsAll)-3)))
	colnames(Haps[["dataHapRaw"]])<- c("Name","Chr","Position",paste(colnames(dataHap1)[-c(1:3)],par1,sep=""),paste(colnames(dataHap1)[-c(1:3)],par2,sep=""))
	Haps[["dataHap"]] <- cbind(dataHap1,matrix(0,nrow(dataHap1),(ncol(dataHap1)-3)))
	colnames(Haps[["dataHap"]])<- colnames(Haps[["dataHapRaw"]])
	Haps[["Htype"]]  <- Htype
	colnames(Haps[["Htype"]]) <- colnames(Htype)	


	Haps

}#end function


