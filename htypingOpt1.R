htypingOpt1 <- function(Gtypes,dataPo,ParScore,Parent1,SibPattern,Family,outPath){

	print("******************************************")
	print("****     Option1Htype analysis...     ****")
	
	ScSexes <- matrix(NA,length(ParScore),2)
	rownames(ScSexes)<-names(ParScore)

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
		warning(paste("On the chr. X  of the father",sum(Father=="AB"),"(",sum(Father=="AB")/sum(Gtypes$Chr=="X"), "%) of the SNP-calls are heterozygous!!"))
		print("These will be treated as NoCalls")
		Father[Father=="AB"] <- "NC"#There could not be htz SNPs on chromosome X of the father
		Gtypes[Gtypes[,"Chr"]=="X",paste("Father_",Family,sep="")] <- Father
	}




	#------------	Phasing chr X	-------------
	HapsX <- chrxhtypingOpt1(Gtypes,Parent1,ScSexes,Family)
	HapXMat <- HapsX[["HapXMat"]]
	GtypesChrX <- HapsX[["GtypesChrX"]]



	#------------	Phasing autosome  -------------
	rownames(Gtypes) <- NULL	
	PhasedPar <- Gtypes[,Parent1]


	Grandfather = paste("Grandfather_",Family,sep="")
	Grandmother = paste("Grandmother_",Family,sep="")
	PhasedPar[(Gtypes[,Grandfather]=="AB" & Gtypes[,Grandmother]=="AA" & Gtypes[,Parent1]=="AB")  |
		(Gtypes[,Grandfather]=="BB" & Gtypes[,Grandmother]=="AB" & Gtypes[,Parent1]=="AB") |
		(Gtypes[,Grandfather]=="BB" & Gtypes[,Grandmother]=="AA" & Gtypes[,Parent1]=="AB") ] <- "BA"

	PhasedPar[(Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="AA" & Gtypes[,Parent1]=="AB") |
		(Gtypes[,Grandfather]=="BB" & Gtypes[,Grandmother]=="BB" & Gtypes[,Parent1]=="AB")   |
		(Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="BB" & Gtypes[,Parent1]=="BB")   | 
		(Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="BB" & Gtypes[,Parent1]=="AA")   | 
		((Gtypes[,Grandfather]=="AA" | Gtypes[,Grandmother]=="AA") & Gtypes[,Parent1]=="BB") |
		((Gtypes[,Grandfather]=="AA" | Gtypes[,Grandmother]=="AB") & Gtypes[,Parent1]=="BB") |
                ((Gtypes[,Grandfather]=="BB" | Gtypes[,Grandmother]=="AA") & Gtypes[,Parent1]=="BB") |
                ((Gtypes[,Grandfather]=="BB" | Gtypes[,Grandmother]=="AA") & Gtypes[,Parent1]=="AA") |
                ((Gtypes[,Grandfather]=="BB" | Gtypes[,Grandmother]=="AB") & Gtypes[,Parent1]=="AA") |
		((Gtypes[,Grandfather]=="BB" | Gtypes[,Grandmother]=="BB") & Gtypes[,Parent1]=="AA") |
		((Gtypes[,Grandfather]=="NC" | Gtypes[,Grandmother]=="NC") & Gtypes[,Parent1]=="AB") ] <- "NC"

##### <------ I would like to simplified it into #######
#PhasedPar_RN  <- c(rownames(Gtypes[(Gtypes$Grandfather=="AB" & Gtypes$Grandmother=="AA" & Gtypes[,Parent1]=="AB") |(Gtypes$Grandfather=="BB" & Gtypes$Grandmother=="AB" & Gtypes[,Parent1]=="AB") |(Gtypes$Grandfather=="BB" & Gtypes$Grandmother=="AA" & Gtypes[,Parent1]=="AB") |(Gtypes$Grandmother=="AB" & Gtypes$Grandfather=="AA" & Gtypes[,Parent1]=="AB") |(Gtypes$Grandmother=="BB" & Gtypes$Grandfather=="AB" & Gtypes[,Parent1]=="AB") | (Gtypes$Grandmother=="BB" & Gtypes$Grandfather=="AA" & Gtypes[,Parent1]=="AB"),]))
#Gtypes[!(rownames(Gtypes) %in% PhasedPar_RN),][,Parent1] <- "NC"

	GtypesBA <- Gtypes[(Gtypes$Grandfather=="AB" & Gtypes$Grandmother=="AA" & Gtypes[,Parent1]=="AB") |(Gtypes$Grandfather=="BB" & Gtypes$Grandmother=="AB" & Gtypes[,Parent1]=="AB") |(Gtypes$Grandfather=="BB" & Gtypes$Grandmother=="AA" & Gtypes[,Parent1]=="AB"),]
	write.table(GtypesBA,paste(outPath,Family,"_GtypeBA.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
	GtypesAB <- Gtypes[(Gtypes$Grandmother=="AB" & Gtypes$Grandfather=="AA" & Gtypes[,Parent1]=="AB") |(Gtypes$Grandmother=="BB" & Gtypes$Grandfather=="AB" & Gtypes[,Parent1]=="AB") | (Gtypes$Grandmother=="BB" & Gtypes$Grandfather=="AA" & Gtypes[,Parent1]=="AB"),]
	write.table(GtypesAB,paste(outPath,Family,"_GtypeAB.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
##### for my own information -----> #########


	Gtypes_filter <- c()
	Gtypes_filter <- rbind(Gtypes_filter,(Gtypes[(Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="AA" & Gtypes[,Parent1]=="AB"),]))
	Gtypes_filter <- rbind(Gtypes_filter,(Gtypes[(Gtypes[,Grandfather]=="BB" & Gtypes[,Grandmother]=="BB" & Gtypes[,Parent1]=="AB"),]))

	Gtypes_filter <- rbind(Gtypes_filter,(Gtypes[(Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="BB" & Gtypes[,Parent1]=="AA"),]))
	Gtypes_filter <- rbind(Gtypes_filter,(Gtypes[(Gtypes[,Grandfather]=="BB" & Gtypes[,Grandmother]=="AB" & Gtypes[,Parent1]=="AA"),]))
	Gtypes_filter <- rbind(Gtypes_filter,(Gtypes[(Gtypes[,Grandfather]=="BB" & Gtypes[,Grandmother]=="AA" & Gtypes[,Parent1]=="AA"),]))

	Gtypes_filter <- rbind(Gtypes_filter,(Gtypes[(Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="BB" & Gtypes[,Parent1]=="BB"),]))
	Gtypes_filter <- rbind(Gtypes_filter,(Gtypes[(Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="AB" & Gtypes[,Parent1]=="BB"),]))
	Gtypes_filter <- rbind(Gtypes_filter,(Gtypes[(Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="AA" & Gtypes[,Parent1]=="BB"),]))

	if (nrow(Gtypes_filter) > 100) {print("warning: check parental genotypes!")}

	write.table(Gtypes_filter,paste(outPath,Family,"_GtypeFiltered.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


	Gtypes[,Parent1] <- PhasedPar
	Gtypes[Gtypes$Chr=="X",] <- GtypesChrX



	write.table(Gtypes,paste(outPath,Family,"_0.75.haplotype",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

	Sibs <- colnames(Gtypes)[grep(SibPattern,colnames(Gtypes))]
	ParHaps <- NULL
	for(sib in Sibs){
		Hap1 <- rep(0,nrow(Gtypes))
		print(sib)

		GtypeChild <- Gtypes[,sib]
	
		#-------------------------------------------------------------------------------------------------------
		## first allele is from GF, second allele is from GM; single cell inheritated 1st allele; assuming GM is sick, second allele is disease allele=> healthy => dark color
		Hap1[(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AB") |	
			(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="BB") |
			(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="BB") | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AB") | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AA") | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AA")] <- 1	


		## first allele is from GF, second allele is from GM; single cell inheritated 2nd allele; assuming GM is sick, second allele is disease allele => disease => light color
		Hap1[(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AA") |	# 
			(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AB") | # 
			(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AA") | # 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="BB") | # 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AB") | # 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="BB")] <- 2	
			 
		Hap1[Gtypes[,Grandfather]=="AB" & Gtypes[,Grandmother]=="AB"] <- 0	#? don't understand why it's needed
 			 
		#----------------------------------------------------------------------------------------------------------

		ParHaps <-cbind(ParHaps,Hap1)			 
			 
	}#end sib loop			 

	if (length(Sibs) ==1) {
		ParHaps <- as.data.frame(ParHaps)
		colnames(ParHaps) <- sib
	} else {
		colnames(ParHaps) <- Sibs
	}


	#------------------------------------------------------------
	Haps <- vector("list",3)
	names(Haps) <- c("dataHapRaw","dataHap","Parents")


	if(sum(grep("Father",Parent1))!=0 ){Father=Gtypes[,Parent1];Mother=Gtypes$Mother}else{Father=Gtypes$Father;Mother=Gtypes[,Parent1]}
	Parents <- cbind(Father,Mother)

	ParHaps[Gtypes$Chr=="X",] <- HapXMat

	ParHapsAll <- cbind(Gtypes[,c("Name","Chr","Position")],ParHaps)
	dataHapRaw <- ParHapsAll
	dataHap1 <- inthapnew1(dataHapRaw,Window,Int)

	if(Parent1==paste("Father_",Family,sep="")){ par1="_Pat"; par2="_Mat" }else{ par1="_Mat"; par2="_Pat" }

	Haps[["dataHapRaw"]] <- cbind(ParHapsAll,matrix(0,nrow(ParHapsAll),(ncol(ParHapsAll)-3)))
	colnames(Haps[["dataHapRaw"]])<- c("Name","Chr","Position",paste(colnames(dataHap1)[-c(1:3)],par1,sep=""),paste(colnames(dataHap1)[-c(1:3)],par2,sep=""))
	Haps[["dataHap"]] <- cbind(dataHap1,matrix(0,nrow(dataHap1),(ncol(dataHap1)-3)))
	colnames(Haps[["dataHap"]])<- colnames(Haps[["dataHapRaw"]])
	Haps[["Parents"]] <- Parents
	Haps

}#end function


