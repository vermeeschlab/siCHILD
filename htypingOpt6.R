htypingOpt6 <- function(Gtypes,dataPo,ParScore,Parent1,SibPattern,Family,outPath){

	print("******************************************")
	print("****     half-sibling analysis...     ****")
	
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
		print(paste("On the chr. X  of the father",sum(Father=="AB"),"(",sum(Father=="AB")/sum(Gtypes$Chr=="X"), "%) of the SNP-calls are heterozygous!!"))
		print("These will be treated as NoCalls")
		Father[Father=="AB"] <- "NC"#There could not be htz SNPs on chromosome X of the father
		Gtypes[Gtypes[,"Chr"]=="X",paste("Father_",Family,sep="")] <- Father
	}




	#------------	Phasing chr X	-------------
	HapsX <- chrxhtypingOpt6(Gtypes,Parent1,ScSexes,Family)
	HapXMat <- HapsX[["HapXMat"]]
	GtypesChrX <- HapsX[["GtypesChrX"]]


	rownames(Gtypes) <- NULL
	PhasedPar <- Gtypes[,Parent1]


print(paste0("~~~~~~",Seed))
	PhasedPar_RN <- c(rownames(Gtypes[Gtypes[,Seed]=="BB" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed]=="AA" & Gtypes[,Parent1]=="AB",]))

	PhasedPar[Gtypes[,Seed]=="BB" & Gtypes[,Parent1]=="AB" ] <- "BA"	# B = disease allele
	PhasedPar[Gtypes[,Seed]=="AA" & Gtypes[,Parent1]=="AB" ] <- "AB"	# A = disease allele
	Gtypes[,Parent1] <- PhasedPar
	Gtypes[!(rownames(Gtypes) %in% PhasedPar_RN),][,Parent1] <- "NC"	# change all the IBS == NA to "NC"

	Gtypes[Gtypes$Chr=="X",] <- GtypesChrX

	write.table(Gtypes,paste(outPath,Family,"_0.75.haplotype",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


	Sibs <- colnames(Gtypes)[c(grep("Sibling",colnames(Gtypes)),grep(SibPattern,colnames(Gtypes)))]
	ParHaps <- NULL
	for(sib in Sibs){
		Hap1 <- rep(0,nrow(Gtypes))
		#print(sib)

		GtypeChild <- Gtypes[,sib]
	
		#-------------------------------------------------------------------------------------------------------
		## first allele is disease allele, single cell inherited second allele => healthy => light color
		Hap1[(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AA") |	
			(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AB") | 
			(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AA") | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="BB") | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AB") |  
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="BB")] <- 2	# 2 = light color


		## first allele is disease allele, single cell inherited first allele => sick => dark color
		Hap1[(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AB") |	
			(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="BB") | 
			(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="BB") | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AB") | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AA") | 
			(Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AA")] <- 1	# 1 = dark color
		#----------------------------------------------------------------------------------------------------------0

		

		if(sib==Sibs[1]){ParHaps <- Hap1;ParHaps <- as.data.frame(ParHaps)}else{ParHaps <-cbind(ParHaps,Hap1)}			 
			 
	}#end sib loop			 

	

	colnames(ParHaps) <- Sibs


	#------------------------------------------------------------

	ParHaps[Gtypes$Chr=="X",] <- HapXMat
	library(signal)

	Haps <- vector("list",3)
	names(Haps) <- c("dataHapRaw","dataHap","Parents")

#	Gtypes <- Gtypes[Gtypes$Name %in% olaps1,]
	

	if(sum(grep("Father",Parent1))!=0 ){Father=Gtypes[,Parent1];Mother=Gtypes$Mother}else{Father=Gtypes$Father;Mother=Gtypes[,Parent1]}
	Parents <- cbind(Father,Mother)

	ParHapsAll <- cbind(Gtypes[,c("Name","Chr","Position")],ParHaps)
	dataHapRaw <- ParHapsAll
	dataHap1<-inthapnew1(dataHapRaw,Window,Int)

	if(Parent1==paste("Father_",Family,sep="")){ par1="_Pat"; par2="_Mat" }else{ par1="_Mat"; par2="_Pat" }

	Haps[["dataHapRaw"]] <- cbind(ParHapsAll,matrix(0,nrow(ParHapsAll),(ncol(ParHapsAll)-3)))
	colnames(Haps[["dataHapRaw"]])<- c("Name","Chr","Position",paste(colnames(dataHap1)[-c(1:3)],par1,sep=""),paste(colnames(dataHap1)[-c(1:3)],par2,sep=""))
	Haps[["dataHap"]] <- cbind(dataHap1,matrix(0,nrow(dataHap1),(ncol(dataHap1)-3)))
	colnames(Haps[["dataHap"]])<- colnames(Haps[["dataHapRaw"]])
	Haps[["Parents"]] <- Parents

	Haps

}#end function


