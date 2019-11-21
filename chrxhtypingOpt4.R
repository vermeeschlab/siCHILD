chrxhtypingOpt <- function(Gtypes,Parent1,ScSexes,Family) {

	GtypesChrX <- Gtypes[Gtypes$Chr=="X",]	
	MotherChrX <- GtypesChrX[,paste("Mother_",Family,sep="")]
	FatherChrX <- GtypesChrX[,paste("Father_",Family,sep="")]
	Seed1ChrX <- GtypesChrX[,Seed1]
	Seed2ChrX <- GtypesChrX[,Seed2]

	HapXMat <- matrix(0,nrow(GtypesChrX),ncol(as.matrix(GtypesChrX[,grep(SibPattern,colnames(GtypesChrX))])))
	colnames(HapXMat)<-paste(colnames(Gtypes)[grep(SibPattern,colnames(Gtypes))])


	if(Parent1==paste("Mother_",Family,sep="")){	
		print("same as autosome!")
		MotherChrX_RN <- rownames(GtypesChrX[(GtypesChrX[,Parent1]=="AA")|(GtypesChrX[,Parent1]=="BB"),])
		if (Condition == "Sick_Sick") {
			MotherChrX_RN <- c(MotherChrX_RN,rownames(GtypesChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]))
			MotherChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"
			MotherChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "AB"

			MotherChrX_RN <- c(MotherChrX_RN,rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",])  )
			MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "AB"
			MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"

		} else if (Condition == "Healthy_Sick") {
			MotherChrX_RN <- c(MotherChrX_RN,rownames(GtypesChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",]))
			MotherChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "AB"
			MotherChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"

			MotherChrX_RN <- c(MotherChrX_RN,rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",])  )
			MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "BA"
			MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "AB"

		} else if (Condition == "Sick_Healthy") {

		MotherChrX_RN <- c(MotherChrX_RN,rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]))
		MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele A
		MotherChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele B

		MotherChrX_RN <- c(MotherChrX_RN,rownames(GtypesChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",])  )
		MotherChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele B
		MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele A

		} else if (Condition == "Healthy_Healthy") {

		MotherChrX_RN <- c(MotherChrX_RN,rownames(GtypesChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",]))
		MotherChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele A
		MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele B

		MotherChrX_RN <- c(MotherChrX_RN,rownames(GtypesChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",])  )
		MotherChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele B
		MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele A
    		} 

		GtypesChrX[,Parent1] <- MotherChrX
		GtypesChrX[!(rownames(GtypesChrX) %in% MotherChrX_RN),][,Parent1] <- "NC"	# change all the IBS == NA to "NC"
	}   




	MotherChrX <- GtypesChrX[,paste("Mother_",Family,sep="")]
	FatherChrX <- GtypesChrX[,paste("Father_",Family,sep="")]
	SeedChrX <- GtypesChrX[,Seed]



	na.omit(ScSexes[ScSexes[,2]=="male",1])

	for(s in as.character(na.omit(ScSexes[ScSexes[,2]=="male",1]))){	
		print(paste(s,"is a male sample but has",sum(GtypesChrX[,s]=="AB"),"heterozygous SNP-calls; these are treated as NoCalls"))			   
		GtypesChrX[GtypesChrX[,s]=="AB",s]<-"NC"
		## healthy
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,s] == "BB") |
			(MotherChrX=="BA" & GtypesChrX[,s] == "AA"),s] <- 2
		## sick
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,s] == "AA") |
			(MotherChrX=="BA" & GtypesChrX[,s] == "BB"),s] <- 1
		#HapXMat[SeedChrX=="NC",s] <- 0						
		print(s)							
	}#end s loop


	for(s in as.character(na.omit(ScSexes[ScSexes[,2]=="female",1]))){
		## healthy
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,s] == "BB") |
			(MotherChrX=="BA" & GtypesChrX[,s] == "AA") |
			(MotherChrX=="AB" & FatherChrX == "AA" & GtypesChrX[,s] == "AB") |
			(MotherChrX=="BA" & FatherChrX == "BB" & GtypesChrX[,s] == "AB"),s] <- 2
		## sick
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,s] == "AA") |
			(MotherChrX=="BA" & GtypesChrX[,s] == "BB") |
			(MotherChrX=="BA" & FatherChrX == "AA" & GtypesChrX[,s] == "AB")|
			(MotherChrX=="AB" & FatherChrX == "BB" & GtypesChrX[,s] == "AB"),s] <- 1
		#HapXMat[SeedChrX=="NC" | FatherChrX=="NC",s] <- 0						
		print(s)							
 	}#end s loop
	
#	HapXMat
	HapsX <- vector("list",2)
	names(HapsX) <- c("HapXMat","GtypesChrX")
	HapsX[["HapXMat"]] <- HapXMat
	HapsX[["GtypesChrX"]] <- GtypesChrX
	HapsX

}#end chrxhtypingOpt1 function

