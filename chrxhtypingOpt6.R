chrxhtypingOpt6 <- function(Gtypes,Parent1,ScSexes,Family) {

	GtypesChrX <- Gtypes[Gtypes$Chr=="X",]	
	MotherChrX <- GtypesChrX[,paste("Mother_",Family,sep="")]
	FatherChrX <- GtypesChrX[,paste("Father_",Family,sep="")]

	SeedChrX <- GtypesChrX[,Seed]


	HapXMat <- matrix(0,nrow(GtypesChrX),ncol(as.matrix(GtypesChrX[,c(grep("Sibling",colnames(Gtypes)),grep(SibPattern,colnames(Gtypes)))])))
	colnames(HapXMat)<-paste(colnames(Gtypes)[c(grep("Sibling",colnames(Gtypes)),grep(SibPattern,colnames(Gtypes)))])


	##### only consider parental mother case #####

	if(Parent1==paste("Mother_",Family,sep="")){	
		print("same as autosome!")
		MotherChrX_RN <- c(rownames(GtypesChrX[GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1]=="AB",]))
		MotherChrX[GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# B = disease allele
		MotherChrX[GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# A = disease allele

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

