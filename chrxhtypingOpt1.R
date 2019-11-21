chrxhtypingOpt1 <- function(Gtypes,Parent1,ScSexes,Family) {

	GtypesChrX <- Gtypes[Gtypes$Chr=="X",]	
	MotherChrX <- GtypesChrX[,paste("Mother_",Family,sep="")]
	FatherChrX <- GtypesChrX[,paste("Father_",Family,sep="")]
	GrandfatherChrX <- GtypesChrX[,paste("Grandfather_",Family,sep="")]
	GrandmotherChrX <- GtypesChrX[,paste("Grandmother_",Family,sep="")]

	HapXMat <- matrix(0,nrow(GtypesChrX),ncol(as.matrix(GtypesChrX[,grep(SibPattern,colnames(GtypesChrX))])))
	colnames(HapXMat) <- paste(colnames(Gtypes)[grep(SibPattern,colnames(Gtypes))])

	if(Parent1==paste("Mother_",Family,sep="")){
		MotherChrX[(GtypesChrX$Grandfather=="BB" & GtypesChrX$Mother=="AB")] <- "BA"
		MotherChrX[(GtypesChrX$Grandfather=="NC" & GtypesChrX$Mother=="AB") |
			(GtypesChrX$Grandfather=="BB" & GtypesChrX$Grandmother=="BB" & GtypesChrX$Mother=="AB") | 
			(GtypesChrX$Grandfather=="AA" & GtypesChrX$Grandmother=="AA" & GtypesChrX$Mother=="AB") ] <-"NC"
		
		GtypesChrX[,Parent1] <- MotherChrX
	}	
   
	na.omit(ScSexes[ScSexes[,2]=="male",1])

	for(s in as.character(na.omit(ScSexes[ScSexes[,2]=="male",1]))){	
		print(paste(s,"is a male sample but has",sum(GtypesChrX[,s]=="AB"),"heterozygous SNP-calls; these are treated as NoCalls"))		
		## disease	   
		GtypesChrX[GtypesChrX[,s]=="AB",s]<-"NC"
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,s] == "BB") |
			(MotherChrX=="BA" & GtypesChrX[,s] == "AA"),s] <- 2
		## healthy
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,s] == "AA") |
			(MotherChrX=="BA" & GtypesChrX[,s] == "BB"),s] <- 1	
		HapXMat[GrandfatherChrX=="NC" | GrandmotherChrX =="NC",s] <- 0		# why its' needed?						
		print(s)							
	}#end s loop


	for(s in as.character(na.omit(ScSexes[ScSexes[,2]=="female",1]))){
		## disease
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,s] == "BB") |
			(MotherChrX=="BA" & GtypesChrX[,s] == "AA") |
			(MotherChrX=="AB" & FatherChrX == "AA" & GtypesChrX[,s] == "AB") |
			(MotherChrX=="BA" & FatherChrX == "BB" & GtypesChrX[,s] == "AB"),s] <- 2
		## healthy
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,s] == "AA") |
			(MotherChrX=="BA" & GtypesChrX[,s] == "BB") |
			(MotherChrX=="BA" & FatherChrX == "AA" & GtypesChrX[,s] == "AB")|
			(MotherChrX=="AB" & FatherChrX == "BB" & GtypesChrX[,s] == "AB"),s] <- 1
		HapXMat[GrandfatherChrX=="NC" | GrandmotherChrX=="NC" | FatherChrX=="NC",s] <- 0	# why its' needed							
		print(s)							
 	}#end s loop
	
	HapsX <- vector("list",2)
	names(HapsX) <- c("HapXMat","GtypesChrX")
	HapsX[["HapXMat"]] <- HapXMat
	HapsX[["GtypesChrX"]] <- GtypesChrX
	HapsX

}#end chrxhtypingOpt1 function

