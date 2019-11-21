chrxhtypingOpt3 <- function(Gtypes,Parent1,ScSexes,Family,olaps0,olaps1,olaps2) {

	GtypesChrX <- Gtypes[Gtypes$Chr=="X",]	
	MotherChrX <- GtypesChrX[,paste("Mother_",Family,sep="")]
	FatherChrX <- GtypesChrX[,paste("Father_",Family,sep="")]
	SeedChrX <- GtypesChrX[,paste(Seed,sep="")]
write.table(GtypesChrX,paste(outPath,"GtypesChrX_before",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	HapXMat <- matrix(0,nrow(GtypesChrX),ncol(as.matrix(GtypesChrX[,grep(SibPattern,colnames(GtypesChrX))])))
	colnames(HapXMat)<-paste(colnames(Gtypes)[grep(SibPattern,colnames(Gtypes))])


MotherChrX_RN <- c()
FatherChrX_RN <- c()
if (Parent1==paste("Mother_",Family,sep="")) {	# Aunt has "AB", but it's useless as AB & AB is not informative
#if((Parent1==paste("Mother_",Family,sep="")) & (Seed == paste("Aunt_",Family,sep=""))) {
#	print("same as autosome!")

	MotherChrX_RN <- c(rownames(GtypesChrX[GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1]=="AB" & GtypesChrX$Name %in% c(olaps1),]),rownames(GtypesChrX[GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1]=="AB" & GtypesChrX$Name %in% c(olaps1),]), rownames(GtypesChrX[((GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1] == "AB")) & GtypesChrX$Name %in% olaps0,]),rownames(GtypesChrX[((GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1] == "AB")) & GtypesChrX$Name %in% olaps0,]))
	MotherChrX[GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1]=="AB" & GtypesChrX$Name %in% olaps1] <- "BA"
	MotherChrX[GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1]=="AB" & GtypesChrX$Name %in% olaps1] <- "AB"	

	MotherChrX[((GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1] == "AB")) & GtypesChrX$Name %in% olaps0] <- "AB"
	MotherChrX[((GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1] == "AB")) & GtypesChrX$Name %in% olaps0] <- "BA"


#} else if((Parent1==paste("Mother_",Family,sep="")) & (Seed == paste("Uncle_",Family,sep=""))) {	# only has IBS0 & IBS1 => PGD016
#	MotherChrX_RN  <- c(rownames(GtypesChrX[(GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1] == "AB"),]), rownames(GtypesChrX[(GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1] == "AB"),]) )

#	MotherChrX[GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1] == "AB" ] <- "BA"	
#	MotherChrX[GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1] == "AB" ] <- "AB"	


}	


  

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
	
	HapXMat

}#end chrxhtypingOpt1 function

