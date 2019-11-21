#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## phasing with uncle/aunt
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phrelatives <- function(Gtypes,Parent1,Family,olaps0,olaps1,olaps2){


###################
#### autosome  ####
###################

#autosomeChr <- Chroms[Chroms!="X"]

	PhasedPar <- Gtypes[,Parent1]
	if (!is.integer0(grep(Seed,colnames(Gtypes)))) {	
		PhasedPar_RN <- c(rownames(Gtypes[Gtypes[,Seed]=="BB" & Gtypes[,Parent1]=="AB" & Gtypes$Name %in% c(olaps1),]),rownames(Gtypes[Gtypes[,Seed]=="AA" & Gtypes[,Parent1]=="AB" & Gtypes$Name %in% c(olaps1),]), rownames(Gtypes[((Gtypes[,Seed]=="BB" & Gtypes[,Parent1] == "AB") | (Gtypes[,Seed]=="AB" & Gtypes[,Parent1] == "BB")) & Gtypes$Name %in% olaps0,]),rownames(Gtypes[((Gtypes[,Seed]=="AA" & Gtypes[,Parent1] == "AB") | (Gtypes[,Seed]=="AB" & Gtypes[,Parent1] == "AA")) & Gtypes$Name %in% olaps0,]))
		PhasedPar[Gtypes[,Seed]=="BB" & Gtypes[,Parent1]=="AB" & Gtypes$Name %in% olaps1] <- "BA"
		PhasedPar[Gtypes[,Seed]=="AA" & Gtypes[,Parent1]=="AB" & Gtypes$Name %in% olaps1] <- "AB"	


#		### this doesn't work
#		PhasedPar[Gtypes[,Seed]=="BB" & Gtypes[,Parent1] %in% c("AB","AA") & Gtypes$Name %in% olaps0] <- "AB"
#		PhasedPar[Gtypes[,Seed]=="AA" & Gtypes[,Parent1] %in% c("AB","BB") & Gtypes$Name %in% olaps0] <- "BA"


		PhasedPar[((Gtypes[,Seed]=="BB" & Gtypes[,Parent1] == "AB") | (Gtypes[,Seed]=="AB" & Gtypes[,Parent1] == "BB")) & Gtypes$Name %in% olaps0] <- "AB"
		PhasedPar[((Gtypes[,Seed]=="AA" & Gtypes[,Parent1] == "AB") | (Gtypes[,Seed]=="AB" & Gtypes[,Parent1] == "AA")) & Gtypes$Name %in% olaps0] <- "BA"



	}

	Gtypes[,Parent1] <- PhasedPar


########################
#### Sex chromosome ####
########################
	GtypesChrX <- Gtypes[Gtypes$Chr=="X",]	
	MotherChrX <- GtypesChrX[,paste("Mother_",Family,sep="")]
	FatherChrX <- GtypesChrX[,paste("Father_",Family,sep="")]
	SeedChrX <- GtypesChrX[,paste(Seed,sep="")]


MotherChrX_RN <- c()
FatherChrX_RN <- c()

if((Parent1==paste("Mother_",Family,sep="")) & (Seed == paste("Aunt_",Family,sep=""))) {
	print("same as autosome!")

	MotherChrX_RN  <- c(rownames(GtypesChrX[(GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1] %in% c("AB","BA")) & GtypesChrX$Name %in% c(olaps1,olaps0),]), rownames(GtypesChrX[(GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1] %in% c("AB","BA")) & GtypesChrX$Name %in% c(olaps1,olaps0),]) )

	MotherChrX[GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1] %in% c("AB","BA") & GtypesChrX$Name %in% c(olaps1)] <- "BA" 
	MotherChrX[GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1] %in% c("AB","BA") & GtypesChrX$Name %in% c(olaps1)] <- "AB"

	MotherChrX[GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1] %in% c("AB","BA") & GtypesChrX$Name %in% c(olaps0)] <- "AB" 
	MotherChrX[GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1] %in% c("AB","BA") & GtypesChrX$Name %in% c(olaps0)] <- "BA"


} else if((Parent1==paste("Mother_",Family,sep="")) & (Seed == paste("Uncle_",Family,sep=""))) {
	MotherChrX_RN  <- c(rownames(GtypesChrX[(GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1] %in% c("AB","BA")),]), rownames(GtypesChrX[(GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1] %in% c("AB","BA")),]) )

	MotherChrX[GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1] %in% c("AB","BA") ] <- "BA"	
	MotherChrX[GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1] %in% c("AB","BA") ] <- "AB"	


} else if((Parent1==paste("Father_",Family,sep="")) & (Seed == paste("Uncle_",Family,sep=""))) {
	FatherChrX_RN  <- c(rownames(GtypesChrX[(GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1] %in% c("BB")),]), rownames(GtypesChrX[(GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1] %in% c("AA")),]) )

#	FatherChrX[GtypesChrX[,Seed]=="BB" & GtypesChrX[,Parent1] %in% c("BB") ] <- "BB"	
#	FatherChrX[GtypesChrX[,Seed]=="AA" & GtypesChrX[,Parent1] %in% c("AA") ] <- "AA"
	
} else if((Parent1==paste("Father_",Family,sep="")) & (Seed == paste("Aunt_",Family,sep=""))) {

	FatherChrX_RN  <- c(rownames(GtypesChrX[(GtypesChrX[,Parent1]=="BB" & GtypesChrX[,Seed] %in% c("AB","BA")),]), rownames(GtypesChrX[(GtypesChrX[,Parent1]=="AA" & GtypesChrX[,Seed] %in% c("AB","BA")),]) )

#	FatherChrX[GtypesChrX[,Parent1]=="BB" & GtypesChrX[,Seed] %in% c("AB","BA") ] <- "BB"	
#	FatherChrX[GtypesChrX[,Parent1]=="AA" & GtypesChrX[,Seed] %in% c("AB","BA") ] <- "AA"	
}	
	   

	if (Parent1==paste("Mother_",Family,sep="")) {
		Gtypes[Gtypes$Chr=="X",Parent1] <- MotherChrX
		GtypesChrX[!(rownames(GtypesChrX) %in% MotherChrX_RN),][,Parent1] <- "NC"	# change all the IBS == NA to "NC"
	} else if (Parent1==paste("Father_",Family,sep="")) {
		Gtypes[Gtypes$Chr=="X",Parent1] <- FatherChrX
		GtypesChrX[!(rownames(GtypesChrX) %in% FatherChrX_RN),][,Parent1] <- "NC"	# change all the IBS == NA to "NC"
	}


if(sum(grep("Father",Parent1))!=0 ){Father=Gtypes[,Parent1];Mother=Gtypes$Mother}else{Father=Gtypes$Father;Mother=Gtypes[,Parent1]}

Parents <- cbind(Father,Mother)

Parents

}#end function
