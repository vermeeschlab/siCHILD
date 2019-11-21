#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## phasing with stepparent: sibling + one of the grandparents %%%%%%
############# This file is no longer in use ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



phstepparents <- function(Gtypes,Parent1,Family,Condition){



###################
#### autosome  ####
###################
rownames(Gtypes) <- NULL

PhasedPar <- Gtypes[,Parent1]

PhasedPar_RN <- c()

    if (Condition == "Sick_Sick") {	# grandparents-sick, sibling-sick
	## IBS1, sick siblings(seed2), first allele is diease allele
	PhasedPar_RN <- c(rownames(Gtypes[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB",]))
	PhasedPar[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB" ] <- "BA"	# common & disease allele B
	PhasedPar[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB" ] <- "AB"	# common & disease allele A

	## IBS0 or IBS1, assume sick sibling(seed2), first allele is disease allele; healthy sibling, second allele is disease allele
	PhasedPar_RN <- c(PhasedPar_RN,rownames(Gtypes[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB",])  )
	PhasedPar[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB" ] <- "AB"	# common & disease allele A
	PhasedPar[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB" ] <- "BA"	# common & disease allele B

#	## 
#	PhasedPar_RN <- c(rownames(Gtypes[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB",])  )
#	PhasedPar[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB" ] <- "AB"	# common & disease allele A
#	PhasedPar[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB" ] <- "BA"	# common & disease allele B

    } else if (Condition == "Healthy_Sick") {	# grandparents-healthy, sibling-sick

	## IBS0, healthy siblings(seed2), first allele is diease allele
	PhasedPar_RN <- c(rownames(Gtypes[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB",]))
	PhasedPar[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB" ] <- "AB"	# disease allele A
	PhasedPar[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB" ] <- "BA"	# disease allele B

	## IBS0 or IBS1, assume sick sibling(seed2), first allele is disease allele; healthy sibling, second allele is disease allele
	PhasedPar_RN <- c(PhasedPar_RN,rownames(Gtypes[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB",])  )
	PhasedPar[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB" ] <- "AB"	# disease allele B
	PhasedPar[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB" ] <- "BA"	# disease allele A

#	## 
#	PhasedPar_RN <- c(rownames(Gtypes[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB",])  )
#	PhasedPar[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB" ] <- "AB"	# disease allele B
#	PhasedPar[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB" ] <- "BA"	# disease allele A

    } else if (Condition == "Sick_Healthy") {	# grandparents-sick, sibling-healthy

	## IBS0, healthy siblings(seed2), first allele is diease allele
	PhasedPar_RN <- c(rownames(Gtypes[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB",]))
	PhasedPar[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB" ] <- "BA"	# disease allele A
	PhasedPar[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB" ] <- "BA"	# disease allele B

	## IBS0 or IBS1, assume sick sibling(seed2), first allele is disease allele; healthy sibling, second allele is disease allele
	PhasedPar_RN <- c(PhasedPar_RN,rownames(Gtypes[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB",])  )
	PhasedPar[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB" ] <- "AB"	# disease allele B
	PhasedPar[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB" ] <- "AB"	# disease allele A

#	## 
#	PhasedPar_RN <- c(rownames(Gtypes[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB",])  )
#	PhasedPar[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB" ] <- "AB"	# disease allele B
#	PhasedPar[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB" ] <- "BA"	# disease allele A

    } else if (Condition == "Healthy_Healthy") {	# grandparents-healthy, sibling-healthy

	## IBS0, healthy siblings(seed2), first allele is diease allele
	PhasedPar_RN <- c(rownames(Gtypes[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB",]))
	PhasedPar[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB" ] <- "AB"	# disease allele A
	PhasedPar[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="BB" & Gtypes[,Parent1]=="AB" ] <- "AB"	# disease allele B

	## IBS0 or IBS1, assume sick sibling(seed2), first allele is disease allele; healthy sibling, second allele is disease allele
	PhasedPar_RN <- c(PhasedPar_RN,rownames(Gtypes[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB",])  )
	PhasedPar[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB" ] <- "BA"	# disease allele B
	PhasedPar[Gtypes[,Seed1]=="AB" & Gtypes[,Seed2]=="AA" & Gtypes[,Parent1]=="AB" ] <- "BA"	# disease allele A

#	## 
#	PhasedPar_RN <- c(rownames(Gtypes[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB",]),rownames(Gtypes[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB",])  )
#	PhasedPar[Gtypes[,Seed1]=="AA" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB" ] <- "AB"	# disease allele B
#	PhasedPar[Gtypes[,Seed1]=="BB" & Gtypes[,Seed2]=="AB" & Gtypes[,Parent1]=="AB" ] <- "BA"	# disease allele A

    } else {
	print("Check")
    }



	Gtypes[,Parent1] <- PhasedPar

	Gtypes[!(rownames(Gtypes) %in% PhasedPar_RN),][,Parent1] <- "NC"	# change all the IBS == NA to "NC"





########################
#### Sex chromosome ####
########################

GtypesChrX <- Gtypes[Gtypes$Chr=="X",]
MotherChrX <- Gtypes[Gtypes$Chr=="X",paste("Mother_",Family,sep="")]	

if(Parent1==paste("Mother_",Family,sep="")){
MotherChrX_RN <- c()

    if (Condition == "Sick_Sick") {	# grandparents-sick, sibling-sick
	## IBS1, sick siblings(seed2), first allele is diease allele
	MotherChrX_RN <- c(rownames(GtypesChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]))
	MotherChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# common & disease allele B
	MotherChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# common & disease allele A

	## IBS0 or IBS1, assume sick sibling(seed2), first allele is disease allele; healthy sibling, second allele is disease allele
	MotherChrX_RN <- c(MotherChrX_RN,rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",])  )
	MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# common & disease allele A
	MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# common & disease allele B

#	## 
#	MotherChrX_RN <- c(rownames(GtypesChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB",])  )
#	MotherChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# common & disease allele A
#	MotherChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# common & disease allele B

    } else if (Condition == "Healthy_Sick") {	# grandparents-healthy, sibling-sick

	## IBS0, healthy siblings(seed2), first allele is diease allele
	MotherChrX_RN <- c(rownames(GtypesChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",]))
	MotherChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele A
	MotherChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele B

	## IBS0 or IBS1, assume sick sibling(seed2), first allele is disease allele; healthy sibling, second allele is disease allele
	MotherChrX_RN <- c(MotherChrX_RN,rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",])  )
	MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele B
	MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele A

#	## 
#	MotherChrX_RN <- c(rownames(GtypesChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB",])  )
#	MotherChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele B
#	MotherChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele A

    } else if (Condition == "Sick_Healthy") {	# grandparents-sick, sibling-healthy

	## IBS0, healthy siblings(seed2), first allele is diease allele
	MotherChrX_RN <- c(rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]))
	MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele A
	MotherChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele B

	## IBS0 or IBS1, assume sick sibling(seed2), first allele is disease allele; healthy sibling, second allele is disease allele
	MotherChrX_RN <- c(MotherChrX_RN,rownames(GtypesChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",])  )
	MotherChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele B
	MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele A

#	## 
#	MotherChrX_RN <- c(rownames(GtypesChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB",])  )
#	MotherChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele B
#	MotherChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele A

    } else if (Condition == "Healthy_Healthy") {	# grandparents-healthy, sibling-healthy

	## IBS0, healthy siblings(seed2), first allele is diease allele
	MotherChrX_RN <- c(rownames(GtypesChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB",]))
	MotherChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele A
	MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="BB" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele B

	## IBS0 or IBS1, assume sick sibling(seed2), first allele is disease allele; healthy sibling, second allele is disease allele
	MotherChrX_RN <- c(MotherChrX_RN,rownames(GtypesChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB",])  )
	MotherChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele B
	MotherChrX[GtypesChrX[,Seed1]=="AB" & GtypesChrX[,Seed2]=="AA" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele A

#	## 
#	MotherChrX_RN <- c(rownames(GtypesChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB",]),rownames(GtypesChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB",])  )
#	MotherChrX[GtypesChrX[,Seed1]=="AA" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB" ] <- "AB"	# disease allele B
#	MotherChrX[GtypesChrX[,Seed1]=="BB" & GtypesChrX[,Seed2]=="AB" & GtypesChrX[,Parent1]=="AB" ] <- "BA"	# disease allele A

    } else {
	print("Check")
    }


	GtypesChrX[,Parent1] <- MotherChrX
	GtypesChrX[!(rownames(GtypesChrX) %in% MotherChrX_RN),][,Parent1] <- "NC"	
}	
  

if(sum(grep("Father",Parent1))!=0 ){Father=Gtypes[,Parent1];Mother=Gtypes$Mother}else{Father=Gtypes$Father;Mother=Gtypes[,Parent1]}

Parents <- cbind(Father,Mother)

Parents

}#end function
