phasebafOpt3 <- function(Gtypes,Htype,dataraw,Family,Parent1,gamma,outPath,olaps0,olaps1,olaps2,IBS1,IBS0, NewSeed,HapCompare){


	print("Computing phased BAF...")

	#	source(paste(siCHILD_DIR,"callpcfBAF.R",sep=""))
	#	source(paste(siCHILD_DIR,"phrelatives.R",sep=""))



	#############################################################################
	#### to correct artifacts of haplotyping, minimize their effects for BAF ####
	#############################################################################
	Inds <- names(ParScore)
	Inds <- gsub(paste0("_",Family),"",Inds)
	colnames(HapCompare)[grep(paste(Inds,collapse="|"),colnames(HapCompare))] <- Inds

	for (ind in 1:length(Inds)) {
		HapCompare[Inds[ind]][(HapCompare[Inds[ind]]) > 3] <- 0
		HapCompare[Inds[ind]][(HapCompare[Inds[ind]]) < 3] <- 0
		HapCompare[Inds[ind]][(HapCompare[Inds[ind]]) == 3] <- 1
	}
	rownames(HapCompare) <- HapCompare$Name
	write.table(HapCompare,paste(outPath,Family,"_HapCompare_after.hap",sep=""),col.names=T,row.names=F,quote=F,sep="\t")




	if(sum(Gtypes[Gtypes[,"Chr"]=="X", paste("Father_",Family,sep="")]=="AB")>=1){
		warning(paste("On the ChrX  of the father",sum(Gtypes[Gtypes[,"Chr"]=="X", paste("Father_",Family,sep="")]=="AB"),"(",sum(Gtypes[Gtypes[,"Chr"]=="X", paste("Father_",Family,sep="")]=="AB")/sum(Gtypes$Chr=="X"), "%) of the SNP-calls are heterozygous!!","These are treated as NoCalls"))
		Gtypes[Gtypes[,"Chr"]=="X" & Gtypes[,paste("Father_",Family,sep="")]=="AB",paste("Father_",Family,sep="")] <- "NC"#There could not be htz SNPs on chromosome X of the father
	}

	plateau <-100
	PhBAF <- vector("list",9)
	names(PhBAF) <- c("P1","P2","M1","M2","P1Seg","P2Seg","M1Seg","M2Seg","pick")



	Parents <- Htype[Htype[,1] %in% c(olaps1,olaps0),c(grep("Father",colnames(Htype)),grep("Mother",colnames(Htype)))]

	#	if(sum(grep("Father",Parent1))!=0 ){Father=Gtypes[,Parent1];Mother=Gtypes$Mother}else{Father=Gtypes$Father;Mother=Gtypes[,Parent1]}
	#	Parents <- cbind(Father,Mother)
	#	write.table(Parents,paste(outPath,Family,Seed,"debugParents.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

	Father <- Parents[,1]
	Mother <- Parents[,2]




	#################################################################################################################
	#### SNPs need to be removed before segment BAF, e.g. Parent1&Seed == "AB"  e.g. Seed == "NC"  e.g. ADI, ADO ####
	#################################################################################################################
	p1col <- grep(strsplit(Parent1,"_")[[1]][1],colnames(dataraw))[1]
	seedcol <- grep(strsplit(Seed,"_")[[1]][1],colnames(dataraw))[1]
	bothHetRN <- dataraw[dataraw[,p1col] == dataraw[,seedcol] & dataraw[,p1col] %in% "AB",1]
	seedNCRN <- dataraw[dataraw[,seedcol] %in% "NC",1]
	ADIO <- read.table(paste(outPath,"ADI_ADO_",Family,".txt",sep=""),sep="\t",header=T,stringsAsFactors=F)[,1]
	ADIOolaps0 <- intersect(ADIO, olaps0)
	seedolap0RN <- read.table(paste(outPath,"seedolap0RN_",Family,".txt",sep=""),sep="\t",header=T,stringsAsFactors=F)[,1]	# Seed Het, Parent1 Hom



	######################
	#### prepare BAFs ####
	######################
	BAFs <- dataraw[,c(1:3,grep("Father.B",colnames(dataraw)),grep("Mother.B",colnames(dataraw)),grep("^E",colnames(dataraw)))]
	BAFs <- BAFs[,c(1:3,grep(".B.Allele.Freq",colnames(BAFs)))]
	BAFs <- (BAFs[BAFs[,1] %in% c(olaps1,olaps0),])
	write.table(BAFs,paste(outPath,"BAFs_before_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	BAFs_org <- BAFs


	#	#############################################################################
	#	#### swop BAF back base on HapComparison ####
	#	#############################################################################
	#	BAFs_embryo <- BAFs[,grep("^E",colnames(BAFs))]
	#	colnames(BAFs_embryo) <- gsub(".B.Allele.Freq",paste0("_",Family),colnames(BAFs_embryo))
	##	BAFs_embryo <- BAFs[,c(grep(paste(Inds,collapse="|"),colnames(BAFs)))]
	#	HapCompare_embryo <- HapCompare[,grep("^E",colnames(HapCompare))]
	##	HapCompare_embryo <- sapply( HapCompare_embryo, as.numeric )

	#	if ((ncol(BAFs_embryo) == ncol(HapCompare_embryo)) & (nrow(BAFs_embryo) == nrow(HapCompare_embryo))) {
	#		BAFs_embryo <- abs(HapCompare_embryo - BAFs_embryo)
	#		colnames(BAFs_embryo) <- gsub("$",".B.Allele.Freq",colnames(BAFs_embryo))
	#	} else {
	#		print("HapCompare_embryo and HapCompare_embryo have different ncol and nrow")
	#	}
	#	BAFs_tmp <- cbind(BAFs[,-grep("^E",colnames(BAFs))],BAFs_embryo)
	#	BAFs <- rbind(BAFs[BAFs$Name %in% olaps1,], BAFs_tmp[BAFs_tmp$Name %in% olaps0,])
	#	BAFs <- BAFs[order(BAFs[,2],BAFs[,3]),]
	#	rownames(BAFs) <- NULL




	###################################
	##### swop BAF  ######
	###################################
	BAFs[Father=="BA" | Mother== "BA",grep(".B.Allele.Freq",colnames(BAFs))] <- 1 - BAFs[Father=="BA" | Mother== "BA",grep(".B.Allele.Freq",colnames(BAFs))]

	#if (trisomyX == "Mother.B") {
	#	BAFstriX <- BAFs_org[BAFs_org$Chr=="X"]
	#	BAFstriX[Mother== "BA",grep(".B.Allele.Freq",colnames(BAFs))] <- 1 - BAFs[Mother== "BA",grep(".B.Allele.Freq",colnames(BAFs))] #ABB = 0.67
	#}

	colnames(BAFs) <- gsub(".B.Allele.Freq","",colnames(BAFs))
	colnames(BAFs)[4:ncol(BAFs)] <- paste(colnames(BAFs)[4:ncol(BAFs)],"_",Family,sep="")
	write.table(BAFs,paste(outPath,"BAFs_after_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)




	Gtypes <- Gtypes[Gtypes$Name %in% c(olaps1,olaps0),]
	`%ni%` = Negate(`%in%`)
	BAFs <- BAFs[BAFs$Name %ni% c(seedNCRN,bothHetRN,ADIOolaps0),]	# BAFs$Name %ni% ADIO[,1]
	BAFs_org <- BAFs_org[BAFs_org$Name %ni% c(seedNCRN,bothHetRN,ADIOolaps0),]
	colnames(BAFs_org) <- colnames(BAFs)

	Parents <- Htype[Htype[,1] %in% c(olaps1,olaps0) & Htype$Name %ni% c(seedNCRN,bothHetRN,ADIOolaps0) ,c(grep("Father",colnames(Htype)),grep("Mother",colnames(Htype)))]
	Father <- Parents[,1]
	Mother <- Parents[,2]




	ToRep <- BAFs
	ToRep[,4:ncol(BAFs)] <- matrix(-1,ncol(BAFs)-3,nrow(BAFs))


	PhBAF[["P1"]] <- BAFs[((Father=="AB"& Mother=="AA") | (Father=="BA"& Mother=="BB")),]
	PhBAF[["P2"]] <- BAFs[((Father=="AB"& Mother=="BB") | (Father=="BA"& Mother=="AA")),]

	PhBAF[["M1"]] <- BAFs[((Mother=="AB"& Father=="AA") | (Mother=="BA"& Father=="BB")),]
	PhBAF[["M2"]] <- BAFs[((Mother=="AB"& Father=="BB") | (Mother=="BA"& Father=="AA")),]


	P1o <- BAFs_org[((Father=="AB"& Mother=="AA") | (Father=="BA"& Mother=="BB")),]
	P2o <- BAFs_org[((Father=="AB"& Mother=="BB") | (Father=="BA"& Mother=="AA")),]

	M1o <- BAFs_org[((Mother=="AB"& Father=="AA") | (Mother=="BA"& Father=="BB")),]
	M2o <- BAFs_org[((Mother=="AB"& Father=="BB") | (Mother=="BA"& Father=="AA")),]

	write.table(P1o,paste(outPath,"P1Org_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(P2o,paste(outPath,"P2Org_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(M1o,paste(outPath,"M1Org_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(M2o,paste(outPath,"M2Org_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)



	if(sum(grep("Father",Parent1))!=0){PhBAF[["M1"]]<-ToRep;PhBAF[["M2"]]<-ToRep;print("Parent1 is the Father")}else{PhBAF[["P1"]]<-ToRep;PhBAF[["P2"]]<-ToRep;print("Parent1 is the Mother")}

	gammaSC <- gamma
	gammaMC <- gamma



	PhBAF[["P1Seg"]] <- PhBAF[["P1"]][0,]
	PhBAF[["P2Seg"]] <- PhBAF[["P2"]][0,]
	PhBAF[["M1Seg"]] <- PhBAF[["M1"]][0,]
	PhBAF[["M2Seg"]] <- PhBAF[["M2"]][0,]


	Seg0 <- PhBAF[["P1"]][0,][,c(1:3,grep("^E",colnames(BAFs)))]
	Seg1 <- Seg0


	print("....P1Seg / IBS01....")
	if (nrow(PhBAF[["P1"]][PhBAF[["P1"]]$Name %in% olaps0,][,c(1:3,grep("^E",colnames(BAFs)))]) > 0) {
		Seg0 <- callpcfBAFtry(PhBAF[["P1"]][PhBAF[["P1"]]$Name %in% olaps0,][,c(1:3,grep("^E",colnames(BAFs)))],gammaSC,gammaMC,plateau)
		write.table(Seg0,paste(outPath,"debugP1Seg0",Family,Seed,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	} else if (nrow(PhBAF[["P1"]][PhBAF[["P1"]]$Name %in% olaps1,][,c(1:3,grep("^E",colnames(BAFs)))]) > 0) {
		Seg1 <- callpcfBAFtry(PhBAF[["P1"]][PhBAF[["P1"]]$Name %in% olaps1,][,c(1:3,grep("^E",colnames(BAFs)))],gammaSC,gammaMC,plateau)
	}
	PhBAF[["P1Seg"]] <- rbind(Seg1,Seg0)


	print("....P2Seg / IBS01....")
	if (nrow(PhBAF[["P2"]][PhBAF[["P2"]]$Name %in% olaps0,][,c(1:3,grep("^E",colnames(BAFs)))]) > 0) {
		Seg0 <- callpcfBAFtry(PhBAF[["P2"]][PhBAF[["P2"]]$Name %in% olaps0,][,c(1:3,grep("^E",colnames(BAFs)))],gammaSC,gammaMC,plateau)
		write.table(Seg0,paste(outPath,"debugP2Seg0",Family,Seed,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	} else if (nrow(PhBAF[["P2"]][PhBAF[["P2"]]$Name %in% olaps1,][,c(1:3,grep("^E",colnames(BAFs)))]) > 0) {
		Seg1 <- callpcfBAFtry(PhBAF[["P2"]][PhBAF[["P2"]]$Name %in% olaps1,][,c(1:3,grep("^E",colnames(BAFs)))],gammaSC,gammaMC,plateau)
	}
	PhBAF[["P2Seg"]] <- rbind(Seg1,Seg0)


	print("....M1Seg / IBS01....")
	if (nrow(PhBAF[["M1"]][PhBAF[["M1"]]$Name %in% olaps0,][,c(1:3,grep("^E",colnames(BAFs)))]) > 0) {
		Seg0 <- callpcfBAFtry(PhBAF[["M1"]][PhBAF[["M1"]]$Name %in% olaps0,][,c(1:3,grep("^E",colnames(BAFs)))],gammaSC,gammaMC,plateau)
		write.table(Seg0,paste(outPath,"debugM1Seg0",Family,Seed,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	} else if (nrow(PhBAF[["M1"]][PhBAF[["M1"]]$Name %in% olaps1,][,c(1:3,grep("^E",colnames(BAFs)))]) > 0) {
		Seg1 <- callpcfBAFtry(PhBAF[["M1"]][PhBAF[["M1"]]$Name %in% olaps1,][,c(1:3,grep("^E",colnames(BAFs)))],gammaSC,gammaMC,plateau)
	}
	PhBAF[["M1Seg"]] <- rbind(Seg1,Seg0)


	print("....M2Seg / IBS01....")
	if (nrow(PhBAF[["M2"]][PhBAF[["M2"]]$Name %in% olaps0,][,c(1:3,grep("^E",colnames(BAFs)))]) > 0) {
		Seg0 <- callpcfBAFtry(PhBAF[["M2"]][PhBAF[["M2"]]$Name %in% olaps0,][,c(1:3,grep("^E",colnames(BAFs)))],gammaSC,gammaMC,plateau)
		write.table(Seg0,paste(outPath,"debugM2Seg0",Family,Seed,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	} else if (nrow(PhBAF[["M2"]][PhBAF[["M2"]]$Name %in% olaps1,][,c(1:3,grep("^E",colnames(BAFs)))]) > 0) {
		Seg1 <- callpcfBAFtry(PhBAF[["M2"]][PhBAF[["M2"]]$Name %in% olaps1,][,c(1:3,grep("^E",colnames(BAFs)))],gammaSC,gammaMC,plateau)
	}
	PhBAF[["M2Seg"]] <- rbind(Seg1,Seg0)




	#png("tmp.png")
	#par(mfrow=c(1,2))
	#plot(PhBAF[["M1"]][PhBAF[["M1"]]$Chr==4 & PhBAF[["M1"]]$Name %in% olaps0,6],ylim=c(0,1))
	#points(PhBAF[["M1Seg"]][PhBAF[["M1Seg"]]$Chr==4 & PhBAF[["M1Seg"]]$Name %in% olaps0,4],col="red")
	#plot(PhBAF[["M2"]][PhBAF[["M2"]]$Chr==4 & PhBAF[["M2"]]$Name %in% olaps0,6],ylim=c(0,1))
	#points(PhBAF[["M2Seg"]][PhBAF[["M2Seg"]]$Chr==4 & PhBAF[["M2Seg"]]$Name %in% olaps0,4],col="red")
	#dev.off()

	##################################################
	#####  output file write into subdirectory  ######
	##################################################
	if (NewSeed == "IBSseed") {
		subDir <- paste0(outPath,"/",NewSeed,"/")
		if (!file.exists(subDir)){
			dir.create(subDir)
		}
	} else {
		subDir <- outPath
	}

	outPath <- subDir

	write.table(PhBAF[["P1"]],paste(outPath,"P1_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(PhBAF[["P1Seg"]],paste(outPath,"P1Seg_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(PhBAF[["P2"]],paste(outPath,"P2_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(PhBAF[["P2Seg"]],paste(outPath,"P2Seg_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)

	write.table(PhBAF[["M1"]],paste(outPath,"M1_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(PhBAF[["M1Seg"]],paste(outPath,"M1Seg_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(PhBAF[["M2"]],paste(outPath,"M2_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(PhBAF[["M2Seg"]],paste(outPath,"M2Seg_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)



	if (length(grep("Father", Parent1)==1) == 1) {
		PhBAF[["pick"]] <- pickBestSingleCellbyChr(Gtypes[,c(1:3,grep("^E",colnames(Gtypes)))], PhBAF[["P1Seg"]],PhBAF[["P2Seg"]],Parent1,exc)
	} else if (length(grep("Mother", Parent1)==1) == 1) {
		PhBAF[["pick"]] <- pickBestSingleCellbyChr(Gtypes[,c(1:3,grep("^E",colnames(Gtypes)))], PhBAF[["M1Seg"]],PhBAF[["M2Seg"]],Parent1,exc)
	}

	PhBAF

}#end function
