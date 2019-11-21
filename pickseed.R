pickseed <- function(pick_seed,ParScore,Gtypes, dataraw, chr){

	if (pick_seed == "NA")  {

		Gtypes[,ncol(Gtypes)+1] <- NA
		colnames(Gtypes)[ncol(Gtypes)]   <- "Sibling"
		dataraw[,ncol(dataraw)+1] <- NA
		colnames(dataraw)[ncol(dataraw)] <- "Sibling.GType"	
		dataraw[,ncol(dataraw)+1] <- NA
		colnames(dataraw)[ncol(dataraw)] <- "Sibling.B.Allele.Freq"	
		dataraw[,ncol(dataraw)+1] <- NA
		colnames(dataraw)[ncol(dataraw)] <- "Sibling.Log.R.Ratio"	

		empty_chr <- NULL
		for(chr in 1:22){
			if (!is.na(PhBAF[["pick"]][chr,2])) {
				Gtypes[Gtypes[,2] == chr,ncol(Gtypes)] <- Gtypes[Gtypes[,2] == chr,PhBAF[["pick"]][chr,2]]
				NewSiblingCol <- paste(strsplit(PhBAF[["pick"]][chr,2],"_")[[1]][1], strsplit(PhBAF[["pick"]][chr,2],"_")[[1]][2], sep="_")	
				dataraw[dataraw[,2] == chr,(ncol(dataraw)-2):ncol(dataraw)] <- dataraw[dataraw[,2] == chr,grep(NewSiblingCol,colnames(dataraw))]
			} else {
				empty_chr <- chr
			}
		}

		#############################################################################################
		##### ChrX has problem, choose the best single cell by call rate, then choose its ChrX ###### 
		#############################################################################################
		#unlist(lapply(ParScore, function(x) max(x["X",5])))
		NewSeed_select <- cbind(t(as.data.frame(lapply(ParScore, function(x) mean(x[1:22,5])))), t(as.data.frame(lapply(ParScore, function(x) sd(x[1:22,5])))))
		colnames(NewSeed_select) <- c("mean","sd")
		NewSeedX <- rownames(NewSeed_select[order(-NewSeed_select[,1],-NewSeed_select[,2]),])[1]
	
		NewSiblingCol <- paste(strsplit(NewSeedX,"_")[[1]][1], strsplit(NewSeedX,"_")[[1]][2], sep="_")

		Gtypes[Gtypes[,2] %in% c("X","Y",empty_chr),ncol(Gtypes)] <- Gtypes[Gtypes[,2] %in% c("X","Y",empty_chr),NewSeedX]
		dataraw[dataraw[,2] %in% c("X","Y",empty_chr),(ncol(dataraw)-2):ncol(dataraw)] <- dataraw[dataraw[,2] %in% c("X","Y",empty_chr),grep(NewSiblingCol,colnames(dataraw))]


	} else if (pick_seed == "new") {	# use by sibling phasing option

		Gtypes[,ncol(Gtypes)+1] <- NA
		colnames(Gtypes)[ncol(Gtypes)]   <- "Sibling"
		dataraw[,ncol(dataraw)+1] <- NA
		colnames(dataraw)[ncol(dataraw)] <- "Sibling.GType"	
		dataraw[,ncol(dataraw)+1] <- NA
		colnames(dataraw)[ncol(dataraw)] <- "Sibling.B.Allele.Freq"	
		dataraw[,ncol(dataraw)+1] <- NA
		colnames(dataraw)[ncol(dataraw)] <- "Sibling.Log.R.Ratio"	

		empty_chr <- 1:22

		NewSeed_select <- cbind(t(as.data.frame(lapply(ParScore, function(x) mean(x[1:22,5])))), t(as.data.frame(lapply(ParScore, function(x) sd(x[1:22,5])))))
		colnames(NewSeed_select) <- c("mean","sd")
		NewSeed_select <- NewSeed_select[ grep("^E",rownames(NewSeed_select)),]

		NewSeedX <- rownames(NewSeed_select[order(-NewSeed_select[,1],-NewSeed_select[,2]),])[1]

		NewSiblingCol <- paste(strsplit(NewSeedX,"_")[[1]][1], strsplit(NewSeedX,"_")[[1]][2], sep="_")

		Gtypes[Gtypes[,2] %in% c("X","Y",empty_chr),ncol(Gtypes)] <- Gtypes[Gtypes[,2] %in% c("X","Y",empty_chr),NewSeedX]
		dataraw[dataraw[,2] %in% c("X","Y",empty_chr),(ncol(dataraw)-2):ncol(dataraw)] <- dataraw[dataraw[,2] %in% c("X","Y",empty_chr),grep(NewSiblingCol,colnames(dataraw))]



	} else {
		Gtypes[,ncol(Gtypes)+1] <- Gtypes[,grep(pick_seed,colnames(Gtypes))]
		colnames(Gtypes)[ncol(Gtypes)] <- "Sibling"
		nc <- ncol(dataraw)
		dataraw[,(nc+1):(nc+3)] <- dataraw[,grep(pick_seed,colnames(dataraw))]
		colnames(dataraw)[(nc+1):(nc+3)] <- c("Sibling.GType","Sibling.B.Allele.Freq","Sibling.Log.R.Ratio")
	}

	datasGtypes[["Gtypes"]] <- Gtypes
	datasGtypes[["dataraw"]] <- dataraw
	datasGtypes

}
