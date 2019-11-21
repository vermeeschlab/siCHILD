
produce_ped <- function(data,Family,outPath){
	GType_name <- colnames(data)[grep("GType",colnames(data))]
GType <- vector("list", length(GType_name))
GType_all <- NULL
for (ind in 1:length(GType_name)) {
	genotype <- (unlist(strsplit(as.character(a[,GType_name[ind]]),"")))
	genotype[genotype=="N"] <- 0
	genotype[genotype=="C"] <- 0
	genotype[genotype=="A"] <- 1
	genotype[genotype=="B"] <- 2
	GType_name[ind] <- gsub(".GType","",GType_name[ind])
	SampleID <- unlist((strsplit(GType_name[ind],"\\.")))[1]
	FamilyID <- Family
	if (GType_name[ind] %in%  c("Father","Mother","Aunt","Uncle")) { PaternalID <- "Grandfather"; MaternalID <- "Grandmother"} else if (length(GType_name[ind][grep("^E",GType_name[ind])]) > 0) { PaternalID <- "Father"; MaternalID <- "Mother"} else {PaternalID <- 0; MaternalID <- 0;}
	if (GType_name[ind] == "ather") { Sex <- 1 } else if (GType_name[ind] == "other") { Sex <- 2 } else { Sex <- 0 }
	Affection <- 0
	GType[[ind]] <-	c(FamilyID,SampleID,PaternalID,MaternalID,Sex,Affection,genotype)
	GType_all <- rbind(GType_all,GType[[ind]])
}
write.table(GType_all,paste(outPath,Family,".ped",sep=""),sep=" ", quote=F, col.names=F, row.names=F)

}



