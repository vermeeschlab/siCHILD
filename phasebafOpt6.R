phasebafOpt6 <- function(Gtypes,dataraw,Family,Parent1,gamma,outPath,Parents){ 

	print("Computing phased BAF...")

	source(paste(siCHILD_DIR,"callpcfBAF.R",sep=""))
	source(paste(siCHILD_DIR,"phstepparents.R",sep=""))


	if(sum(Gtypes[Gtypes[,"Chr"]=="X", paste("Father_",Family,sep="")]=="AB")>=1){
		warning(paste("On the ChrX  of the father",sum(Gtypes[Gtypes[,"Chr"]=="X", paste("Father_",Family,sep="")]=="AB"),"(",sum(Gtypes[Gtypes[,"Chr"]=="X", paste("Father_",Family,sep="")]=="AB")/sum(Gtypes$Chr=="X"), "%) of the SNP-calls are heterozygous!!","These are treated as NoCalls"))
		Gtypes[Gtypes[,"Chr"]=="X" & Gtypes[,paste("Father_",Family,sep="")]=="AB",paste("Father_",Family,sep="")] <- "NC"#There could not be htz SNPs on chromosome X of the father
	}

	plateau <-100
	PhBAF <- vector("list",8)
	names(PhBAF) <- c("P1","P2","M1","M2","P1Seg","P2Seg","M1Seg","M2Seg")

##### don't need to run the script of phstepparents.R anymore ######
#	Parents <- phstepparents(Gtypes,Parent1,Family,Condition)
#	fn <- paste(outPath,"Parents_",Family,".txt",sep="")
#print(fn)
#	write.table(Parents,fn,quote=F,sep="\t",col.names=T,row.names=F)


	Father <- Parents[,1]
	Mother <- Parents[,2]
write.table(Parents,paste(outPath,Family,"debugParents.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
	BAFs <- dataraw[,c(1:3,grep("Sibling",colnames(dataraw)),grep("^E",colnames(dataraw)),grep("Father.B",colnames(dataraw)),grep("Mother.B",colnames(dataraw)))]
	BAFs <- BAFs[,c(1:3,grep(".B.Allele.Freq",colnames(BAFs)))]
	#BAFs[(Father=="BA" | Mother== "BA") & BAFs$Name %in% olaps1,grep(".B.Allele.Freq",colnames(BAFs))] <- 1 - BAFs[(Father=="BA" | Mother== "BA") & BAFs$Name %in% olaps1, grep(".B.Allele.Freq",colnames(BAFs))]
	BAFs[Father=="BA" | Mother== "BA",grep(".B.Allele.Freq",colnames(BAFs))] <- 1 - BAFs[Father=="BA" | Mother== "BA",grep(".B.Allele.Freq",colnames(BAFs))]

	colnames(BAFs) <- gsub(".B.Allele.Freq","",colnames(BAFs))
	colnames(BAFs)[4:ncol(BAFs)] <- paste(colnames(BAFs)[4:ncol(BAFs)],"_",Family,sep="")

	ToRep <- BAFs 
	ToRep[,4:ncol(BAFs)] <- matrix(-1,ncol(BAFs)-3,nrow(BAFs)) 

	BAFs <- BAFs[Gtypes[,Seed]!="NC",]
	Father <- Father[Gtypes[,Seed]!="NC"]
	Mother <- Mother[Gtypes[,Seed]!="NC"]


	PhBAF[["P1"]] <- BAFs[((Father=="AB"& Mother=="AA") | (Father=="BA"& Mother=="BB")),]
	PhBAF[["P2"]] <- BAFs[((Father=="AB"& Mother=="BB") | (Father=="BA"& Mother=="AA")),]
	PhBAF[["M1"]] <- BAFs[((Mother=="AB"& Father=="AA") | (Mother=="BA"& Father=="BB")),]
	PhBAF[["M2"]] <- BAFs[((Mother=="AB"& Father=="BB") | (Mother=="BA"& Father=="AA")),]

	if(sum(grep("Father",Parent1))!=0){PhBAF[["M1"]]<-ToRep;PhBAF[["M2"]]<-ToRep;print("Parent1 is the Father")}else{PhBAF[["P1"]]<-ToRep;PhBAF[["P2"]]<-ToRep;print("Parent1 is the Mother")}

	gammaSC <- gamma
	gammaMC <- gamma


PhBAF[["P1Seg"]] <- PhBAF[["P1"]][0,]
PhBAF[["P2Seg"]] <- PhBAF[["P2"]][0,]
PhBAF[["M1Seg"]] <- PhBAF[["M1"]][0,]
PhBAF[["M2Seg"]] <- PhBAF[["M2"]][0,]

print("....P1Seg....")
print(nrow(PhBAF[["P1"]]))
	PhBAF[["P1Seg"]] <- callpcfBAFtry(PhBAF[["P1"]],gammaSC,gammaMC,plateau)

print("....P2Seg....")
print(nrow(PhBAF[["P2"]]))
	PhBAF[["P2Seg"]] <- callpcfBAFtry(PhBAF[["P2"]],gammaSC,gammaMC,plateau)

print("....M1Seg....")
print(nrow(PhBAF[["M1"]]))
	PhBAF[["M1Seg"]] <- callpcfBAFtry(PhBAF[["M1"]],gammaSC,gammaMC,plateau)

print("....M2Seg....")
print(nrow(PhBAF[["M2"]]))
	PhBAF[["M2Seg"]]<- callpcfBAFtry(PhBAF[["M2"]],gammaSC,gammaMC,plateau)

	write.table(PhBAF[["P1"]],paste(outPath,"P1_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(PhBAF[["P1Seg"]],paste(outPath,"P1Seg_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(PhBAF[["P2"]],paste(outPath,"P2_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(PhBAF[["P2Seg"]],paste(outPath,"P2Seg_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)

	write.table(PhBAF[["M1"]],paste(outPath,"M1_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(PhBAF[["M1Seg"]],paste(outPath,"M1Seg_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(PhBAF[["M2"]],paste(outPath,"M2_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	write.table(PhBAF[["M2Seg"]],paste(outPath,"M2Seg_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)

	PhBAF

}#end function
