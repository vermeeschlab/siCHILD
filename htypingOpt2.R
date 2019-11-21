htypingOpt2 <- function(Gtypes,Window,Int,dataPo,ParScore,SibPattern,Family,outPath){

	Haps <- vector("list",2)
	names(Haps) <- c("dataHapRaw","dataHap")

	SeedCol <- "Sibling"

	Children <- unique(c(grep(SibPattern,colnames(Gtypes)),grep(SeedCol,colnames(Gtypes))))

	
	print(paste0("Children is ", Children))

	Sibs <- as.data.frame(Gtypes[,Children])

	if(length(Children)==1){Sibs<-as.matrix(Sibs);colnames(Sibs)<-colnames(Gtypes)[grep(SibPattern,colnames(Gtypes))]}
	colnames(Sibs) <- colnames(Gtypes)[Children]

	HapsAut <- htypingOpt2Aut(Gtypes,Gtypes[,paste("Father_",Family,sep="")],Gtypes[,paste("Mother_",Family,sep="")],Gtypes[,grep(SeedCol,colnames(Gtypes))],Sibs)


	print("Autosomes are haplotyped")
	HapsChrX <- chrxhtypingOpt2(Gtypes,ParScore,SibPattern)
	print("chrX is haplotyped")

	Haps[["dataHapRaw"]] <- rbind(HapsAut,HapsChrX)
	Haps[["dataHap"]]<-inthapnew1(Haps[["dataHapRaw"]],Window,Int)
	Haps

}#end function

