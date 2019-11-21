DIR <- "/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/"
Family <- "PGD117"
flank = 150

is.even <- function(x) x %% 2 == 0
is.odd <- function(x) x %% 2 != 0
Chroms <- c(1:22,"X")

names <- read.table(paste0(DIR,"PGD117/names.adj"),header=T, stringsAsFactors=F)
double <- read.table("/uz/data/avalok/symbiosys/raw/HiSeqComputed/new/research/kul/joris_vermeesch/gcpu/samples/PGD117_C1/pgd.sichild/current/result/PGD117_C1_Itp2.hap",header=T, stringsAsFactors=F)
GM <- read.table(paste0(DIR,"PGD117_out_SingleGP_Grandmother_v2/PGD117_Itp2.hap"),header=T, stringsAsFactors=F)
GF <- read.table(paste0(DIR,"PGD117_out_SingleGP_Grandfather_v2/PGD117_Itp2.hap"),header=T, stringsAsFactors=F)
dGF <- merge(names,double,by=c("Name","Chr","Position"),all=T)
dGF <- merge(dGF,GF,by=c("Name","Chr","Position"),all=T)
dGF <- merge(dGF,GM,by=c("Name","Chr","Position"),all=T)
dGF <- dGF[,c(1:5,8,9,12,13)]
dGF <- (dGF[!dGF[,2] %in% c(0,"XY","Y"),])
dGF <- (dGF[complete.cases(dGF),])
#[1] 293068
#dGF[rowSums(is.na(dGF[,5:6]))==0,]
#[1] 293068




dGF[dGF[,4] == 2,4] <- 11
dGF[dGF[,4] == 1,4] <- 22
dGF[dGF[,4] == 22,4] <- 2
dGF[dGF[,4] == 11,4] <- 1
dGF[dGF[,5] == 2,5] <- 11
dGF[dGF[,5] == 1,5] <- 22
dGF[dGF[,5] == 22,5] <- 2
dGF[dGF[,5] == 11,5] <- 1
dGF <- dGF[order(dGF$Chr, dGF$Position),]
#dGF <- dGF[rowSums((dGF[,4:9]))!=0,]
rownames(dGF) <- NULL




#################################
####### set up comparison #######
#################################
comp <- as.data.frame(matrix(c(4,4,5,5,6,8,7,9),ncol=2))
comp <- rbind(comp,cbind(comp[,2],comp[,1]))
comp[,3] <- c("E04_GP_GF","E04_GP_GM", "E07_GP_GF", "E07_GP_GM","E04_GF_GP","E04_GM_GP", "E07_GF_GP", "E07_GM_GP")
for (c in 1:nrow(comp)) {
	k <- 0
	RBs <- as.data.frame(matrix(rep(NA,2*(flank*2+1)),ncol=2))
	RBsnp <- dGF[0,]
	for (chr in Chroms) {
		E <- dGF[dGF[,2] == chr,]
		RB <- as.numeric(rownames(E[diff(E[,comp[c,1]]) %in% c(1,-1,2,-2),]))
		RB_begin_end <- RB[c(1,length(RB))]
		RB <- setdiff(RB,RB_begin_end)

		RB_remove <- c()
		RB_real <- c()
		if (length(RB) > 1) {
			for (r in 1:(length(RB)-1)) {
				if (unique(dGF[(RB[r]+1):RB[r+1],comp[c,1]]) == 0) {
					RB_remove <- c(RB_remove,RB[r])
				}
			}

			RB_real <- setdiff(RB,RB_remove)
		} else {
			RB_real <- RB
		}

		if (length(RB_real) > 0) {
			RBsnp <- rbind(RBsnp,dGF[RB_real,])

			for (r in 1:length(RB_real)) {
				k = k+1
				pair_comp <- as.data.frame(matrix(rep(NA,2*(flank*2+1)),ncol=2))
				pair_comp[,1] <- dGF[(RB_real[r]-flank):(RB_real[r]+flank),comp[c,1]]
				pair_comp[,2] <- dGF[(RB_real[r]-flank):(RB_real[r]+flank),comp[c,2]]
				pair_comp[pair_comp[,1] == pair_comp[,2],3] <- 2
				pair_comp[(pair_comp[,1] != pair_comp[,2]) & (pair_comp[,1]==0 | pair_comp[,2]==0),3] <- 1
				pair_comp[(pair_comp[,1] != pair_comp[,2]) & (pair_comp[,1]!=0 & pair_comp[,2]!=0),3] <- 0
				RBs[,k] <- pair_comp[,3]
			}
		}
	}

	fn <- paste0(DIR,Family,"/",comp[c,3],".csv")
	write.table(RBsnp,fn,quote=F,sep="\t",col.names=T,row.names=F)

w1 <- (rowSums(RBs)/ncol(RBs)/2)
error <- qt(0.995,df=length(w1)-1)*sd(w1)/sqrt(length(w1))	# 99% confidence
left <- mean(w1)-error	# 0.9256324
right <- mean(w1)+error	# 0.9486329
#> mean(w1)
#[1] 0.9371326
#> sd(w1)
#[1] 0.07696653
# There is a 99% probability that the true mean is between 0.9256324 and 0.9486329

	fn <- paste0(DIR,Family,"/",comp[c,3],".png")
	png(fn)
	plot(rowSums(RBs)/ncol(RBs)/2,type="n",main=paste0(comp[c,3],"_", Family," with ", ncol(RBs), "RB"),xlab="SNPs flanking RB", ylab="Accuracy",xaxt = "n",ylim=c(0,1))
	axis(1,at=seq(1,(flank*2+1),by=50),labels=(seq(1,(flank*2+1),by=50)-flank-1))
	lines(rowSums(RBs)/ncol(RBs)/2)
	abline(h=0.99,col="red")
	abline(h=0.95,col="green")
	abline(v=(flank+1),col="red")
	text(flank*2+1,1,"0.99",cex=0.7)
	text(flank*2+1,0.94,"0.95",cex=0.7)
	dev.off()
}





