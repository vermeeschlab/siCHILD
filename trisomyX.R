##### part 1 ####

## [jding0@cnpew17 PGD]$ sed -e 's/,/./g' /uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/PGD263/PGD263.txt |cut -f 1-3,5-22 > /uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/PGD263/PGD263.adj

#library(GenABEL)
#g <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/PGD263/PGD263.adj", header=T, stringsAsFactors=F,sep="\t")
#g1 <-g[,6]
#g1[which(g1 == "AB")] <- 1
#g1[which(g1 == "AA")] <- 0
#g1[which(g1 == "BB")] <- 2
#g1[which(g1 == "NC")] <- "NA"
#P1 <- as.numeric(g1)
#g1 <-g[,9]
#g1[which(g1 == "AB")] <- 1
#g1[which(g1 == "AA")] <- 0
#g1[which(g1 == "BB")] <- 2
#g1[which(g1 == "NC")] <- "NA"
#P2 <- as.numeric(g1)
#off1<-generateOffspring(P1,P2)
#off2<-generateOffspring(P1,P2)
#baf <- off1 
#baf[which(baf %in% NA)] <- runif(1,0,1)
#baf[which(baf %in% 0)] <- runif(1,0,0.1)
#baf[which(baf %in% 2)] <- runif(1,0.9,1)
#baf[which(baf %in% 1)] <- runif(1,0.45,0.55)
#baf1<- baf
#baf <- off2
#baf[which(baf %in% NA)] <- runif(1,0,1)
#baf[which(baf %in% 0)] <- runif(1,0,0.1)
#baf[which(baf %in% 2)] <- runif(1,0.9,1)
#baf[which(baf %in% 1)] <- runif(1,0.45,0.55)
#baf2 <- baf

#g[c(ncol(g)+1):c(ncol(g)+6)] <- cbind(baf1,sample(g[,5]),off1,baf2,sample(g[,5]),off2)
#colnames(g)[(ncol(g)-5):ncol(g)] <- c("E01_Bl001.B.Allele.Freq","E01_Bl001.Log.R.Ratio","E01_Bl001.GType","E02_Bl001.B.Allele.Freq","E02_Bl001.Log.R.Ratio","E02_Bl001.GType")
#g[g[,"E01_Bl001.GType"] %in% 1,"E01_Bl001.GType"] <- "AB"
#g[g[,"E01_Bl001.GType"] %in% 2,"E01_Bl001.GType"] <- "BB"
#g[g[,"E01_Bl001.GType"] %in% 0,"E01_Bl001.GType"] <- "AA"
#g[g[,"E01_Bl001.GType"] %in% NA,"E01_Bl001.GType"] <- "NC"
#g[g[,"E02_Bl001.GType"] %in% 1,"E02_Bl001.GType"] <- "AB"
#g[g[,"E02_Bl001.GType"] %in% 2,"E02_Bl001.GType"] <- "BB"
#g[g[,"E02_Bl001.GType"] %in% 0,"E02_Bl001.GType"] <- "AA"
#g[g[,"E02_Bl001.GType"] %in% NA,"E02_Bl001.GType"] <- "NC"
#write.table(g,"/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/PGD263/PGD263.adj",sep="\t",row.names=F,quote=F)





## [jding0@cnpew06 PGD]$ /cm/shared/apps/R/3.2.4/bin/R --save --args PGD263 Grandparents v2 E*_Bl* /uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/ < /home/jding0/siCHILD_toGC/siCHILD_core.R

##    0     1    10    11    12    13    14    15    16    17    18    19     2 
## 1082 23066 13719 14347 14591  8994  9591 10068  9229  9641  7605  6883 24117 
##   20    21    22     3     4     5     6     7     8     9     X    XY     Y 
## 7596  4327  4988 17733 14464 17821 16513 15657 15779 12207 15318   766  2461 


#### part 2 ####
chr <- c(1:22,"X")
pgd <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/PGD263/PGD263_old.adj",sep="\t",header=T,stringsAsFactors=F)
ht <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/PGD263_out_Grandparents_v1/PGD263_0.75.haplotype",sep="\t",stringsAsFactors=F,header=T)

print("E01")
pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "BB" & ht$Mother_PGD263 == "AB"),]$Name,]$E01_Bl001.GType <- "AB"
pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "BB" & ht$Mother_PGD263 == "BA"),]$Name,]$E01_Bl001.GType <- "BB"
pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "AA" & ht$Mother_PGD263 == "BA"),]$Name,]$E01_Bl001.GType <- "AB"
pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "AA" & ht$Mother_PGD263 == "AB"),]$Name,]$E01_Bl001.GType <- "AA"
pgd[pgd$E01_Bl001.GType == "AA",]$E01_Bl001.B.Allele.Freq <- runif(nrow(pgd[pgd$E01_Bl001.GType == "AA",]),0,0.1)
pgd[pgd$E01_Bl001.GType == "BB",]$E01_Bl001.B.Allele.Freq <- runif(nrow(pgd[pgd$E01_Bl001.GType == "BB",]),0.9,1)
pgd[pgd$E01_Bl001.GType == "AB",]$E01_Bl001.B.Allele.Freq <- runif(nrow(pgd[pgd$E01_Bl001.GType == "AB",]),0.45,0.55)
pgd_org <- pgd
pgd_new <- c()

bp <- 0.5	# first 0.5 is H2, second half is H1
for (c in chr) {
	pgd <- pgd_org[pgd_org$Chr == c,][1:floor(nrow(pgd_org[pgd_org$Chr == c,])*bp),]
	pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "BB" & ht$Mother_PGD263 == "AB"),]$Name,]$E01_Bl001.GType <- "BB"
	pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "BB" & ht$Mother_PGD263 == "BA"),]$Name,]$E01_Bl001.GType <- "AB"
	pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "AA" & ht$Mother_PGD263 == "BA"),]$Name,]$E01_Bl001.GType <- "AA"
	pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "AA" & ht$Mother_PGD263 == "AB"),]$Name,]$E01_Bl001.GType <- "AB"
	pgd[pgd$E01_Bl001.GType == "AA",]$E01_Bl001.B.Allele.Freq <- runif(nrow(pgd[pgd$E01_Bl001.GType == "AA",]),0,0.1)
	pgd[pgd$E01_Bl001.GType == "BB",]$E01_Bl001.B.Allele.Freq <- runif(nrow(pgd[pgd$E01_Bl001.GType == "BB",]),0.9,1)
	pgd[pgd$E01_Bl001.GType == "AB",]$E01_Bl001.B.Allele.Freq <- runif(nrow(pgd[pgd$E01_Bl001.GType == "AB",]),0.45,0.55)
	pgd_new <- rbind(pgd_new,pgd)
	pgd_new <- rbind(pgd_new,pgd_org[pgd_org$Chr == c,][((floor(nrow(pgd_org[pgd_org$Chr == c,])*bp)+1):nrow(pgd_org[pgd_org$Chr ==c,])),])
}
pgd <- rbind(pgd_new,pgd_org[pgd_org$Chr =="Y",])
pn <- sample(pgd$Name,length(pgd$Name)*0.2)
pgd[pgd$Name %in% pn,]$E01_Bl001.B.Allele.Freq <- runif(length(pn),0,1)	# introduce 20% noise to BAF
length(unique(pgd$Name))


print("E02")
pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "BB" & ht$Mother_PGD263 == "AB"),]$Name,]$E02_Bl001.GType <- "AB"
pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "BB" & ht$Mother_PGD263 == "BA"),]$Name,]$E02_Bl001.GType <- "BB"
pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "AA" & ht$Mother_PGD263 == "BA"),]$Name,]$E02_Bl001.GType <- "AB"
pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "AA" & ht$Mother_PGD263 == "AB"),]$Name,]$E02_Bl001.GType <- "AA"
pgd[pgd$E02_Bl001.GType == "AA",]$E02_Bl001.B.Allele.Freq <- runif(nrow(pgd[pgd$E02_Bl001.GType == "AA",]),0,0.1)
pgd[pgd$E02_Bl001.GType == "BB",]$E02_Bl001.B.Allele.Freq <- runif(nrow(pgd[pgd$E02_Bl001.GType == "BB",]),0.9,1)
pgd[pgd$E02_Bl001.GType == "AB",]$E02_Bl001.B.Allele.Freq <- runif(nrow(pgd[pgd$E02_Bl001.GType == "AB",]),0.45,0.55)
pgd_org <- pgd
pgd_new <- c()
bp <- 0.75
for (c in chr) {
	pgd <- pgd_org[pgd_org$Chr == c,][1:floor(nrow(pgd_org[pgd_org$Chr == c,])*bp),]
	pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "BB" & ht$Mother_PGD263 == "AB"),]$Name,]$E02_Bl001.GType <- "BB"
	pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "BB" & ht$Mother_PGD263 == "BA"),]$Name,]$E02_Bl001.GType <- "AB"
	pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "AA" & ht$Mother_PGD263 == "BA"),]$Name,]$E02_Bl001.GType <- "AA"
	pgd[pgd$Name %in% ht[(ht$Father_PGD263 == "AA" & ht$Mother_PGD263 == "AB"),]$Name,]$E02_Bl001.GType <- "AB"
	pgd[pgd$E02_Bl001.GType == "AA",]$E02_Bl001.B.Allele.Freq <- runif(nrow(pgd[pgd$E02_Bl001.GType == "AA",]),0,0.1)
	pgd[pgd$E02_Bl001.GType == "BB",]$E02_Bl001.B.Allele.Freq <- runif(nrow(pgd[pgd$E02_Bl001.GType == "BB",]),0.9,1)
	pgd[pgd$E02_Bl001.GType == "AB",]$E02_Bl001.B.Allele.Freq <- runif(nrow(pgd[pgd$E02_Bl001.GType == "AB",]),0.45,0.55)
	pgd_new <- rbind(pgd_new,pgd)
	pgd_new <- rbind(pgd_new,pgd_org[pgd_org$Chr == c,][((floor(nrow(pgd_org[pgd_org$Chr == c,])*bp)+1):nrow(pgd_org[pgd_org$Chr ==c,])),])
}
pgd <- rbind(pgd_new,pgd_org[pgd_org$Chr =="Y",])
pn <- sample(pgd$Name,length(pgd$Name)*0.2)
pgd[pgd$Name %in% pn,]$E01_Bl001.B.Allele.Freq <- runif(length(pn),0,1)	# introduce 20% noise to BAF
length(unique(pgd$Name))

write.table(pgd,"/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/PGD263/PGD263.adj",sep="\t",row.names=F,quote=F)





