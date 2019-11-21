id <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/familyID.txt",header=T,stringsAsFactors=F)


## compare the number of inforamtive snps

### extended family
itp2_file <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/itp2.txt",header=F,stringsAsFactors=F)
itp_file <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/itp.txt",header=F,stringsAsFactors=F)
raw_file <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/raw.txt",header=F,stringsAsFactors=F)
DIR <- "/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/"


#DIRGP <- list.dirs(path = DIR, pattern = "f0", full.names = TRUE, recursive = TRUE)
#DIRext <- list.dirs(path = DIR, pattern = "f4", full.names = TRUE, recursive = TRUE)

nritp2_info_snp <- c()
nritp_info_snp <- c()
nrraw_info_snp <- c()
for (i in 1:nrow(itp2_file)) {
	itp2 <- read.table(paste0(DIR,"/",itp2_file[i,]),header=T,stringsAsFactors=F)
	itp2_info_snp <- itp2[,4:ncol(itp2)]
	itp2_info_snp[itp2_info_snp == 0] <- NA
	nritp2_info_snp[i] <- nrow(itp2_info_snp[rowSums(is.na(itp2_info_snp))!=ncol(itp2_info_snp), ])

	itp <- read.table(paste0(DIR,"/",itp_file[i,]),header=T,stringsAsFactors=F)
	itp_info_snp <- itp[,4:ncol(itp)]
	itp_info_snp[itp_info_snp == 0] <- NA
	nritp_info_snp[i] <- nrow(itp_info_snp[rowSums(is.na(itp_info_snp))!=ncol(itp_info_snp), ])

	raw <- read.table(paste0(DIR,"/",raw_file[i,]),header=T,stringsAsFactors=F)
	raw_info_snp <- raw[,4:ncol(raw)]
	raw_info_snp[raw_info_snp == 0] <- NA
	nrraw_info_snp[i] <- nrow(raw_info_snp[rowSums(is.na(raw_info_snp))!=ncol(raw_info_snp), ])
}

nr_info_snp <- cbind(nrraw_info_snp,nritp_info_snp,nritp2_info_snp)

write.table(nr_info_snp,"family12_nr_info_snp_IBS.csv",quote=F,sep="\t",col.names=T,row.names=F)



### GP
itp2_file <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/itp2_gp.txt",header=F,stringsAsFactors=F)
itp_file <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/itp_gp.txt",header=F,stringsAsFactors=F)
raw_file <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/raw_gp.txt",header=F,stringsAsFactors=F)
DIR <- "/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/"


#DIRGP <- list.dirs(path = DIR, pattern = "f0", full.names = TRUE, recursive = TRUE)
#DIRext <- list.dirs(path = DIR, pattern = "f4", full.names = TRUE, recursive = TRUE)

nritp2_info_snp <- c()
nritp_info_snp <- c()
nrraw_info_snp <- c()
for (i in 1:nrow(itp2_file)) {
	itp2 <- read.table(paste0(DIR,"/",itp2_file[i,]),header=T,stringsAsFactors=F)
	itp2_info_snp <- itp2[,4:ncol(itp2)]
	itp2_info_snp[itp2_info_snp == 0] <- NA
	nritp2_info_snp[i] <- nrow(itp2_info_snp[rowSums(is.na(itp2_info_snp))!=ncol(itp2_info_snp), ])

	itp <- read.table(paste0(DIR,"/",itp_file[i,]),header=T,stringsAsFactors=F)
	itp_info_snp <- itp[,4:ncol(itp)]
	itp_info_snp[itp_info_snp == 0] <- NA
	nritp_info_snp[i] <- nrow(itp_info_snp[rowSums(is.na(itp_info_snp))!=ncol(itp_info_snp), ])

	raw <- read.table(paste0(DIR,"/",raw_file[i,]),header=T,stringsAsFactors=F)
	raw_info_snp <- raw[,4:ncol(raw)]
	raw_info_snp[raw_info_snp == 0] <- NA
	nrraw_info_snp[i] <- nrow(raw_info_snp[rowSums(is.na(raw_info_snp))!=ncol(raw_info_snp), ])
}

nr_info_snp <- cbind(nrraw_info_snp,nritp_info_snp,nritp2_info_snp)

write.table(nr_info_snp,"family12_nr_info_snp_GP.csv",quote=F,sep="\t",col.names=T,row.names=F)


### single GF
itp2_file <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/itp2_gf.txt",header=F,stringsAsFactors=F)
itp_file <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/itp_gf.txt",header=F,stringsAsFactors=F)
raw_file <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/raw_gf.txt",header=F,stringsAsFactors=F)
DIR <- "/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/"


#DIRGP <- list.dirs(path = DIR, pattern = "f0", full.names = TRUE, recursive = TRUE)
#DIRext <- list.dirs(path = DIR, pattern = "f4", full.names = TRUE, recursive = TRUE)

nritp2_info_snp <- c()
nritp_info_snp <- c()
nrraw_info_snp <- c()
for (i in 1:nrow(itp2_file)) {
	itp2 <- read.table(paste0(DIR,"/",itp2_file[i,]),header=T,stringsAsFactors=F)
	itp2_info_snp <- itp2[,4:ncol(itp2)]
	itp2_info_snp[itp2_info_snp == 0] <- NA
	nritp2_info_snp[i] <- nrow(itp2_info_snp[rowSums(is.na(itp2_info_snp))!=ncol(itp2_info_snp), ])

	itp <- read.table(paste0(DIR,"/",itp_file[i,]),header=T,stringsAsFactors=F)
	itp_info_snp <- itp[,4:ncol(itp)]
	itp_info_snp[itp_info_snp == 0] <- NA
	nritp_info_snp[i] <- nrow(itp_info_snp[rowSums(is.na(itp_info_snp))!=ncol(itp_info_snp), ])

	raw <- read.table(paste0(DIR,"/",raw_file[i,]),header=T,stringsAsFactors=F)
	raw_info_snp <- raw[,4:ncol(raw)]
	raw_info_snp[raw_info_snp == 0] <- NA
	nrraw_info_snp[i] <- nrow(raw_info_snp[rowSums(is.na(raw_info_snp))!=ncol(raw_info_snp), ])
}

nr_info_snp <- cbind(nrraw_info_snp,nritp_info_snp,nritp2_info_snp)

write.table(nr_info_snp,"family12_nr_info_snp_GF.csv",quote=F,sep="\t",col.names=T,row.names=F)





### single GM
itp2_file <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/itp2_gm.txt",header=F,stringsAsFactors=F)
itp_file <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/itp_gm.txt",header=F,stringsAsFactors=F)
raw_file <- read.table("/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/raw_gm.txt",header=F,stringsAsFactors=F)
DIR <- "/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/"


#DIRGP <- list.dirs(path = DIR, pattern = "f0", full.names = TRUE, recursive = TRUE)
#DIRext <- list.dirs(path = DIR, pattern = "f4", full.names = TRUE, recursive = TRUE)

nritp2_info_snp <- c()
nritp_info_snp <- c()
nrraw_info_snp <- c()
for (i in 1:nrow(itp2_file)) {
	itp2 <- read.table(paste0(DIR,"/",itp2_file[i,]),header=T,stringsAsFactors=F)
	itp2_info_snp <- itp2[,4:ncol(itp2)]
	itp2_info_snp[itp2_info_snp == 0] <- NA
	nritp2_info_snp[i] <- nrow(itp2_info_snp[rowSums(is.na(itp2_info_snp))!=ncol(itp2_info_snp), ])

	itp <- read.table(paste0(DIR,"/",itp_file[i,]),header=T,stringsAsFactors=F)
	itp_info_snp <- itp[,4:ncol(itp)]
	itp_info_snp[itp_info_snp == 0] <- NA
	nritp_info_snp[i] <- nrow(itp_info_snp[rowSums(is.na(itp_info_snp))!=ncol(itp_info_snp), ])

	raw <- read.table(paste0(DIR,"/",raw_file[i,]),header=T,stringsAsFactors=F)
	raw_info_snp <- raw[,4:ncol(raw)]
	raw_info_snp[raw_info_snp == 0] <- NA
	nrraw_info_snp[i] <- nrow(raw_info_snp[rowSums(is.na(raw_info_snp))!=ncol(raw_info_snp), ])
}

nr_info_snp <- cbind(nrraw_info_snp,nritp_info_snp,nritp2_info_snp)

write.table(nr_info_snp,"family12_nr_info_snp_GM.csv",quote=F,sep="\t",col.names=T,row.names=F)











