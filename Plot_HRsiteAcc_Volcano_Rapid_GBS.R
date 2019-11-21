
rm(list = ls(all = T))

#args <- commandArgs(TRUE)
#File1 <- args[1]	# F1: "PGD_Hapmap_out_Grandparents_f2/PGD_Hapmap_Itp.hap"
#File2 <- args[2]	# F2: "PGD_Hapmap_out_Grandparents_f2/PGD_Hapmap_Raw.hap"

pgdDIRs <- c("PGD_Hapmap_GBSdataHeleen_GM12882MC", "PGD_Hapmap_GBSdataHeleen")
AccInt <- c(250, 80)

dataPath <- "/uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/PGD_HapMap/"
opt <- as.data.frame(matrix(c("PGD_Hapmap_GBSdataHeleen_GM12882MC", "PGD_Hapmap_GBSdataHeleen", "GM12882_MC", "GM12887_MC","GM12887", "GM12882"), ncol = 3))
colnames(opt) <- c("DIR", "Seed", "EmbryoToCheck")



Pars <- matrix(c("Mat","Pat","Mother","Father","Father","Mother"), nrow = 2, ncol = 3)
for (p in 1:nrow(Pars)) {
	par <- Pars[p,1]
	print(par)

	comb_HRsiteAccItp_SCs <- NULL
	comb_HRsiteAccRaw_SCs <- NULL


	for (pgdDIR in pgdDIRs) {
	print(pgdDIR)
    outPath <- paste0(dataPath, pgdDIR, "/")

    File1 <- paste0(pgdDIR, "/result/", pgdDIR, "_Itp2.hap")
    File2 <- paste0(pgdDIR, "/result/", pgdDIR, "_Raw.hap")
    File3 <- paste0(pgdDIR, "/result/", pgdDIR, "_0.75.gtp")
    #Interval for plotting

    #------------------------------------------------


    dataItp <- read.table(paste(dataPath, File1, sep = "") , sep = "\t", header = T, stringsAsFactors = F)
    dataRaw <- read.table(paste(dataPath, File2, sep = "") , sep = "\t", header = T, stringsAsFactors = F)
    gtp <- read.table(paste(dataPath, File3, sep = "") , sep = "\t", header = T, stringsAsFactors = F)
    ID <- read.table(paste(dataPath, "id.csv", sep = "") , sep = "\t", header = T, stringsAsFactors = F)

    SNPChrPos <- dataItp[, 1 : 3]
    rownames(dataItp) <- SNPChrPos[, 1]
    rownames(dataRaw) <- SNPChrPos[, 1]
    Chroms <- as.character(1 : 22)



    dataItp <- dataItp[, c(1 : 3, grep(par, colnames(dataItp)))]
    dataRaw <- dataRaw[, c(1 : 3, grep(par, colnames(dataRaw)))]
    colnames(dataItp)[4 : ncol(dataItp)] <- gsub(paste0("_",pgdDIR,"_", par), "", colnames(dataItp)[4 : ncol(dataItp)])
    colnames(dataRaw)[4 : ncol(dataRaw)] <- gsub(paste0("_",pgdDIR,"_", par), "", colnames(dataRaw)[4 : ncol(dataRaw)])

	for (c in which(colnames(dataItp) %in% ID[, 2])) {
		colnames(dataItp)[c] <- ID[ID[, 2] %in% colnames(dataItp)[c],1]
	}

	for (c in which(colnames(dataRaw) %in% ID[, 2])) {
		colnames(dataRaw)[c] <- ID[ID[, 2] %in% colnames(dataRaw)[c],1]
	}


	dataItp <- dataItp[,c(1:3,grep(as.character(opt[opt[, 1] %in% pgdDIR, 3]),colnames(dataItp)))]
	dataRaw <- dataRaw[,c(1:3,grep(as.character(opt[opt[, 1] %in% pgdDIR, 3]),colnames(dataRaw)))]


    #	MCRef <- as.data.frame(dataItp[,grep(paste0(as.character(opt[opt[,1] %in% pgdDIR,3]),"_MC"),colnames(dataItp))])
    #	colnames(MCRef) <- colnames(dataItp)[grep(paste0(as.character(opt[opt[,1] %in% pgdDIR,3]),"_MC"),colnames(dataItp))]
    for (run in 1 : 1) { # run=1 => accuracy with all snps; run=2 => accuracy with informative snps.
            if (run == 2) {
                rn <- gtp[gtp[, grep(Pars[p,2], colnames(gtp))] == "AB" &
                    gtp[, grep(Pars[p,3], colnames(gtp))] != "AB" &
                    gtp[, grep(Pars[p,3], colnames(gtp))] != "NC", 1]
                dataItp <- dataItp[dataItp$Name %in% rn,]
                dataRaw <- dataRaw[dataRaw$Name %in% rn,]
            }


        MCRef <- as.data.frame(dataItp[, grep(paste(paste0(as.character(opt[opt[, 1] %in% pgdDIR, 3]), "_MC"), collapse = "|"), colnames(dataItp))])
        colnames(MCRef) <- colnames(dataItp)[grep(paste(paste0(as.character(opt[opt[, 1] %in% pgdDIR, 3]), "_MC"), collapse = "|"), colnames(dataItp))]

        dataItp$Chr <- as.character(dataItp$Chr)

        #Determining HR-sites in the MC reference

        for (mc in 1 : ncol(MCRef)) {
            main_text <- c(paste0("phasing with ", as.character(opt[opt[, 1] %in% pgdDIR, 2])[mc], " & checking accuracy of ", as.character(opt[opt[, 1] %in% pgdDIR, 3])[mc], " with all SNPs. "), paste0("phasing with ", as.character(opt[opt[, 1] %in% pgdDIR, 2])[mc], " & checking accuracy of ", as.character(opt[opt[, 1] %in% pgdDIR, 3])[mc], " with informative SNPs. "))

            SCItp <- dataItp[, grep(gsub("MC", "SC", colnames(MCRef)[mc]), colnames(dataItp))]
            SCRaw <- dataRaw[, grep(gsub("MC", "SC", colnames(MCRef)[mc]), colnames(dataRaw))]
            print(colnames(MCRef)[mc])
            HRsites <- NULL
            for (chr in Chroms) {
                dataChr <- as.data.frame(MCRef[dataItp$Chr == chr,])
                SNPChrPosChr <- dataItp[dataItp$Chr == chr, 1 : 3]
                Blk <- rle(dataChr[, mc])
                Blk1 <- as.data.frame(cbind(Blk$values, Blk$lengths, cumsum(Blk$lengths), matrix(0, length(Blk$values), 1)))
                Blk1 <- Blk1[Blk1[, 1] != 0,]
                Blk2 <- Blk1[Blk1[, 2] >= 5,]
                if (is.vector(Blk2)) {print(paste("No HRsite for Chromosome", chr))
                } else if (nrow(Blk2) == 1) {print(paste("One haplotype, No HRsite for Chromosome", chr))
                } else {
                    for (i in 1 : (nrow(Blk2) - 1)) {if (Blk2[i, 2] >= AccInt[run] & Blk2[i + 1, 2] >= AccInt[run]) { Blk2[i, 4] <- 1}}#end i loop
                    HRsite <- cbind(SNPChrPosChr[Blk2[Blk2[, 4] == 1, 3], 2 : 3], Blk2[Blk2[, 4] == 1, 3])
                    rownames(HRsite) <- SNPChrPosChr[Blk2[Blk2[, 4] == 1, 3], 1]
                    HRsites <- rbind(HRsites, HRsite)
                }
            }#end chr loop

            rownames(dataItp) <- as.character(dataItp$Name)
            HRPos <- dataItp[rownames(HRsites),]

            for (i in 1 : nrow(HRsites)) {
                HRPos <- which(rownames(dataItp) == rownames(HRsites)[i])
                HRsiteAcc_Ref <- MCRef[c((HRPos - AccInt[run]) : (HRPos - 1), (HRPos + 1) : (HRPos + AccInt[run])), mc]
                HRsiteItp_SC <- SCItp[c((HRPos - AccInt[run]) : (HRPos - 1), (HRPos + 1) : (HRPos + AccInt[run])), grep("SC", colnames(SCItp))]
                HRsiteRaw_SC <- SCRaw[c((HRPos - AccInt[run]) : (HRPos - 1), (HRPos + 1) : (HRPos + AccInt[run])), grep("SC", colnames(SCItp))]

                AccItp_SCs <- matrix(NA, ncol(HRsiteItp_SC), AccInt[run] * 2)
                AccRaw_SCs <- matrix(NA, ncol(HRsiteItp_SC), AccInt[run] * 2)

                for (sc in 1 : ncol(HRsiteItp_SC)) {
                    AccItp_SCs[sc,] <- HRsiteAcc_Ref == HRsiteItp_SC[, sc] & HRsiteAcc_Ref != 0
                    AccRaw_SCs[sc,] <- HRsiteAcc_Ref == HRsiteRaw_SC[, sc] & HRsiteAcc_Ref != 0
                }#end sc loop

                if (i == 1) {
                    HRsiteAccItp_SCs <- AccItp_SCs
                    HRsiteAccRaw_SCs <- AccRaw_SCs
                } else {
                    HRsiteAccItp_SCs <- rbind(HRsiteAccItp_SCs, AccItp_SCs)
                    HRsiteAccRaw_SCs <- rbind(HRsiteAccRaw_SCs, AccRaw_SCs)
                }
            }#end i loop

			comb_HRsiteAccItp_SCs <- rbind(comb_HRsiteAccItp_SCs,HRsiteAccItp_SCs)
			comb_HRsiteAccRaw_SCs <- rbind(comb_HRsiteAccRaw_SCs,HRsiteAccRaw_SCs)

            AccItpSC <- (apply(HRsiteAccItp_SCs, 2, sum) / nrow(HRsiteAccItp_SCs)) * 100
            AccRawSC <- (apply(HRsiteAccRaw_SCs, 2, sum) / nrow(HRsiteAccItp_SCs)) * 100

            outfile <- paste(outPath, "HRsiteAcc_", par, colnames(MCRef)[mc], "_run", run, ".pdf", sep = "")
            # outfile <- paste(outPath,"HRsiteAcc_",par,colnames(MCRef)[mc],"_run",run,".tiff",sep="")
            # outfile <- paste(outPath,"HRsiteAcc_",par,colnames(MCRef)[mc],"_run",run,".png",sep="")
            print(outfile)
            pdf(outfile, w = 5, h = 5, bg = "transparent", family = "Helvetica", title = "RB plot")
            # tiff(outfile,width = 10, height = 6,units = "in",res=300)
            # png(outfile)
            # par(mfrow=c(1,2))

            plot(c(- AccInt[run] : - 1, 1 : AccInt[run]), matrix(0, 1, AccInt[run] * 2), xlab = "SNP Position", ylab = "HR-detection accuracy (%)", pch = "*", ylim = c(0, 100), frame = FALSE, cex.main = 0.4, main = main_text[run])

            if (run == 2) {
                points(c(- AccInt[run] : - 1, 1 : AccInt[run]), AccItpSC, "h", col = "#ff000080")
            } else {
                points(c(- AccInt[run] : - 1, 1 : AccInt[run]), AccItpSC, "h", col = "#ff000020")
                points(c(- AccInt[run] : - 1, 1 : AccInt[run]), AccRawSC, "h", col = "#ff000080")
            }


            abline(h = 99, col = "grey")
            abline(h = 95, col = "light green")
            #abline(v=500,lty=2,col="grey",lwd=2)
            dev.off()
        }
        }#end par loop
    }



            AccItpSC <- (apply(comb_HRsiteAccItp_SCs, 2, sum) / nrow(comb_HRsiteAccItp_SCs)) * 100
            AccRawSC <- (apply(comb_HRsiteAccRaw_SCs, 2, sum) / nrow(comb_HRsiteAccItp_SCs)) * 100

            outfile <- paste(outPath, "HRsiteAcc_", par, colnames(MCRef)[mc], "_comb_run", run, ".pdf", sep = "")
            # outfile <- paste(outPath,"HRsiteAcc_",par,colnames(MCRef)[mc],"_run",run,".tiff",sep="")
            # outfile <- paste(outPath,"HRsiteAcc_",par,colnames(MCRef)[mc],"_run",run,".png",sep="")
            print(outfile)
            pdf(outfile, w = 5, h = 5, bg = "transparent", family = "Helvetica", title = "RB plot")
            # tiff(outfile,width = 10, height = 6,units = "in",res=300)
            # png(outfile)
            # par(mfrow=c(1,2))

            plot(c(- AccInt[run] : - 1, 1 : AccInt[run]), matrix(0, 1, AccInt[run] * 2), xlab = "SNP Position", ylab = "HR-detection accuracy (%)", pch = "*", ylim = c(0, 100), frame = FALSE, cex.main = 0.4, main = main_text[run])

            if (run == 2) {
                points(c(- AccInt[run] : - 1, 1 : AccInt[run]), AccItpSC, "h", col = "#ff000080")
            } else {
                points(c(- AccInt[run] : - 1, 1 : AccInt[run]), AccItpSC, "h", col = "#ff000020")
                points(c(- AccInt[run] : - 1, 1 : AccInt[run]), AccRawSC, "h", col = "#ff000080")
            }


            abline(h = 99, col = "grey")
            abline(h = 95, col = "light green")
            #abline(v=500,lty=2,col="grey",lwd=2)
            dev.off()



}


# Illumina CytoSNP12-v2.1










