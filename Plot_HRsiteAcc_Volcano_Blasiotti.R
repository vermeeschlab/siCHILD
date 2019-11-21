#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%       				Volcano Like Plot (HR)          			   	  %%%%%%%%%%%%%#
#
#Author: Jia Ding (modified from MZE)
#
#cDate: 13.March.2017, original version 01/01/2012
#mDtae:
#
#(->):  Read in reconstructed haplotypes
#
#(<-): Volcano like plots
#
#Description: This shows accuracy of determined HR-stites for WGA-Haplotypes before and after interpretation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rm(list = ls(all = T))
# args <- commandArgs(TRUE)
dataPath <- "/uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/"
opt_fn <- "PGD_Blasiotti_org/HR.csv"
Chroms <- as.character(1 : 22)
Pars <- matrix(c("Mat", "Pat", "Mother", "Father", "Father", "Mother"), nrow = 2, ncol = 3)
eid <- c("E01", "E02")
AccInt <- c(250, 80)
snps <- c("all", "informative") # run=1 => accuracy with all snps; run=2 => accuracy with informative snps.



ACCItp_all <- list()


opt <- read.table(paste0(dataPath, opt_fn), header = T, stringsAsFactors = F)

for (o in 1 : nrow(opt)) {
    pgdDIR <- paste0(dataPath, opt[o,]$Dir, opt[o,]$subDir)
    print(paste0("----", pgdDIR, "----"))
    File1 <- paste0(pgdDIR, "/", opt[o,]$Dir, "_Itp2.hap")
    File2 <- paste0(pgdDIR, "/", opt[o,]$Dir, "_Raw.hap")
    File3 <- paste0(pgdDIR, "/", opt[o,]$Dir, "_0.75.gtp")

    dataItp <- read.table(File1 , sep = "\t", header = T, stringsAsFactors = F)
    dataRaw <- read.table(File2 , sep = "\t", header = T, stringsAsFactors = F)
    gtp <- read.table(File3, sep = "\t", header = T, stringsAsFactors = F)

    ibs1_fn <- list.files(path = paste0(pgdDIR, "../"), pattern = "olaps1")
    print(ibs1_fn)
    if (! identical(ibs1_fn, character(0))) {
        ibs1 <- read.table(paste0(pgdDIR, "../", ibs1_fn), sep = "\t", header = T, stringsAsFactors = F)
        dataItp <- dataItp[(dataItp)$Name %in% ibs1[, 1],]
        dataRaw <- dataRaw[(dataRaw)$Name %in% ibs1[, 1],]
        gtp <- gtp[(gtp)$Name %in% ibs1[, 1],]
    }

    SNPChrPos <- dataItp[, 1 : 3]
    rownames(dataItp) <- SNPChrPos[, 1]
    rownames(dataRaw) <- SNPChrPos[, 1]


    #exc_col <- c(grep("E02_Bl001",colnames(dataItp)),grep("E02_Bl002",colnames(dataItp)))
    #dataItp <- dataItp[,-exc_col]
    #exc_col <- c(grep("E02_Bl001",colnames(dataRaw)),grep("E02_Bl002",colnames(dataRaw)))
    #dataRaw <- dataItp[,-exc_col]
    #exc_col <- c(grep("E02_Bl001",colnames(gtp)),grep("E02_Bl002",colnames(gtp)))
    #gtp <- gtp[,-exc_col]

    dataItp_org <- dataItp
    dataRaw_org <- dataRaw
    gtp_org <- gtp

    for (p in 1 : nrow(Pars)) {
        par <- Pars[p, 1]
        dataItp <- dataItp_org
        dataRaw <- dataRaw_org
        gtp <- gtp_org
        print(paste0("-------------------",par,"------------------"))


        dataItp <- dataItp[, c(1 : 3, grep(par, colnames(dataItp)))]
        dataRaw <- dataRaw[, c(1 : 3, grep(par, colnames(dataRaw)))]
        # print(paste0("=======",pgdDIR,"=======",nrow(dataItp)))

        for (run in 1 : 1) {
            if (run == 2) {
                rn <- gtp[gtp[, grep(Pars[p, 2], colnames(gtp))] == "AB" &
                    gtp[, grep(Pars[p, 3], colnames(gtp))] != "AB" &
                    gtp[, grep(Pars[p, 3], colnames(gtp))] != "NC", 1]
                dataItp <- dataItp[dataItp$Name %in% rn,]
                dataRaw <- dataRaw[dataRaw$Name %in% rn,]
            }


            for (e in 1 : length(eid)) {
                print(eid[e])
                x <- (c(grep(eid[e], colnames(dataItp)), grep("_Bl01", colnames(dataItp))))
                mc_cn <- x[duplicated(x)]
                MCRef <- as.data.frame(dataItp[, mc_cn])
                x <- (c(grep(eid[e], colnames(dataItp)), grep("_Bl00", colnames(dataItp))))
                sc_cn <- x[duplicated(x)]

                # print(paste0(run,"------",pgdDIR,"------",nrow(dataItp)))
                embryoID <- paste(unlist(strsplit(colnames(MCRef)[1],"_"))[1],unlist(strsplit(colnames(MCRef)[1],"_"))[5],unlist(strsplit(colnames(MCRef)[1],"_"))[6],sep="_")
                main_text <- paste0("checking accuracy of ", embryoID , " with ", snps[run], " SNPs. ")
                print(main_text)



                HRsiteAccItp_SCs_bothSide <- NULL
                HRsiteAccRaw_SCs_bothSide <- NULL
                HR_list <- NULL
                for (mc in 1 : ncol(MCRef)) {
                    #Determining HR-sites in each of MC reference
                    SCItp <- dataItp[, sc_cn]
                    SCRaw <- dataRaw[, sc_cn]

                    HR_left <- NULL
                    HR_right <- NULL

                    for (chr in Chroms) {
                        dataChr <- as.data.frame(MCRef[dataItp$Chr == chr,])
                        SNPChrPosChr <- dataItp[dataItp$Chr == chr, 1 : 3]
                        Blk <- rle(dataChr[, mc])
                        if ((length(Blk$values) == 1 ) && Blk$value == 0) {break}

                        Blk1 <- as.data.frame(cbind(Blk$values, Blk$lengths, cumsum(Blk$lengths), matrix(0, length(Blk$values), 1)))
                        Blk1 <- Blk1[Blk1[, 1] != 0,]

                        ## HR site from left side
                        Blk2 <- Blk1[Blk1[, 2] >= 5,]

                        ## HR site from right side
                        Blk3 <- cbind(Blk2[, 1], (Blk2[, 2] + 1), (Blk2[, 3] + 1), Blk2[, 4])
                        colnames(Blk3) <- colnames(Blk2)
                        rownames(Blk3) <- rownames(Blk2)

                        ## HR site from left side
                        # print("left side")
                        if (is.vector(Blk2)) {print(paste("No HRsite for Chromosome", chr))
                        } else if (nrow(Blk2) == 1) {print(paste("One haplotype, No HRsite for Chromosome", chr))
                        } else {
                            for (i in 1 : (nrow(Blk2) - 1)) {if (Blk2[i, 2] >= AccInt[run] & Blk2[i + 1, 2] >= AccInt[run]) { Blk2[i, 4] <- 1}}#end i loop
                            HRsite <- cbind(SNPChrPosChr[Blk2[Blk2[, 4] == 1, 3], 2 : 3], Blk2[Blk2[, 4] == 1, 3], rep(mc,length(Blk2[Blk2[, 4] == 1, 3])))
                            rownames(HRsite) <- SNPChrPosChr[Blk2[Blk2[, 4] == 1, 3], 1]
                            HR_left <- rbind(HR_left, HRsite)
                        }

                        # ## HR site from right side
                        # print("right side")
                        # if (is.vector(Blk3)) {print(paste("No HRsite for Chromosome", chr))
                        # } else if (nrow(Blk3) == 1) {print(paste("One haplotype, No HRsite for Chromosome", chr))
                        # } else {
                        #     for (i in 1 : (nrow(Blk3) - 1)) {if (Blk3[i, 2] >= AccInt[run] & Blk3[i + 1, 2] >= AccInt[run]) { Blk3[i, 4] <- 1}}#end i loop
                        #     HRsite <- cbind(SNPChrPosChr[Blk3[Blk3[, 4] == 1, 3], 2 : 3], Blk3[Blk3[, 4] == 1, 3])
                        #     rownames(HRsite) <- SNPChrPosChr[Blk3[Blk3[, 4] == 1, 3], 1]
                        #     HR_right <- rbind(HR_right, HRsite)
                        # }
                    } #end chr loop
                    # colnames(HR_left)[3] <- mc
                    HR_list[[mc]] <- rbind(HR_left, HRsite)

                    pos_calc <- c(1,1)
                    for (l in 1:length(HR_list)) {
                        HRsites <- HR_list[[l]]
                        HRsiteAccItp_SCs <- NULL
                        HRsiteAccRaw_SCs <- NULL

                        if (! is.null(HRsites)) {
                            rownames(dataItp) <- as.character(dataItp$Name)
                            HRPos <- NULL   #dataItp[rownames(HRsites),]

                            for (i in 1 : nrow(HRsites)) {
                                HRPos <- which(rownames(dataItp) == rownames(HRsites)[i])
                                HRsiteAcc_Ref <- as.data.frame(MCRef[c((HRPos - AccInt[run]) : (HRPos - pos_calc[l]), (HRPos + pos_calc[l]) : (HRPos + AccInt[run])), mc])
                                HRsiteItp_SC <- as.data.frame(SCItp[c((HRPos - AccInt[run]) : (HRPos - pos_calc[l]), (HRPos + pos_calc[l]) : (HRPos + AccInt[run])), mc])
                                HRsiteRaw_SC <- as.data.frame(SCRaw[c((HRPos - AccInt[run]) : (HRPos - pos_calc[l]), (HRPos + pos_calc[l]) : (HRPos + AccInt[run])), mc])
                                AccItp_SCs <- matrix(NA, ncol(HRsiteItp_SC), AccInt[run] * 2)
                                AccRaw_SCs <- matrix(NA, ncol(HRsiteItp_SC), AccInt[run] * 2)

                                for (sc in 1 : ncol(HRsiteItp_SC)) {
                                    AccItp_SCs[sc,] <- HRsiteAcc_Ref == HRsiteItp_SC[, sc] & HRsiteAcc_Ref != 0
                                    AccRaw_SCs[sc,] <- HRsiteAcc_Ref == HRsiteRaw_SC[, sc] & HRsiteAcc_Ref != 0
                                }#end sc loop
				if (!grepl("E02_Bl001",colnames(SCItp)[mc]) && !grepl("E02_Bl002",colnames(SCItp)[mc])) {
                               		HRsiteAccItp_SCs <- rbind(HRsiteAccItp_SCs, AccItp_SCs)
                                	HRsiteAccRaw_SCs <- rbind(HRsiteAccRaw_SCs, AccRaw_SCs)
				}

                            }#end i loop
                        }
                        HRsiteAccItp_SCs_bothSide <- rbind(HRsiteAccItp_SCs_bothSide,HRsiteAccItp_SCs)
                        HRsiteAccRaw_SCs_bothSide <- rbind(HRsiteAccRaw_SCs_bothSide,HRsiteAccRaw_SCs)
                    } # end HR_list
                } # end mc loop


                if (! is.null(HRsiteAccItp_SCs_bothSide)) {
                    ## calculate accuracy
                    AccItpSC <- (apply(HRsiteAccItp_SCs_bothSide, 2, sum) / nrow(HRsiteAccItp_SCs_bothSide)) * 100
                    AccRawSC <- (apply(HRsiteAccRaw_SCs_bothSide, 2, sum) / nrow(HRsiteAccItp_SCs_bothSide)) * 100

                    outfile <- paste(pgdDIR, "HRsiteAcc_", colnames(MCRef)[mc], "_run", run, ".pdf", sep = "")
                    print(outfile)
                    pdf(outfile, w = 10, h = 5, bg = "transparent", family = "Helvetica", title = "RB plot")

                    plot(c(- AccInt[run] : - 1, 1 : AccInt[run]), matrix(0, 1, AccInt[run] * 2), xlab = "SNP Position", ylab = "HR-detection accuracy (%)", pch = "*", ylim = c(0, 100), frame = FALSE, cex.main = 0.8, main = main_text[run])

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


                    # only plot left side of the volcano plot
                    outfile <- paste(pgdDIR, "HRsiteAcc_", colnames(MCRef)[mc], "_run", run, ".leftHR.pdf", sep = "")
                    print(outfile)
                    pdf(outfile, w = 10, h = 5, bg = "transparent", family = "Helvetica", title = "RB plot")

                    plot(c(- AccInt[run] : - 1, 0), matrix(0, 1, (AccInt[run]+1) ), xlab = "SNP Position", ylab = "HR-detection accuracy (%)", pch = "*", ylim = c(0, 100), frame = FALSE, cex.main = 0.8, main = main_text[run])

                    if (run == 2) {
                        points(c(- AccInt[run] : - 1, 0), AccItpSC[1:(AccInt[1]+1)], "h", col = "#ff000080")
                    } else {
                        points(c(- AccInt[run] : - 1, 0), AccItpSC[1:(AccInt[1]+1)], "h", col = "#ff000020")
                        points(c(- AccInt[run] : - 1, 0), AccRawSC[1:(AccInt[1]+1)], "h", col = "#ff000080")
                    }

                    abline(h = 99, col = "grey")
                    abline(h = 95, col = "light green")
                    #abline(v=500,lty=2,col="grey",lwd=2)
                    dev.off()

		ACCItp_all[[colnames(MCRef)[mc]]] <- AccItpSC
		
                }
            } # end eid loop
        } # end p loop
    } # end Par loop
} # end o loop


library(RColorBrewer)
#all palette available from RColorBrewer
display.brewer.all()
#we will select the first 4 colors in the Set1 palette
cols<-c(brewer.pal(n=9,name="Set1"),brewer.pal(n=8,name="Set2"),brewer.pal(n=7,name="Set3"))

outfile <- paste(pgdDIR, "HRsiteAcc_all_run", run, ".pdf", sep = "")
pdf(outfile, w = 10, h = 5, bg = "transparent", family = "Helvetica", title = "RB plot")
plot(c(- AccInt[run] : - 1, 0), matrix(0, 1, (AccInt[run]+1) ), xlab = "SNP Position", ylab = "HR-detection accuracy (%)", pch = ".", ylim = c(0, 100), frame = FALSE, cex.main = 0.8, main = main_text[run])
for (i in 1:length(ACCItp_all)) {
	acc <- ACCItp_all[[i]]
	points(c(- AccInt[run] : - 1, 0), acc[1:(AccInt[1]+1)], col = cols[i], cex=0.4, pch="*")
}
abline(h = 99, col = "grey")
abline(h = 95, col = "red")

legend("bottomleft",names(ACCItp_all))
dev.off()


