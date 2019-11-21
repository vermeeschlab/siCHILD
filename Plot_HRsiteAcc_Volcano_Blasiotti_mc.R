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
ct <- c("SC", "FC")
AccInt <- c(250, 80)
snps <- c("all", "informative") # run=1 => accuracy with all snps; run=2 => accuracy with informative snps.


opt <- read.table(paste0(dataPath, opt_fn), header = T, stringsAsFactors = F)

ACCItp_all_df <- c()

for (e in 1 : length(eid)) {    # E01 or E02
    for (c in 1 : length(ct)) { # SC or FC compare with MC
        for (p in 1 : nrow(Pars)) { # Pat or Mat
            ACCItp_all <- list()

            embryoID <- paste(eid[e], ct[c], Pars[p], sep = "_")

            for (o in 1 : nrow(opt)) {  # all conditions: GP GM GF ect...
                   print(paste0("----", embryoID, "----"))

	        pgdDIR <- paste0(dataPath, opt[o,]$Dir, opt[o,]$subDir)
                # print(paste0("----", pgdDIR, "----"))
                File1 <- paste0(pgdDIR, "/", opt[o,]$Dir, "_Itp2.hap")
                File2 <- paste0(pgdDIR, "/", opt[o,]$Dir, "_Raw.hap")
                File3 <- paste0(pgdDIR, "/", opt[o,]$Dir, "_0.75.gtp")

                dataItp <- read.table(File1 , sep = "\t", header = T, stringsAsFactors = F)
                dataRaw <- read.table(File2 , sep = "\t", header = T, stringsAsFactors = F)
                gtp <- read.table(File3, sep = "\t", header = T, stringsAsFactors = F)

                ibs1_fn <- list.files(path = paste0(pgdDIR, "../"), pattern = "olaps1")
                # print(ibs1_fn)
                if (! identical(ibs1_fn, character(0))) {
                    ibs1 <- read.table(paste0(pgdDIR, "../", ibs1_fn), sep = "\t", header = T, stringsAsFactors = F)
                    dataItp <- dataItp[(dataItp)$Name %in% ibs1[, 1],]
                    dataRaw <- dataRaw[(dataRaw)$Name %in% ibs1[, 1],]
                    gtp <- gtp[(gtp)$Name %in% ibs1[, 1],]
                }

                SNPChrPos <- dataItp[, 1 : 3]
                rownames(dataItp) <- SNPChrPos[, 1]
                rownames(dataRaw) <- SNPChrPos[, 1]

                dataItp_org <- dataItp
                dataRaw_org <- dataRaw
                gtp_org <- gtp

                par <- Pars[p, 1]
                dataItp <- dataItp_org
                dataRaw <- dataRaw_org
                gtp <- gtp_org
                # print(paste0("-------------------",par,"------------------"))

                dataItp <- dataItp[, c(1 : 3, grep(par, colnames(dataItp)))]
                dataRaw <- dataRaw[, c(1 : 3, grep(par, colnames(dataRaw)))]
                # print(paste0("=======",pgdDIR,"=======",nrow(dataItp)))

                for (run in 1 : 1) {    # run over all the SNPs or only informative SNPs
                    if (run == 2) {
                        rn <- gtp[gtp[, grep(Pars[p, 2], colnames(gtp))] == "AB" &
                            gtp[, grep(Pars[p, 3], colnames(gtp))] != "AB" &
                            gtp[, grep(Pars[p, 3], colnames(gtp))] != "NC", 1]
                        dataItp <- dataItp[dataItp$Name %in% rn,]
                        dataRaw <- dataRaw[dataRaw$Name %in% rn,]
                    }

                    x <- c(grep(eid[e], colnames(dataItp)), grep("_Bl01", colnames(dataItp)))
                    fc_cn <- x[duplicated(x)]
                    fcRef <- as.data.frame(dataItp[, fc_cn])
                    x <- c(grep(eid[e], colnames(dataItp)), grep("_Bl00", colnames(dataItp)))
                    sc_cn <- x[duplicated(x)]
                    x <- grep(paste0(eid[e], "_Bl099"), colnames(dataItp))
                    mc_cn <- x
                    mcRef <- as.data.frame(dataItp[, mc_cn])
                    colnames(mcRef) <- colnames(dataItp)[mc_cn]
                    mc <- ncol(mcRef)


                    main_text <- paste0("checking accuracy of ", embryoID , " with ", snps[run], " SNPs. ")
                    # print(main_text)

                    HRsiteAccItp_SCs_bothSide <- NULL
                    HRsiteAccRaw_SCs_bothSide <- NULL
                    HR_list <- NULL
                    sc_fc <- as.data.frame(matrix(0, 2 * length(sc_cn), nrow = 2, ncol = length(sc_cn)))
                    rownames(sc_fc) <- ct
                    sc_fc[1,] <- sc_cn
                    sc_fc[2,] <- fc_cn

                    #Determining HR-sites in each of mc reference
                    HR_left <- NULL
                    HR_right <- NULL

                    for (chr in Chroms) {
                        dataChr <- as.data.frame(mcRef[dataItp$Chr == chr,])
                        SNPChrPosChr <- dataItp[dataItp$Chr == chr, 1 : 3]
                        Blk <- rle(dataChr[, mc])
                        if (! ((length(Blk$values) == 1) && Blk$value == 0)) {

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
                                for (i in 1 : (nrow(Blk2) - 1)) {if (Blk2[i, 2] >= AccInt[run] & Blk2[i + 1, 2] >= AccInt[run]) { Blk2[i, 4] <- 1 }}#end i loop
                                HRsite <- cbind(SNPChrPosChr[Blk2[Blk2[, 4] == 1, 3], 2 : 3], Blk2[Blk2[, 4] == 1, 3], rep(mc, length(Blk2[Blk2[, 4] == 1, 3])))
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
                        }
                    } #end chr loop
                    # colnames(HR_left)[3] <- fc
                    HR_list[[mc]] <- rbind(HR_left, HR_right)


                    if (is.null(HR_list[[mc]])) {
                        print(paste0("no RBs in this multicell", embryoID, "phasing with the same cell"))
                    } else {

                        pos_calc <- c(1, 1)
                        for (l in 1 : length(HR_list)) {
                            HRsites <- HR_list[[l]]
                            HRsiteAccItp_SCs <- NULL
                            HRsiteAccRaw_SCs <- NULL

                            rownames(dataItp) <- as.character(dataItp$Name)
                            HRPos <- NULL   #dataItp[rownames(HRsites),]

                            for (i in 1 : nrow(HRsites)) {
                                SCItp <- dataItp[, as.numeric(sc_fc[c,])]
                                SCRaw <- dataRaw[, as.numeric(sc_fc[c,])]

                                HRPos <- which(rownames(dataItp) == rownames(HRsites)[i])

                                HRsiteAcc_Ref <- as.data.frame(mcRef[c((HRPos - AccInt[run]) : (HRPos - pos_calc[l]), (HRPos + pos_calc[l]) : (HRPos + AccInt[run])),])
                                HRsiteItp_SC  <- as.data.frame(SCItp[c((HRPos - AccInt[run]) : (HRPos - pos_calc[l]), (HRPos + pos_calc[l]) : (HRPos + AccInt[run])),])
                                HRsiteRaw_SC  <- as.data.frame(SCRaw[c((HRPos - AccInt[run]) : (HRPos - pos_calc[l]), (HRPos + pos_calc[l]) : (HRPos + AccInt[run])),])
                                AccItp_SCs <- matrix(NA, ncol(HRsiteItp_SC), AccInt[run] * 2)
                                AccRaw_SCs <- matrix(NA, ncol(HRsiteItp_SC), AccInt[run] * 2)

                                for (sc in 1 : ncol(HRsiteItp_SC)) {
                                    if (! grepl("E02_Bl001", colnames(HRsiteItp_SC)[sc]) && ! grepl("E02_Bl002", colnames(HRsiteItp_SC)[sc])) {  # because E02_Bl001 & E02_Bl002 are two bad single cells, they are removed.
                                        AccItp_SCs[sc,] <- HRsiteAcc_Ref == HRsiteItp_SC[, sc] & HRsiteAcc_Ref != 0
                                        AccRaw_SCs[sc,] <- HRsiteAcc_Ref == HRsiteRaw_SC[, sc] & HRsiteAcc_Ref != 0
                                    }
                                }#end sc loop

                                AccItp_SCs <- AccItp_SCs[complete.cases(AccItp_SCs),]
                                AccRaw_SCs <- AccRaw_SCs[complete.cases(AccRaw_SCs),]
                                HRsiteAccItp_SCs <- rbind(HRsiteAccItp_SCs, AccItp_SCs)
                                HRsiteAccRaw_SCs <- rbind(HRsiteAccRaw_SCs, AccRaw_SCs)
                            }

                            HRsiteAccItp_SCs_bothSide <- rbind(HRsiteAccItp_SCs_bothSide, HRsiteAccItp_SCs)
                            HRsiteAccRaw_SCs_bothSide <- rbind(HRsiteAccRaw_SCs_bothSide, HRsiteAccRaw_SCs)
                        } # end j


                        if (! is.null(HRsiteAccItp_SCs_bothSide)) {
                            ## calculate accuracy
                            AccItpSC <- (apply(HRsiteAccItp_SCs_bothSide, 2, sum) / nrow(HRsiteAccItp_SCs_bothSide)) * 100
                            AccRawSC <- (apply(HRsiteAccRaw_SCs_bothSide, 2, sum) / nrow(HRsiteAccItp_SCs_bothSide)) * 100

                            # plot left & right side of the volcano plot
                            outfile <- paste(pgdDIR, "HRsiteAcc_", colnames(mcRef), "_run", run, ".pdf", sep = "")
                            # print(outfile)

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
                            outfile <- paste(pgdDIR, "HRsiteAcc_", colnames(mcRef), "_run", run, ".leftHR.pdf", sep = "")
                            # print(outfile)
                            pdf(outfile, w = 10, h = 5, bg = "transparent", family = "Helvetica", title = "RB plot")

                            plot(c(- AccInt[run] : - 1, 0), matrix(0, 1, (AccInt[run] + 1)), xlab = "SNP Position", ylab = "HR-detection accuracy (%)", pch = "*", ylim = c(0, 100), frame = FALSE, cex.main = 0.8, main = main_text[run])

                            if (run == 2) {
                                points(c(- AccInt[run] : - 1, 0), AccItpSC[1 : (AccInt[1] + 1)], "h", col = "#ff000080")
                            } else {
                                points(c(- AccInt[run] : - 1, 0), AccItpSC[1 : (AccInt[1] + 1)], "h", col = "#ff000020")
                                points(c(- AccInt[run] : - 1, 0), AccRawSC[1 : (AccInt[1] + 1)], "h", col = "#ff000080")
                            }

                            abline(h = 99, col = "grey")
                            abline(h = 95, col = "light green")
                            #abline(v=500,lty=2,col="grey",lwd=2)
                            dev.off()

                            ACCItp_all[[gsub("Bl099", rownames(sc_fc)[c], colnames(mcRef))]] <- AccItpSC
			}
                    } # end l loop
                } # end run loop
            } # end o loop



            #########################################
            ## plot embryo E01/E02, pat/mat, sc/fc ##
            #########################################
            library(RColorBrewer)
            #all palette available from RColorBrewer
            display.brewer.all()
            #we will select the first 4 colors in the Set1 palette
            cols <- c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Set2"), brewer.pal(n = 7, name = "Set3"))

            outfile <- paste(pgdDIR, "../../PGD_Blasiotti_HRsite_plot/HRsiteAcc_all_", embryoID, "_run", run, ".pdf", sep = "")
	    
	    print(outfile)
            pdf(outfile, w = 10, h = 5, bg = "transparent", family = "Helvetica", title = "RB plot")
            plot(c(- AccInt[run] : - 1, 0), matrix(0, 1, (AccInt[run] + 1)), xlab = "SNP Position", ylab = "HR-detection accuracy (%)", pch = ".", ylim = c(0, 100), frame = FALSE, cex.main = 0.8)
            for (i in 1 : length(ACCItp_all)) {
                acc <- ACCItp_all[[i]]
                points(c(- AccInt[run] : - 1, 0), acc[1 : (AccInt[1] + 1)], col = cols[i], cex = 0.4, pch = "*")
            }
            abline(h = 99, col = "grey")
            abline(h = 95, col = "red")

            legend(x = - 250, y = 30, names(ACCItp_all), cex = 0.3, border = "white", text.col = cols)
            dev.off()
	    df <- data.frame(matrix(unlist(ACCItp_all), nrow=length(ACCItp_all), byrow=T))
	    rownames(df) <- names(ACCItp_all)
            outfile <- paste(pgdDIR, "../../PGD_Blasiotti_HRsite_plot/HRsiteAcc_all_", embryoID, "_run", run, ".csv", sep = "")
	    print(outfile)
	    write.table(df, outfile, append= T, sep='\t', row.names=T, col.names=T )
ACCItp_all_df <- rbind(ACCItp_all_df,df)
        }   # end Pars
    } # end ct loop
} # end eid loop

outfile <- paste(pgdDIR, "../../PGD_Blasiotti_HRsite_plot/HRsiteAcc_all",".csv", sep = "")
write.table(ACCItp_all_df, outfile, append= T, sep='\t', row.names=T, col.names=T )

rn_aff_sc <- intersect(intersect(grep("Mat",rownames(ACCItp_all_df)),grep("Aff",rownames(ACCItp_all_df))),grep("SC",rownames(ACCItp_all_df)))
rn_aff_fc <- intersect(intersect(grep("Mat",rownames(ACCItp_all_df)),grep("Aff",rownames(ACCItp_all_df))),grep("FC",rownames(ACCItp_all_df)))

rn_aunt_sc <- intersect(c(grep("Aunt",rownames(ACCItp_all_df)),grep("Uncle",rownames(ACCItp_all_df))),grep("SC",rownames(ACCItp_all_df)))
rn_aunt_fc <- intersect(c(grep("Aunt",rownames(ACCItp_all_df)),grep("Uncle",rownames(ACCItp_all_df))),grep("FC",rownames(ACCItp_all_df)))

df_aunt_sc<- as.matrix(ACCItp_all_df[rn_aunt_sc,])
df_aunt_fc<- as.matrix(ACCItp_all_df[rn_aunt_fc,])
df_aff_fc <- as.matrix(ACCItp_all_df[rn_aff_fc,])
df_aff_sc <- as.matrix(ACCItp_all_df[rn_aff_sc,])

library(Rfast)
df_aunt_sc_mean <- colMeans(x=df_aunt_sc)
df_aunt_sc_min <- colMins(x=df_aunt_sc)
df_aunt_sc_max <- colMaxs(x=df_aunt_sc)

df_aunt_fc_mean <- colMeans(x=df_aunt_fc)
df_aunt_fc_min <- colMins(x=df_aunt_fc)
df_aunt_fc_max <- colMaxs(x=df_aunt_fc)

df_aff_sc_mean <- colMeans(x=df_aff_sc)
df_aff_sc_min <- colMins(x=df_aff_sc)
df_aff_sc_max <- colMaxs(x=df_aff_sc)

df_aff_fc_mean <- colMeans(x=df_aff_fc)
df_aff_fc_min <- colMins(x=df_aff_fc)
df_aff_fc_max <- colMaxs(x=df_aff_fc)



df_all <- rbind(
cbind(as.numeric(apply(t(df_aff_fc[,100:250]), 1, min, na.rm=FALSE)),as.numeric(apply(t(df_aff_fc[,100:250]), 1, max, na.rm=FALSE)),as.numeric(apply(t(df_aff_fc[,100:250]), 1, median, na.rm=FALSE)),-150:0,rep("parental_offspring_FC",151)),
cbind(as.numeric(apply(t(df_aff_sc[,100:250]), 1, min, na.rm=FALSE)),as.numeric(apply(t(df_aff_sc[,100:250]), 1, max, na.rm=FALSE)),as.numeric(apply(t(df_aff_sc[,100:250]), 1, median, na.rm=FALSE)),-150:0,rep("parental_offspring_SC",151)),
cbind(as.numeric(apply(t(df_aunt_fc[,100:250]), 1, min, na.rm=FALSE)),as.numeric(apply(t(df_aunt_fc[,100:250]), 1, max, na.rm=FALSE)),as.numeric(apply(t(df_aunt_fc[,100:250]), 1, median, na.rm=FALSE)),-150:0,rep("parental_sibling_FC",151)),
cbind(as.numeric(apply(t(df_aunt_sc[,100:250]), 1, min, na.rm=FALSE)),as.numeric(apply(t(df_aunt_sc[,100:250]), 1, max, na.rm=FALSE)),as.numeric(apply(t(df_aunt_sc[,100:250]), 1, median, na.rm=FALSE)),-150:0,rep("parental_sibling_SC",151))
)
colnames(df_all) <- c("s.min","s.max","s.median","range","type")
df_all <- as.data.frame(df_all)
df_all[,1:4] <- as.data.frame(apply(df_all, 2, as.numeric)[,1:4])


#--Define axis labels:
ylabel <- "HR-detection accuracy by percentage"
xlabel <- "SNP positions"
out_pdf <- paste(pgdDIR, "../../PGD_Blasiotti_HRsite_plot/HRsiteAcc_all",".pdf", sep = "")
out_tiff <- paste(pgdDIR, "../../PGD_Blasiotti_HRsite_plot/HRsiteAcc_all",".tiff", sep = "")	# when x11 is available

require(ggplot2)
df_all <- as.data.frame(df_all)

p <- ggplot(data=df_all, aes(x=range, y=s.median, ymin=s.min, ymax=s.max, fill=type, linetype=type, shape=type)) +
#	geom_hline(yintercept=99) +
	geom_point(aes(shape = factor(type)), size = 0.5, color = c("black")) +
	geom_ribbon(alpha=0.5) +
	xlab(xlabel) +
	ylab(ylabel)
	ggsave(p, file=out_pdf, width=8, height=4.5
#	ggsave(plot=image, file=out_tiff, width=8, height=4.5, dpi=300
)


