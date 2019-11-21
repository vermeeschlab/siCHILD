par(mai = c(0.1, 0.1, 0.1, 0.1)); # margin: bottom, left, top, right

data_dir <- "/uz/data/avalok/symbiosys/gcpi_r_kul_joris_vermeesch/jding/PGD/figure/"
outfile <- paste0(data_dir,"PGD081.tiff")
print(outfile)
tiff(outfile,width = 10, height = 10,units = "in",res = 300)
plot(c(1,100), c(1,50), xlab = "", ylab="", axes = FALSE, xpd=FALSE, col="white")

#
mid=9.5
ymid = 50

#
rect(3.5, 65-ymid, 4.5, 70-ymid, col="dodgerblue4", border=NA)
rect(5, 65-ymid, 6, 70-ymid, col="darksalmon",  border=NA)

#
rect(13, 65-ymid, 14, 70-ymid, col="dodgerblue4", border=NA)
rect(14.5, 65-ymid, 15.5, 70-ymid, col="darksalmon", border=NA)


# horizontal
segments(mid-2.5, 74-ymid, mid+2.5, 74, lwd=1.5)
points(6.25, 74-ymid, pch=15, col="black", cex=2.5) #GF
points(12.75, 74-ymid, pch=21, col="black", cex=2.5)	#GM

# vertical
segments(mid,     64-ymid, mid,     74-ymid, lwd=1.5)
# horizontal
segments(mid-7.5, 64-ymid, mid+12.5, 64-ymid, lwd=1.5)
# vertical
segments(mid-2.5, 64-ymid, mid-2.5, 58.5-ymid, lwd=1.5)
points(2, 57-ymid, pch=21, col="black", cex=2.5)	#Aunt1
segments(mid+2.5, 64-ymid, mid+2.5, 58.5-ymid, lwd=1.5)
points(7, 57-ymid, pch=21, col="black", cex=2.5)	#Aunt2
segments(mid-7.5, 64-ymid, mid-7.5, 58.5-ymid, lwd=1.5)
points(12, 57-ymid, pch=15, col="black", cex=2.5)	#Uncle
segments(mid+7.5, 64-ymid, mid+7.5, 58.5-ymid, lwd=1.5)
points(17, 57-ymid, pch=19, col="black", cex=2.5)	#Mother
#
segments(mid+12.5, 64-ymid, mid+12.5, 58.5-ymid, lwd=1.5)
points(22, 57-ymid, pch=22, col="black", cex=2.5)	#Father
dev.off()



