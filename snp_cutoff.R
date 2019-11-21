#--Define data file directory:
dir <- "http://www.sr.bham.ac.uk/~ajrs/papers/sanderson06"

#--Read in table of data (from Sanderson et al. 2006):
# This refers to a sample of 20 clusters of galaxies with Chandra X-ray data
CC <- read.table(paste(dir, "mean_Tprofile-CC.txt", sep="/"), header=TRUE)
nCC <- read.table(paste(dir, "mean_Tprofile-nCC.txt", sep="/"), header=TRUE)

#--Load extra library:
## if not already installed, then run:
# install.packages("ggplot2")
require(ggplot2)

#--Combine datasets into a single data frame:
CC$type <- "Cool core"
nCC$type <- "Non-cool core"
A <- rbind(CC, nCC)

#--Define axis labels:
xlabel <- as.expression(expression( paste("Radius (", R[500], ")") ))
ylabel <- "Scaled Temperature"

p <- ggplot(data=A, aes(x=r.r500, y=sckT, ymin=sckT.lo, ymax=sckT.up, fill=type, linetype=type)) + 
 geom_line() + 
 geom_ribbon(alpha=0.5) + 
 scale_x_log10() + 
 scale_y_log10() + 
 xlab(xlabel) + 
 ylab(ylabel)

ggsave(p, file="CC-vs_nCC_kT_prof.pdf", width=8, height=4.5)
