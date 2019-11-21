#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%   	      		        			          Function: inthap            	   											     %%%%%%%%%%%%%
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): raw haplotypes
# 
#
#(->): interpreted haplotypes
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testmedfilt <- function(Hap,Window){
	
Chroms <- unique(Hap$Chr)

#===================================
#                        First cycle
	HapIntGenome <- NULL	
	for(chr in Chroms){
	
	HapChr <- Hap[Hap$Chr==chr,4]
	
	#Defining chr-specific window size 
	WindowChr <- round((Window*length(Hap$Chr==chr))/length(Hap$Chr==1))
	if((WindowChr%%2)==0){WindowChr<-WindowChr+1}

	HapIntChr <- medfilt1(HapChr,WindowChr)
	
	HapIntGenome <- c(HapIntGenome,HapIntChr)
	   
	   #print(chr)  
   }#end chr loop
   
   HapIntGenome
	
}#end function
