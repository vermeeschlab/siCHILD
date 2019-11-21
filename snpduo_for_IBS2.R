
#############################
#	upgrade from SNPduo #
#       Author: Jia Ding    #
#############################

library(plyr)



snpduo_for_IBS2 = function( upload, chrom, col1,col2,genomebuild )
{



comparisonVector = 1	#?
comparisonVectorCounter = 1	
comparisonColumns = c(col1,col2)


# Counter

SKIP = 0

# input = read.delim( upload, colClasses="character", comment.char="", nrow=ROW_HOLDER, sep=SEP, skip=SKIP )
input = read.delim( upload, colClasses="character", comment.char="", sep="\t", skip=0 )
input <- input[order(input[,2],input[,1]),]



if ( dim( input )[2] == 1 )
{
	stop( paste( "One column found in uploaded data. Check your file format to be sure this is correct" ) )
}

names( input )[ which( names( input ) == "Chr" ) ] = "Chromosome"
names( input )[ which( names( input ) == "Position" ) ] = "Physical.Position"



# chromList = CHR_LIST
chromList = sort(unique(input[,2]))

pswidth = 11

psheight = 8.5

#genomebuild = "vGRCh37"

cytoband = load_chromosome_features(compiled, genomebuild)



# don't produce .ps & .png files
makepostscript = FALSE

makePNG = FALSE

MODE = "Normal"

BED = TRUE

# Load the library now and do it only once. No need to access it repeatedly.

dyn.load( file.path( compiled, paste( "SNPduoCCodes", .Platform$dynlib.ext, sep="" ) ) )




if ( MODE != "Tabulate" )
{
	write.table( comparisonVector, file=paste( upload, "_Comparisons.txt", sep="" ), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t" )
	
	for ( person1index in 1:( length( comparisonColumns ) - 1 ) )
	{
		ind1 = comparisonColumns[person1index]
		
		for ( person2index in (person1index + 1):length( comparisonColumns ) )
		{
			ind2 = comparisonColumns[person2index]

			if( chrom == "Genome" )
			{	
				whole_genome_plot(cytoband, input, ind1, ind2, savename=upload, pswidth, psheight, comparison=comparisonVector[comparisonVectorCounter], doPostscript=makepostscript, makeBED = BED, doPNG=makePNG, chr.offset=load_chromosome_position_offsets( compiled, genomebuild ) )
			} else if ( chrom == "GenomeByChromosome" )
			{	
				genome_by_chromosome(cytoband, input,ind1, ind2, savename=upload, pswidth, psheight, comparison=comparisonVector[comparisonVectorCounter], doPostscript=makepostscript,chromlist=chromList, makeBED = BED, doPNG=makePNG)
			} else
			{	print(names(input))
				snpduo_single_chromosome( cytoband,input,chrom, ind1, ind2, savename=upload, pswidth, psheight, comparison=comparisonVector[comparisonVectorCounter], doPostscript=makepostscript, makeBED = BED, doPNG=makePNG)

			}
			
			comparisonVectorCounter = comparisonVectorCounter + 1
		}
	}
} else
{
	tabulate_ibs( input, chromList, unique( comparisonColumns ), upload )
}

# Unload at the end to be sure everything closes nicely
dyn.unload( file.path( compiled, paste( "SNPduoCCodes", .Platform$dynlib.ext, sep="" ) ) )

}


