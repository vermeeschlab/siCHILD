#!/uz/data/shortcuts/apps/perl


##############################
# Start platform adjustments #
##############################


my $rowcounts = -1; # Set the row count to -1 so that headers are ignored
my $dataDir = "/uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/PGD322_C2/sichild//current/result/tmp/";
my $platform = "Illumina";
my $fh = "/uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/PGD322_C2/sichild//current/result/tmp/PGD322_C2.adj";

############
# Illumina #
############
if ($platform eq "Illumina")
{
	# Upload the file first
	open (LOCAL, ">${dataDir}/${upload}") or error ( "Cannot make file for upload:$!"); 
	
	# Necessary for windows servers. Greater portability by specifying this
	binmode LOCAL;
	binmode $fh;

	IlluminaFH: while (<$fh>)
	{		
		# Skip commented lines
		if (/^\s*\#/)
		{
			next IlluminaFH;
		}
		if (/^\s*$/)
		{
			next IlluminaFH;
		}
		
		# R doesn't like # so get rid of it elsewhere
		s/\#//g;
		   
		# Make the header something the script will find. Substitute Chr field for Chromosome
		s/^Chr${delimiter}/Chromosome${delimiter}/g;
		s/${delimiter}Chr${delimiter}/${delimiter}Chromosome${delimiter}/g;
		s/${delimiter}Chr\n/${delimiter}Chromosome\n/g;
		

		# Adjust the name of the position column
		s/^Position${delimiter}/Physical.Position${delimiter}/g;
		s/${delimiter}Position${delimiter}/${delimiter}Physical.Position${delimiter}/g;
		s/${delimiter}Position\n/${delimiter}Physical.Position\n/g;
		
		# Get read or .GType suffix on genotype columns
		s/\.GType//g;
		
		print LOCAL $_; # Now write the file
		
		++$rowcounts; # Autoincrement the row count
    }

    close LOCAL; # Close things nicely
    
    WriteRTemplate( $codeDir,$dataDir,$upload,$chrom,$chromList,$rComparisonIndexString,$compiledDir,$rowcounts,$postscriptWidth,$postscriptHeight,$genomeBuild, $totalNumberOfComparisons,$delimiter,$makePostscript, $runmode,$segmentation, 0 );
	# $rComparisonIndexString = "c\($ind1, $ind2\)";
	# if ($segmentation eq "TRUE")




}








##################################
# Subroutine for error reporting #
##################################
sub error
{
	my ($out, $message) = @_;
	
	print $out->header("text/html");
	print ("<html>\n<head>\n<title>Error Page\n</title>\n\n<style type=\"text/css\">\ntable {margin-left: auto; margin-right: auto}\ntr {text-align: center; vertical-align: middle}\nth {text-align: center}\nbody {font-family: verdana, arial, helvetica, sans-serif}\n</style>\n</head>\n<body>\n");
	
	print ("<p>An error has occurred during file processing\n<br>Error message: ${message}\n</p>\n");
	print ("</body>\n</html>");

	exit 1;
}

####################################
# Subroutine for timestamping runs #
####################################
sub timestamp
{
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	
	$year = 1900 + $yearOffset;
	
	if ($minute < 10)
	{
		$minute = "0" . $minute;
	}
	
	if ($second < 10)
	{
		$second = "0" . $second;
	}
	
	$theTime = "$hour:$minute:$second on $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	
	return ($theTime);
}

##############################################
# Subroutine for generating custom R scripts #
##############################################
sub WriteRTemplate
{
	my ( $codeDir, $dataDir, $uploadName, $chromName, $chromListObjects, $RcomparisonVector,
	$compiledDir,$rowCounts,$psWidth,$psHeight,$Build,$ComparisonVector,$Delimiter,$Makeps, $RunMode, $BEDMode, $skipValue ) = @_;
	
	open TEMPLATE, "${codeDir}/snpduo_code_template.R" or error ( $cgi, "Cannot open snpduo_code_template: $!" );
	open R_CODE, ">${dataDir}/${uploadName}.R" or error ( $cgi, "Cannot make custom R code: $!" );
	
	while (<TEMPLATE>)
	{
		s/SKIP_HOLDER/$skipValue/g;
		s/UPLOAD_HOLDER/$uploadName/g;
		s/CHR_HOLDER/$chromName/g;
		s/CHR_LIST/$chromListObjects/g;
		s/IND_HOLDER/$RcomparisonVector/g;
		s/SCRIPT_SOURCE/$codeDir\/SNPduoFunctions\.R/g;
		s/COMPILE_DIR/$compiledDir/g;
		s/ROW_HOLDER/$rowCounts/g;
		s/WIDTH_HOLDER/$psWidth/g;
		s/HEIGHT_HOLDER/$psHeight/g;
		s/BUILD_HOLDER/$Build/g;
		s/COMPARISON_VECTOR/1\:$ComparisonVector/;
		s/SEP_HOLDER/$Delimiter/;
		s/PS_HOLDER/$Makeps/g;
		s/MODE_HOLDER/$RunMode/g;
		s/BED_HOLDER/$BEDMode/g;
		if (CAIRO eq "TRUE")
		{
			s/PNG_HOLDER/TRUE/g;
		}
		else
		{
			s/PNG_HOLDER/FALSE/g;
		}
		
		print R_CODE;
	}
	
	close TEMPLATE;
	close R_CODE;
}

########################################
# Subroutine for creating images using #
# the PerlMagick Module                #
########################################
if ( PERLMAGICK eq "TRUE" ) {
	sub PerlMagickConvertPStoPNG
	{
		
		my ($dataDir,$uploadName,$currentName,$psWidth,$psHeight,$Thumbnail) = @_;
		
		my $imageCMD = Image::Magick->new;
		my $image = $imageCMD->Read("${dataDir}/${uploadName}_${currentName}.ps");
		$image = $imageCMD->Rotate(degrees=>90);
		my $pngwidth = int($psWidth*72);
		my $pngheight = int($psHeight*72);
		$image = $imageCMD->Resize(geometry=>"${pngwidth}x${pngheight}");
		$image = $imageCMD->Write("${dataDir}/${uploadName}_${currentName}.png");
		
		if( $Thumbnail eq "TRUE" )
		{
			$image = $imageCMD->Resize(geometry=>"124x96");
			$image = $imageCMD->Write("${dataDir}/${uploadName}thumb_${currentName}.png");
		}
	}
}

########################################
# Subroutine for creating images using #
# the command-line version of          #
# Image Magick convert                 #
########################################
sub CLForkConvertPStoPNG
{
	my ($dataDir,$uploadName,$currentName,$psWidth,$psHeight,$Thumbnail) = @_;
	
	my $pngwidth = int($psWidth*72);
	my $pngheight = int($psHeight*72);
	
	system( "convert -rotate \"90\" +antialias -resize ${pngwidth}x${pngheight} ${dataDir}/${uploadName}_${currentName}.ps ${dataDir}/${uploadName}_${currentName}.png" );
	
	if ( $Thumbnail eq "TRUE" )
	{
		system( "convert -rotate \"90\" +antialias -resize 124x96 ${dataDir}/${uploadName}_${currentName}.ps ${dataDir}/${uploadName}thumb_${currentName}.png" );
	}
}

###########################################
# Subroutine for checking directory sizes #
###########################################
sub DirectoryCheck
{
	my ($uploadDir, $outputDir) = @_;
	
	if ( $uploadDir eq '/' || $uploadDir eq "/home" || $uploadDir eq "/home/" ) {
		error( "WHOA! uploadDir is a BAD choice! You could frag your whole drive." );
	}
	
	if ( $outputDir eq "/" || $outputDir eq "/home" || $uploadDir eq "/home/" ) {
		error( "WHOA! outputDir is a BAD choice! You could frag important data!!!" );
	}
	
	my $uploadDirSize = `du -c $uploadDir | grep total`;
	$uploadDirSize =~ s/^\s*(\d+)\s*total.*$/$1/e;
	my $outputDirSize = `du -c $outputDir | grep total`;
	$outputDirSize =~ s/^\s*(\d+)\s*total.*$/$1/e;
	
	if ($uploadDirSize > DIR_MAX || $outputDirSize > OUTPUT_MAX)
	{
		DirectoryClean ($uploadDir,$outputDir);
	}
}

############################################
# Subroutine for cleaning full directories #
############################################
sub DirectoryClean
{
	my ($uploadDirectory, $outputDirectory) = @_;
	
	if ( $uploadDir eq '/' || $uploadDir eq "/home" || $uploadDir eq "/home/" ) {
		error( "WHOA! uploadDir is a BAD choice! You could frag your whole drive." );
	}
	
	if ( $outputDir eq "/" || $outputDir eq "/home" || $uploadDir eq "/home/" ) {
		error( "WHOA! outputDir is a BAD choice! You could frag important data!!!" );
	}
	
	my $round = 1;
	my $dirClean = 0;
	
	opendir (DIR, $uploadDirectory);
	my @uploadFiles = readdir(DIR);
	shift (@uploadFiles); # remove "."
	shift (@uploadFiles); # remove ".."
	closedir (DIR);
	
	foreach $file (@uploadFiles)
	{ 
		$file =~ s/^([\w\s\-\_\.]+)$/$1/ || die "Can't match file $file for cleaning: $!\n";
		$file = $uploadDirectory . "/" . $file;
	}
	
	opendir (DIR, $outputDirectory);
	my @outputFiles = readdir(DIR);
	shift (@outputFiles); # remove "."
	shift (@outputFiles); # remove ".."
	closedir (DIR);
	
	foreach $file (@outputFiles)
	{
		$file =~ s/^([\w\s\-\_\.]+)$/$1/ || die "Can't match file $file for cleaning: $!\n";
		$file = $outputDirectory . "/" . $file;
	}
	
	# Combine arrays
	my @allFiles = (@uploadFiles, @outputFiles);
		
	# Loop through cleaning files, multiple times if necessary
	until ($dirClean == 1 || $round > 7)
	{
		# Get each file size
		foreach $file (@allFiles)
		{
			next if ($file eq " ");
			
			my $age = -M $file || error($cgi, "In attempting to clean upload directory could not stat $file");
			
			if ($round == 1 && $age >= 7)
			{
				unlink( $file );
				$file = " ";
			}
			
			elsif ($round == 2 && $age >= 3)
			{
				unlink( $file );
				$file = " ";
			}
			
			elsif ($round == 3 && $age >= 1)
			{
				unlink( $file );
				$file = " ";
			}
			
			elsif ($round == 4 && $age >= 0.50)
			{
				unlink( $file );
				$file = " ";
			}
			
			elsif ($round == 5 && $age >= 0.25)
			{
				unlink( $file );
				$file = " ";
			}
			
			elsif ($round == 6 && $age >= 0.04)
			{
				unlink( $file );
				$file = " ";
			}
			
			else
			{
				unlink( $file );
				$file = " ";
			}
		}
		
				
		my $uploadDirSize = `du -c $uploadDirectory | grep total`;
		$uploadDirSize =~ s/^\s*(\d+)\s*.*$/$1/e or die "Can't understand upload directory size $uploadDirSize: $!\n";
		
		my $outputDirSize = `du -c $outputDirectory | grep total`;
		$outputDirSize =~ s/^\s*(\d+)\s*.*$/$1/e or die "Can't understand output directory size $outputDirSize: $!\n";
		
		if ($uploadDirSize < DIR_MAX && $outputDirSize < DIR_MAX)
		{
			$dirClean = 1;
		}
		
		++$round;
	}
	
	if ($dirClean == 0)
	{
		error($cgi, "Upload or output directory full and script was unable to clean it. Please contact the webmaster.")
		;
	}
}

#####################################################
# Subroutine for correctly sorting chromosome order #
#####################################################
sub SortChromosomeList
{
	#############################################
	# Sort the chromosomes so order makes sense #
	#############################################
	my $mitoChar = "NULL";
	
	my @chromsToSort = @_;
	
	foreach $sorting (@chromsToSort)
	{
		$sorting =~ s/XY/24/g;
		$sorting =~ s/X/23/g;
		$sorting =~ s/Y/25/g;
		
		if ( $sorting =~ /MITO/g ) {
			$sorting =~ s/MITO/26/g;
			$mitoChar = "MITO";
		}
		
		elsif ( $sorting =~ /Mito/g ) {
			$sorting =~ s/Mito/26/g;
			$mitoChar = "Mito";
		}
		
		elsif ( $sorting =~ /MT/g ) {
			$sorting =~ s/MT/26/g;
			$mitoChar = "MT";
		}
		
		elsif ( $sorting =~ /M/g ) {
			$sorting =~ s/M/26/g;
			$mitoChar = "M";
		}
	}
	
	@chromsToSort = sort {$a <=> $b} @chromsToSort;
	
	foreach $sorting (@chromsToSort)
	{
		$sorting =~ s/24/XY/g;
		$sorting =~ s/23/X/g;
		$sorting =~ s/25/Y/g;
		$sorting =~ s/26/${mitoChar}/g;
	}
	
	return( @chromsToSort );
}
