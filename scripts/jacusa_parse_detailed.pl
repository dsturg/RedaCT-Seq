#!/usr/bin/perl

# PROGRAM  : jacusa_parse_detailed.pl
# AUTHOR   : Dave Sturgill

# Usage: perl jacusa_parse_script.pl samplenumber labelstring jacusafile

# Note: This script was written for output from version JACUSA_v2.0.0-S-BETA-1

# eg perl jacusa_parse_script.pl FR-SECONDSTRAND_jacusa_pileup_WT_treated_reps.txt 3 "S1,S2,S3" > FR-SECONDSTRAND_jacusa_pileup_WT_treated_reps_parsed.txt
# 


#use strict;
	
##########################
# Read command line args
##########################

my $filename = $ARGV[0];
chomp $filename;

my $numsamples = $ARGV[1];
chomp $numsamples;

my $labelstring = $ARGV[2];
chomp $labelstring;

#####################################
# Usage / Notes
#####################################



######################
#@@@@@@@@@@@@@@@@@@@@@
######################

# parse label string

my @labels ;
my @labels = split(',', $labelstring);

my $datafile;
my @columns ;
my $prestring ;
my $refbase ;


my @depth ;
my @cov_A ;
my @cov_T ;
my @cov_C ;
my @cov_G ;


my @samplecov ;
my $refcov ;
my $mmcov ;
my $mmbase ;
my $mmflag ;

my $linecount = 0 ;
my $strandbase ;
my $mmperc ;
my $totalcounts ;
my $samplestring ;


print "chr\tstart\tend\tstrand\trefbase\tstrand_refbase\ttotalsum";
foreach (@labels)
{
	#print $_, "\n";
	$samplestring = "$_";
	#print $samplestring, "\n";
	print "\t",$samplestring,"_Depth\t",$samplestring,"_Ref\t",$samplestring,"_MMflag\t",$samplestring,"_MMbase\t",$samplestring,"_MMcounts\t",$samplestring,"_MMperc\t",$samplestring,"_MMconv\t",$samplestring,"_basestring\t",$samplestring,"_CovA\t",$samplestring,"_CovC\t",$samplestring,"_CovG\t",$samplestring,"_CovT" ;
}
print "\n";


open(DATAFILE, "$filename");
while ($datafile = <DATAFILE>) {
	++$linecount;
	if ($linecount > 2) {

		@columns = split(' ', $datafile);
	
		
		$refbase = uc($columns[(7 + $numsamples)]);

		
		
		# Flip base if on - strand
		## NOT NEEDED in Jacusa 2.0.1
		if ($columns[5] =~ /\-/) {
			if ($refbase =~ /A/) {$strandbase = "T"}
			if ($refbase =~ /C/) {$strandbase = "G"}
			if ($refbase =~ /G/) {$strandbase = "C"}
			if ($refbase =~ /T/) {$strandbase = "A"}
			
		} else {
			$strandbase = $refbase ;
		}

		$prestring = "$columns[0]\t$columns[1]\t$columns[2]\t$columns[5]\t$refbase\t$strandbase\t$columns[4]\t";
		#print $prestring, "\t",
		# For each ACTG, get the get the coverage for each substitution
	
		# Get the coverages
		@depth = ();
		@cov_A = () ;
		@cov_T = () ;
		@cov_C = () ;
		@cov_G = () ;
		
		my $CovA ;
		my $CovC ;
		my $CovG ;
		my $CovT ;
		
		print $prestring ;
		my $element ;
		my $basestring = "" ;
		
		#~~~~~~~~~~~~~~~~~~~~
		# Iterate over each sample
		#~~~~~~~~~~~~~~~~~~~~
	
		for ($j = 0 ; $j < $numsamples ; $j++) {
			$basestring = @columns[6+$j] ;
			@samplecov = split(',', @columns[6+$j]);
			$CovA = @samplecov[0] ;
			$CovC = @samplecov[1] ;
			$CovG = @samplecov[2] ;
			$CovT = @samplecov[3] ;
			#print "covT is $CovT\n";
			#sleep 2;
			#@samplecov = split(',', @columns[6+$j]);
		
			$totalcount = 0 ;
			#print "\n";
			foreach $element (@samplecov) {
				$totalcount += $element;
				#print $element, "\n";
			}
			#@@@@@@@@@@@@@@@@@@@@@@@
			# Get dominant mismatch
			#@@@@@@@@@@@@@@@@@@@@@@@
			if ($strandbase =~ /A/) {
				$refcov = @samplecov[0] ;
				if ((@samplecov[1] + @samplecov[2] + @samplecov[3]) < 1) {
					$mmbase = "None" ;
					$mmcov = 0 ;
					$mmflag = "NoMM" ;
				} elsif ((@samplecov[1] > @samplecov[2]) & (@samplecov[1] > @samplecov[3])) {
					$mmbase = "C" ;
					$mmcov = @samplecov[1] ;
					$mmflag = "DomMM" ;
				} elsif ((@samplecov[2] > @samplecov[1]) & (@samplecov[2] > @samplecov[3])) {
					$mmbase = "G" ;
					$mmcov = @samplecov[2] ;
					$mmflag = "DomMM" ;
				} elsif ((@samplecov[3] > @samplecov[1]) & (@samplecov[3] > @samplecov[2])) {
					$mmbase = "T" ;
					$mmcov = @samplecov[3] ;
					$mmflag = "DomMM" ;
				} else {
					$mmbase = "None" ;
					$mmcov = 0 ;
					$mmflag = "NoDomMM" ;
				}
			$mmperc = sprintf("%.3f",($mmcov / $totalcount)*100) ;
			print "$totalcount\t$refcov\t$mmflag\t$mmbase\t$mmcov\t$mmperc\t" ;
			} 		
			if ($strandbase =~ /C/) {
				$refcov = @samplecov[1] ;
				if ((@samplecov[0] + @samplecov[2] + @samplecov[3]) < 1) {
					$mmbase = "None" ;
					$mmcov = 0 ;
					$mmflag = "NoMM" ;
				} elsif ((@samplecov[0] > @samplecov[2]) & (@samplecov[0] > @samplecov[3])) {
					$mmbase = "A" ;
					$mmcov = @samplecov[0] ;
					$mmflag = "DomMM" ;
				} elsif ((@samplecov[2] > @samplecov[0]) & (@samplecov[2] > @samplecov[3])) {
					$mmbase = "G" ;
					$mmcov = @samplecov[2] ;
					$mmflag = "DomMM" ;
				} elsif ((@samplecov[3] > @samplecov[0]) & (@samplecov[3] > @samplecov[2])) {
					$mmbase = "T" ;
					$mmcov = @samplecov[3] ;
					$mmflag = "DomMM" ;
				} else {
					$mmbase = "None" ;
					$mmcov = 0 ;
					$mmflag = "NoDomMM" ;
				}
			$mmperc = sprintf("%.3f",($mmcov / $totalcount)*100) ;
			print "$totalcount\t$refcov\t$mmflag\t$mmbase\t$mmcov\t$mmperc\t" ;
			} 
			if ($strandbase =~ /G/) {
				$refcov = @samplecov[2] ;
				if ((@samplecov[0] + @samplecov[1] + @samplecov[3]) < 1) {
					$mmbase = "None" ;
					$mmcov = 0 ;
					$mmflag = "NoMM" ;
				} elsif ((@samplecov[0] > @samplecov[1]) & (@samplecov[0] > @samplecov[3])) {
					$mmbase = "A" ;
					$mmcov = @samplecov[0] ;
					$mmflag = "DomMM" ;
				} elsif ((@samplecov[1] > @samplecov[0]) & (@samplecov[1] > @samplecov[3])) {
					$mmbase = "C" ;
					$mmcov = @samplecov[1] ;
					$mmflag = "DomMM" ;
				} elsif ((@samplecov[3] > @samplecov[0]) & (@samplecov[3] > @samplecov[1])) {
					$mmbase = "T" ;
					$mmcov = @samplecov[3] ;
					$mmflag = "DomMM" ;
				} else {
					$mmbase = "None" ;
					$mmcov = 0 ;
					$mmflag = "NoDomMM" ;
				}
			$mmperc = sprintf("%.3f",($mmcov / $totalcount)*100) ;
			print "$totalcount\t$refcov\t$mmflag\t$mmbase\t$mmcov\t$mmperc\t" ;
			} 

			if ($strandbase =~ /T/) {
				$refcov = @samplecov[3] ;
				if ((@samplecov[0] + @samplecov[1] + @samplecov[2]) < 1) {
					$mmbase = "None" ;
					$mmcov = 0 ;
					$mmflag = "NoMM" ;
				} elsif ((@samplecov[0] > @samplecov[1]) & (@samplecov[0] > @samplecov[2])) {
					$mmbase = "A" ;
					$mmcov = @samplecov[0] ;
					$mmflag = "DomMM" ;
				} elsif ((@samplecov[1] > @samplecov[0]) & (@samplecov[1] > @samplecov[2])) {
					$mmbase = "C" ;
					$mmcov = @samplecov[1] ;
					$mmflag = "DomMM" ;
				} elsif ((@samplecov[2] > @samplecov[0]) & (@samplecov[2] > @samplecov[1])) {
					$mmbase = "G" ;
					$mmcov = @samplecov[2] ;
					$mmflag = "DomMM" ;
				} else {
					$mmbase = "None" ;
					$mmcov = 0 ;
					$mmflag = "NoDomMM" ;
				}
				
			$mmperc = sprintf("%.3f",($mmcov / $totalcount)*100) ;
			print "$totalcount\t$refcov\t$mmflag\t$mmbase\t$mmcov\t$mmperc\t" ;
			} 
			print "$strandbase>$mmbase\t$basestring\t$CovA\t$CovC\t$CovG\t$CovT\t" ;

		}
				print "\n";
	}				
}
			

close(DATAFILE);





