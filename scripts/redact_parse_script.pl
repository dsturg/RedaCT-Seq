#!/usr/bin/perl

# PROGRAM  :readact_parse_script
# AUTHOR   : Dave Sturgill

# Usage: redact_parse_script.pl filename samplenum

# eg scripts/redact_parse_script.pl sampledata/mpileup_output_chr19_min10.txt 3 > sampledata/mpileup_output_chr19_min10_parsed.txt ;



	
##########################
# Read command line args
##########################

my $filename = $ARGV[0];
chomp $filename;

my $numsamples = $ARGV[1];
chomp $numsamples;


#####################################
# Subroutine for printing fasta
#####################################


sub print_sequence {
	my($sequence, $length) = @_;
	for (my $pos=0 ; $pos < length($sequence) ; $pos += $length) {
		print substr($sequence, $pos, $length), "\n";
	}
}


######################
#@@@@@@@@@@@@@@@@@@@@@
######################

my $datafile;
my @columns ;
my $prestring ;
my $refbase ;


print "chr\tloc\tref\tsub" ;
my $i ;
my $j ;
#for ($i in 0:$numsamples) {
for ($i = 0 ; $i < $numsamples ; $i++) {
	print "\tdepth_",$i,"\tRef_",$i,"\tSub_",$i ; 
}
print "\n" ;

my $startpoint ;
my @depth ;
my @cov_A ;
my @cov_T ;
my @cov_C ;
my @cov_G ;

my $linecount = 0 ;

open(DATAFILE, "$filename");
while ($datafile = <DATAFILE>) {
	++$linecount;
	if ($linecount > 1) {

		@columns = split(' ', $datafile);
	
		#print "$columns[0]\t$columns[1]\t$columns[2]\n";
		$refbase = uc($columns[2]);
		$prestring = "$columns[0]\t$columns[1]\t$refbase";
		#print $prestring, "\t",
		# For each ACTG, get the get the coverage for each substitution
	
		# Get the coverages
		@depth = ();
		@cov_A = () ;
		@cov_T = () ;
		@cov_C = () ;
		@cov_G = () ;
		
		#~~~~~~~~~~~~~~~~~~~~
		# Iterate over each sample
		#~~~~~~~~~~~~~~~~~~~~
	
		for ($j = 0 ; $j < $numsamples ; $j++) {
			$startpoint = 3 + (11 * $j) ; # This changes with each sample
			push(@depth, $columns[$startpoint]) ;
			push(@cov_A,$columns[($startpoint + 1)] + $columns[($startpoint + 5)]) ;
			push(@cov_T,$columns[($startpoint + 2)] + $columns[($startpoint + 6)]) ;
			push(@cov_C,$columns[($startpoint + 3)] + $columns[($startpoint + 7)]) ;
			push(@cov_G,$columns[($startpoint + 4)] + $columns[($startpoint + 8)]) ;		
		
		}
		# The refbase determines what 'substitutions' to report
		if ($refbase =~ /A/) {
			# Print T mismatches
			print $prestring, "\tT";
			for ($j = 0 ; $j < $numsamples ; $j++) {
			  print "\t$depth[$j]\t$cov_A[$j]\t$cov_T[$j]";
			}
			print "\n";
			# Print C mismatches
			print $prestring, "\tC";
			for ($j = 0 ; $j < $numsamples ; $j++) {
			  print "\t$depth[$j]\t$cov_A[$j]\t$cov_C[$j]";
			}
			print "\n";
			# Print G mismatches
			print $prestring, "\tG";
			for ($j = 0 ; $j < $numsamples ; $j++) {
			  print "\t$depth[$j]\t$cov_A[$j]\t$cov_G[$j]";
			}
			print "\n";

		}		
		if ($refbase =~ /T/) {
			# Print A mismatches
			print $prestring, "\tA";
			for ($j = 0 ; $j < $numsamples ; $j++) {
			  print "\t$depth[$j]\t$cov_T[$j]\t$cov_A[$j]";
			}
			print "\n";
			# Print C mismatches
			print $prestring, "\tC";
			for ($j = 0 ; $j < $numsamples ; $j++) {
			  print "\t$depth[$j]\t$cov_T[$j]\t$cov_C[$j]";
			}
			print "\n";
			# Print G mismatches
			print $prestring, "\tG";
			for ($j = 0 ; $j < $numsamples ; $j++) {
			  print "\t$depth[$j]\t$cov_T[$j]\t$cov_G[$j]";
			}
			print "\n";

		}		
		if ($refbase =~ /C/) {
			# Print A mismatches
			print $prestring, "\tA";
			for ($j = 0 ; $j < $numsamples ; $j++) {
			  print "\t$depth[$j]\t$cov_C[$j]\t$cov_A[$j]";
			}
			print "\n";
			# Print T mismatches
			print $prestring, "\tT";
			for ($j = 0 ; $j < $numsamples ; $j++) {
			  print "\t$depth[$j]\t$cov_C[$j]\t$cov_T[$j]";
			}
			print "\n";
			# Print G mismatches
			print $prestring, "\tG";
			for ($j = 0 ; $j < $numsamples ; $j++) {
			  print "\t$depth[$j]\t$cov_C[$j]\t$cov_G[$j]";
			}
			print "\n";

		}		
		if ($refbase =~ /G/) {
			# Print A mismatches
			print $prestring, "\tA";
			for ($j = 0 ; $j < $numsamples ; $j++) {
			  print "\t$depth[$j]\t$cov_G[$j]\t$cov_A[$j]";
			}
			print "\n";
			# Print T mismatches
			print $prestring, "\tT";
			for ($j = 0 ; $j < $numsamples ; $j++) {
			  print "\t$depth[$j]\t$cov_G[$j]\t$cov_T[$j]";
			}
			print "\n";
			# Print C mismatches
			print $prestring, "\tC";
			for ($j = 0 ; $j < $numsamples ; $j++) {
			  print "\t$depth[$j]\t$cov_G[$j]\t$cov_C[$j]";
			}
			print "\n";

		}		
	}				
}
			

close(DATAFILE);




