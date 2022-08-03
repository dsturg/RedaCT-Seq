#!/usr/bin/perl

# PROGRAM  :readact_parse_script
# AUTHOR   : Dave

# Usage: redact_parse_script.pl filename samplenum

# eg perl redact_parse_script.pl mpileup_results_test.txt 3 > mpileup_parsed_test.txt
# perl redact_parse_script.pl SasChen_mpileup_pA_total_merge.txt 3 > SasChen_mpileup_pA_total_merge_parsed.txt

# Note in my example, there are 4 sampels in the mpileup, but I only parse the first 3

#use strict;
	
##########################
# Read command line args
##########################

my $filename = $ARGV[0];
chomp $filename;

my $numsamples = $ARGV[1];
chomp $numsamples;

#####################################
# Usage / Notes
#####################################
=pod
# Attempt to re-code my awkward AWK into a more generalizable perl 

#awk 'BEGIN{i=0}{if (i < 1) { print "chr\tloc\tref\tsub\tdepth_0\tRef_0\tSub_0\tdepth_1\tRef_1\tSub_1"; i++} else {if (toupper($3) ~ /A/) \
#{ print $1"\t"$2"\t"$3"\tT\t"$4"\t"$5+$9"\t"$6+$10"\t"$15"\t"$16+$20"\t"$17+$21"\n"  \
#$1"\t"$2"\t"$3"\tC\t"$4"\t"$5+$9"\t"$7+$11"\t"$15"\t"$16+$20"\t"$18+$22"\n" \
#$1"\t"$2"\t"$3"\tG\t"$4"\t"$5+$9"\t"$8+$12"\t"$15"\t"$16+$20"\t"$19+$23} else {if (toupper($3) ~ /T/) \
#{ print $1"\t"$2"\t"$3"\tA\t"$4"\t"$6+$10"\t"$5+$9"\t"$15"\t"$17+$21"\t"$16+$20"\n" \
#$1"\t"$2"\t"$3"\tC\t"$4"\t"$6+$10"\t"$7+$11"\t"$15"\t"$17+$21"\t"$18+$22"\n" \
#$1"\t"$2"\t"$3"\tG\t"$4"\t"$6+$10"\t"$8+$12"\t"$15"\t"$17+$21"\t"$19+$23} else {if (toupper($3) ~ /C/) \
#{ print $1"\t"$2"\t"$3"\tA\t"$4"\t"$7+$11"\t"$5+$9"\t"$15"\t"$18+$22"\t"$16+$20"\n" \
#$1"\t"$2"\t"$3"\tT\t"$4"\t"$7+$11"\t"$6+$10"\t"$15"\t"$18+$22"\t"$17+$21"\n"  \
#$1"\t"$2"\t"$3"\tG\t"$4"\t"$7+$11"\t"$8+$12"\t"$15"\t"$18+$22"\t"$19+$23} else {if (toupper($3) ~ /G/) \
#{ print $1"\t"$2"\t"$3"\tA\t"$4"\t"$8+$12"\t"$5+$9"\t"$15"\t"$19+$23"\t"$16+$20"\n" \
#$1"\t"$2"\t"$3"\tT\t"$4"\t"$8+$12"\t"$6+$10"\t"$15"\t"$19+$23"\t"$17+$21"\n"  \
#$1"\t"$2"\t"$3"\tC\t"$4"\t"$8+$12"\t"$7+$11"\t"$15"\t"$19+$23"\t"$18+$22} }}}}}' $1 ;


#chr loc ref depth A T C G a t c g Insertion Deletion depth A T C G a t c g Insertion Deletion 
#chrM 1 G 78 0 0 0 77 0 0 0 1 NA NA 2 0 0 0 2 0 0 0 0 NA NA 
#chrM 2 A 325 315 0 0 0 10 0 0 0 NA NA 3 3 0 0 0 0 0 0 0 NA NA 
#chrM 3 T 387 0 377 0 0 0 10 0 0 NA NA 4 0 4 0 0 0 0 0 0 NA NA 
#chrM 4 C 606 0 0 587 0 0 0 19 0 NA NA 5 0 0 5 0 0 0 0 0 NA NA 
=cut



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
			#print "j is $j startpoint is ", $startpoint, "\tadded $columns[$startpoint]\n" ;
			push(@cov_A,$columns[($startpoint + 1)] + $columns[($startpoint + 5)]) ;
			push(@cov_T,$columns[($startpoint + 2)] + $columns[($startpoint + 6)]) ;
			push(@cov_C,$columns[($startpoint + 3)] + $columns[($startpoint + 7)]) ;
			push(@cov_G,$columns[($startpoint + 4)] + $columns[($startpoint + 8)]) ;
			#push(@cov_A,$columns[($startpoint + 1)]) ;
			#push(@cov_T,$columns[($startpoint + 2)]) ;
			#push(@cov_C,$columns[($startpoint + 3)]) ;
			#push(@cov_G,$columns[($startpoint + 4)]) ;
		
		
		}
		#print "Here is my array:" ;
		#print join(", ", @depth);
		#print "Done\n" ;
	   #$j = 0;
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

=pod
chr	loc	ref	depth	A	T	C	G	a	t	c	g	Insertion	Deletion	depth	A	T	C	G	a	t	c	g	Insertion	Deletion	
chrM	1	G	78	0	0	0	77	0	0	0	1	NA	NA	2	0	0	0	2	0	0	0	0	NA	NA	
0       1   2   3   4   5   6   7   8   9   10  11  12  13  14  14
chr loc ref sub depth_0 Ref_0 Sub_0 depth_1 Ref_1 Sub_1
chrM 1 G A 78 78 0 2 2 0
chrM 1 G T 78 78 0 2 2 0


#awk 'BEGIN{i=0}{if (i < 1) { print "chr\tloc\tref\tsub\tdepth_0\tRef_0\tSub_0\tdepth_1\tRef_1\tSub_1"; i++} else {if (toupper($3) ~ /A/) \
#{ print $1"\t"$2"\t"$3"\tT\t"$4"\t"$5+$9"\t"$6+$10"\t"$15"\t"$16+$20"\t"$17+$21"\n"  \
#$1"\t"$2"\t"$3"\tC\t"$4"\t"$5+$9"\t"$7+$11"\t"$15"\t"$16+$20"\t"$18+$22"\n" \
#$1"\t"$2"\t"$3"\tG\t"$4"\t"$5+$9"\t"$8+$12"\t"$15"\t"$16+$20"\t"$19+$23} else {if (toupper($3) ~ /T/) \
#{ print $1"\t"$2"\t"$3"\tA\t"$4"\t"$6+$10"\t"$5+$9"\t"$15"\t"$17+$21"\t"$16+$20"\n" \
#$1"\t"$2"\t"$3"\tC\t"$4"\t"$6+$10"\t"$7+$11"\t"$15"\t"$17+$21"\t"$18+$22"\n" \
#$1"\t"$2"\t"$3"\tG\t"$4"\t"$6+$10"\t"$8+$12"\t"$15"\t"$17+$21"\t"$19+$23} else {if (toupper($3) ~ /C/) \
#{ print $1"\t"$2"\t"$3"\tA\t"$4"\t"$7+$11"\t"$5+$9"\t"$15"\t"$18+$22"\t"$16+$20"\n" \
#$1"\t"$2"\t"$3"\tT\t"$4"\t"$7+$11"\t"$6+$10"\t"$15"\t"$18+$22"\t"$17+$21"\n"  \
#$1"\t"$2"\t"$3"\tG\t"$4"\t"$7+$11"\t"$8+$12"\t"$15"\t"$18+$22"\t"$19+$23} else {if (toupper($3) ~ /G/) \
#{ print $1"\t"$2"\t"$3"\tA\t"$4"\t"$8+$12"\t"$5+$9"\t"$15"\t"$19+$23"\t"$16+$20"\n" \
#$1"\t"$2"\t"$3"\tT\t"$4"\t"$8+$12"\t"$6+$10"\t"$15"\t"$19+$23"\t"$17+$21"\n"  \
#$1"\t"$2"\t"$3"\tC\t"$4"\t"$8+$12"\t"$7+$11"\t"$15"\t"$19+$23"\t"$18+$22} }}}}}' $1 ;

=cut




