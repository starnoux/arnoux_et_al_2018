# This file dissect the vcf file
# Script from Sylvain Glémin (sylvain.glemin@univ-montp2.fr) / C++ version included in the read2snp software 
# J.David 1/03/2012

use warnings;
use strict;

# gives the arguments meaning the vcf file and the out file
my ($vcf ,$output, $seuil) = @ARGV;

if (not defined $vcf) {
  die "Need vcf location, you fool!\n This vcf must be without stars in it.\n E.g. Script.pl /folder/CSFILTER_file.vcf /folder/CSFILTER_file.out";
}
 
if (not defined $output) {
  die "Need a second argument vcf file name + .out\n";
}

open VCF , $vcf or die $! ;

open (SORTIE , ">$output" ) or die $! ;

my $a ;
my $i ;
my %contig_polym ;
my %contig_htz;
my %contig_obs ;
my %contig_alanc ;
my %contig_Inbreed ;

my %cp;

while(<VCF>) {
# on passe les premieres lignes si elles servent à rien
if ($_!~/#/) { 
	chomp($_) ;

	my@line=split/\t/,$_;
	#print "$line[0]\n";
	
	# le nom du contig est en $line[0]
	#print "$line[0] $line[1] $line[3] $line[4]  $line[9]\n"; getc ;
	if (length($line[3])==1 and length($line[4])==1) {
	
	print SORTIE "$line[0],$line[1]";
	for ( $i=9 ; $i<=$#line ; $i++) {
	# print "$line[$i]" ; getc ;
	
	if (length($line[$i])<=5) { print SORTIE ",0,0,0,0"; } 
	else {
	 	my @line2=split/:/,$line[$i];
	    my @lineA=split/\//,$line2[0];
		my @lineN=split/,/,$line2[1];
		
		#print "$lineA[0],$lineA[1],"; 
		#print "$lineN[0],$lineN[1]\n" ;
		
		$cp{"A"}=0 ; $cp{"T"}=0 ;$cp{"C"}=0 ;  $cp{"G"}=0 ;
		$cp{$line[3]}=$lineN[0];
        $cp{$line[4]}=$lineN[1] ;
		print SORTIE ",$cp{'A'},$cp{'C'},$cp{'G'},$cp{'T'}"; 
		}
	 }
	 print SORTIE "\n"; 
	 }
	 }
	 }

close VCF;
close SORTIE;
