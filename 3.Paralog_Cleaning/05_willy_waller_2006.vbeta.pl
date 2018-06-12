# This script dissect the vcf file
# Script from Sylvain Glémin (sylvain.glemin@univ-montp2.fr) / C++ version included in the read2snp software 
# The given threshold value will exclude all heterozygous sites that have less than the threshold value reads over all the individuals. In our case if 'Seuil = 4' then the site must be present in at least 4 alleles of the total number of the population allele number. (Diploid -> total allele over population = individuals \* 2)

use warnings;
use strict;

#  vcf file
my ($vcf ,$output, $thresh) = @ARGV;

if (not defined $vcf) {
  die "Need vcf location, you fool!\n This vcf must be without stars in it.\n E.g. 05_willy_waller_2006.vbeta.pl /folder/file.vcf name_epluche.output.txt 4_or_otherThresholdNumber";
}
 
if (not defined $output) {
  die "Need a second argument with the epluche.txt name\n E.G. \n 05_willy_waller_2006.vbeta.pl /folder/file.vcf name_epluche.output.txt 4_or_otherThresholdNumber";
}

open VCF , $vcf or die $! ;

open (SORTIE , ">$output" ) or die $! ;

my $thresh_x =  $thresh || 4;

my $a ;
my $i ;
my %contig_polym ;
my %contig_htz;
my %contig_obs ;
my %contig_alanc ;
my %contig_Inbreed ;


while(<VCF>) {
# The first lines will be removed if useless
if ($_!~/#/) { 
	chomp($_) ;

	my@line=split/\t/,$_;
	print "$line[0]\n";
	
	# The contig name is in $line[0]
	if (!exists($contig_polym{$line[0]})) 
	 { $contig_polym{$line[0]}=1 ;} 
	 else 
	 {$contig_polym{$line[0]} ++;}
	 
	  my @Inb=split/;/,$line[7];
	 
	 for ($i =0 ; $i<$#Inb ; $i++ ) { if ($Inb[$i] =~"Inbr") 
	                                 {my @Coef=split/=/,$Inb[$i];
									 if(!exists($contig_Inbreed{$line[0]})) {
									        $contig_Inbreed{$line[0]} = $Coef[1];
									 } else
									 {$contig_Inbreed{$line[0]} += $Coef[1];}
									 }
									}
	 for ( $i=9 ; $i<$#line ; $i++) {
	 	my @line2=split/:/,$line[$i];
	    my @line3=split/\//,$line2[0];
		my $c1 ;
		my $c2 ;
		
		if (exists($line2[1])) { my @line4=split/,/,$line2[1];
		                       $c1=$line4[0] ;
							   $c2=$line4[1] ;
		
		# We keep only the individuals for whom the sum of the 2 X is higher tham the threshold thresh_X  	
			
		if ($c1+$c2>=$thresh_x) { if (!exists($contig_obs{$line[0]})) { $contig_obs{$line[0]}=1} 
		                        else { $contig_obs{$line[0]} ++}
								
								if (!exists($contig_alanc{$line[0]})) { $contig_alanc{$line[0]}=0} 
		                        
								$contig_alanc{$line[0]} += ($line3[1]==0)*1 + ($line3[0]==0) *1 ;
								
								if   ($line3[1]!=$line3[0])
								{ if (!exists($contig_htz{$line[0]})) { $contig_htz{$line[0]}=1} 
		                          else { $contig_htz{$line[0]} ++}
								}  
							}	
		}						
		} 
		
		}
	 
	 }


my $n ;
print SORTIE  "contig\tnb_snp\tnb_htz\tnbobs\tnb_all_ref\tInbreedCoef\n"; 

foreach $n (sort keys %contig_polym) { 
if (!exists($contig_htz{$n})) { $contig_htz{$n}=0} 
if (!exists($contig_obs{$n})) { $contig_obs{$n}=0} 
if (!exists($contig_Inbreed{$n})) { $contig_Inbreed{$n}=-999} 
$contig_Inbreed{$n}=$contig_Inbreed{$n}/$contig_polym{$n};

print SORTIE  "$n\t$contig_polym{$n}\t$contig_htz{$n}\t$contig_obs{$n}\t$contig_alanc{$n}\t$contig_htz{$n}/$contig_obs{$n}\t$contig_Inbreed{$n}\n"; }
