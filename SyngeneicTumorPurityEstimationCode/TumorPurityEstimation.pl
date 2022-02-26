########################################################################################################################
#Description: estimating tumor purity of syngeneic tumors using maximal likelihood method
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date:   Aug-2021
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
#
#
#This script implements the maximum likelihood algorithm described in the section 
#"Tumor purity estimation for syngeneic models" of MATERIALS AND METHODS
#
#It requireds three input files:
#	--mouse_strain_mutfre: record the SNP identity and frequency (always one for a nucleotide) in a mouse strain
#	--syngeneic_mutfreq:   record the SNP identity and frequency of somatic mutations in a syngeneic model
#
#	the above two files contains same SNPs
#	
#	--tumor_depth_file:    this is the read depth file for a syngeneic tumor profiled by a deep NGS assay reported in
#	Chen, X., Qian, W., Song, Z., Li, Q.X. & Guo, S. Authentication, characterization and contamination detection of 
#	cell lines, xenografts and organoids by barcode deep NGS sequencing. NAR Genom Bioinform 2, lqaa060 (2020).			
#
#Running environment:
#	Windows/Linus OS with perl installation (tested version: 5.28.1)
########################################################################################################################
#!/usr/bin/perl -w
use strict;

my $MINREAD=20; 		#minimal number of total reads for a SNP site to be considered
my $MEDIANREADS;		#to be determined from data
my $MEDIANSCALING=5.0;	#To prevent the distortion of  likelihood calculation by high-depth SNPs
						#we put a depth cap which is five times of median depth.

if(@ARGV!=3){
	print "usage: program <mouse_strain_mutfre> <syngeneic_mutfreq> <tumor_depth_file>\n";
	print "usage: program EMT6.mutfre EMT6.somatic.mutfre EMT6-1.depth\n";
	exit;
}

#--read in mouse strain mutfre--#
my %ms_mutfreq;
my @ms_snp;
open FH, $ARGV[0] or die;
while(<FH>){
	chomp;
	next if(/^\s*$|^pos/);
	my @s = split(/\s+/);
	for(my $i=2;$i<6;$i++){
		if($s[$i]==0.0){
			$s[$i]=0.001;	#change 0 to 0.001 to prevent underflow;
		}
	}
	$ms_mutfreq{$s[0]}= [$s[2],$s[3],$s[4],$s[5]];	
	push(@ms_snp, $s[0]);
}
close FH;

#--read in syngeneic cell line mutfre--#
my %sc_mutfreq;
my @sc_snp;
open FH, $ARGV[1] or die;
while(<FH>){
	chomp;
	next if(/^\s*$|^pos/);
	my @s = split(/\s+/);
	
	if($s[2] eq 'NA'){
		;
	}else{
		for(my $i=2;$i<6;$i++){
			if($s[$i]==0){
				$s[$i]=0.001;	#change 0 to 0.001 to prevent underflow;
			}
		}
	}
	$sc_mutfreq{$s[0]}= [$s[2],$s[3],$s[4],$s[5]];	
	push(@sc_snp, $s[0]);
}
close FH;

#--to prevent a few SNPs from dominating the likelihood, we need to set a upper bound for max_read--#
my @reads;
my @tumor_snp;
open FH, $ARGV[2] or die;
while(<FH>){
	chomp;
	next if(/^\s*$|^pos|NA/);
	my @s = split(/\s+/);
	
	my $snp = $s[0];
	my $A=$s[2];
	my $T=$s[3];
	my $C=$s[4];
	my $G=$s[5];
	
	my $sum = $A+$T+$C+$G;
	next if($sum<$MINREAD);
	push(@reads, $sum);
	push(@tumor_snp, $snp);
}
close FH;
$MEDIANREADS = median(@reads);

#--get common snp--#
my @common_snp = @{   intersect(intersect(\@ms_snp, \@sc_snp), \@tumor_snp)  };
my %common_snp;
foreach(@common_snp){
	$common_snp{$_} = 1;
}

#--remove SNPs on X-chromosome--#
my %t;
foreach(keys %common_snp){
	if(!/^X|^x/){
		$t{$_}=1;
	}	
}
%common_snp = %t;

#--read in tumor depth and calculate log-likelihood--#
my $min_logLikelihood = -100000000000;
my $min_theta = -1;
my $validSNP=0;
for(my $theta=0.001; $theta<1.0; ){
	open FH, $ARGV[2] or die;
	my $logLikelihood=0.0;
	
	$validSNP=0;
	while(<FH>){
		chomp;
		next if(/^\s*$|^pos|NA/);
		my @s = split(/\s+/);

		my $snp = $s[0];
		my $A=$s[2];
		my $T=$s[3];
		my $C=$s[4];
		my $G=$s[5];
		
		#next if(!isInArray($snp, \@common_snp));
		next if(!exists $common_snp{$snp});
		
		my $sum = $A+$T+$C+$G;
		next if($sum<$MINREAD);
		if($sum>$MEDIANREADS){
			for(my $i=2; $i<6; $i++){
				$s[$i] = $s[$i] * $MEDIANREADS/$sum * $MEDIANSCALING;
			}			
		}

		$validSNP++;
		
		shift(@s);
		shift(@s);
				
		my @ms_mutfre = @{ $ms_mutfreq{$snp} };
		my @sc_mutfre = @{ $sc_mutfreq{$snp} };

		
		#--calculate log-likelihood
		for(my $i=0;$i<@s;$i++){
			my $n = $s[$i];
			next if($n<1);
			my $ms = $ms_mutfre[$i];
			my $sc = $sc_mutfre[$i];
			
			next if($sc eq 'NA');
			$logLikelihood += $n * log( $theta * $sc + (1.0 - $theta) * $ms);			
		}
	}
	close FH;
	
	if($logLikelihood>$min_logLikelihood){
		$min_logLikelihood = $logLikelihood;
		$min_theta = $theta;
	}
	
	$theta=$theta+0.001;
}

print "$ARGV[2]\t$min_theta\t$validSNP\n";



sub median {
    my (@data) = sort { $a <=> $b } @_;
    if ( scalar(@data) % 2 ) {
        return ( $data[ @data / 2 ] );
    } else {
        my ( $upper, $lower );
        $lower = $data[ @data / 2 ];
        $upper = $data[ @data / 2 - 1 ];
        return ( mean( $lower, $upper ) );
    }
}

sub mean {
    my (@data) = @_;
    my $sum;
    foreach (@data) {
        $sum += $_;
    }
    return ( $sum / @data );
}

sub intersect{
	my @a = @{$_[0]};
	my @b = @{$_[1]};
	my $e;
	my @union;
	my @isect;
	my %union;
	my %isect;
	foreach $e (@a) { $union{$e} = 1 }

	foreach $e (@b) {
		if ( $union{$e} ) { $isect{$e} = 1 }
		$union{$e} = 1;
	}
	@union = keys %union;
	@isect = keys %isect;	
	
	return \@isect;
}

sub isInArray{
	my $e = $_[0];
	my @a = @{ $_[1] };
	foreach(@a){
		if($_ eq $e){
			return 1;
		}
	}
	return 0;	
}