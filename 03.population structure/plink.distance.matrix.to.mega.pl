#!/usr/bin/perl
use strict;
die "Usage : perl $0 <mdist.id> <mdist> <individual number> <out prefix>\n" unless @ARGV==4;
open(IN1, "$ARGV[0]") or die "mdist.id file required!\n";
open(IN2, "$ARGV[1]") or die "mdist file required!\n";
open(OUT, ">$ARGV[3]\.meg") or die "permission denied!\n";
print OUT "#mega!\nTitle: Concatenated Files;\n!Format DataType=Distance DataFormat=LowerLeft NTaxa=$ARGV[2];\n!Description\n No. of Taxa : $ARGV[2]\n Gaps/Missing data : Pairwise Deletion\n Codon Positions : 1st+2nd+3nd+Noncoding\n Distance method : Nucleotide: Tamura-Nei  [Pairwise distances]\n d : Estimate\n;\n";
my $individual = 0;
my %hash;
while(<IN1>){
	chomp;
	$individual++;
	my @ind = split/\s+/;
	print OUT "[$individual] #$ind[1]\n";
}
die "something wrong with the individual number!\n" unless $individual eq $ARGV[2];
print OUT "[";
for my $tmp (1..$individual){
	print OUT " $tmp";
}
print OUT " ]\n";
my $index = 0;
while(<IN2>){
	chomp;
	$index++;
	my @dis = split/\s+/;
	print OUT "[$index]";
	print OUT "  $hash{$index}";
	foreach my $tmp (0..($index-2)){
		print OUT " $dis[$tmp]";
	}
	print OUT "\n";
}
close OUT;
close IN2;
