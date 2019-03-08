#!/usr/bin/perl -w
use strict;
die "Usage: perl $0 <block SFS> <output>\n" unless @ARGV == 2;
open(IN, "$ARGV[0]") or die "block SFS missing!\n";
open(OUT, ">$ARGV[1]") or die "permission denied!\n";
my %SFS;
while(<IN>){
	chomp;
	my @entry = split/\s+/;
	foreach my $index (0..$#entry){
		$SFS{$index} += $entry[$index];
	}
	print STDERR "entry:$#entry+1\n";
}
close IN;
my @SFS;
foreach my $index(sort {$a <=> $b} keys %SFS){
	push @SFS, $SFS{$index};
}
print OUT join(" ", @SFS),"\n";
