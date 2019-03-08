#!/usr/bin/perl -w
use strict;
die "Usage:perl $0 <beagle2treemix output> <outputprefix.geno>\n" unless @ARGV == 2;
open (IN, ($ARGV[0] =~ /\.gz$/)? "gzip -dc $ARGV[0] |" : $ARGV[0]) or die "gzip treemix file missing!\n";
open (OUT, ">$ARGV[1].geno") or die "eigenstrat geno file writing failed!\n";
<IN>;
while (<IN>){
	chomp;
	my @tmp = split/\s+/;
	foreach my $tmp_geno (@tmp) {
		if ($tmp_geno eq "2,0") {
			print OUT "2";
		}
		elsif ($tmp_geno eq "1,1") {
			print OUT "1";
		}
		elsif ($tmp_geno eq "0,2") {
			print OUT "0";
		}
		else{
			die "$tmp_geno wrong!\n";
		}
	}
	print OUT "\n";
}
close IN;
close OUT;
