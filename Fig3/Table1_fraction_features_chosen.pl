#!/usr/bin/perl

$x=10*2;
$y=3*2;
$z=$x+$y;

%counts = (); $totcounts = 0;
@heads = split "\n", `grep -A$x 'traincounts' data/fly8mer7mer-1.out | grep -v -e '--' -e traincounts | awk '{if (count++%2==0){\$1=\"\"; print \$0;}}' | tr -d '\"' | tr -d ' '`;
open IN, "grep -A$z 'traincounts' data/fly8mer7mer-*.out | grep -v -e '--' -e traincounts | awk '{if (count++%2==1){\$1=\"\"; print \$0;}}' | tr -d '\"' | tr -d ' ' | ";
open OUT, ">fly-results.out";

print OUT join("\t", @heads), "\n";
while(<IN>){ chomp; $_ =~ s/ //g;
	print OUT $_;
	for $i (1..($x/2 - 1)){ $_ = <IN>; chomp; print OUT "\t$_"; }
	print OUT "\n";
	for $i (1..($y/2)){
		$_ = <IN>; chomp;
		if ($i == 3){ #2 is BIC, 3 is AIC
			$_ =~ s/ |\\n|I\(|\)|==|\\|\(|TRUE//g; #characters to remove
			@a = split /\+/, $_;
			foreach (@a){ $_=join (":", sort { $a cmp $b } (split /:/, $_)); $counts{$_}++; }
			$totcounts++;
		}
	}
}
close IN;
foreach (sort {$counts{$b} <=> $counts{$a}} keys %counts){ print "$_\t", sprintf("%.2f", $counts{$_}/$totcounts), "\n" if $counts{$_}/$totcounts >= 0; }
close OUT;
