#!/usr/bin/perl

use allfxns;

$home = $ENV{DIR};

$species = "fly6_07";
$region = shift;

$maxbl = 7;

@files = <$species/$region/counts1to2mer.REV/*.counts>;

foreach $file (@files){
	open IN, "<$file";
	while(<IN>){ chomp;
		@counts = (split /\s/, $_);
		$kmer = shift @counts;
		vecsum($counts{$kmer}, \@counts);
	}
	close IN;
}

open OUT, ">$species/$region/kmercounts1to2mer.REV.txt";
select OUT;
for($i = 0; $i <= $maxbl; $i = sprintf("%.2f", $i+0.05)){ push(@header, $i); }
print join("\t", "kmer", @header), "\n";
foreach $i (sort {$a cmp $b} keys %counts){ print join("\t", $i, @{$counts{$i}}), "\n"; }
close OUT;
