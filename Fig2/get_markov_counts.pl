#!/usr/bin/perl

use DBI;

$species = "fly6_07";
$region = shift;
$file = shift;

$getConditionals = 1 if $file !~ /6to8mer/; # get if 1-2mers, don't get if 6-8mers
$minkmersize = ($getConditionals) ? 1 : 6;
$code = (split /\//, $file)[-1];
$code =~ s/(\w+).bls/$1/;
if (!$getConditionals) { ($bin) = (split /\s/, `grep -m1 -P '^$code\t' $species/$region/allgenes.bins`)[-1]; }
$maxbl = 7;

open IN, "<$file";
while(<IN>){ chomp;
	($kmer, $bls) = split;
	($klen, $j) = (length($kmer), 0);
	if ($klen == $minkmersize){ ++$pos; $prevbls = 0; }
	$cons = $kmer.$bin;
	$noncons = substr($kmer, 0, $klen-1).lc(substr($kmer,-1)); #add lowercase rather than uppercase if nonconserved
	for ($i = 0; sprintf("%.2f", $i) <= $maxbl; $i += 0.05){
		if ($bls >= $i){ $seq{$cons}[$j]++; }
		elsif ($getConditionals && $prevbls >= $i){ $seq{$noncons}[$j]++; }
		$seq{$noncons}[$j] = 0 if $getConditionals && ! defined $seq{$noncons}[$j];
		$seq{$cons}[$j] = 0 if ! defined $seq{$cons}[$j];
		++$j;
	}
	$prevbls = $bls;
}
close IN;

foreach $i (keys %seq){ print "$i @{$seq{$i}}\n"; }
