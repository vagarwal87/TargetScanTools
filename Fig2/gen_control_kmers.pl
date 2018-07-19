#!/usr/bin/perl

use allfxns;

$species = "fly6_07";
$region = shift;

$dir = "$species/$region";
$order = 2; #markov order + 1
$outdir = "$dir/logprobs.REV";
mkdir "$outdir" if ! -d "$outdir";

open FH, "<$dir/kmercounts1to2mer.REV.txt"; <FH>;
while(<FH>){ chomp;
	@counts = (split /\s/, $_);
	$kmer = shift @counts;
	$counts{$kmer} = [@counts];
}
close FH;

foreach $mer (keys %counts){
	next if uc($mer) ne $mer;

	$kmersize = () = ($mer =~ /A|C|G|T/g);
	$noncons = substr($mer,0,length($mer)-1).lc(substr($mer,-1));

	$numcols = scalar(@{$counts{$mer}});

	for ($i = 0; $i <= $numcols; $i++){
		if (@{$counts{$mer}}[$i] > 0){ #get log of conditional probability
			if ($kmersize == 1){ $kmerprobs{$mer}[$i] = log2( @{$counts{$mer}}[$i]/@{$counts{$mer}}[0] ); }
			else { $kmerprobs{$mer}[$i] = log2( @{$counts{$mer}}[$i]/(@{$counts{$mer}}[$i]+@{$counts{$noncons}}[$i]) ); }
		}
		else{ $kmerprobs{$mer}[$i] = -99; }
	}
}

foreach $len (6..8){
	getkmers($len, '');
}

sub getkmers{
	local ($kmersize, $initframe) = @_;
	$infile = "all_".$kmersize."mers.txt";
	for ($i = 0; $i < $numcols; $i++){
		$bls = sprintf("%.2f", $i*0.05);
		open OUT, ">$outdir/$bls.$kmersize.probs";
		open FH, "<$infile";
		while($seq = <FH>){ chomp $seq;
			$logprob = 0;
			for ($j = 1-$order; $j <= length($seq)-$order; ++$j){
				$pos = max(0, $j);
				$kmer = substr($seq, $pos, $order-($pos-$j));
				$logprob += $kmerprobs{$kmer}[$i];
			}
			$logprob = sprintf("%.5f", $logprob);
			print OUT "$seq $logprob\n";
		}
	}
	close FH;
	close OUT;
}
