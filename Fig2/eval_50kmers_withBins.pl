#!/usr/bin/perl

use allfxns;

$controls = 50;
$numbins = 10;

$species = "fly6_07";
$region = shift;
$sitetype = shift;
$file = shift;

$dir = "$species/$region";
($bls) = ($file =~ /(\d+\.\d\d)/);
open FH, "zcat $dir/kmercounts6to8mer.REV.txt.gz | " or die "no $dir/kmercounts6to8mer.REV.txt.gz"; <FH>; #_binnedtrees
while(<FH>){ chomp;
	@counts = (split /\s/, $_);
	$kmer = shift @counts;
	if($coding){
		($kmer, $thisframe) = ($kmer =~ /(\w+)(\d)/);
		next if $thisframe ne $frame;
	}
	$counts{$kmer} = $counts[round($bls/0.05)];
	$refnum{$kmer} = $counts[0];
	$ratio{$kmer} = $counts{$kmer}/$refnum{$kmer};
}
close FH;

open FH, "<$species/consSeeds/cons$sitetype.txt";
while(<FH>){ chomp;
	if ($sitetype =~ /8m|7m8/){
		$seeds{substr($_, 0, 6)} = 1; $seeds{substr($_, 1, 6)} = 1;
		$seeds{substr($_, 2, 6)} = 1 if $sitetype =~ /8m/ && $species =~ /fly/; #control for 6A1 site
	}
	elsif ($sitetype =~ /7A1/){
		$seeds{substr($_, 0, 6)} = 1;
		$seeds{substr($_, 1, 6)} = 1 if $sitetype =~ /7A1/ && $species =~ /fly/; #control for 6A1 site
	}
	$seeds{$_} = 1;	#get 8-mers
}
close FH;

open FH, "sort -k2,2nr $file | ";	#get sorted log-probabilities of kmer
while(<FH>){ chomp; ($kmers[$i], $probs[$i]) = split; ++$i; } #$probs[$i] = 10**$probs[$i];	default: subtract in log space (to get ratio)
close FH;

for ($i = 0; $i != $#kmers; ++$i){
	for $bin (1..$numbins){
		$refnumsummed{$kmers[$i]} += $counts{$kmers[$i].$bin};
		$ratio{$kmers[$i].$bin} = 1 if ! exists $ratio{$kmers[$i].$bin};
	}
}

for ($i = 0; $i != $#kmers; ++$i){
	next if $seeds{$kmers[$i]} != 1;
	($j, $k, $signal) = ($i-1, $i+1, int($refnumsummed{$kmers[$i]})); #int in case kmer doesnt exist in counts
	$hasPUF = () = ($kmers[$i] =~ /TGT[AG]/);	#count 1 if PUF site exists
	$cgnucs = () = ($kmers[$i] =~ /C|G/g);	#count C + G occurrences in seed
	$cpgnucs = () = ($kmers[$i] =~ /CG/g);	#count CpG occurrences in seed
	@chosen = ();

	for ($ctl = 0; $ctl < $controls; ++$ctl){
		--$j while ( filter($kmers[$j]) && $j >=    0   );
		++$k while ( filter($kmers[$k]) && $k <= $#kmers);
		if($j >= 0 && ($ctl % 2 || $k > $#kmers)){ #use rank-based method to select controls
			for $bin (1..$numbins){ $back[$num][$ctl] += $ratio{$kmers[$j].$bin}*$refnum{$kmers[$i].$bin}; }
			push @chosen, "$kmers[$j]"; --$j;
		}
		elsif ($k < $#kmers || ($j < 0 && $k < $#kmers)) {
			for $bin (1..$numbins){ $back[$num][$ctl] += $ratio{$kmers[$k].$bin}*$refnum{$kmers[$i].$bin}; }
			push @chosen, "$kmers[$k]"; ++$k;
		}
	}
	die "need more controls for $kmers[$i]! $j $k found: ".scalar(@{$back[$num]}). "controls only" if scalar(@{$back[$num]}) < $controls;

	$back = sprintf("%.2f", mean($back[$num]));	#calculate signal above background
	fisher_yates_shuffle($back[$num]);
	$sab = sprintf("%.2f", $signal - $back); #signal above background
	$totsignal += $signal;
	print "$kmers[$i]\t$sab\t$signal-$back\t@chosen\n";
	++$num;
}

for ($i = 0; $i < $controls; ++$i){
	@col = map {$$_[$i]} @back;	#get column $i of 2D signal above background array, representing 50 cohorts
	$ind_back[$i] = sum(\@col);
}

$totback=round(mean(\@ind_back));
(@sab, @s2b)=((),());
for ($i = 0; $i < $controls; ++$i){ $sab[$i]=$totsignal-$ind_back[$i]; }
for ($i = 0; $i < $controls; ++$i){ $s2b[$i]=($totsignal+1)/($ind_back[$i]+1); } #pseudocount added to prevent division by zero errors at high branch cutoffs

$varsab = sprintf("%.10f",variance(@sab)/$controls); #variance in sig above background
$vars2b = sprintf("%.10f",variance(@s2b)/$controls); #variance in sig:background ratio
print "$totsignal\t$totback\t$varsab\t$vars2b\n";

sub filter{
	$kmer = shift;
	$kmerPUF = () = ($kmer =~ /TGT[AG]/);
	$kmercgnucs = () = ($kmer =~ /C|G/g);	# count C + G occurrences in kmer
	$kmercpgnucs = () = ($kmer =~ /CG/g);	# count CpG occurrences in kmer
	return ($seeds{substr($kmer, 0, 6)} 	# should not contain any 6-mers that match a miR
	     || $seeds{substr($kmer, 1, 6)} 	# should not contain any 6-mers that match a miR
	     || $seeds{substr($kmer, 2, 6)} 	# should not contain any 6-mers that match a miR
	     || (substr($kmer, -1, 1) ne "A" 	# A at position 1 filter
	     	   && $sitetype =~ /8m/)
	     || (substr($kmer, -1, 1) ne "A" 	# A at position 1 filter
	     	   && $sitetype =~ /7A1/) #	     	   && $kmerPUF==0 breaks if use for 6A1 as some are purely A/U in fly
	     || ($hasPUF != $kmerPUF			# PUF sites filter
	     	   && $sitetype =~ /8m|7m8/		# impossible constraint for all 7mer-A1, 6mer, or offset-6mer sites
	     	   && $kmercgnucs < 5)
	     || ($hasPUF != $kmerPUF			# PUF sites filter
	     	   && $sitetype =~ /7A1/
	     	   && $kmercgnucs < 4
	     	   && $species =~ /fly/)
	     || ! exists $refnumsummed{$kmer}		# kmer should exist in ref species at least
	     || $kmercgnucs  != $cgnucs		# C + G content should match	     
	     || ($kmercpgnucs != $cpgnucs	# CpG counts should match if working w/ vertebrate clade
	     	   && $species =~ /human|mouse/)
	);
}
