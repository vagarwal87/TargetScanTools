#!/usr/bin/perl

use allfxns;

open IN, "<mirbase-v21-fly.fa";
while(<IN>){
	$mir = substr((split /\s/, $_)[0],1);
	$seq = <IN>; chomp $seq;
	$ismiR{substr($seq, 0, 20)} = 1;
	$seq = substr($seq, 1, 7);
	$miRs{$seq} = $mir if $miRs{$seq} eq '';
}
close IN;

open IN, "zcat GSM280088_AGO1_IP_S2_raw_data.txt.gz |";
while(<IN>){
	chomp;
	($seq, $count) = (split);
	$seed = substr($seq, 1, 7);
	$mirs{$seed} += 1 if $ismiR{substr($seq, 0, 20)} == 1;
}
close IN;

$count=0;
foreach (sort {$mirs{$b} <=> $mirs{$a}} keys %mirs){
	if ($miRs{$_} ne ''){
		$seed = revCom($_);
		$seed2miR{$seed} = $miRs{$_};
		$seed2miR{substr($seed,1)."A"} = $miRs{$_};
		$seedmatch{$seed} = 1;
		$seedmatch{substr($seed,1)."A"} = 1;
		print "$_\t$miRs{$_}\t$mirs{$_}\n";
		++$count;
	}
}
