#!/usr/bin/perl

use allfxns;
use Env; Env::import();

$dist = 30;
%utrBase = %{readFasta($SEQBASE)};
`mkdir $SITES/` unless (-d "$SITES");

open MIRSEQFILE, "<$MIRSEQFILE" or die "miRseq file not opening";

while(<MIRSEQFILE>){ chomp($_);
	($exp, $miR) = split /\t/, $_;
	$seed = revCom(substr($miR, 1, 7)).'A';
	$seven8 = substr($seed,0,7);
	$sevenA = substr($seed,1,7);
	$sixmer = substr($seed,1,6);
	$offset6 = substr($seed,0,6);
	$sixA1 = substr($seed,2,6);

	$outfile = "$SITES/$exp\_flank.out";

	open FILE, "<targetsites/$exp\_sites.out" or die "can't open diff file targetsites/$exp\_sites.out";
	<FILE>;
	`rm -f $outfile` if (-e "$outfile");
	open OUT, ">$outfile" or die "can't open out file";
	while(<FILE>){ chomp($_);
		($utr_id, $fc_val) = (split /\t/, $_)[0,1];
		$seq = $utrBase{$utr_id}; next if $seq eq '';
		$seqlen = length($seq);

		$pos6a1 = index($seq, $sixA1);
		$poso6  = index($seq, $offset6);
		$pos6   = index($seq, $sixmer);
		$pos7m8 = index($seq, $seven8);
		$pos7a1 = index($seq, $sevenA);
		$pos8   = index($seq, $seed);

		$minpos = ($poso6 != -1)? $poso6 : $seqlen;
		if($pos6 <= $minpos && $pos6 != -1){ $minpos = $pos6; }
		if($pos6a1 <= $minpos && $pos6a1 != -1){ $minpos = $pos6a1; }

		while ($minpos != $seqlen) {
			if    (   $pos8 == $minpos ){ $site_type = '8m';  $pos = $pos8;     $begpos = max(0,$pos-$dist); $subseq = substr($seq,$begpos,$pos-$begpos+8+$dist); }
			elsif ( $pos7m8 == $minpos ){ $site_type = '7m8'; $pos = $pos7m8;   $begpos = max(0,$pos-$dist); $subseq = substr($seq,$begpos,$pos-$begpos+8+$dist); }
			elsif ( $pos7a1 == $minpos ){ $site_type = '7A1'; $pos = $pos7a1-1; $begpos = max(0,$pos-$dist); $subseq = substr($seq,$begpos,$pos-$begpos+8+$dist); }
			elsif (  $poso6 == $minpos ){ $site_type = 'o6';  $pos = $poso6;    $begpos = max(0,$pos-$dist); $subseq = substr($seq,$begpos,$pos-$begpos+8+$dist); }
			elsif (   $pos6 == $minpos ){ $site_type = '6m';  $pos = $pos6-1;   $begpos = max(0,$pos-$dist); $subseq = substr($seq,$begpos,$pos-$begpos+8+$dist); }
			else 				    { $site_type = '6A1'; $pos = $pos6a1-2; $begpos = max(0,$pos-$dist); $subseq = substr($seq,$begpos,$pos-$begpos+8+$dist); }

			print OUT join ("\t", $utr_id, sprintf('%0.2f',$fc_val), $pos, length($seq), substr($seq,$begpos,$pos-$begpos), substr($seq,$pos,8), substr($seq,$pos+8,$dist), $site_type), "\n" 
				if $pos >= 14 && length(substr($seq,$pos,8))==8; #outside ribosomal shadow & all 8 bp exist at least

			$pos8   = index($seq, $seed, $minpos+5);
			$pos7m8 = index($seq, $seven8, $minpos+5);
			$pos7a1 = index($seq, $sevenA, $minpos+6);
			$pos6   = index($seq, $sixmer, $minpos+6);
			$pos6a1 = index($seq, $sixA1, $minpos+7);
			$poso6  = index($seq, $offset6, $minpos+5);

			$minpos = ($poso6 != -1)? $poso6 : $seqlen;
			if($pos6 <= $minpos && $pos6 != -1){ $minpos = $pos6; }
			if($pos6a1 <= $minpos && $pos6a1 != -1){ $minpos = $pos6a1; }
		}
	}
	close FILE; close OUT;
}
close MIRSEQFILE;
