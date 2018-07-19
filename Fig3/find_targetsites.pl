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

	$outfile = "$SITES/$exp\_flank.out";

	open FILE, "<diff/$exp.diff" or die "can't open diff file diff/$exp\_sites.out";
	<FILE>;
	`rm -f $outfile` if (-e "$outfile");
	open OUT, ">$outfile" or die "can't open out file";
	while(<FILE>){ chomp($_);
		($utr_id, $fc_val) = (split /\t/, $_)[0,1];
		$seq = $utrBase{$utr_id}; next if $seq eq '';

		$num8 = () = ($seq =~ /$seed/g);
		$num7m8 = () = ($seq =~ /$seven8/g); $num7m8 -= $num8;
		$num7A1 = () = ($seq =~ /$sevenA/g); $num7A1 -= $num8;

		next if ($num8+$num7m8+$num7A1 != 1);

		$pos7m8 = index($seq, $seven8);
		$pos7a1 = index($seq, $sevenA);
		$pos8   = index($seq, $seed);

		if    (   $pos8 != -1 ){ $site_type = '8m';  $pos = $pos8;     $begpos = max(0,$pos-$dist); $subseq = substr($seq,$begpos,$pos-$begpos+8+$dist); }
		elsif ( $pos7m8 != -1 ){ $site_type = '7m8'; $pos = $pos7m8;   $begpos = max(0,$pos-$dist); $subseq = substr($seq,$begpos,$pos-$begpos+8+$dist); }
		elsif ( $pos7a1 != -1 ){ $site_type = '7A1'; $pos = $pos7a1-1; $begpos = max(0,$pos-$dist); $subseq = substr($seq,$begpos,$pos-$begpos+8+$dist); }

		print OUT join ("\t", $utr_id, sprintf('%0.2f',$fc_val), $pos, length($seq), substr($seq,$begpos,$pos-$begpos), substr($seq,$pos,8), substr($seq,$pos+8,$dist), $site_type), "\n" 
			if $pos >= 14 && length(substr($seq,$pos,8))==8; #outside ribosomal shadow (Grimson et al 2007) & all 8 bp exist at least
	}
	close FILE; close OUT;
}
close MIRSEQFILE;
