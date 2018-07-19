#!/usr/bin/perl

use allfxns;
use Env; Env::import();

`mkdir $SITES/` unless (-d "$SITES");
%utr3Base = %{readFasta($SEQBASE)};
%orfBase = %{readFasta($SEQBASEORF)};
%utr5Base = %{readFasta($SEQBASE5UTR)};

open MIRSEQFILE, "<$MIRSEQFILE" or die "miRseq file not opening";
while(<MIRSEQFILE>){
	chomp($_); ($exp, $miR) = split /\t/, $_;

	$seed = revCom(substr($miR, 1, 7)).'A';
	$seven8 = substr($seed,0,7);
	$sevenA = substr($seed,1,7);
	$sixmer = substr($seed,1,6);
	$offset6 = substr($seed,0,6);
	$sixA1 = substr($seed,2,6);
	$fivemer = substr($seed,2,5);
	$seedG = $seven8.'G';
	$seedC = $seven8.'C';
	$seedT = $seven8.'T';

	open FILE, "<diff_plsr_corrected/$exp.diff" or next; #OR diff_uncorrected
	open OUT, ">$SITES/$exp\_sites.out" or die "can't open out file";

	print OUT join ("\t", "utr_id", "fc", "orf_gc", "orf_len", "utr5_gc", "utr5_len", "utr3_gc", "utr3_len", "utr5p8", "utr5p7m8", "utr5p7A1", "utr5p6", "utr5po6", "utr5p6A1", "orf8", "orf7m8", "orf7A1", "orf6", "orfo6", "orf6A1", "utr8C", "utr8G", "utr8T", "utr8", "utr7m8", "utr7A1", "utr6", "utro6", "utr6A1", "utr5"), "\n";

	while(<FILE>){
		chop($_); ++$count;
		($utr_id, $fc_val) = split /\t/;
		$seq = $utr3Base{$utr_id}; next if $seq eq '';
		$seqorf = $orfBase{$utr_id};
		$seq5utr = $utr5Base{$utr_id};

		$num58 = () = ($seq5utr =~ /$seed/g);
		$num57m8 = () = ($seq5utr =~ /$seven8/g); $num57m8 -= $num58;
		$num57A1 = () = ($seq5utr =~ /$sevenA/g); $num57A1 -= $num58;
		$num56  = () = ($seq5utr =~ /$sixmer/g);  $num56 -= ($num57m8+$num57A1+$num58);
		$num5o6 = () = ($seq5utr =~ /$offset6/g); $num5o6 -= ($num57m8+$num58);
		$num56A1 = () = ($seq5utr =~ /$sixA1/g); $num56A1 -= ($num57A1+$num58);

		$numo8 = () = ($seqorf =~ /$seed/g);
		$numo7m8 = () = ($seqorf =~ /$seven8/g); $numo7m8 -= $numo8;
		$numo7A1 = () = ($seqorf =~ /$sevenA/g); $numo7A1 -= $numo8;
		$numo6  = () = ($seqorf =~ /$sixmer/g);  $numo6 -= ($numo7m8+$numo7A1+$numo8);
		$numoo6 = () = ($seqorf =~ /$offset6/g); $numoo6 -= ($numo7m8+$numo8);
		$numo6A1 = () = ($seqorf =~ /$sixA1/g); $numo6A1 -= ($numo7A1+$numo8);

		print OUT join ("\t", $utr_id, sprintf("%.3f", $fc_val), gcContent($seqorf), int(length($seqorf)), gcContent($seq5utr), int(length($seq5utr)), gcContent($seq), length($seq), $num58, $num57m8, $num57A1, $num56, $num5o6, $num56A1, $numo8, $numo7m8, $numo7A1, $numo6, $numoo6, $numo6A1, '');

		$num8 = () = ($seq =~ /$seed/g);
		$num8C = () = ($seq =~ /$seedC/g);
		$num8G = () = ($seq =~ /$seedG/g);
		$num8T  = () = ($seq =~ /$seedT/g);
		$num7m8 = () = ($seq =~ /$seven8/g); $num7m8 -= $num8;
		$num7A1 = () = ($seq =~ /$sevenA/g); $num7A1 -= $num8;
		$num6  = () = ($seq =~ /$sixmer/g);  $num6 -= ($num7m8+$num7A1+$num8);
		$numo6 = () = ($seq =~ /$offset6/g); $numo6 -= ($num7m8+$num8);
		$num6A1 = () = ($seq =~ /$sixA1/g); $num6A1 -= ($num7A1+$num8);
		$num5 = () = ($seq =~ /$fivemer/g); $num5 -= ($num6+$num6A1+$num7m8+$num7A1+$num8);

		print OUT join ("\t", $num8C, $num8G, $num8T, $num8, $num7m8, $num7A1, $num6, $numo6, $num6A1, $num5), "\n";
	}
	close FILE; close OUT;
}
close MIRSEQFILE;
