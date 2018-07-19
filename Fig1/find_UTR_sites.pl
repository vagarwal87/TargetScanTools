#!/usr/bin/perl

use allfxns;
use Env; Env::import();

%utrBase = %{readFasta($SEQBASE)};
`mkdir $SITES/` unless (-d "$SITES");
$file = shift;

open MIRSEQFILE, "<$MIRSEQFILE" or die "miRseq file not opening";
while(<MIRSEQFILE>){
	chomp($_); ($miRname, $miR) = split /\t/, $_;

	$exp = (split /\//, $file)[-1];
	$exp = (split /\./, $exp)[0];

	$seed = revCom(substr($miR, 1, 7)).'A';
	$seven8 = substr($seed,0,7);
	$sevenA = substr($seed,1,7);

	open FILE, "<$file";
	open OUT, ">$SITES/$exp\_$miRname\_sites.out" or die "can't open out file";

	print OUT join ("\t", "utr_id", "fc", "utr8", "utr7m8", "utr7A1"), "\n";

	while(<FILE>){
		chop($_); ++$count;
		($utr_id, $expr_val) = split /\t/;
		$seq = $utrBase{$utr_id}; next if $seq eq '';

		$num8 = () = ($seq =~ /$seed/g);
		$num7m8 = () = ($seq =~ /$seven8/g); $num7m8 -= $num8;
		$num7A1 = () = ($seq =~ /$sevenA/g); $num7A1 -= $num8;

		print OUT join ("\t", $utr_id, sprintf("%.3f", $expr_val), $num8, $num7m8, $num7A1), "\n";
	}
	close FILE; close OUT;
}
close MIRSEQFILE;
