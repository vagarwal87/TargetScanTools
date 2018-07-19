#!/usr/bin/perl

use allfxns;
use Env; Env::import();

$mylen = shift;

print "FC\tMIR_ID\tTYPE";
$winsize = 20-8+1;
for($i = 9; $i <= 21; $i+=1){ for($j=1; $j <= $winsize; $j+=1){ print "\t$i:$j" if $i+$j <= 22; } }
print "\n";

open MIRSEQFILE, "<$MIRSEQFILE" or die "miRseq file not opening";
while(<MIRSEQFILE>){
	chomp; ($exp, $miR) = split /\t/;
	$seven8 = revCom(substr($miR, 1, 7));
	$miR = rev(rev($miR));
	open FILE, "<$SITES/$exp\_flank.out" or die "can't open out file";
	while(<FILE>){ chomp $_;
		($utr_id, $fc_val, $pos, $seqlen, $flank5, $seed, $flank3, $site_type) = split;
		next if $site_type !~ /7A1|7m8|8m/;

		$target = rev(rev($flank5.$seed));
		$target = substr($target,max(0,length($target)-($mylen+8)));

		print join ("\t", $fc_val, $exp, $site_type);
		$x = 0; @a = ();
		for($i = 8; $i <= 20; $i+=1){
			for($j=1; $j <= $winsize; $j+=1){
				$miR1 = substr($miR,0,$i);
				$tmp = `RNAhybrid -c -s 3utr_fly -f 3,6 $target $miR1`;
				$p_score1 = (split /:/, $tmp)[4];
				if($i+$j <= 21){
					$miR2 = substr($miR,0,$i+$j);
					$tmp = `RNAhybrid -c -s 3utr_fly -f 3,6 $target $miR2`;
					$p_score2 = (split /:/, $tmp)[4];
					print "\t", sprintf("%0.2f", $p_score2 - $p_score1);
				}
			}
		}
		print "\n";
	}
	close FILE;
}
close MIRSEQFILE;
