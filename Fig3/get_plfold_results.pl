#!/usr/bin/perl

use allfxns;
use Env; Env::import();

print "FC\tMIR_ID\tTYPE";
$dist = 35;
$winsize = 30;
for($i=-$dist; $i <= $dist+7; $i+=1){ for($j=1; $j <= $winsize; $j+=1){ print "\t$i:$j"; } }
print "\n";

open MIRSEQFILE, "<$MIRSEQFILE" or die "miRseq file not opening";
while(<MIRSEQFILE>){
	chomp; ($exp, $miR) = split /\t/;
	$seven8 = revCom(substr($miR, 1, 7));

	open FILE, "<$SITES/$exp\_flank.out" or die "can't open out file";
	while(<FILE>){ chomp $_;
		($utr_id, $fc_val, $pos, $seqlen, $flank5, $seed, $flank3, $site_type) = split;
		$maxsize = $dist+$winsize+1+7;
		next if $pos < $dist; #remove sites at edges to reduce impact of missing values
		next if $seqlen-$pos <= $maxsize; #remove sites at edges to reduce impact of missing values
		next if $site_type !~ /7A1|7m8|8m/;

		$pos += 1; #+1 for 0/1-based indexing correction
		@lines = split /\n/, `grep -A $maxsize -B $dist -P '^$pos\t' ./rna_probs_90/$utr_id\_lunp | cut -f 2-`;
		print join ("\t", $fc_val, $exp, $site_type), "\t";

		$x = 0; $y = -($dist*2 + 7); @a = ();
		foreach (@lines) {
			$_ =~ s/NA/1/g;
			@vals = split /\t/;
			@vals = @vals[max($y, 0)..min($x, $winsize-1)];
			$idx = max($x-($winsize-1),0);
			while(@vals){
				$a[$idx*($winsize-1)+$x] = pop(@vals);
				$idx++;
			}
			++$x; ++$y;
		}
		print join("\t", @a), "\n";
	}
	close FILE;
}
close MIRSEQFILE;
