#!/usr/bin/perl

use allfxns;

=pod
depends on the phastCons package

id#	 medianBLS bin
1000005 2.800707 4
1000010 3.621932 6
1000014 3.529802 6
1000015 1.319594 1

run with "./make_trees_from_bins.pl fly6_07/five_prime_UTR/allgenes.bins"
or with "./make_trees_from_bins.pl fly6_07/three_prime_UTR/allgenes.bins"

varying the binary codes below in order starting with (1,0,0)
=cut

$home = $ENV{DIR};
$model = "REV"; # JC69|F81|HKY85|HKY85+Gap|REV|SSREV|REV+GC|UNREST|R2|R2S|U2|U2S|R3|R3S|U3|U3S other types of substitution models
$region = "five_prime_UTR"; # five_prime_UTR three_prime_UTR
$species = "fly6_07";
$outdir = "$species";
($genfiles, $merge, $binned) = (1,0,0); #run with the following codes in order: 100, 011, 010+001, 000

$tree = "27way.nh";
$numbins = 5;
$speclist = "dm6,droSim1,droSec1,droYak3,droEre2,droBia2,droSuz1,droAna3,droBip2,droEug2,droEle2,droKik2,droTak2,droRho2,droFic2,droPse3,droPer1,droMir2,droWil2,droVir3,droMoj3,droAlb1,droGri2,musDom2,anoGam1,apiMel4,triCas2";
$toinclude = "dm6,droSim1,droSec1,droYak3,droEre2,droBia2,droSuz1,droAna3,droBip2,droEug2,droEle2,droKik2,droTak2,droRho2,droFic2,droPse3,droPer1,droMir2,droWil2,droVir3,droMoj3,droAlb1,droGri2,musDom2,anoGam1,apiMel4,triCas2";
$init_model = "fly_dm6.mod";

if ($genfiles){
	while(<>){
		($bin, $id) = (split)[-1,0];
		push (@{$ids[$bin]}, "./multifasta_$species/$region/$id.mfa");
	}
}

mkdir $outdir if ! -d "$outdir";
mkdir "$outdir/$region" if ! -d "$outdir/$region";
mkdir "$outdir/$region/trees" if ! -d "$outdir/$region/trees";
$outdir = "$outdir/$region/trees";

foreach $i (1..$numbins){
	print "Bin $i: $#{$ids[$i]} //$region\n";
	if ($genfiles){
		$k = 0;
		for ($j = 0; $j < @{$ids[$i]}; $j+=500){
			++$k;
			$files = (join " ", @{$ids[$i]}[$j..min($j+499, $#{$ids[$i]})]);
			system "bsub -e $home/tmp/sub.err -o $home/tmp/sub.out \"msa_view --out-format SS --aggregate $speclist --seqs $toinclude $files >$outdir/bin$i.$k.ss\""; #-T 3 .3tuple ==> tuple size
		}
	}
	elsif ($merge) {
		if ($binned){
			$files = `ls $outdir/bin$i.*.ss | tr '\n' ' '`;
			system "bsub -e $home/tmp/sub.err -o $home/tmp/sub.out \"msa_view --unordered-ss --in-format SS --out-format SS --aggregate $speclist --seqs $toinclude $files >$outdir/bin$i.ss; rm $files; \"";
		}
		else{
			$files = `ls $outdir/bin*.ss | tr '\n' ' '`;
			system "bsub -e $home/tmp/sub.err -o $home/tmp/sub.out \"msa_view --unordered-ss --in-format SS --out-format SS --aggregate $speclist --seqs $toinclude $files >$outdir/binall.ss; \"";
			last;
		}
	}
	else {
		if ($binned){
			system "bsub -e $home/tmp/sub.err -o $home/tmp/sub.out \"phyloFit -i SS --subst-mod $model --tree $tree -o $outdir/bin$i.$model $outdir/bin$i.ss\""; 
		}
		else{
			system "bsub -e $home/tmp/sub.err -o $home/tmp/sub.out \"phyloFit -i SS --scale-only --subst-mod $model --init-model $init_model -o $outdir/binall.$model $outdir/binall.ss\"";
			last;
		}
	}
}
