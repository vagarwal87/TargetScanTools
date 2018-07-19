#!/usr/bin/perl

use allfxns;
use Env; Env::import();

$dist = 30;

%utrBase = %{readFasta($SEQBASE)};
%utr5Base = %{readFasta2($SEQBASE5UTR)};
%orfBase = %{readFasta2($SEQBASEORF)};

%gene2bin = ();
($species, $region) = ("fly_90%", "three_prime_UTR");

sub readFasta2{
	local $fasta = shift;
	local %fasta = ();
	open DNA, "<$fasta" || die "Could not open fasta file for $fasta\n";
	while ($line = <DNA>){ chomp $line;
		if ($line =~ /^>\s?(\w+)\.?/){ $header = $1; }
		else { $fasta{$header} .= $line; }
	}
	close DNA;
	return \%fasta;
}

###THESE DATA RELY ON GENOME-WIDE ALIGNMENTS THAT ARE NOT PROVIDED WITH THIS PACKAGE..must extract multiple sequence alignments yourself
###comment out code block below, and remove need to compute branch length score later in this script to proceed without branch length scores
open IN, "<../Fig2/$species/$region/allgenes.bins" || die "Could not open fasta file for $fasta\n";
while (<IN>){ chomp;
	@a = split /\t/;
	$gene2bin{$a[0]} = $a[-1];
}
close IN;

print join ("\t", 'FC', 'GENE_ID', 'MIR_ID', 'TYPE',
	'TARNUC10', 'TARNUC9', 'TARNUC8', 'TARNUC1', 'UORF', 'UTR5P8', 'ORF8', 'UTR6', 'UTR6A1', 'UTRO6',
	'SUPP_PAIR', 'LEN_3UTR', 'LEN_5UTR', 'LEN_ORF', 'DIST5P', 'DIST3P', 'MIN_DIST',
	'LOCALAU_30', 'AU_3UTR', 'AU_5UTR', 'AU_ORF', 'PLFOLD', 'HYBRIDSCORE', 'M8SEED', 'BLS', 'SITE', 'MIR'), "\n";

open MIRSEQFILE, "<$MIRSEQFILE" or die "miRseq file not opening";
while(<MIRSEQFILE>){ chomp($_);
	($exp, $miR) = split /\t/;
	$seven8 = revCom(substr($miR, 1, 7));
	open FILE, "<$SITES/$exp\_flank.out" or die "can't open SITES file";
	$count = 1;

	$mirnuc10 = substr($miR, 9, 1);
	$mirnuc9 = substr($miR, 8, 1);
	$mirnuc8 = substr($miR, 7, 1);
	$mirnuc1 = substr($miR, 0, 1);

	$eightmer = revCom(substr($miR, 1, 7)).'A';
	$seven8 = substr($eightmer,0,7);
	$sevenA = substr($eightmer,1,7);
	$sixmer = substr($eightmer,1,6);
	$offset6 = substr($eightmer,0,6);
	$sixA1 = substr($eightmer,2,6);

	while(<FILE>){
		chomp($_); ++$count;
		($utr_id, $fc_val, $pos, $len_3utr, $flank5, $seed, $flank3, $site_type) = split /\t/;

		$seq = $utrBase{$utr_id}; $utr_id1 = (split /\./, $utr_id)[0];
		$seqorf = $orfBase{$utr_id1}; $orflen = length($seqorf);
		$seq5utr = $utr5Base{$utr_id1}; $len_5utr = length($seq5utr);
		die "$utr_id not found" if $seq eq '';

		($dist5p, $dist3p, $min_dist, $target)
			= ($pos+1, $len_3utr-($pos+5), min($pos+1, $len_3utr-($pos+5)), $flank5.$seed);

		$tmp=(split /\t/, `grep -A13 -P '^$dist5p\t' $PLFOLD/$utr_id\_lunp | tail -1 | cut -f 2-`)[24];
		$tmp=0.01 if $tmp eq 'NA';
		$plfold = sprintf("%.3f", log10($tmp)); #take probability of unpairing

		die "ERROR $utr_id" if $plfold eq '';

		$numORF8 = () = ($seqorf =~ /$eightmer/g);
		$numUTR5P8 = () = ($seq5utr =~ /$eightmer/g);
		$num_uORFs  = () = ($seq5utr =~ /ATG/g);

		$num8 = () = ($seq =~ /$eightmer/g);
		$num7m8 = () = ($seq =~ /$seven8/g); $num7m8 -= $num8;
		$num7A1 = () = ($seq =~ /$sevenA/g); $num7A1 -= $num8;
		$num6  = () = ($seq =~ /$sixmer/g);  $num6 -= ($num7m8+$num7A1+$num8);
		$numo6 = () = ($seq =~ /$offset6/g); $numo6 -= ($num7m8+$num8);
		$num6A1 = () = ($seq =~ /$sixA1/g); $num6A1 -= ($num7A1+$num8);

		$tarnuc10 = substr($flank5, -2, 1);
		$tarnuc9 = substr($flank5, -1, 1);
		$tarnuc8 = substr($seed, 0, 1);
		$tarnuc1 = substr($seed, -1, 1);

		$localau = localAU($flank5, $seed, $flank3, $site_type, 30);

		$utr3AU = 1-gcContent($seq);
		$orfAU = 1-gcContent($seqorf);
		$utr5AU = 1-gcContent($seq5utr);

		### comment this out to run script without branch length scores, need to compute multiple sequence alignments otherwise for it to work
		@mybls = split / /, site_bls($seven8, $site_type, $gene2bin{$utr_id}, $pos, "../Fig2/multifasta_$species/$region/$utr_id.mfa");

		$score3p = contextScore($site_type, $localau, $pos, $min_dist, $ta, $sps, $miR, "A".$seq); #add A so ribosome shadow error never occurs
		$rnascore = sprintf("%.2f", rnahybrid(rev(rev($miR)), rev(rev($target))));

		print join ("\t", $fc_val, $utr_id, $exp, $site_type,
			$tarnuc10, $tarnuc9, $tarnuc8, $tarnuc1, $num_uORFs, $numUTR5P8, $numORF8, $num6, $num6A1, $numo6, sprintf("%.2f", $score3p), $len_3utr, $len_5utr, $orflen, $dist5p, $dist3p, $min_dist,  #, $numUTR5P7m8 , $sps, $ta, $ta_orf
			$localau, $utr3AU, $utr5AU, $orfAU, $plfold, $rnascore, $seven8, @mybls, $miR), "\n" if $pos > 14 && index($seven8.'A', $mybls[1]) != -1	; # if $score3p ne ""; # $rnascore\t$pval\t , , $match2 $ddG, sprintf("%.2f", $weightedAU), $dGduplex, $dGopen, $local2bp , $weight 
	}
}
close MIRSEQFILE;

sub contextScore{
	local ($site_type, $localAU, $pos, $min_dist, $ta, $sps, $miR, $seq) = @_;
	local ($contextscore, $score3p, $contextPlus);
	$pos += 2 if $site_type eq "6m" || $site_type eq "7A1";
	$pos += 1 if $site_type eq "7m8" || $site_type eq "8m" || $site_type eq "o6";
	$pos += 3 if $site_type eq "6A1";
	$site_type = "8mer" if $site_type eq "8m";
	$site_type = "6mer" if $site_type eq "6m" || $site_type eq "o6" || $site_type eq "6A1";
	@score = split /\n/, `./ComputeContextScore $site_type $pos $miR $seq | tail -4`;
	$score3p = (split /\s/, $score[2])[-1];
	return $score3p;
}

sub localAU{
	local ($utrUp, $seed, $utrDown, $siteType, $distToCheck) = @_;
	local ($score, $maxAU, $totalUp, $totalDown) = (0,0,0,0);
	$utrUp = substr($utrUp,-$distToCheck) if length($utrUp) > $distToCheck;	#use only 30 bp, consistent with Grimson et al 2007 metric
	$utrDown = substr($utrDown,0,$distToCheck);
	if ($siteType eq '7m8' || $siteType eq '6m' || $siteType eq 'o6'){ $utrDown = substr($seed,-1).$utrDown; }
	@utrUp3to5 = split(//, reverse (uc ($utrUp)));
	@utrDown5to3 = split(//, uc ($utrDown));
	if ($siteType eq '8m'){
		for (my $i = 0; $i <= $#utrUp3to5; $i++){
			$score = 1 / ($i + 1);
			if ($utrUp3to5[$i] =~ /U|T|A/){ $totalUp += $score; }
			$maxAU += $score;
		}
		for ($i = 0; $i <= $#utrDown5to3; $i++){
			$score = 1 / ($i + 2);
			if ($utrDown5to3[$i] =~ /U|T|A/){ $totalDown += $score; }
			$maxAU += $score;
		}
	}
	elsif ($siteType eq '7m8'){
		for ($i = 0; $i <= $#utrUp3to5; $i++){
			$score = 1 / ($i + 1);
			if ($utrUp3to5[$i] =~ /U|T|A/){ $totalUp += $score; }
			$maxAU += $score;
		}
		for ($i = 0; $i <= $#utrDown5to3; $i++){
			$score = 1 / ($i + 1);
			if ($i == 0) { $score = 1 / 2; }
			if ($utrDown5to3[$i] =~ /U|T|A/){ $totalDown += $score; }
			$maxAU += $score;
		}
	}
	elsif ($siteType eq '7A1'){
		for ($i = 0; $i <= $#utrUp3to5; $i++){
			$score = 1 / ($i + 2);
			if ($utrUp3to5[$i] =~ /U|T|A/){ $totalUp += $score; }
			$maxAU += $score;
		}
		for ($i = 0; $i <= $#utrDown5to3; $i++){
			$score = 1 / ($i + 2);
			if ($utrDown5to3[$i] =~ /U|T|A/){ $totalDown += $score; }
			$maxAU += $score;
		}
	}
	elsif ($siteType eq '6m'){
		for ($i = 0; $i <= $#utrUp3to5; $i++){
			$score = 1 / ($i + 2);
			if ($utrUp3to5[$i] =~ /U|T|A/){ $totalUp += $score; }
			$maxAU += $score;
		}
		for ($i = 0; $i <= $#utrDown5to3; $i++){
			$score = 1 / ($i + 1);
			if ($i == 0) { $score = 1 / 2; }
			if ($utrDown5to3[$i] =~ /U|T|A/){ $totalDown += $score; }
			$maxAU += $score;
		}
	}
	elsif ($siteType eq 'o6'){
		for ($i = 0; $i <= $#utrUp3to5; $i++){ #up flank is like 8mer or 7m8
			$score = 1 / ($i + 1);
			if ($utrUp3to5[$i] =~ /U|T|A/){ $totalUp += $score; }
			$maxAU += $score;
		}
		for ($i = 0; $i <= $#utrDown5to3; $i++){ #down flank is like 7A1, treating 1A position as any other A/U pair
			$score = 1 / ($i + 2);
			if ($utrDown5to3[$i] =~ /U|T|A/){ $totalDown += $score; }
			$maxAU += $score;
		}
	}
	elsif ($siteType eq '6A1'){
		for ($i = 0; $i <= $#utrUp3to5; $i++){ #up flank is like 8mer or 7m8
			$score = 1 / ($i + 3);
			if ($utrUp3to5[$i] =~ /U|T|A/){ $totalUp += $score; }
			$maxAU += $score;
		}
		for ($i = 0; $i <= $#utrDown5to3; $i++){ #down flank is like 7A1, treating 1A position as any other A/U pair
			$score = 1 / ($i + 2);
			if ($utrDown5to3[$i] =~ /U|T|A/){ $totalDown += $score; }
			$maxAU += $score;
		}
	}
	else { die "NOT A SITE TYPE\n"; }
	return sprintf('%0.3f', ($totalUp + $totalDown)/$maxAU );
}

sub rev{
	local $seq = shift;
	$seq =~ tr/T/U/;
	$seq = reverse $seq;
	return $seq;
}

sub rnahybrid{
	local ($miR, $target) = @_;
	$target = substr($target,max(0,length($target)-17)); #use upstream length of 9nt as optimal 9+8=17
	$miR = substr($miR,0,17);
	$tmp = `RNAhybrid -c -s 3utr_fly -f 3,6 $target $miR`;
	($p_score, $pval) = (split /:/, $tmp)[4,5];
	$miR = substr($miR,0,12);
	$tmp = `RNAhybrid -c -s 3utr_fly -f 3,6 $target $miR`;
	($p_score2) = (split /:/, $tmp)[4];
	return $p_score - $p_score2;
}

sub site_bls{
	local ($seven8, $kmerlen, $bin, $pos, $file) = @_;
	$pos += 2 if $site_type eq "6m" || $site_type eq "7A1";
	$pos += 1 if $site_type eq "7m8" || $site_type eq "8m" || $site_type eq "o6";
	$pos += 3 if $site_type eq "6A1";
	$kmerlen = 6 if $site_type eq "6m" || $site_type eq "6A1" || $site_type eq "o6";
	$kmerlen = 7 if $site_type eq "7m8" || $site_type eq "7A1";
	$kmerlen = 8 if $site_type eq "8m";

	local $result = `perl site_bls.pl fly6_07 $region $kmerlen $bin $pos $file 2>&-`;
	chomp $result;
	return ($result ne '')? $result : "0\t$seven8";
}
