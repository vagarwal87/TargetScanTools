#!/usr/bin/perl

use allfxns;
use Env; Env::import();

$dist = 30;

%utrBase = %{readFasta($SEQBASE)};
%utr5Base = %{readFasta2($SEQBASE5UTR)};
%orfBase = %{readFasta2($SEQBASEORF)};

%gene2bin = ();
($species, $region) = ("fly6_07", "three_prime_UTR");

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

open IN, "<../Fig2/$species/$region/allgenes.bins" || die "Could not open fasta file for $fasta\n";
while (<IN>){ chomp;
	@a = split /\t/;
	$gene2bin{$a[0]} = $a[-1];
}
close IN;

open IN, "<id2countID.txt" || die "Could not open\n";
while (<IN>){ chomp;
	@a = split /\t/;
	$flyid2countid{$a[0]} = $a[-1];
}
close IN;

print join ("\t", 'FC', 'GENE_ID', 'MIR_ID', 'TYPE', 'UTR5P8', 'ORF8', 'UTR6', 'UTR6A1', 'UTRO6', 'LEN_3UTR', 'LEN_ORF', "HYBRIDSCORE",
	'PLFOLD', 'M8SEED', 'BLS', 'SITE', 'MIR'), "\n";

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
		($utr_id, $FC_val, $pos, $len_3utr, $flank5, $seed, $flank3, $site_type) = split /\t/;

		$seq = $utrBase{$utr_id}; $utr_id1 = (split /\./, $utr_id)[0];
		$seqorf = $orfBase{$utr_id1}; $orflen = length($seqorf);
		$seq5utr = $utr5Base{$utr_id1};

		die "$utr_id not found" if $seq eq '';

		($dist5p, $target) = ($pos+1, $flank5.$seed);

		$plfold = sprintf("%.3f", log10((split /\t/, `grep -A13 -P '^$dist5p\t' $PLFOLD/$utr_id\_lunp | tail -1 | cut -f 2-`)[24])); #take probability of unpairing
		die "ERROR $utr_id" if $plfold eq '';

		$numORF8 = () = ($seqorf =~ /$eightmer/g);
		$numUTR5P8 = () = ($seq5utr =~ /$eightmer/g);

		$num8 = () = ($seq =~ /$eightmer/g);
		$num7m8 = () = ($seq =~ /$seven8/g); $num7m8 -= $num8;
		$num7A1 = () = ($seq =~ /$sevenA/g); $num7A1 -= $num8;
		$num6  = () = ($seq =~ /$sixmer/g);  $num6 -= ($num7m8+$num7A1+$num8);
		$numo6 = () = ($seq =~ /$offset6/g); $numo6 -= ($num7m8+$num8);
		$num6A1 = () = ($seq =~ /$sixA1/g); $num6A1 -= ($num7A1+$num8);

		$consid=$flyid2countid{$utr_id};

		#branch length score calculation requires running pipeline to extract 3' UTR alignments in Fig2
		@mybls = split / /, site_bls($seven8, $site_type, $gene2bin{$consid}, $pos, "/lab/bartel3_ata/agarwal/databases/multifasta_$species/$region/$consid.mfa"); #../Fig2/

		($rnascore) = sprintf("%.2f", rnahybrid(rev(rev($miR)), rev(rev($target))));

		print join ("\t", $FC_val, $utr_id, $exp, $site_type, $numUTR5P8, $numORF8, $num6, $num6A1, $numo6,
		$len_3utr, $orflen, $rnascore, $plfold, $seven8, @mybls, $miR), "\n" if $pos > 14 && index($seven8.'A', $mybls[1]) != -1;
	}
}
close MIRSEQFILE;

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

	local $result = `perl site_bls.pl $species $region $kmerlen $bin $pos $file 2>&-`;
	chomp $result;
	return ($result ne '')? $result : "0\t$seven8";
}
