#!/usr/bin/perl

use allfxns;

%gene2bin = ();

$region = shift;
$dir = "targetsites";
`mkdir $dir` unless (-d $dir);
`mkdir "$dir/$region"` unless (-d "$dir/$region");

$species = "fly6_07";

%utrBase = ($region eq 'three_prime_UTR')? %{readFasta("utr3_seqs.fa")} : %{readFasta("utr5_seqs.fa")};

open IN, "<$species/$region/allgenes.bins" || die "Could not open fasta file for $fasta\n";
while (<IN>){ chomp;
	@a = split /\t/;
	$gene2bin{$a[0]} = $a[-1];
}
close IN;

open MIRSEQFILE, "<cons_seeds_TS7.txt" or die "miRseq file not opening";

while(<MIRSEQFILE>){ chomp($_);
	($exp, $miR) = split /\t/, $_;
	$seed = revCom(substr($miR, 1, 7)).'A';
	$seven8 = substr($seed,0,7);
	$sevenA = substr($seed,1,7);
	$sixmer = substr($seed,1,6);
	$offset6 = substr($seed,0,6);
	$sixA1 = substr($seed,2,6);

	$outfile = "$dir/$region/$exp\_flank.out";
	$infile = ($region eq 'three_prime_UTR')? "utr3_ids.txt" : "utr5_ids.txt";
	open FILE, "<$infile";
	open OUT, ">$outfile" or die "can't open out file";
	while(<FILE>){ chomp($_);
		$utr_id = $_;
		$seq = $utrBase{$utr_id};
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
			if    (   $pos8 == $minpos ){ $site_type = '8m';  $pos = $pos8;     }
			elsif ( $pos7m8 == $minpos ){ $site_type = '7m8'; $pos = $pos7m8;   }
			elsif ( $pos7a1 == $minpos ){ $site_type = '7A1'; $pos = $pos7a1-1; }
			elsif (  $poso6 == $minpos ){ $site_type = 'o6';  $pos = $poso6;    }
			elsif (   $pos6 == $minpos ){ $site_type = '6m';  $pos = $pos6-1;   }
			else 				    { $site_type = '6A1'; $pos = $pos6a1-2; }

			#retain only sites enriched above background at BLS = 1.0 for each region
			if (($region eq "three_prime_UTR" && $site_type =~ /8m|7m8|7A1|o6|6m/) || ($region eq "five_prime_UTR" && $site_type =~ /8m|7m8|7A1/)){
				@mybls = split / /, site_bls($seven8, $site_type, $gene2bin{$utr_id}, $pos, "/lab/bartel3_ata/agarwal/databases/multifasta_$species/$region/$utr_id.mfa");
				print OUT join ("\t", $utr_id, $pos, $site_type, @mybls) if $mybls[0] >= 1.0; #print if bls exceeds 1.0 (where signal - background is maximized)
			}

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

## if 5' UTR --> for each 8mer find 7mer subsequence with maximal BLS
## if 3' UTR --> for each 7mer or 8mer find 6mer subsequence with maximal BLS
sub site_bls{
	local ($seven8, $site_type, $bin, $pos, $file) = @_;
	local @vals1, @vals2 = (), ();
	$pos += 1;
	if ($region eq "three_prime_UTR"){
		$kmerlen = 6;
		@vals1 = split / /, `perl site_bls.pl $region $kmerlen $bin $pos $file 2>&-`;
		$pos += 1;
		@vals2 = split / /, `perl site_bls.pl $region $kmerlen $bin $pos $file 2>&-`;
		$result = join("\t", @vals1) if $site_type eq "o6";
		$result = join("\t", @vals2) if $site_type =~ "6m|7A1";
		if ($site_type =~ "7m8|8m"){ $result = ( max($vals1[0], $vals2[0]) == $vals1[0] ) ? join("\t", @vals1) : join("\t", @vals2); }
	}
	else{
		$kmerlen = 7;
		@vals1 = split / /, `perl site_bls.pl $region $kmerlen $bin $pos $file 2>&-`;
		$pos += 1;
		@vals2 = split / /, `perl site_bls.pl $region $kmerlen $bin $pos $file 2>&-`;
		$result = join("\t", @vals1) if $site_type eq "7m8";
		$result = join("\t", @vals2) if $site_type eq "7A1";
		if ($site_type eq "8m"){ $result = ( max($vals1[0], $vals2[0]) == $vals1[0] ) ? join("\t", @vals1) : join("\t", @vals2); }
	}
	return ($result ne '')? $result : "0\t$seven8\n";
}
