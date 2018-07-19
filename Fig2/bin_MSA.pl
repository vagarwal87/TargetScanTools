#!/usr/bin/perl

#in first try, use with getmedian to determine the appropriate bin for each gene, then use bins later for word counts

use DBI;
use Bio::AlignIO;
use allfxns;

$home = $ENV{DIR};
$species = "fly6_07";
$region = shift;
$kmerlen = shift;
$bin = shift;
$file = shift;

$code = `basename $file .mfa`; chomp $code;

$getmedian = ($bin eq "all" && $kmerlen == 1) ? 1 : 0;
if ($getmedian){ $kmerlen = 1-1; } #use single nucleotide as the unit for median bls calculation
else { $kmerlen--; }

$refspec = "dm6";
$tree = "27way.nh";

if ($bin eq "all"){
	open TREE, "<$tree";
	$treestring = do { local($/); <TREE> }; #slurp up nh tree file
	close TREE;
} else{
	$treestring = `tail -1 $species/$region/trees/bin$bin.REV.mod | tail -c+7`;
}

$mmnewick = parse_tree($treestring);
#$max = getMaxBLS($mmnewick);

$in = Bio::AlignIO->new(-file => $file, -format => "fasta");
$aln = $in->next_aln;
$refseq = $aln->get_seq_by_id($refspec)->seq; $refseq =~ tr/-//d;
$len = length($refseq);
$diff = 1;

for($i=1; $i<=$len; $i+=$diff){
	$pos = $aln->column_from_residue_number($refspec,$i);
	for $klen (max(0,$kmerlen-2)..$kmerlen){ # get 6-8mers or 1-2mers
		next if ($i+$klen > $len) || ($getmedian && $klen != $kmerlen);
		$endpos = $aln->column_from_residue_number($refspec,$i+$klen);
		$base = $aln->get_seq_by_id($refspec)->subseq($pos, $endpos);
		($count, %leaves, %totleaves) = (0, (), ());
		foreach $seq ($aln->each_seq){
			$spec = $seq->id;
			$base = $seq->subseq($pos, $endpos);
			$refbase = $base if ($count == 0);
			$leaves{$spec} = 1 if $base eq $refbase;
			$count++;
		}
		$bls = sprintf("%.3f", BLS($mmnewick, \%leaves));
		push(@blsvals, $bls);
		$refbase =~ tr/-//d;
		print "$refbase $bls\n" if !$getmedian;
	}
}

print join("\t", $code, length($refseq), gcContent($refseq), sprintf("%.3f", median(\@blsvals))), "\n" if $getmedian;
