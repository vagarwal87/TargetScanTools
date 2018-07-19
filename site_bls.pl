#!/usr/bin/perl

use DBI;
use Bio::AlignIO;
use allfxns;

$home = $ENV{DIR};
$species = shift;
$region = shift;
$klen = shift;
$bin = shift;
$seqpos = shift;
$file = shift;

$code = (split /\./, (split /\//, $file)[-1])[0];
$klen--;

$refspec = "dm6";

if ($bin eq "all"){
	open TREE, "<$tree";
	$treestring = do { local($/); <TREE> }; #slurp up nh tree file
	close TREE;
} else{
	$treestring = `tail -1 ../Fig2/$species/$region/trees/bin$bin.REV.mod | tail -c+7`;
}

$mmnewick = parse_tree($treestring);

$in = Bio::AlignIO->new(-file => $file, -format => "fasta");
$aln = $in->next_aln;
$refseq = $aln->get_seq_by_id($refspec)->seq; $refseq =~ tr/-//d;
$len = length($refseq);
$diff = 1;

$i=$seqpos;
$pos = $aln->column_from_residue_number($refspec,$i);
$endpos = $aln->column_from_residue_number($refspec,$i+$klen);
$endloc = $aln->column_from_residue_number($refspec,$len);
$refbase = $aln->get_seq_by_id($refspec)->subseq($pos, $endpos);
$refbase =~ s/[^ACTGN]//g;

(%leaves, %totleaves, $numspec) = ((), (), 0);
foreach $seq ($aln->each_seq){
	$spec = $seq->id;
	$base = $seq->subseq(max($pos-10,1), min($endpos+10,$endloc));
	$base =~ s/[^ACTGN]//g;
	if ($base =~ $refbase){
		$leaves{$spec} = 1;
		$numspec+=1;
		$lastspec = $spec;
	}
}

$bls = 0;
$bls = sprintf("%.3f", BLS($mmnewick, \%leaves));
print "$bls\t$refbase\n";
