#!/usr/bin/perl

use Bio::AlignIO;

$home = $ENV{DIR};
$species = "fly6_07";
$file = shift;

($id = ($file =~ /(\d+).mfa/);

$refspec = "dm6";

$in = Bio::AlignIO->new(-file => $file, -format => "fasta");
$aln = $in->next_aln;
$refseq = $aln->get_seq_by_id($refspec)->seq; $refseq =~ tr/-//d;

print ">$id\n$refseq\n";
