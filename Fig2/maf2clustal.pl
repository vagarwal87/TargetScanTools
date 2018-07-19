#!/usr/bin/perl

use Bio::AlignIO;
use Bio::SimpleAlign;

$file = shift;

$in  = Bio::AlignIO->newFh(-format => "maf", -file => $file);
$out = Bio::AlignIO->newFh(-format => "fasta");

print $out $_ while (<$in>);
