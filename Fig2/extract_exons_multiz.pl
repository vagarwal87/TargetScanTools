#!/usr/bin/perl

# RUN FROM /tmp DIRECTORY

use DBI;
use IO::Uncompress::Gunzip;

$tmp = "/tmp";

$outfile  = pop @ARGV;

($specid, $uniqid, $chr, $str, $starts, $stops) = @ARGV;
@starts = split /,/, $starts;
@stops = split /,/, $stops;

$speclist = "dm6,droSim1,droSec1,droYak3,droEre2,droBia2,droSuz1,droAna3,droBip2,droEug2,droEle2,droKik2,droTak2,droRho2,droFic2,droPse3,droPer1,droMir2,droWil2,droVir3,droMoj3,droAlb1,droGri2,musDom2,anoGam1,apiMel4,triCas2";

for ($i = 0; $i <= $#starts; ++$i) {
	system "mafFrag dm6 multiz27way chr$chr $starts[$i] $stops[$i] $str $tmp/$uniqid.$i.maf";
	system "./maf2clustal.pl $tmp/$uniqid.$i.maf | modfasta.pl >$tmp/$uniqid.$i.fa";
	$files[$i] = "$tmp/$uniqid.$i.fa";
}

@files = reverse @files if $str eq '-';
$files = join " ", @files;

system "msa_view --aggregate $speclist $files >$outfile";
system "rm $tmp/$uniqid.*";
