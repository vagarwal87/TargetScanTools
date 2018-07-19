#!/usr/bin/perl

use DBI;
use allfxns;

$home = $ENV{DIR};
$tmp = "$home/tmp";
$suberr = "$tmp/sub.err";
$subout = "$tmp/sub.out";
$subopts = ''; #"-q short" queue for job submission, 
$region = "five_prime_UTR"; #three_prime_UTR or five_prime_UTR
$specid = "fly6_07";

$dbh = DBI->connect("dbi:mysql:database=username;host=HOSTNAME", "username", "password");
$exo = $dbh->prepare("select uniqid,chr,str,exonids from splicing_graph_$specid where region like \"$region\""); #get all regions
$exo->execute();

$loc = $dbh->prepare("select start,stop from $specid where uniqid = ?"); #get all exons

$jobnum = 10000;
$region =~ s/'//g;

while(($uniqid, $chr, $str, $exonids) = $exo->fetchrow_array){
	$outfile = "multifasta_$specid/$region/$uniqid.mfa";
	mkdir "multifasta_$specid/" if ! -d "multifasta_$specid/";
	mkdir "multifasta_$specid/$region" if ! -d "multifasta_$specid/$region";
	next if (-s $outfile != 0); # next if filesize == 0

	@exons = split /,/, $exonids;
	($starts, $stops, $doit) = ("", "", 1);
	for ($i = 0; $i <= $#exons; ++$i){
		$doit = 0 if ($badids{$exons[$i]} == 1); #skip region if one of the exons in transcript maps redundantly
		$loc->execute($exons[$i]);
		while(($start, $stop) = $loc->fetchrow_array) {
			$start--;
			$starts .= $start.',';
			$stops  .= $stop.',';
		}
	}
	print "$uniqid SKIPPED\n" if !$doit;
	next if !$doit;

	chop $starts; chop $stops;

	$jobfile = "$tmp/job$jobnum.sh";
	open SH, ">>$jobfile" or die "can't open $jobfile";
	print SH "extract_exons_multiz.pl $specid $uniqid $chr $str $starts $stops $outfile;\n";
	close SH;

	$count++;

	if ($count % 100 == 0) {
		bsub($subopts, $suberr, $subout, $jobfile);
		$jobnum++;
	}
}

bsub($subopts, $suberr, $subout, $jobfile);

$loc->finish;
$exo->finish;
$dbh->disconnect;
