#!/usr/bin/perl

use DBI;

$reftable = "fly6_07";
$table = "splicing_graph_$reftable";

$dbh = DBI->connect("dbi:mysql:database=username;host=HOSTNAME", "username", "password");

$exo = $dbh->prepare("select uniqid,exonids from $table");
$exo->execute();

while(($uniqid,$exonid) = $exo->fetchrow_array){
	$exonid = (split /,/, $exonid)[0];
	($region, $chr, $start, $stop, $str, $geneid) = $dbh->selectrow_array("select region,chr,start,stop,str,ids from $reftable where uniqid = \"$exonid\"");
	print "$uniqid\t$region\tchr$chr\t$start\t$stop\t$str\t$geneid\n";
}

$exo->finish;
$dbh->disconnect;
