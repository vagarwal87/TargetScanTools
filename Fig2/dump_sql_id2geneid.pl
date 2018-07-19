#!/usr/bin/perl

use DBI;

$reftable = "fly6_07";
$table = "splicing_graph_$reftable";

$dbh = DBI->connect("dbi:mysql:database=username;host=HOSTNAME", "username", "password");
$ids = $dbh->prepare("select ids from $reftable where uniqid = ?");
$exo = $dbh->prepare("select uniqid,region,exonids from $table order by uniqid"); #get all exons
$exo->execute();

while(($uniqid, $region, $exonids) = $exo->fetchrow_array){
	$exonids = (split /,/, $exonids)[0];
	$ids->execute($exonids);
	$id = $ids->fetchrow_array;
	print "$uniqid\t$region\t$id\n";
}

$ids->finish;
$exo->finish;
$dbh->disconnect;
