#!/usr/bin/perl

use DBI;

$reftable = "fly6_07";
$table = "splicing_graph_$reftable";
$count = 1000000;

$dbh = DBI->connect("dbi:mysql:database=username;host=HOSTNAME", "username", "password");

$stmt = $dbh->do("create table $table (
	uniqid integer unsigned,
	chr varchar(10),
	region varchar(15),
	str char(1),
	exonids text)");

$sth = $dbh->prepare("insert into $table values (?, ?, ?, ?, ?, ?)");
$ids = $dbh->prepare("select uniqid from $reftable where ids = ? and region = ? order by uniqid");
$exo = $dbh->prepare("select uniqid,chr,region,str,ids from $reftable order by uniqid"); #get all exons
$exo->execute();

while(($id, $chr, $region, $str, $geneid) = $exo->fetchrow_array){
	next if exists($seenid{$geneid.$region});
	$seenid{$id} = 1;
	$chain_exons = chain_exons($geneid, $region);
	$sth->execute($count, $chr, $region, $str, $chain_exons);
	++$count;
}

$sth->finish;
$ids->finish;
$exo->finish;
$dbh->disconnect;

sub chain_exons{
	local ($geneid, $region) = @_;
	$ids->execute($geneid, $region);
	$seenid{$geneid.$region} = 1;
	return join ",", map { @$_ } @{ $ids->fetchall_arrayref() };
}
