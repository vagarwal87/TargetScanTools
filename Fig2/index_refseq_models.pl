#!/usr/bin/perl

use DBI;

# can create an index for faster querying afterwards: i.e., "create index fly6_07geneindex on fly6_07All (uniqid);"

$HOME = $ENV{DIR};
$count = 100000;

$table = "fly6_07";

$file = "$table.gtf";

open FILE, "sort -k 1,1 -k 9,9 -k 4,4n $file |";
while(<FILE>){ chomp;
	($chr, $start, $stop, $str, $id) = (split)[0,3,4,6,-1];
	$uniqid{$chr.$str.$start.$stop.$id} = $count;
	$line{$chr.$str.$start.$stop.$id} = $_;
	++$count;
}
close FILE;

open FILE, "sort -k 1,1 -k 9,9 -k 4,4n $file |";
while(<FILE>){
	($chr, $start, $stop, $str, $id) = (split)[0,3,4,6,-1]; $id =~ s/"//g;
	$this = $uniqid{$chr.$str.$start.$stop.$id};
	$exonids{$this} .= "$id," if $this > 0;		#see which ids use this exon in their model
	if ($id eq $previd){
		$prev = $uniqid{$prevchr.$prevstr.$prevstart.$prevstop.$id};
		if ($prev > 0 && $this > 0){
			$nextlink{$prev}{$this} = 1;
			$prevlink{$this}{$prev} = 1;
		}
	}
	($previd, $prevchr, $prevstart, $prevstop, $prevstr) = ($id, $chr, $start, $stop, $str);
}
close FILE;

$dbh = DBI->connect("dbi:mysql:database=username;host=HOSTNAME", "username", "password");

$stmt = $dbh->do("create table $table (
	uniqid integer unsigned,
	chr varchar(10),
	region varchar(15),
	start integer unsigned,
	stop integer unsigned,
	str char(1),
	frame char(1),
	ids varchar(255),
	prevexons varchar(255),
	nextexons varchar(255))");

$sth = $dbh->prepare("insert into $table values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");

foreach $pos (keys %uniqid){
	$val = $uniqid{$pos};
	($chr, $region, $start, $stop, $str, $frame) = (split /\s/, $line{$pos})[0,2,3,4,6,7];
	($next, $prev) = ('',''); $chr =~ s/chr//;
	foreach $i (keys %{ $nextlink{$val} }){ $next .= "$i,"; }
	foreach $i (keys %{ $prevlink{$val} }){ $prev .= "$i,"; }
	chop $next; chop $prev; chop $exonids{$val};
	$sth->execute($val, $chr, $region, $start, $stop, $str, $frame, $exonids{$val}, $prev, $next);
}

$sth->finish;
$dbh->disconnect;
