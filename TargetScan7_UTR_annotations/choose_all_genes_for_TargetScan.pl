#!/usr/bin/perl

use allfxns;

$file = "dmel-all-no-analysis-r6.19.gff.gz";

#map transcript IDs to gene IDs
open IN, "zgrep -P '\tmRNA\t' $file | sort -k1,1 -k4,4n | ";
while(<IN>){
	($id, $parent) = ($_ =~ /ID=(FBtr\d+).*Parent=(FBgn\d+)/);
	$id2parent{$id}=$parent;
}
close IN;

#find most proximal 3'UTR positions for each transcript, most distal 3' UTR position for each gene, and lengths of all annotated 3'UTRs
open IN, "zgrep -P '\tthree_prime_UTR\t' $file | sort -k1,1 -k4,4n | ";
while(<IN>){ chomp;
	($chr, $region, $start, $stop, $str, $last) = (split /\t/)[0,2,3,4,6,-1];
	($id) = ($last =~ /Parent=(FBtr.*)/);
	@ids = split /,/, $id;
	foreach $id (@ids){
		if($str eq '+'){
			if (! exists $utrEndpoint{$id2parent{$id}}{$chr} || $stop > $utrEndpoint{$id2parent{$id}}{$chr}){
				$utrEndpoint{$id2parent{$id}}{$chr} = $stop;
			}
			if (! exists $utrStartpoint{$id} || $start < $utrStartpoint{$id}){
				$utrStartpoint{$id} = $start;
			}
		}
		else{
			if (! exists $utrEndpoint{$id2parent{$id}}{$chr} || $start < $utrEndpoint{$id2parent{$id}}{$chr}){
				$utrEndpoint{$id2parent{$id}}{$chr} = $start;
			}
			if (! exists $utrStartpoint{$id} || $stop > $utrStartpoint{$id}){
				$utrStartpoint{$id} = $stop;
			}
		}
		$lengths{$id}{$region} += ($stop-$start);
		if(defined($prevpos{$id})){ $intronlen3pUTR{$id}+=($start-$prevpos{$id}-1); }
		$prevpos{$id}=$stop;
	}
}
close IN;

open IN, "zcat 3Pseq_pA_allstages.summed.gz | ";
while(<IN>){
	($chr, $start, $stop, $str, $count) = split /\t/;
	$chr = substr($chr, 3);
	$tagcounts{$chr}{$stop}{$str} += $count if $str eq "+";
	$tagcounts{$chr}{$start}{$str} += $count if $str eq "-";
}
close IN;
	
open IN, "<lai_2017_3utr_ends.gff";
while(<IN>){
	($chr, $start, $str) = (split /\t/)[0,3,6];
	$chr = substr($chr,3);
	($parentid) = ($_ =~ /(FBgn\d+)/);
	if($str eq '+'){
		if (! exists $utrEndpoint{$parentid}{$chr} || $start > $utrEndpoint{$parentid}{$chr}){
			$utrEndpoint{$parentid}{$chr} = $start;
		}
	}
	else{
		if (! exists $utrEndpoint{$parentid}{$chr} || $start < $utrEndpoint{$parentid}{$chr}){
			$utrEndpoint{$parentid}{$chr} = $start;
		}
	}
}
close IN;

#choose representative transcript with longest UTR for each unique 3' utr start position
open IN, "zgrep -P '\tthree_prime_UTR\t' $file | ";
while(<IN>){
	($start, $stop, $str, $last) = (split /\t/)[3,4,6,-1];
	$mypos = ($str eq '+') ? $start : $stop;
	($id) = ($last =~ /Parent=(FBtr.*)/);
	@ids = split /,/, $id;
	foreach $id (@ids){
		next if $utrStartpoint{$id} != $mypos; #skip consideration if dealing with a 3' UTR exon that isn't attached to stop codon
		$parent = $id2parent{$id};
		$parent2str{$parent} = $str;
		if (! defined $reptranscript{$parent}{$mypos}){
			$reptranscript{$parent}{$mypos} = $id;
		}
		if ($lengths{$id}{"three_prime_UTR"} > $lengths{$reptranscript{$parent}{$mypos}}{"three_prime_UTR"}){
			$reptranscript{$parent}{$mypos} = $id;
		}
	}
}
close IN;

#find position of most distal stop codon, and corresponding ID for representative transcript harboring it
#compile a list of representative transcript IDs
for $parentid ( keys %reptranscript ) {
	for $mypos ( keys %{ $reptranscript{$parentid} } ) {
		$str = $parent2str{$parentid};
		if($str eq '+'){ 
			if (! exists $distalstop{$parentid} || $mypos > $distalstop{$parentid}){
				$distalstop{$parentid} = $mypos;
				$distalstop_transcript{$parentid} = $reptranscript{$parentid}{$mypos};
			}
		}
		else{
			if (! exists $distalstop{$parentid} || $mypos < $distalstop{$parentid}){
				$distalstop{$parentid} = $mypos;
				$distalstop_transcript{$parentid} = $reptranscript{$parentid}{$mypos};
			}
		}
		$okids{$reptranscript{$parentid}{$mypos}}=1;
	}
}

#store exon boundaries
open IN, "zgrep -P '\texon\t' $file | ";
while(<IN>){
	($chr, $start, $stop, $str) = (split /\t/)[0,3,4,6];
	if($str eq '+'){ $exonpos{$chr}{$start}{$str} = 1; }
	else{ $exonpos{$chr}{$stop}{$str} = 1; }
}
close IN;

#determine 3' UTR boundaries of largest possible region, and terminal 3' UTR positions for representative transcripts
open IN, "zgrep -P '\tthree_prime_UTR\t' $file | ";
while(<IN>){
	($chr, $region, $start, $stop, $str, $last) = (split /\t/)[0,2,3,4,6,-1];
	if($str eq '+'){
		if (! exists $allregions{$chr}{$start} || $stop > $allregions{$chr}{$start}){
			$allregions{$chr}{$start} = $stop;
		}
	}
	else{
		if (! exists $allregions{$chr}{$stop} || $start < $allregions{$chr}{$stop}){
			$allregions{$chr}{$stop} = $start;
		}
	}
	($id) = ($last =~ /Parent=(FBtr.*)/);
	@ids = split /,/, $id;
	foreach $id (@ids){
		if ($okids{$id}){ #if representative transcript for coding gene
			if($str eq '+'){ 
				if (! exists $lastpos{$id} || $stop > $lastpos{$id}){
					$lastpos{$id} = $stop;
				}
			}
			else{
				if (! exists $lastpos{$id} || $start < $lastpos{$id}){
					$lastpos{$id} = $start;
				}
			}
		}
	}
}
close IN;

#build 3'UTR models for representative transcripts
#keep FlyBase annotations if not 3' UTR region
#if most distal 3' UTR exon, link it to most distal FlyBase or Lai 3' UTR annotations
#if not most distal 3' UTR but overlaps last exon, extend it to most distal FlyBase or Lai 3' UTR annotations
#if not most distal 3' UTR but doesn't overlap last exon, search for 3' UTR within the intron using mapped 3'seq reads from Lai et al and our study
open IN, "zgrep -P '\tfive_prime_UTR|three_prime_UTR|CDS\t' $file | ";
while($line = <IN>){
	($chr, $region, $start, $stop, $str, $last) = (split /\t/, $line)[0,2,3,4,6,-1];
	($id) = ($last =~ /Parent=(FBtr.*)/);
	@ids = split /,/, $id;
	foreach $id (@ids){
		if ($okids{$id}){
			@a = split /\t/, $line;
			$a[3] -= 1; #convert to 0-based coordinate system
			$a[4] -= 1; #convert to 0-based coordinate system
			$a[-1] = $id.':'.$id2parent{$id};
			$mostdistal = $utrEndpoint{$id2parent{$id}}{$chr};
			if($a[2] ne "three_prime_UTR" || ($lastpos{$id} != $a[3]+1 && $lastpos{$id} != $a[4]+1)){ $_ = join ("\t", @a)."\n"; }
			elsif ($distalstop_transcript{$id2parent{$id}} ne $id){ #if is last exon of 3' UTR but not transcript with the most distal stop codon
				if ($a[6] eq '+'){ $a[4] = findUTR($id, $a[0], $allregions{$chr}{$start}, $mostdistal, $a[6]); }
				else{ $a[3] = findUTR($id, $a[0], $allregions{$chr}{$stop}, $mostdistal, $a[6]); }
				$_ = join ("\t", @a)."\n";
			}
			else{ #most distal exon of most distal 3' UTR
				if ($a[6] eq '+'){ $a[4] = $mostdistal; }
				else{ $a[3] = $mostdistal; }
				$_ = join ("\t", @a)."\n";
			}
			### initially considered case here if ORFs existed without any annotated 3' UTRs, but in practice including it hardly changes outcome
			print "chr".$_ if $chr !~ /Un|2110000222/;
		}
	}
}
close IN;

sub findUTR{
	local ($id, $chr, $searchpos, $lastpos, $str) = @_;
	local $mypos = $searchpos;
	if ($str eq '+'){ for($i = $searchpos+1; $i <= $lastpos; ++$i){ last if ($exonpos{$chr}{$i}{$str}); $mypos = $i if ($tagcounts{$chr}{$i}{$str} >= 10 || $i == $lastpos) && $i > $mypos; } } #last if hit start of any other exon $stopcodon+(2200+$intronlen3pUTR{$id})
	else{ 		for($i = $searchpos-1; $i >= $lastpos; --$i){ last if ($exonpos{$chr}{$i}{$str}); $mypos = $i if ($tagcounts{$chr}{$i}{$str} >= 10 || $i == $lastpos) && $i < $mypos; } } #$stopcodon-(2200+$intronlen3pUTR{$id})
	return $mypos;
}v
