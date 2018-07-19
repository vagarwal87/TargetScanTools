#!/usr/bin/perl

use allfxns;

$tmp = $ENV{TMP}; #store job output files in this TMP directory environment variable

=pod
		ACCEPTED OPTIONS:			DEFAULTS:
cmd		program to run
exten		extension of output files	out
indir		location of files
instr		input command			<
onefile	output goes to one file		0
outdir	destination of output files	HOME
outstr	output command			>
qsub		use cluster				0 (false)
skipuniq	files to exclude			""
uniq		files to include

		QSUB/BSUB OPTIONS:		DEFAULTS:
cmdperjob	commands per job			1
jobfolder	where to keep tmp job files	HOME
suberr	job output error file		HOME/sub.err"
subopts	qsub options			"-V" (qsub) or "-q bartel" (bsub) <<<<>>>> -l h_vmem=2G
subout	job output file			HOME/sub.out
=cut

$species = "fly6_07";
`mkdir $species` if !(-d "$species");
$normal = 0;

sub options(){
   local ($klen, $region, $sitetype) = @_;
   %opts = (
#	"indir"    =>	"./multifasta_$species/$region",
#	"outdir"   =>	"./$species/$region/bls1to2mer.REV",
#	"uniq"     =>	"mfa",
#	"cmd"	     =>	"./bin_MSA.pl $region 1 all", #for median of entire UTR
##	"cmd"	     =>	"./bin_MSA.pl $region 2", #for binned 1-2mer counts can also use "all" here
##	"cmd"	     =>	"./bin_MSA.pl $region 8", #for conditionals use 2 or 8
#	"exten"    =>	"bls",
#	"cmdperjob"=>	20,

#	"indir"    =>	"./$species/$region/bls1to2mer.REV", #or 6to8mer
#	"outdir"   =>	"./$species/$region/counts1to2mer.REV", #or swap to "6to8mer"
#	"uniq"     =>	"bls",
#	"cmd"	     =>	"./get_markov_counts.pl $region",
#	"exten"    =>	"counts",
#	"cmdperjob"=>	50,

	"indir"    =>	"./$species/$region/logprobs.REV",
	"outdir"   =>	"./$species/$region/s2b.REV",
	"uniq"     =>	"$klen.probs",
	"cmd"	     =>	"./eval_50kmers_withBins.pl $region $sitetype",
	"exten"    =>	"$sitetype.s2b",
	"cmdperjob"=>	1,

	"skipsame" =>	0, #overwrite files if 0 or skip existing files if 1
	"onefile"  =>	0, #only use 1 with "bin_MSA.pl $species $region 1 all", use 0 otherwise
	"bsub"     =>	1, #use bsub or qsub
	"subopts"  =>	"-oo $tmp/jobsub.out" #bsub options
	);
}

if($normal){
	&options('', '', '');
	parallelize(\%opts);
	die("Normal mode");
}

foreach $region (qw/five_prime_UTR three_prime_UTR/){
###	comment out if using for loop
#	&options('', $region, '', '');
#	print $opts{cmd}, "\n";
#	`mkdir $species/$region` if !(-d "$species/$region");
#	parallelize(\%opts);

###	comment out if using main loop above
	foreach $sitetype (qw/8m 7m8 7A1 6m o6 6A1/){
		$klen = 8;
		$klen = 7 if $sitetype =~ /7A1|7m8/;
		$klen = 6 if $sitetype =~ /6m|o6|6A1/;
		&options($klen, $region, $sitetype);
		print $opts{cmd}, "\n";
		parallelize(\%opts);
	}
}
