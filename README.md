<img src="logo.png" width="400">

# TargetScan supplementary tools

This repository is intended to accompany our submission. For more information please refer to:

Agarwal V, Subtelny AO, Thiru P, Ulitsky I, Bartel DP. [Predicting microRNA targeting efficacy in Drosophila](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1504-3).
 _**Genome Biology**_, 19:152. (2018).

The code is released to enhance reproducibility and as a suite of complementary tools to [TargetScan](http://www.targetscan.org/fly_72/) 
in the hope it might help others in the future who work on new datasets.

These tools can be used in a variety of organisms to:

* identify putative miRNA-binding sites in annotated genes and visualize miRNA-mediated repression in miRNA perturbation 
datasets (i.e., transfection, induction, knockdown, or knockout datasets) (as shown in Fig 1)
* perform statistical analyses on the properties of short motif evolution using genome-wide multiple sequence alignments of UTRs
(as shown in Fig 2)
* optimize and compute context features to help train machine learning algorithms (as shown in Fig 3)
* compare the performance of diverse prediction algorithms on test sets (as shown in Fig 4)

For a codebase to compute context scores for flies or other insects while incorporating 
3' UTR isoform information, the code provided in [TargetScan](http://www.targetscan.org/cgi-bin/targetscan/data_download.fly72.cgi) 
is recommended for use instead of this code.

If you find our code or our precomputed fly miRNA target predictions to be helpful for your work, please cite the paper above.

To better understand the methodological details for the evolutionary analyses, the following resources
 describe the original implementation \[1\], extension of parameters to worm and fly \[2\], and the re-implemented pipeline \[3\].

1. Friedman RC, Farh KK, Burge CB, Bartel DP. [Most Mammalian mRNAs Are Conserved Targets of MicroRNAs](http://genome.cshlp.org/content/19/1/92.full.pdf). _**Genome Research**_, 19:92-105 (2009).
2. Jan CH, Friedman RC, Ruby JG, Bartel DP. [Formation, regulation and evolution of Caenorhabditis elegans 3'UTRs](http://bartellab.wi.mit.edu/publication_reprints/Jan_Nature_2011.pdf). _**Nature**_, 469:97-102. (2011).
3. Agarwal V, Bell GW, Nam J, Bartel DP. [Predicting effective microRNA target sites in mammalian mRNAs](http://elifesciences.org/content/4/e05005/). _**eLife**_, 4:e05005, (2015).


# Dependencies for running entire pipeline (mostly optional)
* [RNAhybrid](https://bibiserv2.cebitec.uni-bielefeld.de/rnahybrid?id=rnahybrid_view_download)

* BranchLengthScoring.py script from [MotifMap](http://motifmap.ics.uci.edu/) and 
[newick 1.2](http://www.daimi.au.dk/~mailund/newick/newick-1.2.tar.gz)

* Installation of the [PHAST](http://compgen.cshl.edu/phast/) package

* Local installation of a mysql server. Replace this line in scripts of Fig2 with your server information:  
  $dbh = DBI->connect("dbi:mysql:database=username;host=HOSTNAME", "username", "password");

* The following perl modules:
  * DBI module
  * [BioPerl](http://bioperl.org/), in particular Bio::AlignIO and Bio::SimpleAlign
  * IO::Uncompress::Gunzip

* [UCSC tools](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads) installation, including mafFrag 
and corresponding genome-wide [multiz27way](http://hgdownload.cse.ucsc.edu/goldenPath/dm6/multiz27way/) alignments

* An LSF-based computing cluster that can submit jobs with the "bsub" command

* Local installation of Matlab

* Numerous R libraries listed individually in each R script, including R Bioconductor libraries


# Instructions for use

Not all code may work immediately because some pieces depend on computing environment, and not all intermediate 
files are provided because some are too large. For R code to work properly, please copy the contents of 
.Rprofile in this folder to your local .Rprofile. Exporting the allfxns.pm module to PERL5LIB might also be required.

Users are advised to read the code closely and modify commented pieces as appropriate to acquire 
desired output for your environment. For example, you will need to download all of the additional 
R library and Perl module dependencies for the code to work. This being said, if you find crucial 
files are missing, making the code unusable, or if you identify a major problem in the code, please 
raise a Github issue.

In each Figure's folder, change directories to it and then run the script "bash runme.sh".
Please read this file first as it provides a general overview of relevant commands that were used sequentially 
to pre-process the data and generate the figures.
This script should be able to run on the precomputed data provided in the folder to generate the figures.

# Additional notes

Our naming convention is slightly different in the code than in the paper. In particular, the "HYBRIDSCORE" 
and "PLFOLD" features in the code are equivalent to "3p_energy" and "SA" features in the paper, respectively.
