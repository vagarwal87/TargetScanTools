gunzip utr3_seqs.fa.gz utr5_seqs.fa.gz 

### Get branch length scoring script from MotifMap
#wget http://motifmap.ics.uci.edu/downloads/BranchLengthScoring.py
#wget http://www.daimi.au.dk/~mailund/newick/newick-1.2.tar.gz
#tar xvvzf newick-1.2.tar.gz
#rm newick-1.2.tar.gz
### Test BranchLengthScoring.py is working...might have to sys.path.append() newick module to this script

##### The following scripts were used, in chronological order, to grab multiple sequence alignments from UCSC's multiz27way given a set of UTR annotations
### Requires a local mysql library to run. Must modify scripts according to username/password accordingly.

### For retriving UTR alignments:

### index_refseq_models.pl 				// get all uniq exons from GTF file and ID them
### splicing_graph_sql.pl				// create exon chains for spliced UTRs
### get_gene_slices_multiz.pl				// create multifasta file of gene chains, splicing together exons
### uniqid2geneid.pl					// get gene ids corresponding to uniqids in splicing graph file (optional)

### For computing UTR bins, phylogenetic trees, control kmers, count statistics (** use ezparallel/allfxns to parallelize most commands **):

### bin_MSA.pl						// find median bls for every base in alignment, can use to also find bls for all kmers in alignments (CDS or UTR ok)
### compute_bins.R					// split median bls's for each region into N bins
### make_trees_from_bins.pl 				// make N bins from median bls in region & scale UCSC tree to functional region (3'UTR, CDS, etc) (no ezparallel needed))
### get_markov_counts.pl				// count mono-/di-/tri- nucleotides, & 6-8mers for 3'UTR vs CDS (& for each frame), at each branch length
### ./aggregate_counts.pl three_prime_UTR		// aggregate counts from all files & at each branch length
### ./gen_control_kmers.pl three_prime_UTR	// get logprobs for all 6/7/8-mers (& in each reading frame if CDS)
### eval_50kmers_withBins.pl				// choose 50 kmers with similar expected freqs and calculate sig2back vals

tar xvvzf fly6_07.tar.gz

#Compile results from all branch lengths to generate Figs 2B-E
Rscript Fig2BCDE.R five_prime_UTR
Rscript Fig2BCDE.R three_prime_UTR

### Method to estimate the number of genes targeted at branch length cutoff of 1.0 (maximal sensitivity)
### Takes a while to run and requires database of UTR alignments, so precomputed results are provided
#for x in /lab/bartel3_ata/agarwal/databases/multifasta_fly6_07/five_prime_UTR/*.mfa; do { ./get_seqs.pl $x 2>&- >>utr5_seqs.fa; } done
#for x in /lab/bartel3_ata/agarwal/databases/multifasta_fly6_07/three_prime_UTR/*.mfa; do { ./get_seqs.pl $x 2>&- >>utr3_seqs.fa; } done
#grep -e '>' utr5_seqs.fa | cut -b 2- >utr5_ids.txt
#grep -e '>' utr3_seqs.fa | cut -b 2- >utr3_ids.txt
#./find_targetsites_and_bls.pl five_prime_UTR
#./find_targetsites_and_bls.pl three_prime_UTR
#./dump_sql_id2geneid.pl | uniq | grep -v CDS > id2geneid.txt
Rscript estimate_transcripts_targeted.R

Rscript Fig2F.R
Rscript Fig2G.R

### script that calculates final parameters of equation to fit Pct score from smoothed Branch length scores for each seed family
### requires Matlab to run and takes a little while to finish, the precomputed output is included already
#Rscript calc_Pct_params.R three_prime_UTR

Rscript Fig2H.R feature_table_7_8mers_only.txt
