gunzip plfold_35x30_sitetypes.stats.gz 3pseq_s2_*.gz

###For Fig3A/B the following 3 files were used (processed by Igor Ulitsky)
###3pseq_s2_rawReads.bed -- raw 3P-seq tags mapped to the genome in genome build dm6 coordinates
###3pseq_s2_30nt_window_clusters.bed -- local read counts within 30nt window clustered together and summed, centered at predominant poly(A) junction
###3pseq_s2_30nt_window_cluster_annotations.bed -- annotation of the clustered tags relative to known poly(A) sites with support from RNA-seq data to link tags to genes

Rscript Fig3AB.R

### code to identify target site and surround context from Fasta files of UTRs and ORFs
#export SEQBASE=./dm6_s2_90%.fa #these are subset of UTRs for which >=90% map to single dominant 3p-seq annotated isoform
#export SEQBASEORF=./dm6.flybase.orfs
#export SEQBASE5UTR=./dm6.flybase.5utrs
#export SITES=./targetsites
#export MIRSEQFILE=./microRNAs.txt
#./find_targetsites.pl

### optimize 3p_energy
#./get_3pscore_results.pl 9 >3pscore_8x13_sitetypes_tarlen=9.txt
#mkdir rna_probs_90
#cd rna_probs_90
#cat ../dm6_s2_90%.fa | RNAplfold -L 40 -W 80 -u 30
#cd ..

### optimize structural accessibility
#./get_plfold_results.pl >plfold_35x30_sitetypes.stats

### this script will not work without the associated genome-wide alignments of fly 3' UTRs, but can run 
### if branch length score computation skipped; would have to manually modify to do so
### I provide it's pre-computed output for later steps
#./collect_features.pl >feature_table_7_8mers_only.txt

Rscript Fig3C.R 3pscore_8x13_sitetypes_tarlen=9.txt

Rscript Fig3D.R plfold_35x30_sitetypes.stats

#for i in {1..1000}; do { 
#	bsub -e sub.err -o sub.out "Rscript prepare_datasets.R feature_table_7_8mers_only.txt $i; Rscript regression_models.R $i >data/fly8mer7mer-$i.out";
#} done;

### exact results may vary from paper due to noise from random sampling procedure
### check all data was generated correctly

#./Table1_fraction_features_chosen.pl >Table1_chosenFraction.txt

Rscript FigS3.R

Rscript Fig3E.R
