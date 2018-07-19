tar xvvzf unpackage_me.tar.gz

###### process microarray data for each of 3 test sets, modify script as appropriate
# Rscript get_vals_from_microarray.R

#####requires Vienna RNA package to run
#mkdir rna_probs rna_probs_shortest
#cd rna_probs
#cat ../fly6_07.3utrs.fa | RNAplfold -L 40 -W 80 -u 30
#cd ../rna_probs_shortest
#cat ../fly6_07.3utrs.shortest.fa | RNAplfold -L 40 -W 80 -u 30
#cd ..

######## find predicted miRNA sites in mRNAs
#export SEQBASE=./fly6_07.3utrs.fa
#export SEQBASEORF=./dm6.flybase.orfs
#export SEQBASE5UTR=./dm6.flybase.5utrs
#export SITES=./targetsites
#export MIRSEQFILE=./microRNAs.txt
#export PLFOLD=./rna_probs
#./find_targetsites.pl
#./collect_features.pl >feature_table_knockouts.txt

#export MIRSEQFILE=./microRNAs_transfections.txt
#./find_targetsites.pl
#./collect_features.pl >feature_table_transfections.txt

#export SEQBASE=./fly6_07.3utrs.shortest.fa
#export SEQBASEORF=./dm6.flybase.orfs
#export SEQBASE5UTR=./dm6.flybase.5utrs
#export SITES=./targetsites_shortest
#export MIRSEQFILE=./microRNAs.txt
#export PLFOLD=./rna_probs_shortest
#./find_targetsites.pl
#./collect_features.pl >feature_table_knockouts_shortestUTR.txt

#export MIRSEQFILE=./microRNAs_transfections.txt
#./find_targetsites.pl
#./collect_features.pl >feature_table_transfections_shortestUTR.txt

### find predictions from other algorithms, extract contents of gzipped tarball:
# tar xvvzf pred.tar.gz
### modify commented lines in this script to accomodate knockouts vs transfection predictions summary
# Rscript combine_all_predictions.R pred/*

### modify code as appropriate (instructions in comments) to produce desired output
#Rscript compute_total_context_scores.R {feature_table_knockouts.txt,feature_table_knockouts_shortestUTR.txt,feature_table_transfections.txt,or feature_table_transfections_shortestUTR.txt}

### Figure 4A
Rscript Fig4.R transfections all mean
### Figure 4-S2A
Rscript Fig4.R transfections all median
### Figure 4B
Rscript Fig4.R transfections NoSite mean
### Figure 4-S1A
Rscript Fig4.R transfections 6merOnly mean

### Figure 4C
Rscript Fig4.R knockouts all mean
### Figure 4-S2B
Rscript Fig4.R knockouts all median
### Figure 4D
Rscript Fig4.R knockouts NoSite mean
### Figure 4-S1B
Rscript Fig4.R knockouts 6merOnly mean
