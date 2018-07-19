tar xvvzf unpackage_me.tar.gz

### cuffdiff/ contains the full list of fold changes for all genes in D. melanogaster
### diff_uncorrected/ contains fold changes computed from cuffdiff, with lowly expressed genes filtered out
### diff_plsr_corrected/ contains fold changes corrected by the first round of normalization (Figure 1-FS1B-D)

### code to identify target site counts from Fasta files of UTRs and ORFs
#export SEQBASE=./fly6_07.3utrs.fa
#export SEQBASEORF=./dm6.flybase.orfs
#export SEQBASE5UTR=./dm6.flybase.5utrs
#export SITES=./targetsites_plsr_normalized
#export MIRSEQFILE=./microRNAs.txt
#./site_counts.pl

#echo "normalizing data"
#Rscript normalization.R targetsites_plsr_normalized/*

echo "generating Fig1C data, uncomment different sections of script to evaluate Fig1A,D-E"
Rscript Fig1_cdfs.R targetsites_normalized/miR1_sites.out targetsites_normalized/miR124_sites.out targetsites_normalized/miR263a_sites.out targetsites_normalized/miR4_sites.out targetsites_normalized/miR92a_sites.out targetsites_normalized/miR994_sites.out

### count endegeous miRNA families from Ago1 IP data
#wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM280nnn/GSM280088/suppl/GSM280088_AGO1_IP_S2_raw_data.txt.gz ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
#zgrep -A1 dme mature.fa.gz >mirbase-v21-fly.fa ##changed Us to Ts afterwards
#./miR-counts.pl >fly_miRNA_counts.txt
#rm mature.fa.gz GSM280088_AGO1_IP_S2_raw_data.txt.gz
echo "generating FigS1E plots"
Rscript FigS1E.R

### code to identify target site counts for endogenous miRNAs from Fasta files of UTRs
#export SEQBASE=./fly6_07.3utrs.fa
#export SEQBASEORF=./dm6.flybase.orfs
#export SEQBASE5UTR=./dm6.flybase.5utrs
#export MIRSEQFILE=./endogenous_microRNAs.txt
#export SITES=./targetsites_uncorrected
#for x in diff_uncorrected/miR*.diff; do { ./find_UTR_sites.pl $x; } done
#for x in targetsites_normalized/*; do { X=`basename $x _sites.out`; tail -n +2 $x | cut -f 1-2 >diff_normalized/$X.diff; } done
#export SITES=./targetsites_normalized
#for x in diff_normalized/miR*.diff; do { ./find_UTR_sites.pl $x; } done

echo "generating FigS1FG plots"
Rscript FigS1FG.R targetsites_uncorrected/miR*_*_sites.out; mv Rplots.pdf FigS1F.pdf
Rscript FigS1FG.R targetsites_normalized/miR*_*_sites.out; mv Rplots.pdf FigS1G.pdf
