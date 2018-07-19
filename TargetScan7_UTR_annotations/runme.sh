#######prepare 3p-seq data from tissues/cell lines mapped and counted
#######download FlyBase6.19 annotations
#wget ftp://ftp.flybase.net/releases/FB2017_06/dmel_r6.19/gff/dmel-all-no-analysis-r6.19.gff.gz

#cat bed/* >3Pseq_pA_allstages.bed.gz
#gunzip 3Pseq_pA_allstages.bed.gz
#grep -P '^chr' 3Pseq_pA_allstages.bed | perl -ne '@a = split / /; print join("\t", @a); ' | gzip -c >3Pseq_pA_allstages.bed.gz
#Rscript compute_bedfile_sum.R
#gzip 3Pseq_pA_allstages.summed
#rm 3Pseq_pA_allstages.bed*

#######script to annotate longest 3' UTRs for each unique stop codon by utilizing 3'-seq-revised isoform information (derived from Sanfilippo et al., 2017, Additional File 7)
./choose_all_genes_for_TargetScan.pl | gzip -c >FlyBase-r6.19.forTScan.gff.gz
