library(Biostrings)

UTR5=read.delim("fly6_07/five_prime_UTR/s2b.REV/1.00.8m.s2b", header=F)
UTR5=UTR5[1:92,1:3]
UTR5$sig = as.numeric(unlist(strsplit(as.character(UTR5$V3),'-'))[seq(1,length(UTR5$V3)*2,2)])
UTR5$back = as.numeric(unlist(strsplit(as.character(UTR5$V3),'-'))[seq(2,length(UTR5$V3)*2,2)])
UTR5$s2b = log2(UTR5$sig/UTR5$back)

UTR3=read.delim("fly6_07/three_prime_UTR/s2b.REV/1.00.8m.s2b", header=F)
UTR3=UTR3[1:92,1:3]
UTR3$sig = as.numeric(unlist(strsplit(as.character(UTR3$V3),'-'))[seq(1,length(UTR3$V3)*2,2)])
UTR3$back = as.numeric(unlist(strsplit(as.character(UTR3$V3),'-'))[seq(2,length(UTR3$V3)*2,2)])
UTR3$s2b = log2(UTR3$sig/UTR3$back)

merged=merge(UTR3, UTR5, by=1)
merged[,1]=substring(sapply(merged[,1], function(x) toString(reverseComplement(DNAString(as.character(x))))),2,8)
merged$bilaterian = as.factor(ifelse(merged[,1] %in% c("GAGGTAG","GGAATGT","GTCATGG","AAGGCAC","CCCTGAG","TGGTCCC","ATTGCTT","ATGGCAC","GGACGGA","GATATGT","ACTGGCC","AATACTG","TGTGCGT","AATCTCA","GATTGTC","AGCTGCC","AGCACCA","GGCAAGA","TGCATTG","GGCAGTG","TTGTTCG","GGAAGAC","CTTTGGT","ATTGCAC","TTGGCAC","ACCCGTA","CGGTGGG", "TCGTTGT", "TAAGTAC"), "Bilaterian", "Invertebrate"))
merged$bilaterian=relevel(merged$bilaterian, ref="Invertebrate")
table(merged$bilaterian)

"Spearman correlation"
cor(merged$s2b.x, merged$s2b.y, method='spearman')

pdf("Fig2F.pdf")
par(mar=c(9,9,5,5), mgp = c(5, 1.5, 0))
plot(merged$s2b.x, merged$s2b.y, cex.axis=2, cex.lab=2, bty="n", xlab="8mer signal-to-background ratio, 3' UTR (log2)", 
ylab="8mer signal-to-background ratio, 5' UTR (log2)", pch=19, xlim=c(-0.5, 2.5), ylim=c(-2, 2), cex=1.5, col=ifelse(merged$bilaterian=='Bilaterian','orange','blue'))

boxplot(merged$s2b.x ~ merged$bilaterian, ylim=c(-0.5,2.5), col=c("blue","orange"), main="S:B, 3' UTR", frame = F, cex.axis=1.2, cex.lab=1.5, cex.main=2, las=2, notch=T)
wilcox.test(merged$s2b.x ~ merged$bilaterian, alternative='less')

boxplot(merged$s2b.y ~ merged$bilaterian, ylim=c(-2,2), col=c("blue","orange"), main="S:B, 5' UTR", frame = F, cex.axis=1.2, cex.lab=1.5, cex.main=2, las=2, notch=T)
wilcox.test(merged$s2b.y ~ merged$bilaterian, alternative='less')
