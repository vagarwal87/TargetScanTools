X1=read.delim("Fig2G_fractions.txt", row.names=1)
X1$avg=apply(X1[,3:8],1,mean)
X1$sd=apply(X1[,3:8],1,sd)

X1
"Fraction repressed is computed in ../Fig1/Fig1_cdfs.R for 5' UTR and 3' UTR sites"

col=c("purple", "red", "blue", "cyan", "darkgrey","darkorange","purple", "red", "blue")

pdf("Fig2G.pdf")
par(mar=c(9,9,5,5), mgp = c(5, 1.5, 0))
plot(X1$FracConserved, X1$avg, cex.axis=2, cex.lab=2, bty="n", ylab="Fraction of sites conferring\nmRNA destabilization", 
xlab="Fraction of sites conserved\nabove background (cutoff 1.0)", pch=19, xlim=c(0, 0.6), ylim=c(0, 0.9), 
col = col, las=1, cex=1.5)
error.bar(X1$FracConserved, X1$avg, X1$sd, col=col)
plot(X1$FracRepressed, X1$FracConserved, cex.axis=2, cex.lab=2, bty="n", xlab="Fraction of sites conferring\nmRNA destabilization", 
ylab="Fraction of sites conserved\nabove background (cutoff 1.0)", pch=19, xlim=c(0, 0.6), ylim=c(0, 0.6), 
col = c("purple", "red", "blue", "cyan", "darkgrey","darkorange","purple", "red", "blue"), las=1, cex=1.5)
dev.off()
