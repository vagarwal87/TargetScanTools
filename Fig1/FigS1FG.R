require(graphics)
targets = list()
nontargets = list()

for (i in args){
	miR = strsplit(strsplit(basename(i),"\\.")[[1]][1],'_')[[1]][2:2]

	a = read.table(i, header=T, sep="\t")
	b = a[a$utr8 == 1 | a$utr7m8 == 1 | a$utr7A1 == 1, 2]
	c = a[a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0, 2]

	b = b-median(c)
	c = c-median(c)

	targets[[miR]] = c(targets[[miR]], b)
	nontargets = c(nontargets, c)
}

pdf("Rplots.pdf", width=8.5, height=7)
par(mar=c(7, 7, 7, 13), mgp = c(5, 1, 0))
plot.ecdf(as.numeric(nontargets), xlim=c(-1,1), verticals= TRUE, do.points = FALSE, col="black", ylab="Cumulative fraction", xlab="Fold change (log2)", lwd=2, cex.axis=2, xaxt="n", cex.lab=2, bty="n", col.01line = "white", main="", las=1)

pvals = c()
colors = rainbow(5)
ordered = read.delim("endogenous_microRNAs.txt",header=F)
ordered = as.factor(ordered$V1)
j = 1
for (i in ordered) {
	say(i)
	plot.ecdf(as.numeric(targets[[i]]), verticals= TRUE, do.points = FALSE, add = TRUE, col = colors[j], lwd=2, bty="n", col.01line = "white", las=1)
	pvals[j] = formatC(wilcox.test(as.numeric(targets[[i]]),as.numeric(nontargets),alternative="greater")$p.value, digits = 2, format = 'g')
	j = j+1
}

plot.ecdf(as.numeric(nontargets), verticals= TRUE, do.points = FALSE, col="black", ylab="Cumulative fraction", xlab="Fold change (log2)", lwd=2, cex.axis=2, cex.lab=2, bty="n", col.01line = "white", las=1, add=T)
axis (1, at=c(-1, 0, 1), labels=c(-1, 0, 1), cex.axis=2)
par(xpd=TRUE)
legend("bottomright", inset = c(-0.5, 0), cex=1.7, bg="white", bty="n", legend = c(paste(ordered, pvals), "No site"), text.col = c(colors, "black"))
