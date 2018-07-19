require(graphics)
options(warn=-1)

ball = c()
call = c()
dall = c()
eall = c()
fall = c()
gall = c()
hall = c()

pdf("Fig1.pdf")

for (i in args){
	miR = strsplit(strsplit(i,"\\/")[[1]][2],'_')[[1]][1]
	say("Processing...", miR, "...fraction of functional sites (used for data presented in Fig2G)")
	a = read.delim(i, header=T, sep="\t")

####For Fig1A
#	b = a[a$utr8 == 1 & a$utr8C == 0 & a$utr8G == 0 & a$utr8T == 0, "fc"] # & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0
#	c = a[a$utr8 == 0 & a$utr8C == 1 & a$utr8G == 0 & a$utr8T == 0, "fc"] # & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0
#	d = a[a$utr8 == 0 & a$utr8C == 0 & a$utr8G == 1 & a$utr8T == 0, "fc"] # & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0
#	e = a[a$utr8 == 0 & a$utr8C == 0 & a$utr8G == 0 & a$utr8T == 1, "fc"] # & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0
#	f = a[a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0, "fc"] # & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0 & a$utr5p5 == 0  & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0
#	g = a[a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0, "fc"] # & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0 & a$utr5p5 == 0  & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0
#	h = a[a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0, "fc"] # & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0 & a$utr5p5 == 0  & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0
#	say(length(b),length(c),length(d),length(e),length(f),length(g),length(h))

####For Fig1C (for 3' UTR)
	b = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 1 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
	c = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 1 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
	d = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 1 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
	e = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 1 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
	f = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 1 & a$utr6A1 == 0, "fc"]
	h = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 1, "fc"]
	g = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]

####For Fig1D (for ORF)
#	b = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 1 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
#	c = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 1 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
#	d = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 1 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
#	e = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 1 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
#	f = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 1 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
#	h = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 1 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
#	g = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]

####For Fig1E (for 5' UTR)
#	b = a[a$utr5p8 == 1 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
#	c = a[a$utr5p8 == 0 & a$utr5p7m8 == 1 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
#	d = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 1 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
#	e = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 1 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
#	f = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 1 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
#	h = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 1 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]
#	g = a[a$utr5p8 == 0 & a$utr5p7m8 == 0 & a$utr5p7A1 == 0 & a$utr5p6 == 0 & a$utr5po6 == 0 & a$utr5p6A1 == 0 & a$orf8 == 0 & a$orf7m8 == 0 & a$orf7A1 == 0 & a$orf6 == 0 & a$orfo6 == 0 & a$orf6A1 == 0 & a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0, "fc"]

	par(mar=c(7,7,5,5), mgp = c(5, 1, 0))
	plot.ecdf(b, xlim=c(-2,2), ylim=c(0,1), verticals= TRUE, do.points = FALSE, col="purple", main=paste("Effect of miRNA sites, ", miR),
		ylab="Cumulative fraction", xlab="Fold change (log2)", col.01line = "white", lwd=2, cex.axis=2, cex.lab=2, cex.main=1, bty="n")
	plot.ecdf(c, verticals= TRUE, do.points = FALSE, add = TRUE, col = "red", col.01line = "white", lwd=2, bty="n")
	plot.ecdf(d, verticals= TRUE, do.points = FALSE, add = TRUE, col = "blue", col.01line = "white", lwd=2, bty="n")
	plot.ecdf(e, verticals= TRUE, do.points = FALSE, add = TRUE, col = "cyan", col.01line = "white", lwd=2, bty="n")
	plot.ecdf(f, verticals= TRUE, do.points = FALSE, add = TRUE, col = "grey", col.01line = "white", lwd=2, bty="n")
	plot.ecdf(h, verticals= TRUE, do.points = FALSE, add = TRUE, col = "darkorange", col.01line = "white", lwd=2, bty="n")
	plot.ecdf(g, verticals= TRUE, do.points = FALSE, add = TRUE, col = "black", col.01line = "white", lwd=2, bty="n")

	legend("topleft", bg="white", bty="n", legend = 
	    c( #change alternative="less" if dealing w/ a miRNA knockout or knockdown dataset, for one-sided K-S tests
		paste("8mer P < ", formatC(ks.test(b,g,alternative="greater")$p.value, digits = 2, format = 'g'), " (", length(b), ")", sep=""),
		paste("7mer-m8 P < ", formatC(ks.test(c,g,alternative="greater")$p.value, digits = 2, format = 'g'), " (", length(c), ")", sep=""),
		paste("7mer-A1 P < ", formatC(ks.test(d,g,alternative="greater")$p.value, digits = 2, format = 'g'), " (", length(d), ")", sep=""),
		paste("6mer P < ", formatC(ks.test(e,g,alternative="greater")$p.value, digits = 2, format = 'g'), " (", length(e), ")", sep=""),
		paste("offset-6mer P < ", formatC(ks.test(f,g,alternative="greater")$p.value, digits = 2, format = 'g'), " (", length(f), ")", sep=""),
		paste("6mer-A1 P < ", formatC(ks.test(h,g,alternative="greater")$p.value, digits = 2, format = 'g'), " (", length(h), ")", sep=""),
		paste("No site", " (", length(g), ")", sep="")), text.col = c("purple","red","blue","cyan","grey","darkorange","black"))

	print(paste("8mer:", formatC(ks.test(b,g,alternative="greater")$statistic, digits = 2, format = 'g'), sep=''))
	print(paste("7mer-m8:", formatC(ks.test(c,g,alternative="greater")$statistic, digits = 2, format = 'g'), sep=""))
	print(paste("7mer-A1:", formatC(ks.test(d,g,alternative="greater")$statistic, digits = 2, format = 'g'), sep=""))
	print(paste("6mer:", formatC(ks.test(e,g,alternative="greater")$statistic, digits = 2, format = 'g'), sep=""))
	print(paste("Offset 6mer:", formatC(ks.test(f,g,alternative="greater")$statistic, digits = 2, format = 'g'), sep=""))
	print(paste("6mer-A1:", formatC(ks.test(h,g,alternative="greater")$statistic, digits = 2, format = 'g'), sep=""))

	ball = c(ball, b)
	call = c(call, c)
	dall = c(dall, d)
	eall = c(eall, e)
	fall = c(fall, f)
	gall = c(gall, g)
	hall = c(hall, h)
}

par(mar=c(7,7,5,5), mgp = c(5, 1, 0))
plot.ecdf(ball, xlim=c(-2,2), ylim=c(0,1), verticals= TRUE, do.points = FALSE, col="purple", main="Effect of miRNA sites, pooled",
	ylab="Cumulative fraction", xlab="Fold change (log2)", col.01line = "white", lwd=2, cex.axis=2, cex.lab=2, cex.main=1, bty="n")
plot.ecdf(call, verticals= TRUE, do.points = FALSE, add = TRUE, col = "red", col.01line = "white", lwd=2, bty="n")
plot.ecdf(dall, verticals= TRUE, do.points = FALSE, add = TRUE, col = "blue", col.01line = "white", lwd=2, bty="n")
plot.ecdf(eall, verticals= TRUE, do.points = FALSE, add = TRUE, col = "cyan", col.01line = "white", lwd=2, bty="n")
plot.ecdf(fall, verticals= TRUE, do.points = FALSE, add = TRUE, col = "grey", col.01line = "white", lwd=2, bty="n")
plot.ecdf(hall, verticals= TRUE, do.points = FALSE, add = TRUE, col = "darkorange", col.01line = "white", lwd=2, bty="n")
plot.ecdf(gall, verticals= TRUE, do.points = FALSE, add = TRUE, col = "black", col.01line = "white", lwd=2, bty="n")

legend("topleft", bg="white", bty="n", legend = 
    c(paste("8mer P < ", formatC(ks.test(ball,gall,alternative="greater")$p.value, digits = 2, format = 'g'), " (", length(ball), ")", sep=""),
	paste("7mer-m8 P < ", formatC(ks.test(call,gall,alternative="greater")$p.value, digits = 2, format = 'g'), " (", length(call), ")", sep=""),
	paste("7mer-A1 P < ", formatC(ks.test(dall,gall,alternative="greater")$p.value, digits = 2, format = 'g'), " (", length(dall), ")", sep=""),
	paste("6mer P < ", formatC(ks.test(eall,gall,alternative="greater")$p.value, digits = 2, format = 'g'), " (", length(eall), ")", sep=""),
	paste("offset-6mer P < ", formatC(ks.test(fall,gall,alternative="greater")$p.value, digits = 2, format = 'g'), " (", length(fall), ")", sep=""),
	paste("6mer-A1 P < ", formatC(ks.test(hall,gall,alternative="greater")$p.value, digits = 2, format = 'g'), " (", length(hall), ")", sep=""),
	paste("No site", " (", length(gall), ")", sep="")), text.col = c("purple","red","blue","cyan","grey","darkorange","black"))

#fraction of functional sites (used for data presented in Fig2G)
print("")
"ALL DATA POOLED:"
paste("8mer:", formatC(ks.test(ball,gall,alternative="greater")$statistic, digits = 2, format = 'g'), sep='')
paste("7mer-m8:", formatC(ks.test(call,gall,alternative="greater")$statistic, digits = 2, format = 'g'), sep="")
paste("7mer-A1:", formatC(ks.test(dall,gall,alternative="greater")$statistic, digits = 2, format = 'g'), sep="")
paste("6mer:", formatC(ks.test(eall,gall,alternative="greater")$statistic, digits = 2, format = 'g'), sep="")
paste("Offset 6mer:", formatC(ks.test(fall,gall,alternative="greater")$statistic, digits = 2, format = 'g'), sep="")
paste("6mer-A1:", formatC(ks.test(hall,gall,alternative="greater")$statistic, digits = 2, format = 'g'), sep="")
print("")
"ALL PAIR-WISE COMPARISONS BETWEEN SITE CLASSES (USED TO BUILD TABLE S2)"
sites=c("8mer", "7mer-m8", "7mer-A1", "6mer", "Offset 6mer", "6mer-A1", "No site")
say(c("",sites))
j=1
for (x in list(ball, call, dall, eall, fall, hall)){
	tmp=c()
	for (y in list(ball, call, dall, eall, fall, hall, gall)){
		tmp=c(tmp,formatC(p.adjust(ks.test(x,y,alternative="greater")$p.value, method="bonferroni", n=36), digits = 2, format = 'g'))
	}
	say(c(sites[j], tmp))
	j=j+1
}
