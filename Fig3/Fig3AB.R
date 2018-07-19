pdf("Fig3A.pdf", width=7, height=7)
par(mar=c(14, 14, 10, 10), mgp = c(5, 1, 0))

a = read.table("Fly_3Pseq_sites.txt", header=F, sep="\t")
head(a)
lbls = paste(a[,1], ", ", a[,2], sep='')
say(lbls)
pie(a[,2], labels = lbls,col=c(rainbow(5),"black"), clockwise = T, cex.main = 2, cex=.7) #, init.angle=270 # in Liver Cells

pdf("Fig3B.pdf", width=7, height=9)
par(mar=c(20, 7, 5, 5), mgp = c(5, 1, 0))

b = read.table("novel_annotations.txt", header=F, sep="\t")
lbls = b[,1]

barplot(b[,2], log = "y", yaxt="n", ylim = c(1, 10000), cex.axis=2, cex.lab=2, bty="n", ylab="Number of poly(A) sites", col = "red4", las=2, cex=1.5, names.arg = lbls)
y1 <- floor(log10(range(b[,2])))
pow <- seq(0, 4)
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(2, 10^pow, las=1, cex.axis=1.5)
axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
