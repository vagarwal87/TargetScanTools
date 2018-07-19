pdf("FigS1E.pdf", width=7, height=7)
par(mar=c(5, 9, 8, 8), mgp = c(5, 1, 0))

N = 5 #number of top miRNAS to consider

a = read.table("fly_miRNA_counts.txt", header=F, row.names=1, sep="\t")
b = a[1:N,]
b$V2=as.character(b$V2)
b = rbind(b, c("Other", sum(a[(N+1):nrow(a),2])))
b$V3=as.numeric(b$V3)
b
pie(b[,2], labels = b[,1],col=c(rainbow(N),"black"), clockwise = T, cex.main = 2, cex=1.5)
