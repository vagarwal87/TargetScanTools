library(fields)
source("pcor.R")

pdf("Fig3D.pdf", width=10, height=4.8)
par(mar=c(7, 7, 5, 5), mgp = c(5, 1, 0), cex.axis=2, cex.lab=2, bty="n")

file = args[1]
names = strsplit(basename(args),"\\.")[[1]]

b = read.delim(file)

lowerQ=quantile(b$FC, 0.25)
diff=quantile(b$FC, 0.75)-quantile(b$FC, 0.25)
for(x in unique(b$MIR_ID)){
	a=b$FC[b$MIR_ID == x]
	b$FC[b$MIR_ID == x] = scale(a, center=quantile(a, 0.25),scale=quantile(a, 0.75)-quantile(a, 0.25))*diff + lowerQ
}

b.matrix = model.matrix(FC~TYPE, b)[,2:3]
corr.mat = sapply(4:ncol(b), function(x){ pcor.test(b$FC, -log10(b[,x]), b.matrix, use = c("mat"), method = c("pearson"))[[1]] }) #log10()
dim(corr.mat)=c(30,78)
corr.mat=t(corr.mat)

corr.mat[corr.mat < 0 ] = 0
corr.mat=corr.mat[,1:25]
image.plot(-35:42, 1:ncol(corr.mat), corr.mat, xlab="Position relative to site", ylab="Window size", legend.width=2, legend.mar = 15, zlim=c(0, 0.156)) #5
