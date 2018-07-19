library(fields)
library(lattice)

source("pcor.R")

file = args[1]

pdf("Fig3C.pdf", width=10, height=10)
par(mar=c(7, 7, 5, 5), mgp = c(5, 1, 0), cex.axis=2, cex.lab=2, bty="n")

names = strsplit(basename(args),"\\.")[[1]]

b = read.delim(file)

lowerQ=quantile(b$FC, 0.25)
diff=quantile(b$FC, 0.75)-quantile(b$FC, 0.25)
for(x in unique(b$MIR_ID)){
	a=b$FC[b$MIR_ID == x]
	b$FC[b$MIR_ID == x] = scale(a, center=quantile(a, 0.25),scale=quantile(a, 0.75)-quantile(a, 0.25))*diff + lowerQ
}

b.matrix = model.matrix(FC~TYPE, b)[,2:3]
corr.mat = data.frame( cor = sapply(4:ncol(b), function(x){ pcor.test(b$FC, b[,x], b.matrix, use = c("mat"), method = c("pearson"))[[1]] }) ) #log10()
corr.mat$names = colnames(b)[4:ncol(b)]
corr.mat[corr.mat < 0 ] = 0

corr.mat$startPos = unlist(strsplit(corr.mat$names,"\\."))[seq(1,nrow(corr.mat)*2,2)]
corr.mat$windowSize = as.numeric(unlist(strsplit(corr.mat$names,"\\."))[seq(0,nrow(corr.mat)*2,2)])
corr.mat$startPos = as.numeric(gsub("X","", corr.mat$startPos))
corr.mat=corr.mat[corr.mat$windowSize >= 4, ]

head(corr.mat)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
levelplot(cor ~ startPos + windowSize, data=corr.mat, col.regions=jet.colors, xlim=c(18.5,8.5))
