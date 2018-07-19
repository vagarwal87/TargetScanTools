a = read.delim("fly-results.out")
a = unique(a)

a = na.omit(a)
a=a^2

a = as.data.frame(a)
a$"context"=a[,2]
a[,2]=NULL
colnames(a)[2]="All subsets"
head(a)

median(a[,1])
median(a[,4])
median(a[,ncol(a)])

wilcox.test(a[,1], a[,4], paired = T)

pdf("FigS3.pdf")
par(mar=c(10, 6, 4, 4))
boxplot(a, col="red", main="Performance of Learning Procedures",
	xlab="", ylab=expression(paste("r"^"2", " to held-out data")), frame = F, cex.axis=1.2, cex.lab=1.5, cex.main=2, las=2, notch=T)

files = list.files("targetsites_normalized", "out", full.names=FALSE)
alltargets = c()
allnontargets = c()
allvars = sapply(files, function(file){
	tmp=read.delim(paste("targetsites_normalized", file, sep='/'))
	tmp$fc = (tmp$fc - mean(tmp$fc)) / sd(tmp$fc)
	targetids = tmp[tmp$utr8 > 0 | tmp$utr7m8 > 0 | tmp$utr7A1 > 0, "utr_id"]
	fcs = tmp[tmp$utr_id %in% targetids, "fc"]
	alltargets <<- c(alltargets, fcs)
	vars = sapply(files[files!=file], function(file2){
		tmp2=read.delim(paste("targetsites_normalized", file2, sep='/'))
		tmp2$fc = (tmp2$fc - mean(tmp2$fc)) / sd(tmp2$fc)
		fcs2 = tmp2[tmp2$utr_id %in% targetids & tmp2$utr8 == 0 & tmp2$utr7m8 == 0 & tmp2$utr7A1 == 0, "fc"]
		allnontargets <<- c(allnontargets, fcs2)
		var(fcs2)
	})
	thisvar = var(fcs)
	names(thisvar) = file
	thisvar=data.frame(vars = c(thisvar, vars))
	thisvar[match(files, rownames(thisvar)),]
})

colnames(allvars)=gsub("_sites.out", "", files)
barplot( allvars, las=2, beside=TRUE, col=c("red","blue","darkgrey","orange","purple","cyan"), border=F, ylab = "Variance" )
legend("topright", bg="white", bty="n", legend = colnames(allvars), text.col=c("red","blue","darkgrey","orange","purple","cyan"))
barplot( c(var(alltargets), var(allnontargets)), names.arg = c("predicted targets","predicted non-targets"), las=2, beside=TRUE, col="black", border=F, ylim=c(0,2.6), ylab = "Variance" )
dev.off()

"targets"
var(alltargets)
"non-targets"
var(allnontargets)
