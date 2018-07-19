library(Biobase)
library(GEOquery)
library(gcrma)
library(limma)
library(drosophila2.db)

# load series and platform data from GEO
#GEO <- "GSE25008" #miR34
GEO <- "GSE20202" #miR14
#GEO <- "EMEXP3785" #miR-277

if (!file.exists(GEO)){
	rawfiles <- getGEOSuppFiles(GEO)
	tarfile <- grep("\\.tar$", rownames(rawfiles), value = TRUE)
	say(tarfile)
	untar(tarfile, exdir = GEO)
}

ABfile <- paste(GEO, sep = "/", "AB.Rdata") 
if (!file.exists(ABfile)){
	celfiles <- list.files(GEO, pattern = "\\.CEL.gz", full.names = TRUE)
	print(celfiles)
	AB <- ReadAffy(filenames = celfiles)
	save(AB, file = ABfile)
} else load(ABfile)

EXPfile <- paste(GEO, sep = "/", "EXP.Rdata") 
if (!file.exists(EXPfile)){
	gset <- gcrma(AB, fast = FALSE)
	save(gset, file = EXPfile)
} else load(EXPfile)

#sml <- c("X","G1","X","G0","X","G1","X","G0","X","G1","X","G0","X","G1","X","G0","X","G1","X","G0"); # Day20, miR34 vs WT
sml <- c("G0","G0","G1","G1","X","X"); # miR14 vs WT
#sml <- c("G0","G0","G0","G1","G1","G1");  #miR277 vs WT
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

design <- cbind( WT = as.double(sml =="G0" ), KO = as.double(sml=="G1"))
fit <- lmFit(gset,design,method="robust")
cont.matrix <- makeContrasts(rKO=KO-WT, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, adjust="BH", sort.by="logFC", number=1e10)

#convert probe IDs to Refseq IDs
x <- drosophila2FLYBASE
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
tT$ID=rownames(tT)
head(xx)
head(tT)
tT$ID <- paste(xx[tT$ID])

#remove absent IDs & merge multiple ID entries w/ median value
tT <- tT[tT$ID != "NULL",] #remove IDs without ID
tT <- aggregate(tT[,-ncol(tT)], by=list(tT$ID), median)
colnames(tT)[1]="ID_REF"

head(tT)
plot(tT$AveExpr, tT$logFC, cex=.1)
abline(v=2.5)
abline(v=median(tT$AveExpr),col='blue')

#write.table(tT[tT$AveExpr >= median(tT$AveExpr),c(1,2,3)], file="diff/miR277_medianThresh.diff", row.names=F, col.names=F, sep="\t", quote=F)
write.table(tT[tT$AveExpr >= median(tT$AveExpr),c(1,2,3)], file="diff/miR14_medianThresh.diff", row.names=F, col.names=F, sep="\t", quote=F)
#write.table(tT[tT$AveExpr >= median(tT$AveExpr),c(1,2,3)], file="diff/miR34_Day3_medianThresh.diff", row.names=F, col.names=F, sep="\t", quote=F)
