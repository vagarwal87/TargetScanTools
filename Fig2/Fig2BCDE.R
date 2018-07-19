options(warn=-1)

region = args[1]

files = list.files(path=paste("fly6_07/", region, "/s2b.REV", sep=''), pattern='.s2b', full.names=T)

files = files[sapply(files, FUN=function(file){ file.info(file)$size > 0 })]

sbtable <- do.call("rbind", lapply(files, FUN=function(file){
	tmp = read.table(textConnection(system(paste("tail -1", file), intern=TRUE)));
	tmp$file = file
	tmp
}))

sbtable$bls = sapply(sbtable$file, function(x) { as.numeric(sub(".*?([0-9]+\\.[0-9]+).*", "\\1", x, perl=TRUE)) } )
sbtable$sitetype = sapply(sbtable$file, function(x) { as.character(sub(".*\\.(o?[6-8][m|A]?[1|8]?).*", "\\1", x, perl=TRUE)) } )
sbtable$file = NULL
sbtable = sbtable[order(sbtable$bls) & sbtable$bls <= 3,]

colnames(sbtable)[1:4]=c("sig","back","sab_var","s2b_var")

##########SUBTRACT NESTED SITES#############
sbtable[sbtable$sitetype == "7m8",c("sig","back")] = sbtable[sbtable$sitetype == "7m8",c("sig","back")]-sbtable[sbtable$sitetype == "8m",c("sig","back")]
sbtable[sbtable$sitetype == "7A1",c("sig","back")] = sbtable[sbtable$sitetype == "7A1",c("sig","back")]-sbtable[sbtable$sitetype == "8m",c("sig","back")]
sbtable[sbtable$sitetype == "6m",c("sig","back")] = sbtable[sbtable$sitetype == "6m",c("sig","back")]-sbtable[sbtable$sitetype == "7m8",c("sig","back")]
sbtable[sbtable$sitetype == "6m",c("sig","back")] = sbtable[sbtable$sitetype == "6m",c("sig","back")]-sbtable[sbtable$sitetype == "7A1",c("sig","back")]
sbtable[sbtable$sitetype == "6m",c("sig","back")] = sbtable[sbtable$sitetype == "6m",c("sig","back")]-sbtable[sbtable$sitetype == "8m",c("sig","back")]
sbtable[sbtable$sitetype == "o6",c("sig","back")] = sbtable[sbtable$sitetype == "o6",c("sig","back")]-sbtable[sbtable$sitetype == "7m8",c("sig","back")]
sbtable[sbtable$sitetype == "o6",c("sig","back")] = sbtable[sbtable$sitetype == "o6",c("sig","back")]-sbtable[sbtable$sitetype == "8m",c("sig","back")]
sbtable[sbtable$sitetype == "6A1",c("sig","back")] = sbtable[sbtable$sitetype == "6A1",c("sig","back")]-sbtable[sbtable$sitetype == "7A1",c("sig","back")]
sbtable[sbtable$sitetype == "6A1",c("sig","back")] = sbtable[sbtable$sitetype == "6A1",c("sig","back")]-sbtable[sbtable$sitetype == "8m",c("sig","back")]

sbtable$sab = sbtable$sig-sbtable$back #signal above background
sbtable$s2b = sbtable$sig/sbtable$back
sbtable$sablower = sbtable$sab + qnorm(0.05)*sqrt(sbtable$sab_var)
sbtable$s2blower = sbtable$s2b + qnorm(0.05)*sqrt(sbtable$s2b_var)

sbtable$consPercent = sbtable$sab/sbtable$sig
"Stats used to compute Fig2G fraction conserved at cutoff 1.0 for 5' UTRs and 3' UTRs as well as number of sites conserved at maximal sensitivity; mean and variance used in estimate_transcripts_targeted.R script"
sbtable[sbtable$bls == 1,]

pdf(paste("Fig2", ".", region, ".pdf", sep=''))
par(mar=c(6, 6, 4, 4))
plot(unique(sbtable$bls), sbtable$sab[sbtable$sitetype=="8m"], col="purple", type="l", lwd=3, xlab="Branch-length cutoff",
 ylab="Signal above background", bty = "n", cex.axis=2, cex.lab = 2, ylim=c(0,4000), xlim=c(0,2.5))
lines(unique(sbtable$bls), sbtable$sablower[sbtable$sitetype=="8m"], col="purple", type="l", lty=4, lwd=3)
lines(unique(sbtable$bls), sbtable$sab[sbtable$sitetype=="7m8"], col="red", type="l", lty=1, lwd=3)
lines(unique(sbtable$bls), sbtable$sablower[sbtable$sitetype=="7m8"], col="red", type="l", lty=4, lwd=3)
lines(unique(sbtable$bls), sbtable$sab[sbtable$sitetype=="7A1"], col="blue", type="l", lty=1, lwd=3)
lines(unique(sbtable$bls), sbtable$sablower[sbtable$sitetype=="7A1"], col="blue", type="l", lty=4, lwd=3)
lines(unique(sbtable$bls), sbtable$sab[sbtable$sitetype=="6m"], col="cyan", type="l", lty=1, lwd=3)
lines(unique(sbtable$bls), sbtable$sablower[sbtable$sitetype=="6m"], col="cyan", type="l", lty=4, lwd=3)
lines(unique(sbtable$bls), sbtable$sab[sbtable$sitetype=="o6"], col="grey", type="l", lty=1, lwd=3)
lines(unique(sbtable$bls), sbtable$sablower[sbtable$sitetype=="o6"], col="grey", type="l", lty=4, lwd=3)
lines(unique(sbtable$bls), sbtable$sab[sbtable$sitetype=="6A1"], col="darkorange", type="l", lty=1, lwd=3)
lines(unique(sbtable$bls), sbtable$sablower[sbtable$sitetype=="6A1"], col="darkorange", type="l", lty=4, lwd=3)
legend("topright", bty="n", cex=1.7, c("8mer", "7mer-m8", "7mer-A1", "6mer", "Offset 6mer", "6mer-A1"), text.col = c("purple","red","blue","cyan","grey", "darkorange"))
abline(0,0)

plot(unique(sbtable$bls), sbtable$s2b[sbtable$sitetype=="8m"], col="purple", type="l", lwd=3, xlab="Branch-length cutoff",
 ylab="Signal-to-background ratio", bty = "n", cex.axis=2, cex.lab = 2, log="y", c(0,2.5), ylim=c(0.75,5))
lines(unique(sbtable$bls), sbtable$s2blower[sbtable$sitetype=="8m"], col="purple", type="l", lty=4, lwd=3)
lines(unique(sbtable$bls), sbtable$s2b[sbtable$sitetype=="7m8"], col="red", type="l", lty=1, lwd=3)
lines(unique(sbtable$bls), sbtable$s2blower[sbtable$sitetype=="7m8"], col="red", type="l", lty=4, lwd=3)
lines(unique(sbtable$bls), sbtable$s2b[sbtable$sitetype=="7A1"], col="blue", type="l", lty=1, lwd=3)
lines(unique(sbtable$bls), sbtable$s2blower[sbtable$sitetype=="7A1"], col="blue", type="l", lty=4, lwd=3)
lines(unique(sbtable$bls), sbtable$s2b[sbtable$sitetype=="6m"], col="cyan", type="l", lty=1, lwd=3)
lines(unique(sbtable$bls), sbtable$s2blower[sbtable$sitetype=="6m"], col="cyan", type="l", lty=4, lwd=3)
lines(unique(sbtable$bls), sbtable$s2b[sbtable$sitetype=="o6"], col="grey", type="l", lty=1, lwd=3)
lines(unique(sbtable$bls), sbtable$s2blower[sbtable$sitetype=="o6"], col="grey", type="l", lty=4, lwd=3)
lines(unique(sbtable$bls), sbtable$s2b[sbtable$sitetype=="6A1"], col="darkorange", type="l", lty=1, lwd=3)
lines(unique(sbtable$bls), sbtable$s2blower[sbtable$sitetype=="6A1"], col="darkorange", type="l", lty=4, lwd=3)
axis(2,at=1:10)
abline(0,0)
legend("topleft", bty="n", cex=1.7, c("8mer", "7mer-m8", "7mer-A1", "6mer", "Offset 6mer", "6mer-A1"), text.col = c("purple","red","blue","cyan","grey","darkorange"))
