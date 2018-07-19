library(MASS)
library("gplots")

load("fly6_07.three_prime_UTR.RData")

calcPct <- function(site, bls) {
	if(!(site %in% names(pctparams))) return(NA);
	coefs = unlist(pctparams[[site]])
	max(0, coefs[1] + (coefs[2] / (1 + exp(-coefs[3]*bls + coefs[4]))))
}

preprocessdata <- function(data){
	data$PCT = sapply(1:nrow(data),function(x){ calcPct(data$SITE[x], data$BLS[x]) })
	data[,c("FC","MIR_ID","TYPE","PCT")]
}

file = args[1]

fullset = read.delim(file, colClasses=c("numeric", rep("character",2), rep("factor",5), rep("numeric",19), "character", "numeric", rep("character",2)))
fullset = preprocessdata(fullset)

head(fullset)
ymin=-1.4
ymax=0

pdf("Fig2H.pdf")
plot(c(0,1), c(ymin,ymax), type="n", xlab="Probability of conserved targeting (Pct)", ylab="mRNA change",xaxt="n", yaxt="n", bty="n")
xpos=seq(0,1,0.2)
ypos=round(seq(ymin,ymax,0.2),2)
axis(side=3,at=xpos, cex.axis=1, cex.lab=2)
axis(side=2,at=ypos, cex.axis=1, cex.lab=2)

N=0
cols=c("purple","red","blue")
for (types in list("8m","7m8","7A1")){
	N=N+1
	say(types)
	a = subset(fullset, TYPE %in% types)

	print(quantile(a$PCT, seq(0, 1, len = 6+1)))

	bins=6
	vals=data.frame(PCT=0,FC=0, serr=0)
	a=a[order(a$PCT,decreasing=T),]
	z=0
	for(x in seq(1,nrow(a),nrow(a)/bins)){
		z=z+1
		x=as.integer(x)
		y=min(x+as.integer(nrow(a)/bins), nrow(a))
		vals[z,]=c(mean(a[x:y,"PCT"]), mean(a[x:y,"FC"]), sd(a[x:y,"FC"])/sqrt(length(a[x:y,"FC"])))
	}

	print(summary(lm(FC ~ PCT,data=a)))
	lines(vals$PCT, vals$FC, type='p', col=cols[N], pch=16)
	error.bar(vals$PCT, vals$FC, vals$serr, col=cols[N])
	abline(lm(FC ~ PCT,data=vals),col=cols[N])
}
