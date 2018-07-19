library(MASS)
library("gplots")
library(ROCR)

load("../Fig2/fly6_07.three_prime_UTR.RData")

calcPct <- function(site, bls) {
	if(!(site %in% names(pctparams))) return(0);
	coefs = unlist(pctparams[[site]])
	max(0, coefs[1] + (coefs[2] / (1 + exp(-coefs[3]*bls + coefs[4]))))
}

preprocessdata <- function(data){
	data = data[complete.cases(data),]
	data[, c("LEN_3UTR", "LEN_ORF", "MIN_DIST")] = log10(data[, c("LEN_3UTR", "LEN_ORF", "MIN_DIST")]+1)
	data$PCT = sapply(1:nrow(data),function(x){ calcPct(as.character(data$SITE[x]), data$BLS[x]) })
	data$OTHERSITES = data$ORF8 + data$UTR5P8 + data$UTR6 + data$UTRO6 + data$UTR6A1
	data
}

say("FLY")
#load("data/fly8mer7mer-90-1.RData")
fullset=read.delim("feature_table_7_8mers_only.txt2")
train=fullset
train = preprocessdata(train)

predictions = data.frame()
mycoefs = data.frame()
mystderr = data.frame()
x=1

for (types in list("8m","7m8","7A1")){
	say(types)
	site.train = subset(train, TYPE %in% types)

	numeric = sapply(site.train, is.numeric) & colnames(site.train) != "M8SEED" & colnames(site.train) != "FC"

	say("5th percentile values, values shown in Table 2")
	print(apply(site.train[,numeric], 2, function(x) quantile(x, 0.05)))
	say("95th percentile values, values shown in Table 2")
	print(apply(site.train[,numeric], 2, function(x) quantile(x, 0.95)))

	tmpmean = apply(site.train[,numeric], 2, function(x) quantile(x, 0.05))
	tmpsd = apply(site.train[,numeric], 2, function(x) quantile(x, 0.95)-quantile(x, 0.05))

	site.train[,numeric] = scale(site.train[,numeric], center=tmpmean, scale=tmpsd) # scale both sets in terms of training set points

	site.train[,numeric][site.train[,numeric] > 1] = 1
	site.train[,numeric][site.train[,numeric] < 0] = 0

	context.formula = formula(FC ~ HYBRIDSCORE+PLFOLD+LEN_ORF+LEN_3UTR+OTHERSITES+PCT)

	lm.context = lm(context.formula, data = site.train)
	print(summary(lm.context))

	tmpdata = as.data.frame(coef(lm.context))
	tmpdata2 = as.data.frame(summary(lm.context)$coefficients[, 2])
	if(x==1){
		mycoefs = tmpdata;
		mystderr = tmpdata2;
		colnames(mycoefs) = types;
		colnames(mystderr) = types;
	}
	else {
		mycoefs[,x] = tmpdata[,1];
		mystderr = merge(mystderr,tmpdata2,by=0,all=T);
		rownames(mystderr) = mystderr[,1];
		mystderr[,1] = NULL;
		colnames(mycoefs)[x] = types;
		colnames(mystderr)[x] = types;
	}
	x=x+1
}

pdf("Fig3E2.pdf", width=10, height=8)
par(mar=c(10, 5, 5, 1))
mycoefs
mymat = as.matrix(t(mycoefs))
mystderr = 1.96*mystderr
(mymat2=mystderr[rownames(mycoefs),])
mymat2 = as.matrix(t(mymat2))
colnames(mymat) = c("(Intercept)","3p_energy","SA","Len_ORF","Len_3UTR","Other_sites","PCT")
bar<-barplot(mymat, beside = TRUE, las=2, col=c("purple","red","blue"), border=F, ylim=c(-2, 2))
error.bar(bar,mymat,mymat2, length=0.02,lwd=.2)
legend("bottomright", bg="white", bty="n", legend = c("8mer", "7mer-m8", "7mer-A1"), text.col = c("purple","red","blue"))
dev.off()
