library(MASS)
library("gplots")

load("../Fig2/fly6_07.three_prime_UTR.RData")

calcPct <- function(site, bls) {
	if(!(site %in% names(pctparams))) return(NA);
	coefs = unlist(pctparams[[site]])
	max(0, coefs[1] + (coefs[2] / (1 + exp(-coefs[3]*bls + coefs[4]))))
}

preprocessdata <- function(data){
	data = data[complete.cases(data),]
	data[, c("LEN_3UTR", "LEN_ORF")] = log10(data[, c("LEN_3UTR", "LEN_ORF")]+1) 
	data$OTHERSITES = data$ORF8+data$UTR5P8+data$UTR6+data$UTRO6+data$UTR6A1
	data$PCT = sapply(1:nrow(data),function(x){ calcPct(as.character(data$SITE[x]), data$BLS[x]) })
	data$PCT_unscaled = data$PCT
	data[,c("FC","GENE_ID","MIR_ID","TYPE","BLS","HYBRIDSCORE","LEN_3UTR","LEN_ORF","PLFOLD","M8SEED","SITE","MIR","OTHERSITES","PCT", "PCT_unscaled")]
}

file = args[1]

load("../Fig3/data/fly8mer7mer-90-1.RData")
#load("../Fig3/data/fly6mers-90-1.RData") #use this for 6merOnly model

train=fullset
colnames(train)[1]="FC"
test = read.delim(file, colClasses=c("numeric", rep("character",1), rep("factor",2), rep("numeric",9), rep("character",1), "numeric", rep("character",2)))

train = preprocessdata(train)
test = preprocessdata(test)

predictions = data.frame()
x=1

for (types in c("8m","7m8","7A1")){ #c("6m","6A1","o6") used instead for 6merOnly model
	say(types)
	site.train = subset(train, TYPE %in% types)
	site.test = subset(test, TYPE %in% types)

	if(nrow(site.test)==0) next

	# add  '& colnames(site.train) != "PCT"' if 6merOnly model
	numeric = sapply(site.train, is.numeric) & colnames(site.train) != "FC" & colnames(site.train) != "BLS" & colnames(site.train) != "PCT_unscaled"
	test.numeric = sapply(site.test, is.numeric) & colnames(site.test) != "FC" & colnames(site.train) != "BLS" & colnames(site.test) != "PCT_unscaled"

	tmpmean = apply(site.train[,numeric], 2, function(x) quantile(x, 0.05))
	tmpsd = apply(site.train[,numeric], 2, function(x) quantile(x, 0.95)-quantile(x, 0.05))

	site.train[,numeric] = scale(site.train[,numeric], center=tmpmean, scale=tmpsd) # scale both sets in terms of training set points
	site.test[,test.numeric] = scale(site.test[,test.numeric], center=tmpmean, scale=tmpsd)

	site.train[,numeric][site.train[,numeric] > 1] = 1
	site.train[,numeric][site.train[,numeric] < 0] = 0
	site.test[,numeric][site.test[,numeric] > 1] = 1
	site.test[,numeric][site.test[,numeric] < 0] = 0
	head(site.test)

	context.formula = formula(FC ~ HYBRIDSCORE + PLFOLD + LEN_ORF + LEN_3UTR + OTHERSITES + PCT) # remove  "+ OTHERSITES + PCT" if 6merOnly model

#	# USE FOR KNOCKOUTS SET ONLY
#	lm.context = lm(context.formula, data = site.train)
#	site.test$CONTEXT = predict(lm.context, newdata = site.test)

	## USE FOR TRANSFECTION SET ONLY -- cross validation procedure
	site.test$CONTEXT = rep(0,nrow(site.test))
	for(miR in as.factor(site.train$MIR_ID) ){
		lm.context = lm(context.formula, data = site.train[site.train$MIR_ID != miR,])
		site.test$CONTEXT[site.test$MIR_ID == miR] = predict(lm.context, newdata = site.test[site.test$MIR_ID == miR,])
	}

	if	  (types=="8m")  site.test$CONTEXT[site.test$CONTEXT > -0.03] = -0.03
	else if (types=="7m8") site.test$CONTEXT[site.test$CONTEXT > -0.02] = -0.02
	else if (types=="7A1") site.test$CONTEXT[site.test$CONTEXT > -0.01] = -0.01
	else site.test$CONTEXT[site.test$CONTEXT > 0] = 0

	#conservation thresholds chosen because they correspond to approximate branch length cutoffs in which S:B ratio is 2:1
	site.test$CONS = rep(0,nrow(site.test))
	if	  (types=="8m")  site.test$CONS[site.test$BLS > 1] = 1
	else if (types=="7m8") site.test$CONS[site.test$BLS > 1.6] = 1
	else if (types=="7A1") site.test$CONS[site.test$BLS > 1.6] = 1
	else site.test$CONS[site.test$BLS > 4] = 1

	pred = site.test[,c("GENE_ID", "MIR_ID", "FC", "CONTEXT", "CONS", "PCT_unscaled")]

	if(x==1){ predictions = pred }
	else { predictions = rbind(predictions, pred) }
	x=x+1
}

predictions$NMID = predictions$GENE_ID
predictions$WEIGHT <- NULL
predictions$V4 = paste(predictions$NMID, "|", predictions$MIR_ID, sep='')
pred = predictions[,c("V4","CONTEXT","CONS", "PCT_unscaled")]
mypred = aggregate(pred[,"CONTEXT"],by=list(as.factor(pred$V4)),FUN=sum, na.rm=TRUE)
mypred2 = aggregate(pred[,"CONS"],by=list(as.factor(pred$V4)),FUN=max, na.rm=TRUE)
mypcts = aggregate(pred[,"PCT_unscaled"],by=list(as.factor(pred$V4)),FUN=function(x) { 1-prod(1-x) })
mypred$x = round(mypred$x, 3)
mypcts$x = round(mypcts$x, 3)
mypred = cbind(mypred,mypred2[,2])

write.table(mypred,file="context_8mer7m87A1_transfections_CrossVal.out",quote=F,row.names=F,sep="\t",col.names=F)
#write.table(mypred,file="context_8mer7m87A1_transfections_CrossVal_shortestUTR.out",quote=F,row.names=F,sep="\t",col.names=F)
#write.table(mypred,file="context_8mer7m87A1_knockouts.out",quote=F,row.names=F,sep="\t",col.names=F)
#write.table(mypred,file="context_8mer7m87A1_knockouts_shortestUTR.out",quote=F,row.names=F,sep="\t",col.names=F)
#write.table(mypred,file="context_6mero66A1_knockouts.out",quote=F,row.names=F,sep="\t",col.names=F)
#write.table(mypred,file="context_6mero66A1_knockouts_shortestUTR.out",quote=F,row.names=F,sep="\t",col.names=F)
#write.table(mypred,file="context_6mero66A1_transfections_CrossVal.out",quote=F,row.names=F,sep="\t",col.names=F)
#write.table(mypred,file="context_6mero66A1_transfections_CrossVal_shortestUTR.out",quote=F,row.names=F,sep="\t",col.names=F)

write.table(mypcts,file="aggregate_pct_8mer7m87A1_transfections.out",quote=F,row.names=F,sep="\t",col.names=F)
#write.table(mypcts,file="aggregate_pct_8mer7m87A1_transfections_shortestUTR.out",quote=F,row.names=F,sep="\t",col.names=F)
#write.table(mypcts,file="aggregate_pct_8mer7m87A1_knockouts.out",quote=F,row.names=F,sep="\t",col.names=F)
#write.table(mypcts,file="aggregate_pct_8mer7m87A1_knockouts_shortestUTR.out",quote=F,row.names=F,sep="\t",col.names=F)
