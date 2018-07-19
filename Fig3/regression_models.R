library(leaps)
library(BMA)
library(MASS)
library(glmnet)
library(pls)
library(caret)
library(mboost)
library(modeltools)
library(earth)
library(randomForest)

load("../Fig2/fly6_07.three_prime_UTR.RData")

calcPct <- function(site, bls) {
	if(!(site %in% names(pctparams))) return(0);
	coefs = unlist(pctparams[[site]])
	max(0, coefs[1] + (coefs[2] / (1 + exp(-coefs[3]*bls + coefs[4]))))
}

preprocessdata <- function(data){
	data = data[complete.cases(data),fields]
	data$OTHERSITES = data$ORF8 + data$UTR5P8 + data$UTR6 + data$UTRO6 + data$UTR6A1
	data[, c("LEN_3UTR", "LEN_5UTR", "LEN_ORF", 'DIST5P', "MIN_DIST")] = log10(data[, c("LEN_3UTR", "LEN_5UTR", "LEN_ORF", 'DIST5P', "MIN_DIST")])
	data$PCT = sapply(1:nrow(data),function(x){ calcPct(as.character(data$SITE[x]), data$BLS[x]) })
	data
}

filenum = args[1]

say("FLY")
load(paste("data/fly8mer7mer-90-", filenum, ".RData", sep=""))
fields = c('FC', 'TYPE', 'TARNUC10', 'TARNUC9', 'TARNUC8', 'TARNUC1', 
	'ORF8', 'UTR5P8', 'SUPP_PAIR', 'HYBRIDSCORE', 'UTR6', 'UTRO6', 'UTR6A1',
	'LEN_3UTR', 'LEN_5UTR', 'LEN_ORF', 'DIST5P', 'MIN_DIST', 'LOCALAU_30', 
	'AU_3UTR', 'AU_5UTR', 'AU_ORF', 'PLFOLD', "BLS", "SITE")

train = preprocessdata(train)
test = preprocessdata(test)

site.binary.all = c()
lm.pred = c()
pred.context.plus = c()
pred.context.only = c()
step.pred = c()
glmnet.pred = c()
bicreg.pred = c()
mars.pred = c()

grouped.formula = formula(FC ~ TYPE + SUPP_PAIR + LEN_5UTR + LEN_ORF + LEN_3UTR + DIST5P + MIN_DIST 
	+ LOCALAU_30 + AU_3UTR + AU_5UTR + AU_ORF + PLFOLD + HYBRIDSCORE + OTHERSITES
	+ TARNUC10 + TARNUC9 + PCT)

site.train = train
site.test = test

train.matrix = model.matrix(grouped.formula, site.train)
test.matrix = model.matrix(grouped.formula, site.test)

dependencies = findLinearCombos(train.matrix) #uses QR decomposition to find linear dependencies in matrix (caret package)
train.matrix = train.matrix[, -c(dependencies$remove, 1)]
test.matrix = test.matrix[, -c(dependencies$remove, 1)]

say("traincounts", dim(site.train), "testcounts", dim(site.test))
numeric = sapply(site.train, is.numeric) & colnames(site.train) != "FC" #DO NOT SCALE FC to [0,1] interval
tmpmean = apply(site.train[,numeric], 2, function(x) quantile(x, 0.05))
tmpsd = apply(site.train[,numeric], 2, function(x) quantile(x, 0.95)-quantile(x, 0.05))

site.train[,numeric] = scale(site.train[,numeric], center=tmpmean, scale=tmpsd) # scale both sets in terms of training set points
site.test[,numeric] = scale(site.test[,numeric], center=tmpmean, scale=tmpsd)

site.train[,numeric][site.train[,numeric] > 1] = 1
site.train[,numeric][site.train[,numeric] < 0] = 0
site.test[,numeric][site.test[,numeric] > 1] = 1
site.test[,numeric][site.test[,numeric] < 0] = 0

numeric = apply(train.matrix, 2, function(x){length(unique(x))>2})
tmpmean = apply(train.matrix[,numeric], 2, function(x) quantile(x, 0.05))
tmpsd = apply(train.matrix[,numeric], 2, function(x) quantile(x, 0.95)-quantile(x, 0.05))
train.matrix[,numeric] = scale(train.matrix[,numeric], center=tmpmean, scale=tmpsd) # scale both sets in terms of training set points
test.matrix[,numeric] = scale(test.matrix[,numeric], center=tmpmean, scale=tmpsd)

train.matrix[,numeric][train.matrix[,numeric] > 1] = 1
train.matrix[,numeric][train.matrix[,numeric] < 0] = 0
test.matrix[,numeric][test.matrix[,numeric] > 1] = 1
test.matrix[,numeric][test.matrix[,numeric] < 0] = 0

say("site-only")
lm.context.only = lm(FC ~ TYPE, data = site.train)
pred.context.only = c(pred.context.only, predict(lm.context.only, newdata = site.test))
say(cor(site.test$FC, predict(lm.context.only, newdata = site.test), method="pearson"))

say("context-optimized")
lm.context.only = lm(FC ~ TYPE+LEN_ORF+LEN_3UTR+PCT+PLFOLD+OTHERSITES+HYBRIDSCORE, data = site.train)
pred.context.only = c(pred.context.only, predict(lm.context.only, newdata = site.test))
say(cor(site.test$FC, predict(lm.context.only, newdata = site.test), method="pearson"))

say("fwd-indiv")
sets.fit=regsubsets(grouped.formula, data = site.train, nvmax=15, nbest=1, method="forward", really.big=T) #from leaps package -- nvmax is max subset size // forward or backward or exhaustive
bestfeatures = (rank(summary(sets.fit)$bic)==1) #get index of model with minimum BIC
coefs = coef(sets.fit, which(bestfeatures)) #get coefficients of best model
indices = summary(sets.fit)$which[bestfeatures][-1] #get indices of chosen features of best model
tmptest = cbind(rep(1,nrow(test.matrix)), test.matrix[,indices])
vals = t(coefs %*% t(tmptest))
best.formula = paste(colnames(summary(sets.fit)$outmat)[indices], collapse="+") #extract model with minimum BIC
say(cor(site.test$FC, vals, method="pearson"))

say("stepBIC")
stepbic.fit = stepAIC(lm(FC ~ 1, data = site.train), scope = list( upper = grouped.formula, lower = ~1),
	direction="both", k = log(nrow(site.train)), trace = 0) #k=logN for BIC/ = 2 for AIC, k = log(nrow(site.train))
step.pred = c(step.pred, predict(stepbic.fit, newdata = site.test))
print(as.character(cor(site.test$FC, predict(stepbic.fit, newdata = site.test), method="pearson")))

say("stepAIC")
step.fit = stepAIC(lm(FC ~ 1, data = site.train), scope = list( upper = grouped.formula, lower = ~1),
	direction="both", k = 2, trace = 0) #k=logN for BIC/ = 2 for AIC, k = log(nrow(site.train))
step.pred = c(step.pred, predict(step.fit, newdata = site.test))
print(as.character(cor(site.test$FC, predict(step.fit, newdata = site.test), method="pearson")))

say("Lasso")
glmnet.fit = cv.glmnet(train.matrix, site.train$FC, family="gaussian", nfolds=10, alpha = 1)
print(as.character(cor(site.test$FC, predict(glmnet.fit, type="link", newx=test.matrix), method="pearson")))

say("MARS")
mars.fit = earth(grouped.formula, data = site.train, degree = 1, trace = 0, nk = 500)
say(cor(site.test$FC, predict(mars.fit, newdata = site.test), method="pearson"))

say("RandomForest")
rf.fit = randomForest(train.matrix, site.train$FC, importance = T)
say(cor(site.test$FC, predict(rf.fit, newdata = test.matrix), method="pearson"))

say("PCR")
pcr.fit = pcr(grouped.formula, data = site.train)
pcr.pred = predict(pcr.fit, ncomp = 5, newdata = site.test)
say(cor(site.test$FC, pcr.pred, method="pearson"))

say("PLSR")
plsr.fit = plsr(grouped.formula, data = site.train)
plsr.pred = predict(plsr.fit, ncomp = 5, newdata = site.test)
say(cor(site.test$FC, plsr.pred, method="pearson"))

say("fwd")
say(best.formula)
say("StepBIC")
say(formula(stepbic.fit)[3])
say("StepAIC")
say(formula(step.fit)[3])
