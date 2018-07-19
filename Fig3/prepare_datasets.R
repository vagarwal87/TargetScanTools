file = args[1]
filenum = args[2]
fullset = read.delim(file, colClasses=c("numeric", rep("character",2), rep("factor",5), rep("numeric",19), "character", "numeric", rep("character",2)))

seeds = unique(fullset$M8SEED)

train=data.frame()
test=data.frame()

mu=mean(fullset$FC)
stdev=sd(fullset$FC)

lowerQ=quantile(fullset$FC, 0.25)
diff=quantile(fullset$FC, 0.75)-quantile(fullset$FC, 0.25)

for(x in seeds){
	b=fullset$FC[fullset$M8SEED == x]
	fullset$FC[fullset$M8SEED == x] = scale(b, center=quantile(b, 0.25),scale=quantile(b, 0.75)-quantile(b, 0.25))*diff + lowerQ

	for(types in c("8m","7m8","7A1")){
		all=fullset[fullset$M8SEED==x & fullset$TYPE==types,]
		indices = sample(1:nrow(all), .7*nrow(all), replace=F)
		train = rbind(train, all[indices,])
		test = rbind(test, all[-indices,])
	}
}

save(train, test, fullset, file = paste("data/fly8mer7mer-90-", filenum, ".RData", sep=''))
