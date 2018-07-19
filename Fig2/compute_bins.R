(dir = dirname(args))
a = read.table(args, header=F, sep="\t")

#i.e. fly6_07/five_prime_UTR/bins.txt

a=a[order(a$V1),]
head(a)

numbins = 5

(splits = quantile(a$V4, seq(0, 1, len = numbins+1)))

a$bls.bin <- cut(a$V4, splits, as.character(1:numbins), include.lowest = TRUE)

a$V4 = round(a$V4, 3)

write.table(a, file=paste(dir,"/allgenes.bins",sep=''), quote=F, row.names=F, sep="\t", col.names=F)
