require(graphics)

for (i in args){
	miR = strsplit(strsplit(i,"\\/")[[1]][2],'_')[[1]][1]
	say("Processing...", miR)
	a = read.delim(i, header=T, sep="\t")

	cuff = read.delim(paste("cuffdiff/", miR, "/gene_exp.diff", sep=''), header=T, sep="\t")
	a=merge(a,cuff,by=1,all.x=T)
	g = a[a$utr8 == 0 & a$utr7m8 == 0 & a$utr7A1 == 0 & a$utr6 == 0 & a$utro6 == 0 & a$utr6A1 == 0 & a$utr5p8 == 0 & a$orf8 == 0, ]
#	model = lm(fc~log10(utr5_len+1)+utr5_gc+log10(orf_len+1)+orf_gc+log10(utr3_len+1)+utr3_gc+log10(value_1), data=g)  #results in nearly identical outcome
	model = lm(fc~utr5_len+utr5_gc+orf_len+orf_gc+utr3_len+utr3_gc+log10(value_1), data=g)
	print(summary(model))
	a$fc = a$fc-predict(model, newdata=a)
	write.table(a, file=paste("targetsites_normalized/",miR,"_sites.out", sep=''), quote=F, row.names=F, sep="\t")
}
