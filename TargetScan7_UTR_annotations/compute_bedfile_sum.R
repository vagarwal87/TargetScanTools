a=read.delim(gzfile("3Pseq_pA_allstages.bed.gz"),header=F,sep='\t')
a=a[,c(1:4,6)]
head(a)
b = aggregate(a$V4,by=list(as.factor(a$V1),a$V2,a$V3,as.factor(a$V6)),FUN=sum)
write.table(b,file="3Pseq_pA_allstages.summed",quote=F,sep="\t",col.names=F,row.names=F)
