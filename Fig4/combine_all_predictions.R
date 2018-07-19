b = data.frame()
count = 1

for (file in args){
	say(file)
	a = read.delim(file, header=F)
	a[,3]=as.numeric(as.character(a[,3]))

#	#for knockouts use this
#	a=a[a[,2] %in% c("dme-miR-277","dme-miR-277-3p","dme-miR-34","dme-miR-34-5p","dme-miR-14","dme-miR-14-3p"),]

#	#for transfections use this
	a=a[a[,2] %in% c("dme-miR-1","dme-miR-124","dme-miR-4","dme-miR-92a","dme-miR-263a","dme-miR-994", "dme-miR-1-3p","dme-miR-124-3p","dme-miR-4-3p","dme-miR-92a-3p","dme-miR-263a-5p","dme-miR-994-5p"),]

	a$V2=gsub("-3p","",a$V2)
	a$V2=gsub("-5p","",a$V2)
	a$V4 = paste(a$V1, "|", a$V2, sep='')
	head(a$V4)
	a=a[,c(4,3)]
	colnames(a)[2] = strsplit(basename(file),"\\.")[[1]][1]
	if(count==1){ b = a }
	else b=merge(b,a,by=1,all=T)
	head(b)
	count=count+1
}

##for knockouts use this
#write.table(b,file="all.dme.testset.pred",quote=F,row.names=F,sep="\t",col.names=T)

##for transfections use this
write.table(b,file="all.dme.transfections.pred",quote=F,row.names=F,sep="\t",col.names=T)
