library(RColorBrewer)

whichSites=args[2]

if(args[3]=='median'){
	median=TRUE
	description="median"
} else{
	median=FALSE
	description="mean"
}

if(args[1]=='transfections'){
	transfections=TRUE
	description=paste("transfections", description, sep='-')
	args=c("diff/miR124.diff","diff/miR1.diff","diff/miR263a.diff","diff/miR4.diff","diff/miR92a.diff","diff/miR994.diff")
} else{
	transfections=FALSE
	description=paste("knockouts", description, sep='-')
	args=c("diff/miR277_medianThresh.diff","diff/miR14_medianThresh.diff","diff/miR34_Day20_medianThresh.diff")
}

description=paste(description, whichSites, sep='-')

numDatasets=length(args)
diff = data.frame()
diffseed = data.frame()
diffseedfiltered = data.frame()

for(x in args){
	miR = strsplit(strsplit(x,"\\/")[[1]][2],'\\.')[[1]][1]
	say(miR)
	diffile = read.delim(x, header=F)
#	b=diffile[,2]
#	diffile[,2] = scale(b, center=quantile(b, 0.25),scale=2.5*(quantile(b, 0.75)-quantile(b, 0.25)))
	diffile[,2]=diffile[,2]-median(diffile[,2])
	sitesfile = read.delim(paste("targetsites/",miR,"_sites.out",sep=""))
	diffile2 = diffile[diffile[,1] %in% sitesfile[sitesfile$utr8 == 0 & sitesfile$utr7m8 == 0 & sitesfile$utr7A1 == 0, "utr_id"], ]
	if (whichSites=='NoSite'){ diffile2 = diffile[diffile[,1] %in% sitesfile[sitesfile$utr8 == 0 & sitesfile$utr7m8 == 0 & sitesfile$utr7A1 == 0 & sitesfile$utr6 == 0 & sitesfile$utro6 == 0 & sitesfile$utr6A1 == 0, "utr_id"], ] }
	diffile3 = diffile[diffile[,1] %in% sitesfile[sitesfile$utr8 > 0 | sitesfile$utr7m8 > 0 | sitesfile$utr7A1 > 0, "utr_id"], ]
	colnames(diffile) = c("GENE_ID", "FC")
	colnames(diffile2) = c("GENE_ID", "FC")
	colnames(diffile3) = c("GENE_ID", "FC")
	say(nrow(diffile), nrow(diffile2), nrow(diffile3), median(diffile$FC))
	diffile$MIR_ID = rep(miR, nrow(diffile))
	diffile2$MIR_ID = rep(miR, nrow(diffile2))
	diffile3$MIR_ID = rep(miR, nrow(diffile3))
	diff = rbind(diff, diffile)
	diffseedfiltered = rbind(diffseedfiltered, diffile2)
	diffseed = rbind(diffseed, diffile3)
}

say(median(diff$FC))

preprocessdata <- function(data){
	data$MIR_ID[data$MIR_ID=="miR14_medianThresh"] = 'dme-miR-14'
	data$MIR_ID[data$MIR_ID=="miR34_Day20_medianThresh"] = 'dme-miR-34'
	data$MIR_ID[data$MIR_ID=="miR277_medianThresh"] = 'dme-miR-277'
	data$MIR_ID[data$MIR_ID=="miR124"] = 'dme-miR-124'
	data$MIR_ID[data$MIR_ID=="miR1"] = 'dme-miR-1'
	data$MIR_ID[data$MIR_ID=="miR263a"] = 'dme-miR-263a'
	data$MIR_ID[data$MIR_ID=="miR4"] = 'dme-miR-4'
	data$MIR_ID[data$MIR_ID=="miR994"] = 'dme-miR-994'
	data$MIR_ID[data$MIR_ID=="miR92a"] = 'dme-miR-92a'
	data$GENE_ID = paste(data$GENE_ID, "|", data$MIR_ID, sep='')
	data$MIR_ID = NULL
	data
}

diff=preprocessdata(diff)
diffseed=preprocessdata(diffseed)
diffseedfiltered=preprocessdata(diffseedfiltered)

if(transfections){
	a=read.delim("context_8mer7m87A1_transfections_CrossVal.out", header=F)
	b=read.delim("context_8mer7m87A1_transfections_CrossVal_shortestUTR.out", header=F)
	if (whichSites=='6merOnly'){
		a=read.delim("context_6mero66A1_transfections_CrossVal.out", header=F)
		b=read.delim("context_6mero66A1_transfections_CrossVal_shortestUTR.out", header=F)
	}
} else{
	a=read.delim("context_8mer7m87A1_knockouts.out", header=F)
	b=read.delim("context_8mer7m87A1_knockouts_shortestUTR.out", header=F)
	if (whichSites=='6merOnly'){
		a=read.delim("context_6mero66A1_knockouts.out", header=F)
		b=read.delim("context_6mero66A1_knockouts_shortestUTR.out", header=F)
	}
}
a=a[,1:2]
b=b[,1:3]
a=merge(a,b,by=1,all=TRUE)
a[is.na(a)]=0
if(transfections){
	b=read.delim("aggregate_pct_8mer7m87A1_transfections.out", header=F)
	c=read.delim("aggregate_pct_8mer7m87A1_transfections_shortestUTR.out", header=F)
} else{
	b=read.delim("aggregate_pct_8mer7m87A1_knockouts.out", header=F)
	c=read.delim("aggregate_pct_8mer7m87A1_knockouts_shortestUTR.out", header=F)
}
b=merge(b,c,by=1,all=TRUE)
b[is.na(b)]=0
b[,2]=(b[,2]+b[,3])/2
b[,3]=NULL

a=merge(a,b,by=1,all=TRUE)
colnames(a)=c("GENE_ID", "context_Longest_UTR", "context_Shortest_UTR", "Cons", "aggregate_pct") #"context8m7m87A1", "context8m", "context8m7m8", , "context8m7m87A16m", "context8m7m87A16mo6", "context8m7m87A16mo66a1"
a$context=0.5*a$context_Longest_UTR+0.5*a$context_Shortest_UTR

a$GENE_ID=gsub("miR14_medianThresh", "dme-miR-14", a$GENE_ID)
a$GENE_ID=gsub("miR34_Day20_medianThresh", "dme-miR-34", a$GENE_ID)
a$GENE_ID=gsub("miR277_medianThresh", "dme-miR-277", a$GENE_ID)
a$GENE_ID=gsub("miR1", "dme-miR-1", a$GENE_ID)
a$GENE_ID=gsub("miR2", "dme-miR-2", a$GENE_ID)
a$GENE_ID=gsub("miR4", "dme-miR-4", a$GENE_ID)
a$GENE_ID=gsub("miR9", "dme-miR-9", a$GENE_ID)

if(transfections){ 
	b=read.delim("all.dme.transfections.txt")
} else{
	b=read.delim("all.dme.knockouts.txt")
}
c=merge(b,a,by=1,all=TRUE)
d=merge(diffseedfiltered,c,by=1,all.x=TRUE) #for mRNAs filtered by canonical sites
z=merge(diffseed,c,by=1,all.x=TRUE) #for 7/8mer targets
c=merge(diff,c,by=1,all.x=TRUE) #for all targets

if (whichSites=='6merOnly' | whichSites=='NoSite'){ c=d } #filters out mRNAs with canonical sites
c=c[,c("GENE_ID","FC","ComiR","DianaMicroTCDS","EIMMO","EMBL","MicroCosm","miRSVR","PicTar","PITA","PITAtop","RNAhybrid","RNA22","TargetSpy","TargetScan", "MinoTar", "context", "aggregate_pct")]
z=z[,c("GENE_ID","FC","Cons","ComiR","DianaMicroTCDS","EIMMO","EMBL","MicroCosm","miRSVR","PicTar","PITA","PITAtop","RNAhybrid","RNA22","TargetSpy","TargetScan", "MinoTar", "context", "aggregate_pct")]

c$miRSVR=-c$miRSVR
c$PITA=-c$PITA
c$PITAtop=-c$PITAtop
c$MicroCosm=-c$MicroCosm
c$RNAhybrid=-c$RNAhybrid
c$context=-c$context

pdf(paste("Fig4-", description, ".pdf", sep=''), width=8, height=10)
par(mar=c(12, 7, 5, 5), mgp = c(5, 1, 0))
set.seed(30)
crp.rg <- colorRampPalette(c("red","orange","green","cyan","blue","purple","magenta"))
cols <- sample(crp.rg(ncol(c)-2))

if(transfections){
	ymin=-1.4
	ymax=0.2
	if(median){ ymin=-1.3 }
} else{
	ymin=-0.2
	ymax=1.1
}

plot(c(4,4^6), c(ymin,ymax), type="n", xlab="Average top predictions considered (per miRNA)", ylab="Mean of mean/median mRNA change",log='x', xaxt="n", yaxt="n", bty="n")
xpos=c(30,100,300,1000,3000,10000,30000)
xpos=4^seq(1,6,0.5)
ypos=round(seq(ymin,ymax,0.1),2)
axis(side=3,at=xpos, cex.axis=1, cex.lab=2)
axis(side=2,at=ypos, cex.axis=1, cex.lab=2, las=2)
set.seed(29)
yvals=sapply(xpos, function(x){ quantile(sapply(1:1000, function(y) { median(sample(diff$FC,x*numDatasets)) }), c(0.025, 0.975)) }) #pooled
polygon(c(xpos, rev(xpos)), c(yvals[1,], rev(yvals[2,])), col="lightgrey",border=F)
abline(h = 0)
set.seed(80)
plotchar <- sample(c(15:18,22:25), ncol(c)-2, replace=T)
c[,1]=as.factor(unlist(strsplit(as.character(c[,1]),'\\|'))[seq(2,nrow(c)*2,2)]) #by miRNA

allfcs = list()
yvals = list()
for (N in 3:(ncol(c))){
	if(sum(!is.na(c[,N])) > 10){
		xvals=4^seq(1,min((sum(!is.na(c[,N]))/numDatasets)^(1/4),6),0.5)
		say(colnames(c)[N], sum(!is.na(c[,N])))
		allfcs[[colnames(c)[N]]] = lapply(xvals, function(x){ unlist(by(c[,c(2,N)],c[,1],function(y){
			y[order(y[,2],decreasing=T)[1:x],1]
		})) } ) #save fold changes for given threshold into a list of lists
		if(median){
			yvals[[colnames(c)[N]]]=sapply(xvals, function(x){ mean(unlist(by(c[,c(2,N)],c[,1],function(y){ 
				if(is.na(y[order(y[,2],decreasing=T)[1],2])){ NA; }
				else { median(y[order(y[,2],decreasing=T)[1:x],1]); }
			})), na.rm=T) } ) #compute mean of median by miRNA
		} else{
			yvals[[colnames(c)[N]]]=sapply(xvals, function(x){ mean(unlist(by(c[,c(2,N)],c[,1],function(y){ 
				if(is.na(y[order(y[,2],decreasing=T)[1],2])){ NA; }
				else { mean(y[order(y[,2],decreasing=T)[1:x],1]); }
			})), na.rm=T) } ) #compute mean of mean by miRNA
		}
		lines(xvals, yvals[[colnames(c)[N]]], type='l', col=cols[N-2])
		lines(xvals, yvals[[colnames(c)[N]]], type='p', col=cols[N-2], pch=plotchar[N-2])
	}
}
legend("bottomright", colnames(c[3:(ncol(c))]), bg="white", text.col=cols, col = cols, pch=plotchar, lty = c(1, 1), lwd=c(1,1), cex=0.9, bty="n",ncol=2)

if (whichSites!='NoSite'){ 
	xvals=4^seq(1,6,0.5)
	algnames = names(yvals)
	sapply(1:length(yvals[['context']]), function(x){
		valsAtThresh = c()
		for (name in names(yvals)) {
			valsAtThresh = c(valsAtThresh, yvals[[name]][x])
		}
		say(x, algnames[rank(valsAtThresh)==2])
		if(transfections){
			if(wilcox.test(unlist(allfcs[["context"]][x]), unlist(allfcs[[algnames[rank(valsAtThresh)==2]]][[x]]), alternative="less")$p.value < 0.05) text(xvals[x], yvals[['context']][x]-0.05, '*')
		} else{
			if(wilcox.test(unlist(allfcs[["context"]][x]), unlist(allfcs[[algnames[rank(-1*valsAtThresh)==2]]][[x]]), alternative="greater")$p.value < 0.05) text(xvals[x], yvals[['context']][x]+0.05, '*')
		}
	})
}

#plot mean of mean/median for subset of mRNAs with 7/8mer sites
nrow(z)/numDatasets
if(median){
	MoM=mean(unlist(by(z[,c(2,3)],z[,1],function(y){ median(y[,"FC"]) })))
} else{
	MoM=mean(unlist(by(z[,c(2,3)],z[,1],function(y){ mean(y[,"FC"]) })))
}
abline(h = MoM, lwd=2)

#plot mean of mean/median for subset of mRNAs with conserved 7/8mer sites
z=z[z$Cons==1,]
nrow(z)/numDatasets
if(median){
	MoM=mean(unlist(by(z[,c(2,3)],z[,1],function(y){ median(y[,"FC"]) })))
} else{
	MoM=mean(unlist(by(z[,c(2,3)],z[,1],function(y){ mean(y[,"FC"]) })))
}
abline(h = MoM, lwd=2)
