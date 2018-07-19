id2geneid = read.delim("id2geneid.txt", F)
id2geneid = id2geneid[,c(1,3)]
colnames(id2geneid) = c("uniqid", "geneid")
say("total genes with annotated UTR =", (nGenes = length(unique(id2geneid$geneid))))

getsites = function(region){
	files = list.files(path=paste("targetsites/", region, sep=''), pattern='.out', full.names=T)

	sites <- do.call("rbind", lapply(files, FUN=function(file){
		tmp = read.table(file);
		tmp[,c(1,2,5)]
	}))

	colnames(sites)=c("uniqid","pos","kmer")
	say("total sites for", region, "=", nrow(sites))
	sites = unique(sites) #account for possibility different miRNA families share same 6mer or 7mer at same position
	say("total unique sites for", region, "=", nrow(sites))
	merge(sites, id2geneid, by="uniqid")
}

sites5utr = getsites("five_prime_UTR")
sites3utr = getsites("three_prime_UTR")

set.seed(124)
N = 1000
utr5Vals = round(rnorm(N, 840, 40.44)) #total mean = 325+165+350; total STD = sqrt(726.9624+793.1087+115.3815)
utr3Vals = round(rnorm(N, 12285, 213.81)) #total mean = 2128+4062+2837+2738+520; total STD = sqrt(18022.9226+3672.9610+1715.8166+358.1926+21944.5926)
#### generate 1000 random samples
#### sample ~12,285 3' UTR sites conserved above background and ~840 5' UTR sites, then find # of unique genes hit
genesHit = sapply(1:N, function(x){
	length(unique(c(sites5utr[sample(1:nrow(sites5utr), utr5Vals[x], replace=F),"geneid"],
	sites3utr[sample(1:nrow(sites3utr), utr3Vals[x], replace=F),"geneid"])))
})

"median"
as.integer(median(genesHit))
"lower 5th %, empirical"
as.integer(quantile(genesHit, 0.05))
"upper 95th %, empirical"
as.integer(quantile(genesHit, 0.95))

"mean"
as.integer(mean(genesHit))
as.integer(mean(genesHit))/nGenes
"lower 5th %, estimated from Gaussian"
as.integer(mean(genesHit) + qnorm(0.05)*sd(genesHit))
as.integer(mean(genesHit) + qnorm(0.05)*sd(genesHit))/nGenes
"upper 95th %, estimated from Gaussian"
as.integer(mean(genesHit) + qnorm(0.95)*sd(genesHit))
as.integer(mean(genesHit) + qnorm(0.95)*sd(genesHit))/nGenes
