library(data.table)
library(reshape)
library(robustbase)
library(R.matlab)
library(seqinr)

Matlab$startServer()
matlab <- Matlab(host="localhost")
if (!open(matlab)) throw("Matlab server is not running: waited 30 seconds.")

### This Matlab function was written by Robin Friedman
setFunction(matlab, " \
function [estimates, model]=fitcurve(xdata, ydata) \
model = @expfun; \
start_point = [-0.5 1.5 3 2]; \
estimates = fminsearch(model, start_point); \
    function [sse, FittedCurve]=expfun(params) \
        b_0 = params(1); \
        b_1 = params(2); \
        b_2 = params(3); \
        b_3 = params(4); \
        FittedCurve = max(0,b_0 + b_1.* (1 ./ (1 + exp(-b_2.*xdata+ b_3)))); \
        ErrorVector = FittedCurve - ydata; \
        sse = sum(ErrorVector .^ 2) + 10 .*  (FittedCurve(1)+FittedCurve(2)); \
    end \
end \
");

species = "fly6_07"
region = args[1]
maxbl = 4
mincounts=20
dir = paste(species, "/", region, sep='')

files = list.files(path=paste(dir, "/s2b.REV", sep=''), pattern='.s2b', full.names=T)
files = files[sapply(files, FUN=function(file){ file.info(file)$size > 0 })]

s2btable <- do.call("rbind", lapply(files, FUN=function(file){ 
	tmp = read.table(textConnection(system(paste("head -n -1", file), intern=TRUE)), fill=T, header=F)
	tmp = tmp[,c(1,3)]
	mat = do.call(rbind, strsplit(as.character(tmp$V3), split = '-'))
	class(mat) <- "numeric"
	tmp = cbind(tmp$V1, data.frame(mat))
	colnames(tmp) = c('motif','sig','back')
	tmp$file = file
	tmp
}))

s2btable$bls = sapply(s2btable$file, function(x) { as.numeric(sub(".*?([0-9]+\\.[0-9]+).*", "\\1", x, perl=TRUE)) } )
s2btable$sitetype = sapply(s2btable$file, function(x) { as.character(sub(".*\\.(o?[6-8][m|A]?[1|8]?).*", "\\1", x, perl=TRUE)) } )
s2btable$file = NULL
s2btable = s2btable[order(s2btable$bls),]
head(s2btable,10)

kmers = list()
sig = list()
back = list()
s2b = list()

nested_correction <- function(counts, type1, type2, kmernames){
	for (i in 1:nrow(counts[[type1]])){
		idxs = which(rownames(counts[[type1]][i,])==kmernames)
		counts[[type1]][i,] = counts[[type1]][i,] - colSums(as.data.frame(counts[[type2]][idxs,]))
	}
	counts
}

for (type in c("8m", "7m8", "7A1", "6m" , "o6", "6A1")){
	say(type)
	sig[[type]] = cast(s2btable[s2btable$sitetype==type,], motif ~ bls, value="sig")
	back[[type]] = cast(s2btable[s2btable$sitetype==type,], motif ~ bls, value="back")

	rownames(sig[[type]]) = sig[[type]][,1]
	rownames(back[[type]]) = back[[type]][,1]
	sig[[type]][,1] = NULL
	back[[type]][,1] = NULL

	N = ncol(sig[[type]])
	for (i in 1:(N-1)) sig[[type]][,i] = sig[[type]][,i] - sig[[type]][,i+1]
	for (i in 1:(N-1)) back[[type]][,i] = back[[type]][,i] - back[[type]][,i+1]

	if (type == "8m"){
		kmers[[type]] = data.frame("S7m8" = substring(rownames(sig[[type]]), 1, 7), "S7A1" = substring(rownames(sig[[type]]), 2, 8),
					    "S6m" = substring(rownames(sig[[type]]), 2, 7),  "So6" = substring(rownames(sig[[type]]), 1, 6), "S6A1" = substring(rownames(sig[[type]]), 3, 8))
	}
	else if (type == "7m8"){
		kmers[[type]] = data.frame("S6m"= substring(rownames(sig[[type]]), 2, 7), "So6" = substring(rownames(sig[[type]]), 1, 6))
		sig = nested_correction(sig, type, "8m", kmers[["8m"]]$"S7m8")
		back = nested_correction(back, type, "8m", kmers[["8m"]]$"S7m8")
	}
	else if (type == "7A1"){
		kmers[[type]] = data.frame("S6m"= substring(rownames(sig[[type]]), 1, 6), "S6A1" = substring(rownames(sig[[type]]), 2, 7))
		sig = nested_correction(sig, type, "8m", kmers[["8m"]]$"S7A1")
		back = nested_correction(back, type, "8m", kmers[["8m"]]$"S7A1")
	}
	else if (type == "6m"){
		sig = nested_correction(sig, type, "8m", kmers[["8m"]]$"S6m")
		back = nested_correction(back, type, "8m", kmers[["8m"]]$"S6m")
		sig = nested_correction(sig, type, "7m8", kmers[["7m8"]]$"S6m")
		back = nested_correction(back, type, "7m8", kmers[["7m8"]]$"S6m")
		sig = nested_correction(sig, type, "7A1", kmers[["7A1"]]$"S6m")
		back = nested_correction(back, type, "7A1", kmers[["7A1"]]$"S6m")
	}
	else if (type == "o6"){
		sig = nested_correction(sig, type, "8m", kmers[["8m"]]$"So6")
		back = nested_correction(back, type, "8m", kmers[["8m"]]$"So6")
		sig = nested_correction(sig, type, "7m8", kmers[["7m8"]]$"So6")
		back = nested_correction(back, type, "7m8", kmers[["7m8"]]$"So6")
	}
	else if (type == "6A1"){
		sig = nested_correction(sig, type, "8m", kmers[["8m"]]$"S6A1")
		back = nested_correction(back, type, "8m", kmers[["8m"]]$"S6A1")
		sig = nested_correction(sig, type, "7A1", kmers[["7A1"]]$"S6A1")
		back = nested_correction(back, type, "7A1", kmers[["7A1"]]$"S6A1")
	}
}

for (type in c("8m", "7m8", "7A1", "6m", "o6", "6A1")){
	say(type)

	tmpsig = sig[[type]]
	tmpback = back[[type]]
	for (i in 1:nrow(tmpsig)){
		lowcounts = which(tmpsig[i,] < mincounts)
		for (j in lowcounts){
			breakFlag = FALSE
			for (window in 1:N){ #for smoothing
				tot = sum(tmpsig[i, max(1,j-window):min(N,j+window)])
				if (tot >= mincounts){
					breakFlag = TRUE
					sig[[type]][i,j] = tot
					back[[type]][i,j] = sum(tmpback[i, max(1,j-window):min(N,j+window)])
					break
				}
			}
			if (!breakFlag){
				sig[[type]][i,] = rep(1,N)
				back[[type]][i,] = rep(1,N)
			}
		}
	}
}

save(s2b, sig, back, kmers, file=paste(species, ".", region, ".RData", sep=''))
load(paste(species, ".", region, ".RData", sep=''))

pctparams = list()
Pct = function(coefs, bls) { coefs[1] + (coefs[2] / (1 + exp(-coefs[3]*bls + coefs[4]))) }

# correct Pct to make sure larger sites always have higher Pct than their nested smaller sites
pct_correction <- function(pcts, type1, type2, kmernames){
	for (i in 1:nrow(pcts[[type1]])){
		idxs = which(rownames(pcts[[type2]]) == kmernames[i])
		pcts[[type1]][i,] = apply(rbind(pcts[[type1]][i,], pcts[[type2]][idxs,]), 2, max)
	}
	pcts
}

pdf("Pct_plots_by_seed.pdf")
for (type in c("8m", "7m8", "7A1", "6m", "o6","6A1")){
	say(type)
	s2b[[type]] = (sig[[type]]) / (back[[type]])
	s2b[[type]] = (s2b[[type]]-1) / s2b[[type]]
	if (type == "8m"){
		s2b = pct_correction(s2b, type, "7m8", kmers[["8m"]]$"S7m8")
		s2b = pct_correction(s2b, type, "7A1", kmers[["8m"]]$"S7A1")
		s2b = pct_correction(s2b, type, "6m", kmers[["8m"]]$"S6m")
		s2b = pct_correction(s2b, type, "o6", kmers[["8m"]]$"So6")
		s2b = pct_correction(s2b, type, "6A1", kmers[["8m"]]$"S6A1")
	}
	else if (type == "7m8"){
		s2b = pct_correction(s2b, type, "6m", kmers[["7m8"]]$"S6m")
		s2b = pct_correction(s2b, type, "o6", kmers[["7m8"]]$"So6")
	}
	else if (type == "7A1"){
		s2b = pct_correction(s2b, type, "6m", kmers[["7A1"]]$"S6m")
		s2b = pct_correction(s2b, type, "6A1", kmers[["7A1"]]$"S6A1")
	}
	s2b[[type]][s2b[[type]] < 0] = 0
	s2b[[type]][s2b[[type]] > 1] = 1
	tmps2b = s2b[[type]]
	N=ncol(tmps2b)
	s2b[[type]] = tmps2b
	tmps2b = tmps2b[,1:length(seq(0,maxbl,0.05))]
	allx = as.numeric(colnames(tmps2b))

	for (i in 1:nrow(tmps2b)){
		say(motif <- rownames(tmps2b[i,]))
		minpos = 1
		pos = ncol(tmps2b)
		x1 = as.numeric(colnames(tmps2b)[1:pos])
		y1 = unlist(tmps2b[i,1:pos])

		setVariable(matlab, x1=x1, y1=y1)
		evaluate(matlab, "[x,y]=fitcurve(x1,y1);")
		matlabModel = unlist(getVariable(matlab, "x"))
		y2 = Pct(matlabModel, allx)
		pctparams[[motif]] = data.frame("B_0" = matlabModel[1], "B_1" = matlabModel[2], "B_2" = matlabModel[3], "B_3" = matlabModel[4])

		plot(allx, tmps2b[i,], type = "l", xlim = c(0,maxbl), ylim=c(0,1), main=motif)
		lines(allx, y2, type = "l", col="green")
	}
}

save(s2b, sig, back, kmers, pctparams, file=paste(species, ".", region, ".RData", sep=''))
close(matlab)
