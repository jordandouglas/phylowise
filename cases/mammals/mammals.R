
library(phylowise)
library(randomForest)
library(ape)
library(corrplot)
library(BMA)




# Return a tree where the branch lengths are equal to this branch-metadata value
get.tree.with.lengths = function(tree, yvar="nS"){


	nbranches = nrow(tree$edge)


	# Total subst
	nT = as.numeric(sapply(strsplit(tree$node.comment, ","), function(ele) {
		gsub("[&]SubstitutionSum.species[=]", "", ele[grep("SubstitutionSum.species[=]", ele)])
	}))
	nT[is.na(nT)] = 0


	# Synonymous subst
	nS = as.numeric(sapply(strsplit(tree$node.comment, ","), function(ele) {
		gsub("^SynonymousSubstSum.species[=]", "", ele[grep("^SynonymousSubstSum.species[=]", ele)])
	}))
	nS[is.na(nS)] = 0


	# Non-synonymous subst
	nN = as.numeric(sapply(strsplit(tree$node.comment, ","), function(ele) {
		gsub("^NonSynonymousSubstSum.species[=]", "", ele[grep("^NonSynonymousSubstSum.species[=]", ele)])
	}))
	nN[is.na(nN)] = 0


		
	nT.branch = sapply(1:nbranches, function(edge.nr){
		parent = tree$edge[edge.nr,1]
		child = tree$edge[edge.nr,2]
		nT[child]
	})
	nS.branch = sapply(1:nbranches, function(edge.nr){
		parent = tree$edge[edge.nr,1]
		child = tree$edge[edge.nr,2]
		nS[child]
	})
	nN.branch = sapply(1:nbranches, function(edge.nr){
		parent = tree$edge[edge.nr,1]
		child = tree$edge[edge.nr,2]
		nN[child]
	})



	if (yvar == "nT"){
		tree$edge.length = nT.branch
	}


	if (yvar == "nS"){
		tree$edge.length = nS.branch

	}

	if (yvar == "nN"){
		tree$edge.length = nN.branch
	}


	tree


}



# Read in data
metadata="data/mammal.traits.tsv"
tree.in = "data/beastmap.sample.trees"


# Read in the sample of beastmap trees, takes ~10 seconds
trees.all = phylowise::readBeastTrees(tree.in, burnin=0.1)


# Traits
traits.df = read.table(metadata, sep="\t", header=T)
traits.df$species = traits.df$Species
traits.df$MaxLongevity = traits.df$MaxLongevity * 30.4 # From months to days
traits.df$OffspringPerYear = traits.df$LitterSize * traits.df$LittersPerYear 



# Our window sizes, in units of expected number of subst per site
DIST.MIN=0.005
DIST.MAX=0.1



# Our life history traits
traits.all = c("AdultBodyMassGrams", "SexualMaturityAgeDays", "MaxLongevity", "OffspringPerYear")


ntrees = length(trees.all)






# Put the traits into the same ordering as our tree
time.tree = trees.all[[1]]
traits.tree = sapply(time.tree$tip.label, function(ele) {
	w = which(traits.df$species == ele)
	if (length(w) == 0){
		print(ele)
		rep(NA, length(traits.all))
	}else{
		as.numeric(traits.df[w[1],traits.all])
	}
	
})
if (length(traits.all) > 1){
	traits.tree = t(traits.tree)
}
traits.tree=as.data.frame(traits.tree)
colnames(traits.tree) = traits.all


# Log and standardise the traits
traits.tree.log = traits.tree
for (var in traits.all){
	val = as.numeric(traits.tree[,var])
	traits.tree.log[,var] = log(val)

	# Scale and log
	val = scale(log(val))
	val[is.nan(val)] = NA
	val = as.numeric(val)

	traits.tree[,var] = val
}

ngood = sum(apply(traits.tree, 1, function(row) all(!is.na(row))))
print(ngood)





# We will peform 10 iterations in 10 trees (sampled with replacement) 
nsamples=10
results.df = data.frame(var1=character(0), var2=character(0), p=numeric(0), rho=numeric(0), ppc.prob=numeric(0), ppc.bf=numeric(0))


features.importance.ns = list()
features.importance.nn = list()
for (var in traits.all){
	features.importance.ns[[var]] = numeric(0)
	features.importance.nn[[var]] = numeric(0)
}


features.pvalue.ns = list()
features.pvalue.nn = list()


bma.all.ns = list()
bma.all.nn = list()

for (tree.nr in 1:nsamples){


	cat(paste(tree.nr, "/", nsamples, "\n"))


	# Sample a random tree
	tree.rand = sample(1:ntrees, 1)
	time.tree = trees.all[[tree.rand]]


	# Get the correlation between nN and nS
	subst.tree1 = get.tree.with.lengths(time.tree, yvar="nS")
	subst.tree2 = get.tree.with.lengths(time.tree, yvar="nN")
	n = length(time.tree$tip.label)
	nS.leaves = scale(log(ape::node.depth.edgelength(subst.tree1)[1:n]))
	nN.leaves = scale(log(ape::node.depth.edgelength(subst.tree2)[1:n]))
	pairs.df = phylowise::sampleTaxonPairs(subst.tree=subst.tree1, response=nS.leaves, covariate=nN.leaves, window.tree=time.tree, dist.min=DIST.MIN, dist.max=DIST.MAX, maximise=FALSE, verbose=F) # Will not trait-max to avoid asymmetries between nN and nS being the response


	if (nrow(pairs.df) > 10){
		ppc = phylowise::PPC.test(pairs.df, logY=F, extreme.value.threshold=3) # Do not log the response variable because it is already logged (i.e., set logY=F)
		trait.cor.df2 = data.frame(var1="nS", var2="nN", p=as.numeric(ppc$p.var), rho=as.numeric(ppc$rho.var), ppc.prob=as.numeric(ppc$bma.probs["trait"]), ppc.bf=as.numeric(ppc$bma.bf["trait"]))
		trait.cor.df3 = data.frame(var1="nN", var2="nS", p=as.numeric(ppc$p.var), rho=as.numeric(ppc$rho.var), ppc.prob=as.numeric(ppc$bma.probs["trait"]), ppc.bf=as.numeric(ppc$bma.bf["trait"]))
		results.df = rbind(results.df, trait.cor.df2, trait.cor.df3)
	}
	

	for (yvar in c("nS", "nN")){


		# Get a tree with branch lengths as nN or nS
		subst.tree = get.tree.with.lengths(time.tree, yvar=yvar)


		for (var in traits.all){
			pairs.df = phylowise::sampleTaxonPairs(subst.tree=subst.tree, window.tree=time.tree, covariate=traits.tree[,var], dist.min=DIST.MIN, dist.max=DIST.MAX, verbose=F)
			
			if (nrow(pairs.df) < 10){
				next
			}else{
				ppc = phylowise::PPC.test(pairs.df, extreme.value.threshold=3)
				print(paste0("n=", nrow(pairs.df), " p=", signif(ppc$p.var)))

			

				trait.cor.df2 = data.frame(var1=yvar, var2=var, p=as.numeric(ppc$p.var), rho=as.numeric(ppc$rho.var), ppc.prob=as.numeric(ppc$bma.probs["trait"]), ppc.bf=as.numeric(ppc$bma.bf["trait"]))
				trait.cor.df3 = data.frame(var2=yvar, var1=var, p=as.numeric(ppc$p.var), rho=as.numeric(ppc$rho.var), ppc.prob=as.numeric(ppc$bma.probs["trait"]), ppc.bf=as.numeric(ppc$bma.bf["trait"]))
				results.df = rbind(results.df, trait.cor.df2, trait.cor.df3)

			}
			

		}

	




		# Joint analyses
		pairs.df = phylowise::sampleTaxonPairs(subst.tree=subst.tree, window.tree=time.tree, covariate=traits.tree.log[,traits.all], dist.min=DIST.MIN, dist.max=DIST.MAX, verbose=F)

		if (nrow(pairs.df) < 10){
			next
		}

		ppc = phylowise::PPC.test(pairs.df)
		print(ppc$n)

		if (ppc$n < length(traits.all)+2){
			next
		}

	



		# Get best GLM 
		dat = ppc$data.df




		# BMA
		prior.weight = 0.2
		bma = bic.glm(x = dat[,traits.all], y = dat$distance.response, glm.family=gaussian(), prior.param=prior.weight, OR=1000)
		if (yvar == "nS"){
			bma.all.ns[[tree.nr]]=bma
		}else{
			bma.all.nn[[tree.nr]]=bma
		}
		print(paste(bma$label[1], "p=", signif(bma$postprob[1], 2)))

		




		# Random forest
		forest = randomForest(distance.response ~ ., data = dat)
		imp = importance(forest)
		best = rownames(imp)[which(imp == max(imp))]
		print(best)


		for (var in traits.all){
			if (yvar == "nS"){
				features.importance.ns[[var]] = c(features.importance.ns[[var]], imp[var,])
			}else{
				features.importance.nn[[var]] = c(features.importance.nn[[var]], imp[var,])
			}
		}


		# PGLS between traits
		for (var1 in traits.all){
			for (var2 in traits.all){

				if (var1 == var2){
					next
				}

				x = traits.tree[,var1]
				y = traits.tree[,var2]


				pairs.df = phylowise::sampleTaxonPairs(subst.tree=subst.tree, response=y, window.tree=time.tree, 
						covariate=x, dist.min=DIST.MIN, dist.max=DIST.MAX, maximise=FALSE, verbose=F)

				if (nrow(pairs.df) < 10){
					next
				}

				ppc = phylowise::PPC.test(pairs.df, logY=F)
				p=as.numeric(ppc$p.var)
				rho = ppc$rho.var
				#print(ppc$n)

				trait.cor.df2 = data.frame(var1=var1, var2=var2, p=p, rho=rho, ppc.prob=as.numeric(ppc$bma.probs["trait"]), ppc.bf=as.numeric(ppc$bma.bf["trait"]))
				results.df = rbind(results.df, trait.cor.df2)



			}
		}

	}

}




draw.points = function(x, y, base.col, nlayers=5){

	par(xpd=F)
	cex.all = seq(from=1.6, to=1.1, length=nlayers)
	points(x, y, pch=21, col=NA, bg="white", cex=cex.all[1])
	for (j in 1:nlayers){
		points(x, y, pch=21, col=NA, bg=paste0(base.col, "99"), cex=cex.all[j])
	}


	cex.all = seq(from=1.1, to=0.75, length=nlayers)
	points(x, y, pch=21, col=NA, bg="white", cex=cex.all[1])
	for (j in 1:nlayers){
		points(x, y, pch=21, col=NA, bg=paste0(base.col, "55"), cex=cex.all[j])
	}


	cex.all = seq(from=0.75, to=0.5, length=nlayers)
	points(x, y, pch=21, col=NA, bg="white", cex=cex.all[1])
	for (j in 1:nlayers){
		points(x, y, pch=21, col=NA, bg=paste0(base.col, "11"), cex=cex.all[j])
	}


	#points(x, y, pch=21, col=NA, bg=paste0(base.col, "33"), cex=0.5)
	par(xpd=T)
}




add.label = function(letter, dy=0){


	usr = par("usr")
	letter = toupper(letter)

	# temporarily set coordinates to 0-1
	par(usr = c(0, 1, 0, 1))

	text(-0.27, 1.17 + dy, letter, adj=0, cex=1.3)

	# restore original coordinates
	par(usr = usr)

}



plotPairs2 = function(tree, pairs.df, edge.col="red", show.tip.label=F, edge.width=1, edge.width.pairs=3, label.cex=1){


	# Validation
	for (i in 1:nrow(pairs.df)){
		tip1 = pairs.df[i,"tip1"]
		tip2 = pairs.df[i,"tip2"]

		index1 = which(tree$tip.label == tip1)
		index2 = which(tree$tip.label == tip2)

		if (length(index1) == 0){
			stop(paste("Cannot find taxon", tip1))
		}
		if (length(index2) == 0){
			stop(paste("Cannot find taxon", tip2))
		}

	}


	tree.plot=tree
	if (show.tip.label){
		tree.plot$tip.label = as.character(sapply(tree.plot$tip.label, function(ele) ifelse(any(ele == pairs.df$tip1) | any(ele == pairs.df$tip2), ele, "")) )
	}

	plot(tree.plot, show.tip.label=show.tip.label, direction="down", edge.width=edge.width, cex=label.cex)
	axis(2, las=2)
	tree.plot$tip.label = tree$tip.label


	# Grab node coordinates created by plot.phylo
	lp = get("last_plot.phylo", envir = .PlotPhyloEnv)

	for (i in 1:nrow(pairs.df)){
		tip1 = pairs.df[i,"tip1"]
		tip2 = pairs.df[i,"tip2"]

		index1 = which(tree.plot$tip.label == tip1)
		index2 = which(tree.plot$tip.label == tip2)

		if (length(index1) == 0){
			stop(paste("Cannot find taxon", tip1))
		}
		if (length(index2) == 0){
			stop(paste("Cannot find taxon", tip2))
		}

		# Get the path between these two, excluding their mrca
		path = nodepath(tree, index1, index2)
		mrca = getMRCA(tree, c(index1, index2))
		path = path[path != mrca]
		for (child in path){

			parent = tree$edge[tree$edge[,2] == child,1]

			# Vertical line
			lines(c(lp$xx[child], lp$xx[child]), c(lp$yy[child],lp$yy[parent]), lwd=edge.width.pairs, col=edge.col)

			# Horizontal shoulder
			lines(c(lp$xx[child], lp$xx[parent]), c(lp$yy[parent], lp$yy[parent]), lwd=edge.width.pairs, col=edge.col)
 

		}


	}

}


tidy.name = function(species){

	bits = strsplit(species, "_")[[1]]
	if (bits[2] == "canadensis"){
		paste0(bits[1], " ", bits[2])
	}else{
		paste0(substr(bits[1], 1, 1), ".", bits[2])
	}

}

mammal.cols = c(AdultBodyMassGrams="#DCAE1E", MaxLongevity="#70482A", SexualMaturityAgeDays="#F09A09", OffspringPerYear="#61534E")
var.names = c(AdultBodyMassGrams="Adult body mass", MaxLongevity="Maximum longevity",  SexualMaturityAgeDays="Sexual maturity age", OffspringPerYear="Offspring per year")
labels.short = c("Mass", "Long.", "Matur.", "Offspr.")
YVAR.NICE="nS"
yvar="nS"


# Plot just one 
var="AdultBodyMassGrams"
p.val=1
time.tree = trees.all[[1]]
subst.tree = get.tree.with.lengths(time.tree, yvar="nS")


time.tree$tip.label = sapply(time.tree$tip.label, tidy.name)
subst.tree$tip.label = sapply(subst.tree$tip.label, tidy.name)


# A random sample, will not always be significant but it usually will be
pairs.df = phylowise::sampleTaxonPairs(subst.tree=subst.tree, window.tree=time.tree, covariate=traits.tree[,var], dist.min=DIST.MIN, dist.max=DIST.MAX, verbose=F)
if (nrow(pairs.df) < 10){
	next
}
ppc = phylowise::PPC.test(pairs.df, extreme.value.threshold=3)
p.val = as.numeric(ppc$p.var)
slope = lm(ppc$data.df$distance.response ~ ppc$data.df$trait + 0)$coefficients




png("mammals.png", width=4200, height=3200, res=500)


layout(matrix(c(1,3,3,3,
				1,2,2,2,
				4,5,6,7, 
				8,9,10,11, 
				8,9,10,11),nrow=5, byrow=T), width=c(1,1,1,1), height=c(0.3, 1.7, 2, 1.6, 0.4))
#par(mfrow=c(nrow, ncol))
par(mar=c(3.5, 3.5, 2.8, 1.1))
par(xpd=T)








par(xpd=F)
xlim=c(-20, 20)
ylim=c(-7, 7)
plot(0,0, type="n", xlim=xlim, ylim=ylim, main="Single trial", axes=F, xaxs="i", yaxs="i", xlab="", ylab="")
draw.points(ppc$data.df$trait, ppc$data.df$distance.response, mammal.cols[var])
mtext(paste0("p=", signif(p.val, 2), ", Pearson=", signif(ppc$rho.var, 2)), cex=0.65)

lines(xlim, c(0,0), lty="dashed")
lines(c(0, 0), ylim, lty="dashed")

lines(xlim, xlim*slope, col="#ff000055", lwd=5)
lines(xlim, xlim*slope, col="#ff000099", lwd=3.5)
lines(xlim, xlim*slope, col="#ff0000", lwd=1.5)





par(xpd=T)
axis(1)
axis(2, las=2)
text(0, ylim[1] - (ylim[2]-ylim[1])*0.25, "Δ Adult body mass", cex=1.1)
text(xlim[1] - (xlim[2]-xlim[1])*0.2, 0, paste0("Δ ", YVAR.NICE), cex=1.1, srt=90)
add.label("A")


# Plot a tree
par(mar=c(0.5,3.5,1,0.5))


subst.tree.plot = subst.tree
subst.tree.plot$edge.length = subst.tree.plot$edge.length / 1000
plotPairs2(subst.tree.plot, pairs.df, edge.col=mammal.cols["AdultBodyMassGrams"], edge.width.pairs=1.2, edge.width=0.25, show.tip.label=T, label.cex=0.35)



# Legend above tree
par(xpd=T)
par(mar=c(0,0,0,0))
plot(0,0,type="n",axes=F,main="", xlab="",ylab="", xlim=c(0,1), ylim=c(0,1))
text(0.5, 0.4, "PPC pair sample", font=2, cex=1.2)
text(-0.005, 0.4, "B", adj=0, cex=1.3)



par(mar=c(3.5, 3.5, 2.8, 1.1))


xlab=T
ylab=T
letterz=c("C", "D", "E", "F")
i=1
for (trait in names(var.names)){

	sub.df = results.df[results.df$var1 == trait & results.df$var2 == yvar,]
	n = sum(!is.na(traits.tree.log[,trait]))
	
	plot(sub.df$rho, log(sub.df$p, 10), xlim=c(-1,1), ylim=c(-8,0.5), axes=F, yaxs="i", xaxs="i", main=paste0(var.names[trait], " vs. nS"), xlab="", ylab="", type="n")

	par(xpd=T)
	if (xlab){
		text(0, -9.9, paste0("Pearson correlation with ", YVAR.NICE), cex=1.1)
	}
	
	if (ylab){
		text(-1.4, -3.25, "Log₁₀ p-value", cex=1.1, srt=90)
		ylab=F
	}
	

	col1=paste0(as.character(mammal.cols[trait]), "77")
	col2=paste0(mammal.cols[trait], "22")
	x = sub.df$rho
	y = log(sub.df$p, 10)


	lines(c(-1,1), c(-2,-2), lty="dashed", col="#d3d3d3")
	lines(c(-1,1), c(-4,-4), lty="dashed", col="#d3d3d3")
	lines(c(-1,1), c(-6,-6), lty="dashed", col="#d3d3d3")



	text(0.05, -7.6, "+ve", adj=0, cex=1.2)
	text(-0.05, -7.6, "−ve", adj=1, cex=1.2)

	par(xpd=F)
	lines(c(0,0), c(-8,0), lty="dashed", col="#d3d3d3")


	draw.points(x, y, mammal.cols[trait])

	axis(1)
	axis(2, las=2)


	ymin=1e-6
	bma.p = sub.df$ppc.prob
	bma.p[bma.p>1-ymin] = 1-ymin
	bma.p[bma.p<ymin] = ymin
	bma.p.mean = plogis(mean(qlogis(bma.p)))
	prior.prob=0.5

	z = (bma.p.mean/(1-bma.p.mean)) / (prior.prob / (1-prior.prob)) 
	z = log(z, 10)


	mtext(paste0("log₁₀ K=", signif(z, 3), ", n=", n, " taxa"), cex=0.6)#, ", p2=", signif(ps1, 2), ", p3=", signif(ps2, 2)), cex=0.7)
	par(xpd=T)

	add.label(letterz[i])
	i = i+1


}

par(mar=c(3.5, 3.5, 2.8, 1.1))
par(xpd=T)




# # Correlations
my_colors = paste0(rev(c("#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30")), "88")


correlations = matrix(0, nrow=6, ncol=6)
rownames(correlations) = c("nN", "nS", names(var.names))
colnames(correlations) = rownames(correlations)
pv=correlations
for (var1 in colnames(correlations)){
	for (var2 in colnames(correlations)){

		sub.df = results.df[(results.df$var1 == var1 & results.df$var2 == var2) |
				 (results.df$var1 == var2 & results.df$var2 == var1),]
		rho=1
		p=1
		if (var1 != var2){
			rho = median(sub.df$rho)
			#p = p.hmp(sub.df$p, L=nrow(sub.df))

			ymin=1e-6
			bma.p = sub.df$ppc.prob
			bma.p[bma.p>1-ymin] = 1-ymin
			bma.p[bma.p<ymin] = ymin
			bma.p.mean = plogis(mean(qlogis(bma.p)))
			prior.prob=0.5

			z = (bma.p.mean/(1-bma.p.mean)) / (prior.prob / (1-prior.prob)) 
			z = log(z, 10)
			print(paste(var1, var2, z))

			if (z < 1){
				rho = NA
			}
		}

		correlations[var1,var2] = rho
		pv[var1,var2] = p
	}
}
rownames(correlations) = c("nN", "nS", labels.short)
colnames(correlations) = rownames(correlations)
colnames(pv) = rownames(correlations)
rownames(pv) = rownames(correlations)

corrplot(correlations, col=my_colors, method="ellipse", type="full", diag=TRUE, tl.col="black", cl.ratio = 0.3, addCoef.col = "black", number.cex=0.7, na.label="–")
add.label("G", dy=0.01)




par(mar=c(2.5, 3.5, 2.8, 1.1))
rw=0.4
xlim = c(1-2*rw, length(var.names)+2*rw)
ymax=0.8
plot(0,0, type="n", xlim=xlim, ylim=c(0, ymax), main="Most important predictor", axes=F, xaxs="i", yaxs="i", xlab="", ylab="")


# What was the best in each random forest
nbest.ns = rep(0, length(var.names))
nbest.nn = rep(0, length(var.names))
names(nbest.ns) = names(var.names)
names(nbest.nn) = names(nbest.ns)
for (j in 1:length(features.importance.ns[[1]])){

	imps = numeric(0)
	for (trait in names(var.names)){
		v = features.importance.ns[[trait]][j]
		imps = c(imps, v)
	}
	best = which(imps == max(imps))
	if (length(best) > 1){
		print("tie")
	}
	nbest.ns[best] = nbest.ns[best] + 1/length(best) 

}
nbest.ns = nbest.ns / sum(nbest.ns)
nbest.rename.ns = nbest.ns
names(nbest.rename.ns) = labels.short

for (j in 1:length(features.importance.nn[[1]])){

	imps = numeric(0)
	for (trait in names(var.names)){
		v = features.importance.nn[[trait]][j]
		imps = c(imps, v)
	}
	best = which(imps == max(imps))
	if (length(best) > 1){
		print("tie")
	}
	nbest.nn[best] = nbest.nn[best] + 1/length(best) 

}
nbest.nn = nbest.nn / sum(nbest.nn)
nbest.rename.nn = nbest.nn
names(nbest.rename.nn) = labels.short


i=1
for (trait in names(nbest.ns)){
	#imp = mean(features.importance[[trait]])
	imp2 = nbest.ns[trait]
	imp1 = nbest.nn[trait]
	#bas = bas.all.ns[trait]

	col1 = paste0(as.character(mammal.cols[trait]), "aa")
	col2 = paste0(as.character(mammal.cols[trait]), "66")



	rect(i-rw, 0, i, imp1, col=col1, border=paste0(as.character(mammal.cols[trait]), "cc"), lwd=1)
	rect(i, 0, i+rw, imp2, col=col1, border=paste0(as.character(mammal.cols[trait]), "cc"), lwd=1)
	text(i-rw/2, imp1 + 0.035, "nN", cex=0.8)
	text(i+rw/2, imp2 + 0.035, "nS", cex=0.8)


	i = i + 1
}

axis(1, at=xlim, labels=c("", ""))
axis(1, at=1:4, labels=names(nbest.rename.ns), cex.axis=0.8)
axis(2, las=2)
text(-0.85, ymax/2, "Proportion of random forests", cex=1.1, srt=90)


add.label("H")




# nS
bas.predictors.ns = list()
for (var in names(var.names)){
	bas.predictors.ns[[var]] = numeric(0)
}
bas.model.probs.ns = numeric(0)
names(bas.predictors.ns) = names(var.names)
for (i in 1:length(bma.all.ns)){
	bma = bma.all.ns[[i]]


	prob.total = numeric(0)
	for (var in names(var.names)){
		prob.total[var] = 0
	}

	for (modelNr in 1:length(bma$label)){
		prob = bma$postprob[modelNr] 
		coefficients = strsplit(bma$label[modelNr], ",")[[1]]
		model.str = paste0(coefficients, collapse="+")	
		
		if (is.na(bas.model.probs.ns[model.str])){
			bas.model.probs.ns[model.str] = prob / length(bma.all.ns)
		}else{
			bas.model.probs.ns[model.str] = bas.model.probs.ns[model.str] + prob / length(bma.all.ns)
		}

		for (coeff in coefficients){
			if (coeff == "NULL"){
				next
			}
			prob.total[coeff] = prob.total[coeff] + prob
		}

	}

	for (var in names(var.names)){
		bas.predictors.ns[[var]] = c(bas.predictors.ns[[var]], as.numeric(prob.total[var]))
	}
}




# nN
bas.predictors.nn = list()
for (var in names(var.names)){
	bas.predictors.nn[[var]] = numeric(0)
}
bas.model.probs.nn = numeric(0)
names(bas.predictors.nn) = names(var.names)
for (i in 1:length(bma.all.nn)){
	bma = bma.all.nn[[i]]


	prob.total = numeric(0)
	for (var in names(var.names)){
		prob.total[var] = 0
	}

	for (modelNr in 1:length(bma$label)){
		prob = bma$postprob[modelNr] 
		coefficients = strsplit(bma$label[modelNr], ",")[[1]]
		model.str = paste0(coefficients, collapse="+")	
		
		if (is.na(bas.model.probs.nn[model.str])){
			bas.model.probs.nn[model.str] = prob / length(bma.all.nn)
		}else{
			bas.model.probs.nn[model.str] = bas.model.probs.nn[model.str] + prob / length(bma.all.nn)
		}

		for (coeff in coefficients){
			if (coeff == "NULL"){
				next
			}
			prob.total[coeff] = prob.total[coeff] + prob
		}

	}

	for (var in names(var.names)){
		bas.predictors.nn[[var]] = c(bas.predictors.nn[[var]], as.numeric(prob.total[var]))
	}
}



ymax=0.8
plot(0,0, type="n", xlim=xlim, ylim=c(0, ymax), main="Bayesian model averaging", axes=F, xaxs="i", yaxs="i", xlab="", ylab="")

ymin=1e-6
prior.prob = 0.2

lines(xlim, c(prior.prob,prior.prob), lty="dashed", lwd=0.5, col="#696969")
text(xlim[2]-0.02, 0.23, "Prior", adj=1,  col="#696969", cex=0.8)

i=1
for (trait in names(bas.predictors.nn)){

	if (trait == "NULL"){
		next
	}
	


	bas1 = mean(bas.predictors.nn[[trait]])
	bas2 = mean(bas.predictors.ns[[trait]])
	


	col1 = paste0(as.character(mammal.cols[trait]), "aa")
	col2 = paste0(as.character(mammal.cols[trait]), "66")


	rect(i-rw/2, bas1, i, bas1 + 0.08, col="#ffffff77", border=NA)
	rect(i, bas2, i+rw/2, bas2 + 0.08, col="#ffffff77", border=NA)

	rect(i-rw, 0, i, bas1, col="white", border=NA, lwd=1)
	rect(i, 0, i+rw, bas2, col="white", border=NA, lwd=1)

	rect(i-rw, 0, i, bas1, col=col1, border=paste0(as.character(mammal.cols[trait]), "cc"), lwd=1)
	rect(i, 0, i+rw, bas2, col=col1, border=paste0(as.character(mammal.cols[trait]), "cc"), lwd=1)


	text(i-rw/2, bas1 + 0.035, "nN", cex=0.8)
	text(i+rw/2, bas2 + 0.035, "nS", cex=0.8)


	i = i + 1
}



axis(1, at=xlim, labels=c("", ""))
axis(1, at=1:4, labels=names(nbest.rename.ns), cex.axis=0.8)
axis(2, las=2)
text(-0.85, ymax/2, "Posterior probability", cex=1.1, srt=90)


add.label("I")





get.best.model.str = function(bma.all, bas.model.probs, yvar, rankNr=1){


	map.model = names(sort(bas.model.probs, decreasing=T)[rankNr])
	prob = as.numeric(sort(bas.model.probs, decreasing=T)[rankNr])
	traits = strsplit(map.model, "[+]")[[1]] # none is not handled here
	if (traits[1] == "NULL"){
		traits = character(0)
	}
	coefficients.all = list()
	for (v in c("Intercept", traits)){
		coefficients.all[[v]] = numeric(0)
	}

	for (i in 1:length(bma.all)){
		bma = bma.all[[i]]

		nmodels = length(bma$label)
		
		modelNr = which(bma$label == map.model)

		intercept = bma$mle[modelNr,1]
		slopes = bma$mle[modelNr,traits]

		coefficients.all[["Intercept"]] = c(coefficients.all[["Intercept"]], intercept)
		for (v in traits){
			coefficients.all[[v]] = c(coefficients.all[[v]], as.numeric(slopes[v]))
		}

		
	}


	coeff.median = numeric(0)
	for (var in names(coefficients.all)){
		coeff.median[var] = mean(coefficients.all[[var]], na.rm=T)
	}
	varnames.coeff = c(Intercept="", AdultBodyMassGrams="ΔMass", MaxLongevity="ΔLong", SexualMaturityAgeDays="ΔMatur", OffspringPerYear="ΔOffspr")



	varnames.coeff = as.character(varnames.coeff[names(coeff.median)])

	str = paste(signif(as.numeric(coeff.median),2), varnames.coeff)
	str = paste0(str, collapse=" + ")
	str = paste0("Δ", yvar, " = ", str)
	str = gsub("[+] -", "− ", str)
	str = gsub("-", "− ", str)


	str = paste0(rankNr, ". ", str)
	p1 = paste0("(p=", signif(prob, 2), ")")

	c(str, p1)


}

str1a = get.best.model.str(bma.all.ns, bas.model.probs.ns, yvar="nS", rankNr=1)
str1b = get.best.model.str(bma.all.ns, bas.model.probs.ns, yvar="nS", rankNr=2)
str1c = get.best.model.str(bma.all.ns, bas.model.probs.ns, yvar="nS", rankNr=3)

str2a = get.best.model.str(bma.all.nn, bas.model.probs.nn, yvar="nN", rankNr=1)
str2b = get.best.model.str(bma.all.nn, bas.model.probs.nn, yvar="nN", rankNr=2)
str2c = get.best.model.str(bma.all.nn, bas.model.probs.nn, yvar="nN", rankNr=3)



size1 = 0.8
size2 = 0.7

plot(0,0,type="n", axes=F, xlab="", ylab="", xaxs="i", yaxs="i", main="Best multivariate models", xlim=c(0,1), ylim=c(0,1))

text(0.1, 0.95, str1a[1], cex=size1, adj=0)
text(0.15, 0.88, str1a[2], cex=size2, adj=0, col="#696969")

text(0.1, 0.8, str1b[1], cex=size1, adj=0)
text(0.15, 0.73, str1b[2], cex=size2, adj=0, col="#696969")

text(0.1, 0.65, str1c[1], cex=size1, adj=0)
text(0.15, 0.58, str1c[2], cex=size2, adj=0, col="#696969")





text(0.1, 0.4, str2a[1], cex=size1, adj=0)
text(0.15, 0.33, str2a[2], cex=size2, adj=0, col="#696969")

text(0.1, 0.25, str2b[1], cex=size1, adj=0)
text(0.15, 0.18, str2b[2], cex=size2, adj=0, col="#696969")

text(0.1, 0.1, str2c[1], cex=size1, adj=0)
text(0.15, 0.03, str2c[2], cex=size2, adj=0, col="#696969")

add.label("J")


dev.off()


