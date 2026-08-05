
library(ape)
library(phylowise)


# Parts of this code were modified from 
	# https://datadryad.org/dataset/doi:10.5061/dryad.43mg3 >>> PGLS folder
	# Lanfear, R., Ho, S., Jonathan Davies, T. et al. Taller plants have lower rates of molecular evolution. Nat Commun 4, 1879 (2013). https://doi.org/10.1038/ncomms2836


# load absolute rates and height
# load the trees, Each tree has a column in the datafile, and the first column is the ML estimate
trees = read.tree("data/R8S_trees.txt")
height.data = read.delim("data/raw_data_sister_pairs.csv", sep=",")
h = as.numeric(c(height.data$Height.A, height.data$Height.B))

temp = as.numeric(c(height.data$TEMP.A, height.data$TEMP.B))
uv = as.numeric(c(height.data$UV.A, height.data$UV.B))

Family = c(as.character(height.data$taxA), as.character(height.data$taxB))
traits.df = as.data.frame(cbind(Family, height=h, temp=temp, uv=uv))

traits.df$height = as.numeric(traits.df$height)
traits.df$temp = as.numeric(traits.df$temp)
traits.df$uv = as.numeric(traits.df$uv)



# Load rates
rates = read.csv("data/bootstrap_rates.txt", sep=",", header=TRUE)
rates[-1] = log(rates[-1])



# Merge traits and rates together
data.df = merge(traits.df, rates)



dist.min.all=c(0.005, 0.01, 0.02, 0.05, 0.10)
dist.span.all=c(2, 3, 5, 10, Inf)




trait = "height" # Could also look at temp or uv, which are non-significant



# Record any tips and rows for which we don't have height data
d = as.character(data.df$Family[which(is.na(data.df[,trait]))])
d = c("Ginkgo_biloba", d)
data.df = data.df[!is.na(data.df[,trait]), ]



# There are 1000 bootstrap replicates in total
# We will just run it on 10 iterations to save time
nreps=10
ppc.df = data.frame(treeNr = numeric(0), trait=character(0), nested=numeric(0), maximise=logical(0), ppc.p = numeric(0), ppc.rho = numeric(0), ppc.prob=numeric(0), ppc.bf=numeric(0), npairs=numeric(0))
for(i in 1:nreps){

	cat(paste(i, "\n"))


	# Prepare the data frame for this replicate
	t = trees[[i]]
	df_new = cbind(data.df[i+5], data.df$height, data.df$temp, data.df$uv, data.df$Family) #columns 1 and 2 are Family name and height
	names(df_new)=c('rate', 'height', 'temp', 'uv','Family')
	rownames(df_new) = as.character(df_new$Family)
	

	# Get a vector of trait values
	traits = numeric(0)
	for (taxon in t$tip.label){

		if (any(df_new$Family == taxon)){
			val = df_new[df_new$Family == taxon,trait]
		}else{
			val=NA
		}
		traits = c(traits, val)
	}



	# Build a distance matrix from the subst tree so that we don't need to recompute this on each iteration
	dmat = getDistanceMatrix(t)


	# Run PPC on a range of different window sizes, using both nested and non-nested
	for (dist.min in dist.min.all){
		cat(paste("\twindow-lower:", dist.min, "\n"))
		for (dist.span in dist.span.all){
			dist.max = dist.min * dist.span
			cat(paste("\t\twindow-upper:", dist.max, ", p= "))
			
			for (nested in c(FALSE, TRUE)){
				for (maximise in c(FALSE, TRUE)){


					# Keep sampling until we have at least 3 pairs to compare
					nobs = 0
					nattempts = 1
					while (nobs < 3 & nattempts < 5){

						# Sample pairs
						pairs.df = phylowise::sampleTaxonPairs(t, traits, dist.min=dist.min, dist.max=dist.max, distance.matrix=dmat, nested=nested, maximise=maximise, verbose=F)
						nobs = nrow(pairs.df)
						nattempts = nattempts + 1
					}

					
					if (nobs < 3){
						ppc.df2 = data.frame(treeNr = i, trait=trait, nested=nested, maximise=maximise, dist.min=dist.min, dist.span=dist.span, ppc.p = 1, ppc.rho = 0, ppc.prob=0, ppc.bf=NA,  npairs=nobs)
						ppc.df = rbind(ppc.df, ppc.df2)
						next
					}
					

					result = phylowise::PPC.test(pairs.df, standardise=1, logY=T, extreme.value.threshold=0)
					cat(paste0(signif(result$p, 3), ", "))
					#print(result$p)


					ppc.df2 = data.frame(treeNr = i, trait=trait, nested=nested, maximise=maximise, dist.min=dist.min, dist.span=dist.span, ppc.p = as.numeric(result$p), ppc.rho = as.numeric(result$rho.var),
							ppc.prob=as.numeric(result$bma.probs["trait"]), ppc.bf=as.numeric(result$bma.bf["trait"]), npairs=result$n)
					ppc.df = rbind(ppc.df, ppc.df2)

				}
			}
			cat("\n")
		}
	}



}





# Make a plot of p-value vs Pearson correlation across many replicates of PPC
plot.ppc = function(sub.df){

	plot(sub.df$ppc.rho, log(sub.df$ppc.p, 10), xlim=c(-1,1), ylim=c(-8,0.5), axes=F, yaxs="i", xaxs="i", main="PPC", xlab="Pearson correlation", ylab="Log₁₀ p-value", type="n")


	lines(c(-1,1), c(-2,-2), lty="dashed", col="#d3d3d3")
	lines(c(-1,1), c(-4,-4), lty="dashed", col="#d3d3d3")
	lines(c(-1,1), c(-6,-6), lty="dashed", col="#d3d3d3")
	lines(c(-1,1), c(-8,-8), lty="dashed", col="#d3d3d3")


	text(0.05, -7.6, "+ve", adj=0)
	text(-0.05, -7.6, "–ve", adj=1)

	lines(c(0,0), c(-8,0), lty="dashed", col="#d3d3d3")

	axis(1)
	axis(2, las=2)



	x = sub.df$ppc.rho
	y = log(sub.df$ppc.p, 10)
	draw.points(x, y)

	
	ymin=1e-6
	bma.p = sub.df$ppc.prob
	bma.p[bma.p>1-ymin] = 1-ymin
	bma.p[bma.p<ymin] = ymin
	bma.p.mean = plogis(mean(qlogis(bma.p)))
	prior.prob=0.5

	z = (bma.p.mean/(1-bma.p.mean)) / (prior.prob / (1-prior.prob)) 
	z = log(z, 10)
	z = signif(z, 3)


	mtext(paste0("log₁₀ K = ", z), cex=0.9)

	

}


# Draw points onto the plots
draw.points = function(x, y, nlayers=5){

	# Jungle colour scheme
	colours = c("#a6e3aa", "#6fb574", "#39823e", "#1c5b20", "#072e0a")

	par(xpd=F)
	cex.all = seq(from=1.5, to=1.1, length=nlayers)
	points(x, y, pch=21, col=NA, bg="white", cex=cex.all[1])
	for (j in 1:nlayers){
		points(x, y, pch=21, col=NA, bg=paste0(colours[4], "aa"), cex=cex.all[j])
	}


	cex.all = seq(from=1.1, to=0.75, length=nlayers)
	points(x, y, pch=21, col=NA, bg="white", cex=cex.all[1])
	for (j in 1:nlayers){
		points(x, y, pch=21, col=NA, bg=paste0(colours[3], "aa"), cex=cex.all[j])
	}


	cex.all = seq(from=0.75, to=0.5, length=nlayers)
	points(x, y, pch=21, col=NA, bg="white", cex=cex.all[1])
	for (j in 1:nlayers){
		points(x, y, pch=21, col=NA, bg=paste0(colours[2], "aa"), cex=cex.all[j])
	}


	points(x, y, pch=21, col=NA, bg=paste0(colours[1], "aa"), cex=0.5)
	par(xpd=T)
}






# Draw a ppc plot for a single setting. Each point in this plot is a single replicate of PPC, which will give different p-values and different Pearson correlations
# When there are enough samples, the correlation is usually negative and significant


# Nested and maximised, with a narrow window of (0.01, 0.05)
sub.df1 = ppc.df[ppc.df$nested & ppc.df$maximise & ppc.df$dist.min == 0.01 & ppc.df$dist.span == 5,]
plot.ppc(sub.df1)



# Nested and not-maximised, with a wide window of (0.005, Inf)
sub.df2 = ppc.df[ppc.df$nested & ppc.df$maximise & ppc.df$dist.min == 0.005 & ppc.df$dist.span == Inf,]
plot.ppc(sub.df2)




