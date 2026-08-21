


#' Phylogenetic pairwise contrast test on a series of taxon pairs, using ordinary 
#' least squares and Bayesian model averaging (BMA). This method does not do multivariate
#' regression, but it returns the data frame that can be used for it.
#'
#' @param pairs.df data frame of taxon pairs
#' @param logY should we take the logarithm of the evolutionary distance? default and recommended setting: true
#' @param standardise should we standardise the evolutionary distance? 0 for no, 1 to standardise by sqrt(distance), 
#'						2 to standardise by sqrt(mrca age), 3 to standarise by mrca age. default and recommended setting: 2
#' @param extreme.value.threshold remove any observations more than this many standard deviations away from 
#'									the mean response variable; set to zero for no extreme value removal
#' @param prior.weight trait inclusion prior probbaility for Bayesian model averaging (see BMA::bic.glm)
#' @param OR Occam's window for Bayesian model averaging (see BMA::bic.glm)
#' @param epsilon precision of BMA probabilities for Bayes factor calculation: p < epsilon and p > 1-epsilon
#'				 	will be respectively set to epsilon and 1-epsilon to avoid NaN calculations
#' @return A vector of p-values and Pearson correlations (one element per trait), a data frame of 
#'		   standardised and logged datapoints for doing regression, and posterior probabilities /
#'		   Bayes factors from BMA analyses
#' @examples
#' # Sample a birth-death tree with 100 taxa
#' time.tree <- ape::rphylo(birth=10, death=5, n=100)
#' ntips <- length(time.tree$tip.label)
#'
#' # Simulate traits down the tree under Brownian motion
#' traits1 <- simulateTrait(time.tree)
#' traits1.leaf <- traits1[1:ntips]
#'
#' # Simulate substitutions that have a positive association with traits
#' sim.result <- simulateSubstitutions(time.tree=time.tree, 
#'										beta=2, 
#'										theta=10, 
#'										traits=traits1, 
#'										method="TD", 
#'										number.of.subst=10000)
#' subst.tree <- sim.result$subst.tree.est
#'
#' # Sample taxon pairs within a window of (0.01, 0.2) time units
#' pairs.df <- sampleTaxonPairs(subst.tree=subst.tree, 
#'								covariate=traits1.leaf, 
#'								window.tree=time.tree, 
#'								dist.min=0.01, dist.max=0.2)
#'
#' # Perform linear regression on the pairs
#' ppc <- PPC.test(pairs.df)
#' p.value <- ppc$p.var
#' pearson <- ppc$rho.var
#' inclusion.prob <- ppc$bma.probs
#' bayes.factor <- ppc$bma.bf
#'
#' # Plot the pairs. Should see a positive trend under these parameters
#' data.df <- ppc$data.df
#' plot(data.df$trait, data.df$distance.response, xlab="Trait contrast", ylab="Distance contrast")
#' @export
PPC.test = function(pairs.df, standardise=2, logY=TRUE, extreme.value.threshold=0, prior.weight=0.5, OR=1000, epsilon=1e-6){


	d.all = c(pairs.df$d1, pairs.df$d2)
	if (all(d.all == 0)) {
		warning("All distances are 0")
		result = list(p.var=NA, rho.var=NA, t.values=NA, n=nrow(pairs.df), data.df=NA)
		return(result)
	}

	if (logY & any(d.all < 0)) {
		warning("Some terms are negative! Try with logY=FALSE")
		result = list(p.var=NA, rho.var=NA, t.values=NA, n=nrow(pairs.df), data.df=NA)
		return(result)
	}

	
	if (!logY){
		pseudocount=0
	} else{ 

		pseudocount=0
		if (any(d.all == 0)){
			d.all = d.all[d.all > 0]
			pseudocount = min(d.all)
		}
	}


	# Dont standardise
	if (standardise == 0){
		s = 1
	}

	# Ivan 2021: "Contrasts in dS and contrasts in trait values are standardized by the 
	# 						squared root of the total genetic divergence between the two species in a pair."	
	if (standardise == 1){
		s = sqrt(pairs.df$d1 + pairs.df$d2 + 2*pseudocount)
	}

	# Welch and Waxman
	if (standardise == 2){
		s = sqrt(pairs.df$mrca.dist)
	}


	if (standardise == 3){
		s = pairs.df$mrca.dist
	}


	# Prepare response (distances)
	y = pairs.df$d1 - pairs.df$d2
	if (logY){
		y = log(pairs.df$d1+pseudocount) - log(pairs.df$d2+pseudocount)
	}
	y = y / s



	# Data frame
	data.df = data.frame(tip1=pairs.df$tip1, tip2=pairs.df$tip2, distance.response=y)


	# Traits
	varnames = colnames(pairs.df)
	varnames = varnames[varnames != "tip1" & varnames != "tip2" &varnames != "d1" & varnames != "d2" & varnames != "mrca.dist"]
	varnames = unique(gsub("[.]tip(1|2)$", "", varnames))
	for (var in varnames){
		x = pairs.df[,paste0(var, ".tip1")] - pairs.df[,paste0(var, ".tip2")]
		x = x / s
		data.df[,var] = x
	}

	
	


	# Remove extreme values
	if (extreme.value.threshold > 0){
		yz = abs(scale(data.df$distance.response))
		#xz = abs(scale(x))
		remove = which(yz > extreme.value.threshold) # Remove anything more than a few sd from the mean
		if (length(remove) > 0){
			#cat(paste("Removing", length(remove), "extreme values\n"))
			data.df = data.df[-remove,]
		}
	}




	# One linear model on each trait
	p.var = numeric(0)
	rho.var = numeric(0)
	t.values = numeric(0)
	for (var in varnames){
		result = cor.test(data.df$distance.response, data.df[,var])
		p = result$p.value
		p.var[var] = p
		rho.var[var] = cor(data.df$distance.response, data.df[,var])

	}



	# Univariate BMA
	bma.all = list()
	bma.probs = numeric(0)
	bma.bf = numeric(0)
	for (var in varnames){
		bma =  tryCatch({
			BMA::bic.glm(x = data.frame(trait=data.df[,var]), y = data.df$distance.response, glm.family=gaussian(), prior.param=prior.weight, OR=OR)
		}, error = function(e) {
			warning(paste("BMA error caught:", e))
		  	NA
		})
		bma.all[[var]] = bma
		if (is.na(bma[1])){
			bma.probs[var] = NA
			bma.bf[var] = NA
		}else{
			p =  bma$postprob[bma$label == "trait"]
			bma.probs[var] = p
			bf = (p/(1-p)) / (prior.weight / (1-prior.weight)) 
			bma.bf[var] = log(bf, 10)
		}
	}


	list(p.var=p.var, rho.var=rho.var, n=nrow(pairs.df), data.df=data.df, bma.probs=bma.probs, bma.bf=bma.bf, bma=bma.all)

}





