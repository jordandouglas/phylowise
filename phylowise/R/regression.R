


#' Phylogenetic pairwise contrast test on a series of taxon pairs, using ordinary least squares and Bayesian model averging (BMA).
#' This method does not do multivariate regression, but it returns the data frame that can be used for it
#'
#' @param pairs.df data frame of taxon pairs
#' @param logY should we take the logarithm of the evolutionary distance? default and recommended setting: true
#' @param standardise should we standardise the evolutionary distance? 0 for no, 1 for standardise by distance, 2 for standardise by mrca age. default and recommended setting: 2
#' @param extreme.value.threshold remove any observations more than this many standard deviations away from the mean (distance); set to zero for no extreme value removal
#' @param prior.weight trait inclusion prior probbaility for Bayesian model averaging (see BMA::bic.glm)
#' @param OR Occam's window for Bayesian model averaging (see BMA::bic.glm)
#' @param epsilon precision of BMA probabilities for Bayes factor calculation: p<epsilon and p>1-epsilon will be respectively set to epsilon and 1-epsilon to avoid NaN calculations
#' @return Vector of p-values and Pearson correlations (one element per trait), a data frame of standardised and logged datapoints for doing regression, and posterior probabilities / Bayes factors from BMA analyses
#' @export
PPC.test = function(pairs.df, standardise=2, logY=TRUE, extreme.value.threshold=0, prior.weight=0.5, OR=1000, epsilon=1e-6){


	d.all = c(pairs.df$d1, pairs.df$d2)
	if (all(d.all == 0)) {
		cat(paste("All distances are 0!\n"))
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
			cat(paste("Removing", length(remove), "extreme values\n"))
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
			cat(paste("BMA error:", e, "\n"))
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





