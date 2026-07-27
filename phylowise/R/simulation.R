


#' Simulate a trait down a tree using Brownian motion (BM)
#'
#'
#' @param time.tree a binary rooted tree (phylo object)
#' @param sigma standard deviation for BM
#' @param root.value trait value at the root of the tree
#' @return A vector of traits, one element for each node in the tree, ordered by node number
#' XXX TODO EXAMPLES XXX
#' XXX
#' XXX
#' XXX
#' XXX
#' XXX
#' XXX
#' XXX
#' @export 
simulateTrait = function(time.tree, sigma=1, root.value=0) {


	# Validation
	if (!inherits(time.tree, "phylo")) {
		stop("time.tree must be a phylo object!")
	}

		if (!is.rooted(time.tree) || !is.binary(time.tree)){
		stop("time.tree must be rooted and binary!")
	}


	sigma = sigma[1]
	if (sigma <= 0) {
		stop("sigma must be positive")
	}


	nnodes = time.tree$Nnode


	ntaxa = length(time.tree$tip.label)
  edge = time.tree$edge
  edge.length = time.tree$edge.length
  node_times = ape::node.depth.edgelength(time.tree)
  
  # Find total time from root
  max_time = max(node_times)
  
  # Preallocate trait values
  rootNr = ntaxa+1
  traits = numeric(max(time.tree$edge))
  traits[1:length(traits)] = NA
  traits[rootNr] = root.value  # root node

  for (i in 1:nrow(edge)) {
    parent = edge[i, 1]
    child = edge[i, 2]
    branch_length = edge.length[i]
    
    # Time since root at start of branch
    t_start = node_times[parent]
    t_end = node_times[child]
    
  	# Brownian noise
  	noise = rnorm(1, mean = 0, sd = sigma * sqrt(branch_length))


	 	# Trait at child = trait at parent + mean + noise
 	 	traits[child] = traits[parent] + noise
   
  }

  
  #names(traits) = c(time.tree$tip.label, time.tree$node.label)
  traits

}





#' Simulate a substitution count down each branch of a tree. First, the branch rates are sampled from a correlated, uncorrelated, or OU process that can be dependent on traits. 
#' Then the number of substitutions along each branch is sampled from a Poisson distribution.
#'
#' @param time.tree a binary rooted tree (phylo object)
#' @param nu variance scale of autocorrelated clock or TD
#' @param sigma standard deviation of UCLN
#' @param theta theta term for OU process in TD
#' @param beta effect size of traits on rates (if traits is not NULL)
#' @param traits one trait per node; leave as NULL if rates are conditionally independent of traits
#' @param number.of.subst expected number of substitutions per unit of time
#' @param method clock model may be uncorrelated lognormal (UCLN), autocorrelated lognormal (AC), or trait dependent (TD)
#' @return Two trees, with branch lengths set to either subst. rates or counts, and two vectors, one of true branch rates and one of estimated branch rates (i.e., count divided by time, which can evaluate to zero)
#' XXX TODO EXAMPLES XXX
#' XXX
#' XXX
#' XXX
#' XXX
#' XXX
#' XXX
#' XXX
#' @export 
simulateSubstitutions = function(time.tree, nu=0.5, sigma=0.5, beta=0.5, theta=1, traits=NULL, number.of.subst=1000, method=c("AC", "UCLN", "TD")) {


	# Validation
	if (!inherits(time.tree, "phylo")) {
		stop("time.tree must be a phylo object!")
	}

	if (!is.rooted(time.tree) || !is.binary(time.tree)){
		stop("tree must be rooted and binary!")
	}


	valid.methods = c("UCLN", "AC", "TD")
	method = as.character(method[1])
	if (!(any(method==valid.methods))) {
		stop(paste("Please select one of the following methods:", paste(valid.methods, collapse=", ")))
	}
	if (method == "TD" & is.null(traits)){
		stop("Please provide traits in order to use the TD method")
	}


	sigma = as.numeric(sigma[1])
	if (sigma <= 0) {
		stop("sigma must be positive")
	}


	nu = as.numeric(nu[1])
	if (nu <= 0) {
		stop("nu must be positive")
	}

	number.of.subst = number.of.subst[1]
	if (number.of.subst <= 0) {
		stop("theta must be positive")
	}

	if (!is.null(traits)){

		nnodes = length(time.tree$tip.label) + time.tree$Nnode
		beta = as.numeric(beta[1])
		if (length(traits) != nnodes){
			stop(paste0("traits must be the same length as the node count ", length(traits), "!=", nnodes, "\n"))
		}
		traits = as.numeric(traits)
		if (any(is.na(traits))){
			stop("NA values detected in traits!")
		}

	}


	subst.tree.true = time.tree
	subst.tree.est = time.tree

	ntaxa = length(time.tree$tip.label)
  edge = time.tree$edge
  edge.length = time.tree$edge.length
  node_times = ape::node.depth.edgelength(time.tree)

  
  # Find total time from root
  max_time = max(node_times)
  
  # Preallocate relative rates
  rootNr = ntaxa+1
  log.rates = numeric(max(time.tree$edge))
  log.rates[1:length(log.rates)] = NA
  log.rates[rootNr] = 0  # root node
  rates.est = numeric(length(log.rates))
  
  for (i in 1:nrow(edge)) {
    parent = edge[i, 1]
    child = edge[i, 2]
    branch_length = edge.length[i]


    if (method == "UCLN"){
    	expected = -0.5*sigma^2 # Mean of 1 in real-space
    	variance = sigma^2
    }


    if (method == "AC"){
    	variance = nu*branch_length
    	expected = log.rates[parent] - 0.5*variance
    	
    }

    if (method == "TD"){

    	variance = sigma^2/(2*theta) * (1 - exp(-2*theta*branch_length))
    	p = exp(-theta*branch_length)
    	trait.mean = log(mean(exp(c(traits[child], traits[parent]))))
			mu = beta * trait.mean
			expected = mu*(1-p) + log.rates[parent]*p - 0.5*variance

			#print(paste(mu, expected, variance, p, log.rates[parent], trait.mean))
    }
   

    # Rate at child
    r = rnorm(1, mean=expected, sd=sqrt(variance))
    log.rates[child] = r



    # Turn rates into counts
    num = rpois(1, exp(log.rates[child]) * branch_length * number.of.subst)
    subst.tree.true$edge.length[i] = exp(log.rates[child]) * branch_length # Subst per site
    subst.tree.est$edge.length[i] = num # Number of subst
    rates.est[child] = num / (branch_length * number.of.subst) # Subst per site, estimated with Poisson error


  }



  list(subst.tree.true=subst.tree.true, subst.tree.est=subst.tree.est, node.rates.true=exp(log.rates), node.rates.est=rates.est)

}


