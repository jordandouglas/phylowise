

#' Sample a set of taxon pairs from the tree. These taxa will have non-overlapping edges
#' between their paths and can this be treated as statistically independent. The algorithm  
#' iteratively searches for two random taxa that satisfy sampling requirements until there 
#' are no more valid pairs. Each pair must descend from an MRCA with 
#' dist.min <= tMRCA <= dist.max. c++ is used to speed up the runtime of this code.
#'
#'
#' @param subst.tree a binary rooted tree, with branch lengths in units of change (phylo object)
#' @param response provide a vector as a response trait instead of genetic distances in the tree (optional)
#' @param covariate a data frame of traits at the tips of the tree (rows are taxa, columns are traits);
#'					rows should be in the same order as tips in the tree
#' @param dist.min minimum distance that two tips must be apart from their ancestor
#' 					(mean of both distances), on window.tree
#' @param dist.max maximum distance that two tips must be apart from their ancestor 
#'					(mean of both distances), on window.tree
#' @param nested can taxon pairs be nested with each other?
#' @param maximise should we maximise the trait difference?
#' @param youngest take the youngest pair at each step (and therefore increase the number of pairs)?
#' @param window.tree tree that is the basis for building the distance matrix, if it is not provided; 
#'						should have same taxa as 'tree'
#' @param distance.matrix provide an n x n distance matrix rather than recalculate from scratch
#' @param verbose print some statements along the way
#' @return A data frame, where each row is a pair of taxa. 
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
#'										beta=1, 
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
#' # If iterating this, we can precompute the distance matrix to save time
#' dmat <- getDistanceMatrix(time.tree)
#' for (i in 1:100){
#'     pairs.df <- sampleTaxonPairs(subst.tree=subst.tree, 
#'									covariate=traits1.leaf, 
#'									distance.matrix=dmat, 
#'									dist.min=0.01, 
#'									dist.max=0.2)
#' }
#'
#' # We can also do this with traits-vs-traits rather than rates-vs-traits
#' traits2 <- simulateTrait(time.tree) # Unassociated with traits1
#' traits2.leaf <- traits2[1:ntips]
#' pairs.df <- sampleTaxonPairs(subst.tree=subst.tree, 
#'								response=traits1.leaf, 
#'								covariate=traits2.leaf, 
#'								window.tree=time.tree, 
#'								dist.min=0.01, 
#'								dist.max=0.2)
#' @export
sampleTaxonPairs = function(subst.tree, covariate, dist.min, dist.max=Inf, response=NULL, nested=TRUE, maximise=TRUE, youngest=FALSE, window.tree=subst.tree, distance.matrix=NULL, verbose=FALSE){



	response.is.tree = is.null(response)

	
	# Validation
	if (!inherits(subst.tree, "phylo")) {
		stop("subst.tree must be a phylo object!")
	}
	if (!ape::is.rooted(subst.tree) || !ape::is.binary(subst.tree)){
		stop("subst.tree must be rooted and binary!")
	}



	tips = subst.tree$tip.label
 	n = length(tips)
 	nnodes = n + subst.tree$Nnode

 	if (subst.tree$Nnode != window.tree$Nnode){
 		stop("subst.tree and window.tree must be the same size!")
 	}



 	# Build a tree pointer object in c++ for efficiency
 	tree_ptr = rcpp_buildTree(subst.tree$edge, n+1)

 	if (is.null(nrow(covariate))){
 		covariate = data.frame(trait = as.numeric(covariate))
 	}

 	if (nrow(covariate) < n) {
		stop(paste0("covariate should be at least the same size as the number of tips(", n, ")"))
	}
	ntraits = ncol(covariate)
	covariate = covariate[1:n,]
	if (ntraits == 1){
		covariate = data.frame(trait = as.numeric(covariate))
	}
	

	if (is.null(distance.matrix)){

		if (!inherits(window.tree, "phylo")) {
			stop("window.tree must be a phylo object!")
		}

		if (!ape::is.rooted(window.tree) || !ape::is.binary(window.tree)){
			stop("tree must be rooted and binary!")
		}

		distance.matrix = phylowise::getDistanceMatrix(window.tree)
	}else{
		if (nrow(distance.matrix) != n || ncol(distance.matrix) != n){
			stop(paste0("distance.matrix must have ", n, "rows and columns to match the tree"))
		}
	}


	if (!response.is.tree && length(response) != n){
		stop(paste0("please give one response element for each tip!", length(response), "!=", n ))
	}


	# Filter out all distances under the minimum
	distance.matrix[distance.matrix < dist.min] = Inf
	distance.matrix[distance.matrix > dist.max] = Inf
	diag(distance.matrix) = Inf

	
	leaves.df = data.frame(taxon=tips)
	for (var in colnames(covariate)){
		leaves.df[,var] = as.numeric(covariate[,var])
	}


	# Prepare the response variable
	if (response.is.tree){
		response = ape::node.depth.edgelength(subst.tree)
	}else{
		# response = response
	}
	names(response) = subst.tree$tip.label



	# Filter out all distances with an NA in any row
	na.rows = unique(c(which(apply(leaves.df, 1, anyNA)), which(is.na(response))))
	for (i in na.rows){
	    distance.matrix[i,] = Inf
	    distance.matrix[,i] = Inf
	    if (verbose){
			cat(paste("Removing", leaves.df[i,"taxon"], "because there are NA trait values\n"))
		}
	}

	

	# Sample pairs randomly, but then mark the whole clade as dirty
	max.npairs = floor(n/2)
	pairs.df = data.frame(tip1 = rep("", max.npairs), tip2 = "", mrca.dist=0, d1 = 0, d2 = 0)
	pair.num = 0

	for (var in colnames(covariate)){
		pairs.df[,paste0(var, ".tip1")] = 0
		pairs.df[,paste0(var, ".tip2")] = 0
	}




	# Repeat until everthing is dirty
	while (TRUE){


		candidates.all = which(distance.matrix < Inf, arr.ind=T) # Very slow line to run


		
		if (nrow(candidates.all) == 0) {
			break
		}

		# Get the candidates to choose from
		if (youngest){
			candidates = which(distance.matrix == min(distance.matrix), arr.ind=T)
		}else{
			candidates = candidates.all
		}

		candidates.random = sample(1:nrow(candidates), 1)

		tip.idx1 = candidates[candidates.random,1]
		tip.idx2 = candidates[candidates.random,2]

		
		tip1 = rownames(distance.matrix)[tip.idx1]
		tip2 = colnames(distance.matrix)[tip.idx2]
		mrca = rcpp_getMRCA(tree_ptr, tip.idx1, tip.idx2)


		if (maximise){ 



			# Maximise the trait diff of either subtree
			children = subst.tree$edge[subst.tree$edge[,1] == mrca,2]
			node1 = children[1]
			node2 = children[2]

			desc1 = rcpp_getDescendantsFast(tree_ptr, node1, n)
			desc2 = rcpp_getDescendantsFast(tree_ptr, node2, n)

		
			# Remove those that are all dirty
			desc1 = desc1[sapply(desc1, function(d) !any(pairs.df$tip1 == tips[d]) & !any(pairs.df$tip2 == tips[d]))]
			desc2 = desc2[sapply(desc2, function(d) !any(pairs.df$tip1 == tips[d]) & !any(pairs.df$tip2 == tips[d]))]

			traits1 = leaves.df[desc1,-1]
			traits2 = leaves.df[desc2,-1]

			desc1.labels = tips[desc1]
			desc2.labels = tips[desc2]


			max.by.dist=FALSE
			if (max.by.dist){



				r1 = response[desc1]
				r2 = response[desc2]


				# Absolute difference
				d = outer(r1, r2, function(x, y) {
					abs(x - y)
				})

				for (d1 in names(traits1)){
					for (d2 in names(traits2)){
						if (is.na(distance.matrix[d1,d2]) |  distance.matrix[d1,d2] == Inf){
							d[d1,d2] = 0
						}else{
							mrca.age = distance.matrix[d1,d2] 
							d[d1,d2] = d[d1,d2] / mrca.age
						}
					}
				}


			}
			

			# Single trait
			else if (ntraits == 1){


				names(traits1) = desc1.labels
				names(traits2) = desc2.labels

				# Absolute difference
				d = outer(traits1, traits2, function(x, y) {
					abs(x - y)
				})

				for (d1 in desc1){
					for (d2 in desc2){
						if (is.na(distance.matrix[d1,d2]) |  distance.matrix[d1,d2] == Inf){
							t1 = tips[d1]
							t2 = tips[d2]
							d[t1,t2] = NA
						}
					}
				}

			}else{
	

				# Euclidean distance between each pair of rows/taxa
				d = matrix(0, nrow=nrow(traits1), ncol=nrow(traits2))
				rownames(d) = desc1.labels
				colnames(d) = desc2.labels
				for (t1 in 1:nrow(traits1)){
					for (t2 in 1:nrow(traits2)){

						d1 = desc1.labels[t1]
						d2 = desc2.labels[t2]

						if (any(d1 == pairs.df$tip1) | any(d1 == pairs.df$tip2)){
							d[t1,t2] = NA
							next
						}
						if (any(d2 == pairs.df$tip1) | any(d2 == pairs.df$tip2)){
							d[t1,t2] = NA
							next
						}

						if (is.na(distance.matrix[d1,d2])){
							d[t1,t2] = NA
							next
						}

						vals1 = traits1[t1,]
						vals2 = traits2[t2,]

						dist = sum((as.numeric(vals1) - as.numeric(vals2))^2) # No need to sqrt
						if (is.na(dist) | dist == Inf){
							dist = NA
						}
						d[t1,t2] = dist
					}
				}
			}


			if (all(is.na(d))){

				for (d1 in desc1){
					for (d2 in desc2){
						distance.matrix[d1,d2] = Inf
						distance.matrix[d2,d1] = Inf
					}
				}
				next
			}


			# Find maximal distance pair
			ij = which(d == max(d, na.rm=T), arr.ind = TRUE)

			#print(d)				
			ij = ij[sample(nrow(ij),1),]
			tip1 = rownames(d)[ij[1]]
			tip2 = colnames(d)[ij[2]]


		}



		# This should not happen but the check is here just to be safe
		if (any(pairs.df$tip1 == tip1) | any(pairs.df$tip2 == tip1)){
			#print("Unexpeced dev error 1")
			next
		}
		if (any(pairs.df$tip1 == tip2) | any(pairs.df$tip2 == tip2)){
			#print("Unexpected dev error 2")
			next
		}

		trait1 = leaves.df[leaves.df$taxon==tip1,-1]
		trait2 = leaves.df[leaves.df$taxon==tip2,-1]


		# Distance summed across all trees that have both tips
		if (response.is.tree){
			d1 = as.numeric(response[tip1] - response[mrca])
			d2 = as.numeric(response[tip2] - response[mrca])
		}else{
			d1 = as.numeric(response[tip1])
			d2 = as.numeric(response[tip2])
		}

	

		if (length(d1) == 0 || length(d2) == 0 || is.na(d1) || is.na(d2)) { #} || d1 == 0 || d2 == 0){
			distance.matrix[tip1,tip2] = Inf
			distance.matrix[tip2,tip1] = Inf
			next
		}



		pair.num = pair.num + 1
		pairs.df$tip1[pair.num] = tip1
		pairs.df$tip2[pair.num] = tip2
		pairs.df$mrca.dist[pair.num] = distance.matrix[tip1,tip2]
		pairs.df$d1[pair.num] = d1
		pairs.df$d2[pair.num] = d2


		if (ntraits == 1){
			pairs.df[pair.num,paste0(var, ".tip1")] = as.numeric(trait1)
			pairs.df[pair.num,paste0(var, ".tip2")] = as.numeric(trait2)
		}else{
			for (var in colnames(covariate)){
				pairs.df[pair.num,paste0(var, ".tip1")] = as.numeric(trait1[var])
				pairs.df[pair.num,paste0(var, ".tip2")] = as.numeric(trait2[var])
			}
		}

		

		distance.matrix[tip1,] = Inf
		distance.matrix[,tip1] = Inf
		distance.matrix[tip2,] = Inf
		distance.matrix[,tip2] = Inf


		if (!nested){

			# Mark the whole clade as dirty by setting distances to Inf
			desc = rcpp_getDescendantsFast(tree_ptr, mrca, n)
			distance.matrix[desc,] = Inf
			distance.matrix[,desc] = Inf


			# Also mark every path that crosses a node ancestral to this mrca as dirty
			rootNr = n+1
			path.to.root = rcpp_nodepathFast(tree_ptr, mrca, rootNr)
			for (m in path.to.root){
				m.children.nodes = subst.tree$edge[subst.tree$edge[,1] == m,2]

				left.leaves = rcpp_getDescendantsFast(tree_ptr, m.children.nodes[1], n)
				right.leaves = rcpp_getDescendantsFast(tree_ptr, m.children.nodes[2], n)
				
				distance.matrix[left.leaves,right.leaves] = Inf
				distance.matrix[right.leaves,left.leaves] = Inf

			}


		}else{



			n1 = which(subst.tree$tip.label == tip1)
			n2 = which(subst.tree$tip.label == tip2)

			# Mark the path from these two leaves to the mrca as dirty
			dirty = rcpp_nodepathFast(tree_ptr, n1, n2)

			for (d in dirty){

				# We can no longer have one from inside and one from outside this clade
				desc1 = rcpp_getDescendantsFast(tree_ptr, d, n)
				not.desc1 = rcpp_getNonDescendantsFast(tree_ptr, d, n)
				distance.matrix[desc1, not.desc1] = Inf
				distance.matrix[not.desc1, desc1] = Inf

			}



		}

	
	}


	# This should never happen unless there is a bug
	if (pair.num > 1){
		all.included = table(c(pairs.df$tip1[1:(pair.num)], pairs.df$tip2[1:(pair.num)]))
		if (length(all.included) > 0 && max(all.included) > 1){
			#print(all.included)
			#print(pairs.df)
			stop("Error: detected multiple copies of the same taxon")
		}
	}

	if (verbose){
		cat(paste("Sampled", pair.num, "pairs\n"))
	}


	if (pair.num == 0){
		pairs.df[0,]
	}else{
		pairs.df[1:pair.num,]
	}


	


}


#' Plot a set of sampled pairs onto a tree
#'
#'
#' @param tree the binary rooted tree used to get samples
#' @param pairs.df a data frame of pairs, obtained using phylowise::sampleTaxonPairs
#' @param edge.col the pairs will be highlighted in this colour
#' @param show.tip.label display tip labels on the tree?
#' @param edge.width edge line width of all branches on the tree
#' @param edge.width.pairs edge line width of paired branches 
#' @param label.cex tree tip label font size, if show.tip.label=TRUE  
#' @return No return value, function is called to make a plot
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
#'										beta=1, 
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
#' # Plot the pairs onto the substitution tree
#' plotPairs(subst.tree, pairs.df, show.tip.label=TRUE)
#'
#' # Plot the pairs onto the time tree
#' plotPairs(time.tree, pairs.df, show.tip.label=TRUE, edge.col="#008cba")
#' @export
plotPairs = function(tree, pairs.df, edge.col="red", show.tip.label=FALSE, edge.width=1, edge.width.pairs=3, label.cex=1){



	# Validation
	if (!inherits(tree, "phylo")) {
		stop("subst.tree must be a phylo object!")
	}
	if (!is.rooted(tree) || !is.binary(tree)){
		stop("subst.tree must be rooted and binary!")
	}



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
	lp = get("last_plot.phylo", envir = ape::.PlotPhyloEnv)

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
		path = ape::nodepath(tree, index1, index2)
		mrca = ape::getMRCA(tree, c(index1, index2))
		path = path[path != mrca]
		for (child in path){

			parent = tree.plot$edge[tree.plot$edge[,2] == child,1]

			# Vertical line
			lines(c(lp$xx[child], lp$xx[child]), c(lp$yy[child],lp$yy[parent]), lwd=edge.width.pairs, col=edge.col)

			# Horizontal shoulder
			lines(c(lp$xx[child], lp$xx[parent]), c(lp$yy[parent], lp$yy[parent]), lwd=edge.width.pairs, col=edge.col)
 

		}


	}

}


