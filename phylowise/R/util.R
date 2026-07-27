


#' Build a distance matrix from a tree (time tree or substitution tree)
#'
#'
#' @param tree a tree (phylo object)
#' @return A distance matrix, where dij is the mean distance from leaf i and j to their ancestor 
getDistanceMatrix = function(tree) {
	
 	if (!inherits(tree, "phylo")) {
		cat("Input must be a phylo object!\n")
		return(NA)
	}
  
  tips = tree$tip.label
  n = length(tips)
  
  # Node depths measured from tips
  node_depths_tips = ape::node.depth.edgelength(tree)
  
  dists = matrix(0, nrow=n, ncol=n)
  rownames(dists)=tips
  colnames(dists)=tips
  
  row.num=1
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      mrca_node = getMRCA(tree, c(i, j))
      distance = sum(node_depths_tips[c(i,j)] - node_depths_tips[mrca_node]) / 2
      dists[i,j] = distance
      dists[j,i] = distance
    }
  }
 
  dists

}




