# Mammals

In this case study we will explore associations between substitution rates and life history traits in mammals. This will be based on phylogenies from Douglas and Bromham 2025. Those trees were inferred under a multispecies coalescent tree with branch lengths expressed as number of synonmyous or non-synonymous changes (estimated via the [BeastMap](https://github.com/jordandouglas/BeastMap) package for BEAST 2).  

The `data` folder contains a posterior distribution of these trees + substitutions (downsampled to a small number to improve processing times) and species traits sourced from Pantheria (Jones et al. 2009) and AnAge (Magalhães et al. 2009).
	


## Requirements

We will be using the following R packages

- phylowise
- ape
- randomForest
- BMA
- corrplot

## Running the code


Launch the R terminal and copy the code from `mammals.R`.  A figure will be generated from this analysis, which will be similar to Figure 6 of the main PPC article, with differences coming from random seeds.

This script will by default run just 10 iterations on 10 random trees from the posterior, applying PPC between all combinations of variables (nS, nS, mass, longetivity, age of sexual maturity, and fecundity). This should take a few minutes to complete. In the main article, we performed 1000 iterations.


## References

Douglas, J., & Bromham, L. (2025). Reconstructing substitution histories on phylogenies, with accuracy, precision, and coverage. bioRxiv, 2025-12. 

Jones, K. E., Bielby, J., Cardillo, M., Fritz, S. A., O'Dell, J., Orme, C. D. L., ... & Purvis, A. (2009). PanTHERIA: a species‐level database of life history, ecology, and geography of extant and recently extinct mammals: Ecological Archives E090‐184. Ecology, 90(9), 2648-2648.

De Magalhães, J. P., & Costa, A. J. (2009). A database of vertebrate longevity records and their relation to other life‐history traits. Journal of evolutionary biology, 22(8), 1770-1774.
