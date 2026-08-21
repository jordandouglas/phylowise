# Phylowise


An R package for finding associations between biological traits and molecular rates on a phylogeny.

The method and code are currently in prerelease.

## Description of method

PPC (Phylogenetic Pairwise Contrast) is a phylogenetic comparative method for finding associations between biological traits and molecular evolutionary rates. The method samples pairs from a phylogeny such that each pair has non-overlapping edge paths, and can therefore be treated as statistically independent observations. Linear regression is performed on the pair contrasts. This approach is similar to phylogenetically independent contrasts (PIC) but without reconstructing the traits at internal nodes, and is better suited for finding trait-rate associations than phylogenetic generalised least squares (PGLS).





## Installation

phylowise can be installed through CRAN. In the R terminal, type:


```
install.packages("phylowise")
library(phylowise)

```



## Directories

`phylowise/` - R and C++ code for the phylowise package

`simulation/` - R scripts to test PPC and PGLS using simulation studies (false positives and statistical power)

`cases/` - R scripts for the two biological case studies (heights in flowering plants and life history in mammals)




## References



Jordan Douglas and Lindell Bromham. Searching for patterns in rate of molecular evolution using phylogenetic pairwise contrasts (2026). BiorXiv. [https://doi.org/10.64898/2026.08.13.744736](https://doi.org/10.64898/2026.08.13.744736).
