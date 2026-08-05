# Phylowise


An R package for finding associations between biological traits and molecular rates on a phylogeny.

The method and code are currently in prerelease.

## Description of method

PPC (Phylogenetic Pairwise Contrast) is a phylogenetic comparative method for finding associations between biological traits and molecular evolutionary rates. The method samples pairs from a phylogeny such that each pair has non-overlapping edge paths, and can therefore be treated as statistically independent observations. Linear regression is performed on the pair contrasts. This approach is similar to phylogenetically independent contrasts (PIC) but without reconstructing the traits at internal nodes, and is better suited for finding trait-rate associations than phylogenetic generalised least squares (PGLS).





## Installation



## Directories
```
cases -- R scripts for the two case studies (plants and mammals)
phylowise -- source code for the R package
simulation -- R scripts to reproduce the simulation studies
```


## References




