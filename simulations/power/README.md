# Statistical power

This folder contains the scripts for assessing the statistical power of comparative methods based on simulated datasets.

The simulations consist of three components:

1. Time trees. These are all precomputed and stored in the `../newick` folder.
2. Traits. These are simulated under the tree using Brownian motion.
3. Substitution counts and rates. These are simulated down the time tree with a positive association on traits (controlled by strength beta) under an Ornstein-Uhlenbeck clock model.


The following comparative methods will test for an association between rates and traits. We would hope that power increases with taxon count, substitution count, and beta.

1. Phylogenetic pairwise contrasts (PPC)
    1. Nesting, no maximisation
    2. Nesting,  maximisation
2. Ordinary least squares (OLS)
3. Phylogenetic generalised least squares (PGLS)
    1. Brownian
    2. Pagel
    3. Grafen
    4. Martins
    5. Blomberg
   


## Requirements

We will be using the following R packages

- phylowise
- ape
- nlme
- IDPmisc

## Use

To deploy the validation protocol, make the following call

```
Rscript runPower.R power1.tsv
```
where `power1.tsv` is the output file. You can run multiple replicates in parallel, making sure to save each to a different file. Make sure that all outfiles follow the pattern `power*.tsv` so they can be easily plotted.

Optionally, you can skip OLS/PGLS and only do PPC using
```
Rscript runPower.R power.ppc.1.tsv TRUE
```

or skip PPC and only do OLS/PGLS using
```
Rscript runPower.R power.pgls.1.tsv FALSE TRUE
```



To plot the results, use

```
Rscript plotPower.R
```
which will read in all files that follow the pattern `power*.tsv` and generate a file called `power.png`. A linear model will be fit to each comparative method to show how power responds to beta, taxon/pair count, substitution count.

Ideally, this pipeline should be run for at least 1000 relicates per comparative method before any conclusions can be made. We performed 10,000 in the main study, which took several days across multiple parallel threads. The plotting script will print the number of observations in each panel in the figure. If there are more than 10,000 then some will be randomly discarded. 

