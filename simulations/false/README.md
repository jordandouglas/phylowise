# False positive rates

This folder contains the scripts for assessing the false positive rate (Type I error rate) of comparative methods based on simulated datasets.

The simulations consist of three components:

1. Time trees. These are all precomputed and stored in the `../newick` folder.
2. Traits. These are simulated under the tree using Brownian motion.
3. Substitution counts and rates. These are simulated down the time tree independently of the traits. By default, these are simulated under an autocorrelated clock model, where the log-rate wanders down the tree via Brownian motion.


The following comparative methods will test for an association between rates and traits. Because traits are independent of rates, we would hope for all p-values to be distributed as Uniform(0,1).

1. Phylogenetic pairwise contrasts (PPC)
- No nesting, no maximisation
- No nesting, maximisation
- Nesting, no maximisation
- Nesting,  maximisation
2. Ordinary least squares (OLS)
3. Phylogenetic generalised least squares (PGLS)
- Brownian
- Pagel
- Grafen
- Martins
- Blomberg

These will be run across different error models of trait and rate.
	



## Requirements

We will be using the following R packages

- phylowise
- ape
- nlme

## Use

To deploy the validation protocol, make the following call

```
Rscript runFP.R false1.tsv
```
where `false1.tsv` is the output file. You can run multiple replicates in parallel, making sure to save each to a different file. Make sure that all outfiles follow the pattern `false*.tsv` so they can be easily plotted.

Optionally, you can skip PGLS and only do PPC using
```
Rscript runFP.R false.ppc.1.tsv TRUE
```

or skip PPC and only do PGLS using
```
Rscript runFP.R false.pgls.1.tsv FALSE TRUE
```



To plot the results, use

```
Rscript plotFP.R
```
which will read in all files that follow the pattern `false*.tsv` and generate a file called `FP.png`. By default, we will only consider trials that generated at least 20 observations (taxa or pairs) and will consider a Kolmogorov-Smirnov p-value threshold of 0.01 as significant evidence that the comparative method p-value is not distributed as Uniform(0,1).  

Ideally, this pipeline should be run for at least 1000 relicates per comparative method before any conclusions can be made. We performed 10,000 in the main study, which took several days across multiple parallel threads. The plotting script will print the number of observations in each histogram in the figure. If there are more than 10,000 then some will be randomly discarded. 

