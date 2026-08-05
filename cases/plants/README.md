# Flowering plants

In this case study we will explore the negative association between height and substitution rate in flowering plants, using the dataset and code from Lanfear et al 2013. All of the contents in the `data` folder were downloaded from their [Dryad repository](https://datadryad.org/dataset/doi:10.5061/dryad.43mg3).
	


## Requirements

We will be using the following R packages

- phylowise
- ape


## Running the code


Launch the R terminal and copy the code from `plants.R`. 

This code will first load in the 1000 bootstrapped trees, estimated per-family rates, and trait data (heights, uv, and temperature). 

PPC will be performed across a sample of these trees, exploring a range of window sizes, with and without nesting and trait maximisation. Running this search across 10 trees takes around 30 seconds. 

Then, the resulting p-values and Pearson correlations can be plotted for any given PPC configuration. In most cases, there is a negative and signifcant correlation between height and rate, provided that there were enough taxon pairs sampled by PPC.


## References


Lanfear, R., Ho, S., Jonathan Davies, T. et al. Taller plants have lower rates of molecular evolution. Nat Commun 4, 1879 (2013). https://doi.org/10.1038/ncomms2836
