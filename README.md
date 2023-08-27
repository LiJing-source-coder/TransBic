# TransBic
TransBic is an R toolkit for extracting bucket trend-preserving biclusters for analysis of gene expression data
# Installation
```
$ R
> library("devtools")
> devtools:install_github("LiJing-source-coder/TransBic")
```
Now TransBic can be loaded in R:
```
> library(TransBic)
```
# Input
The input to TransBic is a gene expression matrix with rows being genes and columns being conditions/samples, and values could be measured by DNA microarrays or RNA-seq.
# Usage 
We provide an example to run TransBic as follows.
# Example
```
> library(TransBic)
> data(sim)
> res_TransBic <- TransBic(sim)
```
