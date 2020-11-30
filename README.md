# CoMiRe package

This repository contains the R package CoMiRe that implements the **Convex Mixture Regression** models first introduced in

- Canale, A., Durante, D., and Dunson, D. B., (2018), [Convex Mixture Regression for Quantitative Risk Assessment](https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12917), __Biometrics__, 74, 1331-1340

Install the R package `CoMiRe_VERSION.tar.gz` with `R CMD INSTALL CoMiRe_VERSION.tar.gz` or

```R
install.packages("devtools")
devtools::install_github("tonycanale/CoMiRe/CoMiRe")
```

This packge is an upgrade of the CoMiRe package developed for the [paper in Biometrics](https://onlinelibrary.wiley.com/doi/full/10.1111/biom.12917). The old version reproducing the results of the paper compiled under R version 3.4 is manteined in the github repository [CoMiRe_BiometricsPaper](https://github.com/tonycanale/CoMiRe_BiometricsPaper). If you landed here after clicking on the link in the supplementary materials of the paper, you may want to check [this page](https://github.com/tonycanale/CoMiRe_BiometricsPaper). 

A  [tutorial](Tutorial.md) illustrating the different models and the package usage is also available. 
