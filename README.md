# circFunEnrich
**circFunEnrich** is a shiny app.

It intends to predict the functions of human circRNA through genes with expression correlation. 
It contains three parts of 
* circBank annotation
* canonical function enrichment analysis(CFEA) 
* gene set enrichment analysis(GSEA).

## OS

What my computer is Windows 10, X86_64

## Depends
As a shiny app, you have to install the **R** and **Rstudio**.

* R (>= 3.6.0)

## Install
```
library(devtools)
install_github("gnilihzeux/circFunEnrich")
```
## Run
```
library(circFunEnrich.beta1.1)
launchApp()
```
## Suggestions
Some functions of the app depends on ***clusterProfiler***, so a more memory is needed.

What my ROM size is 8G.
