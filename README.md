# GONetwork

The purpose of GoNetwork is to calculate networks based on GO annotation and to identify clusters that are differentially enriched between conditions based on expression data. Networks created using this program are easily plotted with molecular visualization software platforms such as Cytoscape. 

To install dependencies:
```r
source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")
```

To install GONetwork in R:
```r
devtools::install_github("saralinker/GONetwork")
```

For detailed information on packages read the GONetwork vignette available in vignettes/


