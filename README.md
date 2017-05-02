# GONetwork

GONetwork is an R toolkit that detects statistically significant functional annotation differences between gene lists and has direct application to RNA-seq, ChIP-seq, and GWAS experiments. Through a short workflow GONetwork: (1) calculates a gene network based on functional overlap, (2) clusters genes within that defined functional space, and (3) statistically tests differential enrichment within each cluster to determine functional annotation that is unique to a given condition.

Currently, standard functional analyses are restricted by the requirement of exact overlap of a gene list with a dictionary of predefined terms. This strategy is limited by the lack of complete annotation. Furthermore, it is limited by reduced power in small genes sets which often results in broad, biologically ambiguous terms. Contrary to these current methods that require functional equivalence, GONetwork compares genes via functional similarity. Functional similarity initially broadens the space of terms identified and then narrows in on the highly specific terms within a gene set. This alternative strategy returns a highly robust and informative set of enriched functional terms. Within the terms that are enriched within a complete dataset, GONetwork then tests differential enrichment of each annotation cluster based on the condition from which the gene was identified. GONetwork returns to the user a highly specific set of functional terms that are relevant to the study as well as those that are specific to a given condition. This enables researchers to narrow in on the informative functional groupings present within their dataset.

To install GONetwork in R:
```r
devtools::install_github("saralinker/GONetwork")
```

For detailed information on packages read the GONetwork vignette available in vignettes/


