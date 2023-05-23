# emptyDrops_multiome
Framework for statistically powerful and accurate detection of nuclei-containing droplets in single-cell GEX+ATAC multiome data. The method builds on a cell calling method for droplet-based scRNA data called [EmptyDrops](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y). It can deal with diverse samples (from highly homogeneous to highly heterogeneous) by creating the ATAC and RNA profile of the ambient noise and then testing each droplet for deviations from it.


## Installation

```

# Install development version:
library(devtools)
devtools::install_github("MarioniLab/emptyDrops_multiome", ref="main") 


```

this works
