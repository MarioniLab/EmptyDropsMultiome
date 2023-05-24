# EmptyDrops-multiome
Framework for statistically powerful and accurate detection of nuclei-containing droplets in single-cell GEX+ATAC multiome data. The method builds on a cell calling method for droplet-based scRNA data called [EmptyDrops](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y). It can deal with diverse samples (from highly homogeneous to highly heterogeneous) by creating the ATAC and RNA profile of the ambient noise and then testing each droplet for deviations from it.


## Installation

```

conda create -n eDm_env -c conda-forge r-base=4.3.0 -y
conda activate eDm_env
mamba install -c conda-forge r-devtools
R
> library(devtools)
> if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install("DropletUtils")
> devtools::install_github("MarioniLab/emptyDrops_multiome",
                         ref="main",
                         auth_token = "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"
                         )
# update all, but don't worry about curl package

```



## Vignette in R

```
library(Seurat)
library(eDv3)
sce <- Read10X_h5("/home/jovyan/my-conda-envs/eDv3_test_data/raw_feature_bc_matrix.h5")
rna <- sce[["Gene Expression"]]
atac <- sce[["Peaks"]]
eD.out <- emptydrops_multiome(count_matrix_rna=rna, count_matrix_atac=atac)
print("the number of cells detected is: ")
print(sum(eD.out$FDR_multi<0.001 & ! is.na(eD.out$FDR_multi)))


```

