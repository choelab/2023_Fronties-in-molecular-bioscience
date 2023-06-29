DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)

base::gc()
#base::rm(list = ls()) # Clear the environment
options(warn=-1)

load.lib <- c("Seurat","SingleCellExperiment","dplyr","tidyverse","data.table","scales","RColorBrewer","ggplot2","SingleCellExperiment","devtools","remotes")

if (!requireNamespace(load.lib, quietly = TRUE))
    install.packages(load.lib)

load.lib2 <- c("cowplot","ggpubr","Matrix","viridis","tibble","ggrepel","forcats","readxl","parallel")

if (!requireNamespace(load.lib2, quietly = TRUE))
    install.packages(load.lib2)


sapply(load.lib,suppressPackageStartupMessages(require),character.only = TRUE)

require.lib <- c("pbapply","reshape2","stringr","ggvenn","stringi","BiocManager","RColorBrewer","xlsx","magrittr")

if (!requireNamespace(require.lib, quietly = TRUE))
    install.packages(require.lib)

 bioconductor.lib <- c(#"EnhancedVolcano","STRINGdb",
                     "BiocParallel",
                      'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr',"tradeSeq","ComplexHeatmap","Nebulosa","UCell")

if (!requireNamespace(bioconductor.lib, quietly = TRUE))
     BiocManager::install(bioconductor.lib)

if (!requireNamespace("nichenetr", quietly = TRUE))
    devtools::install_github("saeyslab/nichenetr")

if (!requireNamespace("seurat-wrappers", quietly = TRUE))
     remotes::install_github('satijalab/seurat-wrappers')

if (!requireNamespace("monocle3", quietly = TRUE))
     devtools::install_github('cole-trapnell-lab/monocle3')
