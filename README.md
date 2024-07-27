## BACT

Current BINRES R package: version 1.2 

The R package BINRES is built to implement the nonparametric Bayesian method named BINRES to carry out the region segmentation for a tissue section by integrating all the three types of data generated during the study --- gene expressions, spatial coordinates, and the histology image. BINRES captures the spatial dependence of neighboring spots and does not require a prespecified region number. It also combines the image and the gene expressions whose contribution weights can be flexibly adjusted in a data-adaptive manner. The computationally scalable extension BINRES-fast is developed for large-scale studies. In this package, a partially collapsed Gibbs sampler is carefully designed for Bayesian posterior inference. BINRES can be installed in Windows, Linux, and Mac OS.

For technical details, please refer to our paper currently accepted in *Journal of the American Statistical Association*: Yinqiao Yan and Xiangyu Luo* (2024), "Bayesian integrative region segmentation in spatially resolved transcriptomic studies".

**Note:** The BINRES R package for reproducibility in the paper is version 1.1.



## Prerequisites and Installation

1. R version >= 4.1.3.
2. CRAN package: mvtnorm (>=1.1.3), MCMCpack (>=1.6.3), rBeta2009 (>=1.0), truncnorm (>=1.0.8), stats (>=4.1.3)
3. Install the package BACT.

```R
devtools::install_github("yinqiaoyan/BACT")
```

## Original datasets information



|   Dataset   | Number of cells | Number of genes | Download link |
| :---------: | :-------------: | :-------------: | :-----------: |
|  STARmap*   |      1207       |      1020       |               |
| MERFISH0.04 |      5488       |       155       |               |
| MERFISH0.09 |      5557       |       155       |               |
| MERFISH0.14 |      5926       |       155       |               |
| MERFISH0.19 |      5803       |       155       |               |
| MERFISH0.24 |      5543       |       155       |               |
|  Slide-seq  |      24847      |      18906      |               |



## Usage example 1: STARmap\*

Read data from .h5ad file

```python
python Read_data.py
```



```R
library(BACT)
library(SingleCellExperiment)
library(aricode)
```



```R
coord = read.csv(paste0(RootPath, "coordinates.csv"))
geneData_raw = read.csv(paste0(RootPath, "gene_count_matrix_raw.csv"))
truth_labels = read.csv(paste0(RootPath, "annotation_file.csv"))$Annotation

spotLoc_char = geneData_raw[, 1]
geneData_raw_noName = as.matrix(geneData_raw[, -1])
rownames(geneData_raw_noName) = spotLoc_char
# Dim: cells * genes =  1207 * 1020
colnames(coord) = c("coord_x", "coord_y")
```

Before implementing BACT, we need to preprocess the raw data by `DataPreprocess` function

```R
tmpx = t(geneData_raw_noName)
sce <- SingleCellExperiment(assays=list(counts=tmpx),
                            colData = coord)
## Preprocess the raw data
# norm.type="logNorm"
sce = DataPreprocess(sce, n.PCs=50, norm.type="logNorm", select.hvg=FALSE)
## Get the processed data after conducting PCA
gene_data_pc = t(reducedDim(sce, "PCA"))
```

Run BACT

```R
# Total execution time is about 15 minutes
# on a MacBook Pro with Intel Core i5 CPU at 2GHz and 16GB of RAM.
res_list = BACT(gene_data_pc = gene_data_pc, coord = coord, platform = "sc",
                num_init = 7, num_nei = 6,
                d1=3, R1_elem=0.5,
                a_eta=0, b_eta=1.5, IGkappa=2, IGtau=10, dpAlpha=1, 
                a_beta=1, tau_beta=1, tau0=0.01, tau1=0.05, M0=50, 
                numOfMCMC=6000, burnIn=3000, 
                Is_beta_zero=FALSE, Is_warm_start=TRUE,
                Is_kmeans_use_mean_sd=TRUE,
                Is_print=TRUE, print_gap=500,
                Is_random_seed=TRUE, random_seed=99)
```

Show results

```R
# Execution time
res_list$exeTime
# Posterior mode of consensus clustering C and marker-gene indicator \gamma
clIds_mode = apply(res_list$clIds_mcmc, 2, getmode)
# 95% credible interval for spatial interaction parameter \beta
quantile(res_list$pottsBeta_mcmc, c(0.025, 0.975))
# Compared with true labels
table(clIds_mode, truth_labels)
cat("ARI value:", ARI(clIds_mode, truth_labels))
```

Visualization for BACT and save the figure

```R
tmpc = clIds_mode
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(1,4,7,2,5,3,11,8,6,10,9,13)
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2
plot_color=c("#ff6466", "#ffb610", "#c599f3", "#52c084", "#7b92ce", "#d2d1d0", 
             "#6b1499", "#138320", "#3185eb", "#9d766e", "#b2c7e5", "#a8dc93")
ppdata = data.frame(x = coord[,1], y = coord[,2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 6) +
  theme(panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text = element_text(size = 35),
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm')) +
  scale_color_manual(values=plot_color,
                     labels = paste0("C", 1:12)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
par(ask=FALSE)
ggsave("./starmap_bact.png", pp, width = 18, height = 10, dpi = 100)
```

Visualization for cell type annotation and save the figure

```R
tmpc = truth_labels
plot_color=c("#ff6466", "#ffb610", "#c599f3", "#52c084", "#7b92ce",
             "#d2d1d0", "#6b1499", "#138320", "#3185eb", "#9d766e",
             "#b2c7e5", "#a8dc93", "#f29d99", "#FFE0C1", "#EFCDF7","#d57fbe")
ppdata = data.frame(x = coord[,1], y = coord[,2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 6) +
  theme(panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text = element_text(size = 35),
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm')) +
  scale_color_manual(values=plot_color) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
par(ask=FALSE)
ggsave("./starmap_cell_annotation.png", pp, width = 20, height = 10, dpi = 100)
```

## Usage example 2: MERFISH0.19



```R
library(BACT)
library(SingleCellExperiment)
library(aricode)
```



```R
coord = read.csv(paste0(RootPath, "coordinates.csv"))
geneData_raw = read.csv(paste0(RootPath, "gene_count_matrix_raw.csv"))
truth_labels = read.csv(paste0(RootPath, "adata_obs.csv"))$cell_class

spotLoc_char = geneData_raw[, 1]
geneData_raw_noName = as.matrix(geneData_raw[, -1])
rownames(geneData_raw_noName) = spotLoc_char
# Dim: cells * genes =  5488 * 155
colnames(coord) = c("coord_x", "coord_y")
```

Since the raw data have been normalized, we only need to take log before conducting PCA.

```R
tmpx = t(geneData_raw_noName)
sce <- SingleCellExperiment(assays=list(counts=tmpx),
                            colData = coord)
## Preprocess the raw data
# norm.type="logOnly"
sce = DataPreprocess(sce, n.PCs=50, norm.type="logOnly", select.hvg=FALSE)
## Get the processed data after conducting PCA
gene_data_pc = t(reducedDim(sce, "PCA"))
```

Run BACT

```R
# Total execution time is about 30 minutes
# on a MacBook Pro with Intel Core i5 CPU at 2GHz and 16GB of RAM.
res_list = BACT(gene_data_pc = gene_data_pc, coord = coord, platform = "sc",
                num_init = 8, num_nei = 6,
                d1=3, R1_elem=0.5,
                a_eta=0, b_eta=1.5, IGkappa=2, IGtau=10, dpAlpha=1, 
                a_beta=1, tau_beta=1, tau0=0.01, tau1=0.05, M0=50, 
                numOfMCMC=4000, burnIn=2000, 
                Is_beta_zero=FALSE, Is_warm_start=TRUE,
                Is_kmeans_use_mean_sd=TRUE,
                Is_print=TRUE, print_gap=500,
                Is_random_seed=TRUE, random_seed=55)
```







## Example Code

The following code shows an example (Simulation I in the paper) that runs the main function "BINRES" in our package.

```R
library(BACT)
library(aricode)
library(ggplot2)
# Import example data
# (1) coord: Spatial coordinates
# (2) gene_data_pc: processed gene expression data (after log-normalization and PCA)
# (3) truth_labels: Cell type annotation of all cells
data(example_data)
# Dimension of spatial coordinates
dim(coord)
# Dimension of gene expression data
dim(gene_data_pc)
# Auxiliary functions
getmode <- function(v) {
  uniqv <- unique(v)
  res <- uniqv[which.max(tabulate(match(v, uniqv)))]
  return(res)
}
# --- run BACT ---
# Total execution time is about 1.7 minutes
# on a MacBook Pro with Intel Core i5 CPU at 2GHz and 16GB of RAM.
res_list = BACT(gene_data_pc = gene_data_pc, coord = coord, platform = "sc",
                num_init = 7, num_nei = 6,
                d1=3, R1_elem=0.5,
                a_eta=0, b_eta=1.5, IGkappa=2, IGtau=10, dpAlpha=1,
                a_beta=1, tau_beta=1, tau0=0.01, tau1=0.05, M0=50,
                numOfMCMC=600, burnIn=300,
                Is_beta_zero=FALSE, Is_warm_start=TRUE,
                Is_kmeans_use_mean_sd=TRUE,
                Is_print=TRUE, print_gap=100,
                Is_random_seed=TRUE, random_seed=99)
# Execution time
res_list$exeTime
# Posterior mode of consensus clustering C and marker-gene indicator \gamma
clIds_mode = apply(res_list$clIds_mcmc, 2, getmode)
# 95% credible interval for spatial interaction parameter \beta
quantile(res_list$pottsBeta_mcmc, c(0.025, 0.975))
# Compared with true labels
table(clIds_mode, truth_labels)
cat("ARI value:", ARI(clIds_mode, truth_labels))
# --- Visualization ---
tmpc = clIds_mode
tmpc2 = tmpc
tmpc_vec = sort(unique(tmpc))
tmpc2_vec = c(1,4,7,2,5,3,11,8,6,10,9,13)
for (ii in 1:length(unique(tmpc))) {
  tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
}
tmpc = tmpc2
plot_color=c("#ff6466", "#ffb610", "#c599f3", "#52c084", "#7b92ce", "#d2d1d0",
             "#6b1499", "#138320", "#3185eb", "#9d766e", "#b2c7e5", "#a8dc93")
ppdata = data.frame(x = coord[,1], y = coord[,2], c = tmpc)
pp = ggplot(data = ppdata, aes(x=x, y=y)) +
  geom_point(aes(color=factor(c)), size = 6) +
  theme(panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text = element_text(size = 35),
        legend.title = element_blank(),
        legend.key.size = unit(2, 'cm')) +
  scale_color_manual(values=plot_color,
                     labels = paste0("C", 1:12)) +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
par(ask=FALSE)
ggsave("./starmap_bact.png", pp, width = 18, height = 10, dpi = 100)
```

or you can simply run

```R
library(BACT)
# Total execution time is about 1.7 minutes
# on a MacBook Pro with Intel Core i5 CPU at 2GHz and 16GB of RAM.
example("BACT")
```



## Remarks

- If you have any questions regarding this package, please contact Yinqiao Yan at [yanyinqiao@ruc.edu.cn](mailto:yanyinqiao@ruc.edu.cn).