## BACT

Current BACT R package: version 1.0 

The R package BACT is built to implement the nonparametric Bayesian method named BACT to perform BAyesian Cell Typing by utilizing gene expression information and spatial coordinates of cells. BACT incorporates a nonparametric Potts prior to induce neighboring cells' spatial dependency, and more importantly it can automatically learn the cell type number directly from the data without prespecification. In this package, a partially collapsed Gibbs sampler is carefully designed for Bayesian posterior inference. BACT can be installed in Windows, Linux, and Mac OS.



## Prerequisites and Installation

1. R version >= 4.1.3.
2. CRAN package: mvtnorm (>=1.1.3), MCMCpack (>=1.6.3), rBeta2009 (>=1.0), truncnorm (>=1.0.8), stats (>=4.1.3)
3. Install the package BACT.

```R
devtools::install_github("yinqiaoyan/BACT")
```



## Datasets information

The data description is given in the following table.

| ST Dataset  | Cell number | Gene number |                        Download links                        |
| :---------: | :---------: | :---------: | :----------------------------------------------------------: |
|  STARmap*   |    1,207    |    1,020    | Raw data: http://sdmbench.drai.cn/tcm/download/?file_path=/mnt/JINGD/data/file/sdmbench/db/STARmap_20180505_BY3_1k.h5ad  <br/>Cell type annotation: https://drive.google.com/drive/folders/1I1nxheWlc2RXSdiv24dex3YRaEh780my?usp=sharing |
| MERFISH0.04 |    5,488    |     155     | http://sdmbench.drai.cn/tcm/download/?file_path=/mnt/JINGD/data/file/sdmbench/db/MERFISH_0.04.h5ad |
| MERFISH0.09 |    5,557    |     155     | http://sdmbench.drai.cn/tcm/download/?file_path=/mnt/JINGD/data/file/sdmbench/db/MERFISH_0.09.h5ad |
| MERFISH0.14 |    5,926    |     155     | http://sdmbench.drai.cn/tcm/download/?file_path=/mnt/JINGD/data/file/sdmbench/db/MERFISH_0.14.h5ad |
| MERFISH0.19 |    5,803    |     155     | http://sdmbench.drai.cn/tcm/download/?file_path=/mnt/JINGD/data/file/sdmbench/db/MERFISH_0.19.h5ad |
| MERFISH0.24 |    5,543    |     155     | http://sdmbench.drai.cn/tcm/download/?file_path=/mnt/JINGD/data/file/sdmbench/db/MERFISH_0.24.h5ad |
|  Slide-seq  |   24,847    |   18,906    | https://singlecell.broadinstitute.org/single_cell/study/SCP354/slide-seq-study#study-download. Download the barcode file "Puck\_180430\_1.tar.gz" |



## Example Code

The following code shows an example (Application to the mouse visual cortex STARmap* data in the manuscript) that runs the main function "BACT" in our package.

Import the required R packages.

```R
library(BACT)
library(aricode)
library(ggplot2)
```

 Read the example data stored in this package. The example data includes:

* coord: spatial coordinates
* gene_data_pc: processed gene expression data (after log-normalization and PCA)
* truth_labels: cell type annotation of cells

```R
data(example_data)
dim(coord)
# Dim: 1207 * 2
dim(gene_data_pc)
# Dim: 50 * 1207
```

Define an auxiliary function that is used to compute the posterior mode from the MCMC samples of cell type indicator $\mathbf{C}$.

```R
getmode <- function(v) {
  uniqv <- unique(v)
  res <- uniqv[which.max(tabulate(match(v, uniqv)))]
  return(res)
}
```

Run BACT function for model training. Total execution time is about 1.5 minutes on a MacBook Pro with Intel Core i5 CPU at 2GHz and 16GB of RAM.

```R
res_list = BACT(gene_data_pc = gene_data_pc, coord = coord, platform = "sc",
                num_init = 7, num_nei = 6,
                a_eta=0, b_eta=1.5, IGkappa=2, IGtau=10, dpAlpha=1,
                a_beta=1, tau_beta=1, tau0=0.01, tau1=0.05, M0=50,
                numOfMCMC=600, burnIn=300,
                Is_beta_zero=FALSE, Is_warm_start=TRUE,
                Is_kmeans_use_mean_sd=TRUE,
                Is_print=TRUE, print_gap=100,
                Is_random_seed=TRUE, random_seed=99)
```

Output the total execution time.

```R
res_list$exeTime
# Time difference of 1.473178 mins
```

Compute Posterior mode of cell type indicator $\mathbf{C}$. Compare the estimated cell type labels with the cell type annotation.

```R
clIds_mode = apply(res_list$clIds_mcmc, 2, getmode)

table(clIds_mode, truth_labels)
cat("ARI value:", ARI(clIds_mode, truth_labels))
# ARI value: 0.6223543
```

The following codes visualize the cell typing performance of BACT, and save the figure.

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
par(ask = FALSE)

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

ggsave("./starmap_bact.png", pp, width = 18, height = 10, dpi = 100)
```

The following codes visualize the underlying cell type annotation, and save the figure.

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

ggsave("./starmap_cell_annotation.png", pp, width = 20, height = 10, dpi = 100)
```

Users can simply run the code `example("BACT")` to carry out this example.

```R
library(BACT)
# Total execution time is about 1.5 minutes
# on a MacBook Pro with Intel Core i5 CPU at 2GHz and 16GB of RAM.
example("BACT")
```



## Remarks

- If you have any questions regarding this package, please contact Yinqiao Yan at [yanyinqiao@ruc.edu.cn](mailto:yanyinqiao@ruc.edu.cn).