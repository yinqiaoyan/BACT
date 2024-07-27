#' @useDynLib BACT, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#'
NULL



#' Data: example_data
#'
#' example_data contains coord and gene_data_pc
#' coord: coordinates of the cells
#' gene_data_pc: gene expression data after PCA
#' @name example_data
#' @docType data
#' @usage data(example_data)
NULL



#' Convert beta random variables to weights pi_k in stick-breaking process
#' @param sticks A beta random variable vector
#' @return A vector of weights pi_k
#' @export
SticksToPi = function(sticks) {
  edges = 1.0 - cumprod(1.0 - sticks)
  return(diff(c(0, edges)))
}



#' Determine the finite set {k: pi_k >= u_i} for each spot i
#' @param u Augmented variable for spot i.
#' @param K_max Maximum in the cell type indicators C.
#' @param dpAlpha Concentration parameter in Dirichlet process.
#' @param a_eta Mean of the normal prior for \eqn{\eta}.
#' @param b_eta Standard deviation of the normal prior for \eqn{\eta}.
#' @param dpXi Beta random variables in the stick-breaking process.
#' @param eta_k Mean value matrix for gene expression data. Rows: cell type clustering number. Columns: genes.
#' @return FindFiniteSet_c_R returns an R list including the following information.
#' \item{finiteSet}{vector, the finite set this function aims to find.}
#' \item{dpXi}{vector, beta random variables in the stick-breaking process.}
#' \item{eta_k}{matrix, updated mean value matrix (maybe add some samples drawn from prior).}
#' \item{K_max}{integer, updated value of K_max (after possibly adding some samples).}
#' @export
FindFiniteSet_c_R <- function(u, K_max, dpAlpha, a_eta, b_eta, dpXi, eta_k) {
  m <- 0
  num_count <- 0
  dpXi_size <- length(dpXi)
  mass <- 1
  finiteSet <- numeric(0)

  while (TRUE) {
    if (m > dpXi_size - 1) {
      dpXi <- c(dpXi, rbeta(1, 1, dpAlpha))
    }

    if (m > K_max - 1) {
      eta_k <- rbind(eta_k, rnorm(1, mean = a_eta, sd = b_eta))
      K_max <- K_max + 1
    }

    weight <- mass * dpXi[m + 1]
    if (weight >= u) {
      finiteSet <- c(finiteSet, m)
      num_count <- num_count + 1
    }

    rest <- mass * (1 - dpXi[m + 1])
    if (rest < u) {
      break
    } else {
      mass <- rest
      m <- m + 1
    }
  }

  return(list(finiteSet = finiteSet + 1,
              dpXi = dpXi,
              eta_k = eta_k,
              K_max = K_max))
}



#' Preprocess the data
#'
#' This code is based on the R function "spatialPreprocess" in BayesSpace
#'
#' @param sce SingleCellExperiment to preprocess
#' @param n.PCs Number of principal components to compute. Default is 50.
#' @param norm.type Type of normalization. logOnly: only take logarithm.
#'   logNorm: compute log-transformed normalized expression values from a count matrix.
#' @param select.hvg Select highly variable genes to run PCA upon. Default is FALSE.
#' @param n.HVGs Number of highly variable genes. Required if select.hvg is TRUE. Default is 2000.
#' @param assay.type Name of assay in \code{sce} containing normalized counts.
#'   Leave as "logcounts" unless you explicitly pre-computed a different
#'   normalization and added it to \code{sce} under another assay.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying which
#'   algorithm should be used to perform the PCA. By default, an exact PCA is
#'   performed, as current spatial datasets are generally small (<10,000 spots).
#'   To perform a faster approximate PCA, please specify
#'   \code{FastAutoParam()} and set a random seed to ensure
#'   reproducibility.
#'
#' @return SingleCellExperiment with PCA
#'
#' @export
#' @importFrom scater logNormCounts runPCA
#' @importFrom scran modelGeneVar getTopHVGs
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom BiocSingular ExactParam
DataPreprocess <- function(sce, n.PCs=50, norm.type=c("logOnly", "logNorm"),
                           select.hvg=FALSE, n.HVGs=2000,
                           assay.type="logcounts", BSPARAM=ExactParam()) {

  ## Run PCA on HVGs, log-normalizing if necessary
  if (norm.type == "logOnly") {
    counts <- counts(sce)
    log_counts <- log(counts + 1)
    assays(sce)$logcounts <- log_counts
  } else if (norm.type == "logNorm") {
    sce <- logNormCounts(sce)
  }

  if (select.hvg) {
    dec <- modelGeneVar(sce, assay.type=assay.type)
    top <- getTopHVGs(dec, n=n.HVGs)
    rowData(sce)[["is.HVG"]] <- (rownames(sce) %in% top)
  } else {
    top <- rownames(sce)
  }

  sce <- runPCA(sce, subset_row=top, ncomponents=n.PCs,
                exprs_values=assay.type, BSPARAM=BSPARAM)

  return(sce)
}


