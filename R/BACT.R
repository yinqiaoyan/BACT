#' BACT: nonparametric Bayesian cell typing for single-cell spatial transcriptomics data
#'
#' The function BACT is the model training function in this paper in which PCA
#'     is used to reduce the dimensionality for the normalized gene expressions
#'     in logarithmic scale. Normal priors are assigned to top principal
#'     components.
#'
#' @param gene_data_pc n.PCs*n preprocessed gene expression matrix. Obtained by normalizing ST raw count matrix, taking logarithm, and conducting PCA.
#' @param coord Coordinates dataframe (2 columns). 1st column: first dimension coordinate. 2nd column: second dimension coordinate.
#' @param platform Spatial sequencing platform. Used to determine neighborhood structure (ST = square, Visium = hex, sc = single-cell resolved ST).
#' @param num_init Initial region number. Default is 5.
#' @param num_nei Number of neighbors. Required if platform is "sc". Default is 6.
#' @param d1 Degree of freedom for the inverse Wishart prior of \eqn{\Lambda_k}. Default is 3.
#' @param R1_elem Diagonal element of matrix R1. Default is 0.5.
#' @param a_eta Mean of the normal prior for \eqn{\eta}. Default is 0.
#' @param b_eta Standard deviation of the normal prior for \eqn{\eta}. Default is 1.5.
#' @param IGkappa Shape parameter of the inverse gamma prior for \eqn{\sigma_g}. Default is 2.
#' @param IGtau Scale parameter of the inverse gamma prior for \eqn{\sigma_g}. Default is 10.
#' @param dpAlpha Hyperparameter of the GEM distribution for the stick-breaking prior of \eqn{\pi_k}. That is, \eqn{\xi_i} are drawn from Be(1, dpAlpha). Default is 1.
#' @param a_beta Mean of the normal distribution before truncation for the spatial interaction parameter \eqn{\beta}. Default is 1.
#' @param tau_beta Standard deviation of the normal distribution before truncation for \eqn{\beta}. Default is 1.
#' @param tau0 Standard deviation of the normal distribution before truncation for the proposal distribution of \eqn{\xi_k^*} when k < M0. Default is 0.01.
#' @param tau1 Standard deviation of the normal distribution before truncation for the proposal distribution of \eqn{\beta}. Default is 0.05.
#' @param M0 A relatively large fixed positive integer. Used to determine proposal distribution form of \eqn{\xi_k^*}. Default is 50.
#' @param numOfMCMC Number of MCMC iterations. Default is 4000.
#' @param burnIn Number of iterations in burn-in. After burnIn the posterior samples are used and saved to estimate the unknown parameters. Default is 2000.
#' @param Is_beta_zero Logical; if TRUE, \eqn{\beta} is fixed at zero. Default is FALSE.
#' @param Is_warm_start Logical; if TRUE, warm start steps by KMeans are used to initialize C. Default is FALSE.
#' @param Is_kmeans_use_mean_sd Logical; if TRUE, results by KMeans are used to initialize mean and standard deviation of each cluster. Required if Is_warm_start is TRUE. Default is FALSE.
#' @param Is_print Logical; if TRUE, iteration time information of each update step are printed. Default is TRUE.
#' @param print_gap Length of iteration interval to print the number of iterations. Default is 10.
#' @param Is_random_seed Logical; if TRUE, a random seed is used for reproducibility. Default is TRUE.
#' @param random_seed Random seed. Required if Is_random_seed is TRUE. Default is 30.
#'
#' @return BINRES returns an R list including the following information.
#' \item{clIds_mcmc}{matrix, the posterior samples of integrative region indicators for each spot or cell. Rows: MCMC samples. Columns: n cells.}
#' \item{eta_k_mcmc}{list, each element contains the posterior sample of \eqn{\eta_k} for all clusters in each MCMC iteration.}
#' \item{sigma_g_mcmc}{matrix, the posterior samples of \eqn{\sigma_g} for each gene. Rows: MCMC samples. Columns: PCs.}
#' \item{pottsBeta_mcmc}{vector, the posterior samples of spatial interaction parameter \eqn{\beta}.}
#' \item{dpXi_mcmc}{list, each element contains the posterior sample of \eqn{\xi_k} in each MCMC iteration.}
#' \item{exeTime}{Total execution time of running the code.}
#'
#' @examples
#' library(BACT)
#' library(aricode)
#' library(ggplot2)
#' # Import example data
#' # (1) coord: Spatial coordinates
#' # (2) gene_data_pc: processed gene expression data (after log-normalization and PCA)
#' # (3) truth_labels: Cell type annotation of all cells
#' data(example_data)
#' # Dimension of spatial coordinates
#' dim(coord)
#' # Dimension of gene expression data
#' dim(gene_data_pc)
#' # Auxiliary functions
#' getmode <- function(v) {
#'   uniqv <- unique(v)
#'   res <- uniqv[which.max(tabulate(match(v, uniqv)))]
#'   return(res)
#' }
#' # --- run BACT ---
#' # Total execution time is about 1.5 minutes
#' # on a MacBook Pro with Intel Core i5 CPU at 2GHz and 16GB of RAM.
#' res_list = BACT(gene_data_pc = gene_data_pc, coord = coord, platform = "sc",
#'                 num_init = 7, num_nei = 6,
#'                 d1=3, R1_elem=0.5,
#'                 a_eta=0, b_eta=1.5, IGkappa=2, IGtau=10, dpAlpha=1,
#'                 a_beta=1, tau_beta=1, tau0=0.01, tau1=0.05, M0=50,
#'                 numOfMCMC=600, burnIn=300,
#'                 Is_beta_zero=FALSE, Is_warm_start=TRUE,
#'                 Is_kmeans_use_mean_sd=TRUE,
#'                 Is_print=TRUE, print_gap=100,
#'                 Is_random_seed=TRUE, random_seed=99)
#' # Execution time
#' res_list$exeTime
#' # Posterior mode of consensus clustering C
#' clIds_mode = apply(res_list$clIds_mcmc, 2, getmode)
#' # Compared with true labels
#' table(clIds_mode, truth_labels)
#' cat("ARI value:", ARI(clIds_mode, truth_labels))
#'
#' ## --- Visualization for BACT --- ##
#' tmpc = clIds_mode
#' tmpc2 = tmpc
#' tmpc_vec = sort(unique(tmpc))
#' tmpc2_vec = c(1,4,7,2,5,3,11,8,6,10,9,13)
#' for (ii in 1:length(unique(tmpc))) {
#'   tmpc2[tmpc == tmpc_vec[ii]] = tmpc2_vec[[ii]]
#' }
#' tmpc = tmpc2
#' plot_color=c("#ff6466", "#ffb610", "#c599f3", "#52c084", "#7b92ce", "#d2d1d0",
#'              "#6b1499", "#138320", "#3185eb", "#9d766e", "#b2c7e5", "#a8dc93")
#' par(ask = FALSE)
#'
#' ppdata = data.frame(x = coord[,1], y = coord[,2], c = tmpc)
#' pp = ggplot(data = ppdata, aes(x=x, y=y)) +
#'   geom_point(aes(color=factor(c)), size = 6) +
#'   theme(panel.background = element_blank(),
#'         axis.title.x=element_blank(),
#'         axis.text.x=element_blank(),
#'         axis.ticks.x=element_blank(),
#'         axis.title.y=element_blank(),
#'         axis.text.y=element_blank(),
#'         axis.ticks.y=element_blank(),
#'         legend.text = element_text(size = 35),
#'         legend.title = element_blank(),
#'         legend.key.size = unit(2, 'cm')) +
#'   scale_color_manual(values=plot_color,
#'                      labels = paste0("C", 1:12)) +
#'   guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
#'
#' ggsave("./starmap_bact.png", pp, width = 18, height = 10, dpi = 100)
#'
#' ## --- Visualization for Annotation --- ##
#' tmpc = truth_labels
#' plot_color=c("#ff6466", "#ffb610", "#c599f3", "#52c084", "#7b92ce",
#'              "#d2d1d0", "#6b1499", "#138320", "#3185eb", "#9d766e",
#'              "#b2c7e5", "#a8dc93", "#f29d99", "#FFE0C1", "#EFCDF7","#d57fbe")
#' ppdata = data.frame(x = coord[,1], y = coord[,2], c = tmpc)
#' pp = ggplot(data = ppdata, aes(x=x, y=y)) +
#'   geom_point(aes(color=factor(c)), size = 6) +
#'   theme(panel.background = element_blank(),
#'         axis.title.x=element_blank(),
#'         axis.text.x=element_blank(),
#'         axis.ticks.x=element_blank(),
#'         axis.title.y=element_blank(),
#'         axis.text.y=element_blank(),
#'         axis.ticks.y=element_blank(),
#'         legend.text = element_text(size = 35),
#'         legend.title = element_blank(),
#'         legend.key.size = unit(2, 'cm')) +
#'   scale_color_manual(values=plot_color) +
#'   guides(color = guide_legend(override.aes = list(size = 12), ncol = 2))
#'
#' ggsave("./starmap_cell_annotation.png", pp, width = 20, height = 10, dpi = 100)
#'
#' @references
#' @export
#' @importFrom stats kmeans rnorm runif dbeta dnorm dist sd
#' @importFrom mvtnorm rmvnorm
#' @importFrom MCMCpack riwish rinvgamma
#' @importFrom truncnorm rtruncnorm dtruncnorm
#' @importFrom rBeta2009 rbeta
BACT <- function(gene_data_pc, coord, platform=c("ST", "Visium", "sc"),
                 num_init=5, num_nei = 6,
                 d1=3, R1_elem=0.5,
                 a_eta=0, b_eta=1.5, IGkappa=2, IGtau=10, dpAlpha=1,
                 a_beta=1, tau_beta=1, tau0=0.01, tau1=0.05, M0=50,
                 numOfMCMC=4000, burnIn=2000,
                 Is_beta_zero=FALSE, Is_warm_start=FALSE,
                 Is_kmeans_use_mean_sd=FALSE,
                 Is_print=TRUE, print_gap=10,
                 Is_random_seed=TRUE, random_seed=30) {

  # For single-cell resolved ST data, find neighbors of each cell
  if (platform == "sc") {
    dists = as.matrix(dist(coord))
    neighbor_indices = matrix(nrow = ncol(gene_data_pc), ncol = num_nei)
    for (ii in 1:ncol(gene_data_pc)) {
      dists_vec = dists[ii, ]
      neighbor_indices[ii, ] = order(dists_vec)[2:(num_nei + 1)]
    }
  }

  # transport the gene data
  gene_data_pc = t(as.matrix(gene_data_pc))

  # location coordinates
  X_loc = as.integer(coord[, 1])
  Y_loc = as.integer(coord[, 2])

  # data shapes
  numOfData = nrow(gene_data_pc)
  G = ncol(gene_data_pc)


  ##########################
  ###  Training process  ###
  ##########################
  cat(paste0("=== Initialization ===\n"))

  if (Is_random_seed) set.seed(random_seed)
  if (Is_warm_start) {
    # warm start using KMeans
    kmres_c <- kmeans(x = gene_data_pc, centers = num_init)
    clIds  = kmres_c[["cluster"]]
    K_max = max(clIds)

    if (Is_kmeans_use_mean_sd) {
      eta_k = kmres_c[["centers"]]
      sigma_g = apply(gene_data_pc, 2, sd)
    }

  } else {
    K_max = num_init
    clIds = sample(1:num_init, numOfData, replace = T)
  }


  # ONLY gene expression data
  if (!Is_kmeans_use_mean_sd) {
    eta_k = matrix(rnorm(K_max * G, a_eta, b_eta), nrow = K_max)
    sigma_g = sqrt(rinvgamma(G, shape = IGkappa, scale = IGtau))
  }


  if (!Is_beta_zero) {
    pottsBeta = rtruncnorm(1, a = 0, mean = a_beta, sd = tau_beta)
  } else {
    pottsBeta = 0
  }


  # initialize dpXi
  dpXi <- rbeta(M0, 1, dpAlpha)



  ### store MCMC samples
  clIds_mcmc = matrix(0, nrow = numOfMCMC - burnIn, ncol = numOfData)
  eta_k_mcmc <- vector(mode = "list", length = numOfMCMC - burnIn)
  sigma_g_mcmc = matrix(0, nrow = numOfMCMC - burnIn, ncol = G)
  pottsBeta_mcmc = array(0, dim = numOfMCMC - burnIn)
  dpXi_mcmc <- vector(mode = "list", length = numOfMCMC - burnIn)



  cat(paste0("=== MCMC Iterations ===\n"))

  ssTime = Sys.time()
  for (mcmc in 0:numOfMCMC) {

    if (!Is_beta_zero) {
      # update dpXi and pottsBeta (double MH)
      # [-] proposal of \xi_k: tnorm (k <= M0); Be(1, alpha) (k > M0)
      # [-] proposal of \beta: tnorm

      tmplen = length(dpXi)
      dpXi_proposal = numeric(tmplen)
      dpXi_proposal[1:M0] = rtruncnorm(M0, a = 0, b = 1, mean = dpXi[1:M0], sd = tau0)
      beta_proposal = rtruncnorm(1, a = 0, mean = pottsBeta, sd = tau1)
      if (tmplen > M0) dpXi_proposal[(M0 + 1):tmplen] = rbeta(tmplen - M0, 1, dpAlpha)

      if ( !(beta_proposal < 0 | sum(dpXi_proposal < 0 | dpXi_proposal > 1) > 0) ) {
        # Step 2: Propose new C* from P(C | beta*, {pi_k*})
        aux_c = clIds
        for (i in 1:numOfData){
          if (platform == "sc") {
            tmpNei = neighbor_indices[i, ]
          } else {
            tmpNei = FindNeighbors_Rcpp(i - 1, X_loc, Y_loc, platform) + 1
          }

          uniq_tmpNei_c = unique(aux_c[tmpNei])
          aux_pi = SticksToPi(dpXi_proposal)

          # compute normalizing constant NC_i
          NCi = 1
          for (m in uniq_tmpNei_c) {
            NCi = NCi + aux_pi[m] * ( exp(beta_proposal * sum(aux_c[tmpNei] == m)) - 1 )
          }

          aux_p = aux_pi / NCi
          for (m in uniq_tmpNei_c) {
            aux_p[m] = aux_pi[m] * exp(beta_proposal * sum(aux_c[tmpNei] == m)) / NCi
          }

          # inverse-cdf
          ui = runif(1)
          sumP = sum(aux_p)

          while (ui > sumP) {
            newTau = rbeta(2, 1, dpAlpha)  # two new sampled tau for dpXi and dpXi_proposal
            dpXi = c(dpXi, newTau[1])
            newPi = newTau[2] * prod(1.0 - dpXi_proposal)
            dpXi_proposal = c(dpXi_proposal, newTau[2])
            aux_p = c(aux_p, newPi / NCi)
            sumP = sumP + newPi / NCi
          }

          csum_aux_p = cumsum(aux_p)
          aux_c[i] = sum(ui > csum_aux_p) + 1

        }

        pi_old = SticksToPi(dpXi)
        pi_new = SticksToPi(dpXi_proposal)


        # Step 3: accept proposal with probability min(1,r)
        # log ratio
        # prior
        logfrac1 = sum(dbeta(dpXi_proposal[1:M0], 1, dpAlpha, log = TRUE)) +
          log(dtruncnorm(beta_proposal, a=0, mean = a_beta, sd = tau_beta))
        logfrac2 = sum(dbeta(dpXi[1:M0], 1, dpAlpha, log = TRUE)) +
          log(dtruncnorm(pottsBeta, a=0, mean = a_beta, sd = tau_beta))

        # C
        if (platform == "sc") {
          logfrac1 = logfrac1 +
            sum(log(pi_new[clIds])) + sum(log(pi_old[aux_c])) +
            ComputePottsDist_sc_Rcpp(beta_proposal, clIds, X_loc, Y_loc, neighbor_indices) + ComputePottsDist_sc_Rcpp(pottsBeta, aux_c, X_loc, Y_loc, neighbor_indices)
          logfrac2 = logfrac2 +
            sum(log(pi_old[clIds])) + sum(log(pi_new[aux_c])) +
            ComputePottsDist_sc_Rcpp(pottsBeta, clIds, X_loc, Y_loc, neighbor_indices) + ComputePottsDist_sc_Rcpp(beta_proposal, aux_c, X_loc, Y_loc, neighbor_indices)
        } else {
          logfrac1 = logfrac1 +
            sum(log(pi_new[clIds])) + sum(log(pi_old[aux_c])) +
            ComputePottsDist_Rcpp(beta_proposal, clIds, X_loc, Y_loc, platform) + ComputePottsDist_Rcpp(pottsBeta, aux_c, X_loc, Y_loc, platform)
          logfrac2 = logfrac2 +
            sum(log(pi_old[clIds])) + sum(log(pi_new[aux_c])) +
            ComputePottsDist_Rcpp(pottsBeta, clIds, X_loc, Y_loc, platform) + ComputePottsDist_Rcpp(beta_proposal, aux_c, X_loc, Y_loc, platform)
        }

        # proposal
        logfrac1 = logfrac1 +
          sum(log(dtruncnorm(dpXi[1:M0], a=0, b=1, mean = dpXi_proposal[1:M0], sd = tau0))) +
          log(dtruncnorm(pottsBeta, a=0, mean = beta_proposal, sd = tau1))
        logfrac2 = logfrac2 +
          sum(log(dtruncnorm(dpXi_proposal[1:M0], a=0, b=1, mean = dpXi[1:M0], sd = tau0))) +
          log(dtruncnorm(beta_proposal, a=0, mean = pottsBeta, sd = tau1))

        ratio = logfrac1 - logfrac2
        ratio = exp(ratio)
        prob = min(1, ratio)

        tmpU = runif(1)
        if (tmpU < prob) {
          pottsBeta = beta_proposal
          dpXi = dpXi_proposal
          # print("accept")
        } else {
          # print("reject")
        }

      }
    } else {
      # pottsBeta is fixed at zero
      # update dpXi (same procedure in Dirichlet process)

      dpXi = dpXi[1:max(clIds)]
      for (k in 1:max(clIds)) {
        Nk_dpXi = sum(clIds == k)
        N_larger_k_dpXi = sum(clIds > k)
        dpXi[k] = rbeta(1, 1 + Nk_dpXi, dpAlpha + N_larger_k_dpXi)
      }
    }



    # update c_lw
    Pi = SticksToPi(dpXi)
    # act_Inds = which(gamma_g == 1)
    for (i in 1:numOfData){
      # update u
      u_i = runif(1, min=0, max=Pi[clIds[i]])

      res         = FindFiniteSet_c_R(u_i, K_max, dpAlpha, a_eta, b_eta, dpXi, eta_k)
      finiteSet_c = res$finiteSet
      dpXi        = res$dpXi
      eta_k       = res$eta_k
      K_max       = res$K_max

      set_size <- length(finiteSet_c)
      q_c <- numeric(set_size)
      if (platform == "sc") {
        vec_nei <- neighbor_indices[i, ]
      } else {
        vec_nei <- FindNeighbors_Rcpp(i - 1, X_loc, Y_loc, platform) + 1
      }

      if (set_size == 0) next

      for (kk in 1:set_size) {
        q_c[kk] = sum(dnorm(gene_data_pc[i, ], mean = eta_k[finiteSet_c[kk], ],
                            sd = sigma_g, log = TRUE)) +
          sum(clIds[vec_nei] == finiteSet_c[kk]) * pottsBeta
      }
      q_c = q_c - max(q_c)
      q_c = exp(q_c) / sum(exp(q_c))

      if (set_size > 1) {
        clIds[i] <- sample(finiteSet_c, 1, prob = q_c)
      } else {
        clIds[i] <- finiteSet_c[1]
      }
    }
    K_max = max(clIds)
    eta_k <- eta_k[1:K_max, ]



    # update eta_k
    for (k in 1:K_max){
      tmpInds = which(clIds == k)
      Nk = length(tmpInds)

      if (Nk != 0) {
        tmpSum = if (Nk > 1) colSums(gene_data_pc[tmpInds, ]) else gene_data_pc[tmpInds, ]
        tmpSigma2 = 1 / (Nk / sigma_g^2 + 1 / b_eta^2)
        tmpMu = tmpSigma2 * (tmpSum / sigma_g^2 + a_eta / b_eta^2)
        eta_k[k, ] = rnorm(G, mean = tmpMu, sd = sqrt(tmpSigma2))
      } else {
        eta_k[k, ] = rnorm(G, mean = a_eta, sd = b_eta)
      }
    }



    # update sigma_g
    for (g in 1:G){
      tmpSumSqu = sum((gene_data_pc[, g] - eta_k[clIds, g])^2) / 2
      sigma_g[g] = sqrt(rinvgamma(1, shape = numOfData / 2 + IGkappa, scale = IGtau + tmpSumSqu)) # sigma_g is std
    }



    # Results --------------------------------------------------------------------
    if (mcmc > burnIn) {
      clIds_mcmc[mcmc - burnIn, ] = clIds
      eta_k_mcmc[[mcmc - burnIn]] = eta_k
      sigma_g_mcmc[mcmc - burnIn, ] = sigma_g
      pottsBeta_mcmc[mcmc - burnIn] = pottsBeta
      dpXi_mcmc[[mcmc - burnIn]] = dpXi
    }

    ## Output
    if (Is_print) {
      if (mcmc <= burnIn) {
        if (mcmc==0) cat(" Burn-in:") else if (mcmc/print_gap == floor(mcmc/print_gap)) cat(paste0("+++", mcmc))
      } else {
        if (mcmc==burnIn+1) cat("\n MCMC sampling:") else if (mcmc > burnIn & mcmc/print_gap == floor(mcmc/print_gap)) cat(paste0("...", mcmc))
      }
    }

  }
  eeTime = Sys.time()
  exeTime = eeTime - ssTime

  cat(paste0("\n=== End Train ===\n"))
  print(exeTime)

  # return results list
  res_list = vector("list")
  res_list[["clIds_mcmc"]]      = clIds_mcmc
  res_list[["eta_k_mcmc"]]      = eta_k_mcmc
  res_list[["sigma_g_mcmc"]]    = sigma_g_mcmc
  res_list[["pottsBeta_mcmc"]]  = pottsBeta_mcmc
  res_list[["dpXi_mcmc"]]       = dpXi_mcmc
  res_list[["exeTime"]]         = exeTime

  return(res_list)
}






