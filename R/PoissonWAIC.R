PoissonLogPostDensityOneParamSet <- function(mut_acc_data, params, strain_order = c("WT3610", "polC", "MMR-", "polC_mutL")) {
  mut_acc_data %>%
    arrange(factor(strain, levels = strain_order), context_id) %>%
    (function(df) {
      dpois(x = df$m_sc, lambda = params * df$n_c * df$t_s, log = T)
    })
}

PoissonLogPostDensityOneIt <- function(mut_acc_data, mcmc_samples_mu, it, strain_order = c("WT3610", "polC", "MMR-", "polC_mutL")) {
  # mcmc_samples_mu = mcmc_samples %>% select(starts_with("mu"))
  PoissonLogPostDensityOneParamSet(mut_acc_data = mut_acc_data, mcmc_samples_mu[it, ] %>% as.numeric(), strain_order = strain_order)
}

PoissonLogPostDensity <- function(mut_acc_data, mcmc_samples_mu, strain_order = c("WT3610", "polC", "MMR-", "polC_mutL")) {
  # mcmc_samples_mu = mcmc_samples %>% select(starts_with("mu"))
  res <- matrix(NA, nrow = nrow(mcmc_samples_mu), ncol = nrow(mut_acc_data))
  for (it in 1:nrow(mcmc_samples_mu)) {
    res[it, ] <- PoissonLogPostDensityOneParamSet(mut_acc_data = mut_acc_data, mcmc_samples_mu[it, ] %>% as.numeric(), strain_order = strain_order)
  }
  return(res)
}

pWAIC2 <- function(log_post_density) {
  matrixStats::colVars(log_post_density) %>% sum()
}

lppd <- function(log_post_density) {
  log_post_density %>%
    exp() %>%
    matrixStats::colMeans2() %>%
    log() %>%
    sum()
}

#' Title
#'
#' @param mut_acc_data
#' @param mcmc_samples_mu
#' @param strain_order
#'
#' @return
#' @export
#'
#' @examples
PoissonWAIC_from_mcmc_samples_mu <- function(mut_acc_data, mcmc_samples_mu, strain_order = c("WT3610", "polC", "MMR-", "polC_mutL")) {
  log_post_density <- PoissonLogPostDensity(mut_acc_data = mut_acc_data, mcmc_samples_mu = mcmc_samples_mu, strain_order = strain_order)
  return(lppd(log_post_density) - pWAIC2(log_post_density))
}

CreateConvertStrainName <- function() {
  x = c("WT3610", "polC", "MMR-", "polC_mutL")
  strain_conversion_dict <-  setNames(object = c(x, "wt", "PolC", "MMR", "PolC_MMR"), nm = c("wt", "PolC", "MMR", "PolC_MMR", x))

  ConvertStrainName <- function(nms) {
    sapply(nms, function(nm) {
      strain_conversion_dict[nm]
    })
  }

  return(ConvertStrainName)
}

ConvertStrainName <- CreateConvertStrainName()

#' Title
#'
#' @param mut_acc_data
#' @param mcmc_samples_mu
#' @param strain_order
#'
#' @return
#' @export
#'
#' @examples
PoissonWAIC <- function(mut_acc_data, jags_fit) {
  strain_order <- jags_fit$monitor %>%
    grep("mu_", ., value = T) %>%
    (function(x) {
      x[!grepl(pattern = "prior", x = x)]
    }) %>%
    gsub("mu_", "", .) %>%
    ConvertStrainName()

  mcmc_samples_mu <- jags_fit$mcmc %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    select(starts_with("mu"))

  PoissonWAIC_from_mcmc_samples_mu(mut_acc_data = mut_acc_data, mcmc_samples_mu = mcmc_samples_mu, strain_order = strain_order)
}
# lppd_logsumexp = function(LogPostDensity){
#   LogPostDensity %>%
#     apply(2, matrixStats::logSumExp) %>%
#     (function(mat))
#     sum()
# }
