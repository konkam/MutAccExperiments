ResultSummaryToCsv <- function(res, nm = "") {
  res %>%
    summary() %>%
    as.data.frame() %>%
    cbind(tibble(parameter_name = rownames(.)), .) %>%
    write_csv(nm)
}

SaveChainToCsv <- function(jags_fit, nm = "mcmc_chain.csv") {
  jags_fit %>%
    extract_posterior_samples() %>%
    write_csv(nm)
}

#' Extract posterior samples from a jags fit object for further analysis
#'
#' @param jags_fit
#'
#' @returns The posterior samples as a tibble
#' @export
#'
#' @examples
#' data(minimal_input_data_onestrain)
#' jags_fit <- EstimateMusGCM_onestrain(minimal_input_data_onestrain)
#' extract_posterior_samples(jags_fit)
extract_posterior_samples <- function(jags_fit, type = "all") {
  if (type == "all") {
    extract_posterior_samples_all(jags_fit)
  } else if (type == "prior") {
    extract_posterior_samples_all(jags_fit) %>%
      dplyr::select(contains("prior"), iteration, chain_id)
  } else if (type == "mu") {
    extract_posterior_samples_all(jags_fit) %>%
      dplyr::select(contains("mu"), iteration, chain_id) %>%
      dplyr::select(-contains("prior"), -contains("mean"))
  } else if (type == "hyperparameters" | type == "hyper") {
    extract_posterior_samples_all(jags_fit) %>%
      dplyr::select(contains("sigma"), contains("mean"), contains("theta"), loglikelihood, iteration, chain_id) %>%
      dplyr::select(-contains("prior"))
  } else if (type == "predictive" | type == "pred") {
    extract_posterior_samples_all(jags_fit) %>%
      dplyr::select(contains("pred"), iteration, chain_id)
  } else {
    stop("Invalid type")
  }
}

extract_posterior_samples_all <- function(jags_fit) {
  jags_fit$mcmc %>%
    # lapply(as_tibble) %>%
    # (function(ll){mapply(function(xx, id){mutate(xx, chain_id = id)}, ll, seq_along(ll), SIMPLIFY = F)}) %>%
    purrr::imap(\(x, y) as_tibble(x) %>%
      mutate(chain_id = y) %>%
      rowid_to_column(var = "iteration")) %>%
    bind_rows()
}
