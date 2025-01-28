ResultSummaryToCsv = function(res, nm = ""){
  res %>%
    summary %>%
    as.data.frame %>%
    cbind(tibble(parameter_name = rownames(.)), .) %>%
    write_csv(nm)
}

SaveChainToCsv = function(jags_fit, nm = "mcmc_chain.csv"){
  jags_fit %>%
    extract_posterior_samples %>% 
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
extract_posterior_samples = function(jags_fit){
  jags_fit$mcmc %>%
    # lapply(as_tibble) %>%
    # (function(ll){mapply(function(xx, id){mutate(xx, chain_id = id)}, ll, seq_along(ll), SIMPLIFY = F)}) %>%
    purrr::imap(\(x,y) as_tibble(x) %>% mutate(chain_id = y) %>% rowid_to_column(var = "iteration")) %>% 
    bind_rows()
}