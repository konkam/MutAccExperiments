#' Traceplot for MCMC chains
#'
#' @param fit 
#'
#' @returns
#' @export
#'
#' @examples
traceplot = function(fit){
  if(fit$model_type == "MMRsaturation" | fit$model_type == "GCM_one_strain"){
    fit %>% 
      extract_posterior_samples(type="hyperparameters")  %>% 
      tidyr::pivot_longer(cols = -c(chain_id, iteration), names_to = "parameter", values_to = "value") %>%
      ggplot(aes(x = iteration, y = value, color = factor(chain_id))) +
      theme_bw() +
      facet_wrap(~parameter, scales = "free_y") +
      geom_line() +
      xlab("Iteration") +
      ylab("Parameter value") +
      theme(legend.position = "none") +
      scale_colour_viridis_d()
  }
  else{
    plot(fit)
  }
  
}
