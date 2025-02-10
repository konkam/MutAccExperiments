#' Prior posterior plot for the GCM model
#'
#' @param fit_GCM
#'
#' @returns
#' @export
#'
#' @examples
plot_prior_posterior_GCM <- function(fit_GCM) {
  fit_GCM$mcmc %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    (function(df) {
      params <- c("log10_mean", "log10_sigma")
      prior_params <- params %>%
        lapply(FUN = function(x) paste(x, "_prior", sep = "")) %>%
        Reduce(c, .)
      df %>%
        dplyr::select({{ prior_params }}) %>%
        set_names(nm = params) %>%
        mutate(type = "Prior") %>%
        rbind(df %>%
          dplyr::select({{ params }}) %>%
          mutate(type = "Posterior"))
    }) %>%
    tidyr::pivot_longer(cols = -type, names_to = "parameter", values_to = "value") %>%
    ggplot(aes(x = value, colour = type)) +
    theme_bw() +
    facet_wrap(~parameter, scales = "free") +
    geom_density() +
    scale_color_viridis_d(name = "") +
    xlab("Parameter value") +
    ylab("Density")
}

#' Prior posterior plot for the MMRsaturation model
#'
#' @param fit_GCM
#'
#' @returns
#' @export
#'
#' @examples
plot_prior_posterior_MMRsaturation <- function(fit_MMRsaturation) {
  params <- c("log10_mean", "log10_sigma")
  strains <- c("wt", "MMRminus", "MMRproofminus")

  fit_MMRsaturation$mcmc %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    (function(df) {
      prior_params <- params %>%
        lapply(FUN = function(x) paste(x, "_prior", sep = "")) %>%
        Reduce(c, .)

      lapply(strains, function(strain) {
        posterior_params <- params %>%
          lapply(FUN = function(x) paste(x, strain, sep = "_")) %>%
          Reduce(c, .)

        df %>%
          dplyr::select({{ prior_params }}) %>%
          set_names(nm = posterior_params) %>%
          mutate(type = "Prior") %>%
          rbind(df %>%
            dplyr::select({{ posterior_params }}) %>%
            mutate(type = "Posterior")) %>%
          tidyr::pivot_longer(cols = -type, names_to = "parameter", values_to = "value")
      }) %>%
        bind_rows() %>%
        rbind(df %>%
          dplyr::select(theta4_prior) %>%
          set_names(nm = c("theta4")) %>%
          mutate(type = "Prior") %>%
          rbind(df %>%
            dplyr::select(theta4) %>%
            mutate(type = "Posterior")) %>%
          tidyr::pivot_longer(cols = -type, names_to = "parameter", values_to = "value"))
    }) %>%
    ggplot(aes(x = value, colour = type)) +
    theme_bw() +
    facet_wrap(~parameter, scales = "free", ncol = 2) +
    geom_density() +
    scale_color_viridis_d(name = "") +
    xlab("Parameter value") +
    ylab("Density")
}


#' Prior posterior plot
#'
#' @param fit
#'
#' @returns
#' @export
#'
#' @examples
plot_prior_posterior <- function(fit) {
  if (fit$model_type == "GCM_one_strain") {
    fit %>%
      plot_prior_posterior_GCM()
  } else if (fit$model_type == "MMRsaturation") {
    fit %>%
      plot_prior_posterior_MMRsaturation()
  } else {
    stop("Prior posterior plot not yet implemented for this model")
  }
}
