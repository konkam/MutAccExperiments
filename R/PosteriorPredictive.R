predictive_distr_summary <- function(predicted_sample) {
  column_names <- colnames(predicted_sample)
  lapply(column_names, function(colnam) {
    predicted_sample[[colnam]] %>%
      (function(vv) {
        c(mean(vv), quantile(vv, probs = c(0.5, 0.025, 0.975, 0.25, 0.75))) %>%
          setNames(paste("m_pred_", c("mean", "median", "infCI", "supCI", "infquart", "supquart"), sep = "")) %>%
          t() %>%
          as_tibble() %>%
          mutate(colnam = colnam)
      })
  }) %>%
    bind_rows() %>%
    mutate(
      strain = colnam %>% str_split_i(pattern = "_", i = 2),
      mutation_id = colnam %>% str_split_i(pattern = "\\[", i = 2) %>% gsub("]", "", .) %>% as.numeric()
    )
}


predictive_distr_summary_one_strain <- function(predicted_sample_one_strain) {
  predicted_sample_one_strain %>%
    apply(MARGIN = 2, FUN = function(vv) {
      c(mean(vv), quantile(vv, probs = c(0.5, 0.025, 0.975, 0.25, 0.75))) %>%
        setNames(paste("m_pred_", c("mean", "median", "infCI", "supCI", "infquart", "supquart"), sep = ""))
    }) %>%
    t() %>%
    as_tibble() %>%
    rowid_to_column(var = "mutation_id")
}



#' Title
#'
#' @param input_data
#' @param fit_mu
#'
#' @return
#' @export
#'
#' @examples
posterior_predictive_one_strain <- function(fit_mu_one_strain) {
  fit_mu_one_strain %>%
    extract_posterior_samples(type = "predictive") %>%
    select(-iteration, -chain_id) %>%
    predictive_distr_summary_one_strain() %>%
    left_join(fit_mu_one_strain$input_data)

  # # Replace this by a mu extraction function
  # mu_array <- fit_mu$mcmc %>%
  #   lapply(as_tibble) %>%
  #   bind_rows() %>%
  #   select(starts_with("log10_mu[")) %>%
  #   as.matrix() %>%
  #   (function(x) 10^x)
  #
  # predictive_distribution_summary <- simulate_data_mu_one_strain(fit_mu$input_data, mu_array) %>%
  #   predictive_distr_summary_one_strain()
  #
  # fit_mu$input_data %>%
  #   left_join(predictive_distribution_summary)
}

#' Title
#'
#' @param fit_mu
#'
#' @returns
#' @export
#'
#' @examples
posterior_predictive <- function(fit_mu) {
  fit_mu %>%
    extract_posterior_samples(type = "predictive") %>%
    select(-iteration, -chain_id) %>%
    predictive_distr_summary() %>%
    left_join(fit_mu$input_data)
}

#' Title
#'
#' @param fit_mu
#'
#' @returns
#' @export
#'
#' @examples
plot_posterior_predictive_onestrain <- function(fit_mu) {
  to_plot <- fit_mu %>%
    posterior_predictive_one_strain()

  pl <- ggplot(to_plot, aes(x = mutation_id)) +
    geom_segment(aes(y = m_pred_infCI / (n * t), yend = m_pred_supCI / (n * t), xend = mutation_id, color = "95% CI"), width = 0.3, size = 9, alpha = 0.25) +
    geom_segment(aes(y = m_pred_infquart / (n * t), yend = m_pred_supquart / (n * t), xend = mutation_id, color = "Quartiles"), width = 0.3, size = 9) +
    geom_point(aes(y = m / (n * t)), size = 3, fill = NA) +
    geom_point(aes(y = m / (n * t), color = "Empirical rate estimate"), size = 2.) +
    theme_minimal() +
    labs(
      title = "Posterior Predictive check",
      x = "Mutation Type",
      y = "Mutation rate",
      color = "Legend",
      shape = "Legend"
    ) +
    scale_color_manual(values = c("Empirical rate estimate" = "white", "95% CI" = "orange", "Quartiles" = "orange"))

  if ("mutation_label" %in% names(to_plot)) {
    custom_labels <- to_plot %>%
      select(mutation_id, mutation_label) %>%
      distinct() %>%
      arrange(mutation_id) %>%
      pull(mutation_label)
    pl +
      scale_x_continuous(labels = custom_labels, breaks = seq_along(custom_labels)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    pl
  }
}


plot_posterior_predictive <- function(fit_mu) {
  to_plot <- fit_mu %>%
    posterior_predictive()

  pl <- to_plot %>%
    mutate(strain = gsub("minus", "-", strain)) %>%
    ggplot(aes(x = mutation_id)) +
    facet_wrap(~strain, scales = "free_y", ncol = 1) +
    geom_segment(aes(y = m_pred_infCI / (n * t), yend = m_pred_supCI / (n * t), xend = mutation_id, color = "95% CI"), linewidth = 9, alpha = 0.25) +
    geom_segment(aes(y = m_pred_infquart / (n * t), yend = m_pred_supquart / (n * t), xend = mutation_id, color = "Quartiles"), linewidth = 9) +
    geom_point(aes(y = m / (n * t)), size = 3, fill = NA) +
    geom_point(aes(y = m / (n * t), color = "Empirical rate estimate"), size = 2.) +
    theme_minimal() +
    labs(
      title = "Posterior Predictive check",
      x = "Mutation Type",
      y = "Mutation rate",
      color = "Legend",
      shape = "Legend"
    ) +
    scale_color_manual(values = c("Empirical rate estimate" = "white", "95% CI" = "orange", "Quartiles" = "orange"))

  if ("mutation_label" %in% names(to_plot)) {
    custom_labels <- to_plot %>%
      select(mutation_id, mutation_label) %>%
      distinct() %>%
      arrange(mutation_id) %>%
      pull(mutation_label)
    pl +
      scale_x_continuous(labels = custom_labels, breaks = seq_along(custom_labels)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    pl
  }
}
