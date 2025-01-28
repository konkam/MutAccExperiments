PosteriorPredictiveSaturated <- function(mut_acc_data, posterior_sample) {
  wider_dat <- ConvertToWider(mut_acc_data)
  n_samples <- nrow(posterior_sample)

  predictive_distribution_summary <- wider_dat$context_id %>%
    lapply(FUN = function(context_id) {
      SimulateSaturatedOneContext(context_id, posterior_sample, wider_dat) %>%
        PredictiveDistrSummary()
    }) %>%
    bind_rows()

  # SimulateDataSaturated()
  mut_acc_data %>%
    left_join(predictive_distribution_summary)
}

SimulateSaturatedOneContextOneIt <- function(context_id, posterior_sample, sample_id, wider_dat) {
  theta1 <- posterior_sample[[paste("theta1[", context_id, "]", sep = "")]][sample_id]
  theta2 <- posterior_sample[[paste("theta2[", context_id, "]", sep = "")]][sample_id]
  theta3 <- posterior_sample[[paste("theta3[", context_id, "]", sep = "")]][sample_id]
  theta4 <- posterior_sample[["theta4"]][sample_id]

  SimulateDataSaturated(
    theta1 = theta1, theta2 = theta2, theta3 = theta3, theta4 = theta4,
    n = wider_dat$n_c[context_id],
    t_s = c(
      wider_dat$`t_s_MMR-`[context_id],
      wider_dat$t_s_polC[context_id],
      wider_dat$t_s_polC_mutL[context_id],
      wider_dat$t_s_WT3610[context_id]
    ),
    s = c("MMR-", "polC", "polC_mutL", "WT3610")
  ) %>%
    (function(df) {
      df$m_s %>%
        setNames(paste("m_sc_pred_", paste(df$s), sep = "")) %>%
        t() %>%
        as_tibble()
    }) %>%
    mutate(context_id = context_id)
}

SimulateMuOneContextOneIt <- function(context_id, posterior_sample_mu, sample_id, wider_dat) {
  mu_wt <- posterior_sample_mu[[paste("mu_wt[", context_id, "]", sep = "")]][sample_id]
  mu_PolC <- posterior_sample_mu[[paste("mu_PolC[", context_id, "]", sep = "")]][sample_id]
  mu_MMR <- posterior_sample_mu[[paste("mu_MMR[", context_id, "]", sep = "")]][sample_id]
  mu_PolC_MMR <- posterior_sample_mu[[paste("mu_PolC_MMR[", context_id, "]", sep = "")]][sample_id]

  # SimulateDataSaturated(
  #   theta1 = theta1, theta2 = theta2, theta3 = theta3, theta4 = theta4,
  #   n = wider_dat$n_c[context_id],
  #   t_s = c(
  #     wider_dat$`t_s_MMR-`[context_id],
  #     wider_dat$t_s_polC[context_id],
  #     wider_dat$t_s_polC_mutL[context_id],
  #     wider_dat$t_s_WT3610[context_id]
  #   ),
  #   s = c("MMR-", "polC", "polC_mutL", "WT3610")
  # )
  SimulateDataMu(mu_wt, mu_PolC, mu_MMR, mu_PolC_MMR, n = wider_dat$n_c[context_id], t_s = c(
    wider_dat$`t_s_MMR-`[context_id],
    wider_dat$t_s_polC[context_id],
    wider_dat$t_s_polC_mutL[context_id],
    wider_dat$t_s_WT3610[context_id]
  ),
  s = c("MMR-", "polC", "polC_mutL", "WT3610")) %>%
    (function(df) {
      df$m_s %>%
        setNames(paste("m_sc_pred_", paste(df$s), sep = "")) %>%
        t() %>%
        as_tibble()
    }) %>%
    mutate(context_id = context_id)
}

SimulateSaturatedOneContext <- function(context_id, posterior_sample, wider_dat) {
  1:nrow(posterior_sample) %>%
    lapply(function(it) {
      SimulateSaturatedOneContextOneIt(context_id, posterior_sample, it, wider_dat)
    }) %>%
    bind_rows()
}


PredictiveDistrSummary <- function(predicted_sample) {
  column_names <- colnames(predicted_sample) %>% setdiff("context_id")
  lapply(column_names, function(colnam) {
    predicted_sample[[colnam]] %>%
      (function(vv) {
        c(mean(vv), quantile(vv, probs = c(0.5, 0.025, 0.975, 0.25, 0.75))) %>%
          setNames(paste("m_sc_pred_", c("mean", "median", "infCI", "supCI", "infquart", "supquart"), sep = "")) %>%
          t() %>%
          as_tibble() %>%
          mutate(strain = colnam %>% gsub("m_sc_pred_", "", .))
      })
  }) %>%
    bind_rows() %>%
    mutate(context_id = predicted_sample$context_id %>% unique())
}

predictive_distr_summary_one_strain <- function(predicted_sample_one_strain){
  predicted_sample_one_strain %>%
    apply(MARGIN = 2, FUN = function(vv) {
      c(mean(vv), quantile(vv, probs = c(0.5, 0.025, 0.975, 0.25, 0.75))) %>%
        setNames(paste("m_pred_", c("mean", "median", "infCI", "supCI", "infquart", "supquart"), sep = ""))
    }) %>% 
    t %>% 
    as_tibble() %>% 
    rowid_to_column(var = "mutation_type")
}

# plot_predictive_distr_summary_one_strain <- function(predictive_distr_summary_one_strain){
#   predictive_distr_summary_one_strain %>%
#     ggplot(aes(x = mutation_type)) +
#     geom_segment(aes(xend = mutation_type, y = m_pred_infCI, yend = m_pred_supCI)) +
#     geom_segment(aes(xend = mutation_type, y = m_pred_infquart, yend = m_pred_supquart), alpha = 0.9, colour = "grey", linewidth = 2) +
#     geom_point(aes(y = m_pred_mean), colour = "red") +
#     xlab("Mutation type") +
#     ylab("Number of mutations")
# }

SimulateMuOneContext <- function(context_id, posterior_sample, wider_dat) {
  1:nrow(posterior_sample) %>%
    lapply(function(it) {
      SimulateMuOneContextOneIt(context_id, posterior_sample, it, wider_dat)
    }) %>%
    bind_rows()
}

PlotPosteriorPredictiveCheck <- function(data_with_predictive_distribution_summary) {
  data_with_predictive_distribution_summary %>%
    ggplot(aes(x = context_id)) +
    theme_bw() +
    facet_wrap(~strain, ncol = 1, scales = "free_y") +
    geom_segment(aes(xend = context_id, y = m_sc_pred_infCI, yend = m_sc_pred_supCI)) +
    geom_segment(aes(xend = context_id, y = m_sc_pred_infquart, yend = m_sc_pred_supquart), alpha = 0.9, colour = "grey", linewidth = 2) +
    geom_point(
      data = data_with_predictive_distribution_summary %>%
        gather(type, n_mutations, c(m_sc, m_sc_pred_mean)),
      aes(y = n_mutations, colour = type)
    ) +
    # geom_point(aes(y = m_sc_pred_mean), colour = "red") +
    xlab("Context Id") +
    ylab("Number of mutations") +
    scale_colour_manual(name = "", labels = c("Observed", "Predicted"), values = c("red", "blue"))
}


SimulateMuSaturatedOneContextOneIt <- function(context_id, posterior_sample, sample_id, wider_dat) {
  mu_wt <- posterior_sample[[paste("mu_wt[", context_id, "]", sep = "")]][sample_id]
  mu_PolC_MMR <- posterior_sample[[paste("mu_PolC_MMR[", context_id, "]", sep = "")]][sample_id]
  mu_MMR <- posterior_sample[[paste("mu_MMR[", context_id, "]", sep = "")]][sample_id]
  theta4 <- posterior_sample[["theta4"]][sample_id]
  mu_PolC <- mu_PolC_MMR * (theta4 + (1 - theta4) * mu_wt / mu_MMR) ## mu_PolC is actually in the posterior distribution, so we could do without this line and replace by mu_PolC <- posterior_sample[[paste("mu_PolC[", context_id, "]", sep = "")]][sample_id]

  rpois(n = 4, lambda = c(mu_MMR, mu_PolC, mu_PolC_MMR, mu_wt) * wider_dat$n_c[context_id] * c(
    wider_dat$`t_s_MMR-`[context_id],
    wider_dat$t_s_polC[context_id],
    wider_dat$t_s_polC_mutL[context_id],
    wider_dat$t_s_WT3610[context_id]
  )) %>%
    tibble(
      m_s = .,
      s = c("MMR-", "polC", "polC_mutL", "WT3610")
    ) %>%
    (function(df) {
      df$m_s %>%
        setNames(paste("m_sc_pred_", paste(df$s), sep = "")) %>%
        t() %>%
        as_tibble()
    }) %>%
    mutate(context_id = context_id)
}

SimulateMuSaturatedOneContext <- function(context_id, posterior_sample, wider_dat) {
  1:nrow(posterior_sample) %>%
    lapply(function(it) {
      SimulateMuSaturatedOneContextOneIt(context_id, posterior_sample, it, wider_dat)
    }) %>%
    bind_rows()
}

PosteriorPredictiveMuSaturated <- function(mut_acc_data, posterior_sample) {
  wider_dat <- ConvertToWider(mut_acc_data)
  n_samples <- nrow(posterior_sample)

  predictive_distribution_summary <- wider_dat$context_id %>%
    lapply(FUN = function(context_id) {
      SimulateMuSaturatedOneContext(context_id, posterior_sample, wider_dat) %>%
        PredictiveDistrSummary()
    }) %>%
    bind_rows()

  # SimulateDataSaturated()
  mut_acc_data %>%
    left_join(predictive_distribution_summary)
}

#' Title
#'
#' @param mut_acc_data
#' @param posterior_sample
#'
#' @return
#' @export
#'
#' @examples
PosteriorPredictiveMu <- function(mut_acc_data, posterior_sample) {
  wider_dat <- ConvertToWider(mut_acc_data)
  n_samples <- nrow(posterior_sample)

  predictive_distribution_summary <- wider_dat$context_id %>%
    lapply(FUN = function(context_id) {
      SimulateMuOneContext(context_id, posterior_sample, wider_dat) %>%
        PredictiveDistrSummary()
    }) %>%
    bind_rows()

  # SimulateDataSaturated()
  mut_acc_data %>%
    left_join(predictive_distribution_summary)
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
posterior_predictive_one_strain <- function(input_data, fit_mu) {
  
  #Replace this by a mu extraction function
  mu_array = fit_mu$mcmc %>% 
    lapply(as_tibble) %>% 
    bind_rows() %>%
    select(starts_with("log10_mu[")) %>% 
    as.matrix() %>% 
    (function(x)10^x)
  
  predictive_distribution_summary = simulate_data_mu_one_strain(input_data, mu_array) %>% 
    predictive_distr_summary_one_strain
    
  input_data %>%
    left_join(predictive_distribution_summary)
}
