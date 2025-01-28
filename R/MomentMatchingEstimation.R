MomentMatchingEstimation_differential_correction <- function(n, m_PolC_MMR, t_PolC_MMR, m_MMR, t_MMR, m_PolC, t_PolC, m_wt, t_wt) {
  hat_mu_PolC_MMR <- m_PolC_MMR / (n * t_PolC_MMR)
  theta1 <- hat_mu_PolC_MMR

  hat_mu_MMR <- m_MMR / (n * t_MMR)
  theta2 <- hat_mu_MMR / hat_mu_PolC_MMR

  hat_mu_PolC <- m_PolC / (n * t_PolC)
  hat_mu_wt <- m_wt / (n * t_wt)

  theta3 <- (hat_mu_wt - hat_mu_PolC) / (hat_mu_MMR - hat_mu_PolC_MMR)

  theta4 <- theta1 * theta3 - hat_mu_PolC

  return(c(theta1, theta2, theta3, theta4))
}

MomentMatchingEstimation_saturation <- function(n, m_PolC_MMR, t_PolC_MMR, m_MMR, t_MMR, m_PolC, t_PolC, m_wt, t_wt) {
  hat_mu_PolC_MMR <- m_PolC_MMR / (n * t_PolC_MMR)
  hat_theta1 <- hat_mu_PolC_MMR

  hat_mu_MMR <- m_MMR / (n * t_MMR)
  hat_theta2 <- hat_mu_MMR / hat_mu_PolC_MMR

  hat_mu_PolC <- m_PolC / (n * t_PolC)
  hat_mu_wt <- m_wt / (n * t_wt)

  hat_theta3 <- hat_mu_wt / hat_mu_MMR

  hat_theta4 <- (hat_mu_PolC/hat_mu_PolC_MMR-1)/(1-hat_mu_wt / hat_mu_MMR)

  return(c(hat_theta1, hat_theta2, hat_theta3, hat_theta4))
}

MomentMatchingEstimation_saturation_all_data = function(mut_acc_experiment_data){
  ConvertToWider(mut_acc_experiment_data) %>%
    (function(tb){
      mapply(MomentMatchingEstimation_saturation,
             tb$n_c,
             tb$m_sc_polC_mutL,
             tb$t_s_polC_mutL,
             tb$`m_sc_MMR-`,
             tb$`t_s_MMR-`,
             tb$m_sc_polC,
             tb$t_s_polC,
             tb$m_sc_WT3610,
             tb$t_s_WT3610)
    })

}
