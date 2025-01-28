GCMReparametrisationMuToGammaY <- function(mu_wt, mu_MMR, mu_PolC, mu_PolC_MMR) {
  gamma_1 <- mu_MMR
  gamma_2 <- mu_PolC_MMR - mu_MMR
  y_1 <- mu_wt / mu_MMR
  y_2 <- (mu_PolC - mu_wt) / (mu_PolC_MMR - mu_MMR)
  return(c(gamma_1, gamma_2, y_1, y_2) %>% setNames(c("gamma_1", "gamma_2", "y_1", "y_2")))
}

SaturatedReparametrisationMuToTheta <- function(mu_wt, mu_MMR, mu_PolC_MMR) {
  theta_1 <- mu_PolC_MMR
  theta_2 <- mu_MMR / mu_PolC_MMR
  theta_3 <- mu_wt / mu_MMR
  return(c(theta_1, theta_2, theta_3) %>% setNames(c("theta_1", "theta_2", "theta_3")))
}
