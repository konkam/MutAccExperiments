#' Simulate some data for one context
#'
#' @param theta1
#' @param theta2
#' @param theta3
#' @param theta4
#' @param n
#' @param t_s
#' @param s
#'
#' @return
#' @export
#'
#' @examples
#' theta = MomentMatchingEstimation(n = 21836, m_PolC_MMR = 4, t_PolC_MMR = 2.30e2, m_MMR = 6, t_MMR = 6.45e3, m_PolC = 2, t_PolC = 1.89e3, m_wt = 4, t_wt = 2.49e5)
#' SimulateData(theta1 = 7.96451e-07, theta2 = 0.05348838, theta3 = 6.330917e-02, theta4 = 1.961348e-09, n = 21836, t_s = c(6453, 1895, 230, 248920), s = c("MMR-", "polC", "polC_mutL", "WT3610"))


SimulateDataDCM = function(theta1 = 7.96451e-07, theta2 = 0.05348838, theta3 = 6.330917e-02, theta4 = 1.961348e-09, n = 21836, t_s = c(6453, 1895, 230, 248920), s = c("MMR-", "polC", "polC_mutL", "WT3610")){

  mu_wt = theta1*theta2*theta3-theta4
  mu_PolC = theta1*theta3-theta4
  mu_MMR = theta1*theta2
  mu_PolC_MMR = theta1

  # tibble::tibble(s = s, n = n, t_s = t_s, mu_s = c(mu_MMR, mu_PolC, mu_PolC_MMR, mu_wt)) %>%
    # dplyr::mutate(m_s = rpois(n = length(s), lambda = mu_s*n*t_s)) %>%
    # return()

  SimulateDataMu(mu_wt = mu_wt, mu_PolC = mu_PolC, mu_MMR = mu_MMR, mu_PolC_MMR = mu_PolC_MMR, n = n, t_s = t_s, s = s)

}

SimulateDataSaturated = function(theta1 = 7.96451e-07, theta2 = 0.05348838, theta3 = 6.330917e-02, theta4 = 1.961348e-09, n = 21836, t_s = c(6453, 1895, 230, 248920), s = c("MMR-", "polC", "polC_mutL", "WT3610")){

  mu_PolC_MMR = theta1
  mu_PolC = theta1*(theta4 + (1-theta4)*theta3)
  mu_MMR = theta1*theta2
  mu_wt = theta1*theta2*theta3

  # tibble::tibble(s = s, n = n, t_s = t_s, mu_s = c(mu_MMR, mu_PolC, mu_PolC_MMR, mu_wt)) %>%
  #   dplyr::mutate(m_s = rpois(n = length(s), lambda = mu_s*n*t_s)) %>%
  #   return()

  SimulateDataMu(mu_wt = mu_wt, mu_PolC = mu_PolC, mu_MMR = mu_MMR, mu_PolC_MMR = mu_PolC_MMR, n = n, t_s = t_s, s = s)

}

SimulateDataMu = function(mu_wt, mu_PolC, mu_MMR, mu_PolC_MMR, n = 21836, t_s = c(6453, 1895, 230, 248920), s = c("MMR-", "polC", "polC_mutL", "WT3610")){

  tibble::tibble(s = s, n = n, t_s = t_s, mu_s = c(mu_MMR, mu_PolC, mu_PolC_MMR, mu_wt)) %>%
    dplyr::mutate(m_s = rpois(n = length(s), lambda = mu_s*n*t_s)) %>%
    return()

}

simulate_data_mu_one_strain_one_mutation_type = function(n, t, mu){
  rpois(n = length(mu), lambda = mu*n*t)
}

#' Simulate data for one strain, all mutation types, given posterior samples on mu
#'
#' @param input_data 
#' @param mu_array 
#'
#' @returns
#' @export
#'
#' @examples
simulate_data_mu_one_strain = function(input_data, mu_array){
  mapply(FUN = function(i, n, t){simulate_data_mu_one_strain_one_mutation_type(n, t, mu_array[,i])},
         input_data$mutation_type, 
         input_data$n, 
         input_data$t)
}
