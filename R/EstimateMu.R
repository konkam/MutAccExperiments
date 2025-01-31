# PriorPosteriorMu <- function(jags_fit, context) {
#   mcmc_sample <- jags_fit$mcmc %>%
#     lapply(as_tibble) %>%
#     bind_rows() %>%
#     mutate(mu_prior = 10^log10_mu_prior) %>%
#     select(-log10_mu_prior)
# 
#   parnames <- c("mu_wt", "mu_proofminus", "mu_MMR", "mu_MMRproofminus") %>%
#     paste(., "[", context, "]", sep = "")
# 
#   parnames %>%
#     lapply(function(parname) {
#       mcmc_sample %>%
#         select(parname, mu_prior) %>%
#         gather(type, value) %>%
#         mutate(type = ifelse(type == "mu_prior", yes = "Prior", no = "Posterior")) %>%
#         mutate(strain = parname)
#     }) %>%
#     bind_rows() %>%
#     ggplot(aes(x = value, colour = type)) +
#     theme_bw() +
#     facet_wrap(~strain) +
#     geom_histogram() +
#     scale_x_log10()
# }

Mus_MMRsaturation_model <- "model{

  # Likelihood
  for (i in 1:ncontext){
    m_wt[i] ~ dpois(n[i]*mu_wt[i]*t_wt[i])
    m_proofminus[i] ~ dpois(n[i]*mu_proofminus[i]*t_proofminus[i])
    m_MMRminus[i] ~ dpois(n[i]*mu_MMRminus[i]*t_MMRminus[i])
    m_MMRproofminus[i] ~ dpois(n[i]*mu_MMRproofminus[i]*t_MMRproofminus[i])
    
    loglik_ind[i] <- log(dpois(m_wt[i], n[i]*mu_wt[i]*t_wt[i]) + dpois(m_proofminus[i], n[i]*mu_proofminus[i]*t_proofminus[i]) + dpois(m_MMRminus[i], n[i]*mu_MMRminus[i]*t_MMRminus[i]) + dpois(m_MMRproofminus[i], n[i]*mu_MMRproofminus[i]*t_MMRproofminus[i]))
  }
  
  loglikelihood <- sum(loglik_ind) # Calculate the full log likelihood


  # Reparametrisation

  for (i in 1:ncontext){
    mu_wt[i] = 10^log10_mu_wt[i]
    mu_MMRminus[i] = 10^log10_mu_MMRminus[i]
    mu_MMRproofminus[i] = 10^log10_mu_MMRproofminus[i]

    mu_proofminus[i] = mu_MMRproofminus[i]*(theta4 + (1-theta4)*mu_wt[i]/mu_MMRminus[i])
  }

  # Prior

  for (i in 1:ncontext){
    log10_mu_wt[i] ~ dnorm(log10_mean_wt, tau_wt)
    log10_mu_MMRminus[i] ~ dnorm(log10_mean_MMRminus, tau_MMRminus)
    log10_mu_MMRproofminus[i] ~ dnorm(log10_mean_MMRproofminus, tau_MMRproofminus)
  }

  tau_wt = 1/log10_sigma_wt^2
  tau_MMRminus = 1/log10_sigma_MMRminus^2
  tau_MMRproofminus = 1/log10_sigma_MMRproofminus^2

  log10_mean_wt ~ dnorm(-8, 1/3^2)
  log10_mean_MMRminus ~ dnorm(-8, 1/3^2)
  log10_mean_MMRproofminus ~ dnorm(-8, 1/3^2)
  log10_sigma_wt ~ dt(0, 3, 5)T(0,)
  log10_sigma_MMRminus ~ dt(0, 3, 5)T(0,)
  log10_sigma_MMRproofminus ~ dt(0, 3, 5)T(0,)

  theta4 ~ dbeta(1, 1)

  log10_mean_prior ~ dnorm(-8, 1/3^2)
  log10_sigma_prior ~ dt(0, 3, 5)T(0,)
  mean_prior = 10^log10_mean_prior
  tau_prior = 1/log10_sigma_prior^2
  log10_mu_prior  ~ dnorm(log10_mean_prior, tau_prior)
  theta4_prior ~ dbeta(1, 1)

  # Posterior predictive

  for (i in 1:ncontext){
    m_wt_pred[i] ~ dpois(n[i]*mu_wt[i]*t_wt[i])
    m_proofminus_pred[i] ~ dpois(n[i]*mu_proofminus[i]*t_proofminus[i])
    m_MMRminus_pred[i] ~ dpois(n[i]*mu_MMRminus[i]*t_MMRminus[i])
    m_MMRproofminus_pred[i] ~ dpois(n[i]*mu_MMRproofminus[i]*t_MMRproofminus[i])
  }

}"


#' Title
#'
#' @param mut_acc_experiment_data
#'
#' @return
#' @export
#'
#' @examples
EstimateMusMMRsaturation <- function(mut_acc_experiment_data) {
  wider_data <- ConvertToWider(mut_acc_experiment_data)

  data_jags <- list(
    n = wider_data$n,
    m_MMRproofminus = wider_data$m_MMRproofminus,
    t_MMRproofminus = wider_data$t_MMRproofminus,
    m_MMRminus = wider_data$`m_MMRminus`,
    t_MMRminus = wider_data$`t_MMRminus`,
    m_proofminus = wider_data$m_proofminus,
    t_proofminus = wider_data$t_proofminus,
    m_wt = wider_data$m_wt,
    t_wt = wider_data$t_wt,
    ncontext = length(wider_data$mutation_id)
  )

  model_jags <- rjags::jags.model(
    file = textConnection(Mus_MMRsaturation_model),
    data = data_jags
  )

  res = runjags::autorun.jags(
    model = Mus_MMRsaturation_model,
    monitor = c("mu_wt", "mu_proofminus", "mu_MMRminus", "mu_MMRproofminus", "log10_mean_wt", "log10_mean_MMRminus", "log10_mean_MMRproofminus", "log10_sigma_wt", "log10_sigma_MMRminus", "log10_sigma_MMRproofminus", "theta4", "log10_mean_prior", "log10_sigma_prior", "log10_mu_prior", "theta4_prior", "m_wt_pred", "m_proofminus_pred", "m_MMRminus_pred", "m_MMRproofminus_pred", "loglikelihood"),
    data = data_jags,
    max.time = "1m",
    thin.sample = 4000
  )
  
  res$model_type = "MMRsaturation"
  res$input_data = mut_acc_experiment_data
  return(res)
}


MusGCM_one_strain_model <- "model{

  # Likelihood

  for (i in 1:nmutations){
    m[i] ~ dpois(n[i]*mu[i]*t[i])
    loglik_ind[i] <- log(dpois(m[i], n[i]*mu[i]*t[i]))
  }

  loglikelihood <- sum(loglik_ind) # Calculate the full log likelihood

  # Reparametrisation

  for (i in 1:nmutations){
    mu[i] = 10^log10_mu[i]
  }

  # Prior

  for (i in 1:nmutations){
    log10_mu[i] ~ dnorm(log10_mean, tau)
  }

  tau = 1/log10_sigma^2

  # Prior samples

  log10_mean ~ dnorm(-8, 1/3^2)
  log10_sigma ~ dt(0, 3, 5)T(0,)

  log10_mean_prior ~ dnorm(-8, 1/3^2)
  log10_sigma_prior ~ dt(0, 3, 5)T(0,)
  mean_prior = 10^log10_mean_prior
  tau_prior = 1/log10_sigma^2
  log10_mu_prior  ~ dnorm(log10_mean_prior, tau_prior)

  # Posterior predictive

  for (i in 1:nmutations){
    m_pred[i] ~ dpois(n[i]*mu[i]*t[i])
  }


}"

#' Title
#'
#' @param mut_acc_experiment_data
#'
#' @return
#' @export
#'
#' @examples
#' data(minimal_input_data_onestrain)
#' EstimateMusGCM_onestrain(minimal_input_data_onestrain)
EstimateMusGCM_onestrain <- function(mut_acc_experiment_data) {

  data_jags <- list(
    n = mut_acc_experiment_data$n,
    m = mut_acc_experiment_data$m,
    t = mut_acc_experiment_data$t,
    nmutations = length(mut_acc_experiment_data$mutation_id)
  )

  model_jags <- rjags::jags.model(
    file = textConnection(MusGCM_one_strain_model),
    data = data_jags
  )

  res = runjags::autorun.jags(
    model = MusGCM_one_strain_model,
    monitor = c("m_pred", "log10_mu", "log10_mean", "log10_sigma", "log10_mean_prior", "log10_sigma_prior", "log10_mu_prior", "loglikelihood"),
    data = data_jags,
    max.time = "1m",
    thin.sample = 4000
  )
  res$model_type = "GCM_one_strain"
  res$input_data = mut_acc_experiment_data
  return(res)
}
