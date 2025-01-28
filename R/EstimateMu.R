MusGCM_model <- "model{

  # Likelihood
  for (i in 1:ncontext){
    m_wt[i] ~ dpois(n[i]*mu_wt[i]*t_wt[i])
    m_PolC[i] ~ dpois(n[i]*mu_PolC[i]*t_PolC[i])
    m_MMR[i] ~ dpois(n[i]*mu_MMR[i]*t_MMR[i])
    m_PolC_MMR[i] ~ dpois(n[i]*mu_PolC_MMR[i]*t_PolC_MMR[i])
  }


  # Reparametrisation

  for (i in 1:ncontext){
    mu_wt[i] = 10^log10_mu_wt[i]
    mu_PolC[i] = 10^log10_mu_PolC[i]
    mu_MMR[i] = 10^log10_mu_MMR[i]
    mu_PolC_MMR[i] = 10^log10_mu_PolC_MMR[i]
  }

  # Prior

  for (i in 1:ncontext){
    log10_mu_wt[i] ~ dnorm(log10_mean_wt, tau_wt)
    log10_mu_PolC[i] ~ dnorm(log10_mean_PolC, tau_PolC)
    log10_mu_MMR[i] ~ dnorm(log10_mean_MMR, tau_MMR)
    log10_mu_PolC_MMR[i] ~ dnorm(log10_mean_PolC_MMR, tau_PolC_MMR)
  }

  tau_wt = 1/sigma_wt^2
  tau_PolC = 1/sigma_PolC^2
  tau_MMR = 1/sigma_MMR^2
  tau_PolC_MMR = 1/sigma_PolC_MMR^2

  # mean_wt = 10^log10_mean_wt
  # mean_PolC = 10^log10_mean_PolC
  # mean_MMR = 10^log10_mean_MMR
  # mean_PolC_MMR = 10^log10_mean_PolC_MMR

  log10_mean_wt ~ dnorm(-8, 1/3^2)
  log10_mean_PolC ~ dnorm(-8, 1/3^2)
  log10_mean_MMR ~ dnorm(-8, 1/3^2)
  log10_mean_PolC_MMR ~ dnorm(-8, 1/3^2)
  sigma_wt ~ dt(0, 3, 5)T(0,)
  sigma_PolC ~ dt(0, 3, 5)T(0,)
  sigma_MMR ~ dt(0, 3, 5)T(0,)
  sigma_PolC_MMR ~ dt(0, 3, 5)T(0,)

  log10_mean_prior ~ dnorm(-8, 1/3^2)
  sigma_prior ~ dt(0, 3, 5)T(0,)
  mean_prior = 10^log10_mean_prior
  tau_prior = 1/sigma_prior^2
  log10_mu_prior  ~ dnorm(log10_mean_prior, tau_prior)

}"


#' Title
#'
#' @param mut_acc_experiment_data
#'
#' @return
#' @export
#'
#' @examples
EstimateMusGCM <- function(mut_acc_experiment_data) {
  wider_data <- ConvertToWider(mut_acc_experiment_data)

  data_jags <- list(
    n = wider_data$n_c,
    m_PolC_MMR = wider_data$m_sc_polC_mutL,
    t_PolC_MMR = wider_data$t_s_polC_mutL,
    m_MMR = wider_data$`m_sc_MMR-`,
    t_MMR = wider_data$`t_s_MMR-`,
    m_PolC = wider_data$m_sc_polC,
    t_PolC = wider_data$t_s_polC,
    m_wt = wider_data$m_sc_WT3610,
    t_wt = wider_data$t_s_WT3610,
    ncontext = length(wider_data$context_id)
  )

  model_jags <- rjags::jags.model(
    file = textConnection(MusGCM_model),
    # model_jags = rjags::jags.model(file = textConnection(saturated_model_shrinkage_logit_prior_theta4),
    data = data_jags
  )

  runjags::autorun.jags(
    model = MusGCM_model,
    # runjags::autorun.jags(model = saturated_model_shrinkage_logit_prior_theta4,
    monitor = c("mu_wt", "mu_PolC", "mu_MMR", "mu_PolC_MMR", "log10_mean_wt", "log10_mean_PolC", "log10_mean_MMR", "log10_mean_PolC_MMR", "sigma_wt", "sigma_PolC", "sigma_MMR", "sigma_PolC_MMR", "log10_mean_prior", "sigma_prior", "log10_mu_prior"),
    data = data_jags,
    max.time = "1m",
    thin.sample = 4000
  )

  # plot(aa, vars = c( "theta4", "mu_log10_theta1", "sigma_log10_theta1", "alpha2", "beta2", "alpha3", "beta3"))
}

MusGCMNoPooling <- "model{

  # Likelihood
  for (i in 1:ncontext){
    m_wt[i] ~ dpois(n[i]*mu_wt[i]*t_wt[i])
    m_PolC[i] ~ dpois(n[i]*mu_PolC[i]*t_PolC[i])
    m_MMR[i] ~ dpois(n[i]*mu_MMR[i]*t_MMR[i])
    m_PolC_MMR[i] ~ dpois(n[i]*mu_PolC_MMR[i]*t_PolC_MMR[i])
  }


  # Reparametrisation

  for (i in 1:ncontext){
    mu_wt[i] = 10^log10_mu_wt[i]
    mu_PolC[i] = 10^log10_mu_PolC[i]
    mu_MMR[i] = 10^log10_mu_MMR[i]
    mu_PolC_MMR[i] = 10^log10_mu_PolC_MMR[i]
  }

  # Prior

  for (i in 1:ncontext){
    log10_mu_wt[i] ~ dnorm(-8, tau_wt)
    log10_mu_PolC[i] ~ dnorm(-8, tau_PolC)
    log10_mu_MMR[i] ~ dnorm(-8, tau_MMR)
    log10_mu_PolC_MMR[i] ~ dnorm(-8, tau_PolC_MMR)
  }

  tau_wt = 1/10^2
  tau_PolC = 1/10^2
  tau_MMR = 1/10^2
  tau_PolC_MMR = 1/10^2

  log10_mean_prior ~ dnorm(-8, 1/3^2)
  mean_prior = 10^log10_mean_prior
  tau_prior = 1/10^2
  log10_mu_prior  ~ dnorm(log10_mean_prior, tau_prior)

}"

#' Title
#'
#' @param mut_acc_experiment_data
#'
#' @return
#' @export
#'
#' @examples
EstimateMusGCMNoPooling <- function(mut_acc_experiment_data) {
  wider_data <- ConvertToWider(mut_acc_experiment_data)

  data_jags <- list(
    n = wider_data$n_c,
    m_PolC_MMR = wider_data$m_sc_polC_mutL,
    t_PolC_MMR = wider_data$t_s_polC_mutL,
    m_MMR = wider_data$`m_sc_MMR-`,
    t_MMR = wider_data$`t_s_MMR-`,
    m_PolC = wider_data$m_sc_polC,
    t_PolC = wider_data$t_s_polC,
    m_wt = wider_data$m_sc_WT3610,
    t_wt = wider_data$t_s_WT3610,
    ncontext = length(wider_data$context_id)
  )

  model_jags <- rjags::jags.model(
    file = textConnection(MusGCMNoPooling),
    # model_jags = rjags::jags.model(file = textConnection(saturated_model_shrinkage_logit_prior_theta4),
    data = data_jags
  )

  runjags::autorun.jags(
    model = MusGCMNoPooling,
    # runjags::autorun.jags(model = saturated_model_shrinkage_logit_prior_theta4,
    monitor = c("mu_wt", "mu_PolC", "mu_MMR", "mu_PolC_MMR", "log10_mu_prior"),
    data = data_jags,
    max.time = "1m",
    thin.sample = 4000
  )

  # plot(aa, vars = c( "theta4", "mu_log10_theta1", "sigma_log10_theta1", "alpha2", "beta2", "alpha3", "beta3"))
}


DiagnoseMutFit <- function(jags_fit) {
  jags_fit$mcmc %>%
    lapply(as_tibble) %>%
    bind_rows()
}

PriorPosteriorMu <- function(jags_fit, context) {
  mcmc_sample <- jags_fit$mcmc %>%
    lapply(as_tibble) %>%
    bind_rows() %>%
    mutate(mu_prior = 10^log10_mu_prior) %>%
    select(-log10_mu_prior)

  parnames <- c("mu_wt", "mu_PolC", "mu_MMR", "mu_PolC_MMR") %>%
    paste(., "[", context, "]", sep = "")

  parnames %>%
    lapply(function(parname) {
      mcmc_sample %>%
        select(parname, mu_prior) %>%
        gather(type, value) %>%
        mutate(type = ifelse(type == "mu_prior", yes = "Prior", no = "Posterior")) %>%
        mutate(strain = parname)
    }) %>%
    bind_rows() %>%
    ggplot(aes(x = value, colour = type)) +
    theme_bw() +
    facet_wrap(~strain) +
    geom_histogram() +
    scale_x_log10()
}


Mus_saturated_model <- "model{

  # Likelihood
  for (i in 1:ncontext){
    m_wt[i] ~ dpois(n[i]*mu_wt[i]*t_wt[i])
    m_PolC[i] ~ dpois(n[i]*mu_PolC[i]*t_PolC[i])
    m_MMR[i] ~ dpois(n[i]*mu_MMR[i]*t_MMR[i])
    m_PolC_MMR[i] ~ dpois(n[i]*mu_PolC_MMR[i]*t_PolC_MMR[i])
  }


  # Reparametrisation

  for (i in 1:ncontext){
    mu_wt[i] = 10^log10_mu_wt[i]
    # mu_PolC[i] = 10^log10_mu_PolC[i]
    mu_MMR[i] = 10^log10_mu_MMR[i]
    mu_PolC_MMR[i] = 10^log10_mu_PolC_MMR[i]

    mu_PolC[i] = mu_PolC_MMR[i]*(theta4 + (1-theta4)*mu_wt[i]/mu_MMR[i])
  }

  # Prior

  for (i in 1:ncontext){
    log10_mu_wt[i] ~ dnorm(log10_mean_wt, tau_wt)
    # log10_mu_PolC[i] ~ dnorm(log10_mean_PolC, tau_PolC)
    log10_mu_MMR[i] ~ dnorm(log10_mean_MMR, tau_MMR)
    log10_mu_PolC_MMR[i] ~ dnorm(log10_mean_PolC_MMR, tau_PolC_MMR)
  }

  tau_wt = 1/sigma_wt^2
  # tau_PolC = 1/sigma_PolC^2
  tau_MMR = 1/sigma_MMR^2
  tau_PolC_MMR = 1/sigma_PolC_MMR^2

  # mean_wt = 10^log10_mean_wt
  # mean_PolC = 10^log10_mean_PolC
  # mean_MMR = 10^log10_mean_MMR
  # mean_PolC_MMR = 10^log10_mean_PolC_MMR

  log10_mean_wt ~ dnorm(-8, 1/3^2)
  # log10_mean_PolC ~ dnorm(-8, 1/3^2)
  log10_mean_MMR ~ dnorm(-8, 1/3^2)
  log10_mean_PolC_MMR ~ dnorm(-8, 1/3^2)
  sigma_wt ~ dt(0, 3, 5)T(0,)
  # sigma_PolC ~ dt(0, 3, 5)T(0,)
  sigma_MMR ~ dt(0, 3, 5)T(0,)
  sigma_PolC_MMR ~ dt(0, 3, 5)T(0,)

  theta4 ~ dbeta(1, 1)

  log10_mean_prior ~ dnorm(-8, 1/3^2)
  sigma_prior ~ dt(0, 3, 5)T(0,)
  mean_prior = 10^log10_mean_prior
  tau_prior = 1/sigma_prior^2
  log10_mu_prior  ~ dnorm(log10_mean_prior, tau_prior)

}"


#' Title
#'
#' @param mut_acc_experiment_data
#'
#' @return
#' @export
#'
#' @examples
EstimateMusSaturated <- function(mut_acc_experiment_data) {
  wider_data <- ConvertToWider(mut_acc_experiment_data)

  data_jags <- list(
    n = wider_data$n_c,
    m_PolC_MMR = wider_data$m_sc_polC_mutL,
    t_PolC_MMR = wider_data$t_s_polC_mutL,
    m_MMR = wider_data$`m_sc_MMR-`,
    t_MMR = wider_data$`t_s_MMR-`,
    m_PolC = wider_data$m_sc_polC,
    t_PolC = wider_data$t_s_polC,
    m_wt = wider_data$m_sc_WT3610,
    t_wt = wider_data$t_s_WT3610,
    ncontext = length(wider_data$context_id)
  )

  model_jags <- rjags::jags.model(
    file = textConnection(Mus_saturated_model),
    # model_jags = rjags::jags.model(file = textConnection(saturated_model_shrinkage_logit_prior_theta4),
    data = data_jags
  )

  runjags::autorun.jags(
    model = Mus_saturated_model,
    # runjags::autorun.jags(model = saturated_model_shrinkage_logit_prior_theta4,
    monitor = c("mu_wt", "mu_PolC", "mu_MMR", "mu_PolC_MMR", "log10_mean_wt", "log10_mean_MMR", "log10_mean_PolC_MMR", "sigma_wt", "sigma_MMR", "sigma_PolC_MMR", "theta4", "log10_mean_prior", "sigma_prior", "log10_mu_prior"),
    data = data_jags,
    max.time = "1m",
    thin.sample = 4000
  )

  # plot(aa, vars = c( "theta4", "mu_log10_theta1", "sigma_log10_theta1", "alpha2", "beta2", "alpha3", "beta3"))
}


Mus_saturated_model_PolC_from_MMR <- "model{

  # Likelihood
  for (i in 1:ncontext){
    m_wt[i] ~ dpois(n[i]*mu_wt[i]*t_wt[i])
    m_PolC[i] ~ dpois(n[i]*mu_PolC[i]*t_PolC[i])
    m_MMR[i] ~ dpois(n[i]*mu_MMR[i]*t_MMR[i])
    m_PolC_MMR[i] ~ dpois(n[i]*mu_PolC_MMR[i]*t_PolC_MMR[i])
  }


  # Reparametrisation

  for (i in 1:ncontext){
    mu_wt[i] = 10^log10_mu_wt[i]
    # mu_PolC[i] = 10^log10_mu_PolC[i]
    mu_MMR[i] = 10^log10_mu_MMR[i]
    mu_PolC_MMR[i] = 10^log10_mu_PolC_MMR[i]

    mu_PolC[i] = rho*mu_PolC_MMR[i]
  }

  # Prior

  for (i in 1:ncontext){
    log10_mu_wt[i] ~ dnorm(log10_mean_wt, tau_wt)
    # log10_mu_PolC[i] ~ dnorm(log10_mean_PolC, tau_PolC)
    log10_mu_MMR[i] ~ dnorm(log10_mean_MMR, tau_MMR)
    log10_mu_PolC_MMR[i] ~ dnorm(log10_mean_PolC_MMR, tau_PolC_MMR)
  }

  tau_wt = 1/sigma_wt^2
  # tau_PolC = 1/sigma_PolC^2
  tau_MMR = 1/sigma_MMR^2
  tau_PolC_MMR = 1/sigma_PolC_MMR^2

  # mean_wt = 10^log10_mean_wt
  # mean_PolC = 10^log10_mean_PolC
  # mean_MMR = 10^log10_mean_MMR
  # mean_PolC_MMR = 10^log10_mean_PolC_MMR

  log10_mean_wt ~ dnorm(-8, 1/3^2)
  # log10_mean_PolC ~ dnorm(-8, 1/3^2)
  log10_mean_MMR ~ dnorm(-8, 1/3^2)
  log10_mean_PolC_MMR ~ dnorm(-8, 1/3^2)
  sigma_wt ~ dt(0, 3, 5)T(0,)
  # sigma_PolC ~ dt(0, 3, 5)T(0,)
  sigma_MMR ~ dt(0, 3, 5)T(0,)
  sigma_PolC_MMR ~ dt(0, 3, 5)T(0,)

  rho ~ dt(0, 50, 5)T(0,)

  # Prior
  log10_mean_prior ~ dnorm(-8, 1/3^2)
  sigma_prior ~ dt(0, 3, 5)T(0,)
  mean_prior = 10^log10_mean_prior
  tau_prior = 1/sigma_prior^2
  log10_mu_prior  ~ dnorm(log10_mean_prior, tau_prior)
  rho_prior ~ dt(0, 3, 5)T(0,)

}"


#' Title
#'
#' @param mut_acc_experiment_data
#'
#' @return
#' @export
#'
#' @examples
EstimateMusSaturated_PolC_from_MMR <- function(mut_acc_experiment_data) {
  wider_data <- ConvertToWider(mut_acc_experiment_data)

  data_jags <- list(
    n = wider_data$n_c,
    m_PolC_MMR = wider_data$m_sc_polC_mutL,
    t_PolC_MMR = wider_data$t_s_polC_mutL,
    m_MMR = wider_data$`m_sc_MMR-`,
    t_MMR = wider_data$`t_s_MMR-`,
    m_PolC = wider_data$m_sc_polC,
    t_PolC = wider_data$t_s_polC,
    m_wt = wider_data$m_sc_WT3610,
    t_wt = wider_data$t_s_WT3610,
    ncontext = length(wider_data$context_id)
  )

  model_jags <- rjags::jags.model(
    file = textConnection(Mus_saturated_model_PolC_from_MMR),
    data = data_jags
  )

  runjags::autorun.jags(
    model = Mus_saturated_model_PolC_from_MMR,
    monitor = c("mu_wt", "mu_PolC", "mu_MMR", "mu_PolC_MMR", "log10_mean_wt", "log10_mean_MMR", "log10_mean_PolC_MMR", "sigma_wt", "sigma_MMR", "sigma_PolC_MMR", "rho", "rho_prior", "log10_mean_prior", "sigma_prior", "log10_mu_prior"),
    data = data_jags,
    max.time = "1m",
    thin.sample = 4000
  )

}



Mus_saturated_model_no_PolC <- "model{

  # Likelihood
  for (i in 1:ncontext){
    m_wt[i] ~ dpois(n[i]*mu_wt[i]*t_wt[i])
    m_MMR[i] ~ dpois(n[i]*mu_MMR[i]*t_MMR[i])
    m_PolC_MMR[i] ~ dpois(n[i]*mu_PolC_MMR[i]*t_PolC_MMR[i])
  }


  # Reparametrisation

  for (i in 1:ncontext){
    mu_wt[i] = 10^log10_mu_wt[i]
    mu_MMR[i] = 10^log10_mu_MMR[i]
    mu_PolC_MMR[i] = 10^log10_mu_PolC_MMR[i]
  }

  # Prior

  for (i in 1:ncontext){
    log10_mu_wt[i] ~ dnorm(log10_mean_wt, tau_wt)
    log10_mu_MMR[i] ~ dnorm(log10_mean_MMR, tau_MMR)
    log10_mu_PolC_MMR[i] ~ dnorm(log10_mean_PolC_MMR, tau_PolC_MMR)
  }

  tau_wt = 1/sigma_wt^2
  tau_MMR = 1/sigma_MMR^2
  tau_PolC_MMR = 1/sigma_PolC_MMR^2

  log10_mean_wt ~ dnorm(-8, 1/3^2)
  log10_mean_MMR ~ dnorm(-8, 1/3^2)
  log10_mean_PolC_MMR ~ dnorm(-8, 1/3^2)
  sigma_wt ~ dt(0, 3, 5)T(0,)
  sigma_MMR ~ dt(0, 3, 5)T(0,)
  sigma_PolC_MMR ~ dt(0, 3, 5)T(0,)

  log10_mean_prior ~ dnorm(-8, 1/3^2)
  sigma_prior ~ dt(0, 3, 5)T(0,)
  mean_prior = 10^log10_mean_prior
  tau_prior = 1/sigma_prior^2
  log10_mu_prior  ~ dnorm(log10_mean_prior, tau_prior)

}"


#' Title
#'
#' @param mut_acc_experiment_data
#'
#' @return
#' @export
#'
#' @examples
EstimateMusSaturatedNoPolC <- function(mut_acc_experiment_data) {
  wider_data <- ConvertToWider(mut_acc_experiment_data)

  data_jags <- list(
    n = wider_data$n_c,
    m_PolC_MMR = wider_data$m_sc_polC_mutL,
    t_PolC_MMR = wider_data$t_s_polC_mutL,
    m_MMR = wider_data$`m_sc_MMR-`,
    t_MMR = wider_data$`t_s_MMR-`,
    m_wt = wider_data$m_sc_WT3610,
    t_wt = wider_data$t_s_WT3610,
    ncontext = length(wider_data$context_id)
  )

  model_jags <- rjags::jags.model(
    file = textConnection(Mus_saturated_model_no_PolC),
    data = data_jags
  )

  runjags::autorun.jags(
    model = Mus_saturated_model_no_PolC,
    monitor = c("mu_wt", "mu_MMR", "mu_PolC_MMR", "log10_mean_wt", "log10_mean_MMR", "log10_mean_PolC_MMR", "sigma_wt", "sigma_MMR", "log10_mean_prior", "sigma_prior", "log10_mu_prior"),
    data = data_jags,
    max.time = "1m",
    thin.sample = 4000
  )
}

MusGCM_one_strain_model <- "model{

  # Likelihood

  for (i in 1:nmutations){
    m[i] ~ dpois(n[i]*mu[i]*t[i])
  }


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
EstimateMusGCM_onestrain <- function(mut_acc_experiment_data) {

  data_jags <- list(
    n = mut_acc_experiment_data$n,
    m = mut_acc_experiment_data$m,
    t = mut_acc_experiment_data$t,
    nmutations = length(mut_acc_experiment_data$mutation_type)
  )

  model_jags <- rjags::jags.model(
    file = textConnection(MusGCM_one_strain_model),
    data = data_jags
  )

  res = runjags::autorun.jags(
    model = MusGCM_one_strain_model,
    monitor = c("m_pred", "log10_mu", "log10_mean", "log10_sigma", "log10_mean_prior", "log10_sigma_prior", "log10_mu_prior"),
    data = data_jags,
    max.time = "1m",
    thin.sample = 4000
  )
  res$model_type = "GCM_one_strain"
  return(res)
}
