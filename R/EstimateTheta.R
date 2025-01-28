EstimateTheta = function(n, m_PolC_MMR, t_PolC_MMR, m_MMR, t_MMR, m_PolC, t_PolC, m_wt, t_wt){
  # for debugging
  # n = 21836; m_PolC_MMR = 4; t_PolC_MMR = 2.30e2; m_MMR = 6; t_MMR = 6.45e3; m_PolC = 2; t_PolC = 1.89e3; m_wt = 4; t_wt = 2.49e5

  model_string <- "model{

  # Likelihood
  m_wt ~ dpois(lambda_wt_0)
  m_PolC ~ dpois(n*mu_PolC*t_PolC)
  m_MMR ~ dpois(n*mu_MMR*t_MMR)
  m_PolC_MMR ~ dpois(n*mu_PolC_MMR*t_PolC_MMR)

  # Reparametrisation
  lambda_wt_0 = n*mu_wt0*t_wt
  lambda_wt = n*mu_wt*t_wt
  mu_wt0 = theta1_0*theta2_0*theta3-theta4
  mu_wt = theta1*theta2*theta3-theta4 # This may become negative, it may cause some problems
  mu_PolC = theta1*theta3-theta4
  mu_MMR = theta1*theta2
  mu_PolC_MMR = theta1

  theta1_0 = 7.96451e-07
  theta2_0 = 0.05348838
  theta3 = 6.330917e-02
  theta4 = 1.961348e-09

  theta1 = 10^log10theta1
  #
  # sign_theta4 = 2*sign_theta4_indicator - 1
  # norm_theta4 = 10^log10norm_theta4
  # theta4 = sign_theta4*norm_theta4
  #
  # # Prior
  log10theta1 ~ dunif(-7, -6)
  theta2 ~ dbeta(0.5, 1)
  # theta3 ~ dbeta(0.5, 1)
  # log10norm_theta4 ~ dunif(-9, -8)
  # sign_theta4_indicator ~ dbern(0.99)

}"

  data_jags = list(n = n,
              m_PolC_MMR = m_PolC_MMR,
              t_PolC_MMR = t_PolC_MMR,
              m_MMR = m_MMR,
              t_MMR = t_MMR,
              m_PolC = m_PolC,
              t_PolC = t_PolC,
              m_wt = m_wt,
              t_wt = t_wt)

  model_jags = rjags::jags.model(file = textConnection(model_string),
                                 data = data_jags)

  runjags::autorun.jags(model = model_string, monitor = c("theta1", "theta2", "lambda_wt"), data = data_jags)

}


model_string_param2 <- "model{

  # Likelihood
  m_wt ~ dpois(n*mu_wt*t_wt)
  m_PolC ~ dpois(n*mu_PolC*t_PolC)
  m_MMR ~ dpois(n*mu_MMR*t_MMR)
  m_PolC_MMR ~ dpois(n*mu_PolC_MMR*t_PolC_MMR)

  # Reparametrisation
  mu_wt = theta1*theta2*theta3*(1-tilde_theta4)
  mu_PolC = theta1*theta3*(1-theta2*tilde_theta4)
  mu_MMR = theta1*theta2
  mu_PolC_MMR = theta1

  theta1 = 10^log10theta1
  sign_theta4 = 2*sign_theta4_indicator - 1
  tilde_theta4 = sign_theta4*norm_tilde_theta4
  theta4 = tilde_theta4*theta1*theta3*theta3

  # # Prior
  log10theta1 ~ dunif(-7, -6)
  theta2 ~ dbeta(0.5, 1)
  theta3 ~ dbeta(0.5, 1)
  norm_tilde_theta4 ~ dbeta(1,1)
  sign_theta4_indicator ~ dbern(0.99)

}"

saturated_model_1context <- "model{

  # Likelihood
  m_wt ~ dpois(n*mu_wt*t_wt)
  m_PolC ~ dpois(n*mu_PolC*t_PolC)
  m_MMR ~ dpois(n*mu_MMR*t_MMR)
  m_PolC_MMR ~ dpois(n*mu_PolC_MMR*t_PolC_MMR)

  # Reparametrisation
  mu_wt = theta1*theta2*theta3
  mu_PolC = theta1*(1+theta4*(1-theta3))
  mu_MMR = theta1*theta2
  mu_PolC_MMR = theta1

  theta1 = 10^log10theta1

  # # Prior
  log10theta1 ~ dunif(-7, -6)
  theta2 ~ dbeta(0.5, 1)
  theta3 ~ dbeta(0.5, 1)
  theta4 ~ dbeta(0.5, 1)

}"

saturated_model_no_shrinkage <- "model{

  # Likelihood
  for (i in 1:ncontext){
    m_wt[i] ~ dpois(n[i]*mu_wt[i]*t_wt[i])
    m_PolC[i] ~ dpois(n[i]*mu_PolC[i]*t_PolC[i])
    m_MMR[i] ~ dpois(n[i]*mu_MMR[i]*t_MMR[i])
    m_PolC_MMR[i] ~ dpois(n[i]*mu_PolC_MMR[i]*t_PolC_MMR[i])
  }


  # Reparametrisation
  for (i in 1:ncontext){
    mu_wt[i] = theta1[i]*theta2[i]*theta3[i]
    mu_PolC[i] = theta1[i]*(1+theta4*(1-theta3[i]))
    mu_MMR[i] = theta1[i]*theta2[i]
    mu_PolC_MMR[i] = theta1[i]

    theta1[i] = 10^log10theta1[i]
  }


  # # Prior
  for (i in 1:ncontext){
    log10theta1[i] ~ dunif(-7, -6)
    theta2[i] ~ dbeta(0.5, 1)
    theta3[i] ~ dbeta(0.5, 1)
  }

  theta4 ~ dbeta(0.5, 1)

}"

saturated_model_shrinkage <- "model{

  # Likelihood
  for (i in 1:ncontext){
    m_wt[i] ~ dpois(n[i]*mu_wt[i]*t_wt[i])
    m_PolC[i] ~ dpois(n[i]*mu_PolC[i]*t_PolC[i])
    m_MMR[i] ~ dpois(n[i]*mu_MMR[i]*t_MMR[i])
    m_PolC_MMR[i] ~ dpois(n[i]*mu_PolC_MMR[i]*t_PolC_MMR[i])
  }


  # Reparametrisation
  for (i in 1:ncontext){
    mu_wt[i] = theta1[i]*theta2[i]*theta3[i]
    mu_PolC[i] = theta1[i]*(theta4 + (1-theta4)*theta3[i])
    mu_MMR[i] = theta1[i]*theta2[i]
    mu_PolC_MMR[i] = theta1[i]

    theta1[i] = 10^log10theta1[i]
  }


  # # Prior
  for (i in 1:ncontext){
    log10theta1[i] ~ dnorm(mu_log10_theta1, inv_sigma_sq_log10_theta1)
    theta2[i] ~ dbeta(alpha2, beta2)
    theta3[i] ~ dbeta(alpha3, beta3)
  }

  theta4 ~ dbeta(0.5, 1)
  mu_log10_theta1 ~ dunif(-8, -6)

  inv_sigma_sq_log10_theta1 = 1/sigma_log10_theta1^2
  sigma_log10_theta1 ~ dt(0, 0.5, 2)T(0,)

  alpha2  ~ dunif(0, 100)
  alpha3  ~ dunif(0, 100)

  beta2 = 10^log10beta2
  log10beta2  ~ dnorm(0, 4)
  beta3 = 10^log10beta3
  log10beta3  ~ dnorm(0, 4)
}"

saturated_model_shrinkage_no_PolC <- "model{

  # Likelihood
  for (i in 1:ncontext){
    m_wt[i] ~ dpois(n[i]*mu_wt[i]*t_wt[i])
    m_MMR[i] ~ dpois(n[i]*mu_MMR[i]*t_MMR[i])
    m_PolC_MMR[i] ~ dpois(n[i]*mu_PolC_MMR[i]*t_PolC_MMR[i])
  }


  # Reparametrisation
  for (i in 1:ncontext){
    mu_wt[i] = theta1[i]*theta2[i]*theta3[i]
    mu_MMR[i] = theta1[i]*theta2[i]
    mu_PolC_MMR[i] = theta1[i]

    theta1[i] = 10^log10theta1[i]
  }


  # # Prior
  for (i in 1:ncontext){
    log10theta1[i] ~ dnorm(mu_log10_theta1, inv_sigma_sq_log10_theta1)
    theta2[i] ~ dbeta(alpha2, beta2)
    theta3[i] ~ dbeta(alpha3, beta3)
  }

  mu_log10_theta1 ~ dunif(-8, -6)

  inv_sigma_sq_log10_theta1 = 1/sigma_log10_theta1^2
  sigma_log10_theta1 ~ dt(0, 0.5, 2)T(0,)

  alpha2  ~ dunif(0, 100)
  alpha3  ~ dunif(0, 100)

  beta2 = 10^log10beta2
  log10beta2  ~ dnorm(0, 4)
  beta3 = 10^log10beta3
  log10beta3  ~ dnorm(0, 4)
}"

saturated_model_shrinkage_logit_prior_theta4 <- "model{

  # Likelihood
  for (i in 1:ncontext){
    m_wt[i] ~ dpois(n[i]*mu_wt[i]*t_wt[i])
    m_PolC[i] ~ dpois(n[i]*mu_PolC[i]*t_PolC[i])
    m_MMR[i] ~ dpois(n[i]*mu_MMR[i]*t_MMR[i])
    m_PolC_MMR[i] ~ dpois(n[i]*mu_PolC_MMR[i]*t_PolC_MMR[i])
  }


  # Reparametrisation
  for (i in 1:ncontext){
    mu_wt[i] = theta1[i]*theta2[i]*theta3[i]
    # mu_PolC[i] = theta1[i]*(1+theta4*(1-theta3[i]))
    mu_PolC[i] = theta1[i]*(theta4 + (1-theta4)*theta3[i])
    mu_MMR[i] = theta1[i]*theta2[i]
    mu_PolC_MMR[i] = theta1[i]

    theta1[i] = 10^log10theta1[i]
  }


  # # Prior
  for (i in 1:ncontext){
    log10theta1[i] ~ dnorm(mu_log10_theta1, inv_sigma_sq_log10_theta1)
    theta2[i] ~ dbeta(alpha2, beta2)
    theta3[i] ~ dbeta(alpha3, beta3)
  }

  minus_log_theta4 ~ dt(0, 3, 2)T(0,)
  theta4 = exp(-1*minus_log_theta4)
  mu_log10_theta1 ~ dunif(-8, -6)

  inv_sigma_sq_log10_theta1 = 1/sigma_log10_theta1^2
  sigma_log10_theta1 ~ dt(0, 0.5, 2)T(0,)

  alpha2  ~ dunif(0, 100)
  alpha3  ~ dunif(0, 100)

  beta2 = 10^log10beta2
  log10beta2  ~ dnorm(0, 4)
  beta3 = 10^log10beta3
  log10beta3  ~ dnorm(0, 4)
}"

EstimateTheta_param2 = function(n, m_PolC_MMR, t_PolC_MMR, m_MMR, t_MMR, m_PolC, t_PolC, m_wt, t_wt){
  # for debugging
  # n = 21836; m_PolC_MMR = 4; t_PolC_MMR = 2.30e2; m_MMR = 6; t_MMR = 6.45e3; m_PolC = 2; t_PolC = 1.89e3; m_wt = 4; t_wt = 2.49e5

  data_jags = list(n = n,
                   m_PolC_MMR = m_PolC_MMR,
                   t_PolC_MMR = t_PolC_MMR,
                   m_MMR = m_MMR,
                   t_MMR = t_MMR,
                   m_PolC = m_PolC,
                   t_PolC = t_PolC,
                   m_wt = m_wt,
                   t_wt = t_wt)

  model_jags = rjags::jags.model(file = textConnection(model_string_param2),
                                 data = data_jags)

  runjags::autorun.jags(model = model_string_param2, monitor = c("theta1", "theta2", "theta3", "tilde_theta4", "theta4"), data = data_jags)

}

EstimateTheta_param_saturated_one_context = function(n, m_PolC_MMR, t_PolC_MMR, m_MMR, t_MMR, m_PolC, t_PolC, m_wt, t_wt){
  # for debugging
  # n = 21836; m_PolC_MMR = 4; t_PolC_MMR = 2.30e2; m_MMR = 6; t_MMR = 6.45e3; m_PolC = 2; t_PolC = 1.89e3; m_wt = 4; t_wt = 2.49e5

  data_jags = list(n = n,
                   m_PolC_MMR = m_PolC_MMR,
                   t_PolC_MMR = t_PolC_MMR,
                   m_MMR = m_MMR,
                   t_MMR = t_MMR,
                   m_PolC = m_PolC,
                   t_PolC = t_PolC,
                   m_wt = m_wt,
                   t_wt = t_wt)

  model_jags = rjags::jags.model(file = textConnection(saturated_model_1context),
                                 data = data_jags)

  runjags::autorun.jags(model = saturated_model_1context, monitor = c("theta1", "theta2", "theta3", "theta4"), data = data_jags)

}

EstimateTheta_param_saturated_no_shrinkage = function(mut_acc_experiment_data){

  wider_data = mut_acc_experiment_data %>%
    tidyr::pivot_wider(id_cols = c("context_id", "n_c"),
                names_from = c("strain"),
                values_from = c("m_sc", "t_s"))

  data_jags = list(n = wider_data$n_c,
                   m_PolC_MMR = wider_data$m_sc_polC_mutL,
                   t_PolC_MMR = wider_data$t_s_polC_mutL,
                   m_MMR = wider_data$`m_sc_MMR-`,
                   t_MMR = wider_data$`t_s_MMR-`,
                   m_PolC = wider_data$m_sc_polC,
                   t_PolC = wider_data$t_s_polC,
                   m_wt = wider_data$m_sc_WT3610,
                   t_wt = wider_data$t_s_WT3610,
                   ncontext = length(wider_data$context_id))

  model_jags = rjags::jags.model(file = textConnection(saturated_model_no_shrinkage),
                                 data = data_jags)

  runjags::autorun.jags(model = saturated_model_no_shrinkage, monitor = c("theta1", "theta2", "theta3", "theta4"), data = data_jags)

}

EstimateTheta_param_saturated_shrinkage = function(mut_acc_experiment_data){

  wider_data =  ConvertToWider(mut_acc_experiment_data)

  data_jags = list(n = wider_data$n_c,
                   m_PolC_MMR = wider_data$m_sc_polC_mutL,
                   t_PolC_MMR = wider_data$t_s_polC_mutL,
                   m_MMR = wider_data$`m_sc_MMR-`,
                   t_MMR = wider_data$`t_s_MMR-`,
                   m_PolC = wider_data$m_sc_polC,
                   t_PolC = wider_data$t_s_polC,
                   m_wt = wider_data$m_sc_WT3610,
                   t_wt = wider_data$t_s_WT3610,
                   ncontext = length(wider_data$context_id))

  model_jags = rjags::jags.model(file = textConnection(saturated_model_shrinkage),
  # model_jags = rjags::jags.model(file = textConnection(saturated_model_shrinkage_logit_prior_theta4),
                                 data = data_jags)

  runjags::autorun.jags(model = saturated_model_shrinkage,
  # runjags::autorun.jags(model = saturated_model_shrinkage_logit_prior_theta4,
                                              monitor = c("theta1", "theta2", "theta3", "theta4", "mu_wt", "mu_PolC", "mu_MMR", "mu_PolC_MMR","mu_log10_theta1", "sigma_log10_theta1", "alpha2", "beta2", "alpha3", "beta3"), data = data_jags)

  # plot(aa, vars = c( "theta4", "mu_log10_theta1", "sigma_log10_theta1", "alpha2", "beta2", "alpha3", "beta3"))
}

EstimateThetaSaturatedShrinkageNoPolC = function(mut_acc_experiment_data){

  wider_data =  ConvertToWider(mut_acc_experiment_data)

  data_jags = list(n = wider_data$n_c,
                   m_PolC_MMR = wider_data$m_sc_polC_mutL,
                   t_PolC_MMR = wider_data$t_s_polC_mutL,
                   m_MMR = wider_data$`m_sc_MMR-`,
                   t_MMR = wider_data$`t_s_MMR-`,
                   m_wt = wider_data$m_sc_WT3610,
                   t_wt = wider_data$t_s_WT3610,
                   ncontext = length(wider_data$context_id))

  model_jags = rjags::jags.model(file = textConnection(saturated_model_shrinkage_no_PolC),
                                 data = data_jags)

  runjags::autorun.jags(model = saturated_model_shrinkage_no_PolC,
                        monitor = c("theta1", "theta2", "theta3", "mu_wt", "mu_MMR", "mu_PolC_MMR","mu_log10_theta1", "sigma_log10_theta1", "alpha2", "beta2", "alpha3", "beta3"), data = data_jags)

  # plot(aa, vars = c( "theta4", "mu_log10_theta1", "sigma_log10_theta1", "alpha2", "beta2", "alpha3", "beta3"))
}



PriorPredictiveDistribution = function(n, m_PolC_MMR, t_PolC_MMR, m_MMR, t_MMR, m_PolC, t_PolC, m_wt, t_wt){
  data_jags = list(n = n,
                   t_PolC_MMR = t_PolC_MMR,
                   t_MMR = t_MMR,
                   t_PolC = t_PolC,
                   t_wt = t_wt)

  model_jags = rjags::jags.model(file = textConnection(model_string_param2),
                                 data = data_jags)

  runjags::autorun.jags(model = model_string_param2, monitor = c("theta1", "theta2", "theta3", "tilde_theta4", "theta4", "m_wt", "m_PolC", "m_MMR", "m_PolC_MMR"), data = data_jags)
}
