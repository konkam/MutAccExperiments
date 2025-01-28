Theta4FromLambda = function(lambda, posterior_sample_gamma, mut_acc_data){

  n_c = mut_acc_data %>%
    select(n_c, context_id) %>%
    unique()

  phi_sum_gamma_c_times_nc = posterior_sample_gamma %>%
    select(starts_with("theta1")) %>%
    rowid_to_column(var = "iteration") %>%
    gather(context_id, gamma_c, -iteration) %>%
    mutate(context_id = context_id %>%
             gsub("theta1\\[", "", .) %>%
             gsub("\\]", "", .) %>%
             as.numeric()) %>%
    left_join(n_c) %>%
    mutate(gamma_c_times_n_c = gamma_c*n_c) %>%
    group_by(iteration) %>%
    summarise(sum_gamma_c_n_c = sum(gamma_c_times_n_c)) %>%
    mutate(phi = Vectorize(FUN = PhiGammanFromLambda)(lambda, sum_gamma_c_n_c))


  return(phi_sum_gamma_c_times_nc)

}

PhiGammanFromLambda = function(lambda, gamma_times_n){
  (gamma_times_n * qGammaN(lambda = lambda-1, gamma_times_n = gamma_times_n) - lambda*qGammaN(lambda = lambda, gamma_times_n = gamma_times_n))/gamma_times_n
}

# qGammaN = function(lambda, posterior_sample_gamma){
#   1 - mean(ppois(q = lambda, lambda = posterior_sample_gamma))
# }
qGammaN = function(lambda, gamma_times_n){
  1 - ppois(q = lambda, lambda = gamma_times_n)
}

