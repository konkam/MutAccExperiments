ConvertToWider = function(mut_acc_experiment_data){
  tidyr::pivot_wider(mut_acc_experiment_data, id_cols = c("context_id", "n_c"),
                     names_from = c("strain"),
                     values_from = c("m_sc", "t_s"))
}
