ConvertToWider <- function(mut_acc_experiment_data) {
  tidyr::pivot_wider(mut_acc_experiment_data,
    id_cols = c("mutation_id", "n"),
    names_from = c("strain"),
    values_from = c("m", "t")
  )
}
