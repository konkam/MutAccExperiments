#' Check input data format for the GCM model
#'
#' @param input_data
#'
#' @returns an error or the message "All good!"
#' @export
#'
#' @examples
check_input_format_GCM <- function(input_data) {
  column_names <- colnames(input_data)
  if (!("mutation_id" %in% column_names)) {
    stop('The column "mutation_id" is missing from the input data, but it is needed')
  }
  logical_vector <- map(c("m", "n", "t"), function(x) {
    x %in% column_names
  }) %>% unlist()
  if (!all(logical_vector)) {
    stop('Columns "m", "n", "t" are missing from the input data, but they are needed')
  }
  print("All good!")
}

# minimal_input_data %>% check_input_format
