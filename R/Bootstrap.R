#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Bootstrap sampling
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#'
#' Generate bootstrap samples from original observations by sampling individuals with replacement.
#'
#' @name lir_bootstrap
#' @title Nonparametric bootstrap sampling for animal movement data
#'
#' @param X A list or matrix containing identities of individuals identified per sampling period
#' @param tp A set of observation time
#' @param seed Random seed; the default is NULL
#' @details
#' For more details on this function, please see Efron and Tibshirani (1993).
#'
#' @return The bootstrap samples of animal movement data.
#'
#' @export
#' @rdname lir_bootstrap
#'
#' @references Efron, B., & Tibshirani, R. J. (1986). Bootstrap methods for standard errors,
#' confidence intervals, and other measures of statistical accuracy. Statistical science, 1(1), 54-75.
#'
#' @references Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap. Chapman and Hall Press.
#'
#'
#' @examples
#' # Example
#' # load data
#' data(simulation_A)
#' tp <- simulation_A@tp
#' # if X is a list
#' list_simulation_A <- simulation_A@list_simulation_A
#' bootstrap_sample <- lir_bootstrap(list_simulation_A, tp)
#' # if X is a matrix
#' matrix_simulation_A <- simulation_A@matrix_simulation_A
#' bootstrap_sample <- lir_bootstrap(matrix_simulation_A, tp)
#'

lir_bootstrap <- function(X, tp, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if(is.list(X)){
    matrix_data <- list_to_matrix(X, tp)
    observed_individual <- unlist(X)
  }else{
    matrix_data <- X
    observed_individual <- seq_len(nrow(X))
    colnames(matrix_data) <- tp
    rownames(matrix_data) <- paste('ID', seq_len(nrow(matrix_data)), sep='')
  }
  unique_observed_individual <- observed_individual[!duplicated(observed_individual)]
  len <- length(unique_observed_individual)
  sample_boot <- sample(unique_observed_individual, len, replace = TRUE)
  # update sample
  boot_sample <- matrix_data[match(sample_boot, unique_observed_individual),]
  return(boot_sample)
}

#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Block bootstrap sampling
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#'
#' Observation time is divided into non-overlapping blocks; obtain one bootstrap subsample from
#' each block by sampling the periods with replacement; a complete bootstrap sample is obtained
#' by merging the bootstrap subsamples according to their order in the original time series.
#'
#'
#' @name lar_bootstrap
#' @title Nonparametric block bootstrap sampling for animal social structure data
#'
#' @param X A list or matrix containing identities of individuals and groups affiliated with them per sampling period
#' @param block_list A block list for a series of observation time
#' @param group_id Groups to which individuals belong; if X is a list, please input group_id;
#' if X is a matrix, this parameter can be skipped and and takes the default `NULL` value
#'
#' @param seed Random seed; the default is NULL
#'
#' @details
#' For more details on this function, please see the block bootstrap sampling for time series analysis.
#'
#' @return The bootstrap samples of animal social structure data.
#'
#' @export
#' @rdname lar_bootstrap
#'
#' @references Politis, D. N. (2003). The impact of bootstrap methods on time series analysis.
#' Statistical science, 18, 219-230.
#'
#'
#' @examples
#'
#' # Example
#' # load data
#' data(simulation_D)
#' block_list <- simulation_D@block_list
#' group_id <- simulation_D@group_id
#' # if X is a list
#' list_simulation_D <- simulation_D@list_simulation_D
#' bootstrap_sample <- lar_bootstrap(list_simulation_D, block_list, group_id)
#' # if X is a matrix
#' matrix_simulation_D <- simulation_D@matrix_simulation_D
#' bootstrap_sample <- lar_bootstrap(matrix_simulation_D, block_list)
#'
#'
lar_bootstrap <- function(X, block_list, group_id = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  tp <- unlist(block_list)
  if(is.list(X)){
    if(is.null(group_id)) {
      stop('please input group_id')
    }else{
      matrix_data <- list_to_matrix(X, tp, group_id)
    }
  }else{
    matrix_data <- X
    colnames(matrix_data) <- tp
    rownames(matrix_data) <- paste('ID', seq_len(nrow(matrix_data)), sep='')
  }
  sample_boot <- unlist(lapply(block_list, function(x)sort(sample(x, replace = TRUE))))
  idx <- match(sample_boot, tp)
  # update sample
  boot_sample <- matrix_data[, idx]
  return(boot_sample)
}

