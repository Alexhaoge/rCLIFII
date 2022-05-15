
#' @name lir_bootstrap
#' @title Nonparametric bootstrapping for LIR models
#'observed_individual
#' @param X A list or matrix containing the identities of individuals identified during each sampling period
#' @param tp A set of observed time.
#' @param seed Random seed
#' @details
#' See Efron and Tibshirani (1993) for details on this function.
#'
#' @return The bootstrap samples.
#'
#' @export
#' @rdname Bootstrap
#'
#' @references Efron, B. and Tibshirani, R. (1986). The bootstrap method for
#' standard errors, confidence intervals, and other measures of statistical accuracy.
#' Statistical Science, Vol 1., No. 1, pp 1-35.
#' @references Efron, B. and Tibshirani, R. (1993) An Introduction to the Bootstrap.
#' Chapman and Hall, New York, London.
#'
#'
#' @examples
#'
#' # Example
#' # generate the original data
#' tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
#' N <- 100; n <- 10
#' # X is a list
#' X <- list()
#' for (i in tp){
#' X[[i]] <- sample(1:N, n)
#' }
#' # bootstrap sample data
#' new_sample <- lir_bootstrap(X, tp)
#' # X is a matrix
#' matrix_data <- list_to_matrix(X, tp)
#' new_sample <- lir_bootstrap(X, tp)
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

#' @name lar_bootstrap
#' @title Nonparametric block bootstrapping for LAR models
#'
#' @param X A list or matrix containing the identities of individuals within
#' study area, and the states or status of individuals
#' during each sampling period.
#' @param block_list A block list for a series of observation time. For example,
#' block_list = list(c(1:5), c(51:55), c(101:105), c(501:505), c(601:605)).
#' @param group_id Groups of individuals. If X is a list, please input group_id. If X is a matrix,
#' this parameter can be skipped and takes the default `NULL` value.
#'
#' @param seed Random seed
#'
#' @details
#' See block bootstrap for time series analysis.
#'
#' @return The block bootstrap samples.
#'
#' @export
#' @rdname block_bootstrap
#'
#' @references Politis, D. N. (2003). The impact of bootstrap methods on time series analysis.
#' Statistical science, 18, 219-230.
#'
#'
#' @examples
#'
#' # Example
#' # generate the original data
#' block_list = list(c(1:5), c(51:55), c(101:105), c(501:505), c(601:605))
#' N <- 100; n <- 4; W=10
#' # X is a list
#' tp <- unlist(block_list)
#' tp <- tp - min(tp) + 1
#' tT <- max(tp)
#' len <- length(tp)
#' pop <- seq(N)
#' group_id <- list()
#' X <- list()
#' n <- .data.clean(n, tp)
#' for (i in tp) {
#'   group_id[[i]] <- sample(W, N, replace = TRUE) # Individuals randomly form W groups
#'   g <- sample(W, n[i]) # n groups are randomly selected from W groups
#'   X[[i]] <- pop[which(group_id[[i]] %in% g)] # The members contained in the selected groups
#' }
#' # bootstrap sample data
#' new_sample <- lar_bootstrap(X, block_list, group_id)
#' # X is a matrix
#' matrix_data <- list_to_matrix(X, tp, group_id)
#' new_sample <- lar_bootstrap(matrix_data, block_list)
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

