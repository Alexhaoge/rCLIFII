
#' @title A nonparametric estimator for lagged identification rate
#'
#' @param X A list or matrix containing the identities of individuals identified during each sampling period
#' @param n A vector or a positive integer, representing the number of individuals identified in each sampling period.
#' It indicates the same number of individuals identified in all sampling periods if a positive integer.
#' @param tp A set of observation time
#'
#' @return A list with the following elements:
#'
#' \item{tau}{a set of time lags}
#' \item{R_tau}{A nonparametric estimator of \eqn{\hat{R}(\tau)} was given by Whitehead (2007)}
#' \item{R_m}{\eqn{\sum_{i,j} \{m_{t_i,t_j}|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|} }
#' \item{R_n}{\eqn{\sum_{i,j} \{n_i*n_j|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|}}
#' \item{mij}{number of individuals identified at both time \eqn{t_i} and \eqn{t_j}}
#' \item{nij}{\eqn{n_i}*\eqn{n_j}}
#' \item{tauij}{\eqn{\tau_{ij}=|t_i-t_j|}}
#'
#'
#' @details
#' The lagged identification rate \eqn{R(\tau)} is the probability that an individual
#' in the study area is identified now and is reidentified after a time lag of \eqn{\tau}.
#' A nonparametric estimator of \eqn{R(\tau)} was given by Whitehead (2007):
#' \deqn{\hat{R}(\tau)=\frac{\sum_{i,j} \{m_{t_i,t_j}|\tau_{ij}=\tau\}}{\sum_{i,j}\{n_i*n_j|\tau_{ij}=\tau\}},}
#' where where \eqn{m_{t_i,t_j}} is the number of individuals identified at both
#' time \eqn{t_i} and \eqn{t_j}, and \eqn{n_{t_i}} represents the number of individuals
#' identified at time \eqn{t_i}.
#'
#'
#' @export
#' @rdname lir_nonparametric_estimation
#'
#' @examples
#' # Example
#' tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
#' N <- 100; n <- 10
#' X <- list()
#' for (i in tp){
#' X[[i]] <- sample(1:N, n)
#' }
#'
#' res <- lir_nonparametric_estimation(X, n, tp)
#' res$R_tau; res$tau
#'
#'
lir_nonparametric_estimation <- function(X, n, tp) {
  n <- .data.clean(n, tp)
  tp <- tp-min(tp)+1
  if (is.list(X)){
    data <- list_to_matrix(X, tp)
  }else{
    data <- X
  }
  len <- length(tp)
  tT <- max(tp)
  R_m <- rep(0, tT)
  R_n <- rep(0, tT)
  mij <- c()
  nij <- c()
  tauij <- c()
  k <- 1
  for (i in 1:(len-1)){
    for (j in (i+1):len){
      mij[k] <- sum(data[,j]*data[,i]==1)
      nij[k] <- n[tp[i]]*n[tp[j]]
      tauij[k] <- tp[j] - tp[i]
      R_m[tauij[k]] <- R_m[tauij[k]] + mij[k]
      R_n[tauij[k]] <- n[tp[i]]*n[tp[j]] + R_n[tauij[k]]
      k <- k + 1
    }
  }
  R_tauij <- R_m[tauij]/R_n[tauij]
  tau <- tauij[!duplicated(tauij)]
  R_tau <- R_m[tau]/R_n[tau]
  R_data <- list(R_tau=R_tau, tau=tau, R_m=R_m, R_n=R_n,
                 mij=mij, nij=nij, tauij=tauij)
  return(R_data)
}

#' @title A nonparametric estimator for lagged association rate
#'
#' @param X A list or matrix containing the identities of individuals within
#' study area, and the states or status of individuals during each sampling period.
#' @param tp A set of observation time.
#' @param group_id Groups of individuals. If X is a list, please input group_id. If X is a matrix,
#' this parameter can be skipped and takes the default `NULL` value.
#'
#' @return A list with the following elements:
#' \item{tau}{a set of time lags}
#' \item{gtau}{A nonparametric estimator of \eqn{\hat{g}(\tau)} was given by Whitehead (2007)}
#' \item{g_m}{\eqn{\sum_{i,j} \{A_{t_i,t_j}|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|} }
#' \item{g_n}{\eqn{\sum_{i,j} \{A_i|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|}}
#' \item{Aij}{number of associations of pair of individuals at both time \eqn{t_i} and \eqn{t_j}}
#' \item{Ai}{number of associations of pair of individuals at time \eqn{t_i}}
#' \item{tauij}{\eqn{\tau_{ij}=|t_i-t_j|}}
#'
#'
#' @details
#' The lagged association rate \eqn{g(\tau)} is the probability that if two individuals
#' are associated now, they will still be associated after a time lag of  \eqn{\tau}.
#' A nonparametric estimator of \eqn{g(\tau)} was given by Whitehead (2007):
#' \deqn{\hat{g}(\tau)=\frac{\sum_{i,j} \{A_{t_i,t_j}|\tau_{ij}=\tau\}}{\sum_{i,j}\{A_i|\tau_{ij}=\tau\}},}
#' where where \eqn{A_{t_i,t_j}} is the number of associations of pair of individuals at both
#' time \eqn{t_i} and \eqn{t_j}, and \eqn{A_{t_i}} represents the number of associations
#' of pair of individuals at time \eqn{t_i}.
#'
#'
#' @export
#' @rdname lar_nonparametric_estimation
#'
lar_nonparametric_estimation <- function(X, tp, group_id = NULL) {
  tp <- tp-min(tp)+1
  if(is.list(X)){
    if(is.null(group_id)) {
      stop('please input group_id')
    }else{
      matrix_data <- list_to_matrix(X, tp, group_id)
    }
  }else{
    data <- X
  }
  len <- length(tp)
  tT <- max(tp)
  g_m <- rep(0, tT)
  g_n <- rep(0, tT)
  Aij <- c()
  Ai <- c()
  tauij <- c()
  p <- 1
  for (i in 1:(len-1)) {
    for (j in (i+1):len){
      tauij[p] <- tp[j] - tp[i]
      num <- data[, i]*data[, j]*(data[, i]+data[, j])
      Aij[p] <- sum(choose(table(num[which(num!=0)]), 2))
      if(sum(num!=0)!=0){
        obj <- data[which(num!=0), i]
        s <- 0
        for(g in 1:length(obj)){
          s <- s + sum(data[which(num!=0)[g]:dim(data)[1],i]==obj[g])
        }
        Ai[p] <- s - length(obj)
      }else{
        Ai[p] <- 0
      }
      g_m[tp[j]-tp[i]] <- g_m[tp[j]-tp[i]] + Aij[p]
      g_n[tp[j]-tp[i]] <- g_n[tp[j]-tp[i]] + Ai[p]
      p <- p + 1
    }
  }
  gtauij <- g_m[tauij]/g_n[tauij]
  tau <- unique(tauij)
  g_tau <- g_m[tau]/g_n[tau]
  g_data <- list(tau=tau, g_tau=g_tau, g_m=g_m, g_n=g_n,
                Aij=Aij, Ai=Ai, tauij=tauij)
  return(g_data)
}


