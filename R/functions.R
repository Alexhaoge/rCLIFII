#' @name functions
#' @title Preparation functions
#'
#' @param X Simulation data.
#' @param n The number of sampling individuals (in LIR) or groups (in LAR).
#' @param tp A set of observed time.
#' @param group_id Groups of individuals.
#' @param R_m Sum of mij.
#' @param R_n Sum of nij.
#' @param g_m Sum of Aij.
#' @param g_n Sum of Ai.
#' @param mij The number of individuals identified at both time t_i and t_j.
#' @param nij The number of individuals identified at time t_i.
#' @param Aij The number of observed associated individuals (in pairs) observed at both time t_i and t_j.
#' @param Ai The number of observed associated individuals (in pairs) at time t_i.
#' @param tauij Time lag between time t_i and t_j.
#' @param mtau The maximum allowable lag time.
#' @param np Expectation value of mij (in LIR) or Aij (in LAR).
#'
#' @export
#' @rdname functions
#'
.data.clean <- function(n, tp){
  observed_time <- tp
  time_point_number <- length(observed_time)
  tp <- observed_time-min(observed_time)+1
  tT <- max(tp)
  if (time_point_number == 1) warning("Only one observation?")
  if (length(n) == 1) n <- rep(n, tT)
  if(length(n)==time_point_number) {
    n0 <- rep(0, tT)
    n0[tp] <- n
    n <- n0
  }
  res <- n
  return(res)
}
#' @export
#' @rdname functions
#'
list_to_matrix <- function(X, tp, group_id=NULL){
  observed_time <- tp
  individual_ids <- unique(unlist(X))
  individual_number <- length(individual_ids)
  time_point_number <- length(observed_time)
  data <- matrix(0, individual_number, time_point_number)
  rownames(data) <- individual_ids
  colnames(data) <- observed_time
  if(is.null(group_id)){
    for (i in seq_len(time_point_number)){
      idx <- rownames(data) %in% X[[observed_time[i]]]
      data[idx, i] <- 1
    }
  }else{
    for (i in seq_len(time_point_number)){
      idx <- match(X[[observed_time[i]]], individual_ids)
      data[idx, i] <- group_id[[observed_time[i]]][idx]
      }
  }
  data <- as.matrix(data)
  colnames(data) <- observed_time
  rownames(data) <- paste('ID', round(individual_ids, 0), sep='')
  return(data)
}
#' @export
#' @rdname functions
#'
data.frame_to_matrix <- function(X){
  input_data <- unique(X)
  individual_ids <- unique(input_data$ID)
  input_data$Date <- input_data$Date - min(input_data$Date) + 1
  observed_time <- unique(input_data$Date)
  individual_number <- length(individual_ids)
  time_point_number <- length(observed_time)
  data <- matrix(0, individual_number, time_point_number)
  rownames(data) <- individual_ids
  colnames(data) <- observed_time
  if(ncol(input_data)==2){
    for (i in seq_len(time_point_number)){
      match_ID <- input_data$ID[which(input_data$Date %in% input_data$Date[i])]
      idx <- which(individual_ids %in% match_ID)
      data[idx, i] <- 1
    }
  }else{
    match_ID <- input_data$ID[which(input_data$Date %in% input_data$Date[i])]
    idx <- which(individual_ids %in% match_ID)
    data[idx, i] <- input_data$Group[idx]
  }
  return(data)
}
#' @export
#' @rdname functions
#'
lir_estimation_c <- function(R_m, R_n, mij, nij, tauij, mtau = 1000){

  theta <- lir.model.res('mod3', mij, nij, tauij, mtau)$par
  R_tau_C <- function(tauij){theta[1]*exp(-theta[2]*tauij) + theta[3]}
  unique_tau <- unique(tauij)
  mij <- R_m[tauij][match(unique(tauij), tauij)]
  np <- R_n[tauij][match(unique(tauij), tauij)]*R_tau_C(unique_tau)

  new = g_np(mij, np, unique_tau)
  new_mij = new$new_mij
  new_np = new$new_np

  df = length(new_mij)
  c = (new_mij-new_np)^2 %*% matrix(1/new_np)/(df-3-1)
  c = ifelse(c>1 & c<3, c, 1)
}
#' @export
#' @rdname functions
#'
lar_estimation_c <- function(g_m, g_n, Aij, Ai, tauij, mtau = 1000){

  theta <- lar.model.res('mod6', Aij, Ai, tauij, mtau)$par
  g_tau <- function(tauij){((1-theta[1])*exp(-theta[2]*tauij)+theta[1])*exp(-theta[3]*tauij)}
  unique_tau <- unique(tauij)
  Aij <- g_m[tauij][match(unique(tauij), tauij)]
  np <- g_n[tauij][match(unique(tauij), tauij)]*g_tau(unique_tau)

  new = g_np(Aij, np, unique_tau)
  new_Aij = new$new_mij
  new_np = new$new_np

  df = length(new_Aij)
  c = (new_Aij-new_np)^2 %*% matrix(1/new_np)/(df-3-1)
  c = ifelse(c>1 & c<3, c, 1)
}

#' @export
#' @rdname functions
#'
#'
#'
g_np <- function(mij, np, tauij){
  unique_tau <- unique(tauij)
  num_np <- length(np)
  g <- rep(NA, num_np)
  idx <- seq_len(num_np)
  a = np
  for(i in idx){
    b = idx[which(cumsum(a) > 6)[1]]
    g[idx[which(idx <= b)]] = b
    a = a[-c(which(idx <= b))]
    idx = idx[-c(which(idx <= b))]
  }

  g[(1:num_np)[-which(g>0)]] <- range(g[which(g>0)])[2]
  group_data = data.frame(cbind(mij, np, unique_tau, g))

  new_data <- group_data %>% dplyr::group_by(g)
  new = new_data %>% dplyr::summarise(
    new_mij = sum(mij),
    new_np = sum(np)
  )
  return(new)
}
