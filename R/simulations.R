
#' @name lir_simulation_A/B/C
#' @title Simulate three types of animal movement patterns
#'
#' @param Z A positive integer, representing the population size within a closed area, members of which can be either inside
#' or outside the study area.
#' @param N A positive integer, representing the population size within the study area.
#' @param n A vector or a positive integer, representing the number of individuals identified in each sampling period.
#' It indicates the same number of individuals identified in all sampling periods if a positive integer.
#' @param tp A set of observed time.
#' @param lambda The rate of emigration from a study area.
#' @param mu The rate of reimmigration from a study area.
#' @import ggplot2
#'
#' @details
#'
#' Type A (Closed): a population of N individuals present in the study area throughout,
#' with no immigration, emigration, birth, death.
#'
#' Type B (Permanent emigration): a population of N individuals in the study area with
#' permanent emigration at a rate of lambda per individual per time unit,
#' with departed individuals being replaced 1:1 by new individuals.
#'
#' Type C (Emigration plus reimmigration): a closed population of Z individuals,
#' members of which can be either inside or outside the study area. Individuals
#' in the study area leave the study area at a rate of lambda per individual per
#' time unit and individuals outside the study area reenter it with a probability
#' of mu per individual per time unit.
#'
#' @return Simulation data for Type A-C.
#'
#'
#' @export
#' @rdname lir_simulations
#'
#' @examples
#'
#'
#' # Example
#' # set population size and number of subsampling
#' Z <- 300; N <- 100; n <- 40
#' # set the observation time
#' tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
#' # set parameters
#' lambda_B <- 0.008; lambda_C <- 0.08; mu_C <- 0.04
#'
#' # simulation of Type A
#' simulation_A <- lir_simulation_A(N, n, tp)
#' lir_data <- lir_nonparametric_estimation(simulation_A, n, tp)
#' mij <- lir_data$mij
#' nij <- lir_data$nij
#' tauij <- lir_data$tauij
#' theta <- round(.mod1(mij, nij, tauij)$par, 8)
#' # a plot object
#' graph_data <- data.frame(tau = lir_data$tau, R_tau = lir_data$R_tau)
#' require(ggplot2)
#' ggplot(data = graph_data, aes(tau, R_tau)) +
#' geom_point() +
#'   scale_y_continuous(limits = c(0, 0.02)) +
#'   theme(axis.text = element_text(size = rel(1.2))) +
#'   labs(x = 'Time Lag', y = 'Lagged Identification Rates', title = 'Simulation_A') +
#'   geom_hline(yintercept = theta, linetype = 2, color = 'blue')
#'
#' # simulation of Type B
#' simulation_B <- lir_simulation_B(N, n, lambda_B, tp)
#' lir_data <- lir_nonparametric_estimation(simulation_B, n, tp)
#' mij <- lir_data$mij
#' nij <- lir_data$nij
#' tauij <- lir_data$tauij
#' theta <- round(.mod2(mij, nij, tauij)$par, 8)
#' alpha <- theta[1]; beta <- theta[2]
#' graph_data <- data.frame(tau = lir_data$tau, R_tau = lir_data$R_tau)
#' tT <- max(tp-min(tp))
#' line_data <- data.frame(tau = seq_len(tT), R_tau = alpha*exp(-beta*seq_len(tT)))
#' require(ggplot2)
#' ggplot(data = graph_data, aes(tau, R_tau)) +
#'   geom_point() +
#'   theme(axis.text = element_text(size = rel(1.2))) +
#'   labs(x = 'Time Lag', y = 'Lagged Identification Rates', title = 'Simulation_B') +
#'   geom_line(aes(tau, R_tau), line_data, linetype = 2, color='blue')
#'
#' # simulation of Type C
#' simulation_C <- lir_simulation_C(Z, N, n, lambda_C, mu_C, tp)
#' lir_data <- lir_nonparametric_estimation(simulation_C, n, tp)
#' mij <- lir_data$mij
#' nij <- lir_data$nij
#' tauij <- lir_data$tauij
#' theta <- round(.mod3(mij, nij, tauij)$par, 8)
#' alpha <- theta[3]; beta <- theta[2]; gamma <- theta[1]
#' graph_data <- data.frame(tau = lir_data$tau, R_tau = lir_data$R_tau)
#' tT <- max(tp-min(tp))
#' line_data <- data.frame(tau = seq_len(tT), R_tau = gamma*exp(-beta*seq_len(tT))+alpha)
#' require(ggplot2)
#' ggplot(data = graph_data, aes(tau, R_tau)) +
#'   geom_point() +
#'   theme(axis.text = element_text(size = rel(1.2))) +
#'   labs(x = 'Time Lag', y = 'Lagged Identification Rates', title = 'Simulation_C') +
#'   geom_line(aes(tau, R_tau), line_data, linetype = 2, color='blue')
#'
#'
#'

lir_simulation_A <- function(N, n, tp) {
  tp <- tp-min(tp)+1
  tT <- max(tp)
  n <- .data.clean(n, tp)
  list_data <- lapply(n, function(ni){ sample(1:N, ni) })
  matrix_data <- list_to_matrix(list_data, tp)
  simulation_A <- matrix_data
  return(simulation_A)
}

#' @rdname lir_simulations
#' @export
#'
lir_simulation_B <- function(N, n, lambda, tp) {
  tp <- tp-min(tp)+1
  tT <- max(tp)
  n <- .data.clean(n, tp)
  pop <- seq(N)
  Npop <- list()
  Npop[[1]] <- pop
  for(i in 2:tT){
    unif <- stats::runif(N, 0, 1)
    Npop[[i]] <- (unif<=lambda)*unif + Npop[[i-1]]
  }
  data <- list()
  for (i in tp) {
    data[[i]] <- sample(Npop[[i]], n[i])
  }
  matrix_data <- list_to_matrix(data, tp)
  simulation_B <- matrix_data
  return(simulation_B)
}

#' @rdname lir_simulations
#' @export
#'
#'
lir_simulation_C <- function(Z, N, n, lambda, mu, tp) {
  tp <- tp-min(tp)+1
  tT <- max(tp)
  n <- .data.clean(n, tp)
  pop <- seq(N)
  pop <- seq(N)
  Z1 <- 1:N
  Z2 <- (N+1):(Z)
  ### updata
  popin_new <- list()
  popin_new[[1]] <- Z1
  popout_new <- list()
  popout_new[[1]] <- Z2
  for(i in 2:tT) {
    unif1 <- stats::runif(length(popin_new[[i-1]]))
    unif2 <- stats::runif(length(popout_new[[i-1]]))
    popin_new[[i]] <- c(popin_new[[i-1]][which(unif1>lambda)],popout_new[[i-1]][which(unif2<=mu)])
    popout_new[[i]] <- c(popin_new[[i-1]][which(unif1<=lambda)],popout_new[[i-1]][which(unif2>mu)])
  }
  ### sample
  data <- list()
  for (i in tp) {
    data[[i]] <- sample(popin_new[[i]], n[i])
  }
  matrix_data <- list_to_matrix(data, tp)
  simulation_C <- matrix_data
  return(simulation_C)
}

#' @name lar_simulation_D/E/F
#' @title Simulate three types of animal social structure patterns.
#'
#' @param N Population size within the study area.
#' @param n A vector or a positive integer, representing the number of groups identified in each sampling period.
#' It indicates the same number of groups identified in all sampling periods if a positive integer.
#' @param W Number of groups formed by individuals.
#' @param U Number of permanent social units allocated by individuals.
#' @param tp A set of observed time.
#' @param lambda The probability that an individual changes the group per time unit.
#'
#' @import ggplot2
#'
#' @details
#'
#' Type D (Random): individuals form W groups, with random allocation of individuals
#' to groups at each observation.
#'
#' Type E (Casual acquaintances): individuals form W groups, and change groups with
#' probability lambda per time unit.
#'
#' Type F (Permanent companions plus casual acquaintances): individuals are randomly
#' allocated to U permanent social units, and these are randomly allocated to W groups,
#' with units changing groups with probability lambda per time unit-permanent companions
#' plus casual acquaintances.
#'
#' @return Simulation data for Type D-F.
#'
#'
#' @export
#' @rdname lar_simulations
#'
#' @examples
#'
#'
#' # Example
#' # set population size and number of subsampling
#' N <- 100; n <- 4; W <- 10; U <- 20
#' # set the observation time
#' tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
#' # set parameters
#' lambda_E <- 0.008; lambda_F <- 0.008
#'
#' # simulation of Type D
#' simulation_D <- lar_simulation_D(N, n, W, tp)
#' lar_data <- lar_nonparametric_estimation(simulation_D, tp)
#' Aij <- lar_data$Aij
#' Ai <- lar_data$Ai
#' tauij <- lar_data$tauij
#' theta <- round(.mod4(Aij, Ai, tauij)$par, 8)
#'
#' # a plot object
#' graph_data <- data.frame(tau=lar_data$tau, g_tau=lar_data$g_tau)
#' require(ggplot2)
#' ggplot(data = graph_data, aes(tau, g_tau)) +
#' geom_point() +
#'   theme(axis.text = element_text(size = rel(1.2))) +
#'   labs(x="Time Lag", y="Lagged Association Rates", title='Simulation_D') +
#'   geom_hline(yintercept = theta, linetype = 2, color='blue')
#'
#' # simulation of Type E
#' simulation_E <- lar_simulation_E(N, n, W, lambda_E, tp)
#' lar_data <- lar_nonparametric_estimation(simulation_E, tp)
#' Aij <- lar_data$Aij
#' Ai <- lar_data$Ai
#' tauij <- lar_data$tauij
#' theta <- round(.mod5(Aij, Ai, tauij)$par, 8)
#' alpha <- theta[1]
#' beta <- theta[2]
#' graph_data <- data.frame(tau=lar_data$tau, g_tau=lar_data$g_tau)
#' tT <- max(tp-min(tp))
#' line_data <- data.frame(tau=seq_len(tT), g_tau=(1-alpha)*exp(-beta*seq_len(tT))+alpha)
#' require(ggplot2)
#' ggplot(data = graph_data, aes(tau, g_tau)) +
#'   geom_point() +
#'   theme(axis.text = element_text(size = rel(1.2))) +
#'   labs(x="Time Lag", y="Lagged Association Rates", title='Simulation_E') +
#'   geom_line(aes(tau, g_tau), line_data, linetype = 2, color='blue')
#'
#' # simulation of Type F
#' simulation_F <- lar_simulation_F(N=100, n=4, W=10, U=20, lambda_F, tp)
#' lar_data <- lar_nonparametric_estimation(simulation_F, tp)
#' Aij <- lar_data$Aij
#' Ai <- lar_data$Ai
#' tauij <- lar_data$tauij
#' theta <- round(.mod5(Aij, Ai, tauij)$par, 8)
#' alpha <- theta[1]
#' beta <- theta[2]
#' graph_data <- data.frame(tau=lar_data$tau, g_tau=lar_data$g_tau)
#' tT <- max(tp-min(tp))
#' line_data <- data.frame(tau=seq_len(tT), g_tau=(1-alpha)*exp(-beta*seq_len(tT))+alpha)
#' ggplot(data = graph_data, aes(tau, g_tau)) +
#'   geom_point() +
#'   theme(axis.text = element_text(size = rel(1.2))) +
#'   labs(x="Time Lag", y="Lagged Association Rates", title='Simulation_F') +
#'   geom_line(aes(tau, g_tau), line_data, linetype = 2, color='blue')
#'
lar_simulation_D <- function(N, n, W, tp) {
  tp <- tp-min(tp)+1
  tT <- max(tp)
  len <- length(tp)
  n <- .data.clean(n, tp)
  pop <- seq(N)
  groupid <- list()
  samp.d <- list()
  for (i in tp) {
    groupid[[i]] <- sample(W, N, replace = TRUE) # Individuals randomly form W groups
    g <- sample(W, n[i]) # n groups are randomly selected from W groups
    samp.d[[i]] <- pop[which(groupid[[i]] %in% g)] # The members contained in the selected groups
  }
  ### Generate data matrix
  data <- matrix(0, N, len)
  k <- 1
  for(t in tp){
    data[match(samp.d[[t]], pop), k] <- groupid[[t]][match(samp.d[[t]], pop)]
    k = k + 1
  }
  matrix_data <- data
  colnames(matrix_data) <- tp
  rownames(matrix_data) <- pop
  simulation_D <- matrix_data
  return(simulation_D)
}


#' @rdname lar_simulations
#' @export
#'
#'
lar_simulation_E <- function(N, n, W, lambda, tp) {
  tp <- tp-min(tp)+1
  tT <- max(tp)
  len <- length(tp)
  n <- .data.clean(n, tp)
  pop <- seq(N)
  groupid <- sample(W, N, replace = TRUE)
  groupid_new <- list()
  for( i in 1:tT) {
    groupid_new[[i]] <- groupid
  }
  for( i in 2:tT) {
    gunif <- stats::runif(N)
    for (j in 1:N) {
      groupid_new[[i]][j] <- (gunif[j]>=lambda)*groupid_new[[i-1]][j] +
        (gunif[j]<lambda)*sample(seq(W)[-groupid_new[[i-1]][j]],1)
    }
  }
  samp.e <- list()
  for (i in tp) {
    g <- sample(W, n[i])
    samp.e[[i]] <- pop[which(groupid_new[[i]]%in%g)]
  }
  ### obtain data matrix
  data <- matrix(0, N, len)
  k <- 1
  for(t in tp){
    data[match(samp.e[[t]],pop),k] <- groupid_new[[t]][match(samp.e[[t]],pop)]
    k = k + 1
  }
  matrix_data <- data
  colnames(matrix_data) <- tp
  rownames(matrix_data) <- pop
  simulation_E <- matrix_data
  return(simulation_E)
}


#' @rdname lar_simulations
#' @export
#'
#'
lar_simulation_F <- function(N, n, W, U, lambda, tp) {
  tp <- tp-min(tp)+1
  tT <- max(tp)
  len <- length(tp)
  n <- .data.clean(n, tp)
  pop <- seq(N)
  unitid <- c(sample(U,U), sample(U, N-U, replace = TRUE))  # 100 individuals belong to 20 units
  groupid <- c(sample(W,W), sample(W, U-W, replace = TRUE)) # 20 units belong to 10 groups
  groupid_new <- list()
  for( i in 1:tT) {
    groupid_new[[i]] <- groupid
  }
  for( i in 2:tT) {
    gunif <- stats::runif(U)
    for (j in 1:U) {
      groupid_new[[i]][j]<-(gunif[j]>=lambda)*groupid_new[[i-1]][j] +
        (gunif[j]<lambda)*sample(seq(W)[-groupid_new[[i-1]][j]], 1)
    }
  }
  samp.f <- list()
  for (i in tp) {
    g <- sample(W, n[i])
    samp.f[[i]] <- pop[unitid %in% which(groupid_new[[i]]%in%g)]
  }
  data <- matrix(0, N, len)
  k <- 1
  for(t in tp){
    data[match(samp.f[[t]], pop), k] <- groupid_new[[t]][unitid[match(samp.f[[t]], pop)]]
    k = k + 1
  }
  matrix_data <- data
  colnames(matrix_data) <- tp
  rownames(matrix_data) <- pop
  simulation_F <- matrix_data
  return(simulation_F)
}





