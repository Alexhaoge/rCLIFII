
#' @name lir_model_selection
#' @title Model selection for models of lagged identification rate
#'
#' @param X A list or matrix containing the identities of individuals identified in each sampling period
#' @param model Models of lagged identification rate, model = 'lir_1', 'lir_2', or 'lir_3'.
#' @param n A vector or a positive integer, representing the number of individuals identified in each sampling period.
#' It indicates the same number of individuals identified in all sampling periods if a positive integer.
#' @param tp A set of observed time.
#' @param nboot The number of bootstrap samples desired
#' @param mtau The maximum allowable lag time.
#' @param ncores doParallel.
#' @param seed Random seed.
#' @details
#' See Akaike (1973) for Akaike information criterion (AIC);
#' See Burnhan and Anderson (2002) for Quasi-Akaike information criterion (QAIC);
#' See this paper for composite likelihood information criterion (CLIC).
#'
#' @return The values of model selection criteria.
#'
#'
#' @export
#' @rdname lir_model_selection


lir_model_selection <- function(X, n, tp, model, nboot, mtau = 1000, ncores = 4, seed = NULL){
  B <- as.integer(nboot)
  if (B <= 1){
    stop("nboot must be a positive integer bigger than 1")
  }
  tT <- max(tp-min(tp))
  len <- length(tp)
  n <- .data.clean(n, tp)
  lir_data <- lir_nonparametric_estimation(X, n, tp)
  R_m <- lir_data$R_m
  R_n <- lir_data$R_n
  mij <- lir_data$mij
  nij <- lir_data$nij
  tauij <- lir_data$tauij
  if(model == 'lir_1'){
    mod0 <- 'mod1'
  }
  if(model == 'lir_2'){
    mod0 <- 'mod2'
  }else{
    mod0 <- 'mod3'
  }
  model.H <- lir.model.res(mod0, mij, nij, tauij, mtau)$H
  model.est <- lir.model.res(mod0, mij, nij, tauij, mtau)$par
  model.val <- lir.model.res(mod0, mij, nij, tauij, mtau)$val
  model.K <- lir.model.res(mod0, mij, nij, tauij, mtau)$K
  if(ncores>1){
    cl <- parallel::makeCluster(ncores) # not to overload your computer
    doParallel::registerDoParallel(cl)
    RESULTS = foreach::`%dopar%`(foreach::foreach(i = seq_len(B), .combine = rbind), {
      sampboot <- lir_bootstrap(X, tp, seed = NULL)
      dat <- lir_nonparametric_estimation(sampboot, n, tp)
      mij <- dat$mij
      nij <- dat$nij
      tauij <- dat$tauij
      res <- lir.model.res(model=mod0, mij=mij, nij=nij, tauij=tauij, mtau)
      out <- res$par
      return(out)
    })
    parallel::stopCluster(cl)
  }else{
    RESULTS <- matrix(0, B, model.K)
    for(i in seq_len(B)){
      sampboot <- lir_bootstrap(X, tp, seed)
      dat <- lir_nonparametric_estimation(sampboot, n, tp)
      mij <- dat$mij
      nij <- dat$nij
      tauij <- dat$tauij
      res <- lir.model.res(model=mod0, mij=mij, nij=nij, tauij=tauij, mtau)
      out <- res$par
      RESULTS[i,] <- out
      }
    }
    dimpara <- sum(diag(model.H%*%(stats::var(RESULTS))))
    AIC <- 2*model.val + 2*model.K
    estimation_c <- lir_estimation_c(R_m, R_n, mij, nij, tauij, mtau)
    QAIC <- 2*model.val/estimation_c + 2*model.K
    CLIC <- 2*model.val + dimpara
    res <- data.frame(AIC=AIC, QAIC=QAIC, CLIC=CLIC)
    return(res)
}


#' @name lar_model_selection
#' @title Model selection for models of lagged association rate
#' @param X A list or matrix containing the identities of individuals within
#' study area, and the states or status of individuals during each sampling period.
#' @param model Models of lagged identification rate, model = 'lar_1', 'lar_2', or 'lar_3'.
#' @param block_list A block list for a series of observation time. For example,
#' block_list = list(c(1:5), c(51:55), c(101:105), c(501:505), c(601:605)).
#' @param group_id Groups of individuals. If X is a list, please input group_id. If X is a matrix,
#' this parameter can be skipped and takes the default `NULL` value.
#' @param nboot The number of bootstrap samples desired
#' @param mtau The maximum allowable lag time
#' @param ncores doParallel
#' @param seed Random seed
#' @details
#' See Akaike (1973) for Akaike information criterion (AIC);
#' See Burnhan and Anderson (2002) for Quasi-Akaike information criterion (QAIC);
#' See this paper for composite likelihood information criterion (CLIC).
#'
#' @return The values of model selection criteria.
#'
#'
#' @export
#' @rdname lar_model_selection


lar_model_selection <- function(X, model, block_list, nboot, group_id = NULL, mtau = 1000, ncores = 4, seed = NULL){
  B <- as.integer(nboot)
  if (B <= 1){
    stop("nboot must be a positive integer bigger than 1")
  }
  tp <- unlist(block_list)
  tT <- max(tp-min(tp))
  len <- length(tp)
  lar_data <- lar_nonparametric_estimation(X, tp, group_id)
  g_m <- lar_data$g_m
  g_n <- lar_data$g_n
  Aij <- lar_data$Aij
  Ai <- lar_data$Ai
  tauij <- lar_data$tauij
  if(model == 'lar_1'){
    mod0 <- 'mod4'
  }
  if(model == 'lar_2'){
    mod0 <- 'mod5'
  }else{
    mod0 <- 'mod6'
  }
  model.H <- lar.model.res(mod0, Aij, Ai, tauij, mtau)$H
  model.est <- lar.model.res(mod0, Aij, Ai, tauij, mtau)$par
  model.val <- lar.model.res(mod0, Aij, Ai, tauij, mtau)$val
  model.K <- lar.model.res(mod0, Aij, Ai, tauij, mtau)$K
  if(ncores>1){
    cl <- parallel::makeCluster(ncores) # not to overload your computer
    doParallel::registerDoParallel(cl)
    RESULTS = foreach::`%dopar%`(foreach::foreach(i = seq_len(B), .combine = rbind), {
      sampboot <- lar_bootstrap(X, block_list, group_id, seed)
      dat <- lar_nonparametric_estimation(sampboot, tp, group_id)
      Aij <- dat$Aij
      Ai <- dat$Ai
      tauij <- dat$tauij
      res <- lar.model.res(model=mod0, Aij=Aij, Ai=Ai, tauij=tauij, mtau)
      out <- res$par
      return(out)
    })
    parallel::stopCluster(cl)
  }else{
    RESULTS <- matrix(0, B, model.K)
    for(i in seq_len(B)){
      sampboot <- lar_bootstrap(X, block_list, group_id, seed)
      dat <- lar_nonparametric_estimation(sampboot, tp, group_id)
      Aij <- dat$Aij
      Ai <- dat$Ai
      tauij <- dat$tauij
      res <- lar.model.res(model=mod0, Aij=Aij, Ai=Ai, tauij=tauij, mtau)
      out <- res$par
      RESULTS[i,] <- out
    }
  }
  dimpara <- sum(diag(model.H%*%(stats::var(RESULTS))))
  AIC <- 2*model.val + 2*model.K
  estimation_c <- lar_estimation_c(g_m, g_n, Aij, Ai, tauij, mtau)
  QAIC <- 2*model.val/estimation_c + 2*model.K
  CLIC <- 2*model.val + dimpara
  res <- data.frame(AIC=AIC, QAIC=QAIC, CLIC=CLIC)
  return(res)
}

















