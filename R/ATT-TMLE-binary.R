## TMLE estimation for ATT under univariate binary mediator  ====
#' Function for implementing TMLE for univariate binary mediators.
#' @param a Treatment level at which the average counterfactual outcome is computed
#' @param data A Dataframe contains treatment, mediators, outcome, and measured confounders
#' @param treatment Variable name for the unvariate binary treatment
#' @param mediator Variable name for the binary univariate mediator
#' @param outcome Variable name for the continuous univariate outcome
#' @param covariates Variable name for the measured confounders
#' @param onestep A logical indicator determines whether one-step estimation is executed. When 'onestep=T', the one-step estimation result is provided.
#' Conversely, if 'onestep=F', the result is withheld.
#' @param superlearner A logical indicator determines whether SuperLearner is adopted for estimating the outcome regression, mediator density, and the propensity score.
#' @param crossfit A logical indicator determines whether SuperLearner+Cross-fitting is adopted for estimating the outcome regression, mediator density, and the propensity score.
#' @param K A integer indicating the number of folds for cross-fitting, the default is 5.
#' @param lib Library of algorithms for SuperLearner or SuperLearner+Cross-fitting.
#' @param n.iter The maximum number of iterations performed when iteratively updating the mediator density and propensity score.
#' @param cvg.criteria A numerical value representing the convergence criteria when iteratively updating the mediator density and propensity score. The default value is 0.01, meaning update stops when \eqn{max(|\Phi_M|,|\Phi_A|)<0.01}.
#' @param formulaY Regression formula for the outcome regression of Y on M, A, X. The default is 'Y ~ 1+ M + A + X'.
#' @param linkY_binary The link function used for the logistic regression of Y on M, A, X when Y is binary. The default is the 'logit' link.
#' @param formulaA Regression formula for the propensity score regression of A on X. The default is 'A ~ 1 + X'.
#' @param linkA The link function used for the logistic regression of A on X. The default is the 'logit' link.
#' @param formulaM Regression formula for the mediator density regression of M on A and X. The default is 'M ~ 1 + A + X'. This parameter is only needed when M is a univariate binary mediator.
#' @param linkM_binary The link function used for the logistic regression of M on A and X. The default is the 'logit' link. This parameter is only needed when M is a univariate binary mediator.
#' @param truncate_lower A numeric variable, setting lower bound for the truncated propensity score. The default is 0.
#' @param truncate_upper A numeric variable, setting upper bound for the truncated propensity score. The default is 1.
#' @return a list of initialization of matrices.
#' @examples
#' \donttest{
#' ATT.TMLE.binary(a=1,data=continuousY_binaryM,
#' treatment="A", mediator="M", outcome="Y", covariates="X")
#' }
#' @import SuperLearner
#' @importFrom dplyr %>% mutate select
#' @return Function outputs a list containing TMLE output (and Onestep estimator output if 'onestep=T' is specified):
#' \describe{
#'       \item{\code{estimated_psi}}{The estimated parameter of interest: \eqn{E(Y^a|1-a)}}
#'       \item{\code{lower.ci}}{Lower bound of the 95\% confidence interval for \code{estimated_psi}}
#'       \item{\code{upper.ci}}{Upper bound of the 95\% confidence interval for \code{estimated_psi}}
#'       \item{\code{theta_x}}{\eqn{\int E(Y|M,A,X)p(M|A=a,X)p(A|X) dM dA}}
#'       \item{\code{p.m1.aX}}{\eqn{\int p(M=1|A=a,X)}}
#'       \item{\code{p.a1.X}}{\eqn{p(A=1|X)}}
#'       \item{\code{or_pred}}{\eqn{E(Y|M,A,X)}}
#'       \item{\code{EIF}}{The estimated efficient influence function evaluated at the observed data}
#'       \item{\code{EDstar}}{A vector of the mapping of \code{EIF} in the tangent space of \eqn{Y|M,A,X}, \eqn{M|A,X}, and \eqn{A|X}.}
#'       \item{\code{EDstar.M_vec}}{A vector containing the average value of the mapping of EIF in tangent space \eqn{M|A,X} over iterations. This is useful for checking the convergence behavior of the mediator density. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{EDstar.ps_vec}}{A vector containing the average value of the mapping of EIF in tangent space \eqn{A|X} over iterations. This is useful for checking the convergence behavior of the propensity score. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{eps.M_vec}}{A vector containing the index for submodels of the mediator density over iterations. This is useful for checking the convergence behavior of the mediator density. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{eps.Y_vec}}{A vector containing the index for submodels of the propensity score over iterations. This is useful for checking the convergence behavior of the propensity score. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{iter}}{Number of iterations where convergence is achieved for the iterative update of the mediator density and propensity score.}}
#' @export
#'

ATT.TMLE.binary <- function(a,data,treatment, mediator, outcome, covariates,
                        onestep=TRUE, superlearner=TRUE,crossfit=FALSE,K=5,
                        lib = c("SL.glm","SL.earth","SL.ranger","SL.mean"), n.iter=500, cvg.criteria=0.01,
                        formulaY="Y ~ .", formulaA="A ~ .", formulaM="M~.", linkY_binary="logit", linkA="logit", linkM_binary="logit",
                        truncate_lower=0, truncate_upper=1){

  # attach(data, warn.conflicts=FALSE)

  n <- nrow(data)

  alt <- 1-a # alternative treatment level

  # Variables
  A <- data[,treatment]
  M <- data[,mediator]
  X <- data[,covariates,drop = FALSE]
  Y <- data[,outcome]

  # new data sets
  data_Aa = data.frame(A=a,X)
  data_Aa1 = data.frame(A=1,X)
  data_Aa0 = data.frame(A=0,X)
  data_Aalt = data.frame(A=alt,X)

  data_A1 = data.frame(M, A=1, X)
  data_A0 = data.frame(M, A=0, X)

  data_MAX = data.frame(M=M,A=A,X)
  data_MaX = data.frame(M,A=a,X)
  data_A1M1 = data.frame(M=1,A=1, X)
  data_A1M0 = data.frame(M=0,A=1, X)
  data_A0M1 = data.frame(M=1,A=0, X)
  data_A0M0 = data.frame(M=0,A=0, X)

  ## p(A) ====
  p.a <- mean(A==a)
  p.alt <- mean(A==alt)

  ## outcome regression ====
  Y.binary <- all(Y %in% c(0,1) ) # whether Y is binary
  Y.family <- if (Y.binary){binomial()}else{gaussian()}


  if (crossfit==T){ # use cross fitting

    or_fit <- CV.SuperLearner(Y=Y, X=data_MAX, family = Y.family, V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
    or_pred <- or_fit$SL.predict # E(Y|M,A,X)
    or_pred_alt <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data.frame(M=M, A=alt,X)[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))] # E(Y|M,a',X)

    # E(Y|m,a',X), m = {0,1}
    or_pred_altm1 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data.frame(M=1, A=alt, X)[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
    or_pred_altm0 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data.frame(M=0, A=alt, X)[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]


  }else if (superlearner==T){ # use superlearner

    or_fit <- SuperLearner(Y=Y, X=data_MAX, family = Y.family, SL.library = lib)
    or_pred <- predict(or_fit)[[1]] %>% as.vector() # E(Y|M,A,X)
    or_pred_alt <- predict(or_fit, newdata=data.frame(M=M,A=alt,X))[[1]] %>% as.vector() # E(Y|M,a',X)

    # E(Y|m,a',X), m = {0,1}
    or_pred_altm1 <- predict(or_fit, newdata=data.frame(M=1,A=alt,X))[[1]] %>% as.vector()
    or_pred_altm0 <- predict(or_fit, newdata=data.frame(M=0,A=alt,X))[[1]] %>% as.vector()


  }else{ # use simple regression

    Y.family <- if (all(Y %in% c(0,1) )){binomial(linkY_binary)}else{gaussian()}

    or_fit <- glm(as.formula(formulaY), data=data_MAX, family = Y.family)
    or_pred <- predict(or_fit, type="response") # E(Y|M,A,X)
    or_pred_alt <- predict(or_fit, newdata=data.frame(M=M,A=alt,X), type="response") # E(Y|M,a',X)


    # E(Y|m,a',X), m = {0,1}
    or_pred_altm1 <- predict(or_fit, newdata=data.frame(M=1,A=alt,X), type="response")
    or_pred_altm0 <- predict(or_fit, newdata=data.frame(M=0,A=alt,X), type="response")

  }



  ## propensity score ====
  if (crossfit==T){

    ps_fit <- CV.SuperLearner(Y=A, X=X, family = binomial(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
    p.a1.X <- ps_fit$SL.predict

  }else if (superlearner==T){

    ps_fit <- SuperLearner(Y=A, X=X, family = binomial(), SL.library = lib)
    p.a1.X <- predict(ps_fit, type = "response")[[1]] %>% as.vector()  # p(A=1|X)

  } else {

    if (truncate_lower!=0 | truncate_upper!=1){ # under weak overlapping issue, it's more stable to run A~X via linear regression

      ps_fit <- lm(as.formula(formulaA), data=X)
      p.a1.X <- predict(ps_fit)

      print("Truncation performed.")

    } else { # without weak overlapping issue. Run A~X via logistic regression

      ps_fit <- glm(as.formula(formulaA), data=X,  family = binomial(linkA))
      p.a1.X <- predict(ps_fit, type = "response")  # p(A=1|X)

    }
  }

  # apply truncation to propensity score to deal with weak overlap. Truncated propensity score within the user specified range of [truncate_lower, truncate_upper]: default=[0,1]
  p.a1.X[p.a1.X < truncate_lower] <- truncate_lower
  p.a1.X[p.a1.X > truncate_upper] <- truncate_upper


  p.a.X <- a*p.a1.X + (1-a)*(1-p.a1.X) # p(a|X)
  p.alt.X <- 1-p.a.X # p(a'|X)

  # estimate the density ratio of p(a,X)/p(a',X) using regression instead of the densratio package to make it more accurate and stable
  # p(a'|X)/p(a|X)
  ratio.A.X <- p.alt.X/p.a.X


  ## mediator density p(M|A,X) ====
  if (crossfit==T){

    M_fit <- CV.SuperLearner(Y=M, X=data.frame(A=A,X), family = binomial(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
    p.m1.AX <- M_fit$SL.predict # p(M=1|A,X)

    p.m1.aX <- unlist(lapply(1:K, function(x) predict(M_fit$AllSL[[x]], newdata=data.frame(A=a,X)[M_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) M_fit$folds[[x]])))] # p(M=1|A=a,X)
    p.m0.aX <- 1 - p.m1.aX # p(M=0|A=a,X)

    p.m1.altX <- unlist(lapply(1:K, function(x) predict(M_fit$AllSL[[x]], newdata=data.frame(A=alt,X)[M_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) M_fit$folds[[x]])))] # p(M=1|A=a',X)
    p.m0.altX <- 1 - p.m1.altX # p(M=0|A=a',X)

  }else if (superlearner==T){

    M_fit <- SuperLearner(Y=M, X=data.frame(A=A,X), family = binomial(), SL.library = lib)
    p.m1.AX <- predict(M_fit, type = "response")[[1]] %>% as.vector()  # p(M=1|A,X)

    p.m1.aX <- predict(M_fit, newdata=data.frame(A=a,X), type = "response")[[1]] %>% as.vector() # p(M=1|A=a,X)
    p.m0.aX <- 1 - p.m1.aX # p(M=0|A=a,X)

    p.m1.altX <- predict(M_fit, newdata=data.frame(A=alt,X), type = "response")[[1]] %>% as.vector() # p(M=1|A=a',X)
    p.m0.altX <- 1 - p.m1.altX # p(M=0|A=a',X)

  }else{

    M_fit <- glm(as.formula(formulaM), data=data.frame(A=A,X), family = binomial(linkM_binary))
    p.m1.AX <- predict(M_fit, type = "response")  # p(M=1|A,X)

    p.m1.aX <- predict(M_fit, newdata=data.frame(A=a,X), type = "response") # p(M=1|A=a,X)
    p.m0.aX <- 1 - p.m1.aX # p(M=0|A=a,X)

    p.m1.altX <- predict(M_fit, newdata=data.frame(A=alt,X), type = "response") # p(M=1|A=a',X)
    p.m0.altX <- 1 - p.m1.altX # p(M=0|A=a',X)

  }

  p.M.aX <- M*p.m1.aX + (1-M)*p.m0.aX # p(M|a,X)
  p.M.altX <- M*p.m1.altX + (1-M)*p.m0.altX # p(M|a',X)


  ############################################################### One step estimator ###############################################################
  ## one-step estimator
  if (onestep==T){

    # E(Dstar) for Y|M,A,X
    or_weight <- (A==alt)*{1/p.alt}*p.M.aX*{1/p.M.altX} # I(A=a')/p(a') * p(M|A=a,X)/p(M|A=a',X)
    EDstar_or <- mean(or_weight*(Y-or_pred_alt)) # weight*(Y - E(Y|M,a',X))

    # E(Dstar) for M=1|A,X
    m_weight <- (A==a)*{1/p.alt}*ratio.A.X*(or_pred_altm1-or_pred_altm0) # I(A=a)/p(a') * p(a'|X)/p(a|X) * E(Y|M=1,A=a',X) - E(Y|M=0,A=a',X)

    EDstar_M <- mean(m_weight*(M-p.m1.aX)) # weight*(E(Y|M=1,A=a',X) - E{ E(Y|A=a',M,X) | a, X})

    # at tangent space of A,X
    ax_weight <- (A==alt)*{1/p.alt}
    kappa <- or_pred_altm1*p.m1.aX + or_pred_altm0*p.m0.aX # E{ E(Y|A=a',M,X) | a, X}


    ######################
    # estimate E[Y(a)]
    ######################

    # estimated psi
    estimated_psi = EDstar_or + EDstar_M + mean(ax_weight*kappa) # estimated psi

    # EIF
    EIF <- or_weight*(Y-or_pred_alt)+ #line 1 of the EIF equation
      m_weight*(M - p.m1.aX)+ #line 2 of the EIF equation
      ax_weight*(kappa - estimated_psi) #line 3 of the EIF equation

    # confidence interval
    lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
    upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))

    onestep.out <- list(estimated_psi=estimated_psi, # estimated parameter
                        lower.ci=lower.ci, # lower bound of 95% CI
                        upper.ci=upper.ci, # upper bound of 95% CI
                        p.m1.aX=p.m1.aX,  # estimated M=1|A,X
                        p.a1.X=p.a1.X,  # estimated A=1|X
                        or_pred_alt=or_pred_alt, # estimated E(Y|M,A,X)
                        #
                        EIF=EIF, # EIF
                        EDstar=c(EDstar_or,EDstar_M)) # E(Dstar) for Y|M,A,X and M|A,X, and A,X

  }



  ############################################################### End of One step estimator ###############################################################




  ############################################################### TMLE estimator ###############################################################

  ## TMLE initialize
  eps.M <- 1 # submodel parameter for mediator density p(M|A,X)
  eps.Y <- 1 # submodel parameter for outcome regression E(Y|M,A,X)


  # record values of eps.M and eps.Y over iterations
  eps.M_vec <- vector(mode = "numeric")
  eps.Y_vec <- vector(mode = "numeric")

  # record values of EDstar over iterations
  EDstar.M_vec <- vector(mode = "numeric")
  EDstar.Y_vec <- vector(mode = "numeric")

  # place holder for clever coefficient
  clever_coef.M <- 0 # clever coefficient for mediator density p(M|A,X)
  clever_coef.Y <- 0 # clever coefficient for outcome regression E(Y|M,A,X)

  # cumulative summation of eps*clever_coef over iterations
  clever_sum.M <- 0
  clever_sum.Y <- 0

  # record average EIF in the sub-tangent space corresponding to M and Y
  EDstar_or <- 1 # initial value
  EDstar_M <- 1 # initial value

  # initialize updated nuisances: p(M=1|a,X), E(Y|M,a',X)
  p.m1.aX_updated <- p.m1.aX # p(M=1|A=a,X)
  or_pred_alt_updated <- or_pred_alt # E(Y|M,A=a',X)

  p.M.aX_updated <- p.M.aX # p(M|a,X)
  p.m0.aX_updated <- p.m0.aX # p(M=0|A=a,X)
  or_pred_altm1_updated <- or_pred_altm1 # E(Y|M=1,A=a',X)
  or_pred_altm0_updated <- or_pred_altm0 # E(Y|M=0,A=a',X)


  # initialize iteration
  iter <- 0

  ## if Y is binary, then iteration is needed during TMLE
  if(Y.binary){

    while ( abs(EDstar_M+EDstar_or)>cvg.criteria & iter<n.iter) {

      ######################
      # Update E(Y|M,alt,X)
      ######################

      ## clever coefficient for E(Y|M,alt,X)
      clever_coef.Y <- {1/p.alt}*p.M.aX_updated*{1/p.M.altX} # used for fitting the target regression, and updating E(Y|M,a',X)
      clever_coef.Y.m1 <- {1/p.alt}*p.m1.aX_updated*{1/p.m1.altX} # used for updating E(Y|M=1,a',X)
      clever_coef.Y.m0 <- {1/p.alt}*p.m0.aX_updated*{1/p.m0.altX} # used for updating E(Y|M=0,a',X)
      or_weight <- (A==alt)*clever_coef.Y  # used if calculating EIF

      # offset term for E(Y|M,alt,X). Only needed if Y is binary
      offset.Y <- qlogis(or_pred_alt_updated)

      model.Y <- glm( Y ~ offset(offset.Y)+clever_coef.Y-1, weights=(A==alt)*1, family=binomial() )

      eps.Y = coef(model.Y)
      eps.Y_vec <- c(eps.Y_vec,eps.Y)

      # update outcome regression
      or_pred_alt_updated <- plogis( qlogis(or_pred_alt_updated) + eps.Y*clever_coef.Y)

      # update related nuisances
      or_pred_altm1_updated <- plogis( qlogis(or_pred_altm1_updated) + eps.Y*clever_coef.Y.m1)
      or_pred_altm0_updated <- plogis( qlogis(or_pred_altm0_updated) + eps.Y*clever_coef.Y.m0)

      # E(Dstar) for E(Y|M,alt,X)
      EDstar_or <- mean(or_weight*(Y-or_pred_alt_updated))
      EDstar.Y_vec <- c(EDstar.Y_vec,EDstar_or)

      ######################
      # update p(M|A=a,X)
      ######################

      ## clever coefficient for M|A=a,X
      clever_coef.M = {1/p.alt}*ratio.A.X*(or_pred_altm1_updated-or_pred_altm0_updated)

      # offset term for M|A=a,X
      offset.M <- qlogis(p.m1.aX_updated)


      # derive eps.M
      model.M <- glm( M ~ offset(offset.M)+clever_coef.M-1, weights=(A==a)*1, family=binomial() )

      eps.M <-  coef(model.M)
      eps.M_vec <- c(eps.M_vec,eps.M)

      # update p(M=1|A=a,X)
      p.m1.aX_updated <- plogis(qlogis(p.m1.aX_updated)+eps.M*clever_coef.M)

      # update related nuisances
      p.m0.aX_updated <- 1-p.m1.aX_updated # p(M=0|a,X)
      p.M.aX_updated <- M*p.m1.aX_updated + (1-M)*p.m0.aX_updated # p(M|A=a,X)

      # E(Dstar) for M|A,X
      EDstar_M <- mean((A==a)*clever_coef.M*(M-p.m1.aX_updated))
      EDstar.M_vec <- c(EDstar.M_vec, EDstar_M)

      # estimated psi
      ax_weight <- (A==alt)*{1/p.alt}
      kappa <- or_pred_altm1_updated*p.m1.aX_updated + or_pred_altm0_updated*p.m0.aX_updated

      estimated_psi = mean(ax_weight*kappa)


      # update iteration
      iter <- iter+1

    } # end of while loop



  }




  ## if Y is continuous, then No iteration is needed during TMLE
  if(!Y.binary){

    ######################
    # update p(M|A=a,X)
    ######################

    ## clever coefficient for M|A=a,X
    clever_coef.M = {1/p.alt}*ratio.A.X*(or_pred_altm1-or_pred_altm0)

    # offset term for M|A=a,X
    offset.M <- qlogis(p.m1.aX)

    # derive eps.M
    model.M <- glm( M ~ offset(offset.M)+clever_coef.M-1, weights=(A==a)*1, family=binomial() )

    eps.M <-  coef(model.M)

    # update p(M=1|A=a,X)
    p.m1.aX_updated <- plogis(qlogis(p.m1.aX_updated)+eps.M*clever_coef.M)

    # update related nuisances
    p.m0.aX_updated <- 1-p.m1.aX_updated # p(M=0|a,X)
    p.M.aX_updated <- M*p.m1.aX_updated + (1-M)*p.m0.aX_updated # p(M|A=a,X)

    # calculated E(Dstar) for M|A=a,X
    m_weight <- (A==a)*{1/p.alt}*ratio.A.X*(or_pred_altm1-or_pred_altm0) # I(A=a)/p(a') * p(a'|X)/p(a|X) * E(Y|M=1,A=a',X) - E(Y|M=0,A=a',X)

    EDstar_M <- mean(m_weight*(M-p.m1.aX_updated))

    ######################
    # Update E(Y|M,alt,X)
    ######################

    ## clever coefficient for E(Y|M,alt,X)
    clever_coef.Y <- {1/p.alt}*p.M.aX_updated*{1/p.M.altX} # used if Y is binary
    or_weight <- (A==alt)*clever_coef.Y  # used if Y is continuous

    # offset term for E(Y|M,alt,X)
    offset.Y <- or_pred_alt

    model.Y <- glm( Y ~ offset(offset.Y)+1, weights=or_weight)

    eps.Y = coef(model.Y)

    # update E(Y|M,alt,X)
    or_pred_alt_updated <- or_pred_alt + eps.Y

    # calculated E(Dstar) for Y|M,alt,X
    EDstar_or <- mean(or_weight*(Y-or_pred_alt_updated)) # weight*(Y - E(Y|M,a',X))

    # update related nuisances
    or_pred_altm1_updated <- or_pred_altm1 + eps.Y
    or_pred_altm0_updated <- or_pred_altm0 + eps.Y

    # estimated psi
    ax_weight <- (A==alt)*{1/p.alt}
    kappa <- or_pred_altm1_updated*p.m1.aX_updated + or_pred_altm0_updated*p.m0.aX_updated

    estimated_psi = mean(ax_weight*kappa)


  }




  # EIF
  EIF <- or_weight*(Y-or_pred_alt_updated)+ #line 1 of the EIF equation
    m_weight*(M-p.m1.aX_updated)+ #line 2 of the EIF equation
    ax_weight*(kappa - estimated_psi) #line 3 of the EIF equation


  # confidence interval
  lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
  upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))

  tmle.out <- list(estimated_psi=estimated_psi, # estimated parameter
                   lower.ci=lower.ci, # lower bound of 95% CI
                   upper.ci=upper.ci, # upper bound of 95% CI
                   p.m1.aX=p.m1.aX_updated,  # estimated M=1|A,X
                   p.a1.X=p.a1.X,  # estimated A=1|X
                   or_pred_alt=or_pred_alt_updated, # estimated E(Y|M,A,X)
                   #
                   EIF=EIF, # EIF
                   EDstar=c(EDstar_or,EDstar_M), # EIF
                   #
                   EDstar.M_vec=EDstar.M_vec, # E(Dstar) for M|A,X over iterations
                   EDstar.Y_vec=EDstar.Y_vec, # E(Dstar) for A|X over iterations
                   #
                   eps.Y_vec=eps.Y_vec, # vector of eps2 over iterations
                   eps.M_vec=eps.M_vec, # vector of eps3 over iterations
                   iter=iter) # number of iterations for M|A,X and A|X to converge

  if (onestep==T){return(list(TMLE=tmle.out,Onestep=onestep.out))}else{return(TMLE=tmle.out)}

}
