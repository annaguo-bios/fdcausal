## TMLE estimation for univariate binary mediator  ====
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
#' @param truncate_lower A numeric variable, setting lower bound for the truncated propensity score. The default is 0.01.
#' @param truncate_upper A numeric variable, setting upper bound for the truncated propensity score. The default is 0.99.
#' @return a list of initialization of matrices.
#' @examples
#' \donttest{
#' res <- TMLE.binary(a=1,data=continuousY_binaryM,
#' treatment="A", mediator="M", outcome="Y", covariates="X")
#' }
#' @import SuperLearner
#' @importFrom dplyr %>% mutate select
#' @return Function outputs a list containing TMLE output (and Onestep estimator output if 'onestep=T' is specified):
#' \describe{
#'       \item{\code{estimated_psi}}{The estimated parameter of interest: \eqn{E(Y^a)}}
#'       \item{\code{lower.ci}}{Lower bound of the 95\% confidence interval for \code{estimated_psi}}
#'       \item{\code{upper.ci}}{Upper bound of the 95\% confidence interval for \code{estimated_psi}}
#'       \item{\code{theta_x}}{\eqn{\int E(Y|M,A,X)p(M|A=a,X)p(A|X) dM dA}}
#'       \item{\code{p.m1.aX}}{\eqn{\int p(M=1|A=a,X)}}
#'       \item{\code{p.a1.X}}{\eqn{p(A=1|X)}}
#'       \item{\code{or_pred}}{\eqn{E(Y|M,A,X)}}
#'       \item{\code{EIF}}{The estimated efficient influence function evaluated at the observed data}
#'       \item{\code{EDstar}}{A vector of the mapping of \code{EIF} in the tangent space of \eqn{Y|M,A,X}, \eqn{M|A,X}, and \eqn{A|X}.}
#'       \item{\code{EDstar_M.vec}}{A vector containing the average value of the mapping of EIF in tangent space \eqn{M|A,X} over iterations. This is useful for checking the convergence behavior of the mediator density. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{EDstar_ps.vec}}{A vector containing the average value of the mapping of EIF in tangent space \eqn{A|X} over iterations. This is useful for checking the convergence behavior of the propensity score. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{eps2_vec}}{A vector containing the index for submodels of the mediator density over iterations. This is useful for checking the convergence behavior of the mediator density. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{eps3_vec}}{A vector containing the index for submodels of the propensity score over iterations. This is useful for checking the convergence behavior of the propensity score. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{iter}}{Number of iterations where convergence is achieved for the iterative update of the mediator density and propensity score.}}
#' @export
#'

TMLE.binary <- function(a,data,treatment, mediator, outcome, covariates,
                        onestep=TRUE, superlearner=TRUE,crossfit=FALSE,K=5,
                        lib = c("SL.glm","SL.earth","SL.ranger","SL.mean"), n.iter=500, cvg.criteria=0.01,
                        formulaY="Y ~ .", formulaA="A ~ .", formulaM="M~.", linkY_binary="logit", linkA="logit", linkM_binary="logit",
                        truncate_lower=0.01, truncate_upper=0.99){

  # attach(data, warn.conflicts=FALSE)

  n <- nrow(data)
  # Variables
  A <- data[,treatment]
  M <- data[,mediator]
  X <- data[,covariates,drop = F]
  Y <- data[,outcome]

  # new data sets
  data_Aa = data.frame(A=a,X)
  data_Aa1 = data.frame(A=1,X)
  data_Aa0 = data.frame(A=0,X)

  data_A1 = data.frame(M, A=1, X)
  data_A0 = data.frame(M, A=0, X)

  data_MAX = data.frame(M=M,A=A,X)
  data_MaX = data.frame(M,A=a,X)
  data_A1M1 = data.frame(M=1,A=1, X)
  data_A1M0 = data.frame(M=0,A=1, X)
  data_A0M1 = data.frame(M=1,A=0, X)
  data_A0M0 = data.frame(M=0,A=0, X)


  ## outcome regression
  if (crossfit==T){

    if (all(Y %in% c(0,1))){ # binary outcome

      or_fit <- CV.SuperLearner(Y=Y, X=data_MAX, family = binomial(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
      or_pred <- or_fit$SL.predict
      or_pred_a1 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A1[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
      or_pred_a0 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A0[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]

      or_pred_a1m1 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A1M1[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
      or_pred_a1m0 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A1M0[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
      or_pred_a0m1 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A0M1[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
      or_pred_a0m0 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A0M0[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]

    } else { # continuous outcome

      or_fit <- CV.SuperLearner(Y=Y, X=data_MAX, family = gaussian(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
      or_pred <- or_fit$SL.predict
      or_pred_a1 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A1[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
      or_pred_a0 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A0[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]

      or_pred_a1m1 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A1M1[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
      or_pred_a1m0 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A1M0[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
      or_pred_a0m1 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A0M1[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
      or_pred_a0m0 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A0M0[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]

    }

  } else if (superlearner==T){

    if (all(Y %in% c(0,1))){ # binary outcome

      or_fit <- SuperLearner(Y=Y, X=data_MAX, family = binomial(),SL.library = lib)
      or_pred <- predict(or_fit)[[1]] %>% as.vector()
      or_pred_a1 <- predict(or_fit, newdata=data_A1)[[1]] %>% as.vector()
      or_pred_a0 <- predict(or_fit, newdata=data_A0)[[1]] %>% as.vector()

      or_pred_a1m1 <- predict(or_fit, newdata=data_A1M1)[[1]] %>% as.vector()
      or_pred_a1m0 <- predict(or_fit, newdata=data_A1M0)[[1]] %>% as.vector()
      or_pred_a0m1 <- predict(or_fit, newdata=data_A0M1)[[1]] %>% as.vector()
      or_pred_a0m0 <- predict(or_fit, newdata=data_A0M0)[[1]] %>% as.vector()

    } else { # continuous outcome

      or_fit <- SuperLearner(Y=Y, X=data_MAX, family = gaussian(),SL.library = lib)
      or_pred <- predict(or_fit)[[1]] %>% as.vector()
      or_pred_a1 <- predict(or_fit, newdata=data_A1)[[1]] %>% as.vector()
      or_pred_a0 <- predict(or_fit, newdata=data_A0)[[1]] %>% as.vector()

      or_pred_a1m1 <- predict(or_fit, newdata=data_A1M1)[[1]] %>% as.vector()
      or_pred_a1m0 <- predict(or_fit, newdata=data_A1M0)[[1]] %>% as.vector()
      or_pred_a0m1 <- predict(or_fit, newdata=data_A0M1)[[1]] %>% as.vector()
      or_pred_a0m0 <- predict(or_fit, newdata=data_A0M0)[[1]] %>% as.vector()

    }


  } else { # simple linear regression with user input regression formula: default="Y ~ ."

    if (all(Y %in% c(0,1))){ # binary outcome

      or_fit <- glm(as.formula(formulaY), data=data_MAX, family = binomial(linkY_binary))
      or_pred <- predict(or_fit, type="response")
      or_pred_a1 <- predict(or_fit, newdata=data_A1, type="response")
      or_pred_a0 <- predict(or_fit, newdata=data_A0, type="response")

      or_pred_a1m1 <- predict(or_fit, newdata=data_A1M1, type="response")
      or_pred_a0m1 <- predict(or_fit, newdata=data_A0M1, type="response")
      or_pred_a1m0 <- predict(or_fit, newdata=data_A1M0, type="response")
      or_pred_a0m0 <- predict(or_fit, newdata=data_A0M0, type="response")

    } else { # continuous outcome

      or_fit <- lm(as.formula(formulaY), data=data_MAX)
      or_pred <- predict(or_fit)
      or_pred_a1 <- predict(or_fit, newdata=data_A1)
      or_pred_a0 <- predict(or_fit, newdata=data_A0)

      or_pred_a1m1 <- predict(or_fit, newdata=data_A1M1)
      or_pred_a0m1 <- predict(or_fit, newdata=data_A0M1)
      or_pred_a1m0 <- predict(or_fit, newdata=data_A1M0)
      or_pred_a0m0 <- predict(or_fit, newdata=data_A0M0)

    }

  }


  ## propensity score
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

  p.a.X <- a*p.a1.X + (1-a)*(1-p.a1.X)


  # p(a|X)/p(a'|X)
  # estimate the density ratio of p(a,X)/p(a',X) using regression instead of the densratio package to make it more accurate and stable
  AXratio <- p.a.X/(1-p.a.X)

  # int E(Y|M,A,X) p(A|X)dA
  xi <- or_pred_a1*p.a1.X+or_pred_a0*(1-p.a1.X)

  ################################################
  ######### theta_x=int xi p(M|a,X)dM ############
  ################################################

  #### Fit nuisance models ####

  if (crossfit==T){

    # subset X, xi to rows where A=a
    X.a <- X[A==a, , drop=F]
    xi.a <- xi[A==a]

    # fit model
    theta.fit <- CV.SuperLearner(Y=xi.a, X=X.a, family = gaussian(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)

    # prediction on rows where A=a
    theta.sub.a <- unlist(lapply(1:K, function(x) predict(theta.fit$AllSL[[x]], newdata=X.a[theta.fit$folds[[x]],,drop=F])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) theta.fit$folds[[x]])))]

    # prediction on rows where A=(1-a)
    theta.sub.a_ <- rowMeans( data.frame(lapply(1:K, function(x) predict(theta.fit$AllSL[[x]], newdata=X[A==(1-a),,drop=F])[[1]] %>% as.vector()), row.names = NULL) )


    # merge tow predictions to get prediction on all rows
    theta_x <- vector(mode="numeric", length=n)

    theta_x[A==a] <- theta.sub.a; theta_x[A==(1-a)] <- theta.sub.a_

  }else if (superlearner==T){

    # subset X, xi to rows where A=a
    X.a <- X[A==a, , drop=F]
    xi.a <- xi[A==a]

    theta.fit <- SuperLearner(Y=xi.a, X=X.a, family = gaussian(),SL.library = lib)

    theta_x <- predict(theta.fit, newdata=X)[[1]] %>% as.vector()


  } else { # estimate theta(X) using linear regression

    theta_x <- predict(lm(xi~., data=X,weights = 1*(A==a)))

  }


  ################################################
  ########## eta=int E(Y|M,A,X)p(M|a,X)dM ########
  ################################################


  #### Fit nuisance models ####

  if (crossfit==T){

    # subset X, or_pred_a1, and or_pred_a0 to rows where A=a
    X.a <- X[A==a, , drop=F]
    or_pred_a1.a <- or_pred_a1[A==a]
    or_pred_a0.a <- or_pred_a0[A==a]

    # fit model
    eta.a1.fit <- CV.SuperLearner(Y= or_pred_a1.a, X=X.a, family = gaussian(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
    eta.a0.fit <- CV.SuperLearner(Y= or_pred_a0.a, X=X.a, family = gaussian(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)

    ## prediction using all observations

    # prediction on rows where A=a
    eta.a1.sub.a <- unlist(lapply(1:K, function(x) predict(eta.a1.fit$AllSL[[x]], newdata=X.a[eta.a1.fit$folds[[x]],,drop=F])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) eta.a1.fit$folds[[x]])))]
    eta.a0.sub.a <- unlist(lapply(1:K, function(x) predict(eta.a0.fit$AllSL[[x]], newdata=X.a[eta.a0.fit$folds[[x]],,drop=F])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) eta.a0.fit$folds[[x]])))]

    # prediction on rows where A=(1-a)
    eta.a1.sub.a_ <- rowMeans( data.frame(lapply(1:K, function(x) predict(eta.a1.fit$AllSL[[x]], newdata=X[A==(1-a),,drop=F])[[1]] %>% as.vector()), row.names = NULL) )
    eta.a0.sub.a_ <- rowMeans( data.frame(lapply(1:K, function(x) predict(eta.a0.fit$AllSL[[x]], newdata=X[A==(1-a),,drop=F])[[1]] %>% as.vector()), row.names = NULL) )

    # merge tow predictions to get prediction on all rows
    eta.a1 <- vector(mode="numeric", length=n)
    eta.a0 <- vector(mode="numeric", length=n)

    eta.a1[A==a] <- eta.a1.sub.a; eta.a1[A==(1-a)] <- eta.a1.sub.a_
    eta.a0[A==a] <- eta.a0.sub.a; eta.a0[A==(1-a)] <- eta.a0.sub.a_

  } else if (superlearner==T){

    # subset X, or_pred_a1, and or_pred_a0 to rows where A=a
    X.a <- X[A==a, , drop=F]
    or_pred_a1.a <- or_pred_a1[A==a]
    or_pred_a0.a <- or_pred_a0[A==a]

    # fit the model
    eta.a1.fit <- SuperLearner(Y=or_pred_a1.a, X=X.a, family = gaussian(),SL.library = lib)
    eta.a0.fit <- SuperLearner(Y=or_pred_a0.a, X=X.a, family = gaussian(),SL.library = lib)

    # prediction using all observations
    eta.a1 <- predict(eta.a1.fit, newdata=X)[[1]] %>% as.vector()
    eta.a0 <- predict(eta.a0.fit, newdata=X)[[1]] %>% as.vector()


  } else {

    eta.a1 <- predict(lm(or_pred_a1~., X ,weights = 1*(A==a)))
    eta.a0 <- predict(lm(or_pred_a0~., X ,weights = 1*(A==a)))

  }




  ## M|A,X
  if (crossfit==T){

    M_fit <- CV.SuperLearner(Y=M, X=data.frame(A=A,X), family = binomial(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
    p.m1.AX <- M_fit$SL.predict
    p.m1.aX <- unlist(lapply(1:K, function(x) predict(M_fit$AllSL[[x]], newdata=data_Aa[M_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) M_fit$folds[[x]])))]
    p.m1.a1X <- unlist(lapply(1:K, function(x) predict(M_fit$AllSL[[x]], newdata=data_Aa1[M_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) M_fit$folds[[x]])))]
    p.m1.a0X <- unlist(lapply(1:K, function(x) predict(M_fit$AllSL[[x]], newdata=data_Aa0[M_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) M_fit$folds[[x]])))]

  }else if (superlearner==T){
    M_fit <- SuperLearner(Y=M, X=data.frame(A=A,X), family = binomial(), SL.library = lib)
    p.m1.AX <- predict(M_fit, type = "response")[[1]] %>% as.vector()  # p(M=1|A,X)
    p.m1.aX <- predict(M_fit, newdata=data_Aa, type = "response")[[1]] %>% as.vector() # p(M=1|A=a,X)
    p.m1.a1X <- predict(M_fit, newdata=data_Aa1, type = "response")[[1]] %>% as.vector() # p(M=1|A=1,X)
    p.m1.a0X <- predict(M_fit, newdata=data_Aa0, type = "response")[[1]] %>% as.vector() # p(M=1|A=0,X)

  }else{

    M_fit <- glm(as.formula(formulaM), data=data.frame(A=A,X), family = binomial(linkM_binary))
    p.m1.AX <- predict(M_fit, type = "response")  # p(M=1|A,X)
    p.m1.aX <- predict(M_fit, newdata=data_Aa, type = "response") # p(M=1|A=a,X)
    p.m1.a1X <- predict(M_fit, newdata=data_Aa1, type = "response") # p(M=1|A=a,X)
    p.m1.a0X <- predict(M_fit, newdata=data_Aa0, type = "response") # p(M=1|A=a,X)

  }


  ############################################################### One step estimator ###############################################################
  ## onestep estimator
  if (onestep==T){
    ######################
    # E[Dstar] calculations
    ######################

    # E(Dstar) for E(Y|M,A,X)
    or_weight <- (M*p.m1.aX+(1-M)*(1-p.m1.aX))/(M*p.m1.AX+(1-M)*(1-p.m1.AX))
    EDstar_or <- mean(or_weight*(Y-or_pred))

    #E(Dstar) for M=1|A,X
    m_weight <- (A==a)*{1/p.a.X}*
      (or_pred_a1m1*p.a1.X+ # E(Y|M=1, A=1, X)*p(A=1|X)
         or_pred_a0m1*(1-p.a1.X)- # E(Y|M=1, A=0, X)*p(A=0|X)
         or_pred_a1m0*p.a1.X- # E(Y|M=0, A=1, X)*p(A=1|X)
         or_pred_a0m0*(1-p.a1.X)) # E(Y|M=0, A=0, X)*p(A=0|X)


    EDstar_M <- mean(m_weight*(M-p.m1.aX))

    # E(Dstar) for A=a|X
    ps_weight <- or_pred_a1m1*p.m1.aX+ # E(Y|M=1,A=1,X)*p(M=1|A=a,X)
      or_pred_a1m0*(1-p.m1.aX)- # E(Y|M=0,A=1,X)*p(M=0|A=a,X)
      or_pred_a0m1*p.m1.aX- # E(Y|M=1,A=0,X)*p(M=1|A=a,X)
      or_pred_a0m0*(1-p.m1.aX) # E(Y|M=0,A=0,X)*p(M=0|A=a,X)

    EDstar_ps <- mean(ps_weight*(A-p.a1.X))


    ######################
    # estimate E[Y(a)]
    ######################

    # theta(X) = int E(Y|M,A,X)p(M|A=a,X)p(A|X)dM dA
    theta_x <- (or_pred_a1m1)*p.m1.aX*p.a1.X+ # M=1,A=1
      (or_pred_a0m1)*p.m1.aX*(1-p.a1.X)+ # M=1, A=0
      (or_pred_a1m0)*(1-p.m1.aX)*p.a1.X+ # M=0, A=1
      (or_pred_a0m0)*(1-p.m1.aX)*(1-p.a1.X) # M=0, A=0

    # estimated psi
    estimated_psi = mean(theta_x)+EDstar_or+EDstar_M+EDstar_ps

    # EIF
    EIF <- or_weight*(Y-or_pred)+ #line 1 of the EIF equation
      m_weight*(M-p.m1.aX)+
      ps_weight*(A-p.a1.X)+
      theta_x - estimated_psi #line 4 of the EIF equation

    # confidence interval
    lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
    upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))

    onestep.out <- list(estimated_psi=estimated_psi, # estimated parameter
                        lower.ci=lower.ci, # lower bound of 95% CI
                        upper.ci=upper.ci, # upper bound of 95% CI
                        theta_x=theta_x, # theta(x)
                        p.m1.aX=p.m1.aX,  # estimated M=1|A,X
                        p.a1.X=p.a1.X,  # estimated A=1|X
                        or_pred=or_pred, # estimated E(Y|M,A,X)
                        #
                        EIF=EIF, # EIF
                        EDstar=c(EDstar_or,EDstar_M, EDstar_ps)) # E(Dstar) for Y|M,A,X and M|A,X, and A|X

    print("One step estimator done")

  }



  ############################################################### End of One step estimator ###############################################################




  ############################################################### TMLE estimator ###############################################################

  if (all(Y %in% c(0,1))){ # binary Y


    ## TMLE initialize
    # initialize eps1, eps2, eps3, eps123
    eps1 <- 1 # submodel parameter for outcome regression E(Y|M,A,X)
    eps2 <- 1 # submodel parameter for mediator density p(M|A,X)
    eps3 <- 1 # submodel parameter for propensity score p(A|X)

    # record the values of eps2, eps3, and eps2n3 over iterations

    eps1_vec <- vector(mode = "numeric")
    eps2_vec <- vector(mode = "numeric")
    eps3_vec <- vector(mode = "numeric")

    # place holder for clever coefficient 1 & 2 & 3
    clever_coef1 <- 0 # clever coefficient for outcome regression E(Y|M,A,X)
    clever_coef2 <- 0 # clever coefficient for mediator density p(M|A,X)
    clever_coef3 <- 0 # clever coefficient for propensity score p(A|X)

    # cumulative summation of eps*clever-coefficient
    clever_coef1_add <- 0
    clever_coef2_add <- 0
    clever_coef3_add <- 0

    # record average EIF in the sub-tangent space corresponding to Y|M,A,X, M|A,X and A|X over iterations
    EDstar_Y <- 1 # initial value
    EDstar_M <- 1 # initial value
    EDstar_ps <- 1 # initial value

    EDstar_Y.vec <- vector(mode = "numeric")
    EDstar_M.vec <- vector(mode = "numeric")
    EDstar_ps.vec <- vector(mode = "numeric")

    # initialize E(Y|M,A,X), p(A=1|X), p(M=1|A=a,X)
    or_pred_updated <- or_pred
    # or_pred_a1m1_updated <- or_pred_a1m1
    # or_pred_a1m0_updated <- or_pred_a1m0
    # or_pred_a0m1_updated <- or_pred_a0m1
    # or_pred_a0m0_updated <- or_pred_a0m0

    p.a1.X_updated <- p.a1.X
    p.m1.aX_updated <- p.m1.aX

    # iterations
    iter <- 0

    while (max(abs(EDstar_M), abs(EDstar_ps), abs(EDstar_Y))>cvg.criteria & iter<n.iter) {

      ######################
      # update p(M|A=a,X)
      ######################

      ## clever coefficient for M|A=a,X
      p.a.X_updated <- a*p.a1.X_updated+(1-a)*(1-p.a1.X_updated) # p(A=a|X)

      clever_coef2 = {1/p.a.X_updated}*
        (or_pred_a1m1*p.a1.X_updated+ # E(Y|M=1, A=1, X)*p(A=1|X)
           or_pred_a0m1*(1-p.a1.X_updated)- # E(Y|M=1, A=0, X)*p(A=0|X)
           or_pred_a1m0*p.a1.X_updated- # E(Y|M=0, A=1, X)*p(A=1|X)
           or_pred_a0m0*(1-p.a1.X_updated)) # E(Y|M=0, A=0, X)*p(A=0|X)

      # offset term for M|A=a,X
      offset_M <- qlogis(p.m1.aX) + clever_coef2_add


      # derive eps2
      M_model <- glm(
        M ~ offset(offset_M)+clever_coef2-1, weights=(A==a)*1, family=binomial()
      )

      eps2 <-  coef(M_model)
      eps2_vec <- c(eps2_vec,eps2)

      # print(paste0('iter: ',iter,'. M model done. eps=',eps2))

      # update cumulative summation of eps2*clever coefficient
      clever_coef2_add = clever_coef2_add + eps2*(clever_coef2)

      # update p(M=1|A=a,X)
      p.m1.aX_updated <- plogis(qlogis(p.m1.aX)+clever_coef2_add)

      # E(Dstar) for M|A,X
      EDstar_M <- mean((A==a)*clever_coef2*(M-p.m1.aX_updated))
      EDstar_M.vec <- c(EDstar_M.vec,EDstar_M)

      ######################
      # update p(A=1|X)
      ######################

      # clever coefficient for propensity score
      clever_coef3=or_pred_a1m1*p.m1.aX_updated+ # E(Y|M=1,A=1,X)*p(M=1|A=a,X)
        or_pred_a1m0*(1-p.m1.aX_updated)- # E(Y|M=0,A=1,X)*p(M=0|A=a,X)
        or_pred_a0m1*p.m1.aX_updated- # E(Y|M=1,A=0,X)*p(M=1|A=a,X)
        or_pred_a0m0*(1-p.m1.aX_updated) # E(Y|M=0,A=0,X)*p(M=0|A=a,X)

      # offset term for the propensity score
      # the same as: qlogis(p.a1.X_updated)
      offset_ps <- qlogis(p.a1.X) + clever_coef3_add

      # derive eps3
      ps_model <- glm(
        A ~ offset(offset_ps)+clever_coef3-1, family=binomial(), start=0
      )

      eps3 <- coef(ps_model)
      eps3_vec <- c(eps3_vec,eps3)

      # print(paste0('iter: ',iter,'. A model done. eps=',eps3))

      # update cumulative summation of eps3*clever coefficient
      clever_coef3_add = clever_coef3_add + eps3*(clever_coef3)

      # updated propensity score
      p.a1.X_updated <- plogis(qlogis(p.a1.X)+clever_coef3_add)

      # E(Dstar) for A|X
      EDstar_ps <- mean(clever_coef3*(A-p.a1.X_updated))
      EDstar_ps.vec <- c(EDstar_ps.vec,EDstar_ps)


      ######################
      # update E(Y|M,A,X)
      ######################

      # clever coefficient for propensity score
      p.m1.AX_updated <- (A==a)*p.m1.aX_updated+(1-(A==a))*p.m1.AX

      or_weight <- (M*p.m1.aX_updated+(1-M)*(1-p.m1.aX_updated))/(M*p.m1.AX_updated+(1-M)*(1-p.m1.AX_updated)) # p(M|a,X)/p(M|A,X)
      clever_coef1 <- or_weight

      # offset term for outcome regression
      offset_Y <- qlogis(or_pred) + clever_coef1_add

      # one iteration
      or_model <- glm(
        Y ~ offset(or_pred)+clever_coef1-1, family=binomial(), start=0
      )

      eps1 = coef(or_model)
      eps1_vec <- c(eps1_vec,eps1)

      # print(paste0('iter: ',iter,'. Y model done. eps=',eps1))

      # update cumulative summation of eps1*clever coefficient
      clever_coef1_add = clever_coef1_add + eps1*(clever_coef1)

      # updated outcome regression
      or_pred_updated = plogis(qlogis(or_pred)+clever_coef3_add)

      or_pred_a1m1 = plogis(qlogis(or_pred_a1m1)+eps1*{(a==1)*1+(a==0)*p.m1.aX_updated/p.m1.a1X})
      or_pred_a0m1 = plogis(qlogis(or_pred_a0m1)+eps1*{(a==0)*1+(a==1)*p.m1.aX_updated/p.m1.a0X})
      or_pred_a1m0 = plogis(qlogis(or_pred_a1m0)+eps1*{(a==1)*1+(a==0)*(1-p.m1.aX_updated)/(1-p.m1.a1X)})
      or_pred_a0m0 = plogis(qlogis(or_pred_a0m0)+eps1*{(a==0)*1+(a==1)*(1-p.m1.aX_updated)/(1-p.m1.a0X)})

      # E(Dstar) for Y|M,A,X
      EDstar_Y <- mean(clever_coef1*(Y-or_pred_updated))
      EDstar_Y.vec <- c(EDstar_Y.vec,EDstar_Y)


      iter <- iter+1
    }


    ######################
    # estimate E[Y(a)]
    ######################

    # theta(X) = int E(Y|M,A,X)p(M|A=a,X)p(A|X)dM dA
    theta_x <- or_pred_a1m1*p.m1.aX_updated*p.a1.X_updated+ # M=1,A=1
      or_pred_a0m1*p.m1.aX_updated*(1-p.a1.X_updated)+ # M=1, A=0
      or_pred_a1m0*(1-p.m1.aX_updated)*p.a1.X_updated+ # M=0, A=1
      or_pred_a0m0*(1-p.m1.aX_updated)*(1-p.a1.X_updated) # M=0, A=0

    # estimated psi
    estimated_psi = mean(theta_x)


  } else { # continuous Y


    ## TMLE initialize
    # initialize eps2, eps3, eps2n3
    eps2 <- 1
    eps3 <- 1

    # record the values of eps2, eps3, and eps2n3 over iterations
    eps2_vec <- vector(mode = "numeric")
    eps3_vec <- vector(mode = "numeric")

    # place holder for clever coefficient 2 & 3
    clever_coef2 <- 0
    clever_coef3 <- 0

    # cumulative summation of eps*clever-coefficient
    clever_coef2_add <- 0
    clever_coef3_add <- 0

    # record E(Dstar) for M|A,X and propensity score over iterations
    EDstar_M <- 1 # initial value
    EDstar_ps <- 1 # initial value
    EDstar_M.vec <- vector(mode = "numeric")
    EDstar_ps.vec <- vector(mode = "numeric")

    # initial p(A=1|X), p(M=1|A=a,X)
    p.a1.X_updated <- p.a1.X
    p.m1.aX_updated <- p.m1.aX

    # iterations
    iter <- 0

    while (max(abs(EDstar_M),abs(EDstar_ps))>cvg.criteria & iter<n.iter) {

      ######################
      # update p(M|A=a,X)
      ######################

      ## clever coefficient for M|A=a,X
      p.a.X_updated <- a*p.a1.X_updated+(1-a)*(1-p.a1.X_updated) # p(A=a|X)

      clever_coef2 = {1/p.a.X_updated}*
        (or_pred_a1m1*p.a1.X_updated+ # E(Y|M=1, A=1, X)*p(A=1|X)
           or_pred_a0m1*(1-p.a1.X_updated)- # E(Y|M=1, A=0, X)*p(A=0|X)
           or_pred_a1m0*p.a1.X_updated- # E(Y|M=0, A=1, X)*p(A=1|X)
           or_pred_a0m0*(1-p.a1.X_updated)) # E(Y|M=0, A=0, X)*p(A=0|X)

      # offset term for M|A=a,X
      offset_M <- qlogis(p.m1.aX) + clever_coef2_add


      # derive eps2
      M_model <- glm(
        M ~ offset(offset_M)+clever_coef2-1, weights=(A==a)*1, family=binomial()
      )


      eps2 <-  coef(M_model)
      eps2_vec <- c(eps2_vec,eps2)

      # update cumulative summation of eps2*clever coefficient
      clever_coef2_add = clever_coef2_add + eps2*(clever_coef2)

      # update p(M=1|A=a,X)
      p.m1.aX_updated <- plogis(qlogis(p.m1.aX)+clever_coef2_add)

      # E(Dstar) for M|A,X
      EDstar_M <- mean((A==a)*clever_coef2*(M-p.m1.aX_updated))
      EDstar_M.vec <- c(EDstar_M.vec,EDstar_M)

      ######################
      # update p(A=1|X)
      ######################

      # clever coefficient for propensity score
      clever_coef3=or_pred_a1m1*p.m1.aX_updated+ # E(Y|M=1,A=1,X)*p(M=1|A=a,X)
        or_pred_a1m0*(1-p.m1.aX_updated)- # E(Y|M=0,A=1,X)*p(M=0|A=a,X)
        or_pred_a0m1*p.m1.aX_updated- # E(Y|M=1,A=0,X)*p(M=1|A=a,X)
        or_pred_a0m0*(1-p.m1.aX_updated) # E(Y|M=0,A=0,X)*p(M=0|A=a,X)

      # offset term for the propensity score
      # the same as: qlogis(p.a1.X_updated)
      offset_ps <- qlogis(p.a1.X) + clever_coef3_add


      # derive eps3
      ps_model <- glm(
        A ~ offset(offset_ps)+clever_coef3-1, family=binomial(), start=0
      )


      eps3 <- coef(ps_model)
      eps3_vec <- c(eps3_vec,eps3)

      # update cumulative summation of eps3*clever coefficient
      clever_coef3_add = clever_coef3_add + eps3*(clever_coef3)

      # updated propensity score
      p.a1.X_updated <- plogis(qlogis(p.a1.X)+clever_coef3_add)

      # E(Dstar)-iteration1
      EDstar_ps <- mean(clever_coef3*(A-p.a1.X_updated))
      EDstar_ps.vec <- c(EDstar_ps.vec,EDstar_ps)

      iter <- iter+1
    }

    ######################
    # update E[Y|A,M,X]
    ######################

    p.m1.AX_updated <- (A==a)*p.m1.aX_updated+(1-(A==a))*p.m1.AX
    or_weight <- (M*p.m1.aX_updated+(1-M)*(1-p.m1.aX_updated))/(M*p.m1.AX_updated+(1-M)*(1-p.m1.AX_updated))

    # one iteration
    or_model <- glm(
      Y ~ offset(or_pred)+1, weights = or_weight
    )


    eps1 = coef(or_model)

    # updated outcome regression
    or_pred_updated = or_pred + eps1

    ######################
    # E[Dstar] calculations
    ######################

    # E(Dstar) for E(Y|M,A,X)
    EDstar_Y <- mean(or_weight*(Y-or_pred_updated))

    #E(Dstar) for M=1|A,X
    # since p.a.X_updated might have changed, need to recalculate clever_coef2
    p.a.X_updated <- a*p.a1.X_updated+(1-a)*(1-p.a1.X_updated) # p(A=a|X)

    clever_coef2 <- {1/p.a.X_updated}*
      (or_pred_a1m1*p.a1.X_updated+ # M=1,A=1
         or_pred_a0m1*(1-p.a1.X_updated)- # M=1, A=0
         or_pred_a1m0*p.a1.X_updated- # M=0, A=1
         or_pred_a0m0*(1-p.a1.X_updated)) # M=0, A=0

    EDstar_M <- mean((A==a)*clever_coef2*(M-p.m1.aX_updated))

    # E(Dstar) for A=a|X
    EDstar_ps <- mean(clever_coef3*(A-p.a1.X_updated))


    ######################
    # estimate E[Y(a)]
    ######################

    # theta(X) = int E(Y|M,A,X)p(M|A=a,X)p(A|X)dM dA
    theta_x <- (or_pred_a1m1+eps1)*p.m1.aX_updated*p.a1.X_updated+ # M=1,A=1
      (or_pred_a0m1+eps1)*p.m1.aX_updated*(1-p.a1.X_updated)+ # M=1, A=0
      (or_pred_a1m0+eps1)*(1-p.m1.aX_updated)*p.a1.X_updated+ # M=0, A=1
      (or_pred_a0m0+eps1)*(1-p.m1.aX_updated)*(1-p.a1.X_updated) # M=0, A=0

    # estimated psi
    estimated_psi = mean(theta_x)

  }

  # EIF
  EIF <- or_weight*(Y-or_pred_updated)+ #line 1 of the EIF equation
    I(A==a)*clever_coef2*(M-p.m1.aX_updated)+ #line 2 of the EIF equation
    clever_coef3*(A-p.a1.X_updated)+ #line 3 of the EIF equation
    theta_x - mean(theta_x) #line 4 of the EIF equation


  # confidence interval
  lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
  upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))

  tmle.out <- list(estimated_psi=estimated_psi, # estimated parameter
                   lower.ci=lower.ci, # lower bound of 95% CI
                   upper.ci=upper.ci, # upper bound of 95% CI
                   theta_x=theta_x, # theta(x)
                   p.m1.aX=p.m1.aX_updated,  # estimated M=1|A,X
                   p.a1.X=p.a1.X_updated,  # estimated A=1|X
                   or_pred=or_pred_updated, # estimated E(Y|M,A,X)
                   EIF=EIF, # EIF
                   #
                   EDstar=c(EDstar_or,EDstar_M, EDstar_ps), # E(Dstar) for Y|M,A,X and M|A,X, and A|X
                   EDstar_M.vec=EDstar_M.vec, # E(Dstar) for M|A,X over iterations
                   EDstar_ps.vec=EDstar_ps.vec, # E(Dstar) for A|X over iterations
                   #
                   eps2_vec=eps2_vec, # vector of eps2 over iterations
                   eps3_vec=eps3_vec, # vector of eps3 over iterations
                   iter=iter) # number of iterations for M|A,X and A|X to converge

  if (onestep==T){return(list(TMLE=tmle.out,Onestep=onestep.out))}
  else if (onestep==F){return(TMLE=tmle.out)}
}
