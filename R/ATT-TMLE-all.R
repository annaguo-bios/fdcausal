## The main function that execute all TMLE estimators for ATT ====
#' Function for average counterfactual outcome E(Y(a)|1-a) estimation. Implementing TMLE with all types of mediators.
#' @param a Treatment level at which the average counterfactual outcome is computed
#' @param data A Dataframe contains treatment, mediators, outcome, and measured confounders
#' @param treatment Variable name for the unvariate binary treatment
#' @param mediators Variable name for the continuous univariate mediator
#' @param outcome Variable name for the continuous univariate outcome
#' @param covariates Variable name for the measured confounders
#' @param onestep A logical indicator determines whether one-step estimation is executed. When 'onestep=T', the one-step estimation result is provided. Conversely, if 'onestep=F', the result is withheld.
#' @param mediator.method When M is univariate binary, regression is adopted for estimating \eqn{p(M|A,X)}. Otherwise, four methods for mediator density estimation is provided, namely "bayes", "densratio", "dnorm", and "np".
#' The "bayes" method estimates the density ratio \eqn{p(M|A,X)/p(M|a,X)} by rewriting it as \eqn{(p(a|M,X)/p(A|M,X))/(p(a|X)/p(A|X))}, where p(A|M,X) is then estimated via regression. TMLE estimators using 'bayes' method avoids updating the mediator density p(M|A,X).
#' The "densratio" method estimates \eqn{p(M|A,X)/p(M|a,X)} by rewriting it as \eqn{(p(M,a,X)/p(M,A,X))/(p(a|X)/p(A|X))}, where \eqn{p(M,A,X)/p(M,a,X)} is estimated with the \link[densratio]{densratio} function,
#' and \eqn{p(a|X)/p(A|X)} is estimated via regression. "densratio" method is only applicable if all mediators are numeric or integer valued. TMLE estimators using 'densratio' method avoids updating the mediator density p(M|A,X).
#' The "dnorm" method estimates the mediator density ratio \eqn{p(M|A,X)/p(M|a,X)} assuming p(M|A,X) follows a normal distribution. TMLE estimators using 'dnorm' method avoids updating the mediator density p(M|A,X).
#' TMLE estimators using 'np' method directly update the mediator density p(M|A,X). When np.dnorm=F, the "np" method estimates \eqn{p(M|a,X)/p(M|A,X)} by estimating \eqn{p(M|A,X)} with the \link[np]{npcdens} function.
#' When np.dnorm=T, the "np" method estimates \eqn{p(M|a,X)/p(M|A,X)} assuming normal distribution. Due to the computational burden, "np" method is only available for univariate continuous mediator.
#' @param superlearner A logical indicator determines whether SuperLearner via the \link[SuperLearner]{SuperLearner} function is adopted for estimating the outcome regression, mediator density, and the propensity score.
#' SuperLearner is the preferred option in cases where complex relationships among variables exist, potentially leading to model misspecification issues when using simple linear models.
#' @param crossfit A logical indicator determines whether SuperLearner+Cross-fitting is adopted for estimating the outcome regression, mediator density, and the propensity score.
#' @param K A integer indicating the number of folds for cross-fitting, the default is 5.
#' @param lib Library of algorithms for SuperLearner.
#' @param n.iter The maximum number of iterations performed when iteratively updating the mediator density and propensity score.
#' @param eps A logical indicator determines the stopping criteria used when iteratively updating the mediator density and propensity score. The default is 'eps=T'.
#' When 'eps=T', \eqn{\sqrt{\epsilon_2^2+\epsilon_3^2}} is used, where \eqn{\epsilon_2} and \eqn{\epsilon_3} are the index of the sub-models for mediator density and propensity score. When 'eps=F', \eqn{max(|\Phi_M|,|\Phi_A|)} is used,
#' where \eqn{\Phi_M} and \eqn{\Phi_A} are the mapping of efficient influcence function (EIF) into the tangent space of \eqn{M|A,X} and \eqn{A|X}. In general, adoption of 'eps=F' results in better convergence while it takes longer time.
#' Conversely, adoption of 'eps=T' usually requires less time but can result in algorithm divergence.
#' @param cvg.criteria A numerical value representing the convergence criteria when iteratively updating the mediator density and propensity score.
#'  The default value is 0.01, meaning update stops when stopping criteria < 0.01. The stopping criteria is chosen by the eps.
#' @param formulaY Regression formula for the outcome regression of Y on M, A, X. The default is 'Y ~ 1+ M + A + X'.
#' @param linkY_binary The link function used for the logistic regression of Y on M, A, X when Y is binary. The default is the 'logit' link.
#' @param formulaA Regression formula for the propensity score regression of A on X. The default is 'A ~ 1 + X'.
#' @param linkA The link function used for the logistic regression of A on X. The default is the 'logit' link.
#' @param formulaM Regression formula for the mediator density regression of M on A and X. The default is 'M ~ 1 + A + X'. This parameter is only needed when M is a univariate binary mediator.
#' @param linkM_binary The link function used for the logistic regression of M on A and X. The default is the 'logit' link. This parameter is only needed when M is a univariate binary mediator.
#' @param formula_bayes Regression formula for the regression of A on M and X. The default is 'A ~ 1 + M + X'. This parameter is only needed when mediator.method="bayes".
#' @param link_bayes The link function used for the logistic regression of A on M and X. The default is 'logit' link. This parameter is only needed when mediator.method="bayes".
#' @param truncate_lower A numeric variable, setting lower bound for the truncated propensity score. The default is 0.
#' @param truncate_upper A numeric variable, setting upper bound for the truncated propensity score. The default is 1.
#' @param np.dnorm A logic variable. If np.dnorm=T, p(M|A,X) is directly estimated assuming normal distribution. If np.dnorm=F, p(M|A,X) is directly estimated using the \link[np]{npcdens} function.
#' @param estimator A character string indicating which estimator is to be used. The options are "onestep" and "tmle".
#' @param boundedsubmodelY An indicator for whether the bounded submodel is used for targeting the outcome regression when Z is discrete. The default is FALSE.
#' @return Function outputs a list containing TMLE results (and Onestep results if 'onestep=T' is specified):
#' \describe{
#'       \item{\code{estimated_psi}}{The estimated parameter of interest: \eqn{E(Y^a)}}
#'       \item{\code{lower.ci}}{Lower bound of the 95\% confidence interval for \code{estimated_psi}}
#'       \item{\code{upper.ci}}{Upper bound of the 95\% confidence interval for \code{estimated_psi}}
#'       \item{\code{theta_x}}{\eqn{\int E(Y|M,A,X)p(M|A=a,X)p(A|X) dM dA}}
#'       \item{\code{p.m1.aX}}{\eqn{p(M=1|A=a,X)}}
#'       \item{\code{p.a1.X}}{\eqn{p(A=1|X)}}
#'       \item{\code{or_pred}}{\eqn{E(Y|M,A,X)}}
#'       \item{\code{EIF}}{The estimated efficient influence function evaluated at the observed data}
#'       \item{\code{EDstar}}{A vector of the mapping of \code{EIF} in the tangent space of \eqn{Y|M,A,X}; \eqn{M|A,X}; and \eqn{A|X}.}
#'       \item{\code{EDstar.M_vec}}{A vector containing the average value of the mapping of EIF in tangent space \eqn{M|A,X} over iterations. This is useful for checking the convergence behavior of the mediator density. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{EDstar.Y_vec}}{A vector containing the average value of the mapping of EIF in tangent space \eqn{A|X} over iterations. This is useful for checking the convergence behavior of the propensity score. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{eps.M_vec}}{A vector containing the index for submodels of the mediator density over iterations. This is useful for checking the convergence behavior of the mediator density. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{eps.Y_vec}}{A vector containing the index for submodels of the propensity score over iterations. This is useful for checking the convergence behavior of the propensity score. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{iter}}{Number of iterations where convergence is achieved for the iterative update of the mediator density and propensity score.}}
#' @examples
#' \donttest{
#' # E(Y(1)|0) estimation. For binary outcome Y and binary mediator M.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <- ATT.TMLE.all(a=1,data=binaryY_binaryM,
#' treatment="A", mediators="M", outcome="Y", covariates="X",
#'   linkA="identity")
#'
#' # E(Y(1)|0) estimation. For continuous outcome Y and binary mediator M
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <- ATT.TMLE.all(a=1,data=continuousY_binaryM,
#' treatment="A", mediators="M", outcome="Y", covariates="X",
#'   linkA="identity")
#'
#' # E(Y(1)|0) estimation. For continuous outcome Y and bivariate mediator M.
#' # Using 'densratio' method for mediator density estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <- ATT.TMLE.all(a=1,data=continuousY_bivariateM,
#' treatment="A", mediators=c("M.1","M.2"), outcome="Y", covariates="X",
#'   linkA="identity", mediator.method="densratio")
#'
#' # E(Y(1)|0) estimation. For continuous outcome Y and bivariate mediator M.
#' # Using 'bayes' method for mediator density estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <- ATT.TMLE.all(a=1,data=continuousY_bivariateM,
#' treatment="A", mediators=c("M.1","M.2"), outcome="Y", covariates="X",
#'   linkA="identity", mediator.method="bayes")
#'
#' # E(Y(1)|0) estimation. For continuous outcome Y and bivariate mediator M.
#' # Using 'dnorm' method for mediator density estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <- ATT.TMLE.all(a=1,data=continuousY_bivariateM,
#' treatment="A", mediators=c("M.1","M.2"), outcome="Y", covariates="X",
#'   linkA="identity", mediator.method="dnorm")
#' }
#' @import np densratio SuperLearner mvtnorm stats
#' @importFrom dplyr %>% mutate select
#' @importFrom MASS mvrnorm
#' @export
#'
#'
ATT.TMLE.all <- function(a,data,treatment, mediators, outcome, covariates,
                      mediator.method="bayes", np.dnorm=TRUE, superlearner=FALSE, crossfit=FALSE,K=5,
                     lib = c("SL.glm","SL.earth","SL.ranger","SL.mean"), n.iter=500, eps=TRUE, cvg.criteria=0.01,
                     formulaY="Y ~ .", formulaA="A ~ .", formulaM="M~.", linkY_binary="logit", linkA="logit", linkM_binary="logit",
                     formula_bayes="A ~ .",link_bayes="logit",
                     truncate_lower=0, truncate_upper=1, print.message=TRUE,
                     estimator='onestep', boundedsubmodelY=F){


  # attach(data, warn.conflicts=FALSE)

  n <- nrow(data)

  alt <- 1-a # alternative treatment level

  # Variables
  A <- data[,treatment]
  M <- data[,mediators, drop = F]
  X <- data[,covariates,drop = F]
  Y <- data[,outcome]

  dat_MAX <- data.frame(M,A=A,X)
  dat_MaX <- data.frame(M,A=a,X)

  # new data sets
  data_Aa = data.frame(M, A=A, X)
  data_A1 = data.frame(M, A=1, X)
  data_A0 = data.frame(M, A=0, X)

  #####################################################################################################################
  # DIRECT DENSITY p(M|a,X) ESTIMATION: binary M ~ regression based direct estimation VS continuous M ~ np
  #####################################################################################################################

  if (length(mediators)==1 & all(unlist(M) %in% 0:1)){ # METHOD 1- BINARY: univariate binary mediator

    if(length(estimator)==2){

      # execute the TMLE.binary function
      tmle.binary.res1 <- ATT.TMLE.binary(a=a,data=data,treatment=treatment, mediator=mediators, outcome=outcome, covariates=covariates,
                                         superlearner=superlearner,crossfit=crossfit, K=K, lib=lib,n.iter=n.iter, cvg.criteria=cvg.criteria,
                                         formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
                                         truncate_lower=truncate_lower, truncate_upper=truncate_upper, estimator=estimator[1], boundedsubmodelY=boundedsubmodelY)

      tmle.binary.res2 <- ATT.TMLE.binary(a=a,data=data,treatment=treatment, mediator=mediators, outcome=outcome, covariates=covariates,
                                         superlearner=superlearner,crossfit=crossfit, K=K, lib=lib,n.iter=n.iter, cvg.criteria=cvg.criteria,
                                         formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
                                         truncate_lower=truncate_lower, truncate_upper=truncate_upper, estimator=estimator[2], boundedsubmodelY=boundedsubmodelY)

      tmle.binary.res <- c(tmle.binary.res1,tmle.binary.res2)

    }else{


      # execute the TMLE.binary function
      tmle.binary.res <- ATT.TMLE.binary(a=a,data=data,treatment=treatment, mediator=mediators, outcome=outcome, covariates=covariates,
                                         superlearner=superlearner,crossfit=crossfit, K=K, lib=lib,n.iter=n.iter, cvg.criteria=cvg.criteria,
                                         formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
                                         truncate_lower=truncate_lower, truncate_upper=truncate_upper, estimator=estimator, boundedsubmodelY=boundedsubmodelY)

    }


    return(tmle.binary.res)

  } else if (mediator.method=="np"){ # METHOD 1- NP

    ## Error 0
    # np method only allow univariate mediator
    if (length(mediators)>1){ stop("np method only allow univariate mediator")} else if (is.numeric(unlist(M))) { # univariate continuous mediator

      ## Error1

      # binary outcome require iterative update among propensity score, mediator density, and outcome regression. Results can be very unstable.
      if (all(Y %in% c(0,1))|boundedsubmodelY){ stop("np method under binary outcome is not stable. Try densratio method, bayes method, or dnorm method instead")}

      ## Error2

      if (crossfit==T){stop("Due to its computational burden. np method is not currently supported for crossfit. Set crossfit=F instead.")}

      if (superlearner==T){stop("Due to its computational burden. np method is not currently supported for superlearner. Set superlearner=F instead.")}


      if (np.dnorm){

        # execute the TMLE.np.dnorm function
        tmle.np.out <- ATT.TMLE.np.dnorm(a=a,data=data,treatment=treatment, mediator=mediators, outcome=outcome, covariates=covariates,
                                     n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
                                     formulaY=formulaY, formulaA=formulaA, linkA=linkA,
                                     truncate_lower=truncate_lower, truncate_upper=truncate_upper, estimator=estimator)} else {

        # execute the TMLE.np function
        tmle.np.out <- ATT.TMLE.np(a=a,data=data,treatment=treatment, mediator=mediators, outcome=outcome, covariates=covariates,
                               n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
                               formulaY=formulaY, formulaA=formulaA, linkA=linkA,
                               truncate_lower=truncate_lower, truncate_upper=truncate_upper, estimator=estimator)}

      return(tmle.np.out)

    }


  }


  #####################################################################################################################
  # Mediator ratio based methods
  #####################################################################################################################

  if(length(estimator)==2){

    cat('Implementing',estimator[1],'\n')

    tmle.ratio.out1 <- ATT.TMLE.ratio(a = a,
                                  data = data,
                                  treatment = treatment,
                                  mediators = mediators,
                                  outcome = outcome,
                                  covariates = covariates,
                                  mediator.method = mediator.method,
                                  np.dnorm = np.dnorm,
                                  superlearner = superlearner,
                                  crossfit = crossfit,
                                  K = K,
                                  lib = lib,
                                  n.iter = n.iter,
                                  eps = eps,
                                  cvg.criteria = cvg.criteria,
                                  formulaY = formulaY,
                                  formulaA = formulaA,
                                  formulaM = formulaM,
                                  linkY_binary = linkY_binary,
                                  linkA = linkA,
                                  linkM_binary = linkM_binary,
                                  formula_bayes = formula_bayes,
                                  link_bayes = link_bayes,
                                  truncate_lower = truncate_lower,
                                  truncate_upper = truncate_upper,
                                  estimator = estimator[1],
                                  boundedsubmodelY = boundedsubmodelY)

    cat('Implementing',estimator[2],'\n')
    tmle.ratio.out2 <- ATT.TMLE.ratio(a = a,
                                  data = data,
                                  treatment = treatment,
                                  mediators = mediators,
                                  outcome = outcome,
                                  covariates = covariates,
                                  mediator.method = mediator.method,
                                  np.dnorm = np.dnorm,
                                  superlearner = superlearner,
                                  crossfit = crossfit,
                                  K = K,
                                  lib = lib,
                                  n.iter = n.iter,
                                  eps = eps,
                                  cvg.criteria = cvg.criteria,
                                  formulaY = formulaY,
                                  formulaA = formulaA,
                                  formulaM = formulaM,
                                  linkY_binary = linkY_binary,
                                  linkA = linkA,
                                  linkM_binary = linkM_binary,
                                  formula_bayes = formula_bayes,
                                  link_bayes = link_bayes,
                                  truncate_lower = truncate_lower,
                                  truncate_upper = truncate_upper,
                                  estimator = estimator[2],
                                  boundedsubmodelY = boundedsubmodelY)

    tmle.ratio.out <- c(tmle.ratio.out1, tmle.ratio.out2)


  }else{

    tmle.ratio.out <- ATT.TMLE.ratio(a = a,
                                 data = data,
                                 treatment = treatment,
                                 mediators = mediators,
                                 outcome = outcome,
                                 covariates = covariates,
                                 mediator.method = mediator.method,
                                 np.dnorm = np.dnorm,
                                 superlearner = superlearner,
                                 crossfit = crossfit,
                                 K = K,
                                 lib = lib,
                                 n.iter = n.iter,
                                 eps = eps,
                                 cvg.criteria = cvg.criteria,
                                 formulaY = formulaY,
                                 formulaA = formulaA,
                                 formulaM = formulaM,
                                 linkY_binary = linkY_binary,
                                 linkA = linkA,
                                 linkM_binary = linkM_binary,
                                 formula_bayes = formula_bayes,
                                 link_bayes = link_bayes,
                                 truncate_lower = truncate_lower,
                                 truncate_upper = truncate_upper,
                                 estimator = estimator,
                                 boundedsubmodelY = boundedsubmodelY)

  }

  return(tmle.ratio.out)


}

