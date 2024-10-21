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
#'  onestep=TRUE, linkA="identity")
#'
#' # E(Y(1)|0) estimation. For continuous outcome Y and binary mediator M
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <- ATT.TMLE.all(a=1,data=continuousY_binaryM,
#' treatment="A", mediators="M", outcome="Y", covariates="X",
#'  onestep=TRUE, linkA="identity")
#'
#' # E(Y(1)|0) estimation. For continuous outcome Y and bivariate mediator M.
#' # Using 'densratio' method for mediator density estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <- ATT.TMLE.all(a=1,data=continuousY_bivariateM,
#' treatment="A", mediators=c("M.1","M.2"), outcome="Y", covariates="X",
#'  onestep=TRUE, linkA="identity", mediator.method="densratio")
#'
#' # E(Y(1)|0) estimation. For continuous outcome Y and bivariate mediator M.
#' # Using 'bayes' method for mediator density estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <- ATT.TMLE.all(a=1,data=continuousY_bivariateM,
#' treatment="A", mediators=c("M.1","M.2"), outcome="Y", covariates="X",
#'  onestep=TRUE, linkA="identity", mediator.method="bayes")
#'
#' # E(Y(1)|0) estimation. For continuous outcome Y and bivariate mediator M.
#' # Using 'dnorm' method for mediator density estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <- ATT.TMLE.all(a=1,data=continuousY_bivariateM,
#' treatment="A", mediators=c("M.1","M.2"), outcome="Y", covariates="X",
#'  onestep=TRUE, linkA="identity", mediator.method="dnorm")
#' }
#' @import np densratio SuperLearner mvtnorm stats
#' @importFrom dplyr %>% mutate select
#' @importFrom MASS mvrnorm
#' @export
#'
#'
ATT.TMLE.all <- function(a,data,treatment, mediators, outcome, covariates,
                     onestep=TRUE, mediator.method="bayes", np.dnorm=TRUE, superlearner=FALSE, crossfit=FALSE,K=5,
                     lib = c("SL.glm","SL.earth","SL.ranger","SL.mean"), n.iter=500, eps=TRUE, cvg.criteria=0.01,
                     formulaY="Y ~ .", formulaA="A ~ .", formulaM="M~.", linkY_binary="logit", linkA="logit", linkM_binary="logit",
                     formula_bayes="A ~ .",link_bayes="logit",
                     truncate_lower=0, truncate_upper=1, print.message=TRUE){


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

    # execute the TMLE.binary function
    tmle.binary.res <- ATT.TMLE.binary(a=a,data=data,treatment=treatment, mediator=mediators, outcome=outcome, covariates=covariates, onestep=onestep ,
                                   superlearner=superlearner,crossfit=crossfit, K=K, lib=lib,n.iter=n.iter, cvg.criteria=cvg.criteria,
                                   formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
                                   truncate_lower=truncate_lower, truncate_upper=truncate_upper)
    return(tmle.binary.res)

  } else if (mediator.method=="np"){ # METHOD 1- NP

    ## Error 0
    # np method only allow univariate mediator
    if (length(mediators)>1){ stop("np method only allow univariate mediator")} else if (is.numeric(unlist(M))) { # univariate continuous mediator

      ## Error1

      # binary outcome require iterative update among propensity score, mediator density, and outcome regression. Results can be very unstable.
      if (all(Y %in% c(0,1))){ stop("np method under binary outcome is not stable. Try densratio method, bayes method, or dnorm method instead")}

      ## Error2

      if (crossfit==T){stop("Due to its computational burden. np method is not currently supported for crossfit. Set crossfit=F instead.")}

      if (superlearner==T){stop("Due to its computational burden. np method is not currently supported for superlearner. Set superlearner=F instead.")}


      if (np.dnorm){

        # execute the TMLE.np.dnorm function
        tmle.np.out <- ATT.TMLE.np.dnorm(a=a,data=data,treatment=treatment, mediator=mediators, outcome=outcome, covariates=covariates,
                                     onestep=onestep, n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
                                     formulaY=formulaY, formulaA=formulaA, linkA=linkA,
                                     truncate_lower=truncate_lower, truncate_upper=truncate_upper)} else {

        # execute the TMLE.np function
        tmle.np.out <- ATT.TMLE.np(a=a,data=data,treatment=treatment, mediator=mediators, outcome=outcome, covariates=covariates,
                               onestep=onestep, n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
                               formulaY=formulaY, formulaA=formulaA, linkA=linkA,
                               truncate_lower=truncate_lower, truncate_upper=truncate_upper)}

      return(tmle.np.out)

    }


  }


  ##################################################################
  ## TMLE initialization for sequential regression based estimator
  ##################################################################

  ##############################
  ############### p(A) #########
  ##############################

  p.a <- mean(A==a)
  p.alt <- mean(A==alt)

  ################################################
  ############### OUTCOME REGRESSION #############
  ################################################

  Y.binary <- all(Y %in% c(0,1) ) # whether Y is binary
  Y.family <- if (Y.binary){binomial()}else{gaussian()}


  if (crossfit==T){ # use cross fitting

    or_fit <- CV.SuperLearner(Y=Y, X=data.frame(M,A=A,X), family = Y.family, V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
    or_pred <- or_fit$SL.predict # E(Y|M,A,X)
    or_pred_alt <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data.frame(M, A=alt,X)[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))] # E(Y|M,a',X)


  } else if (superlearner==T){ # use superlearner

    or_fit <- SuperLearner(Y=Y, X=data.frame(M,A=A,X), family = Y.family, SL.library = lib)
    or_pred <- predict(or_fit)[[1]] %>% as.vector() # E(Y|M,A,X)
    or_pred_alt <- predict(or_fit, newdata=data.frame(M,A=alt,X))[[1]] %>% as.vector() # E(Y|M,a',X)


  } else{ # use simple regression

    Y.family <- if (all(Y %in% c(0,1) )){binomial(linkY_binary)}else{gaussian()}

    or_fit <- glm(as.formula(formulaY), data=data.frame(M,A=A,X), family = Y.family)
    or_pred <- predict(or_fit, type="response") # E(Y|M,A,X)
    or_pred_alt <- predict(or_fit, newdata=data.frame(M,A=alt,X), type="response") # E(Y|M,a',X)

  }


  ################################################
  ############### PROPENSITY SCORE ###############
  ################################################

  if (crossfit==T){

    ps_fit <- CV.SuperLearner(Y=A, X=X, family = binomial(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
    p.a1.X <- ps_fit$SL.predict

  } else if (superlearner==T){

    ps_fit <- SuperLearner(Y=A, X=X, family = binomial(), SL.library = lib)
    p.a1.X <- predict(ps_fit, type = "response")[[1]] %>% as.vector()  # p(A=1|X)

  } else {

    if (truncate_lower!=0 | truncate_upper!=1){ # under weak overlapping issue, it's more stable to run A~X via linear regression

      ps_fit <- lm(as.formula(formulaA), data=X)
      p.a1.X <- predict(ps_fit)

      print("Truncation performed.")

      if(print.message){print("ps fit is"); print(summary(ps_fit))}
      if(print.message){print("X is"); print(summary(X))}

    } else { # without weak overlapping issue. Run A~X via logistic regression

      ps_fit <- glm(as.formula(formulaA), data=X,  family = binomial(linkA))
      p.a1.X <- predict(ps_fit, type = "response")  # p(A=1|X)
    }

  }

  # apply truncation to propensity score to deal with weak overlap. Truncated propensity score within the user specified range of [truncate_lower, truncate_upper]: default=[0,1]
  p.a1.X[p.a1.X < truncate_lower] <- truncate_lower
  p.a1.X[p.a1.X > truncate_upper] <- truncate_upper

  if(print.message==T){print('truncate lower is'); print(truncate_lower)}
  if(print.message==T){print('truncate upper is'); print(truncate_upper)}
  if(print.message==T){print('p.a1.X is '); print(summary(p.a1.X))}

  p.a.X <- a*p.a1.X + (1-a)*(1-p.a1.X) # p(a|X)
  p.alt.X <- 1-p.a.X # p(a'|X)

  # estimate the density ratio of p(a,X)/p(a',X) using regression instead of the densratio package to make it more accurate and stable
  # p(a'|X)/p(a|X)
  ratio.A.X <- p.alt.X/p.a.X

  if(print.message==T){print('p.a.X is '); print(summary(p.a.X))}
  if(print.message==T){print('p.alt.X is '); print(summary(p.alt.X))}
  if(print.message==T){print('ratio.A.X is '); print(summary(ratio.A.X))}

  ##############################################################
  ######### kappa=int E(Y|M,a',X) p(M|a,X)dM ###################
  ##############################################################

  if (crossfit==T){ # use cross fitting

    kappa_fit <- CV.SuperLearner(Y=or_pred_alt, X=data.frame(A=A,X), family = Y.family, V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)

    # kappa = int E(Y|M,a',X) p(M|a,X)dM
    kappa <- unlist(lapply(1:K, function(x) predict(kappa_fit$AllSL[[x]], newdata=data.frame(A=a, X)[kappa_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) kappa_fit$folds[[x]])))]


  }else if (superlearner==T){ # use superlearner

    kappa_fit <- SuperLearner(Y=or_pred_alt, X=data.frame(A=A,X), family = Y.family, SL.library = lib)

    # kappa = int E(Y|M,a',X) p(M|a,X)dM
    kappa <- predict(kappa_fit, newdata=data.frame(A=a,X))[[1]] %>% as.vector()

  }else{ # use simple regression

    kappa_fit <- glm(or_pred_alt ~ ., data=data.frame(A=A,X), family = Y.family)

    # kappa = int E(Y|M,a',X) p(M|a,X)dM
    kappa <- predict(kappa_fit, newdata=data.frame(A=a,X), type="response")

  }


  #################################################################################################################
  ## SEQUENTIAL REGRESSION BASED ESTIMATION USING p(M|a,X)/p(M|A,X): densratio vs bayes vs dnorm
  #################################################################################################################

  ################### METHOD 2A: densratio method  ###################
  if (mediator.method=="densratio"){

    # Error3: densratio method doesn't support factor variables

    if (!all(unlist(sapply(data.frame(M,A=A,X), class)) %in% c("integer","numeric"))){stop("densratio method only support numeric/integer variables, try bayes method instead.")}


    # if M,A,X only consists numeric/integer variables: apply density ratio estimation

    Mdata.a <- data[data$A==a,c(covariates,mediators)]
    Mdata.alt <- data[data$A==alt,c(covariates,mediators)]

    densratio.MAX <- densratio(Mdata.a, Mdata.alt)

    MAXratio <- densratio.MAX$compute_density_ratio(data[,c(covariates,mediators)]) # p(M,a,X)/p(M,a',X)

    M.AXratio <- MAXratio*ratio.A.X # recall ratio.A.X=p(a'|X)/p(a|X)

    # METHOD 2: sequential regression + assuming normally distributed M|A,X

  } else if (mediator.method=="dnorm"){

    M.AXratio <- calculate_M_AX_ratio_dnorm(a,M,A,X,ATT=T)

    MAXratio <- M.AXratio/ratio.A.X # [p(M|a,X)/p(M|a',X)]/[p(a'|X)/p(a|X)]


    ################### METHOD 2B: Bayes method ###################
  } else if (mediator.method=="bayes"){

    if (crossfit==T){

      bayes_fit <- CV.SuperLearner(Y=A, X=data.frame(M,X), family = binomial(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)

      # p(A=1|X,M)
      p.a1.XM <- bayes_fit$SL.predict

      #p(A=a|X,M)
      p.a.XM <- a*p.a1.XM+(1-a)*(1-p.a1.XM)

    } else if (superlearner==T){

      bayes_fit <- SuperLearner(Y=A, X=data.frame(M,X), family = binomial(), SL.library = lib)

      # p(A=1|X,M)
      p.a1.XM <- predict(bayes_fit, type = "response")[[1]] %>% as.vector()  # p(A=1|X)

      #p(A=a|X,M)
      p.a.XM <- a*p.a1.XM+(1-a)*(1-p.a1.XM)

    } else {

      # estimate density ratio using bayes rule
      bays_fit <- glm(as.formula(formula_bayes), data=data.frame(X,M), family = binomial(link_bayes))

      # p(A=1|X,M)
      p.a1.XM <- predict(bays_fit, type = "response")

      #p(A=a|X,M)
      p.a.XM <- a*p.a1.XM+(1-a)*(1-p.a1.XM)

    }


    MAXratio <- p.a.XM/(1-p.a.XM)  #p(A=a|X,M)/p(A=a'|X,M)

    M.AXratio <- MAXratio*ratio.A.X # recall ratio.A.X=p(a'|X)/p(a|X)

  } else {

    stop("Invalid mediator method input.")

  }



  ##################################################################
  #################### One-step estimator ##########################
  ##################################################################

  if (onestep==T){

    ######################
    # E[Dstar] calculations
    ######################

    # E(Dstar) for Y|M,A,X
    or_weight <- (A==alt)*{1/p.alt}*M.AXratio # I(A=a')/p(a') * p(M|A=a,X)/p(M|A=a',X)
    EDstar_or <- mean(or_weight*(Y-or_pred_alt)) # weight*(Y - E(Y|M,a',X))

    # E(Dstar) for M|A,X
    EDstar_M <- mean( (A==a)*{1/p.alt}*ratio.A.X*(or_pred_alt - kappa) ) # weight*(E(Y|M=1,A=a',X) - E{ E(Y|A=a',M,X) | a, X})

    # at tangent space of A,X
    ax_weight <- (A==alt)*{1/p.alt}


    ######################
    # estimate E[Y(a)]
    ######################

    # estimated psi
    estimated_psi = EDstar_or + EDstar_M + mean(ax_weight*kappa) # estimated psi

    # EIF
    EIF <- or_weight*(Y-or_pred_alt)+ #line 1 of the EIF equation
      (A==a)*{1/p.alt}*ratio.A.X*(or_pred_alt - kappa)+ #line 2 of the EIF equation
      ax_weight*(kappa - estimated_psi) #line 3 of the EIF equation

    # confidence interval
    lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
    upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))

    onestep.out <- list(estimated_psi=estimated_psi, # estimated parameter
                        lower.ci=lower.ci, # lower bound of 95% CI
                        upper.ci=upper.ci, # upper bound of 95% CI
                        p.a1.X=p.a1.X,  # estimated A=1|X
                        or_pred_alt=or_pred_alt, # estimated E(Y|M,A,X)
                        #
                        EIF=EIF, # EIF
                        EDstar=c(EDstar_or,EDstar_M)) # E(Dstar) for Y|M,A,X and M|A,X, and A,X

    if(print.message==T){print('onestep done')}

  }

  ##############################################################################
  #################### Sequential regression based TMLE ##########################
  ##############################################################################



  ######################
  # Update E(Y|M,alt,X)
  ######################

  ## clever coefficient for E(Y|M,alt,X)
  clever_coef.Y <- {1/p.alt}*M.AXratio # used if Y is binary
  or_weight <- (A==alt)*clever_coef.Y  # used if Y is continuous

  if(print.message==T){print('p.alt is ');print(summary(p.alt))}
  if(print.message==T){print('M.AXratio is '); print(summary(M.AXratio))}
  if(print.message==T){print(summary(clever_coef.Y))}



  # offset term for E(Y|M,alt,X)
  offset.Y <- if(Y.binary){qlogis(or_pred_alt)}else{or_pred_alt}

  if(print.message==T){print(summary(or_weight))}
  if(print.message==T){print(summary(offset.Y))}

  if(Y.binary){model.Y <- glm(Y ~ offset(offset.Y) + clever_coef.Y -1, weights=(A==alt)*1, family=binomial() )}else{model.Y <- glm( Y ~ offset(offset.Y)+1, weights=or_weight, family = gaussian())}

  if(print.message==T){print(summary(model.Y))}


  # the optimizer
  eps.Y = coef(model.Y)

  if(print.message==T){print(paste0('eps.Y=',eps.Y))}

  # update E(Y|M,alt,X)
  or_pred_alt_updated <- if(Y.binary){plogis(offset.Y + eps.Y*clever_coef.Y)}else{or_pred_alt + eps.Y}

  # calculated E(Dstar) for Y|M,alt,X
  EDstar_or <- mean(or_weight*(Y-or_pred_alt_updated)) # weight*(Y - E(Y|M,a',X))

  if(print.message==T){print('Y updated')}

  #---------------------------------
  # update related nuisances - kappa
  #---------------------------------

  if (crossfit==T){ # use cross fitting

    kappa_fit <- CV.SuperLearner(Y=or_pred_alt_updated, X=data.frame(A=A,X), family = Y.family, V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)

    # kappa = int E(Y|M,a',X) p(M|a,X)dM
    kappa <- unlist(lapply(1:K, function(x) predict(kappa_fit$AllSL[[x]], newdata=data.frame(A=a, X)[kappa_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) kappa_fit$folds[[x]])))]


  } else if (superlearner==T){ # use superlearner

    kappa_fit <- SuperLearner(Y=or_pred_alt_updated, X=data.frame(A=A,X), family = Y.family, SL.library = lib)

    # kappa = int E(Y|M,a',X) p(M|a,X)dM
    kappa <- predict(kappa_fit, newdata=data.frame(A=a,X))[[1]] %>% as.vector()

  } else { # use simple regression

    kappa_fit <- glm(or_pred_alt_updated ~ ., data=data.frame(A=A,X), family = Y.family)

    # kappa = int E(Y|M,a',X) p(M|a,X)dM
    kappa <- predict(kappa_fit, newdata=data.frame(A=a,X), type="response")

  }

  if (print.message==T){print('kappa updated upon Y updated')}

  ######################
  # update kappa= int E(Y|M,a',X) p(M|a,X)dM
  ######################

  ## clever coefficient for M|A=a,X
  clever_coef.M = {1/p.alt}*ratio.A.X # used if Y is binary
  m_weight = (A==a)*clever_coef.M # used if Y is continuous

  # offset term for M|A=a,X
  offset.M <- if(Y.binary){qlogis(kappa)}else{kappa}

  # derive eps.M
  model.M <- if(Y.binary){glm( or_pred_alt_updated ~ offset(offset.M)+clever_coef.M-1, weights=(A==a)*1, family=binomial() )}else{glm(or_pred_alt_updated ~ offset(offset.M)+1, weights=m_weight, family=gaussian())}

  # the optimizer
  eps.M <-  coef(model.M)

  # update kappa
  kappa_updated <- if(Y.binary){plogis(offset.M + eps.M*clever_coef.M)}else{kappa + eps.M}

  # calculated E(Dstar) for M|A=a,X
  EDstar_M <- mean(m_weight*(or_pred_alt_updated - kappa_updated))

  if(print.message==T){print('kappa updated')}

  # estimated psi
  ax_weight <- (A==alt)*{1/p.alt}

  estimated_psi = mean(ax_weight*kappa_updated)





  # EIF
  EIF <- or_weight*(Y-or_pred_alt_updated)+ #line 1 of the EIF equation
    m_weight*(or_pred_alt_updated - kappa_updated)+ #line 2 of the EIF equation
    ax_weight*(kappa_updated - estimated_psi) #line 3 of the EIF equation


  # confidence interval
  lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
  upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))

  tmle.out <- list(estimated_psi=estimated_psi, # estimated parameter
                   lower.ci=lower.ci, # lower bound of 95% CI
                   upper.ci=upper.ci, # upper bound of 95% CI
                   M.AXratio=M.AXratio, # estimated p(M|a,X)/p(M|a',X)
                   p.a1.X=p.a1.X,  # estimated A=1|X
                   or_pred_alt=or_pred_alt_updated, # estimated E(Y|M,A,X)
                   #
                   EIF=EIF, # EIF
                   EDstar=c(EDstar_or,EDstar_M) # EIF
                   ) # number of iterations for M|A,X and A|X to converge

  if (onestep==T){return(list(TMLE=tmle.out,Onestep=onestep.out))}else{return(TMLE=tmle.out)}


}

