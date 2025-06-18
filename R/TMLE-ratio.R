## The main function that execute all TMLE estimators ====
#' Function for average counterfactual outcome E(Y(a)) estimation. Implementing TMLE with all types of mediators.
#' @param a Treatment level at which the average counterfactual outcome is computed
#' @param data A Dataframe contains treatment, mediators, outcome, and measured confounders
#' @param treatment Variable name for the unvariate binary treatment
#' @param mediators Variable name for the continuous univariate mediator
#' @param outcome Variable name for the continuous univariate outcome
#' @param covariates Variable name for the measured confounders
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
#' @param truncate_lower A numeric variable, setting lower bound for the truncated propensity score. The default is 0.01.
#' @param truncate_upper A numeric variable, setting upper bound for the truncated propensity score. The default is 0.99.
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
#'       \item{\code{EDstar_M.vec}}{A vector containing the average value of the mapping of EIF in tangent space \eqn{M|A,X} over iterations. This is useful for checking the convergence behavior of the mediator density. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{EDstar_ps.vec}}{A vector containing the average value of the mapping of EIF in tangent space \eqn{A|X} over iterations. This is useful for checking the convergence behavior of the propensity score. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{eps2_vec}}{A vector containing the index for submodels of the mediator density over iterations. This is useful for checking the convergence behavior of the mediator density. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{eps3_vec}}{A vector containing the index for submodels of the propensity score over iterations. This is useful for checking the convergence behavior of the propensity score. It's expected to be close to 0 when convergence is achieved.}
#'       \item{\code{iter}}{Number of iterations where convergence is achieved for the iterative update of the mediator density and propensity score.}}
#' @import np densratio SuperLearner mvtnorm stats
#' @importFrom dplyr %>% mutate select
#' @importFrom MASS mvrnorm
#' @export
#'
#'
TMLE.ratio <- function(a,data,treatment, mediators, outcome, covariates,
                     mediator.method="bayes", np.dnorm=TRUE, superlearner=FALSE, crossfit=FALSE,K=5,
                     lib = c("SL.glm","SL.earth","SL.ranger","SL.mean"), n.iter=500, eps=TRUE, cvg.criteria=0.01,
                     formulaY="Y ~ .", formulaA="A ~ .", formulaM="M~.", linkY_binary="logit", linkA="logit", linkM_binary="logit",
                     formula_bayes="A ~ .",link_bayes="logit",
                     truncate_lower=0.01, truncate_upper=0.99,
                     estimator='onestep', boundedsubmodelY=F){


  # attach(data, warn.conflicts=FALSE)

  n <- nrow(data)

  # Variables
  A <- data[,treatment]
  M <- data[,mediators, drop = F]
  X <- data[,covariates,drop = F]
  Y <- data[,outcome]

  dat_MAX <- data.frame(Y=Y,M,A=A,X)
  dat_MaX <- data.frame(Y=Y,M,A=a,X)

  # new data sets
  data_Aa = data.frame(Y=Y,M, A=A, X)
  data_A1 = data.frame(Y=Y,M, A=1, X)
  data_A0 = data.frame(Y=Y,M, A=0, X)

  ##################################################################
  ## TMLE initialization for sequential regression based estimator
  ##################################################################


  ################################################
  ############### OUTCOME REGRESSION #############
  ################################################

  #### Fit nuisance models ####

  if (crossfit==T){ #### cross fitting + super learner #####

    if (all(Y %in% c(0,1))){ # binary outcome

      or_fit <- CV.SuperLearner(Y=Y, X=dat_MAX, family = binomial(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
      or_pred <- or_fit$SL.predict
      or_pred_a1 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A1[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
      or_pred_a0 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A0[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]

    } else { # continuous outcome

      if (boundedsubmodelY & estimator=='tmle'){
        or_fit <- CV.SuperLearner(Y=(Y-min(Y))/(max(Y)-min(Y)), X=dat_MAX, family = gaussian(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
      }else{
        or_fit <- CV.SuperLearner(Y=Y, X=dat_MAX, family = gaussian(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
      }

      or_pred <- or_fit$SL.predict
      or_pred_a1 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A1[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]
      or_pred_a0 <- unlist(lapply(1:K, function(x) predict(or_fit$AllSL[[x]], newdata=data_A0[or_fit$folds[[x]],])[[1]] %>% as.vector()))[order(unlist(lapply(1:K, function(x) or_fit$folds[[x]])))]

    }

  } else if (superlearner==T){ #### super learner #####

    if (all(Y %in% c(0,1))){ # binary outcome

      or_fit <- SuperLearner(Y=Y, X=dat_MAX, family = binomial(),SL.library = lib)
      or_pred <- predict(or_fit)[[1]] %>% as.vector()
      or_pred_a1 <- predict(or_fit, newdata=data_A1)[[1]] %>% as.vector()
      or_pred_a0 <- predict(or_fit, newdata=data_A0)[[1]] %>% as.vector()

    } else { # continuous outcome

      if(boundedsubmodelY & estimator=='tmle'){
        or_fit <- SuperLearner(Y=(Y-min(Y))/(max(Y)-min(Y)), X=dat_MAX, family = gaussian(),SL.library = lib)
      }else{
        or_fit <- SuperLearner(Y=Y, X=dat_MAX, family = gaussian(),SL.library = lib)
      }

      or_pred <- predict(or_fit)[[1]] %>% as.vector()
      or_pred_a1 <- predict(or_fit, newdata=data_A1)[[1]] %>% as.vector()
      or_pred_a0 <- predict(or_fit, newdata=data_A0)[[1]] %>% as.vector()

    }


  } else { #### simple linear regression with user input regression formula: default="Y ~ ." ####

    if (all(Y %in% c(0,1))){ # binary outcome

      or_fit <- glm(as.formula(formulaY), data=dat_MAX, family = binomial(linkY_binary))
      or_pred <- predict(or_fit, type="response")
      or_pred_a1 <- predict(or_fit, newdata=data_A1, type="response")
      or_pred_a0 <- predict(or_fit, newdata=data_A0, type="response")

    } else { # continuous outcome

      if(boundedsubmodelY & estimator=='tmle'){
        or_fit <- lm(as.formula(formulaY), data=dat_MAX %>% mutate(Y=(Y-min(Y))/(max(Y)-min(Y)) ))
      }else{
        or_fit <- lm(as.formula(formulaY), data=dat_MAX)
      }

      or_pred <- predict(or_fit)
      or_pred_a1 <- predict(or_fit, newdata=data_A1)
      or_pred_a0 <- predict(or_fit, newdata=data_A0)

    }

  }

  ## for TMLE
  # rescaling to avoid numerical issue in TMLE
  if(boundedsubmodelY & 'tmle' %in% estimator & !all(Y %in% c(0,1))){

    boundedsubmodelY.truncate_lower <- 0.001
    boundedsubmodelY.truncate_upper <- 0.999

    or_pred[or_pred>boundedsubmodelY.truncate_upper] <- boundedsubmodelY.truncate_upper; or_pred[or_pred<boundedsubmodelY.truncate_lower] <- boundedsubmodelY.truncate_lower
    or_pred_a1[or_pred_a1>boundedsubmodelY.truncate_upper] <- boundedsubmodelY.truncate_upper; or_pred_a1[or_pred_a1<boundedsubmodelY.truncate_lower] <- boundedsubmodelY.truncate_lower
    or_pred_a0[or_pred_a0>boundedsubmodelY.truncate_upper] <- boundedsubmodelY.truncate_upper; or_pred_a0[or_pred_a0<boundedsubmodelY.truncate_lower] <- boundedsubmodelY.truncate_lower


  }

  ################################################
  ############### PROPENSITY SCORE ###############
  ################################################

  #### Fit nuisance models ####

  if (crossfit==T){ #### cross fitting + super learner #####

    ps_fit <- CV.SuperLearner(Y=A, X=X, family = binomial(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)
    p.a1.X <- ps_fit$SL.predict

  } else if (superlearner==T){ #### super learner #####

    ps_fit <- SuperLearner(Y=A, X=X, family = binomial(), SL.library = lib)
    p.a1.X <- predict(ps_fit, type = "response")[[1]] %>% as.vector()  # p(A=1|X)

  } else { #### simple linear regression with user input regression formula: default="A ~ ." ####

    if (truncate_lower!=0 | truncate_upper!=1){ # under weak overlapping issue, it's more stable to run A~X via linear regression

      ps_fit <- lm(as.formula(formulaA), data=X)
      p.a1.X <- predict(ps_fit)

    } else { # without weak overlapping issue. Run A~X via logistic regression

      ps_fit <- glm(as.formula(formulaA), data=X,  family = binomial(linkA))
      p.a1.X <- predict(ps_fit, type = "response")  # p(A=1|X)

    }
  }


  # apply truncation to propensity score to deal with weak overlap. Truncated propensity score within the user specified range of [truncate_lower, truncate_upper]: default=[0,1]
  if (min(p.a1.X) < truncate_lower){

    p.a1.X[p.a1.X < truncate_lower] <- truncate_lower

    print(paste0("Truncation applied to propensity score p(A=1|X) to deal with weak overlap: truncation lower bound is ",truncate_lower))

  }

  if(max(p.a1.X) > truncate_upper){

    p.a1.X[p.a1.X > truncate_upper] <- truncate_upper

    print(paste0("Truncation applied to propensity score p(A=1|X) to deal with weak overlap: truncation upper bound is ",truncate_upper))

  }




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

  if (crossfit==T){ #### crossfitting + super learner ####

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

  }else if (superlearner==T){ #### super learner ####

    # subset X, xi to rows where A=a
    X.a <- X[A==a, , drop=F]
    xi.a <- xi[A==a]

    theta.fit <- SuperLearner(Y=xi.a, X=X.a, family = gaussian(),SL.library = lib)

    theta_x <- predict(theta.fit, newdata=X)[[1]] %>% as.vector()


  } else { #### estimate theta(X) using linear regression of xi ~ X at (A==a) ####

    theta_x <- predict(lm(xi~., data=X,weights = 1*(A==a)))

  }


  ################################################
  ########## eta=int E(Y|M,A,X)p(M|a,X)dM ########
  ################################################

  #### Fit nuisance models ####

  if (crossfit==T){ #### crossfitting + super learner ####

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

  } else if (superlearner==T){ #### super learner ####

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


  } else { #### estimate eta(X) using linear regression of E(Y|M,1/0,X) ~ X at (A==a) ####

    eta.a1 <- predict(lm(or_pred_a1~., X ,weights = 1*(A==a)))
    eta.a0 <- predict(lm(or_pred_a0~., X ,weights = 1*(A==a)))

  }



  #################################################################################################################
  ## SEQUENTIAL REGRESSION BASED ESTIMATION USING p(M|a,X)/p(M|A,X): densratio vs bayes vs dnorm
  #################################################################################################################


  if (mediator.method=="densratio"){ ################### METHOD 2A: densratio method  ###################

    # Error3: densratio method doesn't support factor variables

    if (!all(unlist(sapply(dat_MAX, class)) %in% c("integer","numeric"))){
      print("densratio method only support numeric/integer variables, try bayes method instead.")

      stop()
    }


    # if M,A,X only consists numeric/integer variables: apply density ratio estimation

    Mdata.a <- data[data$A==a,c(covariates,mediators)]
    Mdata.A <- data[data$A==(1-a),c(covariates,mediators)]

    densratio.MAX <- densratio(Mdata.a, Mdata.A)

    MAXratio <- densratio.MAX$compute_density_ratio(data[,c(covariates,mediators)]) # p(M,a,X)/p(M,1-a,X)

    M.AXratio <- {(A==a)*1+(A==(1-a))*MAXratio}/{(A==a)*1+(A==(1-a))*AXratio} # recall AXratio=p(a|X)/p(a-1|X)


  } else if (mediator.method=="bayes"){ ################### METHOD 2B: Bayes method ###################

    #### Fit nuisance models ####

    if (crossfit==T){

      bayes_fit <- CV.SuperLearner(Y=A, X=data.frame(M,X), family = binomial(), V = K, SL.library = lib, control = list(saveFitLibrary=T),saveAll = T)

      # p(A=1|X,M)
      p.a1.XM <- bayes_fit$SL.predict

    }else if (superlearner==T){

      bayes_fit <- SuperLearner(Y=A, X=data.frame(M,X), family = binomial(), SL.library = lib)

      # p(A=1|X,M)
      p.a1.XM <- predict(bayes_fit, type = "response")[[1]] %>% as.vector()  # p(A=1|X)


    } else {

      # estimate density ratio using bayes rule
      bays_fit <- glm(as.formula(formula_bayes), data=data.frame(X,M), family = binomial(link_bayes))

      # p(A=1|X,M)
      p.a1.XM <- predict(bays_fit, type = "response")

    }


    if(min(p.a1.XM)==0){

      p.a1.XM[p.a1.XM<0.01] <- 0.01

      print('Truncation applied: p(A=1|X,M) has 0 values, replaced with 0.01')

    }

    if(max(p.a1.XM)==1){

      p.a1.XM[p.a1.XM>0.99] <- 0.99

      print('Truncation applied: p(A=1|X,M) has 1 values, replaced with 0.99')

    }

    #p(A=a|X,M)
    p.a.XM <- a*p.a1.XM+(1-a)*(1-p.a1.XM)


    MAXratio <- p.a.XM/(1-p.a.XM)  #p(A=a|X,M)/p(A=a'|X,M)

    M.AXratio <- {(A==a)*1+(A==(1-a))*MAXratio}/{(A==a)*1+(A==(1-a))*AXratio} # recall AXratio=p(a|X)/p(a-1|X)


  } else if (mediator.method=="dnorm"){ # METHOD 2: sequential regression + assuming normally distributed M|A,X

    M.AXratio <- calculate_M_AX_ratio_dnorm(a,M,A,X)

    MAXratio <- M.AXratio*AXratio # [p(M|a,X)/p(M|1-a,X)]*[p(a|X)/p(a-1|X)]



  } else {

    print("Invalid method input.")

  }

  ##################################################################
  #################### One-step estimator ##########################
  ##################################################################

  if ('onestep' %in% estimator){
    ######################
    # E[Dstar] calculations
    ######################

    # E(Dstar) for E(Y|M,A,X)
    EDstar_or <- mean(M.AXratio*(Y-or_pred))

    #E(Dstar) for M=1|A,X
    EDstar_M <- mean({(A==a)/p.a.X}*(xi-theta_x))

    # E(Dstar) for A=a|X
    EDstar_ps <- mean((eta.a1-eta.a0)*(A-p.a1.X))


    ######################
    # estimate E[Y(a)]
    ######################

    # estimated psi
    estimated_psi = mean(theta_x)+EDstar_or+EDstar_M+EDstar_ps

    # EIF
    EIF <- M.AXratio*(Y-or_pred)+ #line 1 of the EIF equation
      {(A==a)/p.a.X}*(xi-theta_x)+ #line 2 of the EIF equation
      (eta.a1-eta.a0)*(A-p.a1.X)+ #line 3 of the EIF equation
      theta_x - mean(theta_x) #line 4 of the EIF equation

    # confidence interval
    lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
    upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))


    onestep.out <- list(estimated_psi=estimated_psi, # estimated parameter
                        theta_x=theta_x, # theta(x)
                        lower.ci=lower.ci, # lower bound of 95% CI
                        upper.ci=upper.ci, # upper bound of 95% CI
                        M.AXratio=M.AXratio, # estimated M|A,X
                        p.a1.X=p.a1.X,  # estimated A=1|X
                        or_pred=or_pred, # estimated E(Y|M,A,X)
                        EIF=EIF, # EIF
                        EDstar=c(EDstar_or,EDstar_M, EDstar_ps) # E(Dstar) for Y|M,A,X and M|A,X, and A|X
    )

  }

  ##############################################################################
  #################### Sequential regression based TMLE ##########################
  ##############################################################################


  if('tmle' %in% estimator){


    if(all(Y %in% c(0,1))|boundedsubmodelY){ ################# binary Y=====

      ## TMLE initialize
      # initialize eps1, eps2, eps3, eps123
      eps1 <- 1 # submodel parameter for outcome regression E(Y|M,A,X)
      eps3 <- 1 # submodel parameter for propensity score p(A|X)

      # record the values of eps2, eps3, and eps2n3 over iterations

      eps1_vec <- vector(mode = "numeric")
      eps3_vec <- vector(mode = "numeric")

      # place holder for clever coefficient 1 & 2 & 3
      clever_coef1 <- 0 # clever coefficient for outcome regression E(Y|M,A,X)
      clever_coef3 <- 0 # clever coefficient for propensity score p(A|X)

      # cumulative summation of eps*clever-coefficient
      clever_coef1_add <- 0
      clever_coef3_add <- 0

      # record average EIF in the sub-tangent space corresponding to Y|M,A,X, M|A,X and A|X over iterations
      EDstar_Y <- 1 # initial value
      EDstar_ps <- 1 # initial value

      EDstar_Y.vec <- vector(mode = "numeric")
      EDstar_ps.vec <- vector(mode = "numeric")

      # iterations
      iter <- 0

      if(boundedsubmodelY){

        cat("boundedsubmodelY is TRUE, rescaling Y to [0,1] \n")

        minY <- min(Y)
        maxY <- max(Y)

        Y <- (Y-minY)/(maxY-minY) # rescale Y to [0,1]

      }


      while (max(abs(EDstar_ps), abs(EDstar_Y))>cvg.criteria & iter<n.iter) {

        ######################
        # update p(A=1|X)
        ######################

        # clever coefficient for propensity score: eta(1,X)-eta(0,X)
        clever_coef3=eta.a1-eta.a0

        # offset term for the propensity score
        offset_ps <- qlogis(p.a1.X)

        # derive eps3
        ps_model <- glm(
          A ~ offset(offset_ps)+clever_coef3-1, family=binomial(), start=0
        )

        eps3 <- coef(ps_model)
        eps3_vec <- c(eps3_vec,eps3)


        # updated propensity score
        p.a1.X <- plogis(qlogis(p.a1.X)+eps3*(clever_coef3))


        # apply truncation to propensity score to deal with weak overlap. Truncated propensity score within the user specified range of [truncate_lower, truncate_upper]: default=[0,1]
        if (min(p.a1.X) < truncate_lower){

          p.a1.X[p.a1.X < truncate_lower] <- truncate_lower

          print(paste0("Truncation applied to propensity score p(A=1|X) to deal with weak overlap: truncation lower bound is ",truncate_lower))

        }

        if(max(p.a1.X) > truncate_upper){

          p.a1.X[p.a1.X > truncate_upper] <- truncate_upper

          print(paste0("Truncation applied to propensity score p(A=1|X) to deal with weak overlap: truncation upper bound is ",truncate_upper))

        }

        # p(A=a|X)
        p.a.X <- a*p.a1.X+(1-a)*(1-p.a1.X)


        AXratio <- p.a.X/(1-p.a.X)


        # update the p(M|a,X)/p(M|A,X) ratio to make it consistent with the updated propensity score

        if (mediator.method %in% c("densratio","bayes")){

          M.AXratio <- {(A==a)*1+(A==(1-a))*MAXratio}/{(A==a)*1+(A==(1-a))*AXratio}

        } else if (mediator.method=="dnorm"){

          MAXratio <- M.AXratio*AXratio
        }

        # E(Dstar) for A|X
        EDstar_ps <- mean(clever_coef3*(A-p.a1.X))
        EDstar_ps.vec <- c(EDstar_ps.vec,EDstar_ps)


        ######################
        # E[Y|M,A,X]
        ######################

        # one iteration
        or_model <- glm(

          Y ~ offset(qlogis(or_pred))+M.AXratio-1, family=binomial(), start=0
        )

        eps1 = coef(or_model)
        eps1_vec <- c(eps1_vec,eps1)

        # updated outcome regression
        or_pred = plogis(qlogis(or_pred)+eps1*M.AXratio)

        or_pred_a1 = plogis(qlogis(or_pred_a1)+eps1*{(a==1)*1+(a==0)*(MAXratio/AXratio)})
        or_pred_a0 = plogis(qlogis(or_pred_a0)+eps1*{(a==0)*1+(a==1)*(MAXratio/AXratio)})

        # xi=int E(Y|M,A,X)p(A|X)dA
        xi <- or_pred_a1*p.a1.X+or_pred_a0*(1-p.a1.X)

        # update eta

        if (crossfit==T){ #### crossfitting + super learner ####

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

        } else if (superlearner==T){ #### super learner ####

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


        } else { #### estimate eta(X) using linear regression of E(Y|M,1/0,X) ~ X at (A==a) ####

          eta.a1 <- predict(lm(or_pred_a1~., X ,weights = 1*(A==a)))
          eta.a0 <- predict(lm(or_pred_a0~., X ,weights = 1*(A==a)))

        }


        # E(Dstar) for Y|M,A,X
        EDstar_Y <- mean(M.AXratio*(Y-or_pred))
        EDstar_Y.vec <- c(EDstar_Y.vec,EDstar_Y)


        iter <- iter+1

      }



    }else{ ######### continuous Y ====

      ######################
      # update p(A=1|X)
      ######################

      # clever coefficient for propensity score: eta(1,X)-eta(0,X)
      clever_coef3=eta.a1-eta.a0

      # offset term for the propensity score
      offset_ps <- qlogis(p.a1.X)

      # derive eps3
      ps_model <- glm(
        A ~ offset(offset_ps)+clever_coef3-1, family=binomial(), start=0
      )

      eps3 <- coef(ps_model)


      # updated propensity score
      p.a1.X <- plogis(qlogis(p.a1.X)+eps3*(clever_coef3))


      # apply truncation to propensity score to deal with weak overlap. Truncated propensity score within the user specified range of [truncate_lower, truncate_upper]: default=[0,1]
      if (min(p.a1.X) < truncate_lower){

        p.a1.X[p.a1.X < truncate_lower] <- truncate_lower

        print(paste0("Truncation applied to propensity score p(A=1|X) to deal with weak overlap: truncation lower bound is ",truncate_lower))

      }

      if(max(p.a1.X) > truncate_upper){

        p.a1.X[p.a1.X > truncate_upper] <- truncate_upper

        print(paste0("Truncation applied to propensity score p(A=1|X) to deal with weak overlap: truncation upper bound is ",truncate_upper))

      }

      # p(A=a|X)
      p.a.X <- a*p.a1.X+(1-a)*(1-p.a1.X)


      AXratio <- p.a.X/(1-p.a.X)

      # update the p(M|a,X)/p(M|A,X) ratio to make it consistent with the updated propensity score

      if (mediator.method %in% c("densratio","bayes")){

        M.AXratio <- {(A==a)*1+(A==(1-a))*MAXratio}/{(A==a)*1+(A==(1-a))*AXratio}

      } else if (mediator.method=="dnorm"){

        MAXratio <- M.AXratio*AXratio
      }

      ######################
      # update E(Y|M,A,X)
      ######################

      # one iteration
      or_model <- glm(
        Y ~ offset(or_pred)+1, weights = M.AXratio
      )

      eps1 <- coef(or_model)

      or_pred <- or_pred+eps1

      # xi=int E(Y|M,A,X)p(A|X)dA
      xi <- (or_pred_a1+eps1)*p.a1.X+(or_pred_a0+eps1)*(1-p.a1.X)

    }






    ######################
    # update theta(X)
    #####################

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

    # derive eps2
    M_model <- lm(xi ~ offset(theta_x)+1, weights=(A==a)/p.a.X)

    eps2 <- coef(M_model)

    # update theta_x
    theta_x <- theta_x+eps2


    ######################
    # E[Dstar] calculations
    ######################

    # E(Dstar) for E(Y|M,A,X)
    EDstar_or <- mean(M.AXratio*(Y-or_pred))

    #E(Dstar) for M=1|A,X
    EDstar_M <- mean({(A==a)/p.a.X}*(xi-theta_x))

    # E(Dstar) for A=a|X
    EDstar_ps <- mean(clever_coef3*(A-p.a1.X))


    ######################
    # estimate E[Y(a)]
    ######################

    # estimated psi
    estimated_psi = mean(theta_x)

    # EIF
    EIF <- M.AXratio*(Y-or_pred)+ #line 1 of the EIF equation
      {(A==a)/p.a.X}*(xi-theta_x)+ #line 2 of the EIF equation
      clever_coef3*(A-p.a1.X)+ #line 3 of the EIF equation
      theta_x - estimated_psi #line 4 of the EIF equation

    # rescale back
    if(boundedsubmodelY){
      estimated_psi <- estimated_psi*(maxY-minY)+minY
      EIF <- EIF*(maxY-minY)
    }

    # confidence interval
    lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
    upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))

    if(boundedsubmodelY|all(Y %in% c(0,1))){
      tmle.out <- list(estimated_psi=estimated_psi, # estimated parameter
                       lower.ci=lower.ci, # lower bound of 95% CI
                       upper.ci=upper.ci, # upper bound of 95% CI
                       theta_x=theta_x, # theta(x)
                       M.AXratio=M.AXratio, # estimated M|A,X
                       p.a1.X=p.a1.X,  # estimated A=1|X
                       or_pred=or_pred, # estimated E(Y|M,A,X)
                       EIF=EIF, # EIF
                       EDstar=c(EDstar_or,EDstar_M, EDstar_ps), # E(Dstar) for Y|M,A,X and M|A,X, and A|X
                       EDstar_ps.vec=EDstar_ps.vec, # E(Dstar) for A|X over iterations
                       EDstar_or.vec=EDstar_Y.vec, # E(Dstar) for A|X over iterations
                       #
                       eps1_vec=eps1_vec, # vector of eps2 over iterations
                       eps3_vec=eps3_vec, # vector of eps3 over iterations
                       iter=iter
      )
    }else{
      tmle.out <- list(estimated_psi=estimated_psi, # estimated parameter
                       lower.ci=lower.ci, # lower bound of 95% CI
                       upper.ci=upper.ci, # upper bound of 95% CI
                       theta_x=theta_x, # theta(x)
                       M.AXratio=M.AXratio, # estimated M|A,X
                       p.a1.X=p.a1.X,  # estimated A=1|X
                       or_pred=or_pred, # estimated E(Y|M,A,X)
                       EIF=EIF, # EIF
                       EDstar=c(EDstar_or,EDstar_M, EDstar_ps) # E(Dstar) for Y|M,A,X and M|A,X, and A|X
      )
    }


  }

  if (estimator=='onestep'){return(list(Onestep=onestep.out))}
  else if (estimator=='tmle'){return(list(TMLE=tmle.out))}

}

