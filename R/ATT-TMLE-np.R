## TMLE np for ATT ====
#' Function for implementing TMLE with the mediator density \eqn{p(M|A,X)} estimated with nonparametric kernel density method: https://cran.r-project.org/web/packages/np/index.html.
#' @param a Treatment level at which the average counterfactual outcome is computed
#' @param data A Dataframe contains treatment, mediators, outcome, and measured confounders
#' @param treatment Variable name for the unvariate binary treatment
#' @param mediator Variable name for the continuous univariate mediator
#' @param outcome Variable name for the continuous univariate outcome
#' @param covariates Variable name for the measured confounders
#' @param onestep A logical indicator determines whether one-step estimation is executed. When 'onestep=T', the one-step estimation result is provided. Conversely, if 'onestep=F', the result is withheld.
#' @param n.iter The maximum number of iterations performed when iteratively updating the mediator density and propensity score.
#' @param eps A logical indicator determines the stopping criteria used when iteratively updating the mediator density and propensity score. The default is 'eps=T'.
#' When 'eps=T', \eqn{\sqrt{\epsilon_2^2+\epsilon_3^2}} is used. When 'eps=F', \eqn{max(|\Phi_M|,|\Phi_A|} is used. In general, adoption of 'eps=F' results in better convergence while it takes longer time.
#' Conversely, adoption of 'eps=T' usually requires less time but can result in algorithm divergence.
#' @param cvg.criteria A numerical value representing the convergence criteria when iteratively updating the mediator density and propensity score. The default value is 0.01, meaning update stops when stopping criteria < 0.01. The stopping criteria is chosen by the eps.
#' @param formulaY Regression formula for the outcome regression of Y on M, A, X. The default is 'Y ~ 1+ M + A + X'.
#' @param formulaA Regression formula for the propensity score regression of A on X. The default is 'A ~ 1 + X'.
#' @param linkA The link function used for the logistic regression of A on X. The default is the 'logit' link.
#' @param truncate_lower A numeric variable, setting lower bound for the truncated propensity score. The default is 0.
#' @param truncate_upper A numeric variable, setting upper bound for the truncated propensity score. The default is 1.
#' @param estimator A character string indicating which estimator is to be used. The default is c('onestep','tmle'), which returns both one-step estimator and TMLE. Other options are "onestep" and "tmle".
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
#' @examples
#' \donttest{
#' ATT.TMLE.np(1,data=continuousY_continuousM,
#' treatment="A", mediator="M", outcome="Y", covariates="X")
#' }
#' @import np
#' @importFrom dplyr %>% mutate select
#' @export

ATT.TMLE.np <- function(a,data,treatment="A", mediator="M", outcome="Y", covariates=c("X1","X2","X3"),
                    n.iter=500, eps=T, cvg.criteria=0.01,
                    formulaY="Y ~ .", formulaA="A ~ .", linkA="logit",
                    truncate_lower=0, truncate_upper=1, estimator='onestep'){

  # attach(data, warn.conflicts=FALSE)

  # Variables
  A <- data[,treatment]
  M <- data[,mediator]
  X <- data[,covariates, drop=F]
  Y <- data[,outcome]
  dat_MAX <- data.frame(M=M,A=A,X)
  dat_MaX <- data.frame(M=M,A=a,X)

  n <- nrow(data) # number of observations

  alt <- 1-a # alternative treatment level

  # covariates is of length 1, rename the covariates to "X"
  if (length(covariates)==1){ covariates = "X" }

  ######################
  ## Nuisance estimation
  ######################

  ## p(A) ====
  p.a <- mean(A==a)
  p.alt <- mean(A==alt)

  ## outcome regression ====

  or_fit <- lm(as.formula(formulaY), data=dat_MAX)
  or_pred_alt <- predict(or_fit, newdata=data.frame(M=M, A=alt,X))

  # updated function for outcome regression
  f.or <- function(M,A,X,eps.Y=0){predict(or_fit, newdata = data.frame(M=M,A=A,X,row.names = NULL))+eps.Y}

  ## M|A,X ====
  # Methods on density estimation:
  # np: https://cran.r-project.org/web/packages/np/np.pdf

  # Construct the formula using paste0
  formula_str <- paste0(paste(mediator, collapse = " + "), " ~ ", paste(c(treatment,covariates), collapse = " + "))

  # Create the formula
  formula <- as.formula(formula_str)
  bw <- npcdensbw(formula=formula, data=dat_MAX)
  M_fit <- npcdens(bws=bw)

  # prediction
  p.M.aX <- predict(M_fit, newdata=data.frame(M=M, A=a, X)) # p(M|A=a,X)
  p.M.altX <- predict(M_fit, newdata=data.frame(M=M, A=alt, X)) # p(M|A=a',X)


  ## propensity score ====

  if (truncate_lower!=0 | truncate_upper!=1){ # under weak overlapping issue, it's more stable to run A~X via linear regression

    ps_fit <- lm(as.formula(formulaA), data=X)
    p.a1.X <- predict(ps_fit)

  } else { # without weak overlapping issue. Run A~X via logistic regression

    ps_fit <- glm(as.formula(formulaA), data=X,  family = binomial(linkA))
    p.a1.X <- predict(ps_fit, type = "response")  # p(A=1|X)

  }


  p.a1.X[p.a1.X < truncate_lower] <- truncate_lower
  p.a1.X[p.a1.X > truncate_upper] <- truncate_upper

  p.a.X <- a*p.a1.X + (1-a)*(1-p.a1.X) # p(a|X)
  p.alt.X <- 1-p.a.X # p(a'|X)

  # estimate the density ratio of p(a,X)/p(a',X) using regression instead of the densratio package to make it more accurate and stable
  # p(a'|X)/p(a|X)
  ratio.A.X <- p.alt.X/p.a.X

  # updated density function for M|A,X

  ## np
  f.m <- function(M,A,X, eps.M=0, clever_coef.M=0){

    # predict p(M|A,X)
    newdata <- data.frame(M=M,A=A,X, row.names = NULL); names(newdata) <- c("M","A", covariates) # new data frame for prediction

    p.M.AX <- predict(M_fit, newdata=newdata) # original p(M|A,X)

    p.M.AX <- p.M.AX*{1+eps.M*clever_coef.M} # post-target p(M|A,X)

    return(p.M.AX)
  }

  # integrand of kappa: int E(Y|M,A=a',X)p(M|A=a,X) dM
  integrand <- function(m,a,x, eps.Y=0 ,eps.M=0,clever_coef.M=0){

    or_pred_alt <- f.or(m,alt,x,eps.Y) # A=a'

    M_pred <- f.m(m,a,x,eps.M,clever_coef.M) # p(M|A=a,X)

    integrand <- or_pred_alt*M_pred # int E(Y|M,A=a',X) p(M|A=a,X) dM

    return(integrand)
  }


  ## initial estimates
  # kappa: int E(Y|m,a',X) p(m|a,X) dM
  kappa <- sapply(1:n, function(i){integrate(integrand,lower = range(M)[1], upper = range(M)[2], a=a,x=X[i,], eps.Y=0, eps.M=0, clever_coef.M=0, rel.tol = 0.01)$value})

  ################## One-step Estimator ##################
  if ('onestep' %in% estimator){

    # E(Dstar) for Y|M,A,X
    or_weight <- (A==alt)*{1/p.alt}*p.M.aX*{1/p.M.altX} # I(A=a')/p(a') * p(M=1|A=a,X)/p(M=1|A=a',X)
    EDstar_or <- mean(or_weight*(Y-or_pred_alt)) # weight*(Y - E(Y|M,a',X))

    # E(Dstar) for M=1|A,X
    m_weight <- {1/p.alt}*ratio.A.X*(or_pred_alt - kappa) # p(a') * p(a'|X)/p(a|X) * { E(Y|M,A=a',X) - int E(Y|M,A=a',X)p(M|A=a,X) dM }

    EDstar_M <- mean((A==a)*m_weight) # I(A=a)/p(a') * p(a'|X)/p(a|X) * { E(Y|M,A=a',X) - int E(Y|M,A=a',X)p(M|A=a,X) dM }

    # at tangent space of A,X
    ax_weight <- (A==alt)*{1/p.alt}


    ######################
    # estimate E[Y(a)]
    ######################

    # estimated psi
    estimated_psi = EDstar_or + EDstar_M + mean(ax_weight*kappa) # estimated psi

    # EIF
    EIF <- or_weight*(Y-or_pred_alt)+ #line 1 of the EIF equation
      (A==a)*m_weight+ #line 2 of the EIF equation
      ax_weight*(kappa - estimated_psi) #line 3 of the EIF equation

    # confidence interval
    lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
    upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))

    onestep.out <- list(estimated_psi=estimated_psi, # estimated parameter
                        lower.ci=lower.ci, # lower bound of 95% CI
                        upper.ci=upper.ci, # upper bound of 95% CI
                        p.M.aX=p.M.aX,  # estimated M=1|A,X
                        p.a1.X=p.a1.X,  # estimated A=1|X
                        or_pred_alt=or_pred_alt, # estimated E(Y|M,A,X)
                        #
                        EIF=EIF, # EIF
                        EDstar=c(EDstar_or,EDstar_M)) # E(Dstar) for Y|M,A,X and M|A,X, and A,X

  }

  if('tmle' %in% estimator){
    ################## TMLE Estimator ##################

    ######################
    # update p(M|A=a,X)
    #####################

    # the gradient of p(M|A=a,X)
    D_M <- {1/p.alt}*ratio.A.X*(or_pred_alt - kappa)

    # loss function for p(M|A=a,X) is -I(A=a)*[log(p(M|A=a,X))+log{1+eps_f*D_M]
    # minimize the loss is the same as minimize -I(A=a)*log{1+eps_f*D_M}
    loss_M <- function (eps.M) {
      return(-1*mean((A==a)*log((1+eps.M*D_M))))
    }

    # range of eps.M that make the density positive
    range_eps.M <- function(eps.M){ range_eps.M <- c(max(max(-1/D_M[D_M>0]),-1e-5),min(min(-1/D_M[D_M<0]),1e-5))}


    # find the optimal eps2 that minimize the loss function
    eps.M <- optimize(loss_M, range_eps.M(D_M), tol=0.0001, maximum=F)$minimum # the domain of eps2 is not R

    # update p(M|A=a,X)
    p.M.aX_updated <- p.M.aX*(1+eps.M*D_M)



    ######################
    # Update E(Y|M,alt,X)
    ######################

    ## clever coefficient for E(Y|M,alt,X)
    clever_coef.Y <- {1/p.alt}*p.M.aX_updated*{1/p.M.altX} # used for fitting the target regression, and updating E(Y|M,a',X)
    or_weight <- (A==alt)*clever_coef.Y  # used if calculating EIF

    # offset term for E(Y|M,alt,X)
    offset.Y <- or_pred_alt

    model.Y <- glm( Y ~ offset(offset.Y)+1, weights=or_weight)

    eps.Y = coef(model.Y)

    # update E(Y|M,alt,X)
    or_pred_alt_updated <- or_pred_alt + eps.Y

    # calculated E(Dstar) for Y|M,alt,X
    EDstar_or <- mean(or_weight*(Y-or_pred_alt_updated)) # weight*(Y - E(Y|M,a',X))

    # estimated psi
    ax_weight <- (A==alt)*{1/p.alt}
    kappa <- sapply(1:n, function(i){integrate(integrand,lower = range(M)[1], upper = range(M)[2], a=a,x=X[i,], eps.Y=eps.Y, eps.M=eps.M, clever_coef.M=D_M[i], rel.tol = 0.01)$value})

    # E(Dstar) for M=1|A,X
    m_weight <- {1/p.alt}*ratio.A.X*(or_pred_alt_updated - kappa) # p(a') * p(a'|X)/p(a|X) * { E(Y|M,A=a',X) - int E(Y|M,A=a',X)p(M|A=a,X) dM }

    EDstar_M <- mean((A==a)*m_weight) # I(A=a)/p(a') * p(a'|X)/p(a|X) * { E(Y|M,A=a',X) - int E(Y|M,A=a',X)p(M|A=a,X) dM }

    estimated_psi = mean(ax_weight*kappa)

    # EIF
    EIF <- or_weight*(Y-or_pred_alt)+ #line 1 of the EIF equation
      (A==a)*{m_weight}+ #line 2 of the EIF equation
      ax_weight*(kappa - estimated_psi) #line 3 of the EIF equation


    # confidence interval
    lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
    upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))

    tmle.out <- list(estimated_psi=estimated_psi, # estimated parameter
                     lower.ci=lower.ci, # lower bound of 95% CI
                     upper.ci=upper.ci, # upper bound of 95% CI
                     p.M.aX=p.M.aX_updated,  # estimated M=1|A,X
                     p.a1.X=p.a1.X,  # estimated A=1|X
                     or_pred_alt=or_pred_alt_updated, # estimated E(Y|M,A,X)
                     #
                     EIF=EIF, # EIF
                     EDstar=c(EDstar_or,EDstar_M) # EIF
    ) # number of iterations for M|A,X and A|X to converge
  }


  if(length(estimator)==2){
    return(list(TMLE=tmle.out,Onestep=onestep.out))
  }else{
    if ('onestep'==estimator){return(list(Onestep=onestep.out))}
    else if ('tmle'==estimator){return(TMLE=tmle.out)}
  }

}
