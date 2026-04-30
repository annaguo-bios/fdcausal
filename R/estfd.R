## usethis namespace: start
#' @importFrom dplyr %>%
#' @importFrom SuperLearner SuperLearner
#' @importFrom np npcdensbw npcdens
#' @importFrom densratio densratio
#' @importFrom stats coef cov lm glm as.formula binomial gaussian integrate optimize plogis predict qlogis
#' @importFrom mvtnorm dmvnorm
## usethis namespace: end
NULL

#' Estimate ATE and ATT via TMLE and/or one-step estimation.
#'
#' Estimates the ATE \eqn{E(Y^1) - E(Y^0)}, ATT \eqn{E(Y^1) - E(Y^0|A=1)},
#' or marginal counterfactual means \eqn{E(Y^a)}, \eqn{E(Y^a|A=1-a)} under the front-door model.
#'
#' @param a Treatment level(s). Use \code{c(1,0)} for ATE/ATT; use \code{1} or
#'   \code{0} for a single counterfactual mean.
#' @param data Data frame containing treatment, mediator(s), outcome, and covariates.
#' @param treatment Name of the binary treatment variable.
#' @param mediators Name of the mediator variable. Can be a single variable or a vector of variable names for multiple mediators.
#' @param outcome Name of the univariate outcome variable.
#' @param covariates Names of the measured confounders.
#' @param estimator Which estimator to use: \code{"onestep"} or \code{"tmle"}.
#' @param ATT Logical. If \code{TRUE}, estimate the ATT; otherwise estimate the ATE. Default \code{FALSE}.
#'
#' @param mediator.method Method for mediator density estimation (ignored when M is binary,
#'   in which case regression is used automatically):
#'   \describe{
#'     \item{\code{"bayes"}}{Estimates the density ratio via Bayes' rule:
#'       \eqn{p(M|A,X)/p(M|a,X) = [p(a|M,X)/p(A|M,X)] / [p(a|X)/p(A|X)]}.
#'       Avoids directly updating \eqn{p(M|A,X)} in constructing TMLE.}
#'     \item{\code{"densratio"}}{Estimates the density ratio using \code{\link[densratio]{densratio}}.
#'       Requires all mediators to be numeric. Avoids directly updating \eqn{p(M|A,X)} in constructing TMLE.}
#'     \item{\code{"dnorm"}}{Assumes \eqn{p(M|A,X)} is normal to estimate the density ratio.
#'       Avoids directly updating \eqn{p(M|A,X)} in constructing TMLE.}
#'     \item{\code{"np"}}{Directly updates \eqn{p(M|A,X)}. Uses \code{\link[np]{npcdens}}
#'       when \code{np.dnorm=FALSE}, or a normal approximation when \code{np.dnorm=TRUE}.
#'       Only available for univariate continuous M due to computational cost. Requires updating \eqn{p(M|A,X)}
#'       at each iteration of TMLE, which can be computationally intensive.}
#'   }
#' @param np.dnorm Logical. Only relevant when \code{mediator.method = "np"}.
#'   If \code{TRUE}, estimate \eqn{p(M|A,X)} assuming normality; if \code{FALSE},
#'   use \code{\link[np]{npcdens}}. Default \code{FALSE}.
#'
#' @param superlearner Logical. If \code{TRUE}, use SuperLearner
#'   (\code{\link[SuperLearner]{SuperLearner}}) for all nuisance estimates
#'   (outcome regression, mediator density, propensity score). Recommended when
#'   relationships among variables are complex.
#' @param crossfit Logical. If \code{TRUE}, combine SuperLearner with cross-fitting.
#' @param K Number of folds for cross-fitting. Default \code{5}.
#' @param lib Learner library passed to SuperLearner.
#'
#' @param n.iter Maximum number of iterations for the iterative update of the
#'   mediator density and propensity score.
#' @param eps Logical. Stopping criterion for TMLE iterative updates. If \code{TRUE}
#'   (default), stop when \eqn{\sqrt{\epsilon_2^2 + \epsilon_3^2}} is small;
#'   if \code{FALSE}, stop when \eqn{\max(|\Phi_M|, |\Phi_A|)} is small.
#'   \code{FALSE} gives better convergence but is slower.
#' @param cvg.criteria Convergence threshold for iterative updates. Default \code{0.01}.
#'
#' @param formulaY Outcome regression formula (\eqn{Y \sim M + A + X}). Default \code{Y ~ 1 + M + A + X}.
#' @param linkY_binary Link function for binary Y. Default \code{"logit"}.
#' @param formulaA Propensity score formula (\eqn{A \sim X}). Default \code{A ~ 1 + X}.
#' @param linkA Link function for the propensity score. Default \code{"logit"}.
#' @param formulaM Mediator regression formula; only used when M is binary. Default \code{M ~ 1 + A + X}.
#' @param linkM_binary Link function for binary M. Default \code{"logit"}.
#' @param formula_bayes Regression formula for \eqn{A \sim M + X}; only used when
#'   \code{mediator.method = "bayes"}. Default \code{A ~ 1 + M + X}.
#' @param link_bayes Link function for the Bayes regression of A on M, X. Default \code{"logit"}.
#' @param truncate_lower Lower truncation bound for the propensity score. Default \code{0.01}.
#' @param truncate_upper Upper truncation bound for the propensity score. Default \code{0.99}.
#' @param boundedsubmodelY Logical. Whether to use a bounded submodel when targeting
#'   the outcome regression for discrete Z. Default \code{FALSE}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{ATE}}{Estimated ATE: \eqn{E(Y^1) - E(Y^0)}.}
#'     \item{\code{estimated_psi}}{Estimated \eqn{E(Y^a)}.}
#'     \item{\code{lower.ci}, \code{upper.ci}}{95\% confidence interval bounds.}
#'     \item{\code{theta_x}}{\eqn{\int E(Y|M,A,X)\,p(M|A=a,X)\,p(A|X)\,dM\,dA}.}
#'     \item{\code{p.m1.aX}}{\eqn{p(M=1|A=a,X)} (binary M only).}
#'     \item{\code{p.a1.X}}{Estimated propensity score \eqn{p(A=1|X)}.}
#'     \item{\code{or_pred}}{Fitted outcome regression \eqn{E(Y|M,A,X)}.}
#'     \item{\code{EIF}}{Efficient influence function evaluated at observed data.}
#'     \item{\code{EDstar}}{EIF projections onto the tangent spaces of \eqn{Y|M,A,X},
#'       \eqn{M|A,X}, and \eqn{A|X}.}
#'     \item{\code{EDstar_M.vec}, \code{EDstar_ps.vec}}{Per-iteration mean EIF projections
#'       for \eqn{M|A,X} and \eqn{A|X}; should approach 0 at convergence.}
#'     \item{\code{eps2_vec}, \code{eps3_vec}}{Per-iteration submodel indices for
#'       \eqn{M|A,X} and \eqn{A|X}; should approach 0 at convergence.}
#'     \item{\code{iter}}{Number of iterations until convergence.}
#'   }
#' @examples
#' \donttest{
#' # ATT estimation. For binary outcome Y and binary mediator M.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <-estfd(a=c(1,0),data=binaryY_binaryM,
#' treatment="A", mediators="M", outcome="Y", covariates="X",
#'   linkA="identity",ATT=TRUE)
#'
#' # ATE estimation. For binary outcome Y and binary mediator M.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <-estfd(a=c(1,0),data=binaryY_binaryM,
#' treatment="A", mediators="M", outcome="Y", covariates="X",
#'   linkA="identity")
#'
#' # ATE estimation. For continuous outcome Y and binary mediator M
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <-estfd(a=c(1,0),data=continuousY_binaryM,
#' treatment="A", mediators="M", outcome="Y",
#' covariates="X",  linkA="identity")
#'
#' # ATE estimation. For continuous outcome Y and binary mediator M.
#' # Data is generated under p(A=1|X) = 0.001 + 0.998X. And X~Uniform(0,1).
#' # Therefore, this dataset suffers from weak overlapping.
#' # Below we apply truncation to the propensity score to truncate it between (0.001, 0.999).
#' res <- estfd(a=c(1,0),data=continuousY_binaryM_weakoverlap,
#' treatment="A", mediators="M", outcome="Y", covariates="X",
#'  linkA="identity", truncate_lower=0.001, truncate_upper=0.999)
#'
#' # ATE estimation. For continuous outcome Y and univariate continuous mediator M.
#' # Using 'np' method for mediator density estimation.
#' # Setting np.dnorm=F, so that mediator density is estimated via the np function.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <-estfd(a=c(1,0),data=continuousY_continuousM,
#' treatment="A", mediators="M", outcome="Y", covariates="X",
#'  linkA="identity", mediator.method="np", np.dnorm=FALSE)
#'
#' # ATE estimation. For continuous outcome Y and univariate continuous mediator M.
#' # Using 'np' method for mediator density estimation.
#' # Setting np.dnorm=T, so that mediator density is estimated assuming normal distribution.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <-estfd(a=c(1,0),data=continuousY_continuousM,
#' treatment="A", mediators="M", outcome="Y", covariates="X",
#'  linkA="identity", mediator.method="np", np.dnorm=TRUE)
#'
#' # ATE estimation. For continuous outcome Y and univariate continuous mediator M.
#' # Using 'densratio' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.001 + 0.998X. And X~Uniform(0,1).
#' # Therefore, this dataset suffers from weak overlapping.
#' # Below we apply truncation to the propensity score to truncate it between (0.001, 0.999).
#' res <-estfd(a=c(1,0),data=continuousY_continuousM_weakoverlap,
#' treatment="A", mediators="M", outcome="Y", covariates="X",
#'  linkA="identity", mediator.method="densratio",
#' truncate_lower=0.001, truncate_upper=0.999)
#'
#' # ATE estimation. For continuous outcome Y and bivariate mediator M.
#' # Using 'densratio' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <-estfd(a=c(1,0),data=continuousY_bivariateM,
#' treatment="A", mediators=c("M.1","M.2"), outcome="Y", covariates="X",
#'   linkA="identity", mediator.method="densratio")
#'
#' # ATE estimation. For continuous outcome Y and bivariate mediator M.
#' # Using 'bayes' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <-estfd(a=c(1,0),data=continuousY_bivariateM,
#' treatment="A", mediators=c("M.1","M.2"), outcome="Y", covariates="X",
#'   linkA="identity", mediator.method="bayes")
#'
#' # ATE estimation. For continuous outcome Y and bivariate mediator M.
#' # Using 'dnorm' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <-estfd(a=c(1,0),data=continuousY_bivariateM,treatment="A",
#' mediators=c("M.1","M.2"), outcome="Y", covariates="X",
#'  linkA="identity", mediator.method="dnorm")
#'
#' # ATE estimation. For continuous outcome Y and quadrivariate mediator M.
#' # Using 'densratio' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <-estfd(a=c(1,0),data=continuousY_quadrivariateM,treatment="A",
#' mediators=c("M.1","M.2","M.3","M.4"), outcome="Y", covariates="X",
#'  linkA="identity", mediator.method="densratio")
#'
#' # ATE estimation. For continuous outcome Y and quadrivariate mediator M.
#' # Using 'bayes' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <-estfd(a=c(1,0),data=continuousY_quadrivariateM,
#' treatment="A", mediators=c("M.1","M.2","M.3","M.4"), outcome="Y", covariates="X",
#'   linkA="identity", mediator.method="bayes")
#'
#' # ATE estimation. For continuous outcome Y and quadrivariate mediator M.
#' # Using 'dnorm' method for mediator density ratio estimation.
#' # Data is generated under p(A=1|X) = 0.3 + 0.2X.
#' # Therefore, setting link for propensity score to be "identity".
#' res <-estfd(a=c(1,0),data=continuousY_quadrivariateM,
#' treatment="A", mediators=c("M.1","M.2","M.3","M.4"), outcome="Y", covariates="X",
#'   linkA="identity", mediator.method="dnorm")
#' }
#'
#' @import np densratio SuperLearner mvtnorm stats
#' @importFrom dplyr %>% mutate select
#' @importFrom MASS mvrnorm
#' @export
#'
#'
estfd <- function(a,data,treatment, mediators, outcome, covariates,
                 mediator.method="bayes", np.dnorm=T, superlearner=F,crossfit=F,K=5,
                 lib = c("SL.glm","SL.earth","SL.ranger","SL.mean"), n.iter=500, eps=T, cvg.criteria=0.01,
                 formulaY="Y ~ .", formulaA="A ~ .", formulaM="M~.", linkY_binary="logit", linkA="logit", linkM_binary="logit",
                 formula_bayes="A ~ .",link_bayes="logit",
                 truncate_lower=0.01, truncate_upper=0.99, ATT = F, estimator='onestep',boundedsubmodelY=F){

  # sample size

  n <- nrow(data)

  if (is.vector(a) & length(a)>2){ ## Invalid input ==

    stop("Invalid input. Enter a=c(1,0) for Average Causal Effect estimation. Enter a=1 or a=0 for average counterfactual outcome estimation at the specified treatment level.")

  }else if (is.vector(a) & length(a)==2){ ## ATE or ATT estimate ==

    if(ATT==T){ ## ATT estimation --

      ## E[Y(a2)| a1]
      out.a0 <- ATT.TMLE.all(a=a[2],data=data,treatment=treatment, mediators=mediators, outcome=outcome, covariates=covariates,
                         mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
                         lib = lib, n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
                         formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
                         formula_bayes=formula_bayes,link_bayes=link_bayes,
                         truncate_lower=truncate_lower, truncate_upper=truncate_upper,estimator=estimator,boundedsubmodelY=boundedsubmodelY)

      ## ATT = E[Y(a1)| a1] - E[Y(a2)| a1]
      Y <- data[,outcome]
      A <- data[,treatment]

      TMLE <- list(estimated_psi = mean(Y[A==a[1]]), EIF = (A==a[1])/mean(A==a[1])*(Y - mean(Y[A==a[1]])))
      Onestep <- list(estimated_psi = mean(Y[A==a[1]]), EIF = (A==a[1])/mean(A==a[1])*(Y - mean(Y[A==a[1]])))

      out.a1 <- list(TMLE=TMLE, Onestep=Onestep)

    }else{ ## ATE estimation --


      ## TMLE estimator

      out.a1 <- TMLE.all(a=a[1],data=data,treatment=treatment, mediators=mediators, outcome=outcome, covariates=covariates,
                         mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
                         lib = lib, n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
                         formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
                         formula_bayes=formula_bayes,link_bayes=link_bayes,
                         truncate_lower=truncate_lower, truncate_upper=truncate_upper,estimator=estimator,boundedsubmodelY=boundedsubmodelY)

      out.a0 <- TMLE.all(a=a[2],data=data,treatment=treatment, mediators=mediators, outcome=outcome, covariates=covariates,
                         mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
                         lib = lib, n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
                         formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
                         formula_bayes=formula_bayes,link_bayes=link_bayes,
                         truncate_lower=truncate_lower, truncate_upper=truncate_upper,estimator=estimator,boundedsubmodelY=boundedsubmodelY)

    } ## end of if-else for ATE/ATT

    estimand <- ifelse(ATT, 'ATT','ATE')

    if('tmle' %in% estimator){

      # run TMLE
      tmle_output_Y1 <- out.a1$TMLE
      tmle_output_Y0 <- out.a0$TMLE

      # estimate E[Y(1)], E[Y(0)], and ATE
      hat_E.Y1 = tmle_output_Y1$estimated_psi
      hat_E.Y0 = tmle_output_Y0$estimated_psi
      hat_ATE = hat_E.Y1 - hat_E.Y0

      # lower CI
      lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((tmle_output_Y1$EIF-tmle_output_Y0$EIF)^2)/n)

      # upper CI
      upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((tmle_output_Y1$EIF-tmle_output_Y0$EIF)^2)/n)

      tmle.out <- list(
                       lower.ci=lower.ci_ATE, # lower bound of 95% CI
                       upper.ci=upper.ci_ATE, # upper bound of 95% CI
                       EIF=tmle_output_Y1$EIF-tmle_output_Y0$EIF # EIF
      )

      tmle.out[[estimand]] <- hat_ATE # ATE estimate

      cat(paste0("TMLE estimated ", estimand ," : ",round(tmle.out[[estimand]],2),"; 95% CI: (",round(tmle.out$lower.ci,2),", ",round(tmle.out$upper.ci,2),")"))

      if(length(estimator)==1){return(list(TMLE=tmle.out,TMLE.Y1=tmle_output_Y1, TMLE.Y0 = tmle_output_Y0))}

    }


    if ('onestep' %in% estimator){

      # run TMLE
      onestep_output_Y1 <- out.a1$Onestep
      onestep_output_Y0 <- out.a0$Onestep

      # estimate E[Y(1)], E[Y(0)], and ATE
      hat_E.Y1 = onestep_output_Y1$estimated_psi
      hat_E.Y0 = onestep_output_Y0$estimated_psi
      hat_ATE = hat_E.Y1 - hat_E.Y0

      # lower CI
      lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((onestep_output_Y1$EIF-onestep_output_Y0$EIF)^2)/n)

      # upper CI
      upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((onestep_output_Y1$EIF-onestep_output_Y0$EIF)^2)/n)

      onestep.out <- list(
                          lower.ci=lower.ci_ATE, # lower bound of 95% CI
                          upper.ci=upper.ci_ATE, # upper bound of 95% CI
                          EIF=onestep_output_Y1$EIF-onestep_output_Y0$EIF # EIF
      )

      onestep.out[[estimand]] <- hat_ATE # ATE estimate

      cat(paste0("Onestep estimated ", estimand ," : ",round(onestep.out[[estimand]],2),"; 95% CI: (",round(onestep.out$lower.ci,2),", ",round(onestep.out$upper.ci,2),")"))

      if(length(estimator)==1){return(list(Onestep=onestep.out,Onestep.Y1=onestep_output_Y1, Onestep.Y0 = onestep_output_Y0))}
    }

  if(length(estimator)==2){return(list(TMLE=tmle.out,Onestep=onestep.out, TMLE.Y1=tmle_output_Y1, TMLE.Y0 = tmle_output_Y0, Onestep.Y1=onestep_output_Y1, Onestep.Y0=onestep_output_Y0))}


  }else if (length(a)==1) { ## E(Y^1) estimate ==

    if(ATT){ # ATT estimation --

      out.a <- ATT.TMLE.all(a=a,data=data,treatment=treatment, mediators=mediators, outcome=outcome, covariates=covariates,
                        mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
                        lib = lib, n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
                        formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
                        formula_bayes=formula_bayes,link_bayes=link_bayes,
                        truncate_lower=truncate_lower, truncate_upper=truncate_upper,estimator=estimator,boundedsubmodelY=boundedsubmodelY)

    }else{ # ATE estimation --


      out.a <- TMLE.all(a=a,data=data,treatment=treatment, mediators=mediators, outcome=outcome, covariates=covariates,
                        mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
                        lib = lib, n.iter=n.iter, eps=eps, cvg.criteria=cvg.criteria,
                        formulaY=formulaY, formulaA=formulaA, formulaM=formulaM, linkY_binary=linkY_binary, linkA=linkA, linkM_binary=linkM_binary,
                        formula_bayes=formula_bayes,link_bayes=link_bayes,
                        truncate_lower=truncate_lower, truncate_upper=truncate_upper,estimator=estimator,boundedsubmodelY=boundedsubmodelY)


    }

    return(out.a)

  } ## end of if-else for input a length

}



