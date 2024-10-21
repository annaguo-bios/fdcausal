## TMLE np-dnorm ====
#' Function for implementing TMLE with the mediator density \eqn{p(M|A,X)} estimated assuming normal distribution.
#' @param a Treatment level at which the average counterfactual outcome is computed
#' @param data A Dataframe contains treatment, mediators, outcome, and measured confounders
#' @param treatment Variable name for the unvariate binary treatment
#' @param mediator Variable name for the continuous univariate mediator
#' @param outcome Variable name for the continuous univariate outcome
#' @param covariates Variable name for the measured confounders
#' @param onestep A logical indicator determines whether one-step estimation is executed. When 'onestep=T', the one-step estimation result is provided. Conversely, if 'onestep=F', the result is withheld.
#' @param n.iter The maximum number of iterations performed when iteratively updating the mediator density and propensity score.
#' @param eps A logical indicator determines the stopping criteria used when iteratively updating the mediator density and propensity score. The default is 'eps=T'.
#' When 'eps=T', \eqn{\sqrt{\epsilon_2^2+\epsilon_3^2}} is used. When 'eps=F', \eqn{max(|\Phi_M|,|\Phi_A|} is used. In general, adoption of 'eps=T' results in better convergence while it takes longer time.
#' Conversely, adoption of 'eps=F' usually requires less time but can result in algorithm divergence.
#' @param cvg.criteria A numerical value representing the convergence criteria when iteratively updating the mediator density and propensity score. The default value is 0.01, meaning update stops when stopping criteria < 0.01. The stopping criteria is chosen by the eps.
#' @param formulaY Regression formula for the outcome regression of Y on M, A, X. The default is 'Y ~ 1+ M + A + X'.
#' @param formulaA Regression formula for the propensity score regression of A on X. The default is 'A ~ 1 + X'.
#' @param linkA The link function used for the logistic regression of A on X. The default is the 'logit' link.
#' @param truncate_lower A numeric variable, setting lower bound for the truncated propensity score. The default is 0.
#' @param truncate_upper A numeric variable, setting upper bound for the truncated propensity score. The default is 1.
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
#' res <- TMLE.np.dnorm(1,data=continuousY_continuousM,
#' treatment="A", mediator="M", outcome="Y", covariates="X")
#' }
#' @import np
#' @importFrom MASS mvrnorm
#' @importFrom dplyr %>% mutate select
#' @export

TMLE.np.dnorm <- function(a,data,treatment="A", mediator="M", outcome="Y", covariates=c("X1","X2","X3"),
                    onestep=T, n.iter=500, eps=T, cvg.criteria=0.01,
                    formulaY="Y ~ .", formulaA="A ~ .", linkA="logit",
                    truncate_lower=0, truncate_upper=1){

  # attach(data, warn.conflicts=FALSE)

  # Variables
  A <- data[,treatment]
  M <- data[,mediator]
  X <- data[,covariates, drop=F]
  Y <- data[,outcome]
  dat_MAX <- data.frame(M=M,A=A,X)
  dat_MaX <- data.frame(M=M,A=a,X)

  n <- nrow(data)

  if (length(covariates)==1){ # covariaets is of length 1, rename the covariate to "X"
    covariaets = "X"
  }

  ######################
  ## store data for TMLE
  ######################
  # initialize eps2, eps3
  eps2 <- 1
  eps3 <- 1

  # record the values of eps2, eps3, and eps2n3 over iterations
  eps2_vec <- c(0)
  eps3_vec <- vector(mode = "numeric")

  # place holder for clever coefficient 3
  clever_coef3 <- 0

  # cumulative summation of eps*clever-coefficient
  clever_coef3_add <- 0

  # matrix to store p(A=1|X) and theta_x over iterations
  p.a1.X.mat <- as.matrix(rep(0.5,n)) # rep(0.5,n) is a random initial place holder
  theta_x.mat <- as.matrix(rep(0.5,n))

  # record E(Dstar) for M|A,X and propensity score over iterations
  EDstar_M <- 1 # initial value
  EDstar_ps <- 1 # initial value
  EDstar_M.vec <- vector(mode = "numeric")
  EDstar_ps.vec <- vector(mode = "numeric")

  # iterations
  iter <- 0

  ######################
  ## TMLE initialization
  ######################

  # new data sets
  data_Aa = data.frame(M=M, A=A, X)
  data_A1 = data.frame(M=M, A=1, X)
  data_A0 = data.frame(M=M, A=0, X)

  ## outcome regression

  or_fit <- lm(as.formula(formulaY), data=dat_MAX)
  or_pred <- predict(or_fit)
  or_pred_a1 <- predict(or_fit, newdata=data_A1)
  or_pred_a0 <- predict(or_fit, newdata=data_A0)

  # updated function for outcome regression
  f.or <- function(M,A,X,eps1=0){predict(or_fit, newdata = data.frame(M=M,A=A,X,row.names = NULL))+eps1}

  ## M|A,X
  # Methods on density estimation: assuming normal distribution for p(M|A,X)

  m.fit <- lm(M~., data=data.frame(A,X))
  m.coef <- m.fit$coefficients # coefficient of the regression
  l.coef <- length(m.coef) # length of the coefficient
  m.sd <- sd(M-predict(m.fit))

  p.M.AX <- dnorm(M,m.coef[1]+m.coef[2]*A+ as.matrix(X) %*% t(m.coef[3:l.coef]),m.sd)
  p.M.aX <- dnorm(M,m.coef[1]+m.coef[2]*a+ as.matrix(X) %*% t(m.coef[3:l.coef]),m.sd)


  ## propensity score


  if (truncate_lower!=0 | truncate_upper!=1){ # under weak overlapping issue, it's more stable to run A~X via linear regression

    ps_fit <- lm(as.formula(formulaA), data=X)
    p.a1.X <- predict(ps_fit)

  } else { # without weak overlapping issue. Run A~X via logistic regression

    ps_fit <- glm(as.formula(formulaA), data=X,  family = binomial(linkA))
    p.a1.X <- predict(ps_fit, type = "response")  # p(A=1|X)
  }


  p.a1.X[p.a1.X < truncate_lower] <- truncate_lower
  p.a1.X[p.a1.X > truncate_upper] <- truncate_upper

  # updated function for propensity score
  f.ps <- function(x,clever_coef3_add){

    # define new data set for prediction
    newdata <- data.frame(x)
    names(newdata) <- covariates # assign column name

    p.a1.X <- predict(ps_fit, newdata = newdata, type = "response")

    p.a1.X[p.a1.X < truncate_lower] <- truncate_lower
    p.a1.X[p.a1.X > truncate_upper] <- truncate_upper

    return(plogis(qlogis(p.a1.X)+clever_coef3_add))
  }


  # cumulative product of (1+eps2*D(M)). A preparation function for f.m
  cumprod.M <- function(m,x,row.indicator){

    # p(A=1|X)
    p.a1.X <- p.a1.X.mat[row.indicator,]

    # initial p(A=a|X)
    p.a.X <- a*p.a1.X+(1-a)*(1-p.a1.X)

    # xi(M,X)
    newdata_A1 <- data.frame(M=m, A=1, x,row.names = NULL); names(newdata_A1) <- c("M","A",covariates) # new data frame for prediction
    newdata_A0 <- data.frame(M=m, A=0, x,row.names = NULL); names(newdata_A0) <- c("M","A",covariates) # new data frame for prediction

    or_pred_a1 <- predict(or_fit, newdata=newdata_A1)
    or_pred_a0 <- predict(or_fit, newdata=newdata_A0)

    xi <- or_pred_a1%*%t(p.a1.X)+or_pred_a0%*%(1-t(p.a1.X))

    # theta_x
    theta_x <- theta_x.mat[row.indicator,]

    # xi-theta_x
    diff <- sweep(xi,2,theta_x,FUN = "-")

    # cumulative product of (1+eps2*D(M))
    return(
      apply(1+sweep(sweep(diff,2,p.a.X,"/"),2,eps2_vec,FUN = "*"), 1, prod)
    )
  }

  # updated density function for M|A,X
  f.m <- function(M,A,X,cumprod.M,row.indicator){

    pM.AX <- dnorm(M,m.coef[1]+m.coef[2]*A+ as.matrix(X) %*% t(m.coef[3:l.coef]),m.sd)

    return(pM.AX*cumprod.M(M,X,row.indicator))
  }

  # integrand of theta_x: int E(Y|M,A,X)p(A|X)p(M|A=a,X)dA
  integrand <- function(m,a,x,f.or,f.ps,f.m,eps1,cumprod.M,clever_coef3_add,row.indicator){
    or_pred1 <- f.or(m,1,x,eps1) # A=1
    or_pred0 <- f.or(m,0,x,eps1) # A=0

    ps_pred <- f.ps(x,clever_coef3_add) # p(A=1|X)

    xi <- or_pred1*ps_pred+or_pred0*(1-ps_pred) # int E(Y|M,A,X) p(A|X)dA
    M_pred <- f.m(m,a,x,cumprod.M,row.indicator) # p(M|A=a,X)
    integrand <- xi*M_pred # int E(Y|M,A,X) p(A|X) p(M|A=a,X)dA
    return(integrand)
  }

  # eta(1,X)-eta(0,X)
  integrand.eta.diff <- function(m,x,a,cumprod.M,row.indicator){
    or_pred1 <- f.or(m,1,x,eps1=0) # A=1
    or_pred0 <- f.or(m,0,x,eps1=0) # A=0

    M_pred <- f.m(m,a,x,cumprod.M,row.indicator) # p(M|A=a,X)
    integrand <- (or_pred1-or_pred0)*M_pred

    return(integrand)
  }

  ## initial estimates

  # initial p(A=a|X)
  p.a.X <- a*p.a1.X+(1-a)*(1-p.a1.X)

  # int E(Y|M,A,X) p(A|X)dA
  xi <- or_pred_a1*p.a1.X+or_pred_a0*(1-p.a1.X)

  # theta_x
  theta_x <- sapply(1:n, function(i){integrate(integrand,lower = range(M)[1], upper = range(M)[2], a=a,x=X[i,],f.or=f.or, f.ps=f.ps, f.m=f.m, eps1=0, cumprod.M=cumprod.M, clever_coef3_add=rep(0,n)[i],row.indicator=i,rel.tol = 0.01)$value})

  ################## One-step Estimator ##################
  if (onestep==T){

    # eta(1,X)-eta(0,X)
    eta_diff <- sapply(1:n, function(i){integrate(integrand.eta.diff,lower = range(M)[1], upper = range(M)[2], a=a,x=X[i,],cumprod.M=cumprod.M, row.indicator=i)$value})


    ######################
    # E[Dstar] calculations
    ######################

    # E(Dstar) for E(Y|M,A,X)
    or_weight <- p.M.aX/p.M.AX
    EDstar_or <- mean(or_weight*(Y-or_pred))

    #E(Dstar) for M=1|A,X
    m_weight <- (A==a)*{1/p.a.X}
    xi <- or_pred_a1*p.a1.X+or_pred_a0*(1-p.a1.X)

    EDstar_M <- mean(m_weight*(xi-theta_x))

    # E(Dstar) for A=a|X
    EDstar_ps <- mean(eta_diff*(A-p.a1.X))


    ######################
    # estimate E[Y(a)]
    ######################

    # estimated psi
    estimated_psi = mean(theta_x)+EDstar_or+EDstar_M+EDstar_ps

    # EIF
    EIF <- or_weight*(Y-or_pred)+ #line 1 of the EIF equation
      m_weight*(xi-theta_x)+
      eta_diff*(A-p.a1.X)+
      theta_x - mean(theta_x) #line 4 of the EIF equation

    # confidence interval
    lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
    upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))


    onestep.out <- list(estimated_psi=estimated_psi, # estimated parameter
                        theta_x=theta_x, # theta(x)
                        lower.ci=lower.ci, # lower bound of 95% CI
                        upper.ci=upper.ci, # upper bound of 95% CI
                        EIF=EIF, # EIF
                        EDstar=c(EDstar_or,EDstar_M, EDstar_ps)) # E(Dstar) for Y|M,A,X and M|A,X, and A|X

  }

  ################## TMLE Estimator ##################
  while ((eps && sqrt(eps2^2+eps3^2)> cvg.criteria) || (!eps && max(abs(EDstar_M),abs(EDstar_ps)) > cvg.criteria)  && iter<n.iter) {

    ######################
    # update p(M|A=a,X)
    #####################

    # the gradient of p(M|A=a,X)
    D_M <- (xi-theta_x)/p.a.X

    # loss function for p(M|A=a,X) is -I(A=a)*[log(p(M|A=a,X))+log{1+eps_f/g_a(X)*(xi(M,X)-theta(X))}]
    # minimize the loss is the same as minimize -I(A=a)*log{1+eps_f/g_a(X)*(xi(M,X)-theta(X))}
    loss_M <- function (eps2) {
      return(-1*mean((A==a)*log((1+eps2*D_M))))
    }

    range_eps2 <- function(eps2){
      range_eps2 <- c(max(max(-1/D_M[D_M>0]),-1e-5),min(min(-1/D_M[D_M<0]),1e-5))
    }


    # find the optimal eps2 that minimize the loss function
    eps2 <- optimize(loss_M, range_eps2(D_M), tol=0.0001, maximum=F)$minimum # the domain of eps2 is not R

    # f.m is updated via updating eps2_vec: f.m=f.m^0 * \prod [1+eps*(xi-theta_x)]
    eps2_vec <- c(eps2_vec,eps2)

    # update matrix that store p(A=1|X) and theta_x over iterations
    p.a1.X.mat <- cbind(p.a1.X.mat,p.a1.X)
    theta_x.mat <- cbind(theta_x.mat,theta_x)

    ######################
    # update p(A=1|X)
    ######################

    # clever coefficient for propensity score: eta(1,X)-eta(0,X)
    clever_coef3=sapply(1:n, function(i){integrate(integrand.eta.diff,lower = range(M)[1], upper = range(M)[2], a=a,x=X[i,],cumprod.M=cumprod.M,row.indicator=i,rel.tol = 0.01)$value})

    # offset term for the propensity score
    offset_ps <- qlogis(p.a1.X)

    # derive eps3
    ps_model <- glm(
      A ~ offset(offset_ps)+clever_coef3-1, family=binomial(), start=0
    )

    eps3 <- coef(ps_model)
    eps3_vec <- c(eps3_vec,eps3)

    # update cumulative summation of eps3*clever coefficient
    clever_coef3_add = clever_coef3_add + eps3*(clever_coef3)

    # updated propensity score
    p.a1.X <- plogis(qlogis(p.a1.X)+eps3*(clever_coef3))

    # p(A=a|X)
    p.a.X <- a*p.a1.X+(1-a)*(1-p.a1.X)


    # update theta_x
    theta_x <- sapply(1:n, function(i){integrate(integrand,lower = range(M)[1], upper = range(M)[2], a=a,x=X[i,],f.or=f.or, f.ps=f.ps, f.m=f.m, eps1=0, cumprod.M=cumprod.M, clever_coef3_add[i],row.indicator=i,rel.tol = 0.01)$value})

    # update xi
    xi <- or_pred_a1*p.a1.X+or_pred_a0*(1-p.a1.X)

    # E(Dstar) for M|A,X: I(A=a)/p(A=a|X)*[xi-(theta_x+eps2*xi)]
    EDstar_M <- mean((A==a)/p.a.X*(xi-theta_x))
    EDstar_M.vec <- c(EDstar_M.vec,EDstar_M)


    # E(Dstar)-iteration1
    EDstar_ps <- mean(clever_coef3*(A-p.a1.X))
    EDstar_ps.vec <- c(EDstar_ps.vec,EDstar_ps)

    iter <- iter+1
    print(paste0("Iteration: ",iter))
  }

  ######################
  # update E[Y|A,M,X]
  ######################
  # or_weight <- p.M.aX/p.M.AX

  p.M.aX_updated <- sapply(1:n, function(i){f.m(M=M[i],A=a,X[i,],cumprod.M,row.indicator=i)})

  or_weight <- p.M.aX_updated/((A==a)*p.M.aX_updated+(1-(A==a))*p.M.AX)

  # one iteration
  or_model <- glm(
    Y ~ offset(or_pred)+1, weights = or_weight
  )

  eps1 = coef(or_model)

  # updated outcome regression
  or_pred = or_pred + eps1

  ######################
  # E[Dstar] calculations
  ######################

  # E(Dstar) for E(Y|M,A,X)
  EDstar_or <- mean(or_weight*(Y-or_pred))

  #E(Dstar) for M=1|A,X
  EDstar_M <- mean((A==a)/p.a.X*(xi-theta_x))

  # E(Dstar) for A=a|X
  EDstar_ps <- mean(clever_coef3*(A-p.a1.X))


  ######################
  # estimate E[Y(a)]
  ######################

  # update theta_x upon updating E(Y|M,A,X): theta(X) = int E(Y|M,A,X)p(M|A=a,X)p(A|X)dM dA
  theta_x <- theta_x+eps1

  # update xi upon updating E(Y|M,A,X): xi(M,X) = int E(Y|M,A,X)p(A|X) dA
  xi <- xi+eps1

  # estimated psi
  estimated_psi = mean(theta_x)

  # EIF
  EIF <- or_weight*(Y-or_pred)+ #line 1 of the EIF equation
    {I(A==a)/p.a.X}*(xi-theta_x)+ #line 2 of the EIF equation
    clever_coef3*(A-p.a1.X)+ #line 3 of the EIF equation
    theta_x - estimated_psi #line 4 of the EIF equation

  # confidence interval
  lower.ci <- estimated_psi-1.96*sqrt(mean(EIF^2)/nrow(data))
  upper.ci <- estimated_psi+1.96*sqrt(mean(EIF^2)/nrow(data))

  tmle.out <- list(estimated_psi=estimated_psi, # estimated parameter
                   lower.ci=lower.ci, # lower bound of 95% CI
                   upper.ci=upper.ci, # upper bound of 95% CI
                   p.M.aX=p.M.aX,  # estimated M|A=a,X
                   p.M.AX=p.M.AX, # estimated M|A,X
                   p.a1.X_updated=p.a1.X,  # estimated A=1|X
                   or_pred_updated=or_pred, # estimated E(Y|M,A,X)
                   EIF=EIF, # EIF
                   #
                   EDstar=c(EDstar_or,EDstar_M, EDstar_ps), # E(Dstar) for Y|M,A,X and M|A,X, and A|X
                   EDstar_M.vec=EDstar_M.vec, # E(Dstar) for M|A,X over iterations
                   EDstar_ps.vec=EDstar_ps.vec, # E(Dstar) for A|X over iterations
                   #
                   eps2_vec=eps2_vec, # vector of eps2 over iterations
                   eps3_vec=eps3_vec, # vector of eps3 over iterations
                   iter=iter) # number of iterations for M|A,X and A|X to converge

  if (onestep==T){return(list(TMLE=tmle.out,Onestep=onestep.out))}else{return(TMLE=tmle.out)}

}
