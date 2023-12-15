#' @title A Sample Data for Analysis. Data generated with the following DGP.
#'
#' @format ## `binaryY_binaryM`
#' A data frame with n=500 rows and 5 columns:
#' \describe{
#'   \item{A}{Binary treatment variable.}
#'   \item{M}{Binary mediator variable.}
#'   \item{Y}{Binary outcome variable.}
#'   \item{X}{Uniformly distributed measure confounder within the range of (0,1).}
#'   \item{U}{Normally distributed unmeasured confounder.}
#' }
#'
#' @examples
#' \donttest{
#' # data generated with the following Data Generating Process (DGP)
#' generate_data <- function(n,parA = c(0.3,0.2), parU=c(1,1,1,0), parM = c(-1,1,1,0), parY = c(1, 1, 1, 0), sd.U=1){ # change the parM to c(-1,1,1,0)
#'
#'  X <- runif(n, 0, 1) # p(X)
#'
#'  A <- rbinom(n, 1, (parA[1] + parA[2]*X)) # p(A|X)
#'
#'  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)
#'
#'  M <- rbinom(n,1,plogis(parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X)) # p(M|A,X)
#'
#'  Y <- rbinom(n, 1, plogis(parY[1]*U + parY[2]*M + parY[3]*X + parY[4]*M*X)) # p(Y|U,M,X)
#'
#'  data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)
#'
#'  # propensity score
#'  ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))
#'
#'  # mediator density ratio: p(M|a,X)/p(M|A,X)
#'  m.ratio.a1 <- dbinom(M,1,plogis(parM[1] + parM[2]*1 + parM[3]*X + parM[4]*1*X))/dbinom(M,1,plogis(parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X))
#'  m.ratio.a0 <- dbinom(M,1,plogis(parM[1] + parM[2]*0 + parM[3]*X + parM[4]*0*X))/dbinom(M,1,plogis(parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X))
#'
#'
#'  return(list(data = data,
#'              parA=parA,
#'              parU=parU,
#'              parM=parM,
#'              parY=parY,
#'              sd.U=sd.U,
#'              ps=ps,
#'              m.ratio.a1=m.ratio.a1,
#'              m.ratio.a0=m.ratio.a0))
#'  }
#'
#' }
#'
#'

#'
"binaryY_binaryM"


#' @title A Sample Data for Analysis
#' @format ## `binaryY_continuousM`
#' A data frame with 500 rows and 5 columns:
#' \describe{
#'   \item{A}{Binary treatment variable}
#'   \item{M}{Normally distributed mediator variable}
#'   \item{Y}{Binary outcome variable}
#'   \item{X}{Uniformly distributed measure confounder within the range of (0,1)}
#'   \item{U}{Normally distributed unmeasured confounder}
#' }
#'
#' @examples
#' \donttest{
#' # data generated with the following Data Generating Process (DGP)
#' generate_data <- function(n,parA = c(0.3,0.2), parU=c(1,1,1,0), parM = c(1,1,1,0), parY = c(1, 1, 1, 0), sd.M=1, sd.U=1){
#'
#'  X <- runif(n, 0, 1) # p(X)
#'
#'  A <- rbinom(n, 1, (parA[1] + parA[2]*X)) # p(A|X)
#'
#'  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)
#'
#'  M <- parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X + rnorm(n,0,sd.M) # p(M|A,X)
#'
#'  Y <- rbinom(n,1,plogis(parY[1]*U + parY[2]*M + parY[3]*X + parY[4]*M*X)) # p(Y|U,M,X)
#'
#'  data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)
#'
#'  # propensity score
#'  ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))
#'
#'  # mediator density ratio: p(M|a,X)/p(M|A,X)
#'  m.ratio.a1 <- dnorm(M,parM[1] + parM[2]*1 + parM[3]*X + parM[4]*1*X,sd.M)/dnorm(M,parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X,sd.M)
#'  m.ratio.a0 <- dnorm(M,parM[1] + parM[2]*0 + parM[3]*X + parM[4]*0*X,sd.M)/dnorm(M,parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X,sd.M)
#'
#'  return(list(data = data,
#'              parA=parA,
#'              parU=parU,
#'              parM=parM,
#'              parY=parY,
#'              sd.U=sd.U,
#'              sd.M=sd.M,
#'              ps=ps,
#'              m.ratio.a1=m.ratio.a1,
#'              m.ratio.a0=m.ratio.a0))
#'  }
#' }
#'
"binaryY_continuousM"


#' @title A Sample Data for Analysis
#' @format ## `binaryY_bivariateM`
#' A data frame with 500 rows and 5 columns:
#' \describe{
#'   \item{A}{Binary treatment variable}
#'   \item{M}{Bivariate Normally distributed mediator variable}
#'   \item{Y}{Binary outcome variable}
#'   \item{X}{Uniformly distributed measure confounder within the range of (0,1)}
#'   \item{U}{Normally distributed unmeasured confounder}
#' }
#'
#' @examples
#' \donttest{
#' # data generated with the following Data Generating Process (DGP)
#' generate_data <- function(n,parA = c(0.3,0.2), parU=c(1,1,1,0), parM = matrix(c(1, 1, 1, 0,-1,-0.5,2,0), nrow = 2,byrow = T), parY = c(1, 1, -0.5,1, 0), sd.U=1){
#'
#'  ########################################################
#'  # M is bivariate normal with mean parameter be
#'  #  1  1.0    1    0
#'  # -1 -0.5    2    0
#'
#'  # and the variance covariance matrix be
#'  # 2 1
#'  # 1 3
#'  ########################################################
#'
#'  X <- runif(n, 0, 1) # p(X)
#'
#'  A <- rbinom(n, 1, (parA[1] + parA[2]*X)) # p(A|X)
#'
#'  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)
#'
#'  M <- cbind(parM[1,1] + parM[1,2]*A + parM[1,3]*X + parM[1,4]*A*X,
#'             parM[2,1] + parM[2,2]*A + parM[2,3]*X + parM[2,4]*A*X)+ mvrnorm(n , mu =c(0,0) , Sigma = matrix(c(2, 1, 1, 3), nrow = 2))
#'
#'  Y <- rbinom(n,1,plogis(parY[1]*U + parY[2]*M[,1]+ parY[3]*M[,2] + parY[4]*X + parY[5]*M[,1]*X)) # p(Y|U,M,X)
#'
#'  data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)
#'
#'  # propensity score
#'  ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))
#'
#'  # mediator density ratio: p(M|a,X)/p(M|A,X)
#'  f.m.ratio.a1 <- function(i){
#'    dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*1 + parM[1,3]*X[i] + parM[1,4]*1*X[i],
#'                                parM[2,1] + parM[2,2]*1 + parM[2,3]*X[i] + parM[2,4]*1*X[i]),sigma=matrix(c(2, 1, 1, 3), nrow = 2))/
#'      dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*A[i] + parM[1,3]*X[i] + parM[1,4]*A[i]*X[i], parM[2,1] + parM[2,2]*A[i] + parM[2,3]*X[i] + parM[2,4]*A[i]*X[i]),sigma=matrix(c(2, 1, 1, 3), nrow = 2))
#'  }
#'  f.m.ratio.a0 <- function(i){
#'    dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*0 + parM[1,3]*X[i] + parM[1,4]*0*X[i],
#'                                parM[2,1] + parM[2,2]*0 + parM[2,3]*X[i] + parM[2,4]*0*X[i]),sigma=matrix(c(2, 1, 1, 3), nrow = 2))/
#'      dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*A[i] + parM[1,3]*X[i] + parM[1,4]*A[i]*X[i], parM[2,1] + parM[2,2]*A[i] + parM[2,3]*X[i] + parM[2,4]*A[i]*X[i]),sigma=matrix(c(2, 1, 1, 3), nrow = 2))
#'  }
#'
#'  m.ratio.a1 <- sapply(1:n, f.m.ratio.a1)
#'  m.ratio.a0 <- sapply(1:n, f.m.ratio.a0)
#'
#'  return(list(data = data,
#'              parA=parA,
#'              parU=parU,
#'              parM=parM,
#'              parY=parY,
#'              sd.U=sd.U,
#'              sigma.M=matrix(c(2, 1, 1, 3), nrow = 2),
#'              ps=ps,
#'              m.ratio.a1=m.ratio.a1,
#'              m.ratio.a0=m.ratio.a0))
#'  }
#'
#' }
#'
"binaryY_bivariateM"


#' @title A Sample Data for Analysis
#' @format ## `binaryY_quadrivariateM`
#' A data frame with 500 rows and 5 columns:
#' \describe{
#'   \item{A}{Binary treatment variable}
#'   \item{M}{Quadrivariate Normally distributed mediator variable}
#'   \item{Y}{Binary outcome variable}
#'   \item{X}{Uniformly distributed measure confounder within the range of (0,1)}
#'   \item{U}{Normally distributed unmeasured confounder}
#' }
#'
#' @examples
#' \donttest{
#' # data generated with the following Data Generating Process (DGP)
#' ## generate binary outcome Y, quadrivariate continuous mediator M, single measured covariate X====
#'  generate_data <- function(n,parA = c(0.3,0.2), parU=c(1,1,1,0), parM = matrix(c(1, 1, 1, 0,-1,-0.5,2,0,-1,2,1,0,1,0.5,-1,0), nrow = 4,byrow = T), parY = c(1, 1, -0.5,1,-1,1, 0), sd.U=1, sd.Y=1){
#'    ########################################################
#'    # M is now quadrivariate with parM be
#'    #  1  1.0    1    0
#'    # -1 -0.5    2    0
#'    # -1  2.     1.   0
#'    #  1 0.5  -1.  0
#'
#'    # and the variance covariance matrix be
#'    # 5 -1 0 2
#'    # -1  6 1 0
#'    # 0 1 4 3
#'    # 2 0 3 7
#'    ########################################################
#'
#'    X <- runif(n, 0, 1) # p(X)
#'
#'    A <- rbinom(n, 1, (parA[1] + parA[2]*X)) # p(A|X)
#'
#'    U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)
#'
#'    M <- cbind(parM[1,1] + parM[1,2]*A + parM[1,3]*X + parM[1,4]*A*X,
#'               parM[2,1] + parM[2,2]*A + parM[2,3]*X + parM[2,4]*A*X,
#'               parM[3,1] + parM[3,2]*A + parM[3,3]*X + parM[3,4]*A*X,
#'               parM[4,1] + parM[4,2]*A + parM[4,3]*X + parM[4,4]*A*X)+ mvrnorm(n , mu =c(0,0,0,0) , Sigma = matrix(c(5,-1,0,2,-1,6,1,0,0,1,4,3,2,0,3,7), nrow = 4))
#'
#'    Y <- rbinom(n,1,plogis(parY[1]*U + parY[2]*M[,1]+ parY[3]*M[,2] + parY[4]*M[,3] + parY[5]*M[,4] + parY[6]*X + parY[7]*M[,1]*X)) # p(Y|U,M,X)
#'
#'    data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)
#'
#'    # propensity score
#'    ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))
#'
#'    # mediator density ratio: p(M|a,X)/p(M|A,X)
#'    f.m.ratio.a1 <- function(i){
#'      dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*1 + parM[1,3]*X[i] + parM[1,4]*1*X[i],
#'                                  parM[2,1] + parM[2,2]*1 + parM[2,3]*X[i] + parM[2,4]*1*X[i],
#'                                  parM[3,1] + parM[3,2]*1 + parM[3,3]*X[i] + parM[3,4]*1*X[i],
#'                                  parM[4,1] + parM[4,2]*1 + parM[4,3]*X[i] + parM[4,4]*1*X[i]),sigma=matrix(c(5,-1,0,2,-1,6,1,0,0,1,4,3,2,0,3,7), nrow = 4))/dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*A[i] + parM[1,3]*X[i] + parM[1,4]*A[i]*X[i],
#'                                                                                                                                                                                         parM[2,1] + parM[2,2]*A[i] + parM[2,3]*X[i] + parM[2,4]*A[i]*X[i],
#'                                                                                                                                                                                         parM[3,1] + parM[3,2]*A[i] + parM[3,3]*X[i] + parM[3,4]*A[i]*X[i],
#'                                                                                                                                                                                         parM[4,1] + parM[4,2]*A[i] + parM[4,3]*X[i] + parM[4,4]*A[i]*X[i]),sigma=matrix(c(5,-1,0,2,-1,6,1,0,0,1,4,3,2,0,3,7), nrow = 4))
#'    }
#'    f.m.ratio.a0 <- function(i){
#'      dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*0 + parM[1,3]*X[i] + parM[1,4]*0*X[i],
#'                                  parM[2,1] + parM[2,2]*0 + parM[2,3]*X[i] + parM[2,4]*0*X[i],
#'                                  parM[3,1] + parM[3,2]*0 + parM[3,3]*X[i] + parM[3,4]*0*X[i],
#'                                  parM[4,1] + parM[4,2]*0 + parM[4,3]*X[i] + parM[4,4]*0*X[i]),sigma=matrix(c(5,-1,0,2,-1,6,1,0,0,1,4,3,2,0,3,7), nrow = 4))/dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*A[i] + parM[1,3]*X[i] + parM[1,4]*A[i]*X[i],
#'                                                                                                                                                                                         parM[2,1] + parM[2,2]*A[i] + parM[2,3]*X[i] + parM[2,4]*A[i]*X[i],
#'                                                                                                                                                                                         parM[3,1] + parM[3,2]*A[i] + parM[3,3]*X[i] + parM[3,4]*A[i]*X[i],
#'                                                                                                                                                                                         parM[4,1] + parM[4,2]*A[i] + parM[4,3]*X[i] + parM[4,4]*A[i]*X[i]),sigma=matrix(c(5,-1,0,2,-1,6,1,0,0,1,4,3,2,0,3,7), nrow = 4))
#'    }
#'
#'    m.ratio.a1 <- sapply(1:n, f.m.ratio.a1)
#'    m.ratio.a0 <- sapply(1:n, f.m.ratio.a0)
#'
#'    return(list(data = data,
#'                parA=parA,
#'                parU=parU,
#'                parM=parM,
#'                parY=parY,
#'                sd.U=sd.U,
#'                sigma.M=matrix(c(2, 1, 1, 3), nrow = 2),
#'                ps=ps,
#'                m.ratio.a1=m.ratio.a1,
#'                m.ratio.a0=m.ratio.a0))
#'  }
#'
#' }
#'
"binaryY_quadrivariateM"







#' @title A Sample Data for Analysis
#' @format ## `continuousY_binaryM`
#' A data frame with 500 rows and 5 columns:
#' \describe{
#'   \item{A}{Binary treatment variable}
#'   \item{M}{Binary mediator variable}
#'   \item{Y}{Normally distributed outcome variable}
#'   \item{X}{Uniformly distributed measure confounder within the range of (0,1)}
#'   \item{U}{Normally distributed unmeasured confounder}
#' }
#'
#' @examples
#' \donttest{
#' # data generated with the following Data Generating Process (DGP)
#' generate_data <- function(n,parA = c(0.3,0.2), parU=c(1,1,1,0), parM = c(-1,1,1,0), parY = c(1, 1, 1, 0), sd.U=1, sd.Y=1){ # change the parM to c(-1,1,1,0)
#'
#'  X <- runif(n, 0, 1) # p(X)
#'
#'  A <- rbinom(n, 1, (parA[1] + parA[2]*X)) # p(A|X)
#'
#'  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)
#'
#'  M <- rbinom(n,1,plogis(parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X)) # p(M|A,X)
#'
#'  Y <- parY[1]*U + parY[2]*M + parY[3]*X + parY[4]*M*X + rnorm(n, 0, sd.Y) # p(Y|U,M,X)
#'
#'  data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)
#'
#'  # propensity score
#'  ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))
#'
#'  # mediator density ratio: p(M|a,X)/p(M|A,X)
#'  m.ratio.a1 <- dbinom(M,1,plogis(parM[1] + parM[2]*1 + parM[3]*X + parM[4]*1*X))/dbinom(M,1,plogis(parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X))
#'  m.ratio.a0 <- dbinom(M,1,plogis(parM[1] + parM[2]*0 + parM[3]*X + parM[4]*0*X))/dbinom(M,1,plogis(parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X))
#'
#'
#'  return(list(data = data,
#'             parA=parA,
#'              parU=parU,
#'              parM=parM,
#'              parY=parY,
#'              sd.U=sd.U,
#'              sd.Y=sd.Y,
#'              ps=ps,
#'              m.ratio.a1=m.ratio.a1,
#'              m.ratio.a0=m.ratio.a0))
#'  }
#' }
"continuousY_binaryM"


#' @title A Sample Data for Analysis
#' @format ## `continuousY_continuousM`
#' A data frame with 500 rows and 5 columns:
#' \describe{
#'   \item{A}{Binary treatment variable}
#'   \item{M}{Normally distributed mediator variable}
#'   \item{Y}{Normally distributed outcome variable}
#'   \item{X}{Uniformly distributed measure confounder within the range of (0,1)}
#'   \item{U}{Normally distributed unmeasured confounder}
#' }
#' @examples
#' \donttest{
#' # data generated with the following Data Generating Process (DGP)
#' generate_data <- function(n,parA = c(0.3,0.2), parU=c(1,1,1,0), parM = c(1,1,1,0), parY = c(1, 1, 1, 0), sd.M=1, sd.U=1, sd.Y=1){
#'
#'X <- runif(n, 0, 1) # p(X)
#'
#'A <- rbinom(n, 1, (parA[1] + parA[2]*X)) # p(A|X)
#'
#'  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)
#'
#'  M <- parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X + rnorm(n,0,sd.M) # p(M|A,X)
#'
#'  Y <- parY[1]*U + parY[2]*M + parY[3]*X + parY[4]*M*X + rnorm(n, 0, sd.Y) # p(Y|U,M,X)
#'
#'  data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)
#'
#'  # propensity score
#'  ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))
#'
#'  # mediator density ratio: p(M|a,X)/p(M|A,X)
#'  m.ratio.a1 <- dnorm(M,parM[1] + parM[2]*1 + parM[3]*X + parM[4]*1*X,sd.M)/dnorm(M,parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X,sd.M)
#'  m.ratio.a0 <- dnorm(M,parM[1] + parM[2]*0 + parM[3]*X + parM[4]*0*X,sd.M)/dnorm(M,parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X,sd.M)
#'
#'  return(list(data = data,
#'            parA=parA,
#'            parU=parU,
#'            parM=parM,
#'            parY=parY,
#'            sd.U=sd.U,
#'            sd.Y=sd.Y,
#'            sd.M=sd.M,
#'            ps=ps,
#'            m.ratio.a1=m.ratio.a1,
#'            m.ratio.a0=m.ratio.a0))
#'  }
#'
#' }
#'
"continuousY_continuousM"


#' @title A Sample Data for Analysis
#' @format ## `continuousY_bivariateM`
#' A data frame with 500 rows and 5 columns:
#' \describe{
#'   \item{A}{Binary treatment variable}
#'   \item{M}{Bivariate Normally distributed mediator variable}
#'   \item{Y}{Normally distributed outcome variable}
#'   \item{X}{Uniformly distributed measure confounder within the range of (0,1)}
#'   \item{U}{Normally distributed unmeasured confounder}
#' }
#'
#' @examples
#' \donttest{
#' # data generated with the following Data Generating Process (DGP)
#' generate_data <- function(n,parA = c(0.3,0.2), parU=c(1,1,1,0), parM = matrix(c(1, 1, 1, 0,-1,-0.5,2,0), nrow = 2,byrow = T), parY = c(1, 1, -0.5,1, 0), sd.U=1, sd.Y=1){
#'
#'  ########################################################
#'  # M is bivariate normal with mean parameter be
#'  #  1  1.0    1    0
#'  # -1 -0.5    2    0
#'
#'  # and the variance covariance matrix be
#'  # 2 1
#'  # 1 3
#'  ########################################################
#'
#'  X <- runif(n, 0, 1) # p(X)
#'
#'  A <- rbinom(n, 1, (parA[1] + parA[2]*X)) # p(A|X)
#'
#'  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)
#'
#'  M <- cbind(parM[1,1] + parM[1,2]*A + parM[1,3]*X + parM[1,4]*A*X,
#'             parM[2,1] + parM[2,2]*A + parM[2,3]*X + parM[2,4]*A*X)+ mvrnorm(n , mu =c(0,0) , Sigma = matrix(c(2, 1, 1, 3), nrow = 2))
#'
#'  Y <- parY[1]*U + parY[2]*M[,1]+ parY[3]*M[,2] + parY[4]*X + parY[5]*M[,1]*X + rnorm(n, 0, sd.Y) # p(Y|U,M,X)
#'
#'  data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)
#'
#'  # propensity score
#'  ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))
#'
#'  # mediator density ratio: p(M|a,X)/p(M|A,X)
#'  f.m.ratio.a1 <- function(i){
#'    dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*1 + parM[1,3]*X[i] + parM[1,4]*1*X[i],
#'                                parM[2,1] + parM[2,2]*1 + parM[2,3]*X[i] + parM[2,4]*1*X[i]),sigma=matrix(c(2, 1, 1, 3), nrow = 2))/
#'      dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*A[i] + parM[1,3]*X[i] + parM[1,4]*A[i]*X[i], parM[2,1] + parM[2,2]*A[i] + parM[2,3]*X[i] + parM[2,4]*A[i]*X[i]),sigma=matrix(c(2, 1, 1, 3), nrow = 2))
#'  }
#'  f.m.ratio.a0 <- function(i){
#'    dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*0 + parM[1,3]*X[i] + parM[1,4]*0*X[i],
#'                                parM[2,1] + parM[2,2]*0 + parM[2,3]*X[i] + parM[2,4]*0*X[i]),sigma=matrix(c(2, 1, 1, 3), nrow = 2))/
#'      dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*A[i] + parM[1,3]*X[i] + parM[1,4]*A[i]*X[i], parM[2,1] + parM[2,2]*A[i] + parM[2,3]*X[i] + parM[2,4]*A[i]*X[i]),sigma=matrix(c(2, 1, 1, 3), nrow = 2))
#'  }
#'
#'  m.ratio.a1 <- sapply(1:n, f.m.ratio.a1)
#'  m.ratio.a0 <- sapply(1:n, f.m.ratio.a0)
#'
#'  return(list(data = data,
#'              parA=parA,
#'              parU=parU,
#'              parM=parM,
#'              parY=parY,
#'              sd.U=sd.U,
#'              sd.Y=sd.Y,
#'              sigma.M=matrix(c(2, 1, 1, 3), nrow = 2),
#'              ps=ps,
#'              m.ratio.a1=m.ratio.a1,
#'              m.ratio.a0=m.ratio.a0))
#'  }
#'
#' }
#'
"continuousY_bivariateM"


#' @title A Sample Data for Analysis
#' @format ## `continuousY_quadrivariateM`
#' A data frame with 500 rows and 5 columns:
#' \describe{
#'   \item{A}{Binary treatment variable}
#'   \item{M}{Quadrivariate Normally distributed mediator variable}
#'   \item{Y}{Normally distributed outcome variable}
#'   \item{X}{Uniformly distributed measure confounder within the range of (0,1)}
#'   \item{U}{Normally distributed unmeasured confounder}
#' }
#'
#' @examples
#' \donttest{
#' # data generated with the following Data Generating Process (DGP)
#' generate_data <- function(n,parA = c(0.3,0.2), parU=c(1,1,1,0), parM = matrix(c(1, 1, 1, 0,-1,-0.5,2,0,-1,2,1,0,1,0.5,-1,0), nrow = 4,byrow = T), parY = c(1, 1, -0.5,1,-1,1, 0), sd.U=1, sd.Y=1){
#'  ########################################################
#'  # M is now quadrivariate with parM be
#'  #  1  1.0    1    0
#'  # -1 -0.5    2    0
#'  # -1  2.     1.   0
#'  #  1 0.5  -1.  0
#'
#'  # and the variance covariance matrix be
#'  # 5 -1 0 2
#'  # -1  6 1 0
#'  # 0 1 4 3
#'  # 2 0 3 7
#'  ########################################################
#'
#'  X <- runif(n, 0, 1) # p(X)
#'
#'  A <- rbinom(n, 1, (parA[1] + parA[2]*X)) # p(A|X)
#'
#'  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)
#'
#'  M <- cbind(parM[1,1] + parM[1,2]*A + parM[1,3]*X + parM[1,4]*A*X,
#'             parM[2,1] + parM[2,2]*A + parM[2,3]*X + parM[2,4]*A*X,
#'             parM[3,1] + parM[3,2]*A + parM[3,3]*X + parM[3,4]*A*X,
#'             parM[4,1] + parM[4,2]*A + parM[4,3]*X + parM[4,4]*A*X)+ mvrnorm(n , mu =c(0,0,0,0) , Sigma = matrix(c(5,-1,0,2,-1,6,1,0,0,1,4,3,2,0,3,7), nrow = 4))
#'
#'  Y <- parY[1]*U + parY[2]*M[,1]+ parY[3]*M[,2] + parY[4]*M[,3] + parY[5]*M[,4] + parY[6]*X + parY[7]*M[,1]*X + rnorm(n, 0, sd.Y) # p(Y|U,M,X)
#'
#'  data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)
#'
#'  # propensity score
#'  ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))
#'
#'  # mediator density ratio: p(M|a,X)/p(M|A,X)
#'  f.m.ratio.a1 <- function(i){
#'    dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*1 + parM[1,3]*X[i] + parM[1,4]*1*X[i],
#'                                parM[2,1] + parM[2,2]*1 + parM[2,3]*X[i] + parM[2,4]*1*X[i],
#'                                parM[3,1] + parM[3,2]*1 + parM[3,3]*X[i] + parM[3,4]*1*X[i],
#'                                parM[4,1] + parM[4,2]*1 + parM[4,3]*X[i] + parM[4,4]*1*X[i]),sigma=matrix(c(5,-1,0,2,-1,6,1,0,0,1,4,3,2,0,3,7), nrow = 4))/dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*A[i] + parM[1,3]*X[i] + parM[1,4]*A[i]*X[i],
#'                                                                                                                                                                                       parM[2,1] + parM[2,2]*A[i] + parM[2,3]*X[i] + parM[2,4]*A[i]*X[i],
#'                                                                                                                                                                                       parM[3,1] + parM[3,2]*A[i] + parM[3,3]*X[i] + parM[3,4]*A[i]*X[i],
#'                                                                                                                                                                                       parM[4,1] + parM[4,2]*A[i] + parM[4,3]*X[i] + parM[4,4]*A[i]*X[i]),sigma=matrix(c(5,-1,0,2,-1,6,1,0,0,1,4,3,2,0,3,7), nrow = 4))
#'  }
#'  f.m.ratio.a0 <- function(i){
#'    dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*0 + parM[1,3]*X[i] + parM[1,4]*0*X[i],
#'                                parM[2,1] + parM[2,2]*0 + parM[2,3]*X[i] + parM[2,4]*0*X[i],
#'                                parM[3,1] + parM[3,2]*0 + parM[3,3]*X[i] + parM[3,4]*0*X[i],
#'                                parM[4,1] + parM[4,2]*0 + parM[4,3]*X[i] + parM[4,4]*0*X[i]),sigma=matrix(c(5,-1,0,2,-1,6,1,0,0,1,4,3,2,0,3,7), nrow = 4))/dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*A[i] + parM[1,3]*X[i] + parM[1,4]*A[i]*X[i],
#'                                                                                                                                                                                       parM[2,1] + parM[2,2]*A[i] + parM[2,3]*X[i] + parM[2,4]*A[i]*X[i],
#'                                                                                                                                                                                       parM[3,1] + parM[3,2]*A[i] + parM[3,3]*X[i] + parM[3,4]*A[i]*X[i],
#'                                                                                                                                                                                       parM[4,1] + parM[4,2]*A[i] + parM[4,3]*X[i] + parM[4,4]*A[i]*X[i]),sigma=matrix(c(5,-1,0,2,-1,6,1,0,0,1,4,3,2,0,3,7), nrow = 4))
#'  }
#'
#'  m.ratio.a1 <- sapply(1:n, f.m.ratio.a1)
#'  m.ratio.a0 <- sapply(1:n, f.m.ratio.a0)
#'
#'  return(list(data = data,
#'              parA=parA,
#'              parU=parU,
#'              parM=parM,
#'              parY=parY,
#'              sd.U=sd.U,
#'              sd.Y=sd.Y,
#'              sigma.M=matrix(c(2, 1, 1, 3), nrow = 2),
#'              ps=ps,
#'              m.ratio.a1=m.ratio.a1,
#'              m.ratio.a0=m.ratio.a0))
#'  }
#'
#' }
#'
"continuousY_quadrivariateM"



#' @title A Sample Data for Analysis
#' @format ## `continuousY_continuousM_10dX`
#' A data frame with 500 rows and 14 columns:
#' \describe{
#'   \item{A}{Binary treatment variable}
#'   \item{M}{Quadrivariate Normally distributed mediator variable}
#'   \item{Y}{Normally distributed outcome variable}
#'   \item{X.i}{10 Uniformly distributed measure confounder within the range of (0,1)}
#'   \item{U}{Normally distributed unmeasured confounder}
#' }
#'
#' @examples
#' \donttest{
#' # data generated with the following Data Generating Process (DGP)
#' generate_data <- function(n,parA = c(0.48, 0.07, 1.00, -1.00, -0.34, -0.12, 0.30, -0.35, 1.00, -0.10, 0.46, # linear terms
#'  0.33, 0.00, 0.45, 0.1, -0.32, -0.08, -0.2, 0.50, 0.50, -0.03)*0.1,# X^2 high order terms
#'
#'  parU=c(-2.0, -1.0, -1.0, 2.0, 3.0, 0.5, 3.0, 2.0, -1.0, 1.0, -3.0, 1.5, # linear terms
#'         -3.0, -2.0,  1.0,  3.0,  1.5), # A*X interactions
#'
#'  parM = c(3.0, 1.5, -1.5, -1.5, -1.0, -2.0, -3.0, -3.0, -1.5, 2.0, 1.5, 3.0, # linear terms
#'           1.5, 2.0, 0.5, 0.5, 3.0, # A*X interactions
#'           -0.2, -0.33, 0.5, 0.3, -0.5)*0.025, # X^2 high order terms
#'
#'  parY = c(1.0, -2.0, -3.0, -1.5, 1.0, 0.5, -2.0, 1.5, -2.0, -3.0, -3.0, -1.5, # linear terms
#'           -1.0, 0.5, 3.0, 1.5, 0.5, # M*X interactions
#'           3, # M^2 high order term
#'           1.0, 1.5, -2.0, 3.0, -1.0), # X^2 high order terms
#'  sd.M=1, sd.U=1, sd.Y=1){
#'
#'    # generate X from uniform distr
#'    X <- replicate(10, runif(n,0,1))
#'
#'    A <- rbinom(n, 1, plogis(parA[1] + rowSums(sweep(X, 2, parA[2:11], "*")) +  rowSums(sweep(X^2,2,parA[12:21],"*")) )) # p(A|X)
#'
#'    U <- parU[1] + parU[2]*A + rowSums(sweep(X, 2, parU[3:12], "*")) + rowSums(diag(A)%*%sweep(X[,1:5],2,parU[13:17],"*")) + rnorm(n,0,sd.U) # p(U|A,X)
#'
#'    M <- parM[1] + parM[2]*A + rowSums(sweep(X, 2, parM[3:12], "*")) + rowSums(diag(A)%*%sweep(X[,1:5],2,parM[13:17],"*")) + rowSums(sweep(X[,6:10]^2,2,parM[18:22],"*"))  + rnorm(n,0,sd.M) # p(M|A,X)
#'
#'    Y <- parY[1]*U + parY[2]*M + rowSums(sweep(X, 2, parY[3:12], "*")) + rowSums(diag(M)%*%sweep(X[,1:5],2,parY[13:17],"*")) + parY[18]*M^2 + rowSums(sweep(X[,6:10]^2,2,parY[19:23],"*")) + rnorm(n, 0, sd.Y) # p(Y|U,M,X)
#'
#'    data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)
#'
#'    # propensity score
#'    ps <- A*plogis(parA[1] + rowSums(sweep(X, 2, parA[2:11], "*")) +  rowSums(sweep(X^2,2,parA[12:21],"*")) )+(1-A)*(1-plogis(parA[1] + rowSums(sweep(X, 2, parA[2:11], "*")) +  rowSums(sweep(X^2,2,parA[12:21],"*")) ))
#'
#'    # mediator density ratio: p(M|a,X)/p(M|A,X)
#'    m.ratio.a1 <- dnorm(M,parM[1] + parM[2]*1 + rowSums(sweep(X, 2, parM[3:12], "*")) + rowSums(diag(rep(1,n))%*%sweep(X[,1:5],2,parM[13:17],"*")) + rowSums(sweep(X[,6:10]^2,2,parM[18:22],"*")) ,sd.M)/
#'      dnorm(M,parM[1] + parM[2]*A + rowSums(sweep(X, 2, parM[3:12], "*")) + rowSums(diag(A)%*%sweep(X[,1:5],2,parM[13:17],"*")) + rowSums(sweep(X[,6:10]^2,2,parM[18:22],"*")) ,sd.M)
#'    m.ratio.a0 <- dnorm(M,parM[1] + parM[2]*0 + rowSums(sweep(X, 2, parM[3:12], "*")) + rowSums(diag(rep(0,n))%*%sweep(X[,1:5],2,parM[13:17],"*")) + rowSums(sweep(X[,6:10]^2,2,parM[18:22],"*")) ,sd.M)/
#'      dnorm(M,parM[1] + parM[2]*A + rowSums(sweep(X, 2, parM[3:12], "*")) + rowSums(diag(A)%*%sweep(X[,1:5],2,parM[13:17],"*")) + rowSums(sweep(X[,6:10]^2,2,parM[18:22],"*")) ,sd.M)
#'
#'    return(list(data = data,
#'                X=X,
#'                parA=parA,
#'                parU=parU,
#'                parM=parM,
#'                parY=parY,
#'                sd.U=sd.U,
#'                sd.Y=sd.Y,
#'                sd.M=sd.M,
#'                ps=ps,
#'                m.ratio.a1=m.ratio.a1,
#'                m.ratio.a0=m.ratio.a0
#'    ))
#'  }
#'
#'
#' }
#'
"continuousY_continuousM_10dX"




#' @title A Sample Data for Testing Out the Performance of TMLE in Dealing with Weak Overlap
#' @format ## `continuousY_binaryM_weakoverlap`
#' A data frame with 500 rows and 5 columns:
#' \describe{
#'   \item{A}{Binary treatment variable}
#'   \item{M}{Binary mediator variable}
#'   \item{Y}{Normally distributed outcome variable}
#'   \item{X}{Uniformly distributed measure confounder within the range of (0,1)}
#'   \item{U}{Normally distributed unmeasured confounder}
#' }
#'
#' @examples
#' \donttest{
#' # data generated with the following Data Generating Process (DGP)
#' generate_data <- function(n,parA = c(0.001,0.998), parU=c(1,1,1,0), parM = c(-1,1,1,0), parY = c(1, 1, 1, 0), sd.U=1, sd.Y=1){
# X <- rbinom(n, 1, 0.5) # p(X)
#'
#'  X <- runif(n, 0, 1) # p(X)
#'
#'  A <- rbinom(n, 1, parA[1] + parA[2]*X) # p(A|X)
#'
#'  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)
#'
#'  M <- rbinom(n,1,plogis(parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X)) # p(M|A,X)
#'
#'  Y <- parY[1]*U + parY[2]*M + parY[3]*X + parY[4]*M*X + rnorm(n, 0, sd.Y) # p(Y|U,M,X)
#'
#'  data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)
#'
#'  # propensity score
#'  ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))
#'
#'  # mediator density ratio: p(M|a,X)/p(M|A,X)
#'  m.ratio.a1 <- dbinom(M,1,plogis(parM[1] + parM[2]*1 + parM[3]*X + parM[4]*1*X))/dbinom(M,1,plogis(parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X))
#'  m.ratio.a0 <- dbinom(M,1,plogis(parM[1] + parM[2]*0 + parM[3]*X + parM[4]*0*X))/dbinom(M,1,plogis(parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X))
#'
#'
#'  return(list(data = data,
#'              parA=parA,
#'              parU=parU,
#'              parM=parM,
#'              parY=parY,
#'              sd.U=sd.U,
#'              sd.Y=sd.Y,
#'              ps=ps,
#'              m.ratio.a1=m.ratio.a1,
#'              m.ratio.a0=m.ratio.a0))
#'  }
#' }
#'
"continuousY_binaryM_weakoverlap"


#' @title A Sample Data for Testing Out the Performance of TMLE in Dealing with Weak Overlap
#' @format ## `continuousY_continuousM_weakoverlap`
#' A data frame with 500 rows and 5 columns:
#' \describe{
#'   \item{A}{Binary treatment variable}
#'   \item{M}{Normally distributed mediator variable}
#'   \item{Y}{Normally distributed outcome variable}
#'   \item{X}{Uniformly distributed measure confounder within the range of (0,1)}
#'   \item{U}{Normally distributed unmeasured confounder}
#' }
#'
#' @examples
#' \donttest{
#' # data generated with the following Data Generating Process (DGP)
#' generate_data <- function(n,parA = c(0.001,0.998), parU=c(1,1,1,0), parM = c(1,1,1,0), parY = c(1, 1, 1, 0), sd.M=1, sd.U=1, sd.Y=1){
#'
#'  X <- runif(n, 0, 1) # p(X)
#'
#'  A <- rbinom(n, 1, (parA[1] + parA[2]*X)) # p(A|X)
#'
#'  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)
#'
#'  M <- parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X + rnorm(n,0,sd.M) # p(M|A,X)
#'
#'  Y <- parY[1]*U + parY[2]*M + parY[3]*X + parY[4]*M*X + rnorm(n, 0, sd.Y) # p(Y|U,M,X)
#'
#'  data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)
#'
#'  # propensity score
#'  ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))
#'
#'  # mediator density ratio: p(M|a,X)/p(M|A,X)
#'  m.ratio.a1 <- dnorm(M,parM[1] + parM[2]*1 + parM[3]*X + parM[4]*1*X,sd.M)/dnorm(M,parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X,sd.M)
#'  m.ratio.a0 <- dnorm(M,parM[1] + parM[2]*0 + parM[3]*X + parM[4]*0*X,sd.M)/dnorm(M,parM[1] + parM[2]*A + parM[3]*X + parM[4]*A*X,sd.M)
#'
#'  return(list(data = data,
#'              parA=parA,
#'              parU=parU,
#'              parM=parM,
#'              parY=parY,
#'              sd.U=sd.U,
#'              sd.Y=sd.Y,
#'              sd.M=sd.M,
#'              ps=ps,
#'              m.ratio.a1=m.ratio.a1,
#'              m.ratio.a0=m.ratio.a0))
#'  }
#' }
#'
"continuousY_continuousM_weakoverlap"


#' @title A Sample Data for Testing Out the Performance of TMLE in Dealing with Weak Overlap
#' @format ## `continuousY_bivariateM_weakoverlap`
#' A data frame with 500 rows and 5 columns:
#' \describe{
#'   \item{A}{Binary treatment variable}
#'   \item{M}{Bivariate Normally distributed mediator variable}
#'   \item{Y}{Normally distributed outcome variable}
#'   \item{X}{Uniformly distributed measure confounder within the range of (0,1)}
#'   \item{U}{Normally distributed unmeasured confounder}
#' }
#'
#' @examples
#' \donttest{
#' # data generated with the following Data Generating Process (DGP)
#' generate_data <- function(n,parA = c(0.001,0.998), parU=c(1,1,1,0), parM = matrix(c(1, 1, 1, 0,-1,-0.5,2,0), nrow = 2,byrow = T), parY = c(1, 1, -0.5,1, 0), sd.U=1, sd.Y=1){
#'
#'  ########################################################
#'  # M is bivariate normal with mean parameter be
#'  #  1  1.0    1    0
#'  # -1 -0.5    2    0
#'
#'  # and the variance covariance matrix be
#'  # 2 1
#'  # 1 3
#'  ########################################################
#'
#'  X <- runif(n, 0, 1) # p(X)
#'
#'  A <- rbinom(n, 1, (parA[1] + parA[2]*X)) # p(A|X)
#'
#'  U <- parU[1] + parU[2]*A + parU[3]*X + parU[4]*A*X + rnorm(n,0,sd.U) # p(U|A,X)
#'
#'  M <- cbind(parM[1,1] + parM[1,2]*A + parM[1,3]*X + parM[1,4]*A*X,
#'             parM[2,1] + parM[2,2]*A + parM[2,3]*X + parM[2,4]*A*X)+ mvrnorm(n , mu =c(0,0) , Sigma = matrix(c(2, 1, 1, 3), nrow = 2))
#'
#'  Y <- parY[1]*U + parY[2]*M[,1]+ parY[3]*M[,2] + parY[4]*X + parY[5]*M[,1]*X + rnorm(n, 0, sd.Y) # p(Y|U,M,X)
#'
#'  data <- data.frame(X=X, U=U, A=A, M=M, Y=Y)
#'
#'  # propensity score
#'  ps <- A*(parA[1] + parA[2]*X)+(1-A)*(1-(parA[1] + parA[2]*X))
#'
#'  # mediator density ratio: p(M|a,X)/p(M|A,X)
#'  f.m.ratio.a1 <- function(i){
#'    dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*1 + parM[1,3]*X[i] + parM[1,4]*1*X[i],
#'                                parM[2,1] + parM[2,2]*1 + parM[2,3]*X[i] + parM[2,4]*1*X[i]),sigma=matrix(c(2, 1, 1, 3), nrow = 2))/dmvnorm(x=c(0,0), mean=cbind(parM[1,1] + parM[1,2]*A[i] + parM[1,3]*X[i] + parM[1,4]*A[i]*X[i],
#'                                                                                                                                                                 parM[2,1] + parM[2,2]*A[i] + parM[2,3]*X[i] + parM[2,4]*A[i]*X[i]),sigma=matrix(c(2, 1, 1, 3), nrow = 2))
#'  }
#'  f.m.ratio.a0 <- function(i){
#'    dmvnorm(x=M[i,], mean=cbind(parM[1,1] + parM[1,2]*0 + parM[1,3]*X[i] + parM[1,4]*0*X[i],
#'                                parM[2,1] + parM[2,2]*0 + parM[2,3]*X[i] + parM[2,4]*0*X[i]),sigma=matrix(c(2, 1, 1, 3), nrow = 2))/dmvnorm(x=c(0,0), mean=cbind(parM[1,1] + parM[1,2]*A[i] + parM[1,3]*X[i] + parM[1,4]*A[i]*X[i],
#'                                                                                                                                                                 parM[2,1] + parM[2,2]*A[i] + parM[2,3]*X[i] + parM[2,4]*A[i]*X[i]),sigma=matrix(c(2, 1, 1, 3), nrow = 2))
#'  }
#'
#'  m.ratio.a1 <- sapply(1:n, f.m.ratio.a1)
#'  m.ratio.a0 <- sapply(1:n, f.m.ratio.a0)
#'
#'  return(list(data = data,
#'              parA=parA,
#'              parU=parU,
#'              parM=parM,
#'              parY=parY,
#'              sd.U=sd.U,
#'              sd.Y=sd.Y,
#'              sigma.M=matrix(c(2, 1, 1, 3), nrow = 2),
#'              ps=ps,
#'              m.ratio.a1=m.ratio.a1,
#'              m.ratio.a0=m.ratio.a0))
#'  }
#'
#' }
#'
"continuousY_bivariateM_weakoverlap"
