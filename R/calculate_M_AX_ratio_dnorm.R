calculate_M_AX_ratio_dnorm <- function(a,M,A,X){ # A is a vector, M and X are data frame

  num_columns <- ncol(M)
  num_rows <- nrow(M)

  # Initialize an empty matrix to store the fitted coefficients
  fit.parM <- matrix(NA, nrow = num_columns, ncol = 2 + ncol(X)) # ncol: 1 for the intercept + 1 for A + ncol(X) for X
                                                                 # currently only support linear models

  # Initialize vectors to store errors
  model_errors <- matrix(NA, nrow = num_rows, ncol = num_columns)

  # fitting glm/ln along each column of M

  for (i in seq_along(M)) {

    if (all(M[, i] %in% c(0,1))) { # binary columns

      # For binary columns, use glm
      model <- glm(M[, i] ~ . , data=data.frame(A,X))

      # Store model errors
      model_errors[, i] <- M[, i] - predict(model, type="response")

    } else { # continuous columns

      # For continuous columns, use lm
      model <- lm(M[, i] ~ . , data=data.frame(A,X))

      # Store model errors
      model_errors[, i] <- M[, i] - predict(model)
    }

    # Store fitted coefficients
    fit.parM[i, ] <- coef(model)


  }

  # Calculate variance-covariance matrix
  varcov <- cov(data.frame(model_errors))

  # Define a function for the ratio calculation
  f.m.ratio.a <- function(j) {
    dmvnorm(
      x = M[j, ],
      mean = rowSums(cbind(fit.parM[,1], # intercept
                   fit.parM[,2]*a, # coeffA*A
                   fit.parM[,3:ncol(fit.parM)] %*% t(X[j, ]))), # coeffX*X
      sigma = varcov
    ) / dmvnorm(
      x = M[j, ],
      mean = rowSums(cbind(fit.parM[,1], # intercept
                   fit.parM[,2]*A[j], # coeffA*A
                   fit.parM[,3:ncol(fit.parM)] %*% t(X[j, ]))), # coeffX*X
      sigma = varcov
    )
  }

  # Apply the function to each row of M
  M.AXratio <- sapply(1:num_rows, f.m.ratio.a)

  return(M.AXratio)
}
