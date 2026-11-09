# Implementations of Sparse Identification of Nonlinear Dynamics

# SINDy algorithms
# x - predictor/feature
# y - response/output
# iter - max number of iteration
# threshold - coefficients below threshold will be set to 0
sindy_alg <- function(x, y, iter=100, threshold=.5, digits=4, ...){
  fitted_coefs <- rep(Inf, ncol(x))
  x <- as.matrix(x)
  y <- as.matrix(y)

  for(i in 1:iter){
    active_features <- which(fitted_coefs != 0)
    if(length(active_features) == 0) stop("no active feature found")

    # fit model with Lasso regression
    # mod_x <- cv.glmnet(x[,active_features, drop=FALSE], y,
    #             family="gaussian",
    #             alpha=1, # use Lasso regression
    #             intercept = FALSE,
    #             ...)

    # threshold coef
    # new_coef <- coef(mod_x, s = "lambda.min")[-1]

    # fit with least square then threshold
    new_coef <- qr.solve(x[,active_features, drop=FALSE], y)

    # when there is no longer a need for threshold -> stop iteration
    if(sum(new_coef!=0 & abs(new_coef) < threshold) == 0) {
      fitted_coefs[active_features] <- new_coef
      break
    }

    new_coef[which(abs(new_coef) < threshold)] <- 0

    fitted_coefs[active_features] <- new_coef
  }

  # build the formula interms of feature
  terms <- colnames(x)[which(fitted_coefs != 0)]
  coefs <- fitted_coefs[which(fitted_coefs != 0)]
  form <- paste0(
    # sign
    if_else(as.numeric(coefs)<0, "- ", "+ "),
    # coef magnitude
    round(abs(coefs), digits=digits),
    # term
    terms,
    sep="",
    collapse = " "
  )

  list(
    "coefs" = fitted_coefs,
    "formula" = form
  )
}
