
#' Log Likelihood Function for Gaussian Family (Internal)
#'
#' @param mle vector, maximum likelihood estimator
#' @param y vector of numeric response values
#' @param nn vector of positive integer weights
#'
#' @return numeric value of log likelihood function
#'
#'
find.loglik.gaussian	<-	function(mle, y, nn) {
  ll	<-	sum(nn * stats::dnorm(y, mean = mle, log = TRUE))
  return(ll)
}


#' Log Likelihood Function for Binomial Family (Internal)
#'
#' @param mle vector, maximum likelihood estimator
#' @param y vector of binary response values
#' @param nn vector of positive integer weights
#'
#' @return numeric value of log likelihood function
#'
#'
find.loglik.binomial	<-	function(mle, y, nn) {
  newy <- y * nn
  if (!isTRUE(all.equal(newy, round(newy)))) {
    warning("y*weights are not all within integer tolerance,
            values of y*weights will be rounded.")
    newy <- round(newy)
  }
  ll <- sum(stats::dbinom(x = newy, size = nn, prob = mle, log = TRUE))
  return(ll)
}


#' Log Likelihood Function for Poisson Family (Internal)
#'
#' @param mle vector, maximum likelihood estimator
#' @param y vector of positive integer response values
#' @param nn vector of positive integer weights
#'
#' @return numeric value of log likelihood function
#'
#'
find.loglik.poisson   <-  function(mle, y, nn) {
  if (!isTRUE(all.equal(y, round(y)))) {
    warning("Response values are not all within integer tolerance,
            responses will be rounded.")
    newy <- round(y)
  }
  ll <- sum(nn * stats::dpois(x = newy, lambda = mle, log = TRUE))
  return(ll)
}


#' Log Likelihood Function for Gamma with Unknown Shape Family (Internal)
#'
#' @param mle vector, maximum likelihood estimator
#' @param y vector of positive response values
#' @param nn vector of positive integer weights
#'
#' @return numeric value of log likelihood function
#'
#'
find.loglik.gamma   <-  function(mle, y, nn) {
  if (!isTRUE(all.equal(y, abs(y)))) {
    stop("Response values are not all positive, regression cannot continue.")
  }
  ll <- sum(nn * stats::dgamma(x = y, shape = mle, log = TRUE))
  return(ll)
}


###############################################################################
###############################################################################


#' Selects and Runs Likelihood for Given Distribution Family (Internal)
#'
#' @param mle vector, maximum likelihood estimator
#' @param y vector of numeric response values
#' @param nn vector of positive integer weights
#' @param family \code{\link{family}} object passed by \code{simfast}
#'
#' @return numeric value of log likelihood function
#'
find.loglik <- function(mle, y, nn, family) {
  if (family == 'gaussian') {
    find.loglik.gaussian(mle, y, nn)
  } else if (family == 'binomial') {
    find.loglik.binomial(mle, y, nn)
  } else if (family == 'poisson') {
    find.loglik.poisson(mle, y, nn)
  } else if (family == 'Gamma') {
    find.loglik.gamma(mle, y, nn)
  } else {
    stop("No likelihood defined: please ensure valid family is selected.")
  }
}
