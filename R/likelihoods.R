
#' Log Likelihood Function for Gaussian Family (Internal)
#'
#' @param mle vector, maximum likelihood estimator
#' @param y vector of numeric response values
#' @param nn vector of positive integer weights
#'
#' @return numeric value of log likelihood function
#'
#' @noRd
#' @keywords internal
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
#' @noRd
#' @keywords internal
#'
#'
find.loglik.binomial	<-	function(mle, y, nn, eps) {
  ll1			<-	nn*y*log(mle)
  ll1[y==0]	<-	0
  ll2			<-	nn*(1-y)*log(1-mle)
  ll2[y==1]	<-	0
  ll			<- 	sum(ll1+ll2)
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
#' @noRd
#' @keywords internal
#'
#'
find.loglik.poisson   <-  function(mle, y, nn, eps) {
  infnt <- mle == 0
  ll <- nn * (y * log(mle) - mle)
  ll[infnt] <- 0
  ll <- sum(ll)
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
#' @noRd
#' @keywords internal
#'
find.loglik.gamma   <-  function(mle, y, nn) {
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
#' @noRd
#' @keywords internal
#'
find.loglik <- function(mle, y, nn, family, eps) {
  if (family == 'gaussian') {
    find.loglik.gaussian(mle, y, nn)
  } else if (family == 'binomial') {
    find.loglik.binomial(mle, y, nn, eps)
  } else if (family == 'poisson') {
    find.loglik.poisson(mle, y, nn, eps)
  } else if (family == 'Gamma') {
    find.loglik.gamma(mle, y, nn)
  } else {
    stop("No likelihood defined: please ensure valid family is selected.")
  }
}
