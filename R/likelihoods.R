
#' Log Likelihood Function for Gaussian Family (Internal)
#'
#' @param mle vector of maximum likelihood estimator
#' @param y vector of numeric response values
#' @param nn vector of positive integer weights
#'
#' @return numeric value of log likelihood function
#'
#'
find.loglik.gaussian	<-	function(mle, y, nn){
  ll	<-	sum(nn*dnorm(y, mean = mle, log = TRUE))
  return(ll)
}


#' Log Likelihood Function for Binomial Family (Internal)
#'
#' @param mle vector of maximum likelihood estimator
#' @param y vector of binary response values
#' @param nn vector of positive integer weights
#'
#' @return numeric value of log likelihood function
#'
#'
find.loglik.binomial	<-	function(mle, y, nn){
  newy <- y * nn
  if (!isTRUE(all(newy == floor(newy)))) {
    warning("y * weights are not all integers, responses will be rounded.")
    newy <- round(newy)
  }
  ll <- sum(dbinom(x = newy, size = nn, prob = mle, log = TRUE))
}


###############################################################################
###############################################################################


#' Selects and Runs Likelihood for Given Distribution Family (Internal)
#'
#' @param mle vector of maximum likelihood estimator
#' @param y vector of numeric response values
#' @param nn vector of positive integer weights
#' @param family character string matching the error distribution or link function
#'
#' @return numeric value of log likelihood function
#'
find.loglik <- function(mle, y, nn, family){
  if (family == 'gaussian'){
    find.loglik.gaussian(mle, y, nn)
  } else if (family %in% c('logit', 'binomial')){
    find.loglik.logit(mle, y, nn)
  } else {
    stop("No likelihood defined: please ensure valid family is selected.")
  }
}
