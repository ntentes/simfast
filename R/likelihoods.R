
#' Log Likelihood Function for Gaussian Family with Indentity Link (Internal)
#'
#' @param mle vector of maximum likelihood estimator
#' @param y vector of numeric response values
#' @param nn vector of positive integer weights
#'
#' @return numeric value of log likelihood function
#'
#'
find.loglik.gaussian	<-	function(mle, y, nn){
  ll	<-	-sum((y-mle)^2)
  return(ll)
}


#' Log Likelihood Function for Binomial Family with Logistic Link (Internal)
#'
#' @param mle vector of maximum likelihood estimator
#' @param y vector of binary response values
#' @param nn vector of positive integer weights
#'
#' @return numeric value of log likelihood function
#'
#'
find.loglik.logit	<-	function(mle, y, nn){
  ll1			<-	nn*y*log(mle)
  ll1[y==0]	<-	0
  ll2			<-	nn*(1-y)*log(1-mle)
  ll2[y==1]	<-	0
  ll			<- 	sum(ll1+ll2)/sum(nn)
  return(ll)
}


###############################################################################
###############################################################################


#' Selects and Runs Likelihood for Given Distribution Family (Internal)
#'
#' @param mle vector of maximum likelihood estimator
#' @param y vector of numeric response values
#' @param nn vector of positive integer weights
#' @param family character string matching the error distribution and link function
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
