#' Stochastic Search Algorithm for Isotonic Single-Index Regression Fitting (Internal)
#'
#' @param x matrix of predictor values
#' @param y vector of response values
#' @param nn vector of positive integer weights
#' @param family \code{\link{family}} object passed by \code{simfast}
#' @param B positive integer, number of index vectors to test
#' @param k positive integer, number of alpha tests per iteration
#' @param kappa0 positive integer, kappa0
#' @param tol numeric, convergence tolerance
#' @param max.iter numeric, maximum number of iterations allowed
#' @param print logical, prints number and tolerance at each iteration step
#'
#'
#' @return list including N(>=1) alpha index estimates in matrix form, N maximum
#'     likelihood estimates of the response in matrix form, and N single index
#'     estimates in matrix form
#'
#' @noRd
#' @keywords internal
#'
search.mle		<-	function(x, y, nn, family='gaussian', B=10000, k=100, kappa0=100,
                        tol=1e-10, max.iter=20, print=TRUE){

  m			<-	length(x[,1])
  d			<-	length(x[1,])

  alpha		<-	matrix(stats::rnorm(d*B), B, d)
  alpha		<-	alpha/apply(alpha,1, function(x) sqrt(sum(x^2)))

  loglik		<-	rep(0, B)
  for(i in 1:B){
    res				<-	find.pMLE(alpha[i,], y , x, nn, family)
    loglik[i]		<-	res$ll
  }
  loc			<-	which.max(loglik)
  alphahat		<-	alpha[loc,]
  llhat		<-	loglik[loc]

  check		<-	10
  iter    <-  0

  for(j in 1:max.iter){

    if(print==TRUE) print(check)
    if(check < tol) break

    if(print==TRUE) print(j)

    locs		<-	getkmax(loglik, k)

    BB			<-	B/k
    alpha2		<-	matrix(0,BB*k, d)
    kappa		<-	kappa0*(j^2)

    for(i in 1:k){

      mm						<-	alpha[locs[i],]
      temp					<-	movMF::rmovMF(BB-1, theta=kappa*mm)
      alpha2[(i-1)*BB+1:BB,]	<-	rbind(mm, temp)

    }

    loglik2		<-	rep(0, B)
    for(i in 1:B){
      res				<-	find.pMLE(alpha2[i,], y , x, nn, family)
      if(is.nan(res$ll)){
      }
      loglik2[i]		<-	res$ll
    }
    llhat2		<-	max(loglik2)
    check		<-	(llhat2-llhat)/abs(llhat)

    if(print==TRUE) print(c(llhat, llhat2, check))

    loglik		<-	loglik2
    alpha		<-	alpha2
    llhat		<-	llhat2
    iter    <- iter + 1

  }

  loc			<-	which(loglik==max(loglik))
  alphahat	<-	alpha[loc,]

  if(length(loc)>1)	{
    alphahat1	<-	alphahat[1,]
    maxll  <-  loglik[loc[1]]
  } else {
    alphahat1	<-	alphahat
    maxll  <-  loglik[loc]
  }
  mle		 <-	 find.pMLE(alphahat1, y , x, nn, family)
  zhat	 <-	 x %*% alphahat1
  mle1   <-  mle$mle
  zhat   <-  as.vector(zhat)
  zord   <-  order(zhat)
  mle1   <-  mle1[order(zord)]
  k      <-  length(loc)

  return(list('alphahat' = alphahat, 'mle' = mle1, 'zhat' = zhat,
              'k' = k, 'maxll' = maxll, 'iter' = iter, 'tol' = check))

}


#' Exact Algorithm for Isotonic Single-Index Regression Fitting
#' with 2-D Predictor Space (Internal)
#'
#' @param x matrix of predictor values
#' @param y vector of reponse values
#' @param nn vector of positive integer weights
#' @param family \code{\link{family}} object passed by \code{simfast}
#'
#' @return list including alpha index estimate in vector form, maximum
#'     likelihood estimate of the response in vector form, and single index
#'     estimate in vector form
#'
#' @noRd
#' @keywords internal
#'
find.mle	<-	function(x, y, nn, family='gaussian'){

  m		<-	length(x[,1])

  pbseq <- round(seq(from = 0, to = 10, by = 10/(m-1)))

  total <- 100
  pb <- utils::txtProgressBar(min = 0, max = total, style = 3)

  theta	<-	matrix(NA,m,m)
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      d			<-	x[i,]-x[j,]
      d			<-	d/sqrt(sum(d^2))
      z			<-	c(-d[2], d[1])
      theta[i,j]	<-	my.inv(z)
    }
    utils::setTxtProgressBar(pb, pbseq[i])
  }

  phi		<-	as.vector(theta)
  phi		<-	phi[!is.na(phi)]
  phi		<-	c(phi, (phi+pi) %% (2*pi))
  phi		<-	sort(phi)
  phi		<-	unique(phi)

  beta		<-	phi+c(phi[2:length(phi)], phi[1]+2*pi)
  beta		<-	(beta/2) %% (2*pi)

  loglik	<-	rep(0, length(beta))

  pbseq <- round(seq(from = 11, to = 90, by = (90-11)/length(beta)))

  for(i in 1:length(beta)){
    alpha		<-	c(cos(beta[i]), sin(beta[i]))
    res			<-	find.pMLE(alpha, y , x, nn, family)
    loglik[i]	<-	res$ll
    utils::setTxtProgressBar(pb, pbseq[i])
  }

  loc			<-	which(loglik==max(loglik))
  thetahat	<-	beta[loc]
  thetahat	<-	thetahat %% (2*pi)
  alphahat	<-	matrix(c(cos(thetahat), sin(thetahat)), ncol=2)

  phat	<-	matrix(0, nrow=length(loc), m)
  ordhat	<-	matrix(0, nrow=length(loc), m)
  llvals  <-  rep(0, length(loc))

  pbseq2 <- round(seq(from = 91, to = 99, by = (99-91)/length(loc)))

  for(i in 1:length(loc)){
    alpha		<-	alphahat[i,]
    mle			<-	find.pMLE(alpha, y , x, nn, family)
    phat[i,]	<-	mle$mle
    ordhat[i,]	<-	mle$ord
    llvals[i]   <-  mle$ll
    utils::setTxtProgressBar(pb, pbseq2[i])
  }

  if(length(loc)>1){
    alphahat1	<-	alphahat[1,]
    zhat		<-	x %*% alphahat1
  }
  else{
    alphahat1	<-	alphahat
    zhat		<-	x %*% t(alphahat1)
  }
  mle	 	 <-	 phat[1,]
  zhat   <-  as.vector(zhat)
  zord   <-  order(zhat)
  mle    <-  mle[order(zord)]
  k      <-  length(loc)
  maxll  <-  max(llvals)

  utils::setTxtProgressBar(pb, 100)
  close(pb)

  return(list('alphahat' = alphahat, 'mle' = mle, 'zhat' = zhat,
              'k' = k, 'maxll' = maxll, 'iter' = 1, tol = 0))

}



