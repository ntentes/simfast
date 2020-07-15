dsq.test <- function(simfastobj, N = 10000, B = 10000) {
  if (class(simfastobj) !=  "simfast") {
    stop("Goodness of fit test only valid with a simfast object.")
  }
  if (is.null(simfastobj$y)) {
    stop("Please refit the simfast model using argument returndata = TRUE")
  }

  ## extract response and weights
  x <- simfastobj$x
  y <- simfastobj$y
  weights <- simfastobj$weights

  ## Remove offset if included
  if (!is.null(simfastobj$offset)) {
    oldweights <- weights
    weights <- oldweights * linkinv(offset)
    comment(weights) <- 'offset'
    oldy <- y
    fam  <- simfastobj$family
    y    <- fam$linkfun(y) - offset
    y    <- fam$linkinv(y)
  }

  ## Extract alpha-hat
  alphahat <- simfastobj$alphahat

  ## Test statistic


  ## Simulate S_d-1
  alphas		<-	matrix(stats::rnorm(d*N), N, d)
  alphas		<-	alphas/apply(alphas, 1, function(x) sqrt(sum(x^2)))

  ## Compute upper quantile of theoretical stat
  Dsq <- NULL

  ## Helper function to calculate D^2 stat
  dsq <- function(a, newresp, wts, preds) {
    a      <- matrix(a, ncol=1)
    ord    <- as.vector(order(preds%*%a))
    isor   <- Iso::pava(newresp[ord], wts[ord])
    lambda <- wts/sum(wts)

    return(sum(lambda[ord] * (isor[ord] - newresp[ord])^2))
  }

  ## Simulate S^d-1
  d         <-  length(simfastobj$alphahat)
  lambdan   <-  weights/sum(weights)

  dsq.stat		<-	rep(0, B)

  for (i in 1:B) {
    Ti          <-  rnorm(n = length(lambdan), mean = 0, sd = sqrt(1/lambdan))
    res				  <-	apply(X = alphas, FUN = dsq, MARGIN = 1, newresp = Ti, wts = weights, preds = x)
    indmin      <-  which.min(res)
    dsq.stat[i] <-  res[indmin]
  }

  empcdf   <- ecdf(dsq.stat)

  pval  <-  1 - empcdf(ll)

}
