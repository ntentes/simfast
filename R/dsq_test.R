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
  yhat <- simfastobj$yhat
  weights <- simfastobj$weights
  ylength <- length(weights)
  offset <- simfastobj$offset
  fam  <- simfastobj$family

  ## Remove offset if included
  if (!is.null(simfastobj$offset)) {
    oldweights <- weights
    weights <- oldweights * fam$linkinv(offset)
    comment(weights) <- 'offset'
    oldy <- y
    oldyhat <- yhat
    y    <- fam$linkfun(y) - offset
    y    <- fam$linkinv(y)
    yhat <- fam$linkfun(yhat) - offset
    yhat <- fam$linkinv(yhat)
  }

  ## Extract alpha-hat
  alphahat <- simfastobj$alphahat

  ## Test statistic
  Dsq_n <- sum(weights * (y - yhat) ^ 2)

  ## Simulate S_d-1
  alphas		<-	matrix(stats::rnorm(d*N), N, d)
  alphas		<-	alphas/apply(alphas, 1, function(x) sqrt(sum(x^2)))

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

  ## calculate dispersion parameter
  famname <- fam$family
  if (famname %in% c("poisson",
                      "binomial")) {
    disp <- 1
  } else {
    disp <- ((y - yhat) ^ 2) / ylength
  }

  ## Compute upper quantile of theoretical stat
  dsq.stat		<-	rep(0, B)

  pb <- utils::txtProgressBar(min = 0, max = 100, style = 3)
  pbseq <- round(seq(from = 0, to = 100, by = 100/(B-1)))

  for (i in 1:B) {
    Ti          <-  rnorm(n = ylength, mean = 0, sd = sqrt(1/lambdan))
    res				  <-	apply(X = alphas, FUN = dsq, MARGIN = 1, newresp = Ti, wts = weights, preds = x)
    indmin      <-  which.min(res)
    dsq.stat[i] <-  res[indmin]

    utils::setTxtProgressBar(pb, pbseq[i])

  }

  close(pb)

  empcdf   <- ecdf(disp * dsq.stat)

  ## what is ll?
  pval  <-  1 - empcdf(ll)

}
