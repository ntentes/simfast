#' Fitting isotonic generalized single-index regression models via maximum likelihood
#'
#' @param x numeric \code{matrix} of observed predictor values, with \code{n}
#'     rows and \code{d} columns
#' @param y numeric \code{vector} of observed response values, of length \code{n}.
#'     For \code{'binomial'} response, vector can be binary values or a vector of
#'     proportions, as long as a proper weight vector (a vector of the denominators
#'     of the proportions, see details) is included. Categorical vectors (character
#'     strings or factors) will automatically be translated into logical vector with
#'     the baseline factor level a 'success' (takes value \code{1}).
#' @param weights optional: vector of positive integer weights, with length
#'     \code{n}. Takes default value \code{NULL} which uses equal weights.
#' @param family character string naming the error distribution to be used in
#'     the model, \code{'gaussian'} by default. Other options include
#'     \code{'binomial'} (with logit link, equivalent to calling \code{'logit'}).
#'     Other options coming soon.
#' @param returndata optional boolean that when \code{TRUE} (the default value)
#'     returns the predictor matrix and response vector in the simfast object.
#' @param method when \code{x} has \code{d=2} columns, method can take \code{'exact'}
#'     argument, which uses an exact omptimization method instead of a stochastic
#'     search. If \code{d} does not equal \code{2}, \code{simfast_m} will give a
#'     warning and automatically continue with a stochastic search (the default
#'     method, \code{method = 'stochastic'}).
#' @param multialpha optional boolean, if \code{TRUE}, will return more than one
#'     \code{alpha} vector if available (see Value section and Details).
#' @param B positive integer, sets number of index vectors to try when maximizing
#'     the likelihood
#' @param k positive integer, less than \code{B}
#' @param kappa0 positive integer, initial value of kappa
#' @param tol numeric, sets tolerance for convergence for \code{method = 'stochastic'}.
#'     Will give value of \code{0} if \code{'exact'} is used.
#' @param max.iter positive integer limiting number of iterations for
#'     \code{method = 'stochastic'}
#'
#' @return an object of class \code{simfast}.
#' @export
#'
#' @examples
simfast_m <- function(x, y, weights = NULL, family = 'gaussian', returndata = TRUE,
                           method = 'stochastic', multialpha = FALSE, B = 10000,
                           k = 100, kappa0 = 100, tol = 1e-10, max.iter = 20) {
  if (NCOL(y) == 1) {
    if (is.character(y)) y <- factor(y)
    if (is.factor(y)) y <- y != levels(y)[1L]
  } else {
    stop("Response 'y' must be a single column.")
  }
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  print(family)
  print(typeof(family))
  famname <- as.character(substitute(family))[1]
  print(famname)
  if (!famname %in% c('gaussian', 'binomial', 'poisson', 'Gamma')){
    stop('Simfast does not support specified family')
  }
  if (is.null(family$family)) {
    print(family)
    stop("Specified family cannot be found.")
  }
  linkinv <- family$linkinv
  if (!is.function(linkinv))
    stop("Cannot find specified link function")
  x <- as.matrix(x)
  if (!is.numeric(x)){
    stop('Predictor matrix must be numeric.')
  }
  if (method == 'exact') {
    if (NCOL(x) == 2) {
      fit <- find.mle(x, y, nn = weights, family = famname)
    } else {
      warning('Predictor dimension does not equal 2: fit will use stochastic method.')
      method <- 'stochastic'
      fit <- search.mle(x, y, nn = weights, family = famname, B = B,
                        k = k, kappa0 = kappa0, tol = tol,
                        max.iter = max.iter, print = FALSE)
    }
  } else {
    fit <- search.mle(x, y, nn = weights, family = famname, B = B,
                      k = k, kappa0 = kappa0, tol = tol,
                      max.iter = max.iter, print = FALSE)
  }
  zhat <- fit$zhat
  lldev <- -2*fit$maxll
  iter <- fit$iter
  tol  <- fit$tol
  yhat <- firstrow(fit$mle)
  indexvals <- firstrow(fit$zhat)
  pdim <- NCOL(x)

  obj <- list('yhat' = yhat, 'indexvals' = indexvals, 'weights' = weights,
              'family' = family, 'link' = linkinv, 'deviance' = lldev,
              'tol' = tol, 'iter' = iter, 'method' = method, 'model' = NULL,
              'intercept' = NULL)

  if (multialpha == TRUE) {
    obj <- append(list("alphahat" = fit$alphahat), obj)
  } else {
    obj <- append(list("alphahat" = firstrow(fit$alphahat)), obj)
  }

  if (returndata == TRUE){
    obj <- append(list('x' = x, 'y' = y), obj)
  }

  class(obj) <- 'simfast'
  return(obj)
}


#' Fitting isotonic generalized single-index regression models via maximum likelihood
#' with formula support
#'
#' @param formula an object of class \code{\link{formula}}, which is a symbolic
#'     description of the model to be fitted. By default, intercepts are NOT
#'     included, so change argument \code{intercept = TRUE} to include one. NOTE: When
#'     including categorical predictors, be sure that you set
#'     \code{options('contrasts')} in your global options to a desired setting. This
#'     function was designed with the \code{'contr.treatment'} option in mind.
#'     For \code{'binomial'} response, the vector can be binary values or a vector of
#'     proportions, as long as a proper weight vector (a vector of the denominators
#'     of the proportions, see details) is provided. Categorical vectors (character
#'     strings or factors) will automatically be translated into logical vector with
#'     the baseline factor level a 'success' (takes value \code{1}).
#' @param data optional data frame (or object coercible by \code{\link{as.data.frame}})
#'     to a data frame) containing the variables in the model. Variables are taken
#'     from \code{environment(formula)} if not found in \code{data}.
#' @param intercept optional boolean, if \code{FALSE} (the default value), then the
#'     model given by the formula does not include an intercept value (even when
#'     including a 1, for example: \code{z ~ 1 + x + y} will only include columns for
#'     \code{x} and \code{y}).
#' @param weights optional: vector of positive integer weights, with length
#'     \code{n}. Takes default value \code{NULL} which uses equal weights.
#' @param family character string naming the error distribution to be used in
#'     the model, \code{'gaussian'} by default. Other options include
#'     \code{'binomial'} (with logit link, equivalent to calling \code{'logit'}).
#'     Other options coming soon.
#' @param returnmodel optional boolean that when \code{TRUE} (the default value)
#'     attaches the \code{\link{model.frame}} object to the simfast object. Leave
#'     as \code{TRUE} to properly use the \code{\link{predict}} function.
#' @param returndata optional boolean that when \code{TRUE} (the default value)
#'     returns the predictor matrix and response vector in the simfast object.
#' @param method when \code{x} has \code{d=2} columns, method can take \code{'exact'}
#'     argument, which uses an exact omptimization method instead of a stochastic
#'     search. If \code{d} does not equal \code{2}, \code{simfast} will give a
#'     warning and automatically continue with a stochastic search (the default
#'     method, \code{method = 'stochastic'}).
#' @param multialpha optional boolean, if \code{TRUE}, will return more than one
#'     \code{indexvec} if available (see Value section and Details).
#' @param B positive integer, sets number of index vectors to try when maximizing
#'     the likelihood
#' @param k positive integer, less than \code{B}
#' @param kappa0 positive integer, initial value of kappa
#' @param tol numeric, sets tolerance for convergence for \code{method = 'stochastic'}.
#'     Will give value of \code{0} if \code{'exact'} is used.
#' @param max.iter positive integer limiting number of iterations for
#'     \code{method = 'stochastic'}
#'
#' @return an object of class \code{simfast}.
#' @export
#'
#' @examples
simfast <- function(formula, data, intercept = FALSE, weights = NULL,
                    family = 'gaussian', returnmodel = TRUE, returndata = TRUE,
                    method = 'stochastic', multialpha = FALSE, B = 10000, k = 100,
                    kappa0 = 100, tol = 1e-10, max.iter = 20){
  if (missing(data)) {
    data <- environment(formula)
  } else if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  mf <- match.call()
  print(mf)
  m <- match(x = c("formula", "data"),
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.omit"
  mf[[1]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())
  mr <- mf
  mm <- attr(mf, which = "terms")
  if (intercept == FALSE){
    attr(mm, "intercept") <- 0
  }
  xm <- stats::model.matrix(object = mm)
  ym <- stats::model.response(mr)
  simfit <- simfast_m(x = xm, y = ym, weights = weights, family = family,
                      returndata = returndata, method = method, multialpha = multialpha,
                      B = B, k = k, kappa0 = kappa0, tol = tol, max.iter = max.iter)
  if (returnmodel == TRUE){
    simfit[['model']] <- mm
  }
  if (intercept == FALSE){
    simfit[['intercept']] <- FALSE
  }
  return(simfit)
}


