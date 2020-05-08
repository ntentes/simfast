#' Fitting isotonic generalized single-index regression models via maximum likelihood
#'
#'
#' Fitting isotonic generalized single-index regression models via maximum likelihood
#' with support for estimating response values with \code{\link{predict}} and plotting
#' values with \code{\link{plot}}. Also includes support for built-in regression families
#' (see Arguments), and formula support through the wrapper \code{\link{simfast}}.
#'
#' @param x numeric \code{matrix} of observed predictor values, with \code{n}
#'     rows and \code{d} columns
#' @param y numeric \code{vector} of observed response values, of length \code{n}.
#'     For \code{'binomial'} response, vector can be binary values or a vector of
#'     proportions, as long as a proper weight vector (a vector of the denominators
#'     of the proportions, see details) is included. Categorical vectors (character
#'     strings or factors) will automatically be translated into logical vector with
#'     the baseline factor level a 'success' (takes value \code{1}).
#' @param weights optional vector of positive integer weights, with length
#'     \code{n}. Takes default value \code{NULL} which uses equal weights.
#' @param family a choice of the error distribution and link function to
#'     be used in the model. This can be a character string naming a family
#'     function, a family function or the result of a call to a family function.
#'     Currently supporting any of \code{\link{gaussian}, \link{binomial},
#'     \link{poisson},} and \code{\link{Gamma}}. The canonical link function is used
#'     by default, but all link functions available for these families are supported.
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
#'
#' @return an object of class \code{simfast}, with the following structure:
#' \describe{
#'  \item{\code{x}}{if \code{returndata = TRUE}, this is the predictor matrix used to
#'       fit the model, otherwise it is \code{NULL}.}
#'  \item{\code{y}}{if \code{returndata = TRUE}, this is the response vector used to
#'       fit the model, otherwise it is \code{NULL}.}
#'  \item{\code{alphahat}}{\code{alpha} value estimated by the model fit, will return a
#'       matrix of multiple vectors if \code{multialpha = TRUE}}
#'  \item{\code{yhat}}{vector of estimated response values}
#'  \item{\code{indexvals}}{vector of estimated single index values, the matrix product
#'      of \code{x} and \code{alphahat}}
#'  \item{\code{weights}}{vector of the integer weights used in the model fit}
#'  \item{\code{family}}{the \code{\link{family}} function provided to \code{simfast_m}}
#'  \item{\code{link}}{the link function proposed by \code{family}}
#'  \item{\code{tol}}{numeric convergence tolerance acheived during fitting with
#'      \code{method = 'stochastic'}. For \code{method = 'exact'}, this is \code{0}.}
#'  \item{\code{iter}}{number of iterations used to acheieve convergence. For
#'      \code{method = 'exact'}, this is \code{1}.}
#'  \item{\code{method}}{\code{method} used for fitting the model}
#'  \item{\code{model}}{returns \code{NULL} when calling \code{simfast_m}}
#'  \item{\code{intercept}}{returns \code{NULL} when calling \code{simfast_m}}
#' }
#'
#'
#' @seealso \code{\link{simfast}} for formula support.
#'
#' @export
#' @md
#'
#' @author
#'     Hanna Jankowski <hkj@@yorku.ca>
#'     Konstantinos Ntentes <kntentes@@yorku.ca> (maintainer)
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
  famname <- family[[1]]
  if (!famname %in% c('gaussian', 'binomial', 'poisson', 'Gamma')) {
    stop('Simfast does not support specified family')
  }
  if (is.null(family$family)) {
    print(family)
    stop("Specified family cannot be found.")
  }
  if (famname == 'binomial') {
    newy <- y * weights
    if (!isTRUE(all.equal(newy, round(newy)))) {
      warning("y*weights are not all within integer tolerance")
    }
  }
  if (famname == 'poisson') {
    if (!isTRUE(all.equal(y, round(y)))) {
      warning("count response values are not all within integer tolerance")
    }
  }
  if (famname == 'Gamma') {
    if (!isTRUE(all.equal(y, abs(y)))) {
      stop("Response values are not all positive, regression cannot continue.")
    }
  }
  linkfun <- family$linkfun
  if (!is.function(linkfun))
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
              'family' = family, 'link' = linkfun,'tol' = tol, 'iter' = iter,
              'method' = method, 'model' = NULL, 'intercept' = NULL)

  if (multialpha == TRUE) {
    obj <- append(list("alphahat" = fit$alphahat), obj)
  } else {
    obj <- append(list("alphahat" = firstrow(fit$alphahat)), obj)
  }

  if (returndata == TRUE) {
    obj <- append(list('x' = x, 'y' = y), obj)
  } else {
    obj <- append(list('x' = NULL, 'y' = NULL), obj)
  }

  class(obj) <- 'simfast'
  return(obj)
}


#' Fitting isotonic generalized single-index regression models via maximum likelihood
#' with formula support
#'
#' Fitting isotonic generalized single-index regression models via maximum likelihood
#' with support for estimating response values with \code{\link{predict}} and plotting
#' values with \code{\link{plot}}. Also includes support for formula objects, data frames,
#' and built-in regression families (see Arguments).
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
#' @param data optional data frame (or object coercible to a data frame by
#'     \code{\link{as.data.frame}}) containing the variables in the model. Variables
#'     are taken from \code{environment(formula)} if not found in \code{data}.
#' @param intercept optional boolean, if \code{FALSE} (the default value), then the
#'     model given by the formula does not include an intercept value (even when
#'     including a 1, for example: \code{z ~ 1 + x + y} will only include columns for
#'     \code{x} and \code{y}).
#' @param weights optional vector of positive integer weights, with length
#'     \code{n}. Takes default value \code{NULL} which uses equal weights.
#' @param family a choice of the error distribution and link function to
#'     be used in the model. This can be a character string naming a family
#'     function, a family function or the result of a call to a family function.
#'     Currently supporting any of \code{\link{gaussian}, \link{binomial},
#'     \link{poisson},} and \code{\link{Gamma}}. The canonical link function is used
#'     by default, but all link functions available for these families are supported.
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
#' @return an object of class \code{simfast}, with the following structure:
#' \describe{
#'  \item{\code{x}}{if \code{returndata = TRUE}, this is the model matrix used to
#'       fit the model, otherwise it is \code{NULL}.}
#'  \item{\code{y}}{if \code{returndata = TRUE}, this is the response vector used to
#'       fit the model, otherwise it is \code{NULL}.}
#'  \item{\code{alphahat}}{\code{alpha} value estimated by the model fit, will return a
#'       matrix of multiple vectors if \code{multialpha = TRUE}}
#'  \item{\code{yhat}}{vector of estimated response values}
#'  \item{\code{indexvals}}{vector of estimated single index values, the matrix product
#'      of \code{x} and \code{alphahat}}
#'  \item{\code{weights}}{vector of the integer weights used in the model fit}
#'  \item{\code{family}}{the \code{\link{family}} function provided to \code{simfast_m}}
#'  \item{\code{link}}{the link function proposed by \code{family}}
#'  \item{\code{tol}}{numeric convergence tolerance acheived during fitting with
#'      \code{method = 'stochastic'}. For \code{method = 'exact'}, this is \code{0}.}
#'  \item{\code{iter}}{number of iterations used to acheieve convergence. For
#'      \code{method = 'exact'}, this is \code{1}.}
#'  \item{\code{method}}{\code{method} used for fitting the model}
#'  \item{\code{model}}{the \code{\link{model.frame}} generated by the formula object
#'      which is used to generate the \code{\link{model.matrix}} and
#'      \code{\link{model.response}} to pass to \code{simfast_m}}
#'  \item{\code{intercept}}{the \code{intercept} rule selected in the argument}
#' }
#'
#'
#' @seealso \code{\link{simfast_m}} for providing model matrices instead of a formula.
#'
#' @export
#' @md
#'
#' @author
#'     Hanna Jankowski \email{hkj@@yorku.ca}
#'     Konstantinos Ntentes \email{kntentes@@yorku.ca} (maintainer)
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
  simfit[['intercept']] <- intercept
  return(simfit)
}


