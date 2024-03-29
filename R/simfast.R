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
#' @param offset numeric vector of model offsets , with length
#'     \code{n}. Takes default value \code{NULL} which uses no offset.
#' @param family a choice of the error distribution and link function to
#'     be used in the model. This can be a character string naming a family
#'     function, a family function or the result of a call to a family function.
#'     Currently supporting any of \code{\link{gaussian}, \link{binomial},
#'     \link{poisson},} and \code{\link{Gamma}}. The canonical link function is used
#'     by default, but all link functions available for these families are supported.
#' @param returndata logical value that when \code{TRUE} (the default value)
#'     returns the predictor matrix and response vector in the simfast object.
#' @param method when \code{x} has \code{d=2} columns, method can take \code{'exact'}
#'     argument, which uses an exact optimization method instead of a stochastic
#'     search. If \code{d} does not equal \code{2}, \code{simfast_m} will give a
#'     warning and automatically continue with a stochastic search (the default
#'     method, \code{method = 'stochastic'}).
#' @param multiout logical value, if \code{TRUE}, will return more than one
#'     \code{alpha} vector and \code{yhat} vector if available, separately from the main
#'     estimate (see Value section and Details).
#' @param B positive integer, sets number of index vectors to try when maximizing
#'     the likelihood
#' @param k positive integer, algorithmic parameter, more info coming, should be less
#'     than \code{B}
#' @param kappa0 positive integer, initial value of kappa, more info coming
#' @param tol numeric, sets tolerance for convergence for \code{method = 'stochastic'}.
#'     Will give value of \code{0} if \code{'exact'} is used.
#' @param max.iter positive integer limiting number of iterations for
#'     \code{method = 'stochastic'}
#'
#' @details
#'
#' For i=1,...,n, let X_i be the d-dimensional covariates and Y_i be the corresponding
#' one-dimensional response.   The isotonic single index model is written as
#'
#' g(mu) = f(a^T x),
#'
#' where x=(x_1,...,x_d)^T, g is a known link function, a is an d x 1 index vector,
#' and f is a nondecreasing function.  The algorithm finds the maximum likelihood
#' estimate of both f and a, assuming that f is an increasing function.  Implementaton
#' details can be found in ADD REFs, where theoretical justification of our estimator
#' (i.e. uniform consistency) is also given.   For the identifiability of isotonic
#' single index models, we refer to REFs.
#'
#' @return an object of class \code{simfast}, with the following structure:
#' \describe{
#'  \item{\code{x}}{if \code{returndata = TRUE}, this is the predictor matrix used to
#'       fit the model, otherwise it is \code{NULL}.}
#'  \item{\code{y}}{if \code{returndata = TRUE}, this is the response vector used to
#'       fit the model, otherwise it is \code{NULL}.}
#'  \item{\code{alphahat}}{\code{alpha} vector estimated by the model fit}
#'  \item{\code{yhat}}{vector of estimated response values}
#'  \item{\code{indexvals}}{vector of estimated single index values, the matrix product
#'      of \code{x} and \code{alphahat}}
#'  \item{\code{weights}}{vector of the integer weights used in the model fit}
#'  \item{\code{family}}{the \code{\link{family}} function provided to \code{simfast_m}}
#'  \item{\code{loglik}}{a numeric value of the log-likelihood at the estimate.}
#'  \item{\code{offset}}{returns \code{NULL} when calling \code{simfast_m}}
#'  \item{\code{tol}}{numeric convergence tolerance acheived during fitting with
#'      \code{method = 'stochastic'}. For \code{method = 'exact'}, this is \code{0}.}
#'  \item{\code{iter}}{number of iterations used to acheieve convergence. For
#'      \code{method = 'exact'}, this is \code{1}.}
#'  \item{\code{method}}{\code{method} used for fitting the model}
#'  \item{\code{model}}{returns \code{NULL} when calling \code{simfast_m}}
#'  \item{\code{intercept}}{returns \code{NULL} when calling \code{simfast_m}}
#'  \item{\code{multialphahat}}{returns all estimated \code{alphahat} vectors
#'      if \code{multiout = TRUE} as a matrix if there is more than one, and
#'      as a vector if there is only one.}
#'  \item{\code{multiyhat}}{returns all estimated \code{yhat} vectors
#'      if \code{multiout = TRUE} as a matrix if there is more than one, and
#'      as a vector if there is only one.}
#' }
#'
#'
#' @seealso \code{\link{simfast}} for formula support (and support for offsets).
#'
#' @export
#' @md
#'
#' @author
#'     Hanna Jankowski: hkj@@yorku.ca
#'     Konstantinos Ntentes: tino.ntentes@@gmail.com (maintainer)
#'
#' @examples
#'
#' # Generate predictor data uniformly on [-4, 4]^2
#'
#' \dontshow{set.seed(100)}
#' d <- 2
#' preds <- matrix(stats::runif(d*500, min = -4, max = 4), 500, d)
#'
#' # Choose true alpha in R^2 with magnitude 1
#'
#' alpha <- c(-2, 0.5)/sqrt(4.25)
#'
#' # Set true ridge function (piecewise)
#'
#' fn <- function(x){
#'   #  RANGE              #  VALUES
#'    (x< -3)             * (x+5)               +
#'    (x< -1)*(x>= -3)    * 2                   +
#'    (x< 1.525)*(x>= -1) * (0.1*(x + 1)^3 + 2) +
#'    (x>=1.525)          * (0.4*x + 3)
#' }
#'
#' # Generate Gaussian response values with the same number of rows
#'
#' indexvals <- preds %*% alpha
#' mu  <- fn(indexvals)
#' y <- stats::rnorm(500, mean = mu)
#'
#' # Set weights
#'
#' weights <- rep(c(5, 2, 3, 5, 2, 2, 1, 5, 3, 1), 50)
#'
#' # Fit simfast object
#'
#' sfobj <- simfast_m(x = preds, y = y, weights = weights)
#'
#' # Predict response values and inspect a subset
#' # Note: simfast only predicts 'response' values
#'
#' predlink <- predict(sfobj)
#' predlink[1:20]
#'
#' # predict response values with new data and inspect
#'
#' newpreds <- matrix(stats::runif(d*200, min = -3, max = 3), 200, d)
#' predresp <- predict(sfobj, newdata = newpreds,
#'                     type = 'response', rule = 1)
#' predresp[1:20]
#'
#' # Plot simfast object Y-Hat vs. indexvals showcasing weights
#' # and compare to the true ridge function fn()
#' \dontshow{devAskNewPage(ask = FALSE)}
#' plot(sfobj, xlim = c(-5, 5))
#' \dontshow{invisible(readline(prompt="Press <Return> to plot the truth."))}
#' curve(fn, lwd = 4, lty = 2, col = rgb(1, 0, 0, 0.6), add=TRUE)
#'
#' # 3-D Interactive plot of observed predictors vs. Y-Hat
#' # Requires library(plotly)
#' \dontshow{invisible(readline(prompt="Press <Return> to see 3-D Plot if library(plotly) is installed."))}
#' if (requireNamespace("plotly", quietly = TRUE)){
#'   plot(sfobj, predictor = TRUE)
#' }
#'
#'
#'
simfast_m <- function(x, y, weights = NULL, offset = NULL, family = 'gaussian',
                      returndata = TRUE, method = 'stochastic', multiout = FALSE,
                      B = 10000, k = 100, kappa0 = 100, tol = 1e-10, max.iter = 20) {
  if (NCOL(y) == 1) {
    if (is.character(y)) y <- factor(y)
    if (is.factor(y)) y <- y != levels(y)[1L]
  } else {
    stop("Response 'y' must be a single column.")
  }
  if (NROW(x) != NROW(y)) {
    stop("predictor matrix 'x' and response 'y' do not have the same number of rows -- fitting cannot continue.")
  }
  if ((!is.null(weights)) & (length(weights) != NROW(x))){
    stop('Weights vector is the wrong length -- fitting cannot continue.')
  }
  if (!is.null(weights) & !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }
  if ((!is.null(offset)) && (length(offset) != NROW(x))){
    stop('Offset vector is the wrong length -- fitting cannot continue')
  }
  if (sum(!complete.cases(x)) > 0) {
    warning("Missing values in predictor matrix 'x' -- rows with missing data will be removed")
    cc <- complete.cases(x)
    x <- x[cc, ]
    y <- y[cc]
    if (!is.null(weights)) {
      weights <- weights[cc]
    }
    if (!is.null(offset)) {
      offset <- offset[cc]
    }
  }
  if (sum(!complete.cases(y)) > 0) {
    warning("Missing values in response vector 'y' -- rows with missing data will be removed")
    cc <- complete.cases(y)
    x <- x[cc, ]
    y <- y[cc]
    if (!is.null(weights)) {
      weights <- weights[cc]
    }
    if (!is.null(offset)) {
      offset <- offset[cc]
    }
  }
  if (!is.null(weights) & (sum(!complete.cases(weights)) > 0)) {
    warning("Missing values in weights vector -- rows with missing weights will be removed")
    cc <- complete.cases(weights)
    x <- x[cc, ]
    y <- y[cc]
    if (!is.null(weights)) {
      weights <- weights[cc]
    }
    if (!is.null(offset)) {
      offset <- offset[cc]
    }
  }
  if (!is.null(offset)) {
    if (sum(!complete.cases(offset)) > 0) {
      warning("Missing values in offset vector -- rows with missing offset values will be removed")
      cc <- complete.cases(offset)
      x <- x[cc, ]
      y <- y[cc]
      if (!is.null(weights)) {
        weights <- weights[cc]
      }
      if (!is.null(offset)) {
        offset <- offset[cc]
      }
    }
  }
  if (!is.null(weights) & any(weights < 0, na.rm = T))
    stop("negative weights not allowed")
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
  linkfun <- family$linkfun
  if (!is.function(linkfun))
    stop("Cannot find specified link function")
  x <- as.matrix(x)
  if (!is.numeric(x)){
    stop('Predictor matrix must be numeric.')
  }
  linkinv <- family$linkinv
  if (!is.null(offset)) {
    if (all(!is.finite(offset))) {
      stop("Offset value is not finite. Please check offset vector.")
    }
    oldweights <- weights
    weights <- weights * linkinv(offset)
    comment(weights) <- 'offset'
    oldy <- y
    y <- linkfun(y) - offset
    y <- linkinv(y)
  }
  if (famname == 'binomial') {
    newy <- y * weights
    if (is.null(comment(weights))) {
      if (!isTRUE(all.equal(newy, round(newy)))) {
        warning("y*weights are not all within integer tolerance")
      }
    }
  }
  if (famname == 'poisson') {
    if (is.null(comment(weights))) {
      if (!isTRUE(all.equal(y, round(y)))) {
        warning("count response values are not all within integer tolerance")
      }
    }
  }
  if (famname == 'Gamma') {
    if (!isTRUE(all.equal(y, abs(y)))) {
      stop("Response values are not all positive, regression cannot continue.")
    }
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
  alphahat <- firstrow(fit$alphahat)
  loglik <- fit$maxll

  obj <- list('alphahat' = alphahat, 'yhat' = yhat, 'indexvals' = indexvals,
              'weights' = weights, 'family' = family, 'loglik' = loglik,
              'offset' = NULL, 'tol' = tol, 'iter' = iter, 'method' = method,
              'model' = NULL, 'intercept' = NULL)

  if (!is.null(offset)) {
    y <- oldy
    yhat <- obj[['yhat']]
    yhat <- linkfun(yhat) + offset
    obj[['yhat']] <- linkinv(yhat)
    obj[['offset']] <- offset
    obj[['weights']] <- oldweights
  }

  if (multiout == TRUE) {
    obj <- append(obj, list("multialphahat" = fit$alphahat, "multiyhat" = fit$mle))
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
#'     included, so change argument \code{intercept = TRUE} to include one. When
#'     including categorical predictors, be sure to set
#'     \code{options('contrasts')} in your global options to a desired setting.
#'     For \code{'binomial'} response, the vector can be binary values or a vector of
#'     proportions, but should include proper weight vector (a vector of the denominators
#'     of the proportions) is provided. Categorical vectors (character strings or
#'     factors) will automatically be translated into a logical vector with
#'     the baseline factor level a 'success' (takes value \code{1}). Poisson responses
#'     can be integer counts or rates, but should include a proper weight vector (a
#'     vector of the denominators of the rates). Offsets can also be specified in the
#'     formula. Note that multiple offsets are combined, and that duplicate
#'     offsets are only counted once.
#' @param data optional data frame (or object coercible to a data frame by
#'     \code{\link{as.data.frame}}) containing the variables in the model. Variables
#'     are taken from \code{environment(formula)} if not found in \code{data}.
#' @param intercept logical value, if \code{FALSE} (the default value), then the
#'     model given by the formula does not include an intercept value (even when
#'     including a 1, for example: \code{z ~ 1 + x + y} will only include columns for
#'     \code{x} and \code{y}).
#' @param weights optional vector of positive integer weights, with length
#'     \code{n}. Takes default value \code{NULL} which uses equal weights.
#' @param offset numeric vector of model offsets , with length
#'     \code{n}. Takes default value \code{NULL} which uses no offset.
#'     If an offset is provided here and in the formula, they are combined.
#' @param family a choice of the error distribution and link function to
#'     be used in the model. This can be a character string naming a family
#'     function, a family function or the result of a call to a family function.
#'     Currently supporting any of \code{\link{gaussian}, \link{binomial},
#'     \link{poisson},} and \code{\link{Gamma}}. The canonical link function is used
#'     by default, but all link functions available for these families are supported.
#' @param returnmodel logical value that when \code{TRUE} (the default value)
#'     attaches the \code{\link{model.frame}} object to the simfast object. Leave
#'     as \code{TRUE} to properly use the \code{\link{predict}} function.
#' @param returndata logical value that when \code{TRUE} (the default value)
#'     returns the predictor matrix and response vector in the simfast object.
#' @param method when \code{x} has \code{d=2} columns, method can take \code{'exact'}
#'     argument, which uses an exact optimization method instead of a stochastic
#'     search. If \code{d} does not equal \code{2}, \code{simfast} will give a
#'     warning and automatically continue with a stochastic search (the default
#'     method, \code{method = 'stochastic'}).
#' @param multiout logical value, if \code{TRUE}, will return more than one
#'     \code{alpha} vector and \code{yhat} vector if available, separately from the main
#'     estimate (see Value section and Details).
#' @param B positive integer, sets number of index vectors to try when maximizing
#'     the likelihood
#' @param k positive integer, algorithmic parameter, more info coming, should be less
#'     than \code{B}
#' @param kappa0 positive integer, initial value of kappa, more info coming
#' @param tol numeric, sets tolerance for convergence for \code{method = 'stochastic'}.
#'     Will give value of \code{0} if \code{'exact'} is used.
#' @param max.iter positive integer limiting number of iterations for
#'     \code{method = 'stochastic'}
#'
#' @details
#'
#'  For i=1,...,n, let X_i be the d-dimensional covariates and Y_i be the corresponding
#' one-dimensional response.   The isotonic single index model is written as
#'
#' g(mu) = f(a^T x),
#'
#' where x=(x_1,...,x_d)^T, g is a known link function, a is an d x 1 index vector,
#' and f is a nondecreasing function.  The algorithm finds the maximum likelihood
#' estimate of both f and a, assuming that f is an increasing function.  Implementaton
#' details can be found in ADD REFs, where theoretical justification of our estimator
#' (i.e. uniform consistency) is also given.   For the identifiability of isotonic
#' single index models, we refer to REFs.
#'
#' @return an object of class \code{simfast}, with the following structure:
#' \describe{
#'  \item{\code{x}}{if \code{returndata = TRUE}, this is the model matrix used to
#'       fit the model, otherwise it is \code{NULL}.}
#'  \item{\code{y}}{if \code{returndata = TRUE}, this is the response vector used to
#'       fit the model, otherwise it is \code{NULL}.}
#'  \item{\code{alphahat}}{\code{alpha} value estimated by the model fit}
#'  \item{\code{yhat}}{vector of estimated response values}
#'  \item{\code{indexvals}}{vector of estimated single index values, the matrix product
#'      of \code{x} and \code{alphahat}}
#'  \item{\code{weights}}{vector of the integer weights used in the model fit}
#'  \item{\code{family}}{the \code{\link{family}} function provided to \code{simfast_m}}
#'  \item{\code{loglik}}{a numeric value of the log-likelihood at the estimate.}
#'  \item{\code{offset}}{a numeric vector specifying the offset provided in the
#'      model formula.}
#'  \item{\code{tol}}{numeric convergence tolerance acheived during fitting with
#'      \code{method = 'stochastic'}. For \code{method = 'exact'}, this is \code{0}.}
#'  \item{\code{iter}}{number of iterations used to acheieve convergence. For
#'      \code{method = 'exact'}, this is \code{1}.}
#'  \item{\code{method}}{\code{method} used for fitting the model}
#'  \item{\code{model}}{the \code{\link{model.frame}} generated by the formula object
#'      which is used to generate the \code{\link{model.matrix}} and
#'      \code{\link{model.response}} to pass to \code{simfast_m}}
#'  \item{\code{intercept}}{the \code{intercept} rule selected in the argument}
#'  \item{\code{multialphahat}}{returns all estimated \code{alphahat} vectors
#'      if \code{multiout = TRUE} as a matrix if there is more than one, and
#'      as a vector if there is only one.}
#'  \item{\code{multiyhat}}{returns all estimated \code{yhat} vectors
#'      if \code{multiout = TRUE} as a matrix if there is more than one, and
#'      as a vector if there is only one.}
#' }
#'
#'
#' @seealso \code{\link{simfast_m}} for providing model matrices instead of a formula, as well
#'      as more examples.
#'
#' @export
#' @md
#'
#' @author
#'     Hanna Jankowski: hkj@@yorku.ca
#'     Konstantinos Ntentes: tino.ntentes@@gmail.com (maintainer)
#'
#' @examples
#'
#' ## Load esophageal cancer dataset
#' esoph <- datasets::esoph
#' str(esoph) # note that three variables are ordered factors
#' esoph$ntotal <- esoph$ncases + esoph$ncontrols #use as offset
#'
#' ## subset the data frame for training
#' set.seed(1) # keep from getting data OOB warning in predict()
#' nobs <- NROW(esoph)
#' ind <- sample(1:nobs, size = round(nobs * 0.8))
#' esophtrain <- esoph[ind, ]
#' esophtest  <- esoph[-ind, ]
#'
#' ## fit a model with formulas, including ordered/regular factors
#' ## and support for offsets. similar syntax to glm()
#' sfobj <- simfast(ncases ~ offset(log(ntotal)) + tobgp + alcgp + agegp,
#'                  data = esophtrain, family = poisson(link = 'log'))
#'
#' glmobj <- glm(ncases ~ offset(log(ntotal)) + tobgp + alcgp + agegp,
#'               data = esophtrain, family = poisson(link = 'log'))
#'
#' ## Plot the relationship of estimated responses vs. index values
#' # Not isotonic because of offset
#' plot(sfobj)
#' # Y-hats adjusted to same scale
#' \dontshow{invisible(readline(prompt="Press <Return> to plot the model with rate response values."))}
#' \dontshow{devAskNewPage(ask = FALSE)}
#' plot(sfobj, offset = FALSE)
#'
#' ## Predictions from simfast and glm rounded to nearest integer
#' sfpred <- round(predict(sfobj, newdata = esophtest))
#' # Note that simfast only predicts 'response' values
#' sfpred
#' glmpred <- round(predict(glmobj, newdata = esophtest, type = 'response'))
#' glmpred
#'
#' ## Compare squared residuals
#' sum((sfpred - esophtest$ncases)^2)   #simfast prediction
#' sum((glmpred - esophtest$ncases)^2)  #glm prediction
#'
#'
simfast <- function(formula, data, intercept = FALSE, weights = NULL, offset = NULL,
                    family = 'gaussian', returnmodel = TRUE, returndata = TRUE,
                    method = 'stochastic', multiout = FALSE, B = 10000, k = 100,
                    kappa0 = 100, tol = 1e-10, max.iter = 20) {
  if (missing(data)) {
    data <- environment(formula)
  } else if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  mf <- match.call()
  m <- match(x = c("formula", "data", "weights", "offset"),
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.omit"
  mf[[1]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())
  mm <- attr(mf, which = "terms")
  if (intercept == FALSE){
    attr(mm, "intercept") <- 0
  }
  xm <- stats::model.matrix(object = mm, data = mf)
  ym <- stats::model.response(mf)
  os <- stats::model.offset(mf)
  if (!is.null(os) && (length(os) != NROW(xm))){
    stop('Invalid offset vector provided. Must be same length as number of observations')
  }
  mw <- stats::model.weights(mf)
  if (!is.null(mw) && (length(mw) != NROW(xm))){
    stop('Invalid weights vector provided. Must be same length as number of observations')
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
  if (!is.null(data) && (NROW(data) != NROW(xm))) {
    rows_dropped <- NROW(data) - NROW(xm)
    warning(paste0(rows_dropped, " rows dropped before fitting -- check for NAs/missing values in data or weights/offset arguments"))
  }
  simfit <- simfast_m(x = xm, y = ym, weights = mw, family = family, offset = os,
                      returndata = returndata, method = method, multiout = multiout,
                      B = B, k = k, kappa0 = kappa0, tol = tol, max.iter = max.iter)
  if (returnmodel == TRUE){
    simfit[['model']] <- mm
  }
  simfit[['intercept']] <- intercept
  return(simfit)
}


