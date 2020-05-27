#' Predict Values from a \code{simfast} object with a matrix (Internal)
#'
#' @param object an object of class \code{simfast}
#' @param data optional numeric matrix of new data
#' @param rule integer that describes how to handle predictor values
#'     outside the range of the original data. By default, the value
#'     is \code{1}, which uses linear extrapolation to estimate outside
#'     values, but can take value \code{2} which provides the value of
#'     the closest edge point.
#' @param fn interpolation function passed from \code{predict.simfast}
#' @param interp extrapolation function passed from \code{predict.simfast}
#'
#' @return a numeric 2 by n matrix with the first row as predictions on the
#'     response scale and the second row a logical that is \code{TRUE} when
#'     the response is out of the bounds of the original data
#'
#' @noRd
#' @keywords internal
#'
mat_pred <- function(object, newdata, rule, fn, interp){
  newdata  <- as.matrix(newdata)
  newdata  <- stats::na.omit(newdata)
  rwnms <- rownames(newdata)
  leftlim  <- min(object$indexvals)
  rightlim <- max(object$indexvals)
  if (NCOL(newdata) != length(firstrow(object$alphahat))) {
    stop("'newdata' must have the same number of columns as object$x")
  }
  newivs   <- as.vector(newdata %*% firstrow(object$alphahat))
  oob      <- (newivs < leftlim | newivs > rightlim)
  if (sum(oob) > 0){
    if (rule == 1) {
      warning("Some predictors are outside bounds of original data, linear extrapolation
              will be used for these values")
    } else {
      warning("Some predictors are outside bounds of original data, the closest edge point
              will be used to predict these values")
    }
  }
  if (!is.numeric(newdata)) {
    stop("'newdata' must be a numeric matrix or object$model must contain a model.frame")
  } else {
    newyhat <- fn(newivs)
    if (rule == 2) {
      names(newyhat) <- rwnms
      return(rbind(newyhat, 'oob' = oob))
    } else {
      newyhat[oob] <- interp(newivs)[oob]
      names(newyhat) <- rwnms
      return(rbind(newyhat, 'oob' = oob))
    }
  }
}




#' Predict values of a \code{simfast} object with a data frame or a matrix
#'
#' Predict values of a \code{simfast} object with a data frame or a matrix.
#' Can also be called by running \code{\link{predict}} on \code{simfast} object.
#' Note that \code{simfast} fits models on the response level, so only
#' 'response' values are predicted.
#'
#' @param object an object of class \code{simfast}
#' @param newdata optional numeric matrix or data frame of new data. If
#'     none is provided, the original data will be used. If a data frame
#'     is provided, predictions cannot be made unless \code{returnmodel = TRUE}
#'     was selected when fitting your \code{simfast} object (this is the default
#'     value). New numeric matrix or data frame must match the structure of the
#'     original data (unused columns in data frames need not be included). Note
#'     that rows with \code{NA} values will be removed.
#' @param rule integer that describes how to handle predictor values
#'     outside the range of the original data. By default, the value
#'     is \code{1}, which uses linear extrapolation to estimate outside
#'     values, but can take value \code{2} which provides the value of
#'     the closest edge point.
#' @param oob logical parameter that, if \code{TRUE}, provides a second row
#'     of logical values that indicates if the data used to predict this
#'     value was out of the bounds of the original data
#' @param ... further arguments passed to \code{predict}
#'
#' @return a numeric vector of the specified prediction values on the response
#'     scale, unless \code{oob = TRUE}, then returns a matrix with the first row
#'     of prediction values on the response scale and a second row of logical
#'     values indicating if data used to predict was out of bounds.
#' @export
#'
#' @author
#'     Hanna Jankowski \email{hkj@@yorku.ca}
#'     Konstantinos Ntentes \email{kntentes@@yorku.ca} (maintainer)
#'
#'
#' @examples
#'
#' # See the example provided in the \code{\link{simfast}} documentation.
#'
predict.simfast <- function(object, newdata, rule = 1, oob = FALSE, ...){
  family <- object$family
  linkinv <- family$linkinv
  linkfun <- family$linkfun
  if (missing(newdata)){
    if (oob) {
      inr <- rep(TRUE, length(object$yhat))
      return(rbind(object$yhat, 'oob' = inr))
    }
    return(object$yhat)
  }
  oldy <- object$yhat
  predfun <- stats::approxfun(x = object$indexvals, y = oldy,
                              method = "linear", rule = 2, ties = mean)
  interp  <- function(ind) {
    x1 <- min(object$indexvals)
    x2 <- max(object$indexvals)
    y1 <- min(oldy)
    y2 <- max(oldy)
    m  <- (y2-y1)/(x2-x1)
    b  <- y2 - x2*m
    ylink <- m * ind + b
    if (family[[1]] %in% c('Gamma', 'poisson')) {
      ylink[ylink < 0] <- 0
    } else if (family[[1]] == 'binomial') {
      ylink[ylink < 0] <- 0
      ylink[ylink > 1] <- 1
    }
    return(ylink)
  }
  if (is.null(object$model)) {
    if (is.matrix(newdata)) {
      predvec <- mat_pred(object = object, newdata = newdata, rule = rule,
                          fn = predfun, interp = interp)
      if (oob) {
        return(predvec)
      }
      return(predvec[1, ])
    } else {
      newdata <- as.matrix(newdata)
      if (!is.numeric(newdata)){
        stop("If model = NULL in the simfast object, then a numeric matrix
           must be provided, not a data frame.")
      } else {
        predvec <- mat_pred(object = object, newdata = newdata, rule = rule,
                            fn = predfun, interp = interp)
        if (oob) {
          return(predvec)
        }
        return(predvec[1, ])
      }
    }
  } else {
    if (!is.data.frame(newdata)){
      newdata <- as.data.frame(newdata)
    }
    mf <- object$model
    mm <- stats::terms(mf)
    if (object$intercept == FALSE) {
      attr(mm, "intercept") <- 0
    }
    newmf <- stats::model.frame(formula = mm, data = newdata)
    newmm <- stats::model.matrix(object = mm, data = newmf)
    if (!is.null(object$offset)) {
      newoffset <- stats::model.offset(newmf)
      oldy <- object$yhat
      osy <- oldy
      osy <- linkfun(osy) - object$offset
      osy <- linkinv(osy)
      pfos <- stats::approxfun(x = object$indexvals, y = osy,
                               method = "linear", rule = 2, ties = mean)
      intos  <- function(ind) {
        x1 <- min(object$indexvals)
        x2 <- max(object$indexvals)
        y1 <- min(osy)
        y2 <- max(osy)
        m  <- (y2-y1)/(x2-x1)
        b  <- y2 - x2*m
        ylink <- m * ind + b
        if (family[[1]] %in% c('Gamma', 'poisson')) {
          ylink[ylink < 0] <- 0
        } else if (family[[1]] == 'binomial') {
          ylink[ylink < 0] <- 0
          ylink[ylink > 1] <- 1
        }
        return(ylink)
      }
      predvec <- mat_pred(object = object, newdata = newmm, rule = rule,
                          fn = pfos, interp = intos)
      yhatos <- linkfun(predvec[1,]) + newoffset
      predvec[1,] <- linkinv(yhatos)
      if (oob) {
        return(predvec)
      }
      return(predvec[1, ])
    } else {
      predvec <- mat_pred(object = object, newdata = newmm, rule = rule,
                          fn = predfun, interp = interp)
      if (oob) {
        return(predvec)
      }
      return(predvec[1, ])
    }
  }
}









