#' Predict Values from a \code{simfast} object with a matrix (Internal)
#'
#' @param object an object of class \code{simfast}
#' @param data optional numeric matrix of new data
#' @param type character string specifying type of prediction, takes
#'     value \code{'link'} by default for predictions on the scale of
#'     the linear predictors, and takes \code{'response'} for predicted
#'     response values.
#' @param rule integer that describes how to handle predictor values
#'     outside the range of the original data. By default, the value
#'     is \code{1}, which uses linear extrapolation to estimate outside
#'     values, but can take value \code{2} which provides the value of
#'     the closest edge point.
#' @param fn interpolation function passed from \code{predict.simfast}
#' @param interp extrapolation function passed from \code{predict.simfast}
#' @param link link function passed from \code{predict.simfast}
#'
#' @return a numeric vector
#'
mat_pred <- function(object, data, type, rule, fn, interp, link){
  newdata  <- as.matrix(newdata)
  leftlim  <- min(object$indexvals)
  rightlim <- max(object$indexvals)
  oob      <- (object$indexvals < leftlim | object$indexvals > rightlim)
  if (NCOL(newdata) != length(firstrow(object$alphahat))) {
    stop("'newdata' must have the same number of columns as object$x")
  }
  newivs   <- newdata %*% firstrow(object$alphahat)
  if (!is.numeric(newdata)) {
    stop("'newdata' must be a numeric matrix or object$model must contain a model.frame")
  } else {
    newyhat <- fn(newivs)
    if (type == 'response') {
      if (rule == 2) {
        return(newyhat)
      } else {
        newyhat[oob] <- interp(newivs[oob])
        return(newyhat)
      }
    } else {
      if (rule == 2) {
        return(link(newyhat))
      } else {
        newyhat[oob] <- interp(newivs[oob])
        return(link(newyhat))
      }
    }
  }
}




#' Predict values of a \code{simfast} object with a data frame or a matrix
#'
#' @param object an object of class \code{simfast}
#' @param newdata optional numeric matrix or data frame of new data. If
#'     none is provided, the original data will be used. If a data frame
#'     is provided, predictions cannot be made unless \code{returnmodel = TRUE}
#'     was selected (this is the default value). New numeric matrix or data
#'     frame must match the structure of the original data.
#' @param type character string specifying type of prediction, takes
#'     value \code{'link'} by default for predictions on the scale of
#'     the linear predictors, and takes \code{'response'} for predicted
#'     response values (with selected link function applied).
#' @param rule integer that describes how to handle predictor values
#'     outside the range of the original data. By default, the value
#'     is \code{1}, which uses linear extrapolation to estimate outside
#'     values, but can take value \code{2} which provides the value of
#'     the closest edge point.
#' @param ... other parameters passed to \code{\link{predict}}
#'
#' @return a numeric vector of the specified prediction values
#' @export
#'
#' @examples
predict.simfast <- function(object, newdata, type = 'link', rule = 1, ...){
  linkfun <- object$link
  if (missing(newdata)){
    if (type == 'response'){
      return(object$yhat)
    } else {
      etahat <- linkfun(object$yhat)
      return(etahat)
    }
  }
  predfun <- stats::approxfun(x = object$indexvals, y = object$yhat,
                              method = "linear", rule = 2, ties = mean)
  interp  <- function(ind) {
    x1 <- min(object$indexvals)
    x2 <- max(object$indexvals)
    y1 <- min(object$yhat)
    y2 <- max(object$yhat)
    m  <- (y2-y1)/(x2-x1)
    b  <- y2 - x2*m
    yhatint <- m * ind + b
    return(yhatint)
  }
  if (is.null(object$model)) {
    if (is.matrix(newdata)){
      predvec <- mat_pred(object = object, data = newdata, type = type, rule = rule,
                        fn = predfun, interp = interp, link = linkfun)
      return(predvec)
    } else {
      newdata <- as.matrix(newdata)
      if (!is.numeric(newdata)){
        stop("If model = NULL in the simfast object, then a numeric matrix
           must be provided, not a data frame.")
      } else {
        predvec <- mat_pred(object = object, data = newdata, type = type, rule = rule,
                            fn = predfun, interp = interp, link = linkfun)
        return(predvec)
      }
    }
  } else {
    mf <- object$model
    mm <- attr(mf, which = "terms")
    if (object$intercept == FALSE) {
      attr(mm, "intercept") <- 0
    }
    newdata <- stats::model.matrix(object = mm)
    predvec <- mat_pred(object = object, data = newdata, type = type, rule = rule,
                        fn = predfun, interp = interp, link = linkfun)
    return(predvec)
  }
}








