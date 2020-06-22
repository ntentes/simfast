#' Plot the estimated responses and index values of a \code{\link{simfast}} object
#'
#'
#' Plot the estimated responses and index values of a \code{\link{simfast}} object.
#' Can also be called by running \code{\link{plot}} on a \code{simfast} object.
#'
#'
#' @param x an object of class \code{\link{simfast}}.
#' @param points logical, if \code{TRUE}, plots points on the estimated ridge function.
#'     Plot will only contain a line representing the estimated ridge function if this
#'     option is changed to \code{FALSE}.
#' @param weights logical, if \code{TRUE}, sizes the points according to their scaled weights.
#'     This only occurs if \code{pts = TRUE}, nothing happens otherwise.
#' @param predictor logical, if \code{TRUE}, produces a \code{plotly} scatterplot
#'     of the estimated responses (\code{yhat}) and the observed predictors \code{x}
#'     (instead of \code{indexvals} on the predictor axes). This is an interactive 3-D
#'     scatterplot when \code{x} is 2 dimensional and a regular scatterplot if \code{x}
#'     is one dimensional. \code{plot.simfast} will throw an error with this option
#'     selected if \code{plotly} is not installed or \code{x} is larger than 2-dimensional.
#' @param offset logical, if \code{FALSE}, adjusts \code{yhat} values by the offset so
#'     that the isotonic relationship is visible. When \code{TRUE}, \code{yhat} values are
#'     left unscaled.
#' @param ... all other arguments passed to \code{\link{plot}}.
#'
#' @export
#'
#' @author
#'     Hanna Jankowski \email{hkj@@yorku.ca}
#'     Konstantinos Ntentes \email{kntentes@@yorku.ca} (maintainer)
#'
#' @examples
#'
#' # See the example provided in the \code{\link{simfast}} documentation.
#'
plot.simfast <- function(x, points = TRUE, weights = TRUE, predictor = FALSE, offset = TRUE, ...){
  indvals <- x$indexvals
  xord <- order(indvals)
  yhat <- x$yhat
  if (!offset) {
    if (!is.null(x$offset)){
      family <- x$family
      linkinv <- family$linkinv
      linkfun <- family$linkfun
      yhat <- linkfun(yhat) - x$offset
      yhat <- linkinv(yhat)
    }
  }
  if (predictor) {
    if (NCOL(x$x) > 2) {
      stop("Predictor dimension >2, cannot plot predictor plot.")
    } else if(!requireNamespace("plotly", quietly = TRUE)) {
      stop("R Package 'plotly' required to plot interactive visualization.
           Please install and re-run.")
    } else if (NCOL(x$x) == 1) {
      fig <- plotly::plot_ly(x = x$x, y= yhat, type = 'scatter', mode = 'markers',
                             color = I('black'), ...)
      if (is.null(dimnames(x$x)[[2]]) | length(dimnames(x$x)[[2]]) != 1){
        fig <- plotly::layout(fig, title = 'Y-Hat vs. Observed Predictors',
                              xaxis = list(title = 'Observed X', zeroline = FALSE),
                              yaxis = list(title = 'Y-Hat', zeroline = FALSE))
      } else {
        fig <- plotly::layout(fig, title = 'Y-Hat vs. Observed Predictors',
                              xaxis = list(title = dimnames(x$x)[[2]], zeroline = FALSE),
                              yaxis = list(title = 'Y-Hat', zeroline = FALSE))
      }
      return(fig)
    } else {
      fig <- plotly::plot_ly(x = x$x[,1], y = x$x[,2], z = yhat,
                             type = 'scatter3d', mode = 'markers',
                             color = I('black'), ...)
      if (is.null(dimnames(x$x)[[2]]) | length(dimnames(x$x)[[2]]) != 2){
        fig <- plotly::layout(fig, title = 'Y-Hat vs. Observed Predictors',
                              scene = list(xaxis = list(title = 'Observed X_1'),
                                           yaxis = list(title = 'Observed X_2'),
                                           zaxis = list(title = 'Y-Hat')))
      } else {
        fig <- plotly::layout(fig, title = 'Y-Hat vs. Observed Predictors',
                              scene = list(xaxis = list(title = dimnames(x$x)[[2]][1]),
                                           yaxis = list(title = dimnames(x$x)[[2]][2]),
                                           zaxis = list(title = 'Y-Hat')))
      }
      return(fig)
    }
  } else if (points) {
    if (weights) {
      wtsize <- x$weights
      wtsize <- scale(wtsize, center = FALSE)
      wtsize <- (wtsize - min(wtsize))*1.5 + 1
    } else {
      wtsize <- rep(1, length(yhat))
    }
    graphics::plot(x = indvals[xord], y = yhat[xord], col = 'darkgrey', lwd = 2,
                   type = 'l', main = "Estimated response values against single index values",
                   xlab = expression(paste("Single Index Values: ", hat(alpha)^T, x)),
                   ylab = '', cex.main = 1, cex.lab = 0.9, ...)
    graphics::title(ylab=expression(paste("Response Estimates: ", hat(Y))),
                    mgp=c(2.3,1,0), cex.lab=0.9)
    graphics::rug(indvals)
    graphics::points(x = indvals, y = yhat, pch = 20,
                     col = grDevices::rgb(0,0,0, alpha = 0.3), cex = wtsize)
  } else {
    graphics::plot(x = indvals[xord], y = yhat[xord], col = 'darkgrey', lwd = 2,
                   type = 'l', main = "Estimated response values against single index values",
                   xlab = expression(paste("Single Index Values: ", hat(alpha)^T, x)),
                   ylab = '', cex.main = 1, cex.lab = 0.9, ...)
    graphics::title(ylab=expression(paste("Response Estimates: ", hat(Y))),
                    mgp=c(2.3,1,0), cex.lab=0.9)
    graphics::rug(indvals)
  }
}
