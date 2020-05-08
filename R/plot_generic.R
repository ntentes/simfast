#' Plot the estimated responses and index values of a \code{\link{simfast}} object
#'
#'
#' @param x an object of class \code{\link{simfast}}.
#' @param points boolean, if \code{TRUE}, plots points on the estimated ridge function.
#'     Plot will only contain a line representing the estimated ridge function if this
#'     option is changed to \code{FALSE}.
#' @param weights boolean, if \code{TRUE}, sizes the points according to their scaled weights.
#'     This only occurs if \code{pts = TRUE}, nothing happens otherwise.
#' @param predictor boolean, if \code{TRUE}, produces a \code{plotly} scatterplot
#'     of the estimated responses (\code{yhat}) and the observed predictors \code{x}
#'     (instead of \code{indexvals} on the predictor axes). This is an interactive 3-D
#'     scatterplot when \code{x} is 2 dimensional and a regular scatterplot if \code{x}
#'     is one dimensional. \code{plot.simfast} will throw an error with this option
#'     selected if \code{plotly} is not installed or \code{x} is larger than 2-dimensional.
#' @param ... all other arguments passed to \code{\link{plot}}.
#'
#' @export
#'
#' @examples
plot.simfast <- function(x, points = TRUE, weights = TRUE, predictor = FALSE, ...){
  xord <- order(x$indexvals)
  if (predictor) {
    if (NCOL(x$x) > 2) {
      stop("Predictor dimension >2, cannot plot predictor plot.")
    } else if(!requireNamespace("plotly", quietly = TRUE)) {
      stop("R Package 'plotly' required to plot interactive visualization.
           Please install and re-run.")
    } else if (NCOL(x$x) == 1) {
      fig <- plotly::plot_ly(x = x$x, y= x$yhat, type = 'scatter', mode = 'markers',
                     color = I('black'))
      fig <- plotly::layout(fig, title = 'Y-Hat vs. Observed Predictors',
                    xaxis = list(title = 'Observed X', zeroline = FALSE),
                    yaxis = list(title = 'Y-Hat', zeroline = FALSE))
      return(fig)
    } else {
      fig <- plotly::plot_ly(x = x$x[,1], y = x$x[,2], z = x$yhat,
                             type = 'scatter3d', mode = 'markers', color = I('black'))
      fig <- plotly::layout(fig, title = 'Y-Hat vs. Observed Predictors',
                            scene = list(xaxis = list(title = 'Observed X_1'),
                                    yaxis = list(title = 'Observed X_2'),
                                    zaxis = list(title = 'Y-Hat')))
      return(fig)
    }
  } else if (points) {
    if (weights) {
      wtsize <- x$weights
      wtsize <- scale(wtsize, center = FALSE)
      wtsize <- wtsize - min(wtsize) + 1
    } else {
      wtsize <- rep(1, length(x$yhat))
    }
    graphics::plot(x = x$indexvals[xord], y = x$yhat[xord], col = 'darkgrey', lwd = 2,
                   type = 'l', main = "Estimated Response Values vs. Single Index Values",
                   xlab = expression(paste("Single Index Values: ", x^T, hat(alpha))),
                   ylab = '', cex.main = 1, cex.lab = 0.9)
    graphics::title(ylab=expression(paste("Response Estimates: ", hat(Y))),
                    mgp=c(2.3,1,0), cex.lab=0.9)
    graphics::rug(x$indexvals)
    graphics::points(x = x$indexvals, y = x$yhat, pch = 20,
                     col = grDevices::rgb(0,0,0, alpha = 0.5), cex = wtsize)
  } else {
    graphics::plot(x = x$indexvals[xord], y = x$yhat[xord], col = 'darkgrey', lwd = 2,
                   type = 'l', main = "Estimated Response Values vs. Single Index Values",
                   xlab = expression(paste("Single Index Values: ", x^T, hat(alpha))),
                   ylab = '', cex.main = 1, cex.lab = 0.9)
    graphics::title(ylab=expression(paste("Response Estimates: ", hat(Y))),
                    mgp=c(2.3,1,0), cex.lab=0.9)
    graphics::rug(x$indexvals)
  }
}
