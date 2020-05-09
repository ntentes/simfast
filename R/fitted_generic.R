#' Extract fitted values of a \code{simfast} object
#'
#' @param object on object of class \code{simfast}
#' @param ... other parameters passed to \code{\link{fitted}}
#'
#' @return a numeric vector of fitted values
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
fitted.simfast <- function (object, ...){
  predict.simfast(object, type = 'response')
}
