#' Find pMLE (Internal)
#'
#' @param alpha vector with magnitude 1, proposed index
#' @param y vector of reponse values
#' @param x matrix of predictor values
#' @param nn vector of positive integer weights
#' @param family \code{\link{family}} object passed by \code{simfast}
#'
#' @return list of numeric vectors
#'
#' @noRd
#' @keywords internal
#'
find.pMLE	<- function(alpha, y, x, nn, family){

	delta			<-	x %*% alpha
	tt				<-	order(delta)
	mle				<-	Iso::pava(y[tt],nn[tt])
	ll				<- 	find.loglik(mle, y[tt], nn[tt], family)

	return(list(ll=ll, ord=tt, mle=mle))
}


#' Finds Angle (Internal)
#'
#' @param z numeric vector of size 2
#'
#' @return numeric, angle
#'
#' @noRd
#' @keywords internal
#'
#'
my.inv	<-	function(z){

	if((z[1]>=0)*(z[2]>=0)) {theta <- asin(z[2])}
	if((z[1]<0)*(z[2]>=0)) {theta <- acos(z[1])}
	if((z[1]<0)*(z[2]<0)) {theta <- 2*pi-acos(z[1])}
	if((z[1]>=0)*(z[2]<0)) {theta <- 2*pi+asin(z[2])}

	return(theta)
}


#' Get Elements Greater Than kth Element (Internal)
#'
#' @param x numeric vector
#' @param k positive integer, location index
#'
#' @return numeric vector of index locations
#'
#' @noRd
#' @keywords internal
#'
getkmax		<-	function(x, k){

	y		<-	sort(x, decreasing=TRUE)
	kmax		<-	y[k]
	locs		<-	which(x >= kmax)

	return(locs)
}


#' Give First Row of Matrix, Return if Vector (Internal)
#'
#' @param x a matrix or vector
#'
#' @return a vector (first row of the matrix)
#'
#' @noRd
#' @keywords internal
#'
firstrow <- function(x){
  if (is.matrix(x)){
    return(x[1,])
  } else {
    return(x)
  }
}

