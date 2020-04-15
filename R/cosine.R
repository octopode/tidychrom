#' Cosine Similarity of Two Vectors
#' Lifted verbatim from Fridolin Wild's LSA package (https://cran.r-project.org/web/packages/lsa/)
#' Serves here as a helper to cosine_spectra().
#'
#' Takes two vector-like objects
#' @param x vector 1
#' @param y vector 2
#' @return cosine similarity score of the two vector (dbl)
#' @keywords cosine similarity vector
#' @export
#' @examples
#' cos <- cosine_spectra(ions %>% filter(scan == 114), ions %>% filter(scan == 114))
#' > 0.94
#'

### cosine.R
###
### 2009-09-14:
###   * added vector vs. matrix comparison
###     to reduce data load when looking for associations
### 2005-11-21:
###   * added lazy calculation:
###     calc only below diagonale; diag = 1; add t(co)
###   * sqrt() over crossprod(x) and crossprod(y)
### 2005-11-09:
###   * bugfix cosvecs
###   * integrated cosvecs into cosine by doing type dependant processing
### 2005-08-26:
###   * rewrote cosvecs function to crossprod
###

cosine <- function( x, y=NULL ) {

    if ( is.matrix(x) && is.null(y) ) {

        co = array(0,c(ncol(x),ncol(x)))
        f = colnames( x )
        dimnames(co) = list(f,f)

        for (i in 2:ncol(x)) {
            for (j in 1:(i-1)) {
                co[i,j] = cosine(x[,i], x[,j])
            }
        }
        co = co + t(co)
        diag(co) = 1

        return (as.matrix(co))

    } else if ( is.vector(x) && is.vector(y) ) {
        return ( crossprod(x,y) / sqrt( crossprod(x)*crossprod(y) ) )
    } else if ( is.vector(x) && is.matrix(y) ) {

		 co = vector(mode='numeric', length=ncol(y))
		 names(co) = colnames(y)
		 for (i in 1:ncol(y)) {
			co[i] = cosine(x,y[,i])
		 }
		 return(co)

	 }	else {
        stop("argument mismatch. Either one matrix or two vectors needed as input.")
    }

}
