#' S3 plotting function for `CVIoutput`objects'
#'
#' @param x a CVIoutput object
#' @param ... additional arguments
#'
#' @import patchwork
#'
#' @export

plot.CVIoutput <- function(x, ...){

  x$ELBO_viz + x$PCA_viz

}
