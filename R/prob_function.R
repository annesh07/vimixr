#' Cummulative summation functions for the probability matrix
#'
#' @param x the probability matrix
#'
#' @return cummulative summations
#' @export
#'
#' @examples
cum_clustprop_fullR <- function(x){
  cs <- Rfast::colCumSums(t(x))              # Cumulative sums for each row
  totals <- cs[nrow(cs),]                     # Total cumulative sum of each row of P1
  result <- totals - cs                # Element-wise matrix operation
  return(Rfast::rowsums(result))              # Row sums
}

cum_clustprop_var_fullR <- function(x){
  cs <- Rfast::colCumSums(t(x))              # Cumulative sums for each row
  totals <- cs[nrow(cs),]                     # Total cumulative sum of each row of P1
  result <- cs * (totals - cs)                # Element-wise matrix operation
  return(Rfast::rowsums(result))              # Row sums
}
