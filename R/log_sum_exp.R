#' Log-sum-exponential computation on the log probability allocation matrix
#'
#' @param Plog log probability allocation matrix
#'
#' @returns per sample ordered log probability allocation matrix
#' 
#' @importFrom Rfast rowMaxs rowsums
#' 
#' @export
log_sum_exp <- function(Plog){
  mx   <- Rfast::rowMaxs(Plog, value = TRUE)   
  Plog <- Plog - mx                              
  Plog <- Plog - log(Rfast::rowsums(exp(Plog))) 
  Plog                                           
}