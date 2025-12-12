#' Root for a0 hyper-parameter for Sparse DPMM
#'
#' @param logP log of probability allocation matrix
#' @param X observed data
#' @param a_min minimum value of a0 for grid search
#' @param a_max maximum value of a0 for grid search
#' @param grid_points number of points for grid search
#'
#' @importFrom stats uniroot
#' 
#' @export
eBa0 <- function(logP, X, 
                        a_min = min(1e-8, 1/ncol(X)),
                        a_max = max(1e+6, ncol(X)),
                        grid_points = min(ncol(X), 10000)) {
  
  N <- nrow(logP)
  T0 <- ncol(logP)
  D <- ncol(X)
  
  P <- exp(logP)
  Rk  <- Rfast::colsums(P)
  Sk <- 0.5*mat_mult(t(P),X^2)
  PSk <- 0.5*Sk^2 - 0.5*Sk*Rk
  
  df_vec <- function(a_vec) {
    n0 <- length(a_vec)
    f0_mat <- outer(a_vec, Rk, "+")
    f0 <- 0.5* D * Rfast::rowsums(sweep((2*f0_mat + 1) / (f0_mat^2), 2, Rk,"*"))
    
    f1_list <- lapply(1:n0, function(i){(PSk + Rk*a_vec[i])/((a_vec[i]+Sk)^2)})
    output <- 0.5*(f0 - sapply(f1_list, sum))
    return(output)
  }
  vll <- function(r_vec){
    n0 <- length(r_vec)
    f_list <- lapply(1:n0, function(i){D*Rk*digamma(r_vec[i] + Rk) -
        Rk*log(r_vec[i] + 0.5*Sk) - (Sk*(r_vec[i] + Rk))/(r_vec[i]+0.5*Sk) -
        D*Rk/(N+1+Rk)})
    val0 <- -0.5*D*N*T0*log(2*pi) + 0.5*sapply(f_list, sum)
    return(val0)
  }
  
  #Interval for root
  a_grid <- 10^(seq(log10(a_min), log10(a_max), length.out = grid_points))
  dfvals_vec <- df_vec(a_grid) 
  
  #vectorized root search
  root_search <- Vectorize(function(i) {
    stats::uniroot(df_vec, interval = c(a_grid[i], a_grid[i + 1]))$root
  }) 
  
  #Detect sign change in df
  idx <- which(dfvals_vec[-length(dfvals_vec)] * dfvals_vec[-1] < 0)
  
  if (length(idx) == 0) {
    diffs <- abs(diff(dfvals_vec))
    p <- floor(log10(max(diffs))) + 1  
    
    threshold <- min(0.01, 1/max(p, 1)) 
    
    d0 <- ((diffs / max(diffs)) < threshold)
    
    idx0 <- length(a_grid) - which(cumprod(rev(d0))==0)[1]
    if (idx0 == 1){
      root <- N #default a0
    } else {root <- a_grid[idx0]} 
  } else {
    roots <- root_search(idx)
    root <- roots[which.max(vll(roots))]
  }
  
  return(root)
}

