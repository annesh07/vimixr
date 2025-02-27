elbo_fixed_diagonal <- function(X, inverts, params){

  L1 <- params[["post_mean_eta"]]
  L2 <- params[["post_precision_scalar_eta"]]
  L20 <- params[["prior_precision_scalar_eta"]]
  inv_C0 <- inverts[["inv_C0"]]
  inv_C00 <- inverts[["inv_C00"]]
  P <- params["P"]
  D <- ncol(X) # dimension of the data


  #the eta's
  L21 <- sweep(L1, 1, L2, "/")
  e20 <- diag(-0.5 * post_mean_eta %*% tcrossprod(inv_C00, L1))
  e21 <- -0.5 * (D * L20)/L2
  e22 <- Mu0 %*% inv_C00 %*% t(post_mean_eta)
  e2 <-  T0 * (-D/2 * log(2 * pi) + D * 0.5 * log(L20) - 0.5 * Mu0 %*% tcrossprod(inv_C00, Mu0)) +
    sum(e20) + sum(e21) + sum(e22)

  #the X's
  diag(X %*% inv_C0 %*% t(X))
  e30 <- sweep(P, 1, -0.5 * diag(X %*% tcrossprod(inv_C0, X)), "*")
  e31 <- P * t(L21 %*% tcrossprod(inv_C0, X))
  e32 <- sweep(P, 2, -0.5 * diag(L21 %*% tcrossprod(inv_C0, L21)), "*")
  e33 <- sweep(P, 2, -0.5 * sum(diag(inv_C0))/L2, "*")
  e3 <- sum(P * (-0.5 * D * log(2 * pi) + 0.5 * determinant(inv_C0, logarithm=TRUE)$modulus))
  + sum(e30) + sum(e31) + sum(e32) + sum(e33)

  #the variational distribution
  e4 <- sum(D * log(L2)/2) - 0.5 * D * T0 * (log(2 * pi) + 1)

  return(c("e2"=e2, "e3"=e3, "me4"=-e4))
}
