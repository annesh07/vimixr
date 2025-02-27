ELBO_C0 <- function(N, D, T0, s1, s2, C00, nu0, V0, X,
                 W1, W2, L1, L2, Plog, nu, V){
  Mu0 <- matrix(c(rep(0,D)), nrow=1)
  #C00 <- diag(D)/L20
  Mu00 <- Mu0%*%solve(C00)
  #C0 <- diag(D)
  W1 <- W1
  W2 <- W2
  L1 <- matrix(L1, nrow = T0, ncol = D)
  L2 <- array(L2, c(D, D, T0))
  Plog <- matrix(Plog, nrow = N, ncol = T0)
  nu <- nu
  V <- matrix(V, nrow = D, ncol = D)

  P0 <- exp(Plog)
  #the alpha
  e0 <- s1*log(s2) - lgamma(s1) + (s1 - 1)*(digamma(W1)-log(W2)) - s2*(W1/W2)

  #the z's
  eni <- colSums(P0)
  vni <- colSums(P0*(1 - P0))
  enj <- rowSums(apply(P0, 1, f0))
  vnj <- rowSums(apply(P0, 1, f1))
  e10 <- lgamma(1 + eni) + 0.5*trigamma(1 + eni)*vni +
    lgamma((W1/W2) + enj) + 0.5*trigamma((W1/W2) + enj)*((W1/(W2^2)) + vnj) -
    lgamma(1 + (W1/W2)+eni+enj) -
    0.5*trigamma(1 + (W1/W2)+eni+enj)*((W1/(W2^2))+vni+vnj)
  e1 <- T0*(digamma(W1) - log(W2)) + sum(e10)

  #the eta's
  # e20 <- rep(NA, T0)
  inv_C00 <- solve(C00)
  L21 <- matrix(0, nrow = T0, ncol = D)
  for (i in 1:T0){
    L21[i,] = L1[i,, drop = FALSE] %*% solve(L2[,,i])
  }
  e20 <- diag(-0.5*L21 %*% inv_C00 %*% t(L21))
  e21 <- apply(L2, 3, function(x){-0.5*sum(diag(inv_C00 %*% solve(x)))})
  e22 <- Mu0 %*% inv_C00 %*% t(L21)
  e2 <-  T0*(-D/2*log(2*pi) - 0.5*determinant(C00, logarithm = TRUE)$modulus -
               0.5*Mu0 %*% inv_C00 %*% t(Mu0)) + sum(e20) + sum(e21) + sum(e22)
  # e2 <- sum(e20)

  #the X's
  #e30 <- matrix(NA, nrow = N, ncol = T0)
  inv_C0 <- nu*V
  e30 <- sweep(P0, 1, -0.5*diag(X %*% inv_C0 %*% t(X)), "*")
  e31 <- P0*t(L21 %*% inv_C0 %*% t(X))
  e32 <- sweep(P0, 2, -0.5*diag(L21 %*% inv_C0 %*% t(L21)), "*")
  e33 <- apply(L2, 3, function(x){-0.5*sum(diag(inv_C0 %*% solve(x)))})
  e34 <- sweep(P0, 2, e33, "*")
  e3 <- sum(P0*(-0.5*D*log(2*pi) + 0.5*(sum(digamma(0.5*(nu + 1 - c(1:D)))) +
        D*log(2) + determinant(V, logarithm = TRUE)$modulus)))
  + sum(e30) + sum(e31) + sum(e32) + sum(e34)

  #the C0
  inv_V0 <- solve(V0)
  e4 <- -0.5*nu0*D*log(2) - 0.25*D*(D-1)*log(pi)
    - sum(lgamma(0.5*(nu0 + 1 - c(1:D))))
    + 0.5*nu0*determinant(inv_V0, logarithm = TRUE)$modulus
    - 0.5*(nu0 + D + 1)*(sum(digamma(0.5*(nu + 1 - c(1:D)))) +
    D*log(2) + determinant(V, logarithm = TRUE)$modulus) -
    0.5*sum(diag(inv_V0 %*% inv_C0))

  #the variationa distributions
  e50 <- W1*log(W2) - lgamma(W1) + (W1-1)*(-log(W2) + digamma(W1)) - W1
  e51 <- sum(exp(Plog)*Plog)
  e52 <- sum(apply(L2, 3,
    function(x){0.5*determinant(x, logarithm = TRUE)$modulus}))
    - 0.5*D*T0*(log(2*pi)+1)
  inv_V <- solve(V)
  e53 <- -0.5*nu*D*log(2) - 0.25*D*(D-1)*log(pi)
    - sum(lgamma(0.5*(nu + 1 - c(1:D))))
    + 0.5*nu*determinant(inv_V, logarithm = TRUE)$modulus
    - 0.5*(nu + D + 1)*(sum(digamma(0.5*(nu + 1 - c(1:D)))) +
    D*log(2) + determinant(V, logarithm = TRUE)$modulus) - 0.5*nu*D
  e5 <- e50 + e51 + e52 + e53

  return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "e4"=e4, "me5"=-e5))
}
