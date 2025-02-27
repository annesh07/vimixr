ELBO_InvWishart0 <- function(N, D, T0, C00, nu0, Lamda0, k0, Lamda, k,
                            X, L1, L2, nu, V, Plog){
  Mu0 <- matrix(c(rep(0,D)), nrow=1)

  k <- k
  Lamda <- matrix(Lamda, nrow = D, ncol = D)
  nu <- nu
  V <- matrix(V, nrow = D, ncol = D)
  # W1 <- W1
  # W2 <- W2
  L1 <- matrix(L1, nrow = T0, ncol = D)
  L2 <- array(L2, c(D, D, T0))
  Plog <- matrix(Plog, nrow = N, ncol = T0)
  P0 <- exp(Plog)

  # #the alpha
  # e0 <- s1*log(s2) - lgamma(s1) + (s1 - 1)*(digamma(W1)-log(W2)) - s2*(W1/W2)

  #the z's
  eni <- colSums(P0)
  vni <- colSums(P0*(1 - P0))
  enj <- rowSums(apply(P0, 1, f0))
  vnj <- rowSums(apply(P0, 1, f1))
  e10 <- lgamma(1 + eni) + 0.5*trigamma(1 + eni)*vni +
    lgamma(0.4 + enj) + 0.5*trigamma(0.4 + enj)*(vnj) -
    lgamma(1 + 0.4 + eni + enj) -
    0.5*trigamma(1 + 0.4 + eni + enj)*(vni+vnj)
  e1 <- T0*log(0.4) + sum(e10)

  #the eta_i's
  inv_C00 <- solve(C00)
  L21 <- matrix(0, nrow = T0, ncol = D)
  for (i in 1:T0){
    L21[i,] = L1[i,, drop = FALSE] %*% solve(L2[,,i])
  }
  e20 <- -0.5*sum(diag(L21 %*% inv_C00 %*% t(L21)))
  e21 <- sum(L21 %*% inv_C00 %*% t(Mu0))
  e2 <- T0*(-0.5*D*log(2*pi) - 0.5*determinant(C00, logarithm = TRUE)$modulus
            - 0.5*Mu0 %*% inv_C00 %*% t(Mu0)) + e20 + e21

  #the x_n's
  EC0_inv <- nu*V
  e30 <- sweep(P0, 1, -0.5*(diag(X %*% EC0_inv %*% t(X))), "*")
  e31 <- P0*t(L21 %*% EC0_inv %*% t(X))
  e32 <- sweep(P0, 2, -0.5*diag(L21 %*% EC0_inv %*% t(L21)), "*")
  e33 <- apply(L2, 3, function(x){-0.5*sum(diag(EC0_inv %*% solve(x)))})
  e34 <- sweep(P0, 2, e33, "*")
  e3 <- sum(P0*(-0.5*D*log(2*pi) + 0.5*(sum(digamma(0.5*(nu + 1 - c(1:D)))) +
                                          D*log(2) + determinant(V, logarithm = TRUE)$modulus))) +
    sum(e30) + sum(e31) + sum(e32) + sum(e34)

  #the V0_inv
  e4 <- -0.5*k0*D*log(2) - 0.5*k0*determinant(Lamda0, logarithm = TRUE)$modulus -
    0.25*D*(D-1)*log(pi) - sum(lgamma(0.5*(k0 + 1 - c(1:D)))) +
    0.5*(k0 - D - 1)*(sum(digamma(0.5*(k + 1 - c(1:D)))) + D*log(2) +
                        determinant(Lamda, logarithm = TRUE)$modulus) -
    0.5*sum(diag(solve(Lamda0) %*% (k*Lamda)))

  #the C0_inv
  e5 <- -0.5*nu0*D*log(2) + 0.5*nu0*(sum(digamma(0.5*(k + 1 - c(1:D)))) + D*log(2) +
                                       determinant(Lamda, logarithm = TRUE)$modulus) -
    0.25*D*(D-1)*log(pi) - sum(lgamma(0.5*(nu0 + 1 - c(1:D)))) +
    0.5*(nu0 + D + 1)*(sum(digamma(0.5*(nu + 1 - c(1:D)))) +
                         D*log(2) + determinant(V, logarithm = TRUE)$modulus) -
    0.5*sum(diag((k*Lamda) %*% EC0_inv))

  #the variationa distributions
  # e60 <- W1*log(W2) - lgamma(W1) + (W1-1)*(digamma(W1) - log(W2)) - W1
  e61 <- sum(exp(Plog)*Plog)
  e62 <- sum(apply(L2, 3,
                   function(x){0.5*determinant(x, logarithm = TRUE)$modulus}))
  - 0.5*D*T0*(log(2*pi)+1)
  e63 <- -0.5*nu*D*log(2) - 0.5*nu*determinant(V, logarithm = TRUE)$modulus -
    0.25*D*(D-1)*log(pi) - sum(lgamma(0.5*(nu + 1 - c(1:D)))) +
    0.5*(nu + D + 1)*(sum(digamma(0.5*(nu + 1 - c(1:D)))) +
                        D*log(2) + determinant(V, logarithm = TRUE)$modulus) -
    0.5*D*nu
  e64 <- -0.5*k*D*log(2) - 0.5*k*determinant(Lamda, logarithm = TRUE)$modulus -
    0.25*D*(D-1)*log(pi) - sum(lgamma(0.5*(k + 1 - c(1:D)))) +
    0.5*(k - D - 1)*(sum(digamma(0.5*(k + 1 - c(1:D)))) + D*log(2) +
                       determinant(Lamda, logarithm = TRUE)$modulus) -
    0.5*D*k
  e6 <- e61 + e62 + e63 + e64

  return(c("e1"=e1, "e2"=e2, "e3"=e3, "e4"=e4, "e5"=e5, "me6"=-e6))
}
