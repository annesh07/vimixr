ELBO_C0_i <- function(N, D, T0, s1, s2, Mu0, nu0, V0, k0, X,
                   W1, W2, L1, nu1, V1, Plog){
  P0 <- exp(Plog)
  RP0 <- colsums(P0)

  V1_inv <- array(apply(V1, 3, function(x){spdinv(x)}), dim = dim(V1))
  E_C0 <- sweep(V1_inv, 3, nu1, "*")
  E_log_C0_1 <- apply(nu1, 2, function(x){sum(digamma(0.5*(x + 1 - c(1:D))))})
  E_log_C0_2 <- apply(V1_inv, 3, function(x){D*log(2) +
      determinant(x, logarithm = TRUE)$modulus})
  E_log_C0 <- matrix((E_log_C0_1 + E_log_C0_2), nrow = 1, ncol = T0)

  L2 <- sweep(V1, 3, nu1*(1/k0 + RP0), "/")

  #the alpha
  e0 <- s1*log(s2) - lgamma(s1) + (s1 - 1)*(digamma(W1)-log(W2)) - s2*(W1/W2)

  #the z's
  eni <- colsums(P0)
  vni <- colsums(P0*(1 - P0))
  enj <- rowsums(apply(P0, 1, f0))
  vnj <- rowsums(apply(P0, 1, f1))
  e10 <- lgamma(1 + eni) + 0.5*trigamma(1 + eni)*vni +
    lgamma((W1/W2) + enj) + 0.5*trigamma((W1/W2) + enj)*((W1/(W2^2)) + vnj) -
    lgamma(1 + (W1/W2)+eni+enj) -
    0.5*trigamma(1 + (W1/W2)+eni+enj)*((W1/(W2^2))+vni+vnj)
  e1 <- T0*(digamma(W1) - log(W2)) + sum(e10)

  #the eta_i's and C0_i's
  e20 <- rep(0, T0)
  e21 <- rep(0, T0)
  for (i in 1:T0){
    e20[i] <- -(1/(2*k0))*L1[i,,drop=FALSE] %*% E_C0[,,i] %*% t(L1[i,,drop=FALSE])
    e21[i] <- (1/k0)*L1[i,,drop=FALSE] %*% E_C0[,,i] %*% t(Mu0)
  }
  e22 <- apply(E_C0, 3, function(x){-(1/(2*k0))*Mu0 %*% x %*% t(Mu0)})
  e23 <- apply(E_C0, 3, function(x){-0.5*sum(t(V0)*x)})
  e2 <- -0.5*D*T0*log(2*pi) + 0.5*(nu0 + D + 2)*sum(E_log_C0) + sum(e20) -
    (1/(2*k0))*D*sum(1/(1/k0 + RP0)) + sum(e21) + sum(e22) - 0.5*T0*nu0*D*log(2) -
    T0*(0.25*D*(D-1)*log(pi) + sum(lgamma(0.5*(nu0 + 1 - c(1:D))))) +
    0.5*T0*nu0*determinant(V0, logarithm = TRUE)$modulus + sum(e23)

  #the data x_n's
  e30 <- matrix(0, nrow = N, ncol = T0)
  e31 <- matrix(0, nrow = N, ncol = T0)
  for (n in 1:N){
    for (i in 1:T0){
      e30[n,i] <- -0.5*P0[n,i]*X[n,,drop=FALSE] %*% E_C0[,,i] %*% t(X[n,,drop=FALSE])
      e31[n,i] <- P0[n,i]*L1[i,,drop=FALSE] %*% E_C0[,,i] %*% t(X[n,,drop=FALSE])
    }
  }
  e3 <- -N*0.5*D*log(2*pi) + sum(RP0*0.5*E_log_C0) + sum(e30) + sum(e31) +
    sum(RP0*e20) - 0.5*D*sum(RP0/(1/k0 + RP0))

  #the variationa distributions
  e40 <- W1*log(W2) - lgamma(W1) + (W1-1)*(-log(W2) + digamma(W1)) - W1
  e41 <- sum(exp(Plog)*Plog)
  e420 <- apply(L2, 3, function(x){-0.5*determinant(x, logarithm = TRUE)$modulus})
  e421 <- apply(nu1, 2, function(x){-0.25*D*(D-1)*log(pi) -
      sum(lgamma(0.5*(x + 1 - c(1:D))))})
  e422 <- apply(V1, 3, function(x){determinant(x, logarithm = TRUE)$modulus})
  e42 <- -0.5*D*T0*log(2*pi) + sum(e420) - 0.5*D*T0 - 0.5*D*log(2)*sum(nu1) +
    sum(e421) + 0.5*sum(nu1*e422) + 0.5*sum((nu1 + D + 1)*E_log_C0) - 0.5*D*sum(nu1)

  e4 <- e40 + e41 + e42

  return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "me4"=-e4))
}
