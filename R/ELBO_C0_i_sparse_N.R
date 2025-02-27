ELBO_C0_i_sparse_N <- function(N, D, T0, s1, s2, Mu0, a0, b0, k0, X,
                             W1, W2, L1, a1, B1, C1, Plog){
  P0 <- exp(Plog)
  RP0 <- colsums(P0)

  E_C0_inv <- array(0, c(D, D, T0))
  E_C0_inv_inv <- array(0, c(D, D, T0))
  for (i in 1:T0){
    E_C0_inv[,,i] <- Diag.fill(C1[,,i], a1[1,i]/B1[i,])
    E_C0_inv_inv[,,i] <- spdinv(E_C0_inv[,,i])
  }
  L2 <- sweep(E_C0_inv_inv, 3, (1/k0 + RP0), "/")

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
  e200 <- rep(0, T0)
  e201 <- e200
  e202 <- e200
  for (i in 1:T0){
    e200[i] <- 0.5*sum(digamma(a1[1,i]) - log(B1[i,]))
    e201[i] <- -(1/(2*k0))*L1[i,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(L1[i,,drop=FALSE])
    e202[i] <- (1/k0)*L1[i,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(Mu0)
  }
  e203 <- apply(E_C0_inv, 3, function(x){-(1/(2*k0))*Mu0 %*% x %*% t(Mu0)})
  e20 <- -0.5*T0*log(2*pi) - 0.5*log(k0) + sum(e200) + sum(e201) -
    (0.5/k0)*D/sum(1/k0 + RP0) + sum(e202) + sum(e203)

  e210 <- rep(0, T0)
  for (i in 1:T0){
    e210[i] <- sum(a0*log(b0) - lgamma(a0) +
                     (a0 - 1)*(digamma(a1[1,i]) - log(B1[i,])) - b0*a1[1,i]/B1[i,])
  }
  e21 <- sum(e210) + sum(-0.5*log(2*pi) - 0.5*(1 + C1[!diag(D)]^2))
  e2 <- e20 + e21

  #the data X
  e30 <- matrix(0, nrow = N, ncol = T0)
  e31 <- matrix(0, nrow = N, ncol = T0)
  for (n in 1:N){
    for (i in 1:T0){
      e30[n,i] <- -0.5*P0[n,i]*X[n,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(X[n,,drop=FALSE])
      e31[n,i] <- P0[n,i]*L1[i,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(X[n,,drop=FALSE])
    }
  }
  e3 <- -N*0.5*D*log(2*pi) + sum(RP0*e200) + sum(e30) + sum(e31) +
    sum(RP0*k0*e201) - 0.5*D*sum(RP0/(1/k0 + RP0))

  #the variationa distributions
  e40 <- W1*log(W2) - lgamma(W1) + (W1-1)*(-log(W2) + digamma(W1)) - W1
  e41 <- sum(exp(Plog)*Plog)
  e420 <- rep(0, T0)
  for (i in 1:T0){
    e420[i] <- sum(a1[1,i]*log(B1[i,]) - lgamma(a1[1,i]) +
                     (a1[1,i] -1)*(digamma(a1[1,i]) - log(B1[i,])) - a1[1,i])
  }
  e42 <- sum(e420) + (D^2 - D)*(-0.5*log(2*pi) - 0.5)
  e4 <- e40 + e41 + e42

  return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "me4"=-e4))
}


