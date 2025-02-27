ELBO_sigma1 <- function(N, D, T0, s1, s2, L20, b1, b2, X, W1, W2, L1, L2, G1, G2,
                       Plog){
  Mu0 <- matrix(c(rep(0,D)), nrow=1)
  C00 <- sweep(diag(D), 1, L20, "/")
  inv_C00 <- solve(C00)

  G1 <- matrix(G1, nrow = 1, ncol = D)
  G2 <- matrix(G2, nrow = 1, ncol = D)
  W1 <- W1
  W2 <- W2
  L1 <- matrix(L1, nrow = T0, ncol = D)
  L2 <- matrix(L2, nrow = T0, ncol = D)
  Plog <- matrix(Plog, nrow = N, ncol = T0)

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
  L21 <- L1/L2
  e20 <- diag(-0.5*L21 %*% inv_C00 %*% t(L21))
  e21 <- - 0.5*sum(sweep((1/L2), 2, L20, "*"))
  e22 <- Mu0 %*% inv_C00 %*% t(L21)
  e2 <-  T0*(-D/2*log(2*pi) + 0.5*sum(log(L20)) - 0.5*Mu0 %*% inv_C00 %*% t(Mu0))
    + sum(e20) + sum(e21) + sum(e22)

  #the X's
  e30 <- sweep(P0, 1, -0.5*rowSums(sweep(X^2, 2, (G1/G2), "*")), "*")
  e31 <- P0*t(sweep(L21, 2, (G1/G2), "*") %*% t(X))
  e32 <- sweep((L21^2 + 1/L2), 2 , (G1/G2), "*")
  e33 <- sweep(P0, 2, -0.5*rowSums(e32), "*")
  e3 <- sum(P0*(-0.5*D*log(2*pi) + 0.5*sum(digamma(G1) - log(G2))))
  + sum(e30) + sum(e31) + sum(e33)

  #the sigma^2
  e4 <- sum(b1*log(b2) - lgamma(b1) + (b1 - 1)*(digamma(G1) - log(G2)) - b2*G1/G2)

  #the variational distributions
  e50 <- W1*log(W2) - lgamma(W1) + (W1-1)*(digamma(W1) - log(W2)) - W1
  e51 <- sum(exp(Plog)*Plog)
  e52 <- sum(log(L2)/2) - 0.5*D*T0*(log(2*pi)+1)
  e53 <- sum(G1*log(G2) - lgamma(G1) + (G1-1)*(digamma(G1) - log(G2)) - G1)
  e5 <- e50 + e51 + e52 + e53

  return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "e4"=e4, "me5"=-e5))
}
