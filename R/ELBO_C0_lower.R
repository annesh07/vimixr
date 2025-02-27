ELBO_C0_lower <- function(N, D, T0, s1, s2, C00, a0, b0, mu0, c0, X,
                          W1, W2, L1, L2, Plog, a1, b1, mu1, c1){
  Mu0 <- matrix(c(rep(0,D)), nrow=1)
  #C00 <- diag(D)/L20
  inv_C00 <- spdinv(C00)
  Mu00 <- mat.mult(Mu0, inv_C00)
  #C0 <- diag(D)
  W1 <- W1
  W2 <- W2
  L1 <- matrix(L1, nrow = T0, ncol = D)
  L2 <- array(L2, c(D, D, T0))
  Plog <- matrix(Plog, nrow = N, ncol = T0)
  a1 <- matrix(a1, nrow = 1, ncol = D)
  b1 <- matrix(b1, nrow = 1, ncol = D)
  mu1 <- matrix(mu1, nrow = 1, ncol = D*(D-1)/2)
  c1 <- matrix(c1, nrow = 1, ncol = D*(D-1)/2)

  mean_lower <- matrix(0, nrow = D, ncol = D)
  mean_lower[lower.tri(mean_lower, diag = FALSE)] <- mu1
  sigma_lower <- matrix(0, nrow = D, ncol = D)
  sigma_lower[lower.tri(sigma_lower, diag = FALSE)] <- c1

  mean_L <- mean_lower + diag(sqrt(1/b1)*sqrt(pi)/beta(a1,0.5))
  diag(sigma_lower) <- (1/b1)*(a1 - (sqrt(pi)/beta(a1,0.5))^2)
  E_C0 <- mat.mult(mean_L, t(mean_L)) + diag(rowsums(sigma_lower))
  #diag(E_C0) <- a1/b1 + rowsums((mean_lower)^2 + sigma_lower)
  #E_C0 <- solve(E_C0)
  P0 <- exp(Plog)
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

  #the eta's
  # e20 <- rep(NA, T0)
  L21 <- matrix(0, nrow = T0, ncol = D)
  for (i in 1:T0){
    L21[i,] = mat.mult(L1[i,, drop = FALSE], L2[,,i])
  }
  e20 <- diag(-0.5*mat.mult(mat.mult(L21, inv_C00), t(L21)))
  e21 <- apply(L2, 3, function(x){-0.5*sum(diag(mat.mult(inv_C00, x)))})
  e22 <- mat.mult(mat.mult(Mu0, inv_C00), t(L21))
  e2 <-  T0*(-D/2*log(2*pi) - 0.5*determinant(C00, logarithm = TRUE)$modulus -
               0.5*mat.mult(mat.mult(Mu0, inv_C00), t(Mu0))) + sum(e20) + sum(e21) + sum(e22)

  #the X's
  #e30 <- matrix(NA, nrow = N, ncol = T0)
  e30 <- eachcol.apply(P0, -0.5*diag(mat.mult(mat.mult(X, E_C0), t(X))),
                       oper = "*")
  e31 <- P0*t(mat.mult(mat.mult(L21, E_C0), t(X)))
  e32 <- eachrow(P0, -0.5*diag(mat.mult(mat.mult(L21, E_C0), t(L21))), oper = "*")
  e33 <- apply(L2, 3, function(x){-0.5*sum(diag(mat.mult(E_C0, x)))})
  e34 <- eachrow(P0, e33, oper = "*")
  e3 <- sum(P0*(-0.5*D*log(2*pi) + 0.5*sum(digamma(a1) - log(b1)))) + sum(e30) +
    sum(e31) + sum(e32) + sum(e34)

  #the C0
  # e40 <- sum(a0*log(b0) - lgamma(a0) + (a0 - 1)*(digamma(a1) - log(b1)) -
  #              a1*b0/b1)
  e40 <- sum(log(2) + a0*log(b0) - lgamma(a0) -
               (2*a0 - 1)*(log(sqrt(b1)*sqrt(pi)/beta(a1,0.5)) -
                             0.5*(b1*(a1-(sqrt(pi)/beta(a1,0.5))^2))/(a1*b1)) -
               a1*b1*b0)
  e41 <- sum(-0.5*log(2*pi*c0) - 0.5*(mu1^2 + c1 - 2*mu0*mu1 + mu0^2)/c0)
  e4 <- e40 + e41

  #the variationa distributions
  e50 <- W1*log(W2) - lgamma(W1) + (W1-1)*(-log(W2) + digamma(W1)) - W1
  e51 <- sum(exp(Plog)*Plog)
  e52 <- sum(apply(L2, 3,
                   function(x){-0.5*determinant(x, logarithm = TRUE)$modulus})) -
    0.5*D*T0*(log(2*pi)+1)
  # e530 <- sum(a1*log(b1) - lgamma(a1) + (a1 - 1)*(digamma(a1) - log(b1)) - a1)
  e530 <- sum(log(2) - a1*log(b1) - lgamma(a1) -
                (2*a1 - 1)*(log(sqrt(b1)*sqrt(pi)/beta(a1,0.5)) -
                              0.5*(b1*(a1-(sqrt(pi)/beta(a1,0.5))^2))/(a1*b1)) -
                a1)
  e531 <- sum(-0.5*(log(2*pi*c1) + 1))
  e5 <- e50 + e51 + e52 + e530 + e531

  return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "e4"=e4, "me5"=-e5))
}
