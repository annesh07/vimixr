ELBO_0 <- function(N, D, T0, s1, s2, C00, C0, X,
                           W1, W2, L1, L2, Plog){
  Mu0 <- matrix(c(rep(0,D)), nrow=1)
  inv_C00 <- spdinv(C00)
  Mu00 <- mat.mult(Mu0, inv_C00)
  E_C0 <- spdinv(C0)
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
  L21 <- matrix(0, nrow = T0, ncol = D)
  for (i in 1:T0){
    L21[i,] = mat.mult(L1[i,, drop = FALSE], L2[,,i])
  }
  e20 <- apply(L21, 1, function(x){-0.5*t(x)%*%inv_C00%*%x})
  e21 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C00)*x)})
  e22 <- mat.mult(mat.mult(Mu0, inv_C00), t(L21))
  e2 <-  T0*(-D/2*log(2*pi) - 0.5*determinant(C00, logarithm = TRUE)$modulus -
               0.5*mat.mult(mat.mult(Mu0, inv_C00), t(Mu0))) +
    sum(e20) + sum(e21) + sum(e22)

  #the X's
  e30 <- eachcol.apply(P0, apply(X, 1, function(x){-0.5*t(x)%*%E_C0%*%x}),
                       oper = "*")
  e31 <- P0*t(mat.mult(mat.mult(L21, E_C0), t(X)))
  e32 <- P0*matrix(apply(L21, 1, function(x){-0.5*t(x)%*%E_C0%*%x}),
                   nrow = N, ncol = T0, byrow = TRUE)
  e33 <- apply(L2, 3, function(x){-0.5*sum(t(E_C0)*x)})
  e34 <- P0*matrix(e33, nrow = N, ncol = T0, byrow=TRUE)
  e3 <- N*(-0.5*D*log(2*pi) + 0.5*determinant(E_C0, logarithm = TRUE)$modulus) +
    sum(e30) + sum(e31) + sum(e32) + sum(e34)

  #the variationa distributions
  e40 <- W1*log(W2) - lgamma(W1) + (W1-1)*(-log(W2) + digamma(W1)) - W1
  e41 <- sum(exp(Plog)*Plog)
  e42 <- apply(L2, 3,function(x){-0.5*determinant(x, logarithm = TRUE)$modulus})

  e4 <- e40 + e41 + (sum(e42) - 0.5*D*T0*(log(2*pi)+1))

  return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "me4"=-e4))
}
