
elbo_fixed_diagonal <- function(X, inverts, params){
  N <- params$N
  D <- params$D
  T0 <- params$T0
  Mu0 <- params$prior_mean_eta
  P <- params$P
  RP <- Rfast::colsums(P)
  L1 <- params[["post_mean_eta"]]
  L2 <- params[["post_precision_scalar_eta"]]
  L20 <- params[["prior_precision_scalar_eta"]]
  inv_C0 <- inverts[["inv_C0"]]
  inv_C00 <- inverts[["inv_C00"]]

  #the eta's
  L21 <- sweep(L1, 1, L2, "/")
  e20 <- -0.5*quadratic_form_diag(L21, inv_C00)
  e21 <- -0.5 * (D * L20)/L2
  e22 <- Rfast::mat.mult(Mu0, Rfast::Tcrossprod(inv_C00, L21))
  e2 <-  T0 * (-D/2 * log(2 * pi) + D * 0.5 * log(L20) -
                 0.5 * Rfast::mat.mult(Mu0, Rfast::Tcrossprod(inv_C00, Mu0))) +
    sum(e20) + sum(e21) + sum(e22)

  #the X's
  e30 <- sweep(P, 1, -0.5 * quadratic_form_diag(X, inv_C0), "*")
  e31 <- P * t(Rfast::mat.mult(L21, Rfast::Tcrossprod(inv_C0, X)))
  e32 <- sweep(P, 2, -0.5 * quadratic_form_diag(L21, inv_C0), "*")
  e33 <- sweep(P, 2, -0.5 * sum(diag(inv_C0))/L2, "*")
  e3 <- sum(P * (-0.5 * D * log(2 * pi) + 0.5 * determinant(inv_C0, logarithm=TRUE)$modulus))
  + sum(e30) + sum(e31) + sum(e32) + sum(e33)

  #the variational distribution
  e4 <- sum(D * log(L2)/2) - 0.5 * D * T0 * (log(2 * pi) + 1)

  return(c("e_mean"=e2, "e_data"=e3, "me_var"=-e4))
}

elbo_varied_diagonal <- function(X, inverts, params){
  N <- params$N
  D <- params$D
  T0 <- params$T0
  Mu0 <- params$prior_mean_eta
  P <- params$P
  RP <- Rfast::colsums(P)
  L1 <- params$post_mean_eta
  b1 <- params$prior_shape_scalar_cov
  b2 <- params$prior_rate_scalar_cov
  G1 <- params$post_shape_scalar_cov
  G2 <- params$post_rate_scalar_cov
  L2 <- params$post_precision_scalar_eta
  L20 <- params$prior_precision_scalar_eta
  inv_C00 <- inverts[["inv_C00"]]    #inverse of C00

  #the eta's
  L21 <- sweep(L1, 1, L2, "/")
  e20 <- -0.5 * quadratic_form_diag(L21, inv_C00)
  e21 <- - 0.5*(D*L20)/L2
  e22 <- Rfast::mat.mult(Mu0, Rfast::Tcrossprod(inv_C00, L21))
  e2 <-  T0*(-D/2*log(2*pi) + D*0.5*log(L20) -
               0.5*Rfast::mat.mult(Mu0, Rfast::Tcrossprod(inv_C00, Mu0))) +
    sum(e20) + sum(e21) + sum(e22)

  #the X's
  e30 <- sweep(P, 1, -0.5*(G1/G2)*diag(Rfast::Tcrossprod(X)), "*")
  e31 <- P*t((G1/G2)*(Rfast::Tcrossprod(L21, X)))
  e32 <- sweep(P, 2, -0.5*(G1/G2)*diag(Rfast::Tcrossprod(L21)), "*")
  e33 <- sweep(P, 2, -0.5*(G1/G2)*D/L2, "*")
  e3 <- sum(P*(-0.5*D*log(2*pi) + 0.5*D*(digamma(G1) - log(G2))))
  + sum(e30) + sum(e31) + sum(e32) + sum(e33)

  #the sigma^2
  e4 <- b1*log(b2) - lgamma(b1) + (b1 - 1)*(digamma(G1) - log(G2)) - b2*G1/G2

  #the variational distributions
  e52 <- sum(D*log(L2)/2) - 0.5*D*T0*(log(2*pi)+1)
  e53 <- G1*log(G2) - lgamma(G1) + (G1-1)*(digamma(G1) - log(G2)) - G1
  e5 <- e52 + e53

  return(c("e_mean"=e2, "e_data"=e3, "e_cov"=e4, "me_var"=-e5))
}

elbo_fixed_full <- function(X, inverts, params){
  N <- params$N
  D <- params$D
  T0 <- params$T0
  Mu0 <- params$prior_mean_eta
  P <- params$P
  RP <- Rfast::colsums(P)
  L1 <- params$post_mean_eta
  C0 <- params$cov_data
  L2 <- params$post_cov_eta
  C00 <- params$prior_cov_eta
  inv_C0 <- inverts[["inv_C0"]]      #inverse of C00
  inv_C00 <- inverts[["inv_C00"]]    #inverse of C0

  #the eta's
  L21 <- matrix(0, nrow = T0, ncol = D)
  for (i in 1:T0){
    L21[i,] = Rfast::mat.mult(L1[i,, drop = FALSE], L2[,,i])
  }
  e20 <- -0.5 * quadratic_form_diag(L21, inv_C00)
  e21 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C00)*x)})
  e22 <- Rfast::mat.mult(Mu0, Rfast::Tcrossprod(inv_C00, L21))
  e2 <-  T0*(-D/2*log(2*pi) - 0.5*determinant(C00, logarithm = TRUE)$modulus -
               0.5*Rfast::mat.mult(Mu0, Rfast::Tcrossprod(inv_C00, Mu0))) +
    sum(e20) + sum(e21) + sum(e22)

  #the X's
  e30 <- sweep(P, 1, -0.5 * quadratic_form_diag(X, inv_C0), "*")
  e31 <- P*t(Rfast::mat.mult(L21, Rfast::Tcrossprod(inv_C0, X)))
  e32 <- sweep(P, 2, -0.5 * quadratic_form_diag(L21, inv_C0), "*")
  e33 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C0)*x)})
  e34 <- P*matrix(e33, nrow = N, ncol = T0, byrow=TRUE)
  e3 <- N*(-0.5*D*log(2*pi) + 0.5*determinant(inv_C0, logarithm = TRUE)$modulus) +
    sum(e30) + sum(e31) + sum(e32) + sum(e34)

  #the variational distribution
  e42 <- apply(L2, 3,function(x){-0.5*determinant(x, logarithm = TRUE)$modulus})
  e4 <- sum(e42) - 0.5*D*T0*(log(2*pi)+1)

  return(c("e_mean"=e2, "e_data"=e3, "me_var"=-e4))
}

elbo_varied_IW_full <- function(X, inverts, params){
  N <- params$N
  D <- params$D
  T0 <- params$T0
  Mu0 <- params$prior_mean_eta
  P <- params$P
  RP <- Rfast::colsums(P)
  L1 <- params$post_mean_eta
  nu0 <- params$prior_df_cov
  V0 <- params$prior_scale_cov
  nu <- params$post_df_cov
  V <- params$post_scale_cov
  L2 <- params$post_cov_eta
  C00 <- params$prior_cov_eta
  inv_C00 <- inverts[["inv_C00"]]    #inverse of C00
  inv_V0 <- inverts[["inv_V0"]]      #inverse of V0

  inv_C0 <- nu*V           #expected inverse of C0; covariance matrix of data

  #the eta's
  L21 <- matrix(0, nrow = T0, ncol = D)
  for (i in 1:T0){
    L21[i,] = Rfast::mat.mult(L1[i,, drop = FALSE], L2[,,i])
  }
  e20 <- -0.5 * quadratic_form_diag(L21, inv_C00)
  e21 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C00) * x)})
  e22 <- Rfast::mat.mult(Mu0, Rfast::Tcrossprod(inv_C00, L21))
  e2 <-  T0*(-D/2*log(2*pi) - 0.5*determinant(C00, logarithm = TRUE)$modulus -
               0.5*Rfast::mat.mult(Mu0, Rfast::Tcrossprod(inv_C00, Mu0))) +
    sum(e20) + sum(e21) + sum(e22)

  #the X's
  e30 <- sweep(P, 1, -0.5 * quadratic_form_diag(X, inv_C0), "*")
  e31 <- P*t(Rfast::mat.mult(L21, Rfast::Tcrossprod(inv_C0, X)))
  e32 <- sweep(P, 2, -0.5 * quadratic_form_diag(L21, inv_C0), "*")
  e33 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C0)*x)})
  e34 <- sweep(P, 2, e33, "*")
  e3 <- sum(P*(-0.5*D*log(2*pi) + 0.5*(sum(digamma(0.5*(nu + 1 - c(1:D)))) +
                                         D*log(2) +
                                         determinant(V, logarithm = TRUE)$modulus)))
  + sum(e30) + sum(e31) + sum(e32) + sum(e34)

  #the C0
  e4 <- -0.5*nu0*D*log(2) - 0.25*D*(D-1)*log(pi) -
    sum(lgamma(0.5*(nu0 + 1 - c(1:D)))) +
    0.5*nu0*determinant(inv_V0, logarithm = TRUE)$modulus -
    0.5*(nu0 + D + 1)*(sum(digamma(0.5*(nu + 1 - c(1:D)))) +
                         D*log(2) + determinant(V, logarithm = TRUE)$modulus) -
    0.5*sum(t(inv_V0) * inv_C0)

  #the variationa distributions
  e52 <- sum(apply(L2, 3,
                   function(x){-0.5*determinant(x, logarithm = TRUE)$modulus}))
  - 0.5*D*T0*(log(2*pi)+1)
  e53 <- -0.5*nu*D*log(2) - 0.25*D*(D-1)*log(pi) -
    sum(lgamma(0.5*(nu + 1 - c(1:D)))) -
    0.5*nu*determinant(V, logarithm = TRUE)$modulus -
    0.5*(nu + D + 1)*(sum(digamma(0.5*(nu + 1 - c(1:D)))) +
                        D*log(2) + determinant(V, logarithm = TRUE)$modulus) -
    0.5*nu*D
  e5 <- e52 + e53

  return(c("e_mean"=e2, "e_data"=e3, "e_cov"=e4, "me_var"=-e5))
}

elbo_varied_decomposed_full <- function(X, inverts, params){
  N <- params$N
  D <- params$D
  T0 <- params$T0
  Mu0 <- params$prior_mean_eta
  P <- params$P
  RP <- Rfast::colsums(P)
  a0 <- params$prior_shape_diag_decomp
  b0 <- params$prior_rate_diag_decomp
  mu0 <- params$prior_mean_offdiag_decomp
  c0 <- params$prior_var_offdiag_decomp
  a1 <- params$post_shape_diag_decomp
  b1 <- params$post_rate_diag_decomp
  mu1 <- params$post_mean_offdiag_decomp
  c1 <- params$post_var_offdiag_decomp
  L2 <- params$post_cov_eta
  C00 <- params$prior_cov_eta
  inv_C00 <- inverts[["inv_C00"]]    #inverse of C00

  mean_lower <- matrix(0, nrow = D, ncol = D) #mean matrix of the decomposed
  mean_lower[lower.tri(mean_lower, diag = FALSE)] <- mu1
  sigma_lower <- matrix(0, nrow = D, ncol = D) #var matrix of the decomposed
  sigma_lower[lower.tri(sigma_lower, diag = FALSE)] <- c1
  mean_L <- mean_lower + diag(sqrt(1/b1)*sqrt(pi)/beta(a1,0.5))
  diag(sigma_lower) <- (1/b1)*(a1 - (sqrt(pi)/beta(a1,0.5))^2)
  #expected inverse of C0; covariance matrix of data
  inv_C0 <- Rfast::Tcrossprod(mean_L) + diag(rowsums(sigma_lower))

  #the eta's
  L21 <- matrix(0, nrow = T0, ncol = D)
  for (i in 1:T0){
    L21[i,] = Rfast::mat.mult(L1[i,, drop = FALSE], L2[,,i])
  }
  e20 <- -0.5 * quadratic_form_diag(L21, inv_C00)
  e21 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C00) * x)})
  e22 <- Rfast::mat.mult(Mu0, Rfast::Tcrossprod(inv_C00, L21))
  e2 <-  T0*(-D/2*log(2*pi) - 0.5*determinant(C00, logarithm = TRUE)$modulus -
               0.5*Rfast::mat.mult(Mu0, Rfast::Tcrossprod(inv_C00, Mu0))) +
    sum(e20) + sum(e21) + sum(e22)

  #the X's
  e30 <- sweep(P, 1, -0.5 * quadratic_form_diag(X, inv_C0), "*")
  e31 <- P*t(Rfast::mat.mult(L21, Rfast::Tcrossprod(inv_C0, X)))
  e32 <- sweep(P, 2, -0.5 * quadratic_form_diag(L21, inv_C0), "*")
  e33 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C0) * x)})
  e34 <- sweep(P, 2, e33, "*")
  e3 <- sum(P*(-0.5*D*log(2*pi) + 0.5*sum(digamma(a1) - log(b1)))) + sum(e30) +
    sum(e31) + sum(e32) + sum(e34)

  #the C0
  e40 <- sum(log(2) + a0*log(b0) - lgamma(a0) -
               (2*a0 - 1)*(log(sqrt(b1)*sqrt(pi)/beta(a1,0.5)) -
                             0.5*(b1*(a1-(sqrt(pi)/beta(a1,0.5))^2))/(a1*b1)) -
               a1*b1*b0)
  e41 <- sum(-0.5*log(2*pi*c0) - 0.5*(mu1^2 + c1 - 2*mu0*mu1 + mu0^2)/c0)
  e4 <- e40 + e41

  #the variationa distributions
  e52 <- sum(apply(L2, 3,
                   function(x){-0.5*determinant(x, logarithm = TRUE)$modulus})) -
    0.5*D*T0*(log(2*pi)+1)
  e530 <- sum(log(2) - a1*log(b1) - lgamma(a1) -
                (2*a1 - 1)*(log(sqrt(b1)*sqrt(pi)/beta(a1,0.5)) -
                              0.5*(b1*(a1-(sqrt(pi)/beta(a1,0.5))^2))/(a1*b1)) -
                a1)
  e531 <- sum(-0.5*(log(2*pi*c1) + 1))
  e5 <- e52 + e530 + e531

  return(c("e_mean"=e2, "e_data"=e3, "e_cov"=e4, "me_var"=-e5))
}

elbo_cs_IW <- function(X, inverts, params){
  N <- params$N
  D <- params$D
  T0 <- params$T0
  Mu0 <- params$prior_mean_eta
  P <- params$P
  RP <- Rfast::colsums(P)
  nu0 <- params$prior_df_cs_cov
  V0 <- params$prior_scale_cs_cov
  nu1 <- params$post_df_cs_cov
  V1 <- params$post_scale_cs_cov
  k0 <- params$scaling_cov_eta

  V1_inv <- array(apply(V1, 3, function(x){spdinv(x)}), dim = dim(V1))
  #expectation of inverse of data covariance matrix
  inv_C0 <- sweep(V1_inv, 3, nu1, "*")
  E_log_C0_1 <- apply(nu1, 2, function(x){sum(digamma(0.5*(x + 1 - c(1:D))))})
  E_log_C0_2 <- apply(V1_inv, 3, function(x){D*log(2) +
      determinant(x, logarithm = TRUE)$modulus})
  #expectation of log-determinant of inverse of data covariance matrix
  E_log_C0 <- matrix((E_log_C0_1 + E_log_C0_2), nrow = 1, ncol = T0)
  #covariance parameter of eta's
  L2 <- sweep(V1, 3, nu1*(1/k0 + RP), "/")

  #the eta_i's and C0_i's
  e20 <- rep(0, T0)
  e21 <- rep(0, T0)
  for (i in 1:T0){
    e20[i] <- -(1/(2*k0))*Rfast::mat.mult(L1[i,,drop=FALSE],
                                          Rfast::Tcrossprod(inv_C0[,,i],
                                                           L1[i,,drop=FALSE]))
    e21[i] <- (1/k0)*Rfast::mat.mult(L1[i,,drop=FALSE],
                                     Rfast::Tcrossprod(inv_C0[,,i], Mu0))
  }
  e22 <- apply(inv_C0, 3, function(x){-(1/(2*k0))*
                                      Rfast::mat.mult(Mu0,
                                                      Rfast::Tcrossprod(x, Mu0))})
  e23 <- apply(inv_C0, 3, function(x){-0.5*sum(t(V0)*x)})
  e2 <- -0.5*D*T0*log(2*pi) + 0.5*(nu0 + D + 2)*sum(E_log_C0) + sum(e20) -
    (1/(2*k0))*D*sum(1/(1/k0 + RP)) + sum(e21) + sum(e22) - 0.5*T0*nu0*D*log(2) -
    T0*(0.25*D*(D-1)*log(pi) + sum(lgamma(0.5*(nu0 + 1 - c(1:D))))) +
    0.5*T0*nu0*determinant(V0, logarithm = TRUE)$modulus + sum(e23)

  #the data x_n's
  e30 <- matrix(0, nrow = N, ncol = T0)
  e31 <- matrix(0, nrow = N, ncol = T0)
  for (n in 1:N){
    for (i in 1:T0){
      e30[n,i] <- -0.5*P[n,i]*Rfast::mat.mult(X[n,,drop=FALSE],
                                              Rfast::Tcrossprod(inv_C0[,,i],
                                                                X[n,,drop=FALSE]))
      e31[n,i] <- P[n,i]*Rfast::mat.mult(L1[i,,drop=FALSE],
                                         Rfast::Tcrossprod(inv_C0[,,i],
                                                           X[n,,drop=FALSE]))
    }
  }
  e3 <- -N*0.5*D*log(2*pi) + sum(RP*0.5*E_log_C0) + sum(e30) + sum(e31) +
    sum(RP*e20) - 0.5*D*sum(RP/(1/k0 + RP))

  #the variationa distributions
  e420 <- apply(L2, 3, function(x){-0.5*determinant(x, logarithm = TRUE)$modulus})
  e421 <- apply(nu1, 2, function(x){-0.25*D*(D-1)*log(pi) -
      sum(lgamma(0.5*(x + 1 - c(1:D))))})
  e422 <- apply(V1, 3, function(x){determinant(x, logarithm = TRUE)$modulus})
  e4 <- -0.5*D*T0*log(2*pi) + sum(e420) - 0.5*D*T0 - 0.5*D*log(2)*sum(nu1) +
    sum(e421) + 0.5*sum(nu1*e422) + 0.5*sum((nu1 + D + 1)*E_log_C0) -
    0.5*D*sum(nu1)

  return(c("e_mean_cov"=e2, "e_data"=e3, "me_var"=-e4))
}

elbo_cs_sparse <- function(X, inverts, params){
  N <- params$N
  D <- params$D
  T0 <- params$T0
  Mu0 <- params$prior_mean_eta
  P <- params$P
  RP <- Rfast::colsums(P)
  a0 <- params$prior_shape_d_cs_cov
  b0 <- params$prior_rate_d_cs_cov
  c0 <- params$prior_var_offd_cs_cov
  a1 <- params$post_shape_d_cs_cov
  B1 <- params$post_rate_d_cs_cov
  C1 <- params$post_var_offd_cs_cov
  k0 <- params$scaling_cov_eta

  #expectation of inverse of C0, data covariance matrix
  inv_C0 <- array(0, c(D, D, T0))
  #inverse of expectation of inverse of C0, data covariance matrix
  E_C0_inv_inv <- array(0, c(D, D, T0))
  for (i in 1:T0){
    inv_C0[,,i] <- Rfast::Diag.matrix(D, a1[1,i]/B1[i,])
    E_C0_inv_inv[,,i] <- Rfast::Diag.matrix(D, B1[i,]/a1[1,i])
  }
  #covariance parameter of eta's
  L2 <- sweep(E_C0_inv_inv, 3, (1/k0 + RP), "/")

  #the eta_i's and C0_i's
  e200 <- rep(0, T0)
  e201 <- e200
  e202 <- e200
  for (i in 1:T0){
    e200[i] <- 0.5*sum(digamma(a1[1,i]) - log(B1[i,]))
    e201[i] <- -(1/(2*k0))*Rfast::mat.mult(L1[i,,drop=FALSE],
                                           Rfast::Tcrossprod(inv_C0[,,i],
                                                             L1[i,,drop=FALSE]))
    e202[i] <- (1/k0)*Rfast::mat.mult(L1[i,,drop=FALSE],
                                      Rfast::Tcrossprod(inv_C0[,,i], Mu0))
  }
  e203 <- apply(inv_C0, 3, function(x){-(1/(2*k0))*
                                       Rfast::mat.mult(Mu0,
                                                       Rfast::Tcrossprod(x, Mu0))})
  e20 <- -0.5*T0*log(2*pi) - 0.5*log(k0) + sum(e200) + sum(e201) -
    (0.5/k0)*D/sum(1/k0 + RP) + sum(e202) + sum(e203)

  e210 <- rep(0, T0)
  for (i in 1:T0){
    e210[i] <- sum(a0*log(b0) - lgamma(a0) +
                     (a0 - 1)*(digamma(a1[1,i]) - log(B1[i,])) - b0*a1[1,i]/B1[i,])
  }
  e21 <- sum(e210) + sum(-log(2*c0) - C1[!diag(D)]/c0)
  e2 <- e20 + e21

  #the data X
  e30 <- matrix(0, nrow = N, ncol = T0)
  e31 <- matrix(0, nrow = N, ncol = T0)
  for (n in 1:N){
    for (i in 1:T0){
      e30[n,i] <- -0.5*P[n,i]*Rfast::mat.mult(X[n,,drop=FALSE],
                                              Rfast::Tcrossprod(inv_C0[,,i],
                                                                X[n,,drop=FALSE]))
      e31[n,i] <- P[n,i]*Rfast::mat.mult(L1[i,,drop=FALSE],
                                         Rfast::Tcrossprod(inv_C0[,,i],
                                                           X[n,,drop=FALSE]))
    }
  }
  e3 <- -N*0.5*D*log(2*pi) + sum(RP*e200) + sum(e30) + sum(e31) +
    sum(RP*k0*e201) - 0.5*D*sum(RP/(1/k0 + RP))

  #the variationa distributions
  e420 <- rep(0, T0)
  for (i in 1:T0){
    e420[i] <- sum(a1[1,i]*log(B1[i,]) - lgamma(a1[1,i]) +
                     (a1[1,i] -1)*(digamma(a1[1,i]) - log(B1[i,])) - a1[1,i])
  }
  e4 <- sum(e420) + sum(-log(2*C1[!diag(D)]) - 1)

  return(c("e_mean_cov"=e2, "e_data"=e3, "me_var"=-e4))
}

elbo_cs_offd_normal <- function(X, inverts, params){
  N <- params$N
  D <- params$D
  T0 <- params$T0
  Mu0 <- params$prior_mean_eta
  P <- params$P
  RP <- Rfast::colsums(P)
  a0 <- params$prior_shape_d_cs_cov
  b0 <- params$prior_rate_d_cs_cov
  a1 <- params$post_shape_d_cs_cov
  B1 <- params$post_rate_d_cs_cov
  C1 <- params$post_mean_offd_cs_cov
  k0 <- params$scaling_cov_eta

  #expectation of inverse of data covariance matrix
  inv_C0 <- array(0, c(D, D, T0))
  #inverse of expectation of inverse of data covariance matrix
  E_C0_inv_inv <- array(0, c(D, D, T0))
  for (i in 1:T0){
    inv_C0[,,i] <- Rfast::Diag.fill(C1[,,i], a1[1,i]/B1[i,])
    E_C0_inv_inv[,,i] <- Rfast::spdinv(inv_C0[,,i])
  }
  #covariance parameter of eta's
  L2 <- sweep(E_C0_inv_inv, 3, (1/k0 + RP), "/")

  #the eta_i's and C0_i's
  e200 <- rep(0, T0)
  e201 <- e200
  e202 <- e200
  for (i in 1:T0){
    e200[i] <- 0.5*sum(digamma(a1[1,i]) - log(B1[i,]))
    e201[i] <- -(1/(2*k0))*Rfast::mat.mult(L1[i,,drop=FALSE],
                                           Rfast::Tcrossprod(inv_C0[,,i],
                                                             L1[i,,drop=FALSE]))
    e202[i] <- (1/k0)*Rfast::mat.mult(L1[i,,drop=FALSE],
                                      Rfast::Tcrossprod(inv_C0[,,i], Mu0))
  }
  e203 <- apply(inv_C0, 3, function(x){-(1/(2*k0))*
                                       Rfast::mat.mult(Mu0,
                                                       Rfast::Tcrossprod(x, Mu0))})
  e20 <- -0.5*T0*log(2*pi) - 0.5*log(k0) + sum(e200) + sum(e201) -
    (0.5/k0)*D/sum(1/k0 + RP) + sum(e202) + sum(e203)

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
      e30[n,i] <- -0.5*P[n,i]*Rfast::mat.mult(X[n,,drop=FALSE],
                                              Rfast::Tcrossprod(inv_C0[,,i],
                                                                X[n,,drop=FALSE]))
      e31[n,i] <- P[n,i]*Rfast::mat.mult(L1[i,,drop=FALSE],
                                         Rfast::Tcrossprod(inv_C0[,,i],
                                                           X[n,,drop=FALSE]))
    }
  }
  e3 <- -N*0.5*D*log(2*pi) + sum(RP*e200) + sum(e30) + sum(e31) +
    sum(RP*k0*e201) - 0.5*D*sum(RP/(1/k0 + RP))

  #the variationa distributions
  e420 <- rep(0, T0)
  for (i in 1:T0){
    e420[i] <- sum(a1[1,i]*log(B1[i,]) - lgamma(a1[1,i]) +
                     (a1[1,i] -1)*(digamma(a1[1,i]) - log(B1[i,])) - a1[1,i])
  }
  e4 <- sum(e420) + (D^2 - D)*(-0.5*log(2*pi) - 0.5)

  return(c("e_mean_cov"=e2, "e_data"=e3, "me_var"=-e4))
}
