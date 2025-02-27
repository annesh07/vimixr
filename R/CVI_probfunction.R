CVI_probfunction <- function(Plog, type, params){
  P <- exp(Plog)
  RP <- colsums(P)

  #probability matrix update based on latent allocations
  P2 <- matrix(0, N, T0)
  for (n in 1:N){
    #update of the n^th vector is done by considering all the vecors except
    #the n^th one
    P1 <- P[-n,]

    p_eni <- colSums(P1)
    p_vni <- colSums(P1*(1-P1))
    p_enj <- rowSums(apply(P1, 1, f0))
    p_vnj <- rowSums(apply(P1, 1, f1))

    P20 <- log(1 + p_eni) - p_vni/((1 + p_eni)^2) - log(1 + p_eni + p_enj +
                                                          (W1/W2)) +
      (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)

    P21 <- log((W1/W2) + p_enj) - (p_vnj + (W1/(W2^2)))/(((W1/W2) + p_enj)^2) -
      log(1 + p_eni + p_enj + (W1/W2)) +
      (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)
    P22 <- c(0, cumsum(P21)[1:(T0-1)])

    P2[n,] <- P20 + P22
  }

  if (type == "fixed-scalar"){
    C0 <- params$cov_data
    L2 <- params$post_precision_scalar_eta
    L20 <- params$prior_precision_scalar_eta

    inv_C0 <- spdinv(C0) #inverse of C0

    L21 <- sweep(L1, 1, L2, "/")
    P230 <- L21 %*% inv_C0 %*% t(X)
    P231 <- diag(-0.5*L21 %*% inv_C0 %*% t(L21))
    P232 <- - 0.5*sum(diag(inv_C0))/L2
    P233 <- diag(- 0.5*X %*% inv_C0 %*% t(X))
    P_const <- -0.5*(D*log(2*pi) + determinant(C0, logarithm = TRUE)$modulus)

    Plog <- P2 + t(P230) + P_const
    Plog <- sweep(Plog, 2, P231+P232, "+")
    Plog <- sweep(Plog, 1, P233, "+")
    #log-sum-exp trick
    Plog <- t(apply(Plog, 1, function(x){x - max(x) - log(sum(exp(x - max(x))))}))

    return(c("Plog"=Plog))

  } else if (type == "fixed-matrix"){
    C0 <- params$cov_data
    L2 <- params$post_cov_eta
    C00 <- params$prior_cov_eta

    inv_C0 <- spdinv(C0)      #inverse of C0

    L21 <- matrix(0, nrow = T0, ncol = D)
    for (i in 1:T0){
      L21[i,] = L1[i,, drop = FALSE] %*% L2[,,i]
    }

    P230 <- L21 %*% inv_C0 %*% t(X)
    P231 <- apply(L21, 1, function(x){-0.5*t(x)%*%inv_C0%*%x})
    P232 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C0)*x)})
    P233 <- apply(X, 1, function(x){-0.5*t(x)%*%inv_C0%*%x})
    P_const <- -0.5*(D*log(2*pi) - determinant(inv_C0, logarithm = TRUE)$modulus)

    return(c("P230"=P230, "P231"=P231, "P232"=P232, "P233"=P233,
             "P_const"=P_const))

  } else if (type == "varied-scalar"){
    b1 <- params$prior_shape_scalar_cov
    b2 <- params$prior_rate_scalar_cov
    G1 <- params$post_shape_scalar_cov
    G2 <- params$post_rate_scalar_cov
    L2 <- params$post_precision_scalar_eta
    L20 <- params$prior_precision_scalar_eta

    L21 <- sweep(L1, 1, L2, "/")
    P230 <- (G1/G2)*L21 %*% t(X)
    P231 <- diag(-0.5*(G1/G2)*L21 %*% t(L21))
    P232 <- - 0.5*D*(G1/G2)/L2
    P233 <- (G1/G2)*diag(X %*% t(X))
    P_const <- - 0.5*(D*log(2*pi) - D*(digamma(G1) - log(G2)))

    return(c("P230"=P230, "P231"=P231, "P232"=P232, "P233"=P233,
             "P_const"=P_const))

  } else if (type == "varied-matrix-IW"){
    nu0 <- params$prior_df_cov
    V0 <- params$prior_scale_cov
    nu <- params$post_df_cov
    V <- params$post_scale_cov
    L2 <- params$post_cov_eta
    C00 <- params$prior_cov_eta

    inv_C0 <- nu*V           #expected inverse of C0; covariance matrix of data

    L21 <- matrix(0, nrow = T0, ncol = D)
    for (i in 1:T0){
      L21[i,] = L1[i,, drop = FALSE] %*% L2[,,i]
    }
    P230 <- L21 %*% inv_C0 %*% t(X)
    P231 <- diag(-0.5*L21 %*% inv_C0 %*% t(L21))
    P232 <- apply(L2, 3, function(x){-0.5*sum(diag(inv_C0 %*% x))})
    P233 <- diag(X %*% inv_C0 %*% t(X))
    P_const <- -0.5*(D*log(2*pi) - (sum(digamma(0.5*(nu + 1 - c(1:D)))) +
                                      D*log(2) +
                                      determinant(V, logarithm = TRUE)$modulus))

    return(c("P230"=P230, "P231"=P231, "P232"=P232, "P233"=P233,
             "P_const"=P_const))

  } else if (type == "varied-matrix-decomposed"){
    a0 <- params$prior_scale_diag_decomp
    b0 <- params$prior_rate_diag_decomp
    mu0 <- params$prior_mean_offdiag_decomp
    c0 <- params$prior_var_offdiag_decomp
    a1 <- params$post_scale_diag_decomp
    b1 <- params$post_rate_diag_decomp
    mu1 <- params$post_mean_offdiag_decomp
    c1 <- params$post_var_offdiag_decomp
    L2 <- params$post_cov_eta
    C00 <- params$prior_cov_eta

    mean_lower <- matrix(0, nrow = D, ncol = D) #mean matrix of the decomposed
    mean_lower[lower.tri(mean_lower, diag = FALSE)] <- mu1
    sigma_lower <- matrix(0, nrow = D, ncol = D) #var matrix of the decomposed
    sigma_lower[lower.tri(sigma_lower, diag = FALSE)] <- c1
    mean_L <- mean_lower + diag(sqrt(1/b1)*sqrt(pi)/beta(a1,0.5))
    diag(sigma_lower) <- (1/b1)*(a1 - (sqrt(pi)/beta(a1,0.5))^2)
    #expected inverse of C0; covariance matrix of data
    inv_C0 <- mean_L %*% t(mean_L) + diag(rowsums(sigma_lower))

    L21 <- matrix(0, nrow = T0, ncol = D)
    for (i in 1:T0){
      L21[i,] = L1[i,, drop = FALSE] %*% L2[,,i]
    }
    P230 <- L21 %*% inv_C0 %*% t(X)
    P231 <- diag(-0.5*L21 %*% inv_C0 %*% t(L21))
    P232 <- apply(L2, 3, function(x){-0.5*sum(diag(inv_C0 %*% x))})
    P233 <- diag(X %*% inv_C0 %*% t(X))
    P_const <- -0.5*(D*log(2*pi) - sum(digamma(a1) - log(b1)))

    return(c("P230"=P230, "P231"=P231, "P232"=P232, "P233"=P233,
             "P_const"=P_const))

  } else if (type == "cluster-specific-IW"){
    nu0 <- params$prior_df_cs_cov
    V0 <- params$prior_scale_cs_cov
    nu1 <- params$post_df_cs_cov
    V1 <- params$post_scale_cs_cov
    k0 <- params$scaling_cov_eta

    V1_inv <- array(apply(V1, 3, function(x){spdinv(x)}), dim = dim(V1))
    #expectation of inverse of data covariance matrix
    E_C0 <- sweep(V1_inv, 3, nu1, "*")
    E_log_C0_1 <- apply(nu1, 2, function(x){sum(digamma(0.5*(x + 1 - c(1:D))))})
    E_log_C0_2 <- apply(V1_inv, 3, function(x){D*log(2) +
        determinant(x, logarithm = TRUE)$modulus})
    #expectation of log-determinant of inverse of data covariance matrix
    E_log_C0 <- matrix((E_log_C0_1 + E_log_C0_2), nrow = 1, ncol = T0)

    #updating the latent probability values
    P230 <- matrix(0, nrow = N, ncol = T0)
    P231 <- matrix(0, nrow = N, ncol = T0)
    for (n in 1:N){
      for (i in 1:T0){
        P230[n,i] <- -0.5*X[n,,drop=FALSE] %*% E_C0[,,i] %*% t(X[n,,drop=FALSE])
        P231[n,i] <- L1[i,,drop=FALSE] %*% E_C0[,,i] %*% t(X[n,,drop=FALSE])
      }
    }
    P232 <- matrix(0, nrow = 1, ncol = T0)
    for (i in 1:T0){
      P232[1,i] <- -0.5*L1[i,,drop=FALSE] %*% E_C0[,,i] %*% t(L1[i,,drop=FALSE])
    }
    P233 <- 0.5*E_log_C0
    P_const <- -0.5*D*(log(2*pi) + 1/(1/k0 + RP))

    return(c("P230"=P230, "P231"=P231, "P232"=P232, "P233"=P233,
             "P_const"=P_const))

  } else if (type == "cluster-specific-sparse"){
    a0 <- params$prior_scale_d_cs_cov
    b0 <- params$prior_rate_d_cs_cov
    c0 <- params$prior_var_offd_cs_cov
    a1 <- params$post_scale_d_cs_cov
    B1 <- params$post_rate_d_cs_cov
    C1 <- params$post_var_offd_cs_cov
    k0 <- params$scaling_cov_eta

    #expectation of inverse of C0, data covariance matrix
    E_C0_inv <- array(0, c(D, D, T0))
    for (i in 1:T0){
      E_C0_inv[,,i] <- Diag.matrix(D, a1[1,i]/B1[i,])
    }

    #updating the latent probability values
    P230 <- matrix(0, nrow = N, ncol = T0)
    P231 <- matrix(0, nrow = N, ncol = T0)
    for (n in 1:N){
      for (i in 1:T0){
        P230[n,i] <- -0.5*X[n,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(X[n,,drop=FALSE])
        P231[n,i] <- L1[i,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(X[n,,drop=FALSE])
      }
    }
    P232 <- matrix(0, nrow = 1, ncol = T0)
    P233 <- P232
    for (i in 1:T0){
      P232[1,i] <- -0.5*L1[i,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(L1[i,,drop=FALSE])
      P233[1,i] <- 0.5*sum(digamma(a1[1,i]) - log(B1[i,]))
    }
    P_const <- -0.5*D*(log(2*pi) + 1/(1/k0 + RP))

    return(c("P230"=P230, "P231"=P231, "P232"=P232, "P233"=P233,
             "P_const"=P_const))

  } else if (type == "cluster-specific-offdNormal"){
    a0 <- params$prior_scale_d_cs_cov
    b0 <- params$prior_rate_d_cs_cov
    a1 <- params$post_scale_d_cs_cov
    B1 <- params$post_rate_d_cs_cov
    C1 <- params$post_mean_offd_cs_cov
    k0 <- params$scaling_cov_eta

    #expectation of inverse of data covariance matrix
    E_C0_inv <- array(0, c(D, D, T0))
    for (i in 1:T0){
      E_C0_inv[,,i] <- Diag.fill(C1[,,i], a1[1,i]/B1[i,])
    }

    #updating the latent probability values
    P230 <- matrix(0, nrow = N, ncol = T0)
    P231 <- matrix(0, nrow = N, ncol = T0)
    for (n in 1:N){
      for (i in 1:T0){
        P230[n,i] <- -0.5*X[n,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(X[n,,drop=FALSE])
        P231[n,i] <- L1[i,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(X[n,,drop=FALSE])
      }
    }
    P232 <- matrix(0, nrow = 1, ncol = T0)
    P233 <- P232
    for (i in 1:T0){
      P232[1,i] <- -0.5*L1[i,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(L1[i,,drop=FALSE])
      P233[1,i] <- 0.5*sum(digamma(a1[1,i]) - log(B1[i,]))
    }
    P_const <- -0.5*D*(log(2*pi) + 1/(1/k0 + RP))

    return(c("P230"=P230, "P231"=P231, "P232"=P232, "P233"=P233,
             "P_const"=P_const))
  }
}
