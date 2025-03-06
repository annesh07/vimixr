CVI_eta_covfunction <- function(fixed_variance = FALSE,
                                covariance_type = "diagonal",
                                cluster_specific_covariance = TRUE,
                                variance_prior_type = c("IW", "decomposed",
                                                        "sparse",
                                                        "off-diagonal normal"),
                                X, inverts, params){
  L1 <- params$post_mean_eta     #posterior mean parameter for DP mean parameters
  P <- params$P                  #allocation probability matrix
  D <- ncol(X)                   #dimension of the data
  N <- nrow(X)                   #samples of the data

  RP <- Rfast::colsums(P)
  CP <- Rfast::rowsums(P)

  if(covariance_type == "diagonal") {

    if(fixed_variance) {
      L20 <- params$prior_precision_scalar_eta
      L2 <- params$post_precision_scalar_eta

      inv_C0 <- inverts[["inv_C0"]] #inverse of data covariance
      inv_C00 <- inverts[["inv_C00"]] #inverse of covariance of DP mean parameters
      Mu00 <- Mu0%*%inv_C00

      #updated parameters of eta's
      for (i in 1:T0){
        L1[i,] <- Mu00 + crossprod(P[, i, drop=FALSE], X) %*% inv_C0
        L2[i, 1] <- L20 + sum(P[, i])
      }
      params$post_mean_eta <- L1
      params$post_precision_scalar_eta <- L2
      out <- params

    } else {

      b1 <- params$prior_shape_scalar_cov
      b2 <- params$prior_rate_scalar_cov
      G1 <- params$post_shape_scalar_cov
      G2 <- params$post_rate_scalar_cov
      L2 <- params$post_precision_scalar_eta
      L20 <- params$prior_precision_scalar_eta

      inv_C00 <- inverts[["inv_C00"]] #inverse of covariance of DP mean parameters
      Mu00 <- Mu0%*%inv_C00

      for (i in 1:T0){
        L1[i,] <- Mu00 + (G1/G2)*crossprod(P[, i, drop=FALSE], X)
        L2[i, 1] <- L20 + (G1/G2)*sum(P[, i])
      }
      L21 <- sweep(L1, 1, L2, "/")

      #updated parameters of the scalar multiple of the data covariance matrix
      G1 <- b1 + 0.5*D*sum(P)
      G20 <- sweep(P, 1, 0.5*diag(tcrossprod(X)), "*")
      G21 <- P*t(- tcrossprod(L21, X))
      G22 <- sweep(P, 2, 0.5*diag(tcrossprod(L21)), "*")
      G23 <- sweep(P, 2, 0.5*D/L2, "*")
      G2 <- b2 + sum(G20) + sum(G21) + sum(G22) + sum(G23)

      params$post_mean_eta <- L1
      params$post_precision_scalar_eta <- L2
      params$post_shape_scalar_cov <- G1
      params$post_rate_scalar_cov <- G2
      out <- params

    }

  } else if(covariance_type == "full") {

    if(fixed_variance) {
      L2 <- params$post_cov_eta

      inv_C0 <- inverts[["inv_C0"]]      #inverse of C0
      inv_C00 <- inverts[["inv_C00"]] #inverse of covariance of DP mean parameters
      Mu00 <- Mu0%*%inv_C00

      #updated parameters of eta's
      for (i in 1:T0){
        L1[i,] <- Mu00 + crossprod(P[, i, drop=FALSE], X) %*% inv_C0
        L2[,,i] <- Rfast::spdinv(inv_C00 + sum(P[, i])*(inv_C0))
      }

      params$post_mean_eta <- L1
      params$post_cov_eta <- L2
      out <- params

    } else {
      if(!cluster_specific_covariance) {
        if(variance_prior_type == "IW"){

          nu0 <- params$prior_df_cov
          V0 <- params$prior_scale_cov
          nu <- params$post_df_cov
          V <- params$post_scale_cov
          L2 <- params$post_cov_eta

          inv_C0 <- nu*V  #expected inverse of C0; covariance matrix of data
          inv_V0 <- inverts[["inv_V0"]] #inverse of prior scale matrix of C0_1
          inv_C00 <- inverts[["inv_C00"]] #inverse of covariance of DP mean parameters
          Mu00 <- Mu0%*%inv_C00

          #updated parameters of eta's
          L21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            L1[i,] <- Mu00 + crossprod(P[, i, drop=FALSE], X) %*% inv_C0
            L2[,,i] <- Rfast::spdinv(inv_C00 + sum(P[, i])*inv_C0)
            L21[i,] = L1[i,, drop = FALSE] %*% L2[,,i]
          }

          #updated parameters of C0
          nu <- nu0 + sum(P)
          V1 <- inv_V0
          for (n in 1:N){
            V1 <- V1 + CP[n]*crossprod(X[n,, drop = FALSE])
          }
          for (n in 1:N){
            for (i in 1:T0){
              V1 <- V1 - P[n,i]*2*crossprod(X[n,, drop = FALSE],
                                            L21[i,, drop = FALSE])
            }
          }
          for (i in 1:T0){
            V1 <- V1 + RP[i]*(crossprod(L21[i,, drop = FALSE])) +
                                  Rfast::spdinv(L2[,,i])
          }
          V <- Rfast::spdinv(V1)

          params$post_mean_eta <- L1
          params$post_df_cov <- nu
          params$post_scale_cov <- V
          params$post_cov_eta <- L2
          out <- params

        } else if (variance_prior_type == "decomposed"){

          a0 <- params$prior_scale_diag_decomp
          b0 <- params$prior_rate_diag_decomp
          mu0 <- params$prior_mean_offdiag_decomp
          c0 <- params$prior_var_offdiag_decomp
          a1 <- params$post_scale_diag_decomp
          b1 <- params$post_rate_diag_decomp
          mu1 <- params$post_mean_offdiag_decomp
          c1 <- params$post_var_offdiag_decomp
          L2 <- params$post_cov_eta

          inv_C00 <- inverts[["inv_C00"]] #inverse of covariance of DP mean parameters
          Mu00 <- Mu0%*%inv_C00

          L21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            L21[i,] = L1[i,, drop = FALSE] %*% L2[,,i]
          }
          .
          #updates of C0
          a1 <- rep(a0, D) + 0.5*sum(P)

          b21 <- matrix(0, nrow = T0, ncol = D)
          colPf0 <- colsums(P)
          for (i in 1:T0){
            b20 <- eachrow(X, L21[i,], oper = "-")
            b21[i,] <- colsums(P[,i]*(b20^2)) + colPf0[i]*diag(L2[,,i])
          }
          b21 <- colsums(b21)

          b1 <- (b0 + 0.5*b21)

          c10 <- 1/(1/c0 + b21)
          c10 <- c10[-1]
          c1 <- rep(c10, times = seq_along(c10))
          sigma_lower <- matrix(0, nrow = D, ncol = D)
          sigma_lower[lower.tri(sigma_lower, diag = FALSE)] <- c1

          diag_L <- sqrt(1/b1)*sqrt(pi)/beta(a1,0.5)
          lowerL <- diag(diag_L)
          for (k in 2:D){
            mu10 <- eachcol.apply(X[, 1:(k-1), drop=FALSE],
                                  rowsums(P)*X[,k, drop=FALSE], oper = "*")

            mu20 <- rep(0, (k-1))
            for (n in 1:N){
              for (i in 1:T0){
                mu20 <- mu20 +
                  P[n,i]*(L21[i, 1:(k-1), drop=FALSE]*X[n,k] +
                            X[n, 1:(k-1), drop=FALSE]*L21[i,k])
              }
            }
            mu30 <- rep(0, (k-1))
            for (i in 1:T0){
              mu30 <- mu30 +
                RP[i]*(L2[,,i][k, 1:(k-1)] + L21[i, 1:(k-1), drop=FALSE]*L21[i, k])
            }

            lower_L0 <- lowerL[1:(k-1), 1:(k-1), drop = FALSE]
            muf0 <- eachcol.apply(lower_L0, (mu10 - mu20 + mu30), oper = "*")
            muf <- (mu0/c0 - muf0)/sigma_lower[k, 1:(k-1)]

            lowerL[k,] <- c(muf, diag_L[k], rep(0, (D - (length(muf)+1))))
          }
          mu1 <- lowerL[lower.tri(lowerL, diag = FALSE)]

          a1 <- matrix(a1, nrow = 1)
          b1 <- matrix(b1, nrow = 1)
          mu1 <- matrix(mu1, nrow = 1)
          c1 <- matrix(c1, nrow = 1)

          mean_lower <- matrix(0, nrow = D, ncol = D)
          mean_lower[lower.tri(mean_lower, diag = FALSE)] <- mu1
          sigma_lower <- matrix(0, nrow = D, ncol = D)
          sigma_lower[lower.tri(sigma_lower, diag = FALSE)] <- c1

          mean_L <- mean_lower + diag(as.vector(sqrt(1/b1)*sqrt(pi)/beta(a1,0.5)))
          diag(sigma_lower) <- (1/b1)*(a1 - (sqrt(pi)/beta(a1,0.5))^2)
          inv_C0 <- tcrossprod(mean_L) + diag(rowsums(sigma_lower))


          #update for eta_i's
          for (i in 1:T0){
            L1[i,] <- Mu00 + crossprod(P[, i, drop=FALSE], X) %*% inv_C0
            L2[,,i] <- Rfast::spdinv(inv_C00 + sum(P[, i])*(inv_C0))
          }

          params$post_mean_eta <- L1
          params$post_cov_eta <- L2
          params$post_scale_diag_decomp <- a1
          params$post_rate_diag_decomp <- b1
          params$post_mean_offdiag_decomp <- mu1
          params$post_var_offdiag_decomp <- c1
          out <- params

        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is FALSE")
        }

      }else{
        if(variance_prior_type == "IW"){

          nu0 <- params$prior_df_cs_cov
          V0 <- params$prior_scale_cs_cov
          nu1 <- params$post_df_cs_cov
          V1 <- params$post_scale_cs_cov
          k0 <- params$scaling_cov_eta

          nu1 <- nu0 + matrix(RP, nrow=1)
          for (i in 1:T0){
            V1[,,i] <- matrix(0, D, D)
            for (n in 1:N){
              V1[,,i] <- V1[,,i] + P[n,i]*crossprod(X[n,,drop=FALSE])
            }
            V1[,,i] <- V0 + (1/k0)*crossprod(Mu0) + V1[,,i] + diag(1e-6, D)
            L1[i,] <- (Mu0/k0 + colsums(sweep(X, 1, P[,i], "*")))/(1/k0 + RP[i])
          }

          params$post_df_cs_cov <- nu1
          params$post_scale_cs_cov <- V1
          params$post_mean_eta <- L1
          out <- params

        } else if (variance_prior_type == "sparse"){

          a0 <- params$prior_scale_d_cs_cov
          b0 <- params$prior_rate_d_cs_cov
          c0 <- params$prior_var_offd_cs_cov
          a1 <- params$post_scale_d_cs_cov
          B1 <- params$post_rate_d_cs_cov
          C1 <- params$post_var_offd_cs_cov
          k0 <- params$scaling_cov_eta

          a1 <- matrix(a0 + RP, nrow = 1, ncol = T0)
          for (i in 1:T0){
            B1[i,] <- b0 + Rfast::colsums(sweep(X^2, 1, P[,i], "*"))
            C01 <- 1/c0 + 0.5*abs(crossprod(sweep(X, 1, P[,i], "*"), X))
            C1[,,i] <- Diag.fill(1/C01, rep(0, D))
            L1[i,] <- (Mu0/k0 +
                         Rfast::colsums(sweep(X, 1, P[,i], "*")))/(1/k0 + RP[i])
          }

          params$post_scale_d_cs_cov <- a1
          params$post_rate_d_cs_cov <- B1
          params$post_var_offd_cs_cov <- C1
          params$post_mean_eta <- L1
          out <- params

        } else if (variance_prior_type == "off-diagonal normal"){

          a0 <- params$prior_scale_d_cs_cov
          b0 <- params$prior_rate_d_cs_cov
          a1 <- params$post_scale_d_cs_cov
          B1 <- params$post_rate_d_cs_cov
          C1 <- params$post_mean_offd_cs_cov
          k0 <- params$scaling_cov_eta

          a1 <- matrix(a0 + RP, nrow = 1, ncol = T0)
          for (i in 1:T0){
            B1[i,] <- b0 + Rfast::colsums(sweep(X^2, 1, P[,i], "*"))
            C01 <- 0.5*(crossprod(sweep(X, 1, P[,i], "*"), X))
            C1[,,i] <- Diag.fill(C01, rep(0, D))
            L1[i,] <- (Mu0/k0 +
                         Rfast::colsums(sweep(X, 1, P[,i], "*")))/(1/k0 + RP[i])
          }

          params$post_scale_d_cs_cov <- a1
          params$post_rate_d_cs_cov <- B1
          params$post_mean_offd_cs_cov <- C1
          params$post_mean_eta <- L1
          out <- params

        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is TRUE")
        }
      }

    } else {
      stop("covariance_type can only be either 'diagonal' or 'full'.")
    }
  }

  return(out)
}
