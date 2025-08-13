#' Update of the variational parameters
#'
#' @param fixed_variance whether the covariance is fixed or estimated.
#' Default is \code{FALSE} which means it is estimated.
#' @param covariance_type The assumed type of the covariance matrix.
#' Can be either \code{"diagonal"} if it is the identify multiplied by a scalar,
#' or \code{"full"} for a fully unspecified covariance matrix.
#' @param cluster_specific_covariance whether the the covariance is shared across
#' estimated clusters or is cluster specific. Default is \code{TRUE} which means it is cluster specific.
#' @param variance_prior_type character string specifying the type of prior distribution
#' for the covariance when cluster_specific_covariance is \code{TRUE}.
#' Can be either \code{"IW"} or \code{"decomposed"} if \code{cluster_specific_covariance} is \code{FALSE},
#' and can be either \code{"IW"}, \code{"sparse"} or \code{"off-diagonal normal"} otherwise.
#' @param X the data matrix
#' @param inverts a list of inverses
#' @param params a list of required arguments
#'
#' @return Updated parameters
#'
#' @importFrom Rfast rowsums colsums spdinv eachrow eachcol.apply Diag.fill
#' Diag.matrix
#'
#' @export
#'
CVI_update_function <- function(fixed_variance = FALSE,
                                covariance_type = "diagonal",
                                cluster_specific_covariance = TRUE,
                                variance_prior_type = c("IW", "decomposed",
                                                        "sparse",
                                                        "off-diagonal normal"),
                                X, inverts, params){
  N <- params$N
  D <- params$D
  T0 <- params$T0
  s1 <- params$prior_shape_alpha #shape parameter for alpha prior
  s2 <- params$prior_rate_alpha  #rate parameter for alpha prior
  W1 <- params$post_shape_alpha  #shape parameter for alpha posterior
  W2 <- params$post_rate_alpha   #rate parameter for alpha posterior
  Mu0 <- params$prior_mean_eta   #prior mean for DP mean parameters; vector
  L1 <- params$post_mean_eta     #posterior mean for DP mean parameters; matrix
  P <- params$P                  #allocation probability matrix
  
  RP <- Rfast::colsums(P)
  #probability matrix update based on latent allocations
  P2 <- matrix(0, N, T0)
  for (n in 1:N){
    #update of the n^th vector is done by considering all the vecors except
    #the n^th one
    P1 <- P[-n,]

    p_eni <- Rfast::colsums(P1) #cluster proportions
    p_vni <- Rfast::colsums(P1*(1-P1)) #cluster proportion variance
    p_enj <- cum_clustprop(P1) #cummulative cluster proportion
    p_vnj <- cum_clustprop_var(P1) #cummulative cluster proportion variance



    P20 <- log(1 + p_eni) - p_vni/((1 + p_eni)^2) - log(1 + p_eni + p_enj +
                                                          (W1/W2)) +
      (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)

    P21 <- log((W1/W2) + p_enj) - (p_vnj + (W1/(W2^2)))/(((W1/W2) + p_enj)^2) -
      log(1 + p_eni + p_enj + (W1/W2)) +
      (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)
    P22 <- c(0, cumsum(P21[1:(T0-1)]))

    P2[n,] <- P20 + P22
  }

  if(covariance_type == "diagonal") {

    if(fixed_variance) {
      L20 <- params$prior_precision_scalar_eta
      L2 <- params$post_precision_scalar_eta

      inv_C0 <- inverts[["inv_C0"]] #inverse of C0
      inv_C00 <- inverts[["inv_C00"]] #inverse of covariance of DP mean parameters
      Mu00 <- mat_mult(Mu0, inv_C00)

      L21 <- L1/c(L2)
      P230 <- mat_mult_t(X, inv_C0, L21)
      P231 <- - 0.5*quadratic_form_diag(L21, inv_C0)
      P232 <- - 0.5*sum(diag(inv_C0))/L2
      P233 <- - 0.5*quadratic_form_diag(X, inv_C0)
      P_const <- -0.5*(D*log(2*pi) - determinant(inv_C0, logarithm = TRUE)$modulus)
      #log probability matrix update
      Plog <- P2 + P230 + matrix((P_const + P231 + P232), nrow = N, ncol = T0,
                                 byrow = TRUE) + P233
      #log-sum-exp trick
      Plog <- t(apply(Plog, 1, function(x){
        mx <- max(x)
        x - mx - log(sum(exp(x - mx)))
      }))
      P <- exp(Plog)
      RP <- Rfast::colsums(P)
      ord <- order(RP, decreasing = TRUE)
      Plog <- Plog[,ord]
      P <- P[,ord]
      RP <- Rfast::colsums(P)

      #updated parameters of eta's
      for (i in 1:T0){
        L1[i,] <- Mu00 + t_mat_mult(P[, i, drop=FALSE], X, inv_C0)
        L2[i, 1] <- L20 + sum(P[, i])
      }

      params$post_precision_scalar_eta <- L2
      params$post_mean_eta <- L1
      params$log_prob_matrix <- Plog
      params$P <- P
      out <- params

    } else {

      b1 <- params$prior_shape_scalar_cov
      b2 <- params$prior_rate_scalar_cov
      G1 <- params$post_shape_scalar_cov
      G2 <- params$post_rate_scalar_cov
      L2 <- params$post_precision_scalar_eta
      L20 <- params$prior_precision_scalar_eta

      inv_C00 <- inverts[["inv_C00"]] #inverse of covariance of DP mean parameters
      Mu00 <- mat_mult(Mu0, inv_C00)

      L21 <- L1/c(L2)
      P230 <- (G1/G2)*mat_mult(X, t(L21))
      P231 <- -0.5*(G1/G2)*Rfast::rowsums(L21^2)
      P232 <- - 0.5*D*(G1/G2)/L2
      P233 <- -0.5*(G1/G2)*Rfast::rowsums(X^2)
      P_const <- - 0.5*(D*log(2*pi) - D*(digamma(G1) - log(G2)))
      #log probability matrix update
      Plog <- P2 + P230 + matrix((P_const + P231 + P232), nrow = N, ncol = T0,
                                 byrow = TRUE) + P233
      #log-sum-exp trick
      Plog <- t(apply(Plog, 1, function(x){
        mx <- max(x)
        x - mx - log(sum(exp(x - mx)))
      }))
      P <- exp(Plog)
      RP <- Rfast::colsums(P)
      ord <- order(RP, decreasing = TRUE)
      Plog <- Plog[,ord]
      P <- P[,ord]
      RP <- Rfast::colsums(P)

      #update of parameters of eta's
      for (i in 1:T0){
        L1[i,] <- Mu00 + (G1/G2)*Rfast::colsums(X*P[, i])
        L2[i, 1] <- L20 + (G1/G2)*sum(P[, i])
      }
      L21 <- L1/c(L2)

      #updated parameters of the scalar multiple of the data covariance matrix
      G1 <- b1 + 0.5*D*sum(P)
      G20 <- P * c(0.5*Rfast::rowsums(X^2))
      G21 <- P*(-mat_mult(X, t(L21)))
      G22 <- Rfast::eachrow(P, 0.5*Rfast::rowsums(L21^2), oper = "*")
      G23 <- Rfast::eachrow(P, 0.5*D/L2, oper = "*")
      G2 <- b2 + sum(G20) + sum(G21) + sum(G22) + sum(G23)

      params$post_precision_scalar_eta <- L2
      params$post_shape_scalar_cov <- G1
      params$post_rate_scalar_cov <- G2
      params$post_mean_eta <- L1
      params$log_prob_matrix <- Plog
      params$P <- P
      out <- params

    }

  } else if(covariance_type == "full") {

    if(fixed_variance) {
      L2 <- params$post_cov_eta

      inv_C0 <- inverts[["inv_C0"]]      #inverse of C0
      inv_C00 <- inverts[["inv_C00"]] #inverse of covariance of DP mean parameters
      Mu00 <- mat_mult(Mu0, inv_C00)

      L21 <- matrix(0, nrow = T0, ncol = D)
      for (i in 1:T0){
        L21[i,] = mat_mult(L1[i,, drop = FALSE], L2[,,i])
      }

      P230 <- mat_mult_t(X, inv_C0, L21)
      P231 <- - 0.5*quadratic_form_diag(L21, inv_C0)
      P232 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C0)*x)})
      P233 <- - 0.5*quadratic_form_diag(X, inv_C0)
      P_const <- -0.5*(D*log(2*pi) - determinant(inv_C0, logarithm = TRUE)$modulus)
      #log probability matrix update
      Plog <- P2 + P230 + matrix((P_const + P231 + P232), nrow = N, ncol = T0,
                                 byrow = TRUE) + P233
      #log-sum-exp trick
      Plog <- t(apply(Plog, 1, function(x){
        mx <- max(x)
        x - mx - log(sum(exp(x - mx)))
      }))
      P <- exp(Plog)
      RP <- Rfast::colsums(P)
      ord <- order(RP, decreasing = TRUE)
      Plog <- Plog[,ord]
      P <- P[,ord]
      RP <- Rfast::colsums(P)

      #updated parameters of eta's
      for (i in 1:T0){
        L1[i,] <- Mu00 + t_mat_mult(P[, i, drop=FALSE], X, inv_C0)
        L2[,,i] <- Rfast::spdinv(inv_C00 + sum(P[, i])*(inv_C0))
      }

      params$post_cov_eta <- L2
      params$post_mean_eta <- L1
      params$log_prob_matrix <- Plog
      params$P <- P
      out <- params

    } else {
      if(!cluster_specific_covariance) {
        if(variance_prior_type == "IW"){

          nu0 <- params$prior_df_cov
          V0 <- params$prior_scale_cov
          nu <- params$post_df_cov
          V <- params$post_scale_cov
          L2 <- params$post_cov_eta

          inv_C0 <- nu*V           #expected inverse of C0; covariance matrix of data
          inv_V0 <- inverts[["inv_V0"]] #inverse of prior scale matrix of C0_1
          inv_C00 <- inverts[["inv_C00"]] #inverse of covariance of DP mean parameters
          Mu00 <- mat_mult(Mu0, inv_C00)

          L21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            L21[i,] = mat_mult(L1[i,, drop = FALSE], L2[,,i])
          }
          P230 <- mat_mult_t(X, inv_C0, L21)
          P231 <- - 0.5*quadratic_form_diag(L21, inv_C0)
          P232 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C0)*x)})
          P233 <- - 0.5*quadratic_form_diag(X, inv_C0)
          P_const <- -0.5*(D*log(2*pi) - (sum(digamma(0.5*(nu + 1 - c(1:D)))) +
                                            D*log(2) +
                                            determinant(V, logarithm = TRUE)$modulus))
          #log probability matrix update
          Plog <- P2 + P230 + matrix((P_const + P231 + P232), nrow = N, ncol = T0,
                                     byrow = TRUE) + P233
          #log-sum-exp trick
          Plog <- t(apply(Plog, 1, function(x){
            mx <- max(x)
            x - mx - log(sum(exp(x - mx)))
          }))
          P <- exp(Plog)
          RP <- Rfast::colsums(P)
          ord <- order(RP, decreasing = TRUE)
          Plog <- Plog[,ord]
          P <- P[,ord]
          RP <- Rfast::colsums(P)
          CP <- Rfast::rowsums(P)

          #updated parameters of eta's
          L21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            L1[i,] <- Mu00 + t_mat_mult(P[, i, drop=FALSE], X, inv_C0)
            L2[,,i] <- Rfast::spdinv(inv_C00 + sum(P[, i])*inv_C0)
            L21[i,] = mat_mult(L1[i,, drop = FALSE], L2[,,i])
          }

          #updated parameters of C0
          nu <- nu0 + sum(P)
          V10 <- t_mat_mult(X, diag(CP), X)
          V20 <- 2*t_mat_mult(X, P, L21)
          V30 <- t_mat_mult(L21, diag(RP), L21)
          #missing solve(L2[,,i])
          inv_L2 <- apply(L2, 3, function(x){Rfast::spdinv(x)})
          V40 <- sweep_3D(inv_L2, RP, c(D, D, T0))
          V1 <- inv_V0 + V10 - V20 + V30 + rowSums(V40, dims = 2)
          V <- Rfast::spdinv(V1)

          params$post_df_cov <- nu
          params$post_scale_cov <- V
          params$post_cov_eta <- L2
          params$post_mean_eta <- L1
          params$log_prob_matrix <- Plog
          params$P <- P
          out <- params

        } else if (variance_prior_type == "decomposed"){

          a0 <- params$prior_shape_diag_decomp
          b0 <- params$prior_rate_diag_decomp
          mu0 <- params$prior_mean_offdiag_decomp
          c0 <- params$prior_var_offdiag_decomp
          a1 <- params$post_shape_diag_decomp
          b1 <- params$post_rate_diag_decomp
          mu1 <- params$post_mean_offdiag_decomp
          c1 <- params$post_var_offdiag_decomp
          L2 <- params$post_cov_eta

          inv_C00 <- inverts[["inv_C00"]] #inverse of covariance of DP mean parameters
          Mu00 <- mat_mult(Mu0, inv_C00)

          mean_lower <- matrix(0, nrow = D, ncol = D) #mean matrix of the decomposed
          mean_lower[lower.tri(mean_lower, diag = FALSE)] <- mu1
          sigma_lower <- matrix(0, nrow = D, ncol = D) #var matrix of the decomposed
          sigma_lower[lower.tri(sigma_lower, diag = FALSE)] <- c1
          mean_L <- mean_lower + diag(sqrt(1/b1)*sqrt(pi)/beta(a1,0.5))
          diag(sigma_lower) <- (1/b1)*(a1 - (sqrt(pi)/beta(a1,0.5))^2)
          #expected inverse of C0; covariance matrix of data
          inv_C0 <- mat_mult(mean_L, t(mean_L)) +
            diag(Rfast::rowsums(sigma_lower))

          L21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            L21[i,] = mat_mult(L1[i,, drop = FALSE], L2[,,i])
          }
          P230 <- mat_mult_t(X, inv_C0, L21)
          P231 <- - 0.5*quadratic_form_diag(L21, inv_C0)
          P232 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C0)*x)})
          P233 <- - 0.5*quadratic_form_diag(X, inv_C0)
          P_const <- -0.5*(D*log(2*pi) - sum(digamma(a1) - log(b1)))
          #log probability matrix update
          Plog <- P2 + P230 + matrix((P_const + P231 + P232), nrow = N, ncol = T0,
                                     byrow = TRUE) + P233
          #log-sum-exp trick
          Plog <- t(apply(Plog, 1, function(x){
            mx <- max(x)
            x - mx - log(sum(exp(x - mx)))
          }))
          P <- exp(Plog)
          RP <- Rfast::colsums(P)
          ord <- order(RP, decreasing = TRUE)
          Plog <- Plog[,ord]
          P <- P[,ord]
          RP <- Rfast::colsums(P)

          #update for eta_i's
          L21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            L1[i,] <- Mu00 + t_mat_mult(P[, i, drop=FALSE], X ,inv_C0)
            L2[,,i] <- Rfast::spdinv(inv_C00 + sum(P[, i])*(inv_C0))
            L21[i,] = mat_mult(L1[i,, drop = FALSE], L2[,,i])
          }


          #updates of C0
          a1 <- rep(a0, D) + 0.5*sum(P)

          b21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            b20 <- sweep(X, 2, L21[i,], "-")
            b21[i,] <- Rfast::colsums(P[,i]*(b20^2)) + RP[i]*diag(L2[,,i])
          }
          b21 <- Rfast::colsums(b21)

          b1 <- (b0 + 0.5*b21)

          c10 <- 1/(1/c0 + b21)
          c10 <- c10[-1]
          c1 <- rep(c10, times = seq_along(c10))
          sigma_lower <- matrix(0, nrow = D, ncol = D)
          sigma_lower[lower.tri(sigma_lower, diag = FALSE)] <- c1

          diag_L <- sqrt(1/b1)*sqrt(pi)/beta(a1,0.5)
          lowerL <- diag(diag_L)
          for (k in 2:D){
            mu10 <- Rfast::eachcol.apply(X[, 1:(k-1), drop=FALSE],
                                         Rfast::rowsums(P)*X[,k, drop=FALSE],
                                          oper = "*")

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
            muf0 <- Rfast::eachcol.apply(lower_L0, (mu10 - mu20 + mu30), oper = "*")
            muf <- (mu0/c0 - muf0)/sigma_lower[k, 1:(k-1)]
            lowerL[k,] <- c(muf, diag_L[k], rep(0, (D - (length(muf)+1))))
          }
          mu1 <- lowerL[lower.tri(lowerL, diag = FALSE)]
          a1 <- matrix(a1, nrow = 1)
          b1 <- matrix(b1, nrow = 1)
          mu1 <- matrix(mu1, nrow = 1)
          c1 <- matrix(c1, nrow = 1)

          params$post_cov_eta <- L2
          params$post_shape_diag_decomp <- a1
          params$post_rate_diag_decomp <- b1
          params$post_mean_offdiag_decomp <- mu1
          params$post_var_offdiag_decomp <- c1
          params$post_mean_eta <- L1
          params$log_prob_matrix <- Plog
          params$P <- P
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
          
          RP <- Rfast::colsums(P)
          
          nu1 <- nu0 + matrix(RP, nrow=1)
          for (i in 1:T0){
            V10 <- t_mat_mult(X, diag(P[,i]), X)
            V1[,,i] <- V0 + (1/k0)*matrix(Mu0, nrow = D, ncol = D,
                                          byrow = TRUE)*c(Mu0) +
              V10 + diag(1e-6, D)
            L1[i,] <- (Mu0/k0 + Rfast::eachcol.apply(X, P[,i], oper = "*"))/(1/k0 + RP[i])
          }
          
          V1_inv <- array(apply(V1, 3, function(x){spdinv(x)}), dim = dim(V1))
          #expectation of inverse of data covariance matrix
          inv_C0 <- sweep_3D(V1_inv, nu1, c(D, D, T0))
          E_log_C0_1 <- apply(nu1, 2, function(x){sum(digamma(0.5*(x + 1 - c(1:D))))})
          E_log_C0_2 <- apply(V1_inv, 3, function(x){D*log(2) +
              determinant(x, logarithm = TRUE)$modulus})
          #expectation of log-determinant of inverse of data covariance matrix
          E_log_C0 <- matrix((E_log_C0_1 + E_log_C0_2), nrow = 1, ncol = T0)

          #updating the latent probability values
          P230 <- matrix(0, nrow = N, ncol = T0)
          P233 <- matrix(0, nrow = N, ncol = T0)
          P231 <- matrix(0, nrow = 1, ncol = T0)
          for (i in 1:T0){
            temp <- inv_C0[,,i]
            P233[,i] <- -0.5*quadratic_form_diag(X, temp)
            P230[,i] <- mat_mult_t(L1[i,,drop=FALSE], temp, X)
            P231[1,i] <- -0.5*mat_mult_t(L1[i,,drop=FALSE], temp,
                                         L1[i,,drop=FALSE])
          }
          P232 <- 0.5*E_log_C0
          P_const <- -0.5*D*(log(2*pi) + 1/(1/k0 + RP))
          #log probability matrix update
          Plog <- P2 + P230 + matrix((P_const + P231 + P232), nrow = N, ncol = T0,
                                     byrow = TRUE) + P233
          #log-sum-exp trick
          Plog <- t(apply(Plog, 1, function(x){
            mx <- max(x)
            x - mx - log(sum(exp(x - mx)))
          }))
          P <- exp(Plog)
          RP <- Rfast::colsums(P)
          ord <- order(RP, decreasing = TRUE)
          Plog <- Plog[,ord]
          P <- P[,ord]

          params$post_df_cs_cov <- nu1
          params$post_scale_cs_cov <- V1
          params$post_mean_eta <- L1
          params$log_prob_matrix <- Plog
          params$P <- P
          out <- params

        } else if (variance_prior_type == "sparse"){

          a0 <- params$prior_shape_d_cs_cov
          b0 <- params$prior_rate_d_cs_cov
          c0 <- params$prior_var_offd_cs_cov
          a1 <- params$post_shape_d_cs_cov
          B1 <- params$post_rate_d_cs_cov
          C1 <- params$post_var_offd_cs_cov
          k0 <- params$scaling_cov_eta
          
          RP <- Rfast::colsums(P)
          a1 <- matrix(a0 + RP, nrow = 1, ncol = T0)
          for (i in 1:T0){
            B1[i,] <- b0[i,] + 0.5*Rfast::eachcol.apply(X^2, P[,i], oper = "*") +
              (0.5/k0)*(Mu0^2)
            C01 <- 1/c0 + abs(t_mat_mult(X, diag(P[,i]), X) + 
                                (1/k0)*mat_mult(t(Mu0), Mu0))
            C1[,,i] <- Rfast::Diag.fill(1/C01, rep(0, D))
            L1[i,] <- (Mu0/k0 +
                         Rfast::eachcol.apply(X, P[,i], oper = "*"))/(1/k0 + RP[i])
          }
          #expectation of inverse of C0, data covariance matrix
          inv_C0 <- array(0, c(D, D, T0))

          #updating the latent probability values
          P230 <- matrix(0, nrow = N, ncol = T0)
          P233 <- matrix(0, nrow = N, ncol = T0)
          P231 <- matrix(0, nrow = 1, ncol = T0)
          P232 <- matrix(0, nrow = 1, ncol = T0)
          for (i in 1:T0){
            inv_C0[,,i] <- temp <- Rfast::Diag.matrix(D, a1[1,i]/B1[i,])
            P233[,i] <- -0.5*quadratic_form_diag(X, temp)
            P230[,i] <- mat_mult_t(L1[i,,drop=FALSE], temp, X)
            P231[1,i] <- -0.5*mat_mult_t(L1[i,,drop=FALSE], temp,
                                              L1[i,,drop=FALSE])
            P232[1,i] <- 0.5*sum(digamma(a1[1,i]) - log(B1[i,]))
          }
          P_const <- -0.5*D*(1/(1/k0 + RP))
          #log probability matrix update
          Plog <- P2 + P230 + matrix((P_const + P231 + P232), nrow = N, ncol = T0,
                                     byrow = TRUE) + P233

          #log-sum-exp trick
          Plog <- t(apply(Plog, 1, function(x){
            mx <- max(x)
            x - mx - log(sum(exp(x - mx)))
          }))
          P <- exp(Plog)
          RP <- Rfast::colsums(P)
          ord <- order(RP, decreasing = TRUE)
          Plog <- Plog[,ord]
          P <- P[,ord]
          L1 <- L1[ord,]
          a1 <- a1[1,ord, drop = FALSE]
          B1 <- B1[ord,]
          C1 <- C1[,,ord]
          
          params$post_shape_d_cs_cov <- a1
          params$post_rate_d_cs_cov <- B1
          params$post_mean_eta <- L1
          params$log_prob_matrix <- Plog
          params$P <- P
          out <- params

        } else if (variance_prior_type == "off-diagonal normal"){

          a0 <- params$prior_shape_d_cs_cov
          b0 <- params$prior_rate_d_cs_cov
          c0 <- params$prior_var_offd_cs_cov
          a1 <- params$post_shape_d_cs_cov
          B1 <- params$post_rate_d_cs_cov
          C1 <- params$post_mean_offd_cs_cov
          k0 <- params$scaling_cov_eta
          
          RP <- Rfast::colsums(P)
          
          a1 <- matrix(a0 + RP, nrow = 1, ncol = T0)
          for (i in 1:T0){
            B1[i,] <- b0 + Rfast::eachcol.apply(X^2, P[,i], oper = "*")
            C01 <- -0.5*c0*t_mat_mult(X, diag(P[,i]), X)
            C1[,,i] <- Rfast::Diag.fill(C01, rep(0, D))
            L1[i,] <- (Mu0/k0 +
                         Rfast::eachcol.apply(X, P[,i], oper = "*"))/(1/k0 + RP[i])
          }
          
          #expectation of inverse of data covariance matrix
          inv_C0 <- array(0, c(D, D, T0))

          #updating the latent probability values
          P230 <- matrix(0, nrow = N, ncol = T0)
          P233 <- matrix(0, nrow = N, ncol = T0)
          P231 <- matrix(0, nrow = 1, ncol = T0)
          P232 <- P231
          for (i in 1:T0){
            inv_C0[,,i] <- temp <- Rfast::Diag.fill(C1[,,i], a1[1,i]/B1[i,])
            P233[,i] <- -0.5*quadratic_form_diag(X, temp)
            P230[,i] <- mat_mult_t(L1[i,,drop=FALSE], temp, X)
            P231[1,i] <- -0.5*mat_mult_t(L1[i,,drop=FALSE], temp,
                                         L1[i,,drop=FALSE])
            P232[1,i] <- 0.5*sum(digamma(a1[1,i]) - log(B1[i,]))
          }
          P_const <- -0.5*D*(log(2*pi) + 1/(1/k0 + RP))
          #log probability matrix update
          Plog <- P2 + P230 + matrix((P_const + P231 + P232), nrow = N, ncol = T0,
                                     byrow = TRUE) + P233
          #log-sum-exp trick
          Plog <- t(apply(Plog, 1, function(x){
            mx <- max(x)
            x - mx - log(sum(exp(x - mx)))
          }))
          P <- exp(Plog)
          RP <- Rfast::colsums(P)
          ord <- order(RP, decreasing = TRUE)
          Plog <- Plog[,ord]
          P <- P[,ord]
          
          params$post_shape_d_cs_cov <- a1
          params$post_rate_d_cs_cov <- B1
          params$post_mean_offd_cs_cov <- c0*C1
          params$post_mean_eta <- L1
          params$log_prob_matrix <- Plog
          params$P <- P
          out <- params

        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed'
               when 'cluster_specific_covariance' is TRUE")
        }
      }
    }
  } else {
    stop("covariance_type can only be either 'diagonal' or 'full'.")
  }
  params$P <- P
  
  #update for concentration parameter alpha
  C0sum <- Rfast::colsums(P)
  index <- which(C0sum >= 1)
  l0 <- max(index)
  #update of the shape parameter of alpha
  W1 <- s1 + l0 - 1
  #update of the rate parameter of alpha
  alpha0 <- l0/log(N)
  if (l0 > 1){
    a_eni <- Rfast::colsums(P[,1:l0])
    a_vni <- Rfast::colsums(P[,1:l0]*(1-P[,1:l0]))
    a_enj <- cum_clustprop(P[,1:l0])
    a_vnj <- cum_clustprop_var(P[,1:l0])
    W20 <- log(alpha0 + a_eni[1:(l0 - 1)] + a_enj[1:(l0 - 1)]) -
      0.5*(a_vni[1:(l0 - 1)] +
             a_vnj[1:(l0 - 1)])/((alpha0 + a_eni[1:(l0 - 1)] +
                                    a_enj[1:(l0 - 1)])^2) -
      log(alpha0 + a_enj[1:(l0 - 1)]) +
      0.5*a_vnj[1:(l0 - 1)]/((alpha0 + a_enj[1:(l0 - 1)])^2)
    W21 <- log(alpha0 + a_eni[l0]) - 0.5*a_vni[l0]/((alpha0 + a_eni[l0])^2) -
      log(alpha0)
    W2f <- sum(W20) + W21
  } else {
    W2f <- log(alpha0 + sum(P[,l0])) -
      0.5*sum(P[,l0]*(1-P[,l0]))/((alpha0 + sum(P[,l0]))^2) - log(alpha0)
  }
  W2 <- s2  + W2f
  
  params$post_shape_alpha <- W1
  params$post_rate_alpha <- W2

  return(params)
}
