#' General ELBO function
#'
#' @param fixed_variance whether the covariance is fixed or estimated.
#' Default is \code{FALSE} which means it is estimated.
#' @param covariance_type The assumed type of the covariance matrix. Can be either \code{"diagonal"} if it is the identify multiplied by a scalar,
#' or \code{"full"} for a fully unspecified covariance matrix.
#' @param cluster_specific_covariance whether the the covariance is shared across
#' estimated clusters or is cluster specific. Default is \code{TRUE} which means it is cluster specific.
#' @param variance_prior_type character string specifying the type of prior distribution
#' for the covariance when cluster_specific_covariance is \code{TRUE}.
#' Can be either \code{"IW"} or \code{"decomposed} if \code{cluster_specific_covariance} is \code{FALSE},
#' and can be either \code{"IW"}, \code{"sparse"} or \code{"off-diagonal normal"} otherwise.
#' @param params a list
#'
#' @returns
#'
#' @importFrom Rfast rowsums colsums
#'
#' @export
#'
#' @examples
ELBO_function <- function(fixed_variance=FALSE, covariance_type="diagonal",
                          cluster_specific_covariance=TRUE,
                          variance_prior_type=c("IW", "decomposed", "sparse", "off-diagonal normal"),
                          X,
                          inverts,
                          params){


  s1 <- params[["s1"]] #shape1 parameter for alpha prior
  s2 <- params[["s2"]] #shape2 parameter for alpha prior
  W1 <- params[["W1"]] #posterior shape1 parameter for alpha
  W2 <- params[["W2"]] #posterior shape2 parameter for alpha

  Plog <- params[["Plog"]] #log of posterior probability allocation matrix
  P <- params[["P"]]

  #expectation of alpha prior
  e_alpha <- s1*log(s2) - lgamma(s1) + (s1 - 1)*(digamma(W1) - log(W2)) - s2 * (W1 / W2)

  #expectation of latent probability allocation prior
  cp <- Rfast::colsums(P) #cluster proportions
  v_cp <- Rfast::colsums(P*(1 - P)) #variance of the cluster proportions
  ccp <- Rfast::rowsums(apply(P, 1, f0)) #cumulative cluster proportions
  v_ccp <- Rfast::rowsums(apply(P, 1, f1))  #variance of the cumulative cluster proportions
  e_indiv_alloc <- lgamma(1 + cp) + 0.5 * trigamma(1 + cp) * v_cp +
    lgamma(W1 / W2 + enj) + 0.5 * trigamma(W1 / W2 + enj) * ((W1 / W2^2) + v_ccp) -
    lgamma(1 + W1 / W2 + cp + enj) -
    0.5*trigamma(1 + W1 / W2 + cp + enj) * (W1 / W2^2 + v_cp + v_ccp)
  e_alloc <- T0 * (digamma(W1) - log(W2)) + sum(e_indiv_alloc)

  #Variational expectation of alpha & latent allocations
  e_alpha_post <- W1 * log(W2) - lgamma(W1) + (W1 - 1)*(-log(W2) + digamma(W1)) - W1
  e_alloc_post <- sum(exp(Plog) * Plog)

  if(covariance_type == "diagonal") {

    if(fixed_variance) {

      fixed_variance_elbo <- elbo_fixed_diagonal(X, inverts, params)
      fixed_variance_elbo$me4 <- fixed_variance_elbo$me4 - e_alpha_post - e_alloc_post
      out <- c("e1" = e_alpha + e_alloc,
               fixed_variance_elbo)

    } else {
      b1 <- params$prior_shape_scalar_cov
      b2 <- params$prior_rate_scalar_cov
      G1 <- params$post_shape_scalar_cov
      G2 <- params$post_rate_scalar_cov
      L2 <- params$post_precision_scalar_eta
      L20 <- params$prior_precision_scalar_eta

      C00 <- diag(D)/L20       #C00; covariance matrix for eta_i's
      inv_C00 <- spdinv(C00)    #inverse of C00

      #the eta's
      L21 <- sweep(L1, 1, L2, "/")
      e20 <- diag(-0.5*L21 %*% inv_C00 %*% t(L21))
      e21 <- - 0.5*(D*L20)/L2
      e22 <- Mu0 %*% inv_C00 %*% t(L21)
      e2 <-  T0*(-D/2*log(2*pi) + D*0.5*log(L20)- 0.5*Mu0 %*% inv_C00 %*% t(Mu0)) +
        sum(e20) + sum(e21) + sum(e22)

      #the X's
      e30 <- sweep(P, 1, -0.5*(G1/G2)*diag(X %*% t(X)), "*")
      e31 <- P*t((G1/G2)*(L21 %*% t(X)))
      e32 <- sweep(P, 2, -0.5*(G1/G2)*diag(L21 %*% t(L21)), "*")
      e33 <- sweep(P, 2, -0.5*(G1/G2)*D/L2, "*")
      e3 <- sum(P*(-0.5*D*log(2*pi) + 0.5*D*(digamma(G1) - log(G2))))
      + sum(e30) + sum(e31) + sum(e32) + sum(e33)

      #the sigma^2
      e4 <- b1*log(b2) - lgamma(b1) + (b1 - 1)*(digamma(G1) - log(G2)) - b2*G1/G2

      #the variationa distributions
      e52 <- sum(D*log(L2)/2) - 0.5*D*T0*(log(2*pi)+1)
      e53 <- G1*log(G2) - lgamma(G1) + (G1-1)*(digamma(G1) - log(G2)) - G1
      e5 <- e52 + e53 + ev

      return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "e4"=e4, "me5"=-e5))
    }

  } else if(covariance_type == "full") {

    if(fixed_variance) {

      C0 <- params$cov_data
      L2 <- params$post_cov_eta
      C00 <- params$prior_cov_eta

      inv_C00 <- spdinv(C00)    #inverse of C00
      inv_C0 <- spdinv(C0)      #inverse of C0

      #the eta's
      L21 <- matrix(0, nrow = T0, ncol = D)
      for (i in 1:T0){
        L21[i,] = L1[i,, drop = FALSE] %*% L2[,,i]
      }
      e20 <- apply(L21, 1, function(x){-0.5*t(x)%*%inv_C00%*%x})
      e21 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C00)*x)})
      e22 <- Mu0 %*% inv_C00 %*% t(L21)
      e2 <-  T0*(-D/2*log(2*pi) - 0.5*determinant(C00, logarithm = TRUE)$modulus -
                   0.5*mat.mult(mat.mult(Mu0, inv_C00), t(Mu0))) +
        sum(e20) + sum(e21) + sum(e22)

      #the X's
      e30 <- eachcol.apply(P, apply(X, 1, function(x){-0.5*t(x)%*%inv_C0%*%x}),
                           oper = "*")
      e31 <- P*t(L21 %*% inv_C0 %*% t(X))
      e32 <- P*matrix(apply(L21, 1, function(x){-0.5*t(x)%*%inv_C0%*%x}),
                      nrow = N, ncol = T0, byrow = TRUE)
      e33 <- apply(L2, 3, function(x){-0.5*sum(t(inv_C0)*x)})
      e34 <- P*matrix(e33, nrow = N, ncol = T0, byrow=TRUE)
      e3 <- N*(-0.5*D*log(2*pi) + 0.5*determinant(inv_C0, logarithm = TRUE)$modulus) +
        sum(e30) + sum(e31) + sum(e32) + sum(e34)

      #the variationa distribution
      e42 <- apply(L2, 3,function(x){-0.5*determinant(x, logarithm = TRUE)$modulus})
      e4 <- sum(e42) - 0.5*D*T0*(log(2*pi)+1) + ev

      return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "me4"=-e4))

    } else {
      if(!cluster_specific_covariance) {
        if(variance_prior_type == "IW"){
          nu0 <- params$prior_df_cov
          V0 <- params$prior_scale_cov
          nu <- params$post_df_cov
          V <- params$post_scale_cov
          L2 <- params$post_cov_eta
          C00 <- params$prior_cov_eta

          inv_C00 <- spdinv(C00)    #inverse of C00
          inv_V0 <- spdinv(V0)      #inverse of V0
          inv_C0 <- nu*V           #expected inverse of C0; covariance matrix of data

          #the eta's
          L21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            L21[i,] = L1[i,, drop = FALSE] %*% L2[,,i]
          }
          e20 <- diag(-0.5*L21 %*% inv_C00 %*% t(L21))
          e21 <- apply(L2, 3, function(x){-0.5*sum(diag(inv_C00 %*% x))})
          e22 <- Mu0 %*% inv_C00 %*% t(L21)
          e2 <-  T0*(-D/2*log(2*pi) - 0.5*determinant(C00, logarithm = TRUE)$modulus -
                       0.5*Mu0 %*% inv_C00 %*% t(Mu0)) + sum(e20) + sum(e21) + sum(e22)

          #the X's
          e30 <- sweep(P, 1, -0.5*diag(X %*% inv_C0 %*% t(X)), "*")
          e31 <- P*t(L21 %*% inv_C0 %*% t(X))
          e32 <- sweep(P, 2, -0.5*diag(L21 %*% inv_C0 %*% t(L21)), "*")
          e33 <- apply(L2, 3, function(x){-0.5*sum(diag(inv_C0 %*% x))})
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
            0.5*sum(diag(inv_V0 %*% inv_C0))

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
          e5 <- e52 + e53 + ev

          return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "e4"=e4, "me5"=-e5))
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
          C00 <- params$prior_cov_eta

          inv_C00 <- spdinv(C00)    #inverse of C00
          mean_lower <- matrix(0, nrow = D, ncol = D) #mean matrix of the decomposed
          mean_lower[lower.tri(mean_lower, diag = FALSE)] <- mu1
          sigma_lower <- matrix(0, nrow = D, ncol = D) #var matrix of the decomposed
          sigma_lower[lower.tri(sigma_lower, diag = FALSE)] <- c1
          mean_L <- mean_lower + diag(sqrt(1/b1)*sqrt(pi)/beta(a1,0.5))
          diag(sigma_lower) <- (1/b1)*(a1 - (sqrt(pi)/beta(a1,0.5))^2)
          #expected inverse of C0; covariance matrix of data
          inv_C0 <- mean_L %*% t(mean_L) + diag(rowsums(sigma_lower))

          #the eta's
          L21 <- matrix(0, nrow = T0, ncol = D)
          for (i in 1:T0){
            L21[i,] = L1[i,, drop = FALSE] %*% L2[,,i]
          }
          e20 <- diag(-0.5*L21 %*% inv_C00 %*% t(L21))
          e21 <- apply(L2, 3, function(x){-0.5*sum(diag(inv_C00 %*% x))})
          e22 <- Mu0 %*% inv_C00 %*% t(L21)
          e2 <-  T0*(-D/2*log(2*pi) - 0.5*determinant(C00, logarithm = TRUE)$modulus -
                       0.5*Mu0 %*% inv_C00 %*% t(Mu0)) + sum(e20) + sum(e21) + sum(e22)

          #the X's
          e30 <- eachcol.apply(P, -0.5*diag(X %*% inv_C0 %*% t(X)), oper = "*")
          e31 <- P*t(L21%*%inv_C0%*%t(X))
          e32 <- eachrow(P, -0.5*diag(L21%*%inv_C0%*%t(L21)), oper = "*")
          e33 <- apply(L2, 3, function(x){-0.5*sum(diag(inv_C0 %*% x))})
          e34 <- eachrow(P, e33, oper = "*")
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
          e5 <- e52 + e530 + e531 + ev

          return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "e4"=e4, "me5"=-e5))
        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed' when 'cluster_specific_covariance' is FALSE")
        }
      }else{
        if(variance_prior_type == "IW"){
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
          #covariance parameter of eta's
          L2 <- sweep(V1, 3, nu1*(1/k0 + RP), "/")

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
            (1/(2*k0))*D*sum(1/(1/k0 + RP)) + sum(e21) + sum(e22) - 0.5*T0*nu0*D*log(2) -
            T0*(0.25*D*(D-1)*log(pi) + sum(lgamma(0.5*(nu0 + 1 - c(1:D))))) +
            0.5*T0*nu0*determinant(V0, logarithm = TRUE)$modulus + sum(e23)

          #the data x_n's
          e30 <- matrix(0, nrow = N, ncol = T0)
          e31 <- matrix(0, nrow = N, ncol = T0)
          for (n in 1:N){
            for (i in 1:T0){
              e30[n,i] <- -0.5*P[n,i]*X[n,,drop=FALSE] %*% E_C0[,,i] %*% t(X[n,,drop=FALSE])
              e31[n,i] <- P[n,i]*L1[i,,drop=FALSE] %*% E_C0[,,i] %*% t(X[n,,drop=FALSE])
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
            0.5*D*sum(nu1) + ev

          return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "me4"=-e4))

        } else if (variance_prior_type == "sparse"){

          a0 <- params$prior_scale_d_cs_cov
          b0 <- params$prior_rate_d_cs_cov
          c0 <- params$prior_var_offd_cs_cov
          a1 <- params$post_scale_d_cs_cov
          B1 <- params$post_rate_d_cs_cov
          C1 <- params$post_var_offd_cs_cov
          k0 <- params$scaling_cov_eta

          #expectation of inverse of C0, data covariance matrix
          E_C0_inv <- array(0, c(D, D, T0))
          #inverse of expectation of inverse of C0, data covariance matrix
          E_C0_inv_inv <- array(0, c(D, D, T0))
          for (i in 1:T0){
            E_C0_inv[,,i] <- Diag.matrix(D, a1[1,i]/B1[i,])
            E_C0_inv_inv[,,i] <- Diag.matrix(D, B1[i,]/a1[1,i])
          }
          #covariance parameter of eta's
          L2 <- sweep(E_C0_inv_inv, 3, (1/k0 + RP), "/")

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
              e30[n,i] <- -0.5*P[n,i]*X[n,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(X[n,,drop=FALSE])
              e31[n,i] <- P[n,i]*L1[i,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(X[n,,drop=FALSE])
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
          e4 <- sum(e420) + sum(-log(2*C1[!diag(D)]) - 1) + ev

          return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "me4"=-e4))

        } else if (variance_prior_type == "off-diagonal normal"){
          a0 <- params$prior_scale_d_cs_cov
          b0 <- params$prior_rate_d_cs_cov
          a1 <- params$post_scale_d_cs_cov
          B1 <- params$post_rate_d_cs_cov
          C1 <- params$post_mean_offd_cs_cov
          k0 <- params$scaling_cov_eta

          #expectation of inverse of data covariance matrix
          E_C0_inv <- array(0, c(D, D, T0))
          #inverse of expectation of inverse of data covariance matrix
          E_C0_inv_inv <- array(0, c(D, D, T0))
          for (i in 1:T0){
            E_C0_inv[,,i] <- Diag.fill(C1[,,i], a1[1,i]/B1[i,])
            E_C0_inv_inv[,,i] <- spdinv(E_C0_inv[,,i])
          }
          #covariance parameter of eta's
          L2 <- sweep(E_C0_inv_inv, 3, (1/k0 + RP), "/")

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
              e30[n,i] <- -0.5*P[n,i]*X[n,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(X[n,,drop=FALSE])
              e31[n,i] <- P[n,i]*L1[i,,drop=FALSE] %*% E_C0_inv[,,i] %*% t(X[n,,drop=FALSE])
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
          e4 <- sum(e420) + (D^2 - D)*(-0.5*log(2*pi) - 0.5) + ev

          return(c("e0"=e0, "e1"=e1, "e2"=e2, "e3"=e3, "me4"=-e4))
        } else {
          stop("'variance_prior_type' can only be either 'IW' or 'decomposed' when 'cluster_specific_covariance' is TRUE")
        }
      }


    } else {
      stop("covariance_type can only be either 'diagonal' or 'full'.")
    }
  }

  return("ELBO" = out)

}


