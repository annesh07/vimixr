CVI_C0_i <- function(N, D, T0, s1, s2, Mu0, nu0, V0, k0, X,
                    W1, W2, L1, nu1, V1, Plog, maxit){
  library("Rfast")
  #store the output of ELBO function for every iteration of updates
  f <- list()
  f[[1]] <- ELBO_C0_i(N, D, T0, s1, s2, Mu0, nu0, V0, k0, X,
                     W1, W2, L1, nu1, V1, Plog)
  for (m in 1:maxit){
    P0 <- exp(Plog)
    RP0 <- colsums(P0)

    #update for eta_i's and C0_i's
    nu1 <- nu0 + matrix(RP0, nrow=1)
    for (i in 1:T0){
      V1[,,i] <- matrix(0, D, D)
      for (n in 1:N){
        V1[,,i] <- V1[,,i] + P0[n,i]*t(X[n,,drop=FALSE])%*%X[n,,drop=FALSE]
      }
      V1[,,i] <- V0 + (1/k0)*t(Mu0)%*%Mu0 + V1[,,i] + diag(1e-6, D)
      L1[i,] <- (Mu0/k0 + colsums(sweep(X, 1, P0[,i], "*")))/(1/k0 + RP0[i])
    }

    V1_inv <- array(apply(V1, 3, function(x){spdinv(x)}), dim = dim(V1))
    E_C0 <- sweep(V1_inv, 3, nu1, "*")
    E_log_C0_1 <- apply(nu1, 2, function(x){sum(digamma(0.5*(x + 1 - c(1:D))))})
    E_log_C0_2 <- apply(V1_inv, 3, function(x){D*log(2) +
        determinant(x, logarithm = TRUE)$modulus})
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
    for (n in 1:N){
      #update of the n^th vector is done by considering all the vecors except
      #the n^th one
      P1 <- P0[-n,]

      p_eni <- colsums(P1)
      p_vni <- colsums(P1*(1-P1))
      p_enj <- rowsums(apply(P1, 1, f0))
      p_vnj <- rowsums(apply(P1, 1, f1))

      P20 <- log(1 + p_eni) - p_vni/((1 + p_eni)^2) -
        log(1 + p_eni + p_enj + (W1/W2)) +
        (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)

      P21 <- log((W1/W2) + p_enj) -
        (p_vnj + (W1/(W2^2)))/(((W1/W2) + p_enj)^2) -
        log(1 + p_eni + p_enj + (W1/W2)) +
        (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)
      P22 <- c(0, cumsum(P21)[1:(T0-1)])

      P2 <- P20 + P22 + P230[n,] + P231[n,] + P232 -
        0.5*D*(log(2*pi) + 1/(1/k0 + RP0)) + 0.5*E_log_C0
      #log-sum-exp trick
      p0 <- max(P2)
      Plog[n,] <- P2 - p0 - log(sum(exp(P2 - p0)))
    }
    Pf0 <- exp(Plog)
    RPf0 <- colsums(Pf0)

    Csum <- colsums(Pf0)
    ord <- order(Csum, decreasing = TRUE)
    Plog0 <- Plog[,ord]
    Pf0 <- exp(Plog0)
    C0sum <- colsums(Pf0)
    index <- which((C0sum) >= 1)
    l0 <- max(index)
    #update of the shape parameter of alpha
    W1 <- s1 + l0 - 1
    #update of the rate parameter of alpha
    alpha0 <- l0/log(N)
    if (l0 > 1){
      a_eni <- colsums(Pf0[,1:l0])
      a_vni <- colsums(Pf0[,1:l0]*(1-Pf0[,1:l0]))
      a_enj <- rowsums(apply(Pf0[,1:l0], 1, f0))
      a_vnj <- rowsums(apply(Pf0[,1:l0], 1, f1))
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
      W2f <- log(alpha0 + sum(Pf0[,l0])) -
        0.5*sum(Pf0[,l0]*(1-Pf0[,l0]))/((alpha0 + sum(Pf0[,l0]))^2) - log(alpha0)
    }
    W2 <- s2  + W2f

    f[[m+1]] <- ELBO_C0_i(N, D, T0, s1, s2, Mu0, nu0, V0, k0, X,
                         W1, W2, L1, nu1, V1, Plog)
    if (abs(sum(f[[m]]) - sum(f[[m + 1]])) < 0.000001 ){
      break
    }
    message("outer loop: ", m,"\n", f[[m + 1]], '\n', sep="")
  }
  alpha0 <- W1/W2
  clustering <- apply(Plog, MARGIN = 1, FUN=which.max)
  clust <- table(clustering)
  clustnum <- length(unique(clustering))
  posterior <- list("alpha"=alpha0, "Clusters"=clustnum, "Proportions"=clust,
                    "Clustering" = Plog)
  optimisation <- list("ELBO" = f)

  output <-  list("posterior" = posterior, "optimisation" = optimisation)
  class(output) <- "CVIoutput"

  return(output)
}
