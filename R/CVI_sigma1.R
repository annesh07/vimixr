CVI_sigma1 <- function(N, D, T0, s1, s2, L20, b1, b2, X, W1, W2, L1, L2, G1, G2,
                      Plog, maxit){
  #the mean vector of the parameters eta_i
  Mu0 <- matrix(c(rep(0,D)), nrow=1)
  #the covariance matrix of the parameters eta_i
  C00 <- sweep(diag(D), 1, L20, "/")
  Mu00 <- Mu0%*%solve(C00)
  #store the output of ELBO function for every iteration of updates
  f <- list()
  f[[1]] <- ELBO_sigma1(N, D, T0, s1, s2, L20, b1, b2, X, W1, W2, L1, L2, G1, G2,
                       Plog)

  for (m in 1:maxit){
    #updating the latent probability values
    P0 <- exp(Plog)
    #different updates for i = 1, i = {2, ..., T0-1} and i = T0
    L21 <- L1/L2
    P230 <- sweep(L21, 2, (G1/G2), "*") %*% t(X)
    P231 <- rowSums(sweep((L21^2 + 1/L2), 2 , (G1/G2), "*"))
    P232 <- rowSums(sweep(X^2, 2, (G1/G2), "*"))
    for (n in 1:N){
      #update of the n^th vector is done by considering all the vecors except
      #the n^th one
      P1 <- P0[-n,]

      p_eni <- colSums(P1)
      p_vni <- colSums(P1*(1-P1))
      p_enj <- rowSums(apply(P1, 1, f0))
      p_vnj <- rowSums(apply(P1, 1, f1))

      P20 <- log(1 + p_eni) - p_vni/((1 + p_eni)^2) - log(1 + p_eni + p_enj +
        (W1/W2)) + (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)

      P21 <- log((W1/W2) + p_enj) - (p_vnj + (W1/(W2^2)))/(((W1/W2) + p_enj)^2) -
        log(1 + p_eni + p_enj + (W1/W2)) +
        (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)
      P22 <- c(0, cumsum(P21)[1:(T0-1)])

      P2 <- P20 + P22 + P230[,n] + P231 -
        0.5*(P232[n] + D*log(2*pi) - sum(digamma(G1) - log(G2)))
      #log-sum-exp trick
      p0 <- max(P2)
      Plog[n,] <- P2 - p0 - log(sum(exp(P2 - p0)))
    }
    Pf0 <- exp(Plog)

    for (i in 1:T0){
      #update of the 1st parameter vector of eta_i
      L1[i,] <- Mu00 + (G1/G2)*t(Pf0[, i, drop=FALSE])%*%X
      #update of the 2nd parameter value of eta_i
      L2[i,] <- L20 + (G1/G2)*sum(Pf0[, i])
    }

    L21 <- L1/L2
    #update of the shape parameter of sigma^2
    G1 <- b1 + 0.5*sum(Pf0)
    G20 <- 0.5*colSums(sweep(X^2, 1, rowSums(Pf0), "*"))
    #G21 <- -diag(t(L21) %*% t(Pf0) %*% X)
    G21 <- rep(0, D)
    for (n in 1:N){
      for (i in 1:T0){
        G21 <- G21 + Pf0[n,i]*(L21[i,, drop=FALSE] * X[n,, drop=FALSE])
      }
    }
    G22 <- 0.5*colSums(sweep((L21^2 + 1/L2), 1, colSums(Pf0), "*"))
    G2 <- b2 + G20 - G21 + G22

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

    # PP[[m]] = Pf0
    # alpha0[m] <- W1/W2
    # sigma0[[m]] <- G1/G2
    # c1[[m]] <- table(apply(log(Pf0), MARGIN = 1, FUN=which.max))
    # cl[m] <- length(unique(c1[[m]]))
    f[[m+1]] <- ELBO_sigma1(N, D, T0, s1, s2, L20, b1, b2, X, W1, W2, L1, L2, G1,
                           G2, Plog)
    if (abs(sum(f[[m]]) - sum(f[[m + 1]])) < 0.000001){
      break
    }
    message("outer loop: ", m,"\n", f[[m + 1]], '\n', sep="")
  }

  alpha0 <- W1/W2
  sigma0 <- G1/G2
  clustering <- apply(Plog, MARGIN = 1, FUN=which.max)
  clust <- table(clustering)
  clustnum <- length(unique(clustering))

  posterior <- list("alpha"=alpha0, "sigma^2"=sigma0, "Clusters"= clustnum,
                    "Proportions"=clust, "Clustering" = Plog)
  optimisation <- list("ELBO" = f)

  output <-  list("posterior" = posterior, "optimisation" = optimisation)
  class(output) <- "CVIoutput"

  return(output)
}
