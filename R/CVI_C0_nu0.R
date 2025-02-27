CVI_C0_nu0 <- function(N, D, T0, s1, s2, C00, b10, b20, V0, X,
                   W1, W2, L1, L2, Plog, b2, V, maxit){
  #C0 <- diag(D)
  #inv_C0 <- solve(C0)
  #d_inv_C0 <- diag(inv_C0)
  #the mean vector of the parameters eta_i
  Mu0 <- matrix(c(rep(0,D)), nrow=1)
  #the covariance matrix of the parameters eta_i
  #C00 <- diag(D)/L20
  Mu00 <- Mu0%*%solve(C00)
  #store the output of ELBO function for every iteration of updates
  f <- list()
  f[[1]] <- ELBO_C0_nu0(N, D, T0, s1, s2, C00, b10, b20, V0, X,
                    W1, W2, L1, L2, Plog, b2, V)

  for (m in 1:maxit){
    nu <- b10/(b2 + 0.5*D*(log(2) + log(0.5) -
                             determinant(solve(V0), logarithm = TRUE)$modulus))
    nu0 <- b10/(b20 + 0.5*D*(log(2) + log(0.5) -
                               determinant(solve(V0), logarithm = TRUE)$modulus))
    inv_C0 <- nu*V
    #updating the latent probability values
    P0 <- exp(Plog)
    #different updates for i = 1, i = {2, ..., T0-1} and i = T0
    L21 <- matrix(0, nrow = T0, ncol = D)
    for (i in 1:T0){
      L21[i,] = L1[i,, drop = FALSE] %*% solve(L2[,,i])
    }
    P230 <- L21 %*% inv_C0 %*% t(X)
    P231 <- diag(-0.5*L21 %*% inv_C0 %*% t(L21))
    P232 <- apply(L2, 3, function(x){-0.5*sum(diag(inv_C0 %*% solve(x)))})
    P233 <- diag(X %*% inv_C0 %*% t(X))
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

      P2 <- P20 + P22 + P230[,n] + P231 + sum(P232) -
        0.5*(P233[n] + D*log(2*pi) - (sum(digamma(0.5*(nu + 1 - c(1:D)))) +
                                        D*log(2) + determinant(V, logarithm = TRUE)$modulus))
      #log-sum-exp trick
      p0 <- max(P2)
      Plog[n,] <- P2 - p0 - log(sum(exp(P2 - p0)))
    }
    #labelling the probability matrix so that non-zero cluster allocations
    #are present in the beginning of the matrix
    P00 <- exp(Plog)
    Csum <- colSums(P00)
    index <- which(Csum > 1)
    l0 <- length(index)
    for (l in 1:l0){
      Plog[, c(l, index[l])] <- Plog[, c(index[l], l)]
    }
    #final updated and labelled probability matrix
    Pf0 <- exp(Plog)

    for (i in 1:T0){
      #update of the 1st parameter vector of eta_i, including C0 matrix
      L1[i,] <- Mu00 + (t(Pf0[, i, drop=FALSE]) %*% X) %*% inv_C0
      #update of the 2nd parameter value of eta_i, just fr calculation purpose,
      #C0 is not included; but here if implies that C0 must be a scalar times
      #identity matrix
      L2[,,i] <- solve(C00) + sum(Pf0[, i])*inv_C0
    }
    L21 <- matrix(0, nrow = T0, ncol = D)
    for (i in 1:T0){
      L21[i,] = L1[i,, drop = FALSE] %*% solve(L2[,,i])
    }

    #update of the shape parameter of alpha
    W1 <- s1 + l0 - 1
    #update of the rate parameter of alpha
    a0 <- l0/log(N)
    a_eni <- colSums(Pf0[,1:l0])
    a_vni <- colSums(Pf0[,1:l0]*(1-Pf0[,1:l0]))
    a_enj <- rowSums(apply(Pf0[,1:l0], 1, f0))
    a_vnj <- rowSums(apply(Pf0[,1:l0], 1, f1))
    W20 <- log(a0 + a_eni[1:(l0 - 1)] + a_enj[1:(l0 - 1)]) -
      0.5*(a_vni[1:(l0 - 1)] + a_vnj[1:(l0 - 1)])/((a0 + a_eni[1:(l0 - 1)]
                                                    + a_enj[1:(l0 - 1)])^2) - log(a0 + a_enj[1:(l0 - 1)]) +
      0.5*a_vnj[1:(l0 - 1)]/((a0 + a_enj[1:(l0 - 1)])^2)
    W21 <- log(a0 + a_eni[l0]) - 0.5*a_vni[l0]/((a0 + a_eni[l0])^2) -
      log(a0)
    W2 <- s2  + sum(W20) + W21

    #updates of nu0
    b2 <- b20 + 0.5*D*(log(2) + log(0.5) -
                         determinant(solve(V0), logarithm = TRUE)$modulus) -
      0.5*(sum(digamma(0.5*(nu + 1 - c(1:D)))) + D*log(2) +
             determinant(V, logarithm = TRUE)$modulus)

    #updates of C0
    V1 <- solve(V0)
    RPf0 <- rowSums(Pf0)
    CPf0 <- colSums(Pf0)
    for (n in 1:N){
      V1 <- V1 + (RPf0[n])*(t(X[n,, drop = FALSE]) %*% X[n,, drop = FALSE])
    }
    for (n in 1:N){
      for (i in 1:T0){
        V1 <- V1 - (Pf0[n,i])*2*(t(X[n,, drop = FALSE]) %*% L21[i,, drop = FALSE])
      }
    }
    for (i in 1:T0){
      V1 <- V1 + (CPf0[i])*(t(L21[i,, drop = FALSE]) %*% L21[i,, drop = FALSE] +
                              solve(L2[,,i]))
    }
    V <- solve(V1)

    f[[m+1]] <- ELBO_C0_nu0(N, D, T0, s1, s2, C00, b10, b20, V0, X,
                        W1, W2, L1, L2, Plog, b2, V)
    if (abs(sum(f[[m]]) - sum(f[[m + 1]])) < 0.000001){
      break
    }
    message("outer loop: ", m,"\n", f[[m + 1]], '\n', sep="")
  }

  alpha0 <- W1/W2
  clustering <- apply(Plog, MARGIN = 1, FUN=which.max)
  clust <- table(clustering)
  clustnum <- length(unique(clustering))
  Scale <- V
  dof <- b10/(b2 + 0.5*D*(log(2) + log(0.5) -
                            determinant(solve(V0), logarithm = TRUE)$modulus))

  posterior <- list("alpha"=alpha0, "Clusters"=clustnum, "Proportions"=clust,
                    "Clustering" = Plog, "Scale matrix" = Scale,
                    "degrees of freedom" = dof)
  optimisation <- list("ELBO" = f)

  output <-  list("posterior" = posterior, "optimisation" = optimisation)
  class(output) <- "CVIoutput"

  return(output)
}
