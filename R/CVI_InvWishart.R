CVI_InvWishart <- function(N, D, T0, s1, s2, C00, nu0, Lamda0, k0, Lamda, k,
                           X, W1, W2, L1, L2, nu, V, Plog, maxit){
  #the mean vector of the parameters eta_i
  Mu0 <- matrix(c(rep(0,D)), nrow=1)
  inv_C00 <- solve(C00)
  Mu00 <- Mu0 %*% inv_C00
  #store the output of ELBO function for every iteration of updates
  f <- list()
  f[[1]] <- ELBO_InvWishart(N, D, T0, s1, s2, C00, nu0, Lamda0, k0, Lamda, k,
                            X, W1, W2, L1, L2, nu, V, Plog)

  for (m in 1:maxit){
    L21 <- matrix(0, nrow = T0, ncol = D)
    for (i in 1:T0){
      L21[i,] = L1[i,, drop = FALSE] %*% solve(L2[,,i])
    }
    EC0_inv <- nu*V
    #updating the latent probability values
    P0 <- exp(Plog)
    P230 <- diag(X %*% EC0_inv %*% t(X))
    P231 <- L21 %*% EC0_inv %*% t(X)
    P232 <- -0.5*diag(L21 %*% EC0_inv %*% t(L21))
    P233 <- apply(L2, 3, function(x){-0.5*sum(diag(EC0_inv %*% solve(x)))})
    P234 <- sum(digamma(0.5*(nu + 1 - c(1:D)))) +
      D*log(2) + determinant(V, logarithm = TRUE)$modulus
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

      P2 <- P20 + P22 + P231[,n] + P232 + P233 - 0.5*(P230[n] + D*log(2*pi)
                                                      - P234)
      #log-sum-exp trick
      p0 <- max(P2)
      Plog[n,] <- P2 - p0 - log(sum(exp(P2 - p0)))
    }
    #labelling the probability matrix so that non-zero cluster allocations
    #are present in the beginning of the matrix
    P00 <- exp(Plog)
    Csum <- colSums(P00)
    index <- which(Csum > 0.000001)
    l0 <- length(index)
    for (l in 1:l0){
      Plog[, c(l, index[l])] <- Plog[, c(index[l], l)]
    }
    #final updated and labelled probability matrix
    Pf0 <- exp(Plog)

    for (i in 1:T0){
      #update of the 1st parameter vector of eta_i
      L1[i,] <- Mu00 + (t(Pf0[, i, drop=FALSE])%*%X) %*% EC0_inv
      #update of the 2nd parameter matrix of eta_i
      L2[,,i] <- inv_C00 + sum(Pf0[, i])*EC0_inv
    }

    L21 <- matrix(0, nrow = T0, ncol = D)
    for (i in 1:T0){
      L21[i,] <- L1[i,, drop = FALSE] %*% solve(L2[,,i])
    }

    #update of degree of freedom of V0_inv
    k <- k0 + nu0
    #update of scalar matrix of V0_inv
    Lamda <- solve(solve(Lamda0) + (nu*V))

    #update of degree of freedom of C0_inv
    nu <- nu0 + 0.5*sum(Pf0)
    #update of matrix of C0_inv
    RPf0 <- rowSums(Pf0)
    CPf0 <- colSums(Pf0)
    V00 <- matrix(0, nrow = D, ncol = D)
    for (n in 1:N){
      V00 <- V00 + (RPf0[n])*(t(X[n,, drop = FALSE]) %*% X[n,, drop = FALSE])
    }
    for (n in 1:N){
      for (i in 1:T0){
        V00 <- V00 - (Pf0[n,i])*2*(t(X[n,, drop = FALSE]) %*% L21[i,, drop = FALSE])
      }
    }
    for (i in 1:T0){
      V00 <- V00 + (CPf0[i])*(t(L21[i,, drop = FALSE]) %*% L21[i,, drop = FALSE])
    }
    for (i in 1:T0){
      V00 <- V00 + CPf0[i]*solve(L2[,,i])
    }
    V <- solve((k*Lamda) + V00)

    #update of the shape parameter of alpha
    #a0 <- l0/log(N)
    W1 <- s1 + l0 - 1
    #update of the rate parameter of alpha
    #a0 <- l0/(log(N)^(l0+log(N)))
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

    f[[m+1]] <- ELBO_InvWishart(N, D, T0, s1, s2, C00, nu0, Lamda0, k0, Lamda, k,
                                X, W1, W2, L1, L2, nu, V, Plog)
    if (abs(sum(f[[m]]) - sum(f[[m + 1]])) < 0.000001){
      break
    }
    message("outer loop: ", m,"\n", f[[m + 1]], '\n', sep="")
  }
  alpha0 <- W1/W2
  dof <- nu
  sigma0 <- nu*V
  clustering <- apply(Plog, MARGIN = 1, FUN=which.max)
  clust <- table(clustering)
  clustnum <- length(unique(clustering))

  posterior <- list("alpha"=alpha0, "degrees of freedom"=dof,
                    "sigma^2"=sigma0, "Clusters"= clustnum,
                    "Proportions"=clust, "Clustering" = Plog)
  optimisation <- list("ELBO" = f)

  output <-  list("posterior" = posterior, "optimisation" = optimisation)
  class(output) <- "CVIoutput"

  return(output)
}
