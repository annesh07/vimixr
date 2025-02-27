CVI_C0_lower <- function(N, D, T0, s1, s2, C00, a0, b0, mu0, c0, X,
                         W1, W2, L1, L2, Plog, a1, b1, mu1, c1, maxit){
  library("Rfast")
  #the mean vector of the parameters eta_i
  Mu0 <- matrix(c(rep(0,D)), nrow=1)
  #the covariance matrix of the parameters eta_i
  inv_C00 <- spdinv(C00)
  Mu00 <- mat.mult(Mu0, inv_C00)
  #store the output of ELBO function for every iteration of updates
  f <- list()
  f[[1]] <- ELBO_C0_lower(N, D, T0, s1, s2, C00, a0, b0, mu0, c0, X,
                          W1, W2, L1, L2, Plog, a1, b1, mu1, c1)
  for (m in 1:maxit){
    L21 <- matrix(0, nrow = T0, ncol = D)
    for (i in 1:T0){
      L21[i,] = mat.mult(L1[i,, drop = FALSE], L2[,,i])
    }
    Pf0 <- exp(Plog)
    #updates of C0
    a1 <- rep(a0, D) + 0.5*sum(Pf0)

    b21 <- matrix(0, nrow = T0, ncol = D)
    colPf0 <- colsums(Pf0)
    for (i in 1:T0){
      b20 <- eachrow(X, L21[i,], oper = "-")
      b21[i,] <- colsums(Pf0[,i]*(b20^2)) + colPf0[i]*diag(L2[,,i])
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
                            rowsums(Pf0)*X[,k, drop=FALSE], oper = "*")

      mu20 <- rep(0, (k-1))
      for (n in 1:N){
        for (i in 1:T0){
          mu20 <- mu20 +
            Pf0[n,i]*(L21[i, 1:(k-1), drop=FALSE]*X[n,k] +
                        X[n, 1:(k-1), drop=FALSE]*L21[i,k])
        }
      }
      col_Pf0 <- colsums(Pf0)
      mu30 <- rep(0, (k-1))
      for (i in 1:T0){
        mu30 <- mu30 +
          col_Pf0[i]*(L2[,,i][k, 1:(k-1)] + L21[i, 1:(k-1), drop=FALSE]*L21[i, k])
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
    E_C0 <- mat.mult(mean_L, t(mean_L)) + diag(rowsums(sigma_lower))
    #diag(E_C0) <- a1/b1 + rowsums((mean_lower)^2 + sigma_lower)
    #E_C0 <- solve(E_C0)

    #update for eta_i's
    for (i in 1:T0){
      #update of the 1st parameter vector of eta_i, including C0 matrix
      L1[i,] <- Mu00 + mat.mult(mat.mult(t(Pf0[, i, drop=FALSE]), X), E_C0)
      #update of the 2nd parameter value of eta_i, just fr calculation purpose,
      #C0 is not included; but here if implies that C0 must be a scalar times
      #identity matrix
      L2[,,i] <- spdinv(inv_C00 + sum(Pf0[, i])*(E_C0))
    }
    L21 <- matrix(0, nrow = T0, ncol = D)
    for (i in 1:T0){
      L21[i,] = mat.mult(L1[i,, drop = FALSE], L2[,,i])
    }

    #updating the latent probability values
    P0 <- exp(Plog)
    P230 <- mat.mult(mat.mult(L21, E_C0), t(X))
    P231 <- diag(-0.5*mat.mult(mat.mult(L21, E_C0), t(L21)))
    P232 <- apply(L2, 3, function(x){-0.5*sum(diag(mat.mult(E_C0, x)))})
    P233 <- diag(mat.mult(mat.mult(X, E_C0), t(X)))
    for (n in 1:N){
      #update of the n^th vector is done by considering all the vecors except
      #the n^th one
      P1 <- P0[-n,]

      p_eni <- colsums(P1)
      p_vni <- colsums(P1*(1-P1))
      p_enj <- rowsums(apply(P1, 1, f0))
      p_vnj <- rowsums(apply(P1, 1, f1))

      P20 <- log(1 + p_eni) - p_vni/((1 + p_eni)^2) - log(1 + p_eni + p_enj +
                                                            (W1/W2)) + (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)

      P21 <- log((W1/W2) + p_enj) - (p_vnj + (W1/(W2^2)))/(((W1/W2) + p_enj)^2) -
        log(1 + p_eni + p_enj + (W1/W2)) +
        (p_vni + p_vnj + (W1/(W2^2)))/((1 + p_eni + p_enj + (W1/W2))^2)
      P22 <- c(0, cumsum(P21)[1:(T0-1)])

      P2 <- P20 + P22 + P230[,n] + P231 + (P232) -
        0.5*(P233[n] + D*log(2*pi) - sum(digamma(a1) - log(b1)))
      #log-sum-exp trick
      p0 <- max(P2)
      Plog[n,] <- P2 - p0 - log(sum(exp(P2 - p0)))
    }
    #labelling the probability matrix so that non-zero cluster allocations
    #are present in the beginning of the matrix
    Csum <- colsums(exp(Plog))
    ord <- order(Csum, decreasing = TRUE)
    Plog <- Plog[,ord]
    L1 <- L1[ord,]
    L2 <- L2[,,ord]
    Pf0 <- exp(Plog)
    C0sum <- colsums(Pf0)
    index <- which(C0sum > 1)
    l0 <- length(index)
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
        0.5*(a_vni[1:(l0 - 1)] + a_vnj[1:(l0 - 1)])/((alpha0 + a_eni[1:(l0 - 1)]
                                                      + a_enj[1:(l0 - 1)])^2) - log(alpha0 + a_enj[1:(l0 - 1)]) +
        0.5*a_vnj[1:(l0 - 1)]/((alpha0 + a_enj[1:(l0 - 1)])^2)
      W21 <- log(alpha0 + a_eni[l0]) - 0.5*a_vni[l0]/((alpha0 + a_eni[l0])^2) -
        log(alpha0)
      W2f <- sum(W20) + W21
    } else {
      W2f <- log(alpha0 + sum(Pf0[,l0])) -
        0.5*sum(Pf0[,l0]*(1-Pf0[,l0]))/((alpha0 + sum(Pf0[,l0]))^2) - log(alpha0)
    }

    W2 <- s2  + W2f

    f[[m+1]] <- ELBO_C0_lower(N, D, T0, s1, s2, C00, a0, b0, mu0, c0, X,
                              W1, W2, L1, L2, Plog, a1, b1, mu1, c1)
    if (abs(sum(f[[m]]) - sum(f[[m + 1]])) < 0.000001 ){
      break
    }
    message("outer loop: ", m,"\n", f[[m + 1]], '\n', sep="")
  }
  alpha <- W1/W2
  clustering <- apply(Plog, MARGIN = 1, FUN=which.max)
  clust <- table(clustering)
  clustnum <- length(unique(clustering))
  mean_lower <- matrix(0, nrow = D, ncol = D)
  mean_lower[lower.tri(mean_lower, diag = FALSE)] <- mu1
  sigma_lower <- matrix(0, nrow = D, ncol = D)
  sigma_lower[lower.tri(sigma_lower, diag = FALSE)] <- c1
  mean_L <- mean_lower + diag(sqrt(1/b1)*sqrt(pi)/beta(a1,0.5))
  diag(sigma_lower) <- (1/b1)*(a1 - (sqrt(pi)/beta(a1,0.5))^2)
  E_C0 <- mat.mult(mean_L, t(mean_L)) + diag(rowsums(sigma_lower))
  #diag(E_C0) <- a1/b1 + rowsums((mean_lower)^2 + sigma_lower)
  #E_C0 <- solve(E_C0)

  posterior <- list("alpha"=alpha, "Clusters"=clustnum, "Proportions"=clust,
                    "Clustering" = Plog, "Covariance" = E_C0)
  optimisation <- list("ELBO" = f)

  output <-  list("posterior" = posterior, "optimisation" = optimisation)
  class(output) <- "CVIoutput"

  return(output)
}
