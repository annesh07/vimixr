#optimisation of C0 parameters
#approximated by IW (V_1, nu)
#so C0_1 has W(V, nu)
#initial values, a0, b0, k0, Lambda0, truncated from 0.01 to D*N
D <- 20
N <- 100
a0 <- 0.001
b0 <- 0.001
k0 <- 2
Lambda0 <- diag(D)
prob_function <- function(nu, V){
  #V <- solve(V_1)
  inv_L <- solve(Lambda0)
  k1 <- k0/2
  k2 <- (D - 1)*k1 + k1*(0.5*D*(D+1) - 1) + D*k1*(k1 - 1)
  a1 <- a0 + D*k1
  b1 <- b0 + 0.5*(determinant((V + inv_L), logarithm = TRUE)$modulus -
                    determinant(V, logarithm = TRUE)$modulus)
  psi_nu <- sum(digamma(0.5*(nu + 1 - c(1:D))))
  F_00 <- pgamma(0.01, shape = a0, rate = b0)
  F_01 <- pgamma(D*N, shape = a0, rate = b0)
  F_10 <- pgamma(0.01, shape = a1, rate = b1)
  F_11 <- pgamma(D*N, shape = a1, rate = b1)

  # logf <- 0.5*(D+1)*(psi_nu + D*log(2) + determinant(V, logarithm = TRUE)$modulus) -
  #   0.25*D*(D-1)*log(pi) - sum(lgamma(0.5*(k0 + 1 - c(1:D)))) -
  #   0.5*k0*(psi_nu + D*log(2) + determinant((diag(D) + Lambda0%*%V), logarithm = TRUE)$modulus) +
  #   a0*log(b0) - lgamma(a0) - log(F_01 - F_00) -
  #   a1*log(b1) + lgamma(a1) + log(F_11 - F_10) + k2*log(D*N/0.01)

  logf <- 0.5*(D+1)*(psi_nu + D*log(2) + determinant(V, logarithm = TRUE)$modulus)-
    0.5*k0*(psi_nu + D*log(2) + determinant((diag(D) + Lambda0%*%V), logarithm = TRUE)$modulus)-
    a1*log(b1)

  return(logf)
}

var_function <- function(nu, V){
  #V <- solve(V_1)
  psi_nu <- sum(digamma(0.5*(nu + 1 - c(1:D))))
  gamma_nu <- 0.25*D*(D-1)*log(pi) + sum(lgamma(0.5*(nu + 1 - c(1:D))))

  logq <- -0.5*D*nu*log(2) - gamma_nu - 0.5*nu*determinant(V, logarithm = TRUE)$modulus +
    0.5*(nu+D+1)*(psi_nu + D*log(2) + determinant(V, logarithm = TRUE)$modulus) -
    0.5*D*nu

  return(logq)
}

ELBO <- function(params){
  LM <- matrix(params[-length(params)], nrow = D)
  V <- LM %*% t(LM)
  nu <- params[length(params)]

  e0 <- prob_function(nu, V) - var_function(nu, V) -
    sum(diag(LM)^2 + 10000/abs(diag(LM)))
  return(e0)
}

Bound <- matrix(0, nrow = D, ncol = D)
Bound[upper.tri(diag(D), diag=TRUE)] <- 1
Lbound <- Bound*(-round((max(cov(X)))^(1/D)))
Ubound <- Bound*(round((max(cov(X)))^(1/D)))
library("GA")
#param_init <- c(10,10,0,10,4)
set.seed(22112024)
Results <- ga(type = "real-valued", fitness = ELBO,
              lower = c(Lbound, 2),
              upper = c(Ubound, max(N,D)+2), maxiter = 10000,
              popSize = N*D, elitism = 100, run = 100,
              pcrossover = 0.8, pmutation = 0.1)
LM0 <- matrix(Results@solution[-length(Results@solution)], nrow = D)
nu0 <- Results@solution[length(Results@solution)]
V0 <- LM0%*%t(LM0)
det_V0 <- det(V0)
if (det_V0 < 1){
  V0 <- V0*(det_V0)
} else {
  V0 <- V0/(det_V0)
}
