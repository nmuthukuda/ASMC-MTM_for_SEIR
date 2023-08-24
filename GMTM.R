getwd()
setwd("C:\\SFU\\Research\\MTM")
set.seed(444)

# MTM Approach

# Solving I(t) DE trajectory using ODE solver
library(deSolve)
library(mvtnorm)
library(MASS)
library(coda)
library(tictoc)


# initial (state) values for SIR model
N <- 1000
z0 <- c(S = N-6, I = 6, R = 0)

# vector of time steps
times <- 0:80

# vector of parameters used in SIR model
params <- c(beta = 0.3, gamma = 0.1)

SIR <- function(t, z, params) {
  with(as.list(c(params, z)), {
    dS <- -beta * S * I / N 
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I 
    list(c(dS, dI, dR))
  })
}

r <-  rk4(z0, times, SIR, params)
plot(r)

I_values <- r[, "I"]

data <- rnbinom(81, I_values, 0.54)
plot(times, data, main = 'Observations', xlab = 'time')

library(MCMCpack)
library(truncnorm)
library(numDeriv)

# Define the target distribution as a function of all unknown parameters
td <- function(x) {
  beta <- x[1] # transmission rate
  gamma <- x[2] # recovery rate
  phi <- x[3] # dispersion parameter of the negative binomial model
  S0 <- x[4] # initial number of susceptible individuals
  I0 <- x[5] # initial number of infected individuals
  R0 <- x[6] # initial number of recovered individuals
  
  x <- c(beta, gamma, phi, S0, I0, R0)
  
  params2 <- c(beta, gamma)
  
  SIR_model <- function(t, z, params2) {
    with(as.list(c(params2, z)), {
      dS <- -beta * S * I / N 
      dI <- beta * S * I / N - gamma * I
      dR <- gamma * I 
      list(c(dS, dI, dR))
    })
  }
  
  z_0 <- c(S = S0, I = I0, R = R0)
  out <-  rk4(z_0, times, SIR_model, params2)
  
  I_values <- out[, "I"]
  
  # Calculate the log-likelihood
  loglikelihood <- sum(dnbinom(data, I_values, phi, log = TRUE))
  
  #print(rbind(data, I_values))
  log_prior_beta <- dunif(beta, min = 0, max = 1, log = TRUE)
  log_prior_gamma <- dunif(gamma, min = 0, max = 1, log = TRUE)
  log_prior_phi <- dunif(phi, min = 0, max = 1, log = TRUE)
  log_prior_S0 <- dnorm(S0, 994, 1, log = TRUE)
  log_prior_I0 <- dnorm(I0, 6, 1, log = TRUE)
  log_prior_R0 <- dnorm(R0, 0, 0.01, log = TRUE)
  
  # Return the product of the mixture and prior distributions
  
  return(loglikelihood + log_prior_beta + log_prior_gamma + 
           log_prior_phi  + log_prior_S0  + log_prior_I0 + log_prior_R0)
}

n_iter <- 7500

k <- 5
d <- 6
# current parameters
xt <- c(0.3, 0.1, 0.5, 994, 6, 0)

# Initialize empty matrix to store samples
samples_matrix <- matrix(nrow = n_iter, ncol = d)

tic()
# Run the Metropolis-Hastings algorithm
for (i in 1:n_iter) {
  
  # Propose sets of new values
  y <- matrix(NA, k, d)
  for (j in 1:k){
    y[j,] <- c(rnorm(1, xt[1],0.01), rnorm(1, xt[2],0.01), rnorm(1, xt[3], 0.01), 
               rtruncnorm(1, 0, 1000, xt[4], 1), rtruncnorm(1, 0, 1000, xt[5], 1),
               rtruncnorm(1, 0, 1000, xt[6], 0.1))
  }
  
  # Compute the weights w(yj, x_t) for each trial proposal yj
  logp <- c()
  for (j in 1:k){
    logp[j] <- td(y[j,])
  }
  
  p <- exp(logp)
  
  # Step 2: Select Y among the trial set with probability proportional to weights
  selected_index <- sample(1:k, 1, prob = p)  # Select index based on weights
  
  Y <- y[selected_index,]  # Select Y based on the selected index
  
  # Draw x1*, ..., x(k-1)* from the distribution T(Y, .)
  X_star <- matrix(NA, k-1, d)
  
  for (j in 1:k-1){
    X_star[j,] <- c(rnorm(1, Y[1],0.01), rnorm(1, Y[2],0.01), rnorm(1, Y[3], 0.01), 
                    rtruncnorm(1, 0, 1000, Y[4], 1), rtruncnorm(1, 0, 1000, Y[5], 1), 
                    rtruncnorm(1, 0,1000 , Y[6], 0.1))
  }
  
  x_star <- rbind(X_star, xt)
  
  # Compute the weights w(x_starj,y) for each x_starj
  logpx <- c()
  for (j in 1:k){
    logpx[j] <- td(x_star[j,]) 
  }
  
  px <- exp(logpx)
  
  # calculating the acceptance probability
  a <- sum(p) / sum(px)
  rg <- min(1, a)
  
  
  # Decide whether to accept or reject the proposal
  if (!is.na(rg) && runif(1) < rg) {
    xt <- Y
  }
  
  # Store the current parameters
  samples_matrix[i, ] <- xt
  
} 
toc()
rg

sigma <- cov(samples_matrix)

samples_beta <- data.frame(Iterations = 1: n_iter, transmission_rate = samples_matrix[,1])
samples_gamma <- data.frame(Iterations = 1: n_iter, recovery_rate = samples_matrix[,2])
samples_phi <- data.frame(Iterations = 1: n_iter, disp_para = samples_matrix[,3])
samples_S0 <- data.frame(Iterations = 1 : n_iter, init_S = samples_matrix[,4])
samples_I0 <- data.frame(Iterations = 1 : n_iter, init_I = samples_matrix[,5])
samples_R0 <- data.frame(Iterations = 1 : n_iter, init_R = samples_matrix[,6])

par(mfrow = c(2, 2)) 
plot(samples_beta, type = 'l', main = 'Trace plot of beta')

plot(samples_gamma, type = 'l', main = 'Trace plot of gamma')

plot(samples_phi, type = 'l', main = 'Trace plot of phi')

plot(samples_S0, type = 'l', main = 'Trace plot of S0')

plot(samples_I0, type = 'l', main = 'Trace plot of I0')

plot(samples_R0, type = 'l', main = 'Trace plot of R0')

par(mfrow = c(3, 2)) 

hist(samples_matrix[,1], freq = FALSE, main = 'Histogram of beta', xlab = 'beta')

hist(samples_matrix[,2], freq = FALSE, main = 'Histogram of gamma', xlab = 'gamma')

hist(samples_matrix[,3], freq = FALSE, main = 'Histogram of phi', xlab = 'phi')

hist(samples_matrix[,4], freq = FALSE, main = 'Histogram of S0', xlab = 'S0')

hist(samples_matrix[,5], freq = FALSE, main = 'Histogram of I0', xlab = 'I0')

hist(samples_matrix[,6], freq = FALSE, main = 'Histogram of R0', xlab = 'R0')

(beta_approx <- mean(samples_beta$transmission_rate))
(gamma_approx <- mean(samples_gamma$recovery_rate))
(phi_approx <- mean(samples_phi$disp_para))
(S0_approx <- mean(samples_S0$init_S))
(I0_approx <- mean(samples_I0$init_I))
(R0_approx <- mean(samples_R0$init_R))
(R <- beta_approx/gamma_approx)

(ess1 <- effectiveSize(samples_matrix))

# GMTM
Tp1 <- function(x, y){
  b * dmvnorm(x, y, 2.38^2 * sigma/d, log = TRUE) + 
    (1 - b) * dmvnorm(x, y, 0.1^2 * diag(d)/d, log = TRUE)
}

s <- function(x) {
  grad(td, x)
}

D <- function(x){
  hessian(td, x)
}

n_iter1 <- 1500

b <- 0.988

# Initialize empty matrix to store samples
samples_matrix1 <- matrix(nrow = n_iter1, ncol = d)

tic()
# Run the Metropolis-Hastings algorithm
for (i in 1:n_iter1) {
  
  # Propose sets of new values
  y1 <- matrix(NA, k, d)
  for (j in 1:k){
    y1[j, ] <- b * rmvnorm(1, xt, 2.38^2 * sigma/d) +
      (1 - b) * rmvnorm(1, xt, 0.1^2 * diag(d)/d)
    }
  
  # Compute the value of A(y,x_t)
  #s1 <- s(xt)
  #D1 <- D(xt)
  logA1 <- c()
  #A1 <- c()
  for (j in 1:k){
    #logA1[j] <- t(s1) %*% (y1[j,]- xt) + (0.5 * t(y1[j, ]- xt) %*% D1 %*% (y1[j, ] - xt))
    logA1[j] <- t(s(xt)) %*% (y1[j,]- xt) + (0.5 * t(y1[j, ]- xt) %*% D(xt) %*% (y1[j, ] - xt))
    #A1[j] <- exp(logA1[j])
  }
  
  logpy1 <- c()
  logT1 <- c()
  for (j in 1:k){
    logT1[j] <- Tp1(xt, y1[j, ])
    #logpy1[j] <- logA1[j] - logT1[j] - log(sum(A1[j])/sum(exp(logT1[j])))
    logpy1[j] <- logA1[j] - logT1[j] - log(sum(exp(logA1[j]))/sum(exp(logT1[j])))
  }

  # Step 2: Select Y among the trial set with probability proportional to weights
  selected_index1 <- sample(1:k, 1, prob = exp(logpy1))  # Select index based on weights
  
  Y1 <- y1[selected_index1,]  # Select Y based on the selected index
  
  # Draw x1*, ..., x(k-1)* from the distribution T(Y, .)
  X_star1 <- matrix(NA, k-1, 6)
  
  for (j in 1:k-1){
    X_star1[j, ] <- b * rmvnorm(1, Y1, 2.38^2 * sigma/d) +
      (1 - b) * rmvnorm(1, Y1, 0.1^2 * diag(d)/d)
  }
  
  x_star1 <- rbind(X_star1, xt)
  
  #s2 <- s(Y1)
  #D2 <- D(Y1)
  logA2 <- c()
  #A2 <- c()
  for (j in 1:k){
    logA2[j] <- t(s(Y1)) %*% (x_star1[j,]- Y1) + (0.5 * t(x_star1[j, ]- Y1) %*% D(Y1) %*% (x_star1[j, ] - Y1))
    #A2[j] <- exp(logA2[j])
  }
  
  logpxt1 <- c()
  logT2 <- c()
  for (j in 1:k){
    logT2[j] <- Tp1(Y1, x_star1[j, ])
    logpxt1[j] <- logA2[j] - logT2[j] - log(sum(exp(logA2[j]))/sum(exp(logT2[j])))
  }
  
  # calculating the acceptance probability
  a1 <- exp(td(Y1) + Tp1(Y1, xt) + logpxt1 - td(xt) - Tp1(xt, Y1) - logpy1)
  
  rg1 <- min(1, a1)
  # Decide whether to accept or reject the proposal
  if (!is.na(rg1) && runif(1) < rg1) {
    xt <- Y1
  }
  
  # Store the current parameters
  samples_matrix1[i, ] <- xt
} 

toc()
rg1

samples_beta1 <- data.frame(Iterations = 1: n_iter1, transmission_rate = samples_matrix1[,1])
samples_gamma1 <- data.frame(Iterations = 1: n_iter1, recovery_rate = samples_matrix1[,2])
samples_phi1 <- data.frame(Iterations = 1: n_iter1, disp_para = samples_matrix1[,3])
samples_S01 <- data.frame(Iterations = 1 : n_iter1, init_S = samples_matrix1[,4])
samples_I01 <- data.frame(Iterations = 1 : n_iter1, init_I = samples_matrix1[,5])
samples_R01 <- data.frame(Iterations = 1 : n_iter1, init_R = samples_matrix1[,6])

par(mfrow = c(2, 2)) 
plot(samples_beta1, type = 'l', main = 'Trace plot of beta', col = 'darkred')

plot(samples_gamma1, type = 'l', main = 'Trace plot of gamma', col = 'darkred')

plot(samples_phi1, type = 'l', main = 'Trace plot of phi', col = 'darkred')

plot(samples_S01, type = 'l', main = 'Trace plot of S0', col = 'darkred')

plot(samples_I01, type = 'l', main = 'Trace plot of I0', col = 'darkred')

plot(samples_R01, type = 'l', main = 'Trace plot of R0', col = 'darkred')

par(mfrow = c(3, 2)) 

hist(samples_matrix1[,1], freq = FALSE, main = 'Histogram of beta', xlab = 'beta', col = "pink")

hist(samples_matrix1[,2], freq = FALSE, main = 'Histogram of gamma', xlab = 'gamma', col = "pink")

hist(samples_matrix1[,3], freq = FALSE, main = 'Histogram of phi', xlab = 'phi', col = "pink")

hist(samples_matrix1[,4], freq = FALSE, main = 'Histogram of S0', xlab = 'S0', col = "pink")

hist(samples_matrix1[,5], freq = FALSE, main = 'Histogram of I0', xlab = 'I0', col = "pink")

hist(samples_matrix1[,6], freq = FALSE, main = 'Histogram of R0', xlab = 'R0', col = "pink")

(beta1_approx <- mean(samples_beta1$transmission_rate))
(gamma1_approx <- mean(samples_gamma1$recovery_rate))
(phi1_approx <- mean(samples_phi1$disp_para))
(S01_approx <- mean(samples_S01$init_S))
(I01_approx <- mean(samples_I01$init_I))
(R01_approx <- mean(samples_R01$init_R))
(R1 <- beta1_approx/gamma1_approx)

(ess2 <- effectiveSize(samples_matrix1))