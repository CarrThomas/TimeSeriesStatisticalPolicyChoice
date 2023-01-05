# Toy model
rm(list = ls())

library(mvtnorm)
library(invgamma)

# specify model
beta <- 0.96 # discount factor
lambda <- 0.62 # weight on output

theta <- 0.13 # coefficient on pi in the Phillips Curve
rho <- 0.52 # persistance of output deviations
sigma_e <- 1.3 # standard deviation of the shock

# get prior as the regression result from T periods
T_prior <- 30
y <- rep(0, T_prior)
y_lag <- rep(0, T_prior)
pi <- rep(0, T_prior)

for (t in 1:T_prior){
  
  e <- rnorm(1, 0, sigma_e)
  pi[t] <- - runif(1, 0, 1) * y_lag[t] 
  y[t] <- rho * y_lag[t] + theta * pi[t] + e
  
  if (t < T_prior){
    y_lag[t + 1] <- y[t]
  }
  
}

data <- data.frame(y = y, y_lag = y_lag, pi = pi)
prior <- lm(y ~ y_lag + pi - 1, data = data)

P_0 <- t(cbind(y_lag, pi)) %*% cbind(y_lag, pi)
Bhat_0 <- prior$coefficients
s_0 <- sum(prior$residuals ^ 2)
v_0 <- prior$df.residual

# main recursion
T_main <- 10000
y <- rep(0, T_main)
y_lag <- c(y_lag[T_prior], rep(0, T_main - 1))
pi <- rep(0, T_main)

P <- vector(mode ="list", length = T_main)
P[[1]] <- P_0  
Bhat <- vector(mode ="list", length = T_main)
Bhat[[1]] <- Bhat_0
s <- rep(0, T_main)
s[1] <- s_0
v <- rep(0, T_main)
v[1] <- v_0
e <- rep(0, T_main)

for (t in 1:T_main){
  
  # draw parameters from the prior
  sigma_draw <- rinvgamma(1, shape = v[t] / 2, rate = s[t] / 2)
  B_draw <- rmvnorm(1, Bhat[[t]], sigma_draw * solve(P[[t]])) 
  
  # get optimal policy
  rho_draw <- B_draw[1]
  theta_draw <- B_draw[2]
  
  a <- -2 * beta * theta_draw ^ 2
  b <- beta * rho_draw ^ 2 + lambda * beta * theta_draw ^ 2 - 1
  c <- lambda
  
  root_one <- (- b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a)
  root_two <- (- b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a)
  
  xi <- max(root_one, root_two)
  
  coef <- - beta * theta_draw * xi * rho_draw / (1 + beta * xi * theta_draw ^ 2)
  
  # simulate
  # pi[t] <- coef * y_lag[t]
  pi[t] <- rnorm(1)
  e[t] <- rnorm(1, 0, sigma_e)
  y[t] <- rho * y_lag[t] + theta * pi[t] + e[t]
  
  if (t < T_main){
    
    y_lag[t + 1] <- y[t]
    
    # update belief
    P[[t + 1]] <- P[[t]] + 
      matrix(c(y_lag[t], pi[t]), 2, 1) %*% matrix(c(y_lag[t], pi[t]), 1, 2)
    Bhat[[t + 1]] <- solve(P[[t + 1]]) %*% (P[[t]] %*% as.matrix(Bhat[[t]]) + 
      matrix(c(y_lag[t], pi[t]), 2, 1) * y[t])
    s[t + 1] <- s[t] + y[t] ^ 2 + 
      matrix(Bhat[[t]], 1, 2) %*% 
      P[[t]] %*% 
      matrix(Bhat[[t]], 2, 1) - 
      matrix(Bhat[[t + 1]], 1, 2) %*% 
      P[[t + 1]] %*% 
      matrix(Bhat[[t + 1]], 2, 1)
    v[t + 1] <- v[t] + 1
    
  }
}

rho_est <- unlist(Bhat)[seq(from = 1, by = 2, length.out = T_main)]
theta_est <-  unlist(Bhat)[seq(from = 2, by = 2, length.out = T_main)]
sigma_est <- s / (v - 2)

plot(rho_est)
lines(1:T_main, rep(rho, T_main))
plot(theta_est)
lines(1:T_main, rep(theta, T_main))
plot(sigma_est)
lines(1:T_main, rep(sigma_e ^ 2, T_main))

# optimal policy
a <- -2 * beta * theta ^ 2
b <- beta * rho ^ 2 + lambda * beta * theta ^ 2 - 1
c <- lambda

root_one <- (- b + sqrt(b ^ 2 - 4 * a * c)) / (2 * a)
root_two <- (- b - sqrt(b ^ 2 - 4 * a * c)) / (2 * a)

xi <- max(root_one, root_two)

coef <- - beta * theta * xi * rho / (1 + beta * xi * theta ^ 2)

pi_opt <- coef * y_lag

plot(pi, xlab = "t", ylab = expression(pi))
points(pi_opt, col = 2)
legend("topright", bty = "n", legend = c("Thompson Sampling", "Optimal"),  
       pch = c(1, 1), col = c(1, 2))
