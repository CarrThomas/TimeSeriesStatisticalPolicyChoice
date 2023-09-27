################################# FRONT MATTER ################################# 

# Simulation exercise based on Wieland (2000) "Monetary Policy, Parameter 
# Uncertainty and Optimal Learning"

# Simulates a Thompson Sampling central bank and a passive learning/adaptive
# policy/anticipated utility central bank for 100 periods, then builds plots.

rm(list = ls())
gc()

library(ggplot2)
library(tidyr)
library(mvtnorm)

########################## THOMPSON SAMPLING FUNCTION ##########################

run_Thompson <- function(w, pi_star, B0, B1, sigma_e, b_prior, V_prior, e, T){
  
  i_hist <- rep(0, T)
  pi_hist <- rep(0, T)
  b_hist <- matrix(0, 2, T + 1) 
  V_hist <- array(0, dim = c(2, 2, T + 1))
  
  b_hist[, 1] <- b_prior
  V_hist[, , 1] <- V_prior 
  
  b_old <- b_prior
  V_old <- V_prior
  
  for (t in 1:T){
    
    # draw parameter values from the belief
    B_draw <- rmvnorm(1, mean = b_old, sigma = V_old)
    
    # optimal i given the drww
    i <- B_draw[2] * (pi_star - B_draw[1]) / (w + B_draw[2] ^ 2) 
    
    # inflation
    pi <- B0 + B1 * i + e[t]
    
    # update
    V_new <- solve(solve(V_old) + rbind(c(1, i), c(i, i ^ 2)) / (sigma_e ^ 2))
    b_new <- V_new %*% (solve(V_old) %*% b_old + 
                          matrix(c(1, i), 2, 1) * pi / (sigma_e ^ 2))                                                      
    
    # hist results
    i_hist[t] <- i
    pi_hist[t] <- pi
    b_hist[, t + 1] <- b_new
    V_hist[, , t + 1] <- V_new
    
    # next period
    b_old <- b_new
    V_old <- V_new
    
  }
  
  return(list(i = i_hist, pi = pi_hist, b = b_hist, V = V_hist))
  
}

########################## PASSIVE LEARNING FUNCTION ###########################

run_passive <- function(w, pi_star, B0, B1, sigma_e, b_prior, V_prior, e, T){
  
  i_hist <- rep(0, T)
  pi_hist <- rep(0, T)
  b_hist <- matrix(0, 2, T + 1) 
  V_hist <- array(0, dim = c(2, 2, T + 1))
  
  b_hist[, 1] <- b_prior
  V_hist[, , 1] <- V_prior 
  
  b_old <- b_prior
  V_old <- V_prior
  
  for (t in 1:T){
    
    # optimal i given beliefs
    i <- (b_old[2] * (pi_star - b_old[1]) - V_old[1, 2]) / 
      (w + b_old[2] ^ 2 + V_old[2, 2])
    
    # inflation
    pi <- B0 + B1 * i + e[t]
    
    # update
    V_new <- solve(solve(V_old) + rbind(c(1, i), c(i, i ^ 2)) / (sigma_e ^ 2))
    b_new <- V_new %*% (solve(V_old) %*% b_old + 
                          matrix(c(1, i), 2, 1) * pi / (sigma_e ^ 2))                                                      
    
    # store results
    i_hist[t] <- i
    pi_hist[t] <- pi
    b_hist[, t + 1] <- b_new
    V_hist[, , t + 1] <- V_new
    
    # next period
    b_old <- b_new
    V_old <- V_new
    
  }
  
  return(list(i = i_hist, pi = pi_hist, b = b_hist, V = V_hist))
}

################################ SET PARAMETERS ################################

T <- 100 # number of periods

w <- 0.14 # weight on the interest rate
pi_star <- 0 # interest rate target
B0 <- 9 # true intercept
B1 <- - 0.7 # true slope
sigma_e <- 1 # variance of the shocks

# oracle policy
i_star <- B1 * (pi_star - B0) / (w + B1 ^ 2)
opt_pi <- B0 + B1 * i_star

# prior
b_prior <- c(6, -0.8)
V_prior <- rbind(c(13.1, -1.54), 
                 c(-1.54, 0.25))

############################### MANY SIMULATIONS ###############################

set.seed(2147)
J <- 100
Thompson <- vector(mode = "list", length = J)
passive <- vector(mode = "list", length = J)

for (j in 1:J){
  
  # draw a sequence of shocks
  e <- rnorm(T, 0, sigma_e)
  
  # Thompson Sampling
  Thompson[[j]] <- run_Thompson(w, pi_star, B0, B1, sigma_e, b_prior, V_prior,
                                e, T)
  
  # passive learning
  passive[[j]] <- run_passive(w, pi_star, B0, B1, sigma_e, b_prior, V_prior,
                              e, T)
  
}

############################## SPECIFIC RUN PLOTS ##############################

# select a particular run
j <- 5

# Comparison of chosen interest rates (Panel a)
temp_data <- data.frame(t = 1:T, 
                        i_star = rep(i_star, T),
                        Thompson = Thompson[[j]]$i,
                        passive = passive[[j]]$i)

temp_data %>% ggplot() +
  geom_point(aes(x = t, y = Thompson, col = "Thompson"), shape = 19) +
  geom_point(aes(x = t, y = passive, col = "passive"), shape = 15) +
  geom_line(aes(x = t, y = i_star, col = "i star"), linetype = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("Thompson" = "black", 
                                "passive" = "gray",
                                "i star" = "black"), 
                     labels = c("Thompson Sampling",
                                "Anticipated Utility",
                                expression('i'^"*"))) +
  labs(x = "t", y = expression('i'[t])) +
   guides(color = guide_legend(override.aes = list(
    linetype = c("blank", "blank", "dashed"),
    shape = c(19, 15, NA),
    size = c(5, 5, 1)
  ), nrow = 3, byrow = TRUE)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.2),
        legend.text = element_text(size=15),
        legend.text.align = 0)

# Comparison of estimates mean estimates for B0
temp_data <- data.frame(t = 1:(T + 1), 
                        B0 = rep(B0, T + 1),
                        Thompson = Thompson[[j]]$b[1, ],
                        passive = passive[[j]]$b[1, ])

temp_data %>% ggplot() +
  geom_point(aes(x = t, y = Thompson, col = "Thompson")) +
  geom_point(aes(x = t, y = passive, col = "passive")) +
  geom_line(aes(x = t, y = B0, col = "B0"), linetype = "dashed") +
  theme_classic() + 
  scale_color_manual(values = c("Thompson" = "black", 
                                "passive" = "gray",
                                "B0" = "black"), 
                     labels = c("Thompson Sampling",
                                "Passive Learning",
                                expression(beta[0]))) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.2),
        legend.text.align = 0) +
  labs(x = "t", y = "") +
  guides(color = guide_legend(override.aes = list(
    linetype = c("blank", "blank", "dashed"),
    shape = c(19, 15, NA)
  ), nrow = 3, byrow = TRUE))

# Comparison of mean estimates for B1
temp_data <- data.frame(t = 1:(T + 1), 
                        B1 = rep(B1, T + 1),
                        Thompson = Thompson[[j]]$b[2, ],
                        passive = passive[[j]]$b[2, ])

temp_data %>% ggplot() +
  geom_point(aes(x = t, y = Thompson, col = "Thompson")) +
  geom_point(aes(x = t, y = passive, col = "passive")) +
  geom_line(aes(x = t, y = B1, col = "B1"), linetype = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("Thompson" = "black", 
                                "passive" = "gray",
                                "i star" = "black"), 
                     labels = c("Thompson Sampling",
                                "Passive Learning",
                                expression(beta[1]))) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.text.align = 0) +
  labs(x = "t", y = "") +
  guides(color = guide_legend(override.aes = list(
    linetype = c("blank", "blank", "dashed"),
    shape = c(19, 15, NA)
  ), nrow = 3, byrow = TRUE))

# Comparison of element [1, 1] of the covariance matrix (panel B)
temp_data <- data.frame(t = 1:(T + 1), 
                        Thompson = Thompson[[j]]$V[1, 1, ],
                        passive = passive[[j]]$V[1, 1, ])

temp_data %>% ggplot() +
  geom_point(aes(x = t, y = Thompson, col = "Thompson")) +
  geom_point(aes(x = t, y = passive, col = "passive")) +
  theme_classic() +
  scale_color_manual(values = c("Thompson" = "black", 
                                "passive" = "gray"), 
                     labels = c("Thompson Sampling",
                                "Passive Learning")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.9),
        legend.text = element_text(size=15),
        legend.text.align = 0) +
  labs(x = "t", y = "") +
  guides(color = guide_legend(override.aes = list(
    linetype = c("blank", "blank"),
    shape = c(19, 15),
    size = c(5, 5)
  ), nrow = 2, byrow = TRUE))


# Comparison of element [2, 2] of the covariance matrix
temp_data <- data.frame(t = 1:(T + 1), 
                        Thompson = Thompson[[j]]$V[2, 2, ],
                        passive = passive[[j]]$V[2, 2, ])

temp_data %>% ggplot() +
  geom_point(aes(x = t, y = Thompson, col = "Thompson")) +
  geom_point(aes(x = t, y = passive, col = "passive")) +
  theme_classic() +
  scale_color_manual(values = c("Thompson" = "black", 
                                "passive" = "gray"), 
                     labels = c("Thompson Sampling",
                                "Passive Learning")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.text.align = 0) +
  labs(x = "t", y = "") +
  guides(color = guide_legend(override.aes = list(
    linetype = c("blank", "blank"),
    shape = c(19, 15)
  ), nrow = 2, byrow = TRUE))

# Comparison of element [1, 2] of the covariance matrix
temp_data <- data.frame(t = 1:(T + 1), 
                        Thompson = Thompson[[j]]$V[1, 2, ],
                        passive = passive[[j]]$V[1, 2, ])

temp_data %>% ggplot() +
  geom_point(aes(x = t, y = Thompson, col = "Thompson")) +
  geom_point(aes(x = t, y = passive, col = "passive")) +
  theme_classic() +
  scale_color_manual(values = c("Thompson" = "black", 
                                "passive" = "gray"), 
                     labels = c("Thompson Sampling",
                                "Passive Learning")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.2),
        legend.text.align = 0) +
  labs(x = "t", y = "") +
  guides(color = guide_legend(override.aes = list(
    linetype = c("blank", "blank"),
    shape = c(19, 15)
  ), nrow = 2, byrow = TRUE))

########################### SUMMARY STATISTICS PLOTS ###########################

#### Thompson Sampling

# mean and variance of i
Thompson_i <- sapply(sapply(Thompson, "[", "i"), "[")
Thompson_mean_i <- rowMeans(Thompson_i)
Thompson_sd_i <- apply(Thompson_i, 1, sd)

# Plots for summary statistics
temp_data <- data.frame(t = 1:T, 
                        i_star = rep(i_star, T), 
                        mean = Thompson_mean_i, 
                        upper = Thompson_mean_i + Thompson_sd_i, 
                        lower = Thompson_mean_i - Thompson_sd_i)

temp_data %>% ggplot() +
  geom_line(aes(x = t, y = i_star, col = "i star"), lty = "dashed") +
  geom_line(aes(x = t, y = mean, col = "mean i")) +
  geom_ribbon(aes(x = t, ymin = lower, ymax = upper, col = "ci"), fill = "gray", alpha = 0.1, show.legend = F) +
  theme_classic() +
  scale_color_manual(breaks = c("i star", "mean i"),
                     values = c(rep("black", 2), "gray"), 
                     labels = c(expression('i'^"*"), 
                                expression(bar("i")))) +
  scale_fill_manual(breaks = "ci",
                    values = "gray",
                    labels = expression(bar("i") %+-% sigma[i])) +
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.2),
        legend.text = element_text(size=15),
        legend.text.align = 0) +
  labs(x = "t", y = "") +
  guides(color = guide_legend(override.aes = list(
    linetype = c("dashed", "solid"),
    shape = c(NA, NA),
    fill = c(NA, NA)
  ),nrow = 2, byrow = TRUE))

#### Passive Learning

# mean and variance of i
passive_i <- sapply(sapply(passive, "[", "i"), "[")
passive_mean_i <- rowMeans(passive_i)
passive_sd_i <- apply(passive_i, 1, sd)

# Plots for summary statistics
temp_data <- data.frame(t = 1:T, 
                        i_star = rep(i_star, T), 
                        mean = passive_mean_i, 
                        upper = passive_mean_i + passive_sd_i, 
                        lower = passive_mean_i - passive_sd_i)

temp_data %>% ggplot() +
  geom_line(aes(x = t, y = i_star, col = "i star"), lty = "dashed") +
  geom_line(aes(x = t, y = mean, col = "mean i")) +
  geom_ribbon(aes(x = t, ymin = lower, ymax = upper, col = "ci"), fill = "gray", alpha = 0.1, show.legend = F) +
  theme_classic() +
  scale_color_manual(breaks = c("i star", "mean i"),
                     values = c(rep("black", 2), "gray"), 
                     labels = c(expression('i'^"*"), 
                                expression(bar("i")))) +
  scale_fill_manual(breaks = "ci",
                    values = "gray",
                    labels = expression(bar("i") %+-% sigma[i])) +
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.2),
        legend.text = element_text(size=15),
        legend.text.align = 0) +
  labs(x = "t", y = "") +
  guides(color = guide_legend(override.aes = list(
    linetype = c("dashed", "solid"),
    shape = c(NA, NA),
    fill = c(NA, NA)
  ),nrow = 2, byrow = TRUE))

# # mean and variance of b0 and b1, Thompson Sampling
# Thompson_b0 <- sapply(sapply(Thompson, "[", "b"), "[", i = 1, j = 1:(T + 1))
# Thompson_mean_b0 <- rowMeans(Thompson_b0)
# Thompson_sd_b0 <- apply(Thompson_b0, 1, sd)
# 
# Thompson_b1 <- sapply(sapply(Thompson, "[", "b"), "[", i = 2, j = 1:(T + 1))
# Thompson_mean_b1 <- rowMeans(Thompson_b1)
# Thompson_sd_b1 <- apply(Thompson_b1, 1, sd)
# 
# # mean and variance of b0 and b1, passive learning
# passive_b0 <- sapply(sapply(passive, "[", "b"), "[", i = 1, j = 1:(T + 1))
# passive_mean_b0 <- rowMeans(passive_b0)
# passive_sd_b0 <- apply(passive_b0, 1, sd)
# 
# passive_b1 <- sapply(sapply(passive, "[", "b"), "[", i = 2, j = 1:(T + 1))
# passive_mean_b1 <- rowMeans(passive_b1)
# passive_sd_b1 <- apply(passive_b1, 1, sd)