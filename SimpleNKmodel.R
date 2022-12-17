library(ggplot2)
library(tidyr)
library(ggpubr)

# Shock responses in two cases

lambda <- 2.11
kappa <- 1.63
sigma <- 3.12
beta <- 0.98

T <- 20

u <- rep(0, T)
u[2] <- 1


# Solution with discretion
pi <- lambda * u / (lambda + kappa ^ 2)
y <- - kappa * u / (lambda + kappa ^ 2)
i <- - y / sigma

discretion <- list(pi = pi, y = y, i = i)

# Solution with commitment
phi <- 2 / ((1 + beta + kappa ^ 2 / lambda) + 
              sqrt((1 + beta + kappa ^ 2 / lambda) ^ 2 - 4 * beta))
for (t in 2:T){
  
  pi[t] <- phi * (pi[t-1] + u[t] - u[t - 1])
  y[t] <- phi * y[t - 1] - kappa * phi * u[t] / lambda
  i[t] <- (phi - 1) * y[t] / sigma + phi * pi[t] - phi * u[t]
  
  
} 

commitment <- list(pi = pi, y = y, i = i)

# Figures

# pi
make_plot <- function(T, discretion, commitment, title){
  
  data.frame(t = 0:(T-1), discretion = discretion, commitment = commitment) %>% 
    ggplot() +
    geom_line(aes(x = t, y = discretion, col = "discretion")) + 
    geom_line(aes(x = t, y = commitment, col = "commitment")) + 
    labs(x = "t", y = "", title = title) +  
    theme_bw() + theme(legend.position = "none")
}

pi_plot <- make_plot(T, discretion$pi, commitment$pi, "Inflation")
y_plot <- make_plot(T, discretion$y, commitment$y, "Output Gap")
i_plot <- make_plot(T, discretion$i, commitment$i, "Interest Rate")

ggarrange(plotlist = list(pi_plot, y_plot, i_plot), nrow = 2, ncol = 2)

