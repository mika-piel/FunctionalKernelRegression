# Script for plotting the covariables from Chapter 4, Figures 4.2, 4.3 and 4.4

# Load necessary packages
library(dplyr)
library(tidyr)
library(ggplot2)

# Load functions and routines needed to conduct the simulation study
# Please save the scripts below in the same folder as the simulation scripts
source("fct_plot.R")

# Function for creating evaluations of X_S
curve_S <- function(t_j, a_i, b_i, c_i, noise = 0)
{
  x_tj <-
    a_i * cos(2 * t_j) + b_i * sin(4 * t_j) + c_i * (t_j ^ 2 - pi * t_j + (2 /
                                                                             9) * pi ^ 2)
  return(x_tj)
}

# Function for creating evaluations of X_R
curve_R <- function(t_j, a_i, b_i, noise = 0){
  x_tj <- a_i * sin(4*(b_i - t_j)) + b_i + noise
  return(x_tj)
}


# Set random seed for reproducibility purposes
set.seed(121)

# Sample size
n <- 150
# Number of evaluation points
points <- 100
# Creating the evaluation grid
t <- seq(0, pi, length = points)

# Random vectors for the X_S
as <- runif(n, 0, 1)
bs <- runif(n, 0, 1)
cs <- runif(n, 0, 1)

# Random vectors for the X_R
tr <- seq(0, 1, length=points)
a1 <- rnorm(ceiling(n/2), -2, 0.5)
a2 <- rnorm(floor(n/2), 3, 2)
ar <- c(a1, a2)
br <- rnorm(n, 0, 3)


# Initialise an empty matrix for storing the discretised curves 
matS <- matrix(NA, nrow = n, ncol = points)
matR <- matrix(NA, nrow = n, ncol = points)

# Computing the functional data (of type X_S and X_R)
for (i in 1:n) {
  for (j in 1:points) {
    noise_gaussian <- rnorm(1, mean=0, sd=0.5)
    matS[i, j] <- curve_S(t[j], as[i], bs[i], cs[i])
    matR[i,j] <- curve_R(tr[j], ar[i], br[i], noise_gaussian)
  }
}
  
# Plot 150 curves of type X_S
dfS <- as.data.frame(matS)
plotS <- fct_plot(dfS)
show(plotS)

# Compute and plot the n Brownian Motions
# Variance of the increments
sig <- 0.01
matL <- matrix(rnorm(n = n * (length(t) - 1), mean = 0, sd = sqrt(sig)), n, length(t) - 
                  1)
matL <- cbind(rep(0, n), t(apply(matL, 1, cumsum)))

dfL <- as.data.frame(matL)
plotB <- fct_plot(dfL, "Brownian Motions")
show(plotB)

# Plot 150 curves of type X_R
dfR <- as.data.frame(matR)
plotR <- fct_plot(dfR)
show(plotR)