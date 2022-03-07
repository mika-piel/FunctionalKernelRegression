# Script for plotting estimated small ball probabilities against the bandwidths h
# Results correspond to Figures 4.5, 4.6 and 4.7

# Load necessary packages
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(splines)
library(fda)
library(gridExtra)

# Load functions and routines needed for this script
# Please save the scripts below in the same folder as this script
source("auxiliary_functions.R")

# Functions for creating the covariables X_S and X_R
curve_smooth <- function(t_j, a_i, b_i, c_i, noise = 0){
  x_tj <- a_i*cos(2*t_j) + b_i*sin(4*t_j) + c_i*(t_j^2 - pi*t_j + (2/9)*pi^2)
  return(x_tj)
}

curve_rough <- function(t_j, a_i, b_i, noise = 0){
  x_tj <- a_i * sin(4*(b_i - t_j)) + b_i + noise
  return(x_tj)
}

# Set the random seed
set.seed(121)

# Fix sample size and number of evaluation points
n <- 150
points <- 100

# Initialise empty matrices
mat_smooth <- matrix(NA, nrow = n, ncol = points)
mat_rough <- matrix(NA, nrow=n, ncol = points)

# Set the variance of the Brownian motions' increments
sig <- 0.01
# Evaluation grid
t <- seq(0, 1, length=points)

# Compute n Brownian motions
X <- matrix(rnorm(n = n * (length(t) - 1), mean = 0, sd = sqrt(sig)), n, length(t) - 
              1)
X <- cbind(rep(0, n), t(apply(X, 1, cumsum)))

# Initialise random vectors for the covariables X_S and X_R
ts <- seq(0, pi, length=points)
as <- runif(n,0,1)
bs <- runif(n,0,1)
cs <- runif(n,0,1)

tr <- seq(0, 1, length=points)
a1 <- rnorm(ceiling(n/2), -2, 0.5)
a2 <- rnorm(floor(n/2), 3, 2)
ar <- c(a1, a2)
br <- rnorm(n, 0, 3)

# Compute the curves
for (i in 1:n){
  for(j in 1:points) {
    noise_gaussian <- rnorm(1, mean=0, sd=0.5)
    mat_smooth[i,j] <- curve_smooth(ts[j], as[i], bs[i], cs[i])
    mat_rough[i,j] <- curve_rough(tr[j], ar[i], br[i], noise_gaussian)
  } 
}

# Bandwidth vector
h <- seq(0.0001, 1000, length=500)
# Create semi-metric matrices for each curve set (based on derivatives)
dS_deriv <- d_deriv(mat_smooth, mat_smooth, q=2, nknot=20, c(0,1))
dR_deriv <- d_deriv(mat_rough, mat_rough, q=2, nknot=20, c(0,1))
dB_deriv <- d_deriv(X, X, q=2, nknot=20, c(0,1))

# Initialise empty vectors
sb1_deriv <- c()
sb2_deriv <- c()
sb3_deriv <- c()

# Compute the estimated small ball probabilities
for(i in 1:length(h)){
  sb1_deriv[i] <- sbprob(dS_deriv, 1, h[i])
  sb2_deriv[i] <- sbprob(dB_deriv, 1, h[i])
  sb3_deriv[i] <- sbprob(dR_deriv, 1, h[i])
}

# Save the results in a data.frame to be able to use ggplot
dfS_plot1 <- data.frame(h, sb1_deriv)
dfB_plot1 <- data.frame(h, sb2_deriv)
dfR_plot1 <- data.frame(h, sb3_deriv)

# Plot for each covariable
plotS1 <- ggplot(data=dfS_plot1, aes(x=h, y=sb1_deriv)) +
  geom_line(show.legend=FALSE) +
  xlab("Bandwidths") +
  ylab("Estimated Small Ball Probability") +
  ggtitle("Semi-metric based on derivatives") +
  theme(plot.title = element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plotB1 <- ggplot(data=dfB_plot1, aes(x=h, y=sb2_deriv)) +
  geom_line(show.legend=FALSE) +
  xlab("Bandwidths") +
  ylab("Estimated Small Ball Probability") +
  ggtitle("Semi-metric based on derivatives") +
  theme(plot.title = element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plotR1 <- ggplot(data=dfR_plot1, aes(x=h, y=sb3_deriv)) +
  geom_line(show.legend=FALSE) +
  xlab("Bandwidths") +
  ylab("Estimated Small Ball Probability") +
  ggtitle("Semi-metric based on derivatives") +
  theme(plot.title = element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Follow the same procedure as above
# Create semi-metric matrices for each curve set (based on principal components)
dS_pca <- d_pca(mat_smooth, mat_smooth, q=2)
dR_pca <- d_pca(mat_rough, mat_rough, q=2)
dB_pca <- d_pca(X, X, q=2)

# Initialise vectors
sb1_pca <- c()
sb2_pca <- c()
sb3_pca <- c()

# Compute the small ball probabilities
for(i in 1:length(h)){
  sb1_pca[i] <- sbprob(dS_pca, 1, h[i])
  sb2_pca[i] <- sbprob(dB_pca, 1, h[i])
  sb3_pca[i] <- sbprob(dR_pca, 1, h[i])
}

# Save in a data.frame
dfS_plot2 <- data.frame(h, sb1_pca)
dfB_plot2 <- data.frame(h, sb2_pca)
dfR_plot2 <- data.frame(h, sb3_pca)

# Plot for each covariable
plotS2 <- ggplot(data=dfS_plot2, aes(x=h, y=sb1_pca)) +
  geom_line(show.legend=FALSE) +
  xlab("Bandwidths") +
  ylab("Estimated Small Ball Probability") +
  ggtitle("Semi-metric based on principal components") +
  theme(plot.title = element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plotB2 <- ggplot(data=dfB_plot2, aes(x=h, y=sb2_pca)) +
  geom_line(show.legend=FALSE) +
  xlab("Bandwidths") +
  ylab("Estimated Small Ball Probability") +
  ggtitle("Semi-metric based on principal components") +
  theme(plot.title = element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plotR2 <- ggplot(data=dfR_plot2, aes(x=h, y=sb3_pca)) +
  geom_line(show.legend=FALSE) +
  xlab("Bandwidths") +
  ylab("Estimated Small Ball Probability") +
  ggtitle("Semi-metric based on principal components") +
  theme(plot.title = element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Arrange the respective plots
# For the smooth covariables X_S
grid.arrange(plotS1, plotS2, nrow=2, ncol=1)

# For the Brownian Motions
grid.arrange(plotB1, plotB2, nrow=2, ncol=1)

# For the covariables X_R
grid.arrange(plotR1, plotR2, nrow=2, ncol=1)
