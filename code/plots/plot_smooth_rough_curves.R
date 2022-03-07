# Plot smooth and rough curves for comparison, see Figure 4.1

# Load necessary packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

# Load functions and routines needed for this script
# Please save the script below in the same folder as this script
source("fct_plot.R")

# Function for the covariables X_S
curve_S <- function(t_j, a_i, b_i, c_i, noise = 0)
{
  x_tj <- a_i*cos(2*t_j) + b_i*sin(4*t_j) + c_i*(t_j^2 - pi*t_j + (2/9)*pi^2)
  return(x_tj)
}

# Function for the covariables X_R
curve_R <- function(t_j, a_i, b_i, noise = 0){
  x_tj <- a_i * sin(4*(b_i - t_j)) + b_i + noise
  return(x_tj)
}

# Set the random seed
set.seed(121)

# Set sample size and number of evaluation points
n <- 25
points <- 100

# Initialise empty matrices
mat_S <- matrix(NA, nrow = n, ncol = points)
mat_R <- matrix(NA, nrow=n, ncol = points)

# Evaluation grid and random vectors for the X_S
ts <- seq(0,pi,length=points)
as <- runif(n,0,1)
bs <- runif(n,0,1)
cs <- runif(n,0,1)

# Evaluation grid and random vectors for the X_R
tr <- seq(0, 1, length=points)
a1 <- rnorm(ceiling(n/2), -2, 0.5)
a2 <- rnorm(floor(n/2), 3, 2)
ar <- c(a1, a2)
br <- rnorm(n, 0, 3)

# Compute the curves
for (i in 1:n){
  for(j in 1:points) {
    noise_gaussian <- rnorm(1, mean=0, sd=0.5)
    mat_S[i,j] <- curve_S(ts[j], as[i], bs[i], cs[i])
    mat_R[i,j] <- curve_R(tr[j], ar[i], br[i], noise_gaussian)
  } 
}

# Save the data in a data.frame for ggplot
df_S <- as.data.frame(mat_S)
df_R <- as.data.frame(mat_R)
                      
# Plot the functional data with fct_plot
plot1 <- fct_plot(df_S, "S Curves")
plot2 <- fct_plot(df_R, "R Curves")

# Show plots next to each other
grid.arrange(plot1, plot2, ncol=2)
