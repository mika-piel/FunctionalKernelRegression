# ASPE Plot for varying sample sizes for the global cross-validation method
# with the third regression operator and covariables X_R
# Particularly large sample sizes and their effect on the performance of the GCV estimator will be tested
# The plot corresponds to Figure 4.8

# Loading necessary packages
library(ggplot2)
library(tidyverse)
library(splines)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# Load functions and routines needed to conduct the simulation study
# Please save the scripts below in the same folder as the simulation scripts
source("auxiliary_functions.R")
source("functional_kernel_regression_GCV.R")

curve_R <- function(t_j, a_i, b_i, noise = 0){
  x_tj <- a_i * sin(4*(b_i - t_j)) + b_i + noise
  return(x_tj)
}

points <- 100
n_study <- 50
t <- seq(0, 1, length=points)

GCV_n_mse <- c()
GCV_averaged <- c()


integral3 <- c()
n_seq <- seq(50, 2000, length.out=40)
for(n in n_seq){
set.seed(121)
index <- 1
threshold <- 0.7 * n
train_index <- 1:ceiling(threshold)
test_index <- (ceiling(threshold) + 1):n

# Repeat that procedure randomly for n_study times - Monte Carlo Study
while(index <= n_study){
  
  a1 <- rnorm(ceiling(n/2), -2, 0.5)
  a2 <- rnorm(floor(n/2), 3, 2)
  ar <- c(a1, a2)
  br <- rnorm(n, 0, 3)
  # Initialise an empty matrix for storing the discretised curves 
  mat <- matrix(NA, nrow = n, ncol = points)
  # Computing curves X_s as well as the their values for each regression function
  # Looping over every curve (i) and every evaluation point (j)
  for (i in 1:n) {
    for (j in 1:points) {
      noise_gaussian <- rnorm(1, mean=0, sd=0.5)
      mat[i, j] <- curve_R(t[j], ar[i], br[i], noise_gaussian)
    }
    # Initialise the random functions as functions
    fun3 <- function(t) {log(1+abs(ar[i]*sin(4*(br[i]-t)) + br[i] + noise_gaussian))}
    # Calculating the values of the regression function
    integral3[i] <- integrate(fun3,0,pi, subdivisions = 250)$value
  }
  # Store values integral 3 as numerical array
  r <- array(as.numeric(unlist(integral3)), dim=c(1,n,1))
  # Create random normal distributed errors
  errors <- rnorm(n, 0, 1)
  # Create the responses corresponding to the non-parametric regression model
  y_r <- r + errors
  
  # Using the index vectors to divide the data into training and test set
  train_curves <- mat[train_index,]
  y_train <- y_r[train_index]
  test_curves <- mat[test_index,]
  y_test <- y_r[test_index]
  
  # Performing functional kernel regression (with global cross validation)
  GCV_result_n <-
    functional_kernel_regression_GCV(
      y_train,
      train_curves,
      test_curves,
      kernel = 'triangle',
      d = 'pca',
      q = 2
    )
  
  
  GCV_n_mse[index] <- pred_error(y_test, GCV_result_n$y_predicted)
  
  index <- index + 1
}

GCV_averaged[which(n_seq==n, arr.ind=TRUE)] <- sum(GCV_n_mse)/n_study
}

print(GCV_averaged)

plot_df <- as.data.frame(cbind(n_seq, GCV_averaged))
plot <- ggplot(data=plot_df, aes(x=n_seq, y=GCV_averaged)) +
  geom_point(show.legend=FALSE) +
  ggtitle("ASPE for varying sample sizes (GCV)") +
  xlab("Sample Size") +
  ylab("ASPE") +
  theme(plot.title = element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
show(plot)
