# Simulations to test the effect of changing the semi-metrics (for the GCV and AE methods)
# Using the regression operator 3 and covariables X_S
# Results correspond to Table 4.6

# Loading packages
library(splines)

# Load functions and routines needed to conduct the simulation study
# Please save the scripts below in the same folder as the simulation scripts
source("auxiliary_functions.R")
source("functional_kernel_regression_GCV.R")
source("functional_kernel_regression_LCV.R")
source("functional_kernel_regression_adaptive_estimation.R")

# Function for creating evaluations of the smooth covariable
curve_S <- function(t_j, a_i, b_i, c_i, noise = 0)
{
  x_tj <-
    a_i * cos(2 * t_j) + b_i * sin(4 * t_j) + c_i * (t_j ^ 2 - pi * t_j + (2 /
                                                                             9) * pi ^ 2)
  return(x_tj)
}


# Set random seed for reproducibility purposes
set.seed(121)

# Set number of simulation iterations
n_study <- 50

# Fix sample size and evaluation points
n <- 150
points <- 100
# Create evaluation grid
t <- seq(0, pi, length = points)

# Initialising vectors
GCV_s3_L2_ASPE <- c()
GCV_s3_deriv_ASPE <- c()
GCV_s3_pca_ASPE <- c()
AE_s3_L2_ASPE <- c()
AE_s3_deriv_ASPE <- c()
AE_s3_pca_ASPE <- c()
integral3 <- c()

# Set training and test set threshold
threshold <- 0.7 * n
# Training set index
train_index <- 1:ceiling(threshold)
# Test set index
test_index <- (ceiling(threshold) + 1):n

# Iteration index
index <- 1

# Repeat the simulation for n_study times
while(index <= n_study){
  a <- runif(n, 0, 1)
  b <- runif(n, 0, 1)
  c <- runif(n, 0, 1)
  mat <- matrix(NA, nrow = n, ncol = points)
  
  # Creating curves X_s with regression function III
  for (i in 1:n) {
    for (j in 1:points) {
      mat[i, j] <- curve_S(t[j], a[i], b[i], c[i])
    }
    fun3 <- function(t) {log(1+abs(a[i]*cos(2*t) + b[i] * sin(4*t) + c[i]*(t^2 - pi*t + (2/9)*pi^2)))}
    integral3[i] <- integrate(fun3, 0, pi, subdivisions=250)$value
  }
  rs3 <- array(as.numeric(unlist(integral3)), dim=c(1,n,1))
  errors <- rnorm(n, 0, 1)
  y_s3 <- rs3 + errors
  
  # Using the index vectors to divide the data into training and test set
  train_curves <- mat[train_index,]
  y_s3_train <- y_s3[train_index]
  test_curves <- mat[test_index,]
  y_s3_test <- y_s3[test_index]
  
  # Performing functional kernel regression (with global cross validation)
  # using the L2-norm
  GCV_result_s3_L2 <-
    functional_kernel_regression_GCV(
      y_s3_train,
      train_curves,
      test_curves,
      kernel = 'triangle',
      d = 'L2'
    )
  
  # Performing functional kernel regression (with global cross validation)
  # using the semi-metric based on derivatives
  GCV_result_s3_deriv <-
    functional_kernel_regression_GCV(
      y_s3_train,
      train_curves,
      test_curves,
      kernel = 'triangle',
      d = 'deriv',
      q = 2,
      nknot = 20,
      c(0, pi)
    )
  
  # Performing functional kernel regression (with global cross validation)
  # using the semi-metric based on principal components
  GCV_result_s3_pca <-
    functional_kernel_regression_GCV(
      y_s3_train,
      train_curves,
      test_curves,
      kernel = 'triangle',
      d = 'pca',
      q = 2
    )
  
  # Performing functional kernel regression with adaptively estimated bandwidth
  # using the L2-norm
  AE_result_s3_L2 <-
    adaptive_functional_regression(
      y_s3_train,
      train_curves,
      test_curves,
      kernel = "triangle",
      d = "L2"
    )
  
  # Performing functional kernel regression with adaptively estimated bandwidth
  # using the semi-metric based on derivatives
  AE_result_s3_deriv <-
    adaptive_functional_regression(
      y_s3_train,
      train_curves,
      test_curves,
      kernel = "triangle",
      d = "deriv",
      q = 2,
      nknot = 20,
      c(0,pi)
    )
  
  # Performing functional kernel regression with adaptively estimated bandwidth
  # using the semi-metric based on principal components
  AE_result_s3_pca <-
    adaptive_functional_regression(
      y_s3_train,
      train_curves,
      test_curves,
      kernel = "triangle",
      d = "pca",
      q = 2
    )
  
  # Storing the ASPE in these vectors for each iteration
  GCV_s3_L2_ASPE[index] <- pred_error(y_s3_test, GCV_result_s3_L2$y_predicted)
  GCV_s3_deriv_ASPE[index] <- pred_error(y_s3_test, GCV_result_s3_deriv$y_predicted)
  GCV_s3_pca_ASPE[index] <- pred_error(y_s3_test, GCV_result_s3_pca$y_predicted)
  AE_s3_L2_ASPE[index] <- pred_error(y_s3_test, AE_result_s3_L2$y_predicted)
  AE_s3_deriv_ASPE[index] <- pred_error(y_s3_test, AE_result_s3_deriv$y_predicted)
  AE_s3_pca_ASPE[index] <- pred_error(y_s3_test, AE_result_s3_pca$y_predicted)
  
  index <- index + 1
}

# Average the ASPE over all simulations up to the third decimal point
GCV_s3_L2_ASPE_averaged <- round(sum(GCV_s3_L2_ASPE)/n_study,3)
GCV_s3_deriv_ASPE_averaged <- round(sum(GCV_s3_deriv_ASPE)/n_study,3)
GCV_s3_pca_ASPE_averaged <- round(sum(GCV_s3_pca_ASPE)/n_study,3)
AE_s3_L2_ASPE_averaged <- round(sum(AE_s3_L2_ASPE)/n_study,3)
AE_s3_deriv_ASPE_averaged <- round(sum(AE_s3_deriv_ASPE)/n_study,3)
AE_s3_pca_ASPE_averaged <- round(sum(AE_s3_pca_ASPE)/n_study,3)

print(paste("ASPE for GCV for the L2-metric as well as semi-metrics based on derivatives and principal components: ", GCV_s3_L2_ASPE_averaged,GCV_s3_deriv_ASPE_averaged, GCV_s3_pca_ASPE_averaged))
print(paste("ASPE for AE for the L2-metric as well as semi-metrics based on derivatives and principal components: ", AE_s3_L2_ASPE_averaged, AE_s3_deriv_ASPE_averaged, AE_s3_pca_ASPE_averaged))
