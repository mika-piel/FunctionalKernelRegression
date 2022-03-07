# Simulations to test the effect of changing the kernel (for the GCV and AE methods)
# Using the regression operator III and covariables X_S.
# Results correspond to Table 4.5

# Loading necessary packages
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

# Set random seed
set.seed(121)

# Set number of simulation iterations
n_study <- 50

# Fix sample size and evaluation points
n <- 150
points <- 100
# Create evaluation grid
t <- seq(0, pi, length = points)

# Initialise vectors
GCV_s3_box_ASPE <- c()
GCV_s3_quadratic_ASPE <- c()
GCV_s3_triangle_ASPE <- c()
AE_s3_box_ASPE <- c()
AE_s3_quadratic_ASPE <- c()
AE_s3_triangle_ASPE <- c()
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
  
  # For the smooth curves X_s three standard normal distributed r.v. are needed
  a <- runif(n, 0, 1)
  b <- runif(n, 0, 1)
  c <- runif(n, 0, 1)
  
  # Initialise an empty matrix for storing the discretised curves 
  mat <- matrix(NA, nrow = n, ncol = points)
  
  # Creating curves X_s with regression function III
  for (i in 1:n) {
    for (j in 1:points) {
      mat[i, j] <- curve_S(t[j], a[i], b[i], c[i])
    }
    # Create the random covariable function
    fun3 <- function(t) {log(1+abs(a[i]*cos(2*t) + b[i] * sin(4*t) + c[i]*(t^2 - pi*t + (2/9)*pi^2)))}
    # Calculate the regression operator 3 for the random covariable function fun3
    integral3[i] <- integrate(fun3,0,pi, subdivisions=250)$value
  }
  
  # Store values of intregral 3 as a numerical array
  rs3 <- array(as.numeric(unlist(integral3)), dim=c(1,n,1))
  # Create random normal distributed errors
  errors <- rnorm(n, 0, 1)
  # Create the responses corresponding to the non-parametric regression model 3
  y_s3 <- rs3 + errors
  
  # Using the index vectors to divide the data into training and test set
  train_curves <- mat[train_index,]
  y_s3_train <- y_s3[train_index]
  test_curves <- mat[test_index,]
  y_s3_test <- y_s3[test_index]
  
  # Performing functional kernel regression with global cross-validated bandwidth
  # Using the asymmetrical box kernel
  GCV_result_s3_box <-
    functional_kernel_regression_GCV(
      y_s3_train,
      train_curves,
      test_curves,
      kernel = 'box',
      d = 'deriv',
      q = 2,
      nknot = 20,
      c(0, pi)
    )
  
  # Performing functional kernel regression with global cross-validated bandwidth
  # Using the asymmetrical quadratic kernel
  GCV_result_s3_quadratic <-
    functional_kernel_regression_GCV(
      y_s3_train,
      train_curves,
      test_curves,
      kernel = 'quadratic',
      d = 'deriv',
      q = 2,
      nknot = 20,
      c(0, pi)
    )
  
  # Performing functional kernel regression with global cross-validated bandwidth
  # Using the asymmetrical triangle kernel
  GCV_result_s3_triangle <-
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
  
  # Performing functional kernel regression with adaptively estimated bandwidth
  # Using the asymmetrical box kernel
  AE_result_s3_box <-
    adaptive_functional_regression(
      y_s3_train,
      train_curves,
      test_curves,
      kernel = "box",
      d = "deriv",
      q = 2,
      nknot = 20,
      c(0, pi)
    )
  
  # Performing functional kernel regression with adaptively estimated bandwidth
  # Using the asymmetrical quadratic kernel
  AE_result_s3_quadratic <-
    adaptive_functional_regression(
      y_s3_train,
      train_curves,
      test_curves,
      kernel = "quadratic",
      d = "deriv",
      q = 2,
      nknot = 20,
      c(0, pi)
    )
  
  # Performing functional kernel regression with adaptively estimated bandwidth
  # Using the asymmetrical triangle kernel
  AE_result_s3_triangle <-
    adaptive_functional_regression(
      y_s3_train,
      train_curves,
      test_curves,
      kernel = "triangle",
      d = "deriv",
      q = 2,
      nknot = 20,
      c(0, pi)
    )
  
  # For each repetition of the simulation store the respective averaged squared
  # prediction errors of each method and chosen kernel in these vectors
  GCV_s3_box_ASPE[index] <- pred_error(y_s3_test, GCV_result_s3_box$y_predicted)
  GCV_s3_quadratic_ASPE[index] <- pred_error(y_s3_test, GCV_result_s3_quadratic$y_predicted)
  GCV_s3_triangle_ASPE[index] <- pred_error(y_s3_test, GCV_result_s3_triangle$y_predicted)
  AE_s3_box_ASPE[index] <- pred_error(y_s3_test, AE_result_s3_box$y_predicted)
  AE_s3_quadratic_ASPE[index] <- pred_error(y_s3_test, AE_result_s3_quadratic$y_predicted)
  AE_s3_triangle_ASPE[index] <- pred_error(y_s3_test, AE_result_s3_triangle$y_predicted)
  
  # Increment the while-loop index
  index <- index + 1
}

# Average the ASPE over all simulations up to the third decimal point
GCV_s3_box_ASPE_averaged <- round(sum(GCV_s3_box_ASPE)/n_study,3)
GCV_s3_quadratic_ASPE_averaged <- round(sum(GCV_s3_quadratic_ASPE)/n_study,3)
GCV_s3_triangle_ASPE_averaged <- round(sum(GCV_s3_triangle_ASPE)/n_study,3)
AE_s3_box_ASPE_averaged <- round(sum(AE_s3_box_ASPE)/n_study,3)
AE_s3_quadratic_ASPE_averaged <- round(sum(AE_s3_quadratic_ASPE)/n_study,3)
AE_s3_triangle_ASPE_averaged <- round(sum(AE_s3_triangle_ASPE)/n_study,3)

# Print the averaged ASPEs for the two methods for all three choices of kernel function
print(paste("ASPE for GCV for kernels box, quadratic and triangle: ", GCV_s3_box_ASPE_averaged,GCV_s3_quadratic_ASPE_averaged, GCV_s3_triangle_ASPE_averaged))
print(paste("ASPE for AE for kernels box, quadratic and triangle: ", AE_s3_box_ASPE_averaged, AE_s3_quadratic_ASPE_averaged, AE_s3_triangle_ASPE_averaged))