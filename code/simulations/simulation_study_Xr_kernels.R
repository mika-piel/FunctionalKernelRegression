# Simulations to test the effect of changing the kernel (for the GCV and AE methods)
# Using the regression function 3 and covariables X_R.
# Results correspond to Table 4.13

# Loading necessary packages
library(splines)

# Load functions and routines needed to conduct the simulation study
# Please save the scripts below in the same folder as the simulation scripts
source("auxiliary_functions.R")
source("functional_kernel_regression_GCV.R")
source("functional_kernel_regression_LCV.R")
source("functional_kernel_regression_adaptive_estimation.R")

# Function for the covariables X_R
curve_R <- function(t_j, a_i, b_i, noise = 0){
  x_tj <- a_i * sin(4*(b_i - t_j)) + b_i + noise
  return(x_tj)
}

# Number of iterations of the study
n_study <- 50
# Set random seed
set.seed(121)

# Fix sample size and evaluation points
n <- 150
points <- 100

# Create evaluation grid
t <- seq(0, 1, length = points)

# Initialising vectors
GCV_r3_box_ASPE <- c()
GCV_r3_quadratic_ASPE <- c()
GCV_r3_triangle_ASPE <- c()
AE_r3_box_ASPE <- c()
AE_r3_quadratic_ASPE <- c()
AE_r3_triangle_ASPE <- c()
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
  a1 <- rnorm(ceiling(n/2), -3, 0.5)
  a2 <- rnorm(floor(n/2), 4, 3)
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
    
    # Initialise the random function as an R function
    fun3 <- function(t) {log(1+abs(ar[i]*sin(4*(br[i]-t)) + br[i] + noise_gaussian))}
    
    # Calculating the values of the regression function
    integral3[i] <- integrate(fun3,0,pi, subdivisions = 250)$value
  }
  # Store values of intregral 3 as a numerical array
  rr3 <- array(as.numeric(unlist(integral3)), dim=c(1,n,1))
  # Create random normal distributed errors
  errors <- rnorm(n, 0, 1)
  # Create the responses corresponding to the non-parametric regression model 3
  y_r3 <- rr3 + errors
  
  # Using the index vectors to divide the data into training and test set
  train_curves <- mat[train_index,]
  y_r3_train <- y_r3[train_index]
  test_curves <- mat[test_index,]
  y_r3_test <- y_r3[test_index]
  
  # Performing functional kernel regression with global cross-validated bandwidth
  # Using the asymmetrical box kernel
  GCV_result_r3_box <-
    functional_kernel_regression_GCV(
      y_r3_train,
      train_curves,
      test_curves,
      kernel = 'box',
      d = 'pca',
      q = 2
    )
  
  # Performing functional kernel regression with global cross-validated bandwidth
  # Using the asymmetrical quadratic kernel
  GCV_result_r3_quadratic <-
    functional_kernel_regression_GCV(
      y_r3_train,
      train_curves,
      test_curves,
      kernel = 'quadratic',
      d = 'pca',
      q = 2
    )
  
  # Performing functional kernel regression with global cross-validated bandwidth
  # Using the asymmetrical triangle kernel
  GCV_result_r3_triangle <-
    functional_kernel_regression_GCV(
      y_r3_train,
      train_curves,
      test_curves,
      kernel = 'triangle',
      d = 'pca',
      q = 2
    )
  
  # Performing functional kernel regression with adaptively estimated bandwidth
  # Using the asymmetrical box kernel
  AE_result_r3_box <-
    adaptive_functional_regression(
      y_r3_train,
      train_curves,
      test_curves,
      kernel = "box",
      d = "pca",
      q = 2
    )
  
  # Performing functional kernel regression with adaptively estimated bandwidth
  # Using the asymmetrical quadratic kernel
  AE_result_r3_quadratic <-
    adaptive_functional_regression(
      y_r3_train,
      train_curves,
      test_curves,
      kernel = "quadratic",
      d = "pca",
      q = 2
    )
  
  # Performing functional kernel regression with adaptively estimated bandwidth
  # Using the asymmetrical triangle kernel
  AE_result_r3_triangle <-
    adaptive_functional_regression(
      y_r3_train,
      train_curves,
      test_curves,
      kernel = "triangle",
      d = "pca",
      q = 2
    )
  
  # For each repetition of the simulation store the respective Mean Squared Errors
  # of each method and chosen kernel
  GCV_r3_box_ASPE[index] <- pred_error(y_r3_test, GCV_result_r3_box$y_predicted)
  GCV_r3_quadratic_ASPE[index] <- pred_error(y_r3_test, GCV_result_r3_quadratic$y_predicted)
  GCV_r3_triangle_ASPE[index] <- pred_error(y_r3_test, GCV_result_r3_triangle$y_predicted)
  AE_r3_box_ASPE[index] <- pred_error(y_r3_test, AE_result_r3_box$y_predicted)
  AE_r3_quadratic_ASPE[index] <- pred_error(y_r3_test, AE_result_r3_quadratic$y_predicted)
  AE_r3_triangle_ASPE[index] <- pred_error(y_r3_test, AE_result_r3_triangle$y_predicted)
  
  # Increment the while-loop index
  index <- index + 1
}

# Average the ASPE over each simulation cycle up to the third decimal point
GCV_r3_box_ASPE_averaged <- round(sum(GCV_r3_box_ASPE)/n_study,3)
GCV_r3_quadratic_ASPE_averaged <- round(sum(GCV_r3_quadratic_ASPE)/n_study,3)
GCV_r3_triangle_ASPE_averaged <- round(sum(GCV_r3_triangle_ASPE)/n_study,3)
AE_r3_box_ASPE_averaged <- round(sum(AE_r3_box_ASPE)/n_study,3)
AE_r3_quadratic_ASPE_averaged <- round(sum(AE_r3_quadratic_ASPE)/n_study,3)
AE_r3_triangle_ASPE_averaged <- round(sum(AE_r3_triangle_ASPE)/n_study,3)

# Print the averaged ASPEs for the two methods for all three choices of kernel function
print(paste("ASPE for GCV for kernels box, quadratic and triangle: ", GCV_r3_box_ASPE_averaged,GCV_r3_quadratic_ASPE_averaged, GCV_r3_triangle_ASPE_averaged))
print(paste("ASPE for AE for kernels box, quadratic and triangle: ", AE_r3_box_ASPE_averaged, AE_r3_quadratic_ASPE_averaged, AE_r3_triangle_ASPE_averaged))
