# Script for the simulation study of covariables X_R and varying sample sizes
# and the two regression operators r_II and r_III
# Results correspond to Table 4.11 and Table 4.12

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
# Number of evaluation points
points <- 100
# Creating the evaluation grid over [0,1]
t <- seq(0, 1, length = points)

# Conduct this simulation for varying sample sizes
for (n in c(50, 150, 300)) {
  print(paste("Conducting the study for sample size", n))
  # Set random seed for reproducibility purposes
  set.seed(121)
  # Initialising vectors
  GCV_r2_ASPE <- c()
  LCVB_r2_ASPE <- c()
  LCVT_r2_ASPE <- c()
  AE_r2_ASPE <- c()
  
  GCV_r3_ASPE <- c()
  LCVB_r3_ASPE <- c()
  LCVT_r3_ASPE <- c()
  AE_r3_ASPE <- c()
  
  integral2 <- c()
  integral3 <- c()
  
  # Set training and test set threshold
  threshold <- 0.7 * n
  # Training set index
  train_index <- 1:ceiling(threshold)
  # Test set index
  test_index <- (ceiling(threshold) + 1):n
  
  # Iteration index
  index <- 1
  
  # Repeat the simulation randomly for n_study times
  while (index <= n_study) {
    # Creating random variables needed for the covariables X_R
    a1 <- rnorm(ceiling(n / 2),-2, 0.5)
    a2 <- rnorm(floor(n / 2), 3, 2)
    ar <- c(a1, a2)
    br <- rnorm(n, 0, 3)
    
    # Initialise an empty matrix for storing the discretised curves
    mat <- matrix(NA, nrow = n, ncol = points)
    
    # Computing curves X_s as well as the their values for each regression function
    # Looping over every curve (i) and every evaluation point (j)
    for (i in 1:n) {
      for (j in 1:points) {
        # Gaussian noise at point j for the i-th curve
        noise_gaussian <- rnorm(1, mean = 0, sd = 0.5)
        mat[i, j] <- curve_R(t[j], ar[i], br[i], noise_gaussian)
      }
      
      # Initialise the random functions as functions
      fun2 <-
        function(t) {
          1 / (1 + (ar[i] * sin(4 * (br[i] - t)) + br[i] + noise_gaussian) ^ 2)
        }
      fun3 <-
        function(t) {
          log(1 + abs(ar[i] * sin(4 * (br[i] - t)) + br[i] + noise_gaussian))
        }
      # Calculating the values of the regression functions
      integral2[i] <- integrate(fun2, 0, pi, subdivisions = 250)$value
      integral3[i] <- integrate(fun3, 0, pi, subdivisions = 250)$value
    }
    # Store values of integral2 and intregral 3 as numerical arrays
    rr2 <- array(as.numeric(unlist(integral2)), dim = c(1, n, 1))
    rr3 <- array(as.numeric(unlist(integral3)), dim = c(1, n, 1))
    
    # Create random normal distributed errors as used in the regression model
    errors <- rnorm(n, 0, 1)
    
    # Create the responses corresponding to the non-parametric regression model
    y_r2 <- rr2 + errors
    y_r3 <- rr3 + errors
    
    # Using the index vectors to divide the data into training and test set
    train_curves <- mat[train_index, ]
    y_r2_train <- y_r2[train_index]
    y_r3_train <- y_r3[train_index]
    test_curves <- mat[test_index, ]
    y_r2_test <- y_r2[test_index]
    y_r3_test <- y_r3[test_index]
    
    # For Responses created with regression function 2
    GCV_result_r2 <-
      functional_kernel_regression_GCV(
        y_r2_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "pca",
        q = 2
      )
    
    LCVB_result_r2 <-
      functional_kernel_regression_LCV(
        y_r2_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "pca",
        q = 2,
        tri_weight = FALSE
      )
    
    LCVT_result_r2 <-
      functional_kernel_regression_LCV(
        y_r2_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "pca",
        q = 2,
        tri_weight = TRUE
      )
    
    
    AE_result_r2 <-
      adaptive_functional_regression(
        y_r2_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "pca",
        q = 2
      )
    
    # For Responses created with regression function 3
    GCV_result_r3 <-
      functional_kernel_regression_GCV(
        y_r3_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "pca",
        q = 2
      )
    LCVB_result_r3 <-
      functional_kernel_regression_LCV(
        y_r3_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "pca",
        q = 2,
        tri_weight = FALSE
      )
    
    LCVT_result_r3 <-
      functional_kernel_regression_LCV(
        y_r3_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "pca",
        q = 2,
        tri_weight = TRUE
      )
    
    AE_result_r3 <-
      adaptive_functional_regression(
        y_r3_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "pca",
        q = 2
      )
    
    # For each repetition of the simulation store the respective Mean Squared Errors
    # of each method and data in a vector
    GCV_r2_ASPE[index] <-
      pred_error(y_r2_test, GCV_result_r2$y_predicted)
    LCVB_r2_ASPE[index] <-
      pred_error(y_r2_test, LCVB_result_r2$y_predicted)
    LCVT_r2_ASPE[index] <-
      pred_error(y_r2_test, LCVT_result_r2$y_predicted)
    AE_r2_ASPE[index] <-
      pred_error(y_r2_test, AE_result_r2$y_predicted)
    
    GCV_r3_ASPE[index] <-
      pred_error(y_r3_test, GCV_result_r3$y_predicted)
    LCVB_r3_ASPE[index] <-
      pred_error(y_r3_test, LCVB_result_r3$y_predicted)
    LCVT_r3_ASPE[index] <-
      pred_error(y_r3_test, LCVT_result_r3$y_predicted)
    AE_r3_ASPE[index] <-
      pred_error(y_r3_test, AE_result_r3$y_predicted)
    
    # Increment the while-loop index
    index <- index + 1
  }
  
  # Average the ASPE over each simulation cycle up to the third decimal point
  GCV_r2_ASPE_averaged <- round(sum(GCV_r2_ASPE) / n_study, 3)
  LCVB_r2_ASPE_averaged <- round(sum(LCVB_r2_ASPE) / n_study, 3)
  LCVT_r2_ASPE_averaged <- round(sum(LCVT_r2_ASPE) / n_study, 3)
  AE_r2_ASPE_averaged <- round(sum(AE_r2_ASPE) / n_study, 3)
  
  GCV_r3_ASPE_averaged <- round(sum(GCV_r3_ASPE) / n_study, 3)
  LCVB_r3_ASPE_averaged <- round(sum(LCVB_r3_ASPE) / n_study, 3)
  LCVT_r3_ASPE_averaged <- round(sum(LCVT_r3_ASPE) / n_study, 3)
  AE_r3_ASPE_averaged <- round(sum(AE_r3_ASPE) / n_study, 3)
  
  # Print the averaged ASPEs for each method as well as for each regression function
  print(
    paste(
      "Average ASPE for the second regression function II: ",
      GCV_r2_ASPE_averaged,
      LCVB_r2_ASPE_averaged,
      LCVT_r2_ASPE_averaged,
      AE_r2_ASPE_averaged
    )
  )
  print(
    paste(
      "Average ASPE for the third regression function III: ",
      GCV_r3_ASPE_averaged,
      LCVB_r3_ASPE_averaged,
      LCVT_r3_ASPE_averaged,
      AE_r3_ASPE_averaged
    )
  )
}