# Script for the simulation study of covariables X_S and varying sample sizes
# and the three regression operators r_I, r_II and r_III
# Corresponding tables are Table 4.2, Table 4.3 and Table 4.4


# Loading necessary packages
library(splines)

# Load functions and routines needed to conduct the simulation study
# Please save the scripts below in the same folder as the simulation scripts
source("auxiliary_functions.R")
source("functional_kernel_regression_GCV.R")
source("functional_kernel_regression_LCV.R")
source("functional_kernel_regression_adaptive_estimation.R")

# Function for creating discretised curves (evaluations) of type X_S
curve_S <- function(t_j, a_i, b_i, c_i)
{
  x_tj <-
    a_i * cos(2 * t_j) + b_i * sin(4 * t_j) + c_i * (t_j ^ 2 - pi * t_j + (2 /
                                                                             9) * pi ^ 2)
  return(x_tj)
}

# Number of iterations of the study
n_study <- 50
# Number of evaluation points
points <- 100
# Creating the evaluation grid over [0,pi]
t <- seq(0, pi, length = points)

# For sample sizes 50, 150 and 300 the simulation study will be conducted
for (n in c(50, 150, 300)) {
  print(paste("Conducting the study for sample size", n))
  # Set random seed for reproducibility purposes
  set.seed(121)
  
  # Initialise vectors needed for computations
  GCV_s1_ASPE <- c()
  LCVB_s1_ASPE <- c()
  LCVT_s1_ASPE <- c()
  AE_s1_ASPE <- c()
  
  GCV_s2_ASPE <- c()
  LCVB_s2_ASPE <- c()
  LCVT_s2_ASPE <- c()
  AE_s2_ASPE <- c()
  
  GCV_s3_ASPE <- c()
  LCVB_s3_ASPE <- c()
  LCVT_s3_ASPE <- c()
  AE_s3_ASPE <- c()
  
  rs1 <- c()
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
  
  # Repeat the following procedure for n_study times
  while (index <= n_study) {
    # For the smooth curves X_s three standard normal distributed r.v. are needed
    a <- runif(n, 0, 1)
    b <- runif(n, 0, 1)
    c <- runif(n, 0, 1)
    
    # Initialise an empty matrix for storing the discretised curves
    mat <- matrix(NA, nrow = n, ncol = points)
    
    # Computing curves X_s as well as the their values for each regression function
    # Looping over every curve (i) and every evaluation point (j)
    for (i in 1:n) {
      for (j in 1:points) {
        # This quantity will be contains the curves (discretised)
        mat[i, j] <- curve_S(t[j], a[i], b[i], c[i])
      }
      
      # Define each random curve as a function of t in R
      fun2 <-
        function(t) {
          1 / (1 + (a[i] * cos(2 * t) + b[i] * sin(4 * t) + c[i] * (t ^ 2 - pi * t + (2 /
                                                                                        9) * pi ^ 2) ^ 2) ^ 2)
        }
      fun3 <-
        function(t) {
          log(1 + abs(a[i] * cos(2 * t) + b[i] * sin(4 * t) + c[i] * (t ^ 2 - pi *
                                                                        t + (2 / 9) * pi ^ 2)))
        }
      
      # Calculating the values of the regression operators
      rs1[i] <- 10 * (a[i] ^ 2 - b[i] ^ 2) + c[i]
      integral2[i] <-
        integrate(fun2, 0, pi, subdivisions = 250)$value
      integral3[i] <-
        integrate(fun3, 0, pi, subdivisions = 250)$value
    }
    
    # Store values of integral2 and integral 3 as numerical arrays
    rs2 <- array(as.numeric(unlist(integral2)), dim = c(1, n, 1))
    rs3 <- array(as.numeric(unlist(integral3)), dim = c(1, n, 1))
    
    # Create random normal distributed errors as used in the regression model
    errors <- rnorm(n, 0, 1)
    
    # Create the responses corresponding to the non-parametric regression model
    y_s1 <- rs1 + errors
    y_s2 <- rs2 + errors
    y_s3 <- rs3 + errors
    
    # Using the index vectors to divide the data into training and test set
    train_curves <- mat[train_index,]
    y_s1_train <- y_s1[train_index]
    y_s2_train <- y_s2[train_index]
    y_s3_train <- y_s3[train_index]
    test_curves <- mat[test_index,]
    y_s1_test <- y_s1[test_index]
    y_s2_test <- y_s2[test_index]
    y_s3_test <- y_s3[test_index]
    
    # For responses created with regression function I
    # Performing functional kernel regression with global cross-validated bandwidth
    GCV_result_s1 <-
      functional_kernel_regression_GCV(
        y_s1_train,
        train_curves,
        test_curves,
        kernel = 'triangle',
        d = 'deriv',
        q = 2,
        nknot = 20,
        c(0, pi)
      )
    
    # Performing functional kernel regression with local cross-validated bandwidth
    # box weights
    LCVB_result_s1 <-
      functional_kernel_regression_LCV(
        y_s1_train,
        train_curves,
        test_curves,
        kernel = 'triangle',
        d = 'deriv',
        q = 2,
        nknot = 20,
        c(0, pi),
        tri_weight = FALSE
      )
    
    # Performing functional kernel regression with Local cross-validated bandwidth
    # with triangle weights
    LCVT_result_s1 <-
      functional_kernel_regression_LCV(
        y_s1_train,
        train_curves,
        test_curves,
        kernel = 'triangle',
        d = 'deriv',
        q = 2,
        nknot = 20,
        c(0, pi),
        tri_weight = TRUE
      )
    
    # Performing functional kernel regression with adaptively estimated bandwidth
    AE_result_s1 <-
      adaptive_functional_regression(
        y_s1_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "deriv",
        q = 2,
        nknot = 20,
        c(0, pi)
      )
    
    # For responses created with regression function II
    # Performing functional kernel regression with global cross-validated bandwidth
    GCV_result_s2 <-
      functional_kernel_regression_GCV(
        y_s2_train,
        train_curves,
        test_curves,
        kernel = 'triangle',
        d = 'deriv',
        q = 2,
        nknot = 20,
        c(0, pi)
      )
    
    # Performing functional kernel regression with local cross-validated bandwidth
    # box weights
    LCVB_result_s2 <-
      functional_kernel_regression_LCV(
        y_s2_train,
        train_curves,
        test_curves,
        kernel = 'triangle',
        d = 'deriv',
        q = 2,
        nknot = 20,
        c(0, pi),
        tri_weight = FALSE
      )
    
    # Performing functional kernel regression with Local cross-validated bandwidth
    # with triangle weights
    LCVT_result_s2 <-
      functional_kernel_regression_LCV(
        y_s2_train,
        train_curves,
        test_curves,
        kernel = 'triangle',
        d = 'deriv',
        q = 2,
        nknot = 20,
        c(0, pi),
        tri_weight = TRUE
      )
    
    # Performing functional kernel regression with adaptively estimated bandwidth
    AE_result_s2 <-
      adaptive_functional_regression(
        y_s2_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "deriv",
        q = 2,
        nknot = 20,
        c(0, pi)
      )
    
    # For responses created with regression function III
    # Performing functional kernel regression with global cross-validated bandwidth
    GCV_result_s3 <-
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
    
    # Performing functional kernel regression with local cross-validated bandwidth
    # box weights
    LCVB_result_s3 <-
      functional_kernel_regression_LCV(
        y_s3_train,
        train_curves,
        test_curves,
        kernel = 'triangle',
        d = 'deriv',
        q = 2,
        nknot = 20,
        c(0, pi),
        tri_weight = FALSE
      )
    
    # Performing functional kernel regression with Local cross-validated bandwidth
    # with triangle weights
    LCVT_result_s3 <-
      functional_kernel_regression_LCV(
        y_s3_train,
        train_curves,
        test_curves,
        kernel = 'triangle',
        d = 'deriv',
        q = 2,
        nknot = 20,
        c(0, pi),
        tri_weight = TRUE
      )
    
    # Performing functional kernel regression with adaptively estimated bandwidth
    AE_result_s3 <-
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
    
    # For each repetition of the simulation store the respective ASPEs
    # of each method and data in a vector
    GCV_s1_ASPE[index] <-
      pred_error(y_s1_test, GCV_result_s1$y_predicted)
    LCVB_s1_ASPE[index] <-
      pred_error(y_s1_test, LCVB_result_s1$y_predicted)
    LCVT_s1_ASPE[index] <-
      pred_error(y_s1_test, LCVT_result_s1$y_predicted)
    AE_s1_ASPE[index] <-
      pred_error(y_s1_test, AE_result_s1$y_predicted)
    
    GCV_s2_ASPE[index] <-
      pred_error(y_s2_test, GCV_result_s2$y_predicted)
    LCVB_s2_ASPE[index] <-
      pred_error(y_s2_test, LCVB_result_s2$y_predicted)
    LCVT_s2_ASPE[index] <-
      pred_error(y_s2_test, LCVT_result_s2$y_predicted)
    AE_s2_ASPE[index] <-
      pred_error(y_s2_test, AE_result_s2$y_predicted)
    
    GCV_s3_ASPE[index] <-
      pred_error(y_s3_test, GCV_result_s3$y_predicted)
    LCVB_s3_ASPE[index] <-
      pred_error(y_s3_test, LCVB_result_s3$y_predicted)
    LCVT_s3_ASPE[index] <-
      pred_error(y_s3_test, LCVT_result_s3$y_predicted)
    AE_s3_ASPE[index] <-
      pred_error(y_s3_test, AE_result_s3$y_predicted)
    
    # Increment the while-loop index
    index <- index + 1
  }
  
  
  # Average the ASPE over all simulations up to the third decimal point
  GCV_s1_ASPE_averaged <- round(sum(GCV_s1_ASPE) / n_study, 3)
  LCVB_s1_ASPE_averaged <- round(sum(LCVB_s1_ASPE) / n_study, 3)
  LCVT_s1_ASPE_averaged <- round(sum(LCVT_s1_ASPE) / n_study, 3)
  AE_s1_ASPE_averaged <- round(sum(AE_s1_ASPE) / n_study, 3)
  
  GCV_s2_ASPE_averaged <- round(sum(GCV_s2_ASPE) / n_study, 3)
  LCVB_s2_ASPE_averaged <- round(sum(LCVB_s2_ASPE) / n_study, 3)
  LCVT_s2_ASPE_averaged <- round(sum(LCVT_s2_ASPE) / n_study, 3)
  AE_s2_ASPE_averaged <- round(sum(AE_s2_ASPE) / n_study, 3)
  
  GCV_s3_ASPE_averaged <- round(sum(GCV_s3_ASPE) / n_study, 3)
  LCVB_s3_ASPE_averaged <- round(sum(LCVB_s3_ASPE) / n_study, 3)
  LCVT_s3_ASPE_averaged <- round(sum(LCVT_s3_ASPE) / n_study, 3)
  AE_s3_ASPE_averaged <- round(sum(AE_s3_ASPE) / n_study, 3)
  
  # Print the averaged ASPEs for each method as well as for each regression operator
  print(
    paste(
      "Average ASPEs for the first regression operator I: ",
      GCV_s1_ASPE_averaged,
      LCVB_s1_ASPE_averaged,
      LCVT_s1_ASPE_averaged,
      AE_s1_ASPE_averaged
    )
  )
  print(
    paste(
      "Average ASPEs for the second regression operator II: ",
      GCV_s2_ASPE_averaged,
      LCVB_s2_ASPE_averaged,
      LCVT_s2_ASPE_averaged,
      AE_s2_ASPE_averaged
    )
  )
  print(
    paste(
      "Average ASPEs for the third regression operator III: ",
      GCV_s3_ASPE_averaged,
      LCVB_s3_ASPE_averaged,
      LCVT_s3_ASPE_averaged,
      AE_s3_ASPE_averaged
    )
  )
}
