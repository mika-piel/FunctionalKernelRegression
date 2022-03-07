# Script for the simulation study of Brownian Motions as covariables
# and varying sample sizes and the two regression operators r_II and r_III
# Corresponding tables are Table 4.7 and Table 4.8

# Loading necessary packages
library(splines)

# Load functions and routines needed to conduct the simulation study
# Please save the scripts below in the same folder as the simulation scripts
source("auxiliary_functions.R")
source("functional_kernel_regression_GCV.R")
source("functional_kernel_regression_LCV.R")
source("functional_kernel_regression_adaptive_estimation.R")

# Set random seed for reproducibility purposes
set.seed(121)

# Set number of simulation iterations
n_study <- 50

# Number of evaluation points
points <- 100
# Create evaluation grid
t <- seq(0, 1, length = points)

# Set the variance of Brownian Motions
sig <- 0.01

for (n in c(50, 150, 300)) {
  print(paste("Conducting the study for sample size", n))
  
  # Initialise vectors needed for computations
  GCV_l2_ASPE <- c()
  LCVB_l2_ASPE <- c()
  LCVT_l2_ASPE <- c()
  AE_l2_ASPE <- c()
  
  GCV_l3_ASPE <- c()
  LCVB_l3_ASPE <- c()
  LCVT_l3_ASPE <- c()
  AE_l3_ASPE <- c()
  
  integral2 <- c()
  integral3 <- c()
  
  # Set training and test set threshold
  threshold <- 0.7 * n
  # Training set index
  train_index <- 1:ceiling(threshold)
  # Test set index
  test_index <- (ceiling(threshold) + 1):n
  index <- 1
  
  # Repeat the simulations n_study times
  while (index <= n_study){
    # Creating n Brownian Motions
    matXl <-
      matrix(rnorm(
        n = n * (length(t) - 1),
        mean = 0,
        sd = sqrt(sig)
      ), n, length(t) -
        1)
    matXl <- cbind(rep(0, n), t(apply(matXl, 1, cumsum)))
    
    # Compute the integrals occuring in the regression operators II and III
    for (i in 1:n) {
      integral2[i] <-
        integrate(
          approxfun(t, 1 / (1 + matXl[i, ] ^ 2), n = points),
          0,
          1,
          subdivisions = 250,
          stop.on.error = FALSE
        )$value
      integral3[i] <-
        integrate(
          approxfun(t, log(1 + abs(matXl[i, ])), n = points),
          0,
          1,
          subdivisions = 250,
          stop.on.error = FALSE
        )$value
    }
    # Store integral values as numerical arrays
    rl2 <- array(as.numeric(unlist(integral2)), dim = c(1, n, 1))
    rl3 <- array(as.numeric(unlist(integral3)), dim = c(1, n, 1))
    
    # Normal distributed errors
    errors <- rnorm(n, 0, 1)
    
    # Create the responses through the non-parametric regression model
    y_l2 <- rl2 + errors
    y_l3 <- rl3 + errors
    
    # Using the index vectors to divide the data into training and test set
    train_curves <- matXl[train_index, ]
    y_l2_train <- y_l2[train_index]
    y_l3_train <- y_l3[train_index]
    test_curves <- matXl[test_index, ]
    y_l2_test <- y_l2[test_index]
    y_l3_test <- y_l3[test_index]
    
    # Performing functional kernel regression (with global cross validation)
    GCV_result_l2 <-
      functional_kernel_regression_GCV(
        y_l2_train,
        train_curves,
        test_curves,
        kernel = 'triangle',
        d = "deriv",
        q = 2,
        nknot = 20,
        c(0, 1)
      )
    
    # Local cross validation with box weights kernel regression
    LCVB_result_l2 <-
      functional_kernel_regression_LCV(
        y_l2_train,
        train_curves,
        test_curves,
        kernel = 'triangle',
        d = "deriv",
        q = 2,
        nknot = 20,
        c(0, 1),
        tri_weight = FALSE
      )
    
    # Local cross-validated with triangle weights kernel regression
    LCVT_result_l2 <-
      functional_kernel_regression_LCV(
        y_l2_train,
        train_curves,
        test_curves,
        kernel = 'box',
        d = "deriv",
        q = 2,
        nknot = 20,
        c(0, 1),
        tri_weight = TRUE
      )
    
    # Functional kernel regression with adaptively estimated bandwidths
    AE_result_l2 <-
      adaptive_functional_regression(
        y_l2_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "deriv",
        q = 2,
        nknot = 20,
        c(0, 1)
      )
    
    # Performing functional kernel regression (with global cross validation)
    GCV_result_l3 <-
      functional_kernel_regression_GCV(
        y_l3_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "deriv",
        q = 2,
        nknot = 20,
        c(0, 1)
      )
    
    # Local cross validation (with the local cross validation criterion LCV)
    LCVB_result_l3 <-
      functional_kernel_regression_LCV(
        y_l3_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "deriv",
        q = 2,
        nknot = 20,
        c(0, 1),
        tri_weight = FALSE
      )
    
    # Local cross-validated with triangle weights kernel regression
    LCVT_result_l3 <-
      functional_kernel_regression_LCV(
        y_l3_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "deriv",
        q = 2,
        nknot = 20,
        c(0, 1),
        tri_weight = TRUE
      )
    
    # Functional kernel regression with adaptively estimated bandwidths
    AE_result_l3 <-
      adaptive_functional_regression(
        y_l3_train,
        train_curves,
        test_curves,
        kernel = "triangle",
        d = "deriv",
        q = 2,
        nknot = 20,
        c(0, 1)
      )
    
    # Store the averaged squared prediction errors in respective vectors for
    # each simulation iteration
    GCV_l2_ASPE[index] <-
      pred_error(y_l2_test, GCV_result_l2$y_predicted)
    LCVB_l2_ASPE[index] <-
      pred_error(y_l2_test, LCVB_result_l2$y_predicted)
    LCVT_l2_ASPE[index] <-
      pred_error(y_l2_test, LCVT_result_l2$y_predicted)
    AE_l2_ASPE[index] <-
      pred_error(y_l2_test, AE_result_l2$y_predicted)
    
    GCV_l3_ASPE[index] <-
      pred_error(y_l3_test, GCV_result_l3$y_predicted)
    LCVB_l3_ASPE[index] <-
      pred_error(y_l3_test, LCVB_result_l3$y_predicted)
    LCVT_l3_ASPE[index] <-
      pred_error(y_l3_test, LCVT_result_l3$y_predicted)
    AE_l3_ASPE[index] <-
      pred_error(y_l3_test, AE_result_l3$y_predicted)
    
    # Increase the iteration index
    index <- index + 1
  }
  
  # Average the ASPE over all iterations
  GCV_l2_ASPE_averaged <- round(sum(GCV_l2_ASPE) / n_study, 3)
  LCVB_l2_ASPE_averaged <- round(sum(LCVB_l2_ASPE) / n_study, 3)
  LCVT_l2_ASPE_averaged <- round(sum(LCVT_l2_ASPE) / n_study, 3)
  AE_l2_ASPE_averaged <- round(sum(AE_l2_ASPE) / n_study, 3)
  
  GCV_l3_ASPE_averaged <- round(sum(GCV_l3_ASPE) / n_study, 3)
  LCVB_l3_ASPE_averaged <- round(sum(LCVB_l3_ASPE) / n_study, 3)
  LCVT_l3_ASPE_averaged <- round(sum(LCVT_l3_ASPE) / n_study, 3)
  AE_l3_ASPE_averaged <- round(sum(AE_l3_ASPE) / n_study, 3)
  
  # Print the respective ASPEs of each method for the two regression operators
  print(
    paste(
      "Average ASPEs for the second regression operator II: ",
      GCV_l2_ASPE_averaged,
      LCVB_l2_ASPE_averaged,
      LCVT_l2_ASPE_averaged,
      AE_l2_ASPE_averaged
    )
  )
  print(
    paste(
      "Average ASPEs for the third regression operator III: ",
      GCV_l3_ASPE_averaged,
      LCVB_l3_ASPE_averaged,
      LCVT_l3_ASPE_averaged,
      AE_l3_ASPE_averaged
    )
  )
}