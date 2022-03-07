# Simulations to test the effect of changing the kernel (for the GCV and AE methods)
# Using the regression operator 3 and Brownian Motions as covariables
# Results correspond to Table 4.9

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

# Fix sample size and evaluation points
n <- 150
points <- 100

# Create evaluation grid
t <- seq(0, 1, length = points)

# Set the variance of Brownian Motions
sig <- 0.01

# Initialise vectors needed for computations
GCV_l3_box_ASPE <- c()
GCV_l3_quadratic_ASPE <- c()
GCV_l3_triangle_ASPE <- c()
AE_l3_box_ASPE <- c()
AE_l3_quadratic_ASPE <- c()
AE_l3_triangle_ASPE <- c()
integral3 <- c()

# Set training and test set threshold
threshold <- 0.7 * n
# Training set index
train_index <- 1:ceiling(threshold)
# Test set index
test_index <- (ceiling(threshold) + 1):n

# Iteration index
index <- 1

# Repeat the simulations for n_study times
while (index <= n_study) {
  # Creating n Brownian Motions
  matXl <-
    matrix(rnorm(
      n = n * (length(t) - 1),
      mean = 0,
      sd = sqrt(sig)
    ), n, length(t) -
      1)
  matXl <- cbind(rep(0, n), t(apply(matXl, 1, cumsum)))
  # Compute the integral as used in regression operator III
  for (i in 1:n) {
    integral3[i] <-
      integrate(
        approxfun(t, log(1 + abs(matXl[i, ])), n = points),
        0,
        1,
        subdivisions = 250,
        stop.on.error = FALSE
      )$value
  }
  rl3 <- array(as.numeric(unlist(integral3)), dim = c(1, n, 1))
  # Create normal distributed errors
  errors <- rnorm(n, 0, 1)
  # Create responses through the non-parametric regression model
  y_l3 <- rl3 + errors
  
  # Using the index vectors to divide the data into training and test set
  train_curves <- matXl[train_index, ]
  y_l3_train <- y_l3[train_index]
  test_curves <- matXl[test_index, ]
  y_l3_test <- y_l3[test_index]
  
  # Performing functional kernel regression with global cross-validated bandwidth
  # Using the asymmetrical box kernel
  GCV_result_l3_box <-
    functional_kernel_regression_GCV(
      y_l3_train,
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
  GCV_result_l3_quadratic <-
    functional_kernel_regression_GCV(
      y_l3_train,
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
  GCV_result_l3_triangle <-
    functional_kernel_regression_GCV(
      y_l3_train,
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
  AE_result_l3_box <-
    adaptive_functional_regression(
      y_l3_train,
      train_curves,
      test_curves,
      kernel = "box",
      d = 'deriv',
      q = 2,
      nknot = 20,
      c(0, pi)
    )
  
  # Performing functional kernel regression with adaptively estimated bandwidth
  # Using the asymmetrical quadratic kernel
  AE_result_l3_quadratic <-
    adaptive_functional_regression(
      y_l3_train,
      train_curves,
      test_curves,
      kernel = "quadratic",
      d = 'deriv',
      q = 2,
      nknot = 20,
      c(0, pi)
    )
  
  # Performing functional kernel regression with adaptively estimated bandwidth
  # Using the asymmetrical triangle kernel
  AE_result_l3_triangle <-
    adaptive_functional_regression(
      y_l3_train,
      train_curves,
      test_curves,
      kernel = "triangle",
      d = 'deriv',
      q = 2,
      nknot = 20,
      c(0, pi)
    )
  
  # For each repetition of the simulation store the respective averaged squared
  # prediction errors of each method and chosen kernel in these vectors
  GCV_l3_box_ASPE[index] <-
    pred_error(y_l3_test, GCV_result_l3_box$y_predicted)
  GCV_l3_quadratic_ASPE[index] <-
    pred_error(y_l3_test, GCV_result_l3_quadratic$y_predicted)
  GCV_l3_triangle_ASPE[index] <-
    pred_error(y_l3_test, GCV_result_l3_triangle$y_predicted)
  AE_l3_box_ASPE[index] <-
    pred_error(y_l3_test, AE_result_l3_box$y_predicted)
  AE_l3_quadratic_ASPE[index] <-
    pred_error(y_l3_test, AE_result_l3_quadratic$y_predicted)
  AE_l3_triangle_ASPE[index] <-
    pred_error(y_l3_test, AE_result_l3_triangle$y_predicted)
  
  # Increase the iteration index
  index <- index + 1
}

# Average the ASPE all simulations up to the third decimal point
GCV_l3_box_ASPE_averaged <- round(sum(GCV_l3_box_ASPE) / n_study, 3)
GCV_l3_quadratic_ASPE_averaged <-
  round(sum(GCV_l3_quadratic_ASPE) / n_study, 3)
GCV_l3_triangle_ASPE_averaged <-
  round(sum(GCV_l3_triangle_ASPE) / n_study, 3)
AE_l3_box_ASPE_averaged <- round(sum(AE_l3_box_ASPE) / n_study, 3)
AE_l3_quadratic_ASPE_averaged <-
  round(sum(AE_l3_quadratic_ASPE) / n_study, 3)
AE_l3_triangle_ASPE_averaged <-
  round(sum(AE_l3_triangle_ASPE) / n_study, 3)

print(
  paste(
    "ASPE for GCV for kernels box, quadratic and triangle: ",
    GCV_l3_box_ASPE_averaged,
    GCV_l3_quadratic_ASPE_averaged,
    GCV_l3_triangle_ASPE_averaged
  )
)
print(
  paste(
    "ASPE for AE for kernels box, quadratic and triangle: ",
    AE_l3_box_ASPE_averaged,
    AE_l3_quadratic_ASPE_averaged,
    AE_l3_triangle_ASPE_averaged
  )
)
