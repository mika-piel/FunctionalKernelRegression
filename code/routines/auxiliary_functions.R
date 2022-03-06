# This script contains supplementary functions needed for Functional Kernel
# Regression and corresponding bandwidth selection methods

ind <- function(x, a, b, closed=FALSE){
  # Indicator function for x in the interval [a,b] or [a,b) depending on the
  # boolean "closed"
  b0 <- x<a
  if(closed==FALSE){
  b1 <- x >=b
  }
  else{
    b1 <- x>b
  }
  bool <- as.logical((1-b0)*(1-b1))
  x[b0] <- 0
  x[b1] <- 0
  x[bool] <- 1
  return(x)
}


box <- function(x)
{
  # Asymmetrical Box kernel function
  out <- ind(x,0,1)
  return(out)
}

quadratic <- function(x)
{
  # Asymmetrical Quadratic kernel function
  out <- (6/4)*(1 - x^2)*ind(x, 0,1)
  return(out)
}

triangle <- function(x)
{
  # Asymmetrical Triangle kernel function
  out <- 2*(1-x)*ind(x,0,1)
  return(out)
}


sbprob <- function(d_curves, i, h)
{
  # Function to calculate the empirical/estimated Small Ball Probability
  # d_curves: Matrix of semimetric evaluations (d(X_i, X_j)) for some dataset containing
  #           functional data (i.e. discretisation of curves)
  # i: Index for which column the probability vector wants to be computed
  # h: Fixed bandwidth for which the Small Ball Probability shall be estimated
  
  # Set n as the number of columns of d_curves
  n <- dim(d_curves)[2]
  # Create a vector containing the empirical probability of for the given bandwidth h
  pvec <- ind(d_curves/h, 0, 1, closed=TRUE)[,i]
  # Sum the probability vector
  prob <- 1/n * sum(pvec)
  # prob[is.nan(prob)] <- 0
  # Return the estimated Small Ball Probability
  return(prob)
}

H <- function(d_curves, h_length){
  # This function returns local bandwidth sets (columns) for each observed 
  # function X_i (rows) stored in a bandwidth matrix H_curves
  # Reasoning from Chagny and Roche (2016)
  # d_curves: Matrix of semimetric evaluations for some dataset containing
  #           functional data (i.e. discretisation of curves)
  # h_length: Determines the cardinality of the bandwidth sets stored in each row
  
  # Initialising the empty Bandwidth matrix
  H <- matrix(NA, dim(d_curves)[2], h_length)
  # Find the maximum of each column of d_curves for the upper bandwidth bound
  h_max <- apply(d_curves, 2, max)
  # Set count index for the while loop
  k_max <- 1
  # Set n as the number of columns of d_curves
  n <- dim(d_curves)[2]
  # Iterate over each column
  for(index in 1:n){
    # Increase k_max as long as the empirical quantile of the Small Ball Probabilities
    # do not fall below log(n)/n (for each column)
    while (sbprob(d_curves, index,  h_max[index]/k_max) >= log(n)/n)
    {
      if (k_max == 100){
        break
      }
      else{
        k_max <- k_max + 1
      }
    }
    # With the previously computed k_max the minimal bandwidth to satisfy the while
    # loop condition can be set
    h_min <- h_max[index] / k_max
    # Then, create a bandwidth vector (set) with length h_length ranging from h_min
    # to h_max
    H <- rbind(H, seq(h_min, h_max[index], length=h_length))
  }
  # Delete the initially created empty bandwidth matrix which was used for correct
  # storing purposes
  H_curves <- H[-(1:n),]
  # Return the bandwidth matrix
  return(H_curves)
}


sigma_hat <- function(d_curves, y, kernel, H_n)
{
  # Function to estimate the variance of the errors given in the regression model
  # Reasoning from Chagny and Roche (2016)
  # d_curves: Matrix of semimetric evaluations (d(X_i, X_j)) for some dataset containing
  #           functional data (i.e. discretisation of curves)
  # y: True responses of the regression model given as a vector
  # kernel: Name of the kernel which is used for computation
  # H_n: Bandwidth set used for the Kernel Regression
  
  # Save y as a vector (if not already the case)
  y <- as.vector(y)
  # Get the corresponding kernel function
  K <- get(kernel)
  # Set h_min as the minimum of the bandwidth set
  h_min <- min(H_n)
  # Calculate the Functional Kernel weights for the minimal bandwidth
  func_kernel_weight <- K(d_curves/h_min)
  # In the adaptive bandwidth estimation procedure the function H ensures that
  # the minimal bandwidth of H_n ensures a non-zero Kernel Regression denominator
  dnom <- apply(func_kernel_weight, 2, sum)
  # Compute the numerator of the Functional Kernel Estimator
  numerator <- func_kernel_weight * y
  # Combining the previous results yields the predicted responses
  y_hat_min <- apply(numerator, 2, sum)/dnom
  # The squared error is the deviation of y_hat_min from the actual responses y
  sq_error <- (y - y_hat_min)^2
  # Sum each entry of sq_error and divide by the number of responses to gain
  # the estimated noise or error variance
  est_noise <- sum(sq_error)/length(y)
  # Return the estimated error variance
  return(est_noise)
}


# The following three function concerning the (semi-)metrics were taken from
# the accompanying material of Ferraty and Vieu (2006) and slightly adapted.
# The code can be found here: http://www.math.univ-toulouse.fr/staph/npfda/
  
d_L2 = function(curves1, curves2)
{
  # # From http://www.math.univ-toulouse.fr/staph/npfda/ by Ferraty and Vieu (2006)
  # Computes the approximated L2-distance between two curves stored in two 
  # datasets
  #    "curves1" matrix contains a first set of curves stored row by row
  #    "curves2" matrix contains a second set of curves stored row by row
  # Returns a  matrix containing the L2-distances between all entry pairs
  
  d = matrix(0, nrow(curves1), nrow(curves2))
  for(i in 1:nrow(curves2)){
    delta = t(curves1) - curves2[i, ]
    d[, i] = apply(delta^2, 2, sum) / (ncol(curves2) - 1)
  }
  return(sqrt(d))
}

d_deriv <- function(curves1, curves2, q, nknot, grid_range)
{
  # From http://www.math.univ-toulouse.fr/staph/npfda/ by Ferraty and Vieu (2006)
  # Computes a semimetric between curves based on their derivatives.
  #    "curves1" matrix containing a first set of curves stored row by row
  #    "curves2" matrix containing a second set of curves stored row by row
  #    "q" order of derivation
  #    "nknot" number of interior knots (needed for defining the B-spline basis)
  #    "grid_range" vector of length 2 containing the range of the grid at 
  #                 which the curve are evaluated (i.e. range of the 
  #                 discretization)
  # Returns a "semimetric" matrix containing the semimetric computed 
  # between the curves lying to the first sample and the curves lying  
  # to the second one.

  if(is.vector(curves1)) curves1 <- as.matrix(t(curves1))
  if(is.vector(curves2)) curves2 <- as.matrix(t(curves2))
  testfordim <- sum(dim(curves1)==dim(curves2))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(curves1==curves2)!=prod(dim(curves1))

  # B-spline approximation of the curves containing in DATASET :
  # -----------------------------------------------------------
  # "knot" and "x" allow to define the B-spline basis
  # "coef.mat1[, i]" corresponds to the B-spline expansion
  # of the discretized curve contained in DATASET[i, ]. 
  # The B-spline approximation of the curve contained in "curves1[i, ]" 
  # is given by "Bspline %*% coef.mat1[, i]"
  
  p <- ncol(curves1)
  a <- grid_range[1]
  b <- grid_range[2]
  x <- seq(a, b, length = p)
  order.Bspline <- q + 3
  nknotmax <- (p - order.Bspline - 1)%/%2
  if(nknot > nknotmax){
    stop(paste("give a number nknot smaller than ", nknotmax, " for avoiding ill-conditioned matrix"))
  }
  Knot <- seq(a, b, length = nknot + 2)[ - c(1, nknot + 2)]
  delta <- sort(c(rep(c(a, b), order.Bspline), Knot))
  Bspline <- splineDesign(delta, x, order.Bspline)
  Cmat <- crossprod(Bspline)
  Dmat1 <- crossprod(Bspline, t(curves1))
  coef.mat1 <- symsolve(Cmat, Dmat1)
  
  # Numerical integration by the Gauss method :
  # -------------------------------------------
  # The objects ending by "gauss" allow us to compute numerically  
  # integrals by means the "Gauss method" (lx.gauss=6 ==> the computation 
  # of the integral is exact for polynom of degree less or equal to 11).
  
  point.gauss <- c(-0.9324695142, -0.6612093865, -0.2386191861, 
                   0.2386191861, 0.6612093865, 0.9324695142)
  weight.gauss <- c(0.1713244924, 0.360761573, 0.4679139346, 0.4679139346,0.360761573, 0.1713244924)
  x.gauss <- 0.5 * ((b + a) + (b - a) * point.gauss)
  lx.gauss <- length(x.gauss)
  Bspline.deriv <- splineDesign(delta, x.gauss, order.Bspline, rep(q, lx.gauss))
  H <- t(Bspline.deriv) %*% (Bspline.deriv * (weight.gauss * 0.5 * (b - a)))
  eigH <- eigen(H, sym = T)
  eigH$values[eigH$values < 0] <- 0
  Hhalf <- t(eigH$vectors %*% (t(eigH$vectors) * sqrt(eigH$values)))
  coef1 <- t(Hhalf %*% coef.mat1)
  if(twodatasets){
    Dmat2 <- crossprod(Bspline, t(curves2))
    coef.mat2 <- symsolve(Cmat, Dmat2)
    coef2 <- t(Hhalf %*% coef.mat2)
  } else {
    coef2 <- coef1
  }
  d <- 0
  nbasis <- nrow(H)
  for(f in 1:nbasis)
    d <- d + outer(coef1[, f], coef2[, f], "-")^2
  return(sqrt(d))
}

d_pca <- function(curves1, curves2, q)
{
  # From http://www.math.univ-toulouse.fr/staph/npfda/ by Ferraty and Vieu (2006)
  # Computes between curves a pca-type semimetric based on the
  # functional principal components analysis method.
  #    "curves1" matrix containing a first set of curves stored row by row
  #    "curves2" matrix containing a second set of curves stored row by row
  #    "q" the retained number of principal components
  # Returns a "semimetric" matrix containing the semimetric computed 
  # between the curves lying to the first sample and the curves lying  
  # to the second one.
  
  if(is.vector(curves1)) curves1 <- as.matrix(t(curves1))
  if(is.vector(curves2)) curves2 <- as.matrix(t(curves2))
  testfordim <- sum(dim(curves1)==dim(curves2))==2
  twodatasets <- T
  if(testfordim) twodatasets <- sum(curves1==curves2)!=prod(dim(curves1))
  qmax <- ncol(curves1)
  if(q > qmax) stop(paste("give a integer q smaller than ", qmax))
  n <- nrow(curves1)
  cov <- t(curves1) %*% curves1/n
  eigenvecs <- eigen(cov, sym = T)$vectors[, 1:q]
  comp1 <- curves1 %*% eigenvecs
  if(twodatasets) {
    comp2 <- curves2 %*% eigenvecs
  }
  else {
    comp2 <- comp1
  }
  d <- 0
  for(k in 1:q)
    d <- d + outer(comp1[, k], comp2[, k], "-")^2
  return(sqrt(d))
}

