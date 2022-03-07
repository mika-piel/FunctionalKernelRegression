functional_kernel_regression_GCV <- function(y, curves, pred_curves, ..., d = "deriv", kernel = "triangle", h_length = 50)
{
  # This function performs functional kernel regression with global
  # cross-validated bandwidths automatically selected as presented in
  # Rachdi & Vieu (2006)
  
  # y: Training set responses (as numerical vector)
  # curves: Training set curves (stored in a numerical matrix)
  # pred_curves: Test set curves (stored in a numerical matrix)
  # d: Choice of (semi-)metric (either "deriv", "pca" or "L2")
  # kernel: Choice of kernel (either "box", "quadratic" or "triangle")
  # h_length: Number of bandwidths used for the global cv procedure
  # ...: Additional arguments partly needed for the (semi-)metrics
  
  # Get the corresponding kernel function from the auxiliary_functions script
  K <- get(kernel)
  
  d <- get(paste("d_", d, sep = ""))
  d_curves <- d(curves, curves, ...)
  # Find row containing the maximal entry
  max_index <- which(d_curves == max(d_curves), arr.ind=TRUE)
  dh <- as.matrix(d_curves[max_index[1],])
  H_n <- H(dh, h_length)
  
  # Initialising count indices i and j
  i <- 1
  j <- 1
  # Initialising the global cross-validation criterion
  GCV <- 0
  # Save bandwidths which are not zero and satisfy that no entry of the summed
  # denominator contains value zero in a separate vector
  H_n_nonzero <- 0
  # Starting with the smallest bandwidth h in H_n to calculate the kernel weights
  h <- H_n[1]
  # Repeat the following procedure for every bandwidth
  while(i < length(H_n)) {
    # First, incrementally increase the count index
    i <- i + 1
    # The kernel function is applied pointwise for each entry of the d_curves matrix
    func_kernel_weight <- K(d_curves/h)
    # Here, the diagonal entries of the evaluated kernel matrix are set to zero
    # to exclude the i-th entry in the calculation of the i-th cross validated
    # kernel estimate (leave-one-curve-out-estimator)
    diag(func_kernel_weight) <- 0
    # Row-wise summation of the functional kernel weights to get the denominator
    # for each i
    dnom <- apply(func_kernel_weight, 2, sum)
    # The following calculation will only be conducted if the denominator is not zero
    if(sum(dnom == 0)==0){
      numerator <- func_kernel_weight * y
      y_hat <- apply(numerator, 2, sum)/dnom
      # Bandwidths which are not zero and for which no entry of the dnom matrix
      # is zero are saved in a separate array
      H_n_nonzero[j] <- h
      # The global cross validation criterion
      GCV[j] <- sum((y_hat - y)^2)
      j <- j + 1
    }
    # Previous procedure will be conducted again for the next bandwidths
    h <- H_n[i]
  }
  # Choosing the bandwidth which minimises the global cross-validation criterion
  # Order(GCV)[1] returns the index for which the GCV criterion is the smallest
  index_GCV <- order(GCV)[1]
  # Estimating the functional kernel estimator with the global cross validated
  # bandwidth h_GCV
  h_GCV <- H_n_nonzero[index_GCV]
  func_kernel_weight_GCV <- K(d_curves/h_GCV)
  dnom_GCV <- apply(func_kernel_weight_GCV, 2, sum)
  numerator_GCV <- func_kernel_weight_GCV * y
  # Computing the estimated responses of the training set
  y_hat_GCV <- apply(numerator_GCV, 2, sum)/ dnom_GCV

  # If pred_curves are different than curves (i.e. when pred_curves is the test set
  # and curves is the training set)
  prediction <- FALSE
  if (!identical(curves, pred_curves)){
    if(ncol(curves) == ncol(pred_curves)){
      prediction <- TRUE
    }
  }
  
  # Perform prediction on the test set
  if(prediction){
    # The semi-metric matrix for the interaction between train and test set curves
    d_test <- d(curves, pred_curves, ...)
    # Calculate the functional kernel weights with the global cross-validated
    # bandwidth h_GCV
    kernel_weight <- K(d_test/h_GCV)
    dnom_pred <- apply(kernel_weight, 2, sum)
    # If the cross-validated bandwidth also can be applied to the test set
    # continue with computing the quantities of the functional kernel regression
    if(sum(dnom_pred==0) == 0){
      numerator_pred <- kernel_weight * y
      y_hat_pred <- apply(numerator_pred, 2, sum)/dnom_pred
      y_hat_pred[is.nan(y_hat_pred)] <- 0
      return(list(y_estimated = y_hat_GCV, 
                  y_predicted = y_hat_pred, h_GCV = h_GCV, 
                  H_n_GCV = H_n_nonzero))
    }
    # If the denominator contains values zero (on the test set) 
    # the bandwidths from the set will used in the order of an ascending GCV
    # criterion
    else {
      index <- 2
      # As long as the denominator is zero for some value repeat the following
      # procedure
      while(sum(dnom_pred==0)!=0){
        index_GCV <- order(GCV)[index]
        h_GCV <- H_n_nonzero[index_GCV]
        kernel_weight <- K(d_test/h_GCV)
        dnom_pred <- apply(kernel_weight, 2, sum)
        # It is possible that the test data inherits particularly different
        # characteristics and the bandwidth set H_n may not be suitable for
        # prediction. Hence, if for the largest bandwidth in H_n the denominator
        # condition is still not satisfied, even larger bandwidths will be tried
        if (index == length(H_n_nonzero) & (sum(dnom_pred==0) != 0)){
          while(sum(dnom_pred==0)!=0){
              h_GCV <- H_set[i, h_length]*1.05
              H_set[i, h_length] <- H_set[i, h_length]*1.05
              func_kernel_weight <- K(d_test[i,]/h_LCV_pred[i])
              func_kernel_weight_pred[i,] <- func_kernel_weight
          
            dnom_pred <- apply(func_kernel_weight_pred, 2, sum)
          }
          # If the condition is satisfied, stop the while-loop
          break
        }
        index <- index + 1
      }
      # Calculate the estimator with the final kernel weights
      numerator_pred <- kernel_weight * y
      y_hat_pred <- apply(numerator_pred, 2, sum)/dnom_pred
      # Return list of interesting quantities (especially predicted responses)
      return(list(y_estimated = y_hat_GCV, 
                  y_predicted = y_hat_pred, h_GCV = h_GCV))
    }
  }
  else {
    # Return list with estimated training responses and chosen bandwidth
    return(list(y_estimated = y_hat_GCV, h_GCV = h_GCV))
  }
}

