functional_kernel_regression_LCV <- function(y, curves, pred_curves, ..., d = "deriv", kernel = "triangle", h_length=20, tri_weight=TRUE)
{
  # This function performs functional kernel regression with
  # local cross-validated bandwidths automatically selected
  # as presented in Benhenni et. al. (2007)
  
  # y: Training set responses (as numerical vector)
  # curves: Training set curves (stored in a numerical matrix)
  # pred_curves: Test set curves (stored in a numerical matrix)
  # d: Choice of (semi-)metric (either "deriv", "pca" or "L2")
  # kernel: Choice of kernel (either "box", "quadratic" or "triangle")
  # h_length: Length of each row of the bandwidth matri
  # tri_weight: Boolean for the choice of weights used in the local cross-validation
  #             criterion (if TRUE, triangle weights will be used,
  #             if FALSE, box weights (constant weights) will be used.)
  # ...: Additional arguments partly needed for the (semi-)metrics
  
  # Get the corresponding kernel function from another script
  K <- get(kernel)
  # Set n as the number of given train responses
  n <- length(y)
  # Get the corresponding semi-metric
  d <- get(paste("d_", d, sep = ""))
  # Compute the semi-metric matrix of the training curves "curves"
  d_curves <- d(curves, curves, ...)
  # Compute the bandwidth matrix with the function H (see Chagny and Roche(2016))
  H_set <- H(d_curves, h_length)
  
  # Initialising the local cross-validation criterion matrix and vector for storage
  LCV <- matrix(NA, length(y), h_length)
  h_LCV <- 0
  
  # Starting with the smallest bandwidth h in H_n to calculate the kernel weights
  # for every h and every observed curve
  for(i in 1:n) {
    for(j in 1:h_length){
      h <- H_set[i,j]
      # The kernel function is applied pointwise to each entry of the d_curves matrix
      func_kernel_weight <- K(d_curves/h)
      # Here, the diagonal entries of the evaluated kernel matrix are set to zero
      # to exclude the i-th entry in the calculation of the i-th cross validated
      # kernel estimate (leave-one-curve-out-estimator)
      diag(func_kernel_weight) <- 0
      # Column-wise summation of the functional kernel weights to get the denominator
      # for each i
      dnom <- apply(func_kernel_weight, 2, sum)
      # The following calculation will only be conducted if the denominator is not zero
      if(sum(dnom == 0)==0){
        # Compute the numerator and the estimated responses
        numerator <- func_kernel_weight * y
        y_hat <- apply(numerator, 2, sum)/dnom
        # If the boolean tri_weight is TRUE use triangle weights for the criterion
        if(tri_weight==TRUE){
          local_weights <- triangle(d_curves/h)
        }
        else{
        # Else use box weights for computation of the criterion
        local_weights <- box(d_curves/h)
        }
        LCV[i,j] <- sum((y_hat - y)^2*local_weights[,i])
      }
      else{
        # If the denominator is zero set the corresponding entry of the LCV matrix
        # to infinity
        LCV[i,j] <- .Machine$double.xmax
      }
    }
  }
  
  # Store the functional kernel weights of in this vector
  func_kernel_weight_LCV <- c()
  
  # For each curve the functional kernel weights will be calculated
  for(m in 1:n){
    # Order the m-th row by ascending values, i.e. the smallest entry is given
    # on the first position
    index_LCV <- order(LCV[m,])[1]
    # Estimating the functional kernel estimator with the local cross validated
    # bandwidths h_LCV at each curve
    h_LCV[m] <- H_set[m, index_LCV]
    func_kernel_weight <- K(d_curves[m,]/h_LCV[m])
    func_kernel_weight_LCV <- rbind(func_kernel_weight_LCV, func_kernel_weight)
  }
  
  # Further calculations for the Functional Kernel Estimator
  dnom_LCV <- apply(func_kernel_weight_LCV, 2, sum)
  numerator_LCV <- func_kernel_weight_LCV * y
  # Estimated responses of the training set
  y_hat_LCV <- apply(numerator_LCV, 2, sum)/ dnom_LCV
  
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
    h_LCV_pred <- h_LCV
    # The semi-metric matrix for the interaction between train and test set curves
    d_test <- d(curves, pred_curves, ...)
    func_kernel_weight_pred <- c()
    # Calculate the functional kernel weights with the selected bandwidths
    for(m in 1:n){
      func_kernel_weight <- K(d_test[m,]/h_LCV_pred[m])
      func_kernel_weight_pred <- rbind(func_kernel_weight_pred, func_kernel_weight)
    }
    dnom_pred <- apply(func_kernel_weight_pred, 2, sum)
    # Set index
    index <- 2
    # As long as the denominator is zero larger bandwidths have to be chosen
    while(sum(dnom_pred==0) != 0){
      # Adjust bandwidths only where the denominator is zero
      bool_vec <- which((dnom_pred==0)==TRUE)
      for (i in bool_vec){
        index_LCV <- order(LCV[i,])[index]
        h_LCV_pred[i] <- H_set[i, index_LCV]
        func_kernel_weight <- K(d_test[i,]/h_LCV_pred[i])
        func_kernel_weight_pred[i, ] <- func_kernel_weight
      }
      dnom_pred <- apply(func_kernel_weight_pred, 2, sum)
      # It is possible that the test data inherits particularly different
      # characteristics and the bandwidth set H_n may not be suitable for
      # prediction. Hence, if for the largest bandwidth in H_n the denominator
      # condition is still not satisfied, even larger bandwidths will be tried
      if (index == h_length & (sum(dnom_pred==0) != 0)){
        while(sum(dnom_pred==0)!=0){
          bool_vec <- which((dnom_pred==0)==TRUE)
          for(i in bool_vec){
            h_LCV_pred[i] <- H_set[i, h_length]*1.05
            H_set[i, h_length] <- H_set[i, h_length]*1.05
            func_kernel_weight <- K(d_test[i,]/h_LCV_pred[i])
            func_kernel_weight_pred[i,] <- func_kernel_weight
          }
          dnom_pred <- apply(func_kernel_weight_pred, 2, sum)
        }
        # If the condition is satisfied, stop the while-loop
        break
      }
      index <- index + 1
    }
      # Calculate the predicted responses with the final functional kernel weights
      numerator_pred <- func_kernel_weight_pred * y
      y_hat_pred <- apply(numerator_pred, 2, sum)/dnom_pred
      y_hat_pred[is.nan(y_hat_pred)] <- 0
      # Return several quantities of interest (especially predicted responses)
      return(list(y_estimated = y_hat_LCV, 
                  y_predicted = y_hat_pred, h_local_pred = h_LCV_pred))
  }
  else {
    # Return the estimated training responses and the locally chosen bandwidths
    return(list(y_estimated = y_hat_LCV, h_local = h_LCV))
  }
}
