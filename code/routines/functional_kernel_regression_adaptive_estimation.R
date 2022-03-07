adaptive_functional_regression <- function(y, curves, pred_curves, ...,  d = "L2", kernel = "triangle", h_length = 20){
  # This function performs functional kernel regression with adaptive bandwidth
  # estimation as presented in Chagny & Roche (2016)
  
  # y: Training set responses (as numerical vector)
  # curves: Training set curves (stored in a numerical matrix)
  # pred_curves: Test set curves (stored in a numerical matrix)
  # d: Choice of (semi-)metric (either "deriv", "pca" or "L2")
  # kernel: Choice of kernel (either "box", "quadratic" or "triangle")
  # h_length: Length of each row of the bandwidth matrix
  # ...: Additional arguments partly needed for the (semi-)metrics
  
  # Fixed constant in the estimated variance term (see Chagny & Roche Appendix A Section 1.1)
  kappa <- 0.1
  # Get the corresponding kernel function
  K <- get(kernel)
  # Set n as the number of given train responses
  n <- length(y)
  # Get the corresponding semi-metric
  d <- get(paste("d_", d, sep = ""))
  # Compute the semi-metric matrix of the training curves "curves"
  d_curves <- d(curves, curves, ...)
  # Compute the bandwidth matrix with the function H
  H_set <- H(d_curves, h_length)
  # Fix the maximal bandwidth of the smallest bandwidths
  # contained in the first column of H_set, to ensure good behaviour
  # Otherwise, data sparsity may lead to unsatisfying results for the denominator
  h_error_min <- max(H_set[,1])
  # Estimate the error variance (needed in the selection criterion)
  est_error <- sigma_hat(d_curves, y, kernel, h_error_min)
  # If the error variance is estimated to be zero, set the estimated error variance
  # to be very small (but not zero) to maintain an influence of the V_hat term
  if (est_error == 0){
    est_error <- 0.00001
  }

  # Initialise empty matrices and vectors
  V_hat <- matrix(NA, nrow = n, ncol = h_length)
  A_hat <- matrix(NA, nrow = n, ncol = h_length)
  h_AE <- 0
  a_hat <- c()
  # First, calculate the estimated variance term V_hat
  for(i in 1:n){
    for (j in 1:h_length){
      # Compute the estimated Small Ball Probabilities and V_hat for all h in the
      # bandwidth set
      h <- H_set [i,j]
      sb <- sbprob(d_curves, i, h)
      # If the estimated small ball probability is zero it is ensured that the
      # criterion is not minimal for that h
      if(sb != 0)
        {
      V_hat[i,j] <- kappa*est_error*log(n)/(n*sb)
        }
      else{
        V_hat[i,j] <- .Machine$double.xmax
        }
      }
      
      res <- c()
      # As the procedure needs recursive computations two loops are needed
      # Computing the matrix A_hat
      for (j in 1:h_length){
      for (j2 in 1:h_length){
        h <- H_set[i, j]
        h_j <- H_set[i,j2]
        func_kw_j <- K(d_curves/h_j)
        dnom_j <- apply(func_kw_j, 2, sum)
        num_j <- func_kw_j * y
        y_hat_j <- (apply(num_j, 2, sum)/dnom_j)[i]
        max_h <- max(h, h_j)
        if (max_h != h_j){
          func_kw_max_h <- K(d_curves/max_h)
          dnom_max_h <- apply(func_kw_max_h, 2, sum)
          num_max_h <- func_kw_max_h * y
          y_hat_max_h <- (apply(num_max_h, 2, sum)/dnom_max_h)[i]
          res[j2] <-  abs((y_hat_j - y_hat_max_h)^2 - V_hat[i,j2])
        }
        else{
          y_hat_max_h <- y_hat_j
          res[j2] <- abs(V_hat[i,j2])
        }
      }
        # Storing each value in vector a_hat
        a_hat <- c(a_hat, max(res))
        }
  }
  # Transforming the vector a_hat into the matrix A_hat
  A_hat <- matrix(a_hat, nrow = n, byrow=TRUE)
  eval_mat <- A_hat + V_hat
  # Now, the best bandwidths can be selected for each observed curve (local selection)
  func_kernel_weight_AE <- c()
  for(m in 1:n){
    # For each entry the optimal bandwidth corresponding to the evaluation matrix
    # eval_mat are chosen
    index_AE <- order(eval_mat[m,])[1]
    h_AE[m] <- H_set[m, index_AE]
    func_kernel_weight <- K(d_curves[m,]/h_AE[m])
    func_kernel_weight_AE <- rbind(func_kernel_weight_AE, func_kernel_weight)
  }
  dnom_AE <- apply(func_kernel_weight_AE, 2, sum)
  numerator_AE <- func_kernel_weight_AE * y
  # Estimated responses for the training set
  y_hat_AE <- apply(numerator_AE, 2, sum)/ dnom_AE
  
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
    h_AE_pred <- h_AE
    # The semi-metric matrix for the interaction between train and test set curves
    d_test <- d(curves, pred_curves, ...)
    func_kernel_weight_pred <- c()
    # Calculate the functional kernel weights with the selected bandwidths
    for(m in 1:n){
      func_kernel_weight <- K(d_test[m,]/h_AE_pred[m])
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
        index_AE <- order(eval_mat[i,])[index]
        h_AE_pred[i] <- H_set[i, index_AE]
        func_kernel_weight <- K(d_test[i,]/h_AE_pred[i])
        func_kernel_weight_pred[i,] <- func_kernel_weight
      }
      dnom_pred <- apply(func_kernel_weight_pred, 2, sum)
      dnom_pred[is.nan(dnom_pred)] <- 0
      # It is possible that the test data inherits particularly different
      # characteristics and the bandwidth set H_n may not be suitable for
      # prediction. Hence, if for the largest bandwidth in H_n the denominator
      # condition is still not satisfied, even larger bandwidths will be tried
      if (index == h_length & (sum(dnom_pred==0) != 0)){
        while(sum(dnom_pred==0)!=0){
          bool_vec <- which((dnom_pred==0)==TRUE)
          for(i in bool_vec){
            h_AE_pred[i] <- H_set[i, h_length]*1.05
            H_set[i, h_length] <- H_set[i, h_length]*1.05
            func_kernel_weight <- K(d_test[i,]/h_AE_pred[i])
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
    return(list(y_estimated = y_hat_AE, 
                y_predicted = y_hat_pred, h_local_pred = h_AE_pred, est_variance = est_error, A_hat = A_hat, V_hat = V_hat))

  }
  
  else {
    # Return several quantities of interest without predicted responses
    return(list(y_estimated = y_hat_AE, h_local = h_AE, est_variance = est_error))
  }
}

