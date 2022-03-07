# This script applies functional kernel regression and corresponding bandwidth
# selection methods to a real data set given by the fda.usc package (aemet)
# Object of the regression analysis are averaged yearly wind speeds while
# the covariables are yearly temperature curves (with daily measurements)

# Load necessary packages
library(fda.usc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(gridExtra)

# Load functions and routines needed to perform the regression analysis of the data
# Please save the scripts below in the same folder as this script
source("auxiliary_functions.R")
source("functional_kernel_regression_GCV.R")
source("functional_kernel_regression_LCV.R")
source("functional_kernel_regression_adaptive_estimation.R")
source("fct_plot.R")
source("result_plot.R")

# Set random seed for the random separation of training and test data
set.seed(112213)

# Load the data (from fda.usc)
data(aemet)
# Functional data (covariables)
temp_curves <- as.data.frame(aemet$temp$data)
# Average the wind speed to gain a scalar response
y <- (1/365) * apply(aemet$wind.speed$data, 1, sum)

# Set atrributes and names of the data to NULL
attributes(temp_curves)$dimnames <- NULL
names(y) <- NULL

# Get sample index vector
smple <- sample.int(73, 73)

# Get train and test indices
train <- smple[1:51]
test <- smple[52:73]

# Separate the data into training and test data
curves_train <- as.matrix(temp_curves[train,])
curves_test <- as.matrix(temp_curves[test,])
y_train <- y[train]
y_test <- y[test]

# Save both training and test data into data.frames
df_train <- as.data.frame(curves_train)
df_test <- as.data.frame(curves_test)

# Set colour scheme
n_colours_train <- dim(df_train)[1]
n_colours_test <- dim(df_test)[1]
getPalette <- colorRampPalette(brewer.pal(8, "Accent"))

# Plot the functional train data (Figure 5.1)
plot_train <- df_train %>% 
  tibble::rownames_to_column() %>% 
  gather(reading, value, -rowname) %>% 
  group_by(rowname) %>% 
  mutate(x=1:n()) %>% 
  ggplot(aes(x=x, y=value)) +
  geom_line(aes(color=rowname)) +
  scale_colour_manual(values = getPalette(n_colours_train))+
  xlab("Day of Year")+
  guides(color="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Mean Temperature (Training set)")

# Plot the functional test data (Figure 5.2)
plot_test <- df_test %>% 
  tibble::rownames_to_column() %>% 
  gather(reading, value, -rowname) %>% 
  group_by(rowname) %>% 
  mutate(x=1:n()) %>% 
  ggplot(aes(x=x, y=value)) +
  geom_line(aes(color=rowname)) +
  scale_colour_manual(values = getPalette(n_colours_test))+
  xlab("Day of Year")+
  guides(color="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Mean Temperature (Test set)")
grid.arrange(plot_train, plot_test, nrow=2)

# Again, set rownames to NULL
rownames(curves_train) <- NULL
rownames(curves_test) <- NULL

# Now, plot the estimated small ball probabilities of the training data against
# bandwidths
h <- seq(0.0001, 1000, length=500)
dT_deriv <- d_deriv(curves_train, curves_train, q=2, nknot=20, c(0,1))

# Initialise empty vector
s_deriv <- c()

# Compute the estimated small ball probabilities (for semi-metric d.deriv)
for(i in 1:length(h)){
  s_deriv[i] <- sbprob(dT_deriv, 1, h[i])
}

# Save data in a data.frame
dfT_plot1 <- data.frame(h, s_deriv)

# Plot the estimated probabilities against the bandwidths (for semi-metric d.deriv)
plot1 <- ggplot(data=dfT_plot1, aes(x=h, y=s_deriv)) +
  geom_line(show.legend=FALSE) +
  xlab("Bandwidths") +
  ylab("Estimated Small Ball Probability") +
  ggtitle("Semi-metric based on derivatives") +
  theme(plot.title = element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

dT_pca <- d_pca(curves_train, curves_train, q=2)

s_pca <- c()

# Compute the estimated small ball probabilities (for semi-metric d.pca)
for(i in 1:length(h)){
  s_pca[i] <- sbprob(dT_pca, 1, h[i])
}

# Save data in a data.frame
dfT_plot2 <- data.frame(h, s_pca)

# Plot the estimated probabilities against the bandwidths (for semi-metric d.pca)
plot2 <- ggplot(data=dfT_plot2, aes(x=h, y=s_pca)) +
  geom_line(show.legend=FALSE) +
  xlab("Bandwidths") +
  ylab("Estimated Small Ball Probability") +
  ggtitle("Semi-metric based on principal components") +
  theme(plot.title = element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Show both plots next to each other (Figure 5.3)
grid.arrange(plot1, plot2, nrow=2, ncol=1)

# Now, perform functional kernel regression and predict the true responses of
# the test set

# Functional kernel regression with global cross-validated bandwidth
GCV_result_tc <-
  functional_kernel_regression_GCV(
    y_train,
    curves_train,
    curves_test,
    kernel = 'quadratic',
    d = 'pca',
    q = 2
  )

# Functional kernel regression with local cross-validated bandwidth
# using box weights
LCVB_result_tc <-
  functional_kernel_regression_LCV(
    y_train,
    curves_train,
    curves_test,
    kernel = 'quadratic',
    d = 'pca',
    q = 2,
    tri_weight = FALSE
  )

# Functional kernel regression with local cross-validated bandwidth
# using triangle weights
LCVT_result_tc <-
  functional_kernel_regression_LCV(
    y_train,
    curves_train,
    curves_test,
    kernel = 'quadratic',
    d = 'pca',
    q = 2,
    tri_weight = TRUE
  )

# Functional kernel regression with adaptively estimated bandwidth
AE_result_tc <-
  adaptive_functional_regression(
    y_train,
    curves_train,
    curves_test,
    kernel = 'quadratic',
    d = 'pca',
    q = 2
  )

# Calculating the Averaged Squared Prediction Errors
GCV_tc_mse <- pred_error(y_test, GCV_result_tc$y_predicted)
LCVB_tc_mse <- pred_error(y_test, LCVB_result_tc$y_predicted)
LCVT_tc_mse <- pred_error(y_test, LCVT_result_tc$y_predicted)
AE_tc_mse <- pred_error(y_test, AE_result_tc$y_predicted)

# Print the errors fpr each method
print(paste("Average ASPEs for predicting the wind speed given the respective temperature curves: ", GCV_tc_mse, LCVB_tc_mse, LCVT_tc_mse, AE_tc_mse))

# Plot true responses against predicted responses for each bandwidth method
plot_predicted_GCV_tc <- res_plot(y_test, GCV_result_tc$y_predicted, GCV_tc_mse, "GCV")
plot_predicted_LCVB_tc <- res_plot(y_test, LCVB_result_tc$y_predicted, LCVB_tc_mse, "LCVB")
plot_predicted_LCVT_tc <- res_plot(y_test, LCVT_result_tc$y_predicted, LCVT_tc_mse, "LCVT")
plot_predicted_AE_tc <- res_plot(y_test, AE_result_tc$y_predicted, AE_tc_mse, "AE")
grid.arrange(plot_predicted_GCV_tc, plot_predicted_LCVB_tc,  plot_predicted_LCVT_tc, plot_predicted_AE_tc, nrow=2, ncol=2)
