# Creating the Nadaraya-Watson estimator plot (univariate case)
# Plot the estimator for two bandwidths, see Figure 2.2

# Set random seed for reproducibility purposes
set.seed(121)

# Set sample size
n <- 150
# Evaluation grid
x <- seq(0, 1, length.out = n)
# Create normal distributed error vector
eps <- rnorm(n,0,0.1)
# Compute the responses through the univariate non-parametric regression model
y <- sin(2*x)- exp(x^3) + eps

# Plot the estimators and observed data
plot(x,y,pch = 20,xlab ="X", ylab = "Y",cex = 0.6,cex.axis = 1, cex.lab =1.2, cex.main = 1.2)
fit004 <- ksmooth(x,y, bandwidth=0.04)
lines(fit004$x,fit004$y, col = "cyan4")
fit04 <- ksmooth(x,y, bandwidth = 0.4)
lines(fit04$x, fit04$y, col = "indianred2")
legend("bottomleft", cex = 1,  legend = c("0.04", "0.4"), col = c("cyan4", "indianred2"), lty =1, title = "Bandwidths")
