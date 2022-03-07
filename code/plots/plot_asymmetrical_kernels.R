# Here the asymmetrical box and quadratic kernel are plotted, Figure 2.1

# Load necessary packages
require(gridExtra)
require(ggplot2)

# Load functions and routines needed for this script
# Please save the scripts below in the same folder as this script
source("auxiliary_functions.R")

# Number of evaluation points
points <- 500
# Create evaluation grid
x <- seq(-1, 2, length.out = points)
# Apply the kernel functions to x
y_b <- box(x)
y_q <- quadratic(x)
y_t <- triangle(x)
# Save data in a data.frame
d_b <- data.frame(x, y_b)
d_q <- data.frame(x, y_q)
d_t <- data.frame(x,y_t)

# Plot the kernel functions
plot1 <- ggplot(d_b, aes(x=x, y=y_b)) + 
  geom_line() +
  ylim(0,2) +
  theme(axis.title.x = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("(a)") + ylab("")
                                      
plot2 <- ggplot(d_q, aes(x=x, y=y_q)) +
  geom_line() + 
  ylim(0,2) + 
  theme(axis.title.x = element_text(size = 14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("(b)") + ylab("")

plot3 <- ggplot(d_t, aes(x=x, y=y_t)) +
  geom_line() + 
  ylim(0,2) + 
  theme(axis.title.x = element_text(size = 14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("(c)") + ylab("")

# Arrange the three plots in one plot
grid.arrange(plot1, plot2, plot3, ncol=3)

