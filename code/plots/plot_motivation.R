# Script for figure 1.1 used in the motivation
# Excerpt of the data used in Chapter 5

# Load necessary packages
library(fda.usc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Set random seed 
set.seed(12)

# Load data (from fda.usc)
data(aemet)

# Save temperature data in a data.frame
df <- as.data.frame(aemet$temp$data)
# Take a random sample of 10 curves
df_sample <- df[sample(nrow(df), 10),]

# Set colour scheme
n_colours <- dim(df_sample)[1]
getPalette <- colorRampPalette(brewer.pal(8, "Accent"))

# Plot the functional data with ggplot
df_sample %>% 
  tibble::rownames_to_column() %>% 
  gather(reading, value, -rowname) %>% 
  group_by(rowname) %>% 
  mutate(x=1:n()) %>% 
  ggplot(aes(x=x, y=value)) +
  geom_line(aes(color=rowname)) +
  scale_colour_manual(values = getPalette(n_colours))+
  xlab("Day of Year")+
  guides(color="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Mean Temperature")
