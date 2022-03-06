fct_plot <- function(df, ylabel = "Curves"){
  # Function to plot functional data
  # df: Dataframe which contains the functional data (i.e. the curves)
  # ylabel: Label of the y-axis (if not specified, generically chosen as "curves")
  
  # Set a fixed colour scheme
  n_colours <- dim(df)[1]
  getPalette <- colorRampPalette(brewer.pal(8, "Accent"))
  # Plot the curves via ggplot and specified theme
  # Naturally, as the data is discrete a curve approximation for each observation
  # will be plotted via geom_line
  plot_curves <- df %>% 
    tibble::rownames_to_column() %>% 
    gather(reading, value, -rowname) %>% 
    group_by(rowname) %>% 
    mutate(x=1:n()) %>% 
    ggplot(aes(x=x, y=value)) +
    geom_line(aes(color=rowname)) +
    scale_colour_manual(values = getPalette(n_colours)) +
    guides(color="none") + 
    ylab(ylabel)+
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.x=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}
