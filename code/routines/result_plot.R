res_plot <- function(response, predicted, ASPE, label){
  # Function to plot true responses against predicted/estimated responses
  # while showing the ASPE in the plot title
  # response: the true response vector
  # predicted: the predicted/estimated response vector
  # ASPE: Computed Averaged Squared Prediction Error
  # label: String for naming the method corresponding to the ASPE
  plot_df <- as.data.frame(cbind(response, predicted))
  colnames(plot_df) <- c("True", "Predicted")
  plot <- ggplot(data=plot_df, aes(x=True, y=Predicted)) +
    geom_point(show.legend=FALSE) +
    geom_abline() +
    ggtitle(paste0("ASPE (",label, ") =", ASPE)) +
    xlab("True Response") +
    ylab("Predicted Response") +
    theme(plot.title = element_text(hjust=0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  return(plot)
}

