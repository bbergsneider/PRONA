#' Plot bar chart of severity
#'
#' This function takes in a dataframe of symptom
#' severities (formatted as required by PRONA)
#' and returns a bar chart of the severity for each
#' symptom.
#'
#' This function can also be used to plot frequency
#' data (simply change the y_label parameter)
#'
#' @param df Dataframe of symptom severities (formatted as required by PRONA)
#' @param color Color of the barplot (default green)
#' @param y_label Y-axis label
#' @return A ggplot2 object of the severity barplot
#' @export

plot_severity <- function(df, color = "#4caf4a", y_label = "severity"){
    df <- dplyr::select(df, -ID)
    tmp_df <- df
    tmp_df["comm"] <- 1
    long_df <- reshape2::melt(tmp_df, id = "comm")
    severity_plot <- ggplot2::ggplot(long_df, ggplot2::aes(x=reorder(variable, value), y=value)) + ggplot2::geom_boxplot(fill=color, na.rm=TRUE)  + ggplot2::ylab(y_label) + ggplot2::xlab('') + ggplot2::theme_minimal() + ggplot2::coord_flip() + ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 12))
    return(severity_plot)
}

#' Plot bar chart of frequency
#'
#' This is the same as plot_severity, just with different labels

# plot_frequency <- function(){

# }