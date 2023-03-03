#' Plot bar chart of severity
#'
#' This function takes in a dataframe of symptom
#' severities (formatted as required by PRONA)
#' and returns a bar chart of the severity for each
#' symptom.
#'
#'
#' @param df Dataframe of symptom severities (formatted as required by PRONA)
#' @param color Color of the barplot (default green)
#' @param y_label Y-axis label
#' @return A ggplot2 object of the severity barplot
#' @export

plot_severity <- function(df, color = "#4caf4a", y_label = "severity"){
    check_data_format(df) # Check if df is formatted properly
    df <- dplyr::select(df, -ID)
    tmp_df <- df
    tmp_df["comm"] <- 1
    long_df <- reshape2::melt(tmp_df, id = "comm")
    severity_plot <- ggplot2::ggplot(long_df, ggplot2::aes(x=reorder(variable, value), y=value)) +
        ggplot2::geom_boxplot(fill = color, na.rm = TRUE) +
        ggplot2::ylab(y_label) +
        ggplot2::xlab("") +
        ggplot2::theme_minimal() +
        ggplot2::coord_flip() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 12))
    return(severity_plot)
}

#' Plot bar chart of frequency
#'
#' Same function as plot_severity, just with different default y-axis label
#' It is provided to make it clear to users that PRONA can be used to analyze
#' both severity and frequency data.
#'
#' @param df Dataframe of symptom frequencies (formatted as required by PRONA)
#' @param color Color of the barplot (default green)
#' @param y_label Y-axis label
#' @return A ggplot2 object of the frequency barplot
#' @export

plot_frequency <- function(df, color = "#4caf4a", y_label = "frequency"){
    plot_severity(df, color, y_label)
}

#' Plot bar chart of occurrence
#'
#' This function takes in a dataframe of symptom severities or frequencies
#' (formatted as required by PRONA) and returns a bar chart of the percent
#' of times each symptom occurs across the entire cohort.
#'
#' @param df Dataframe of symptom severities or frequencies (PRONA format)
#' @param color Color of the barplot (default blue)
#' @param y_label Y-axis label
#' @param cutoff Severity/Frequency cutoff to define occurrence at (default: 1)
#' @return A ggplot2 object of the occurrence barplot
#' @export

plot_occurrence <- function(df, color = "royalblue", y_label = "occurrence", cutoff = 1) {
    check_data_format(df) # Check if df is formatted properly
    df <- dplyr::select(df, -ID)
    symptom_frequencies <- round(colSums(df >= cutoff, na.rm = TRUE) / nrow(df) * 100, 1)
    frequency_df <- data.frame(symptom = names(symptom_frequencies), frequency = unname(symptom_frequencies))

    frequency_plot <- ggplot2::ggplot(data=frequency_df, ggplot2::aes(x = reorder(symptom, frequency), y = frequency)) +
        ggplot2::geom_bar(stat = "identity", color = 'black', fill = color, alpha=0.75, width=0.8) +
        ggplot2::ylab("% occurrence") +
        ggplot2::xlab("") +
        ggplot2::ylim(0, 100) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", frequency)),hjust=-0.3, size=4) +
        ggplot2::theme_minimal() +
        ggplot2::coord_flip() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 12))

    return(frequency_plot)
}