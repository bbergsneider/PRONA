#' Check normality of every variable in a symptom severity/frequency dataframe
#'
#' This function runs every symptom through the Shapiro-Wilk normality tests
#' and prints the results. A p-value < 0.05 in the Shapiro-Wilk normality test
#' suggests that the variable is not normally distributed. This function
#' is helpful in knowing whether or not to use a non-paranormal transformation
#' when constructing a GGM.
#'
#' @param df Dataframe representing symptom severity/frequency (PRONA format)
#' @export

check_normality <- function(df) {
    check_data_format(df) # Check if df is formatted properly
    df <- dplyr::select(df, -ID) # Get rid of ID column
    apply(df, 2, stats::shapiro.test) # Perform Shapiro-Wilk normality test
}



#' Plot density distribution of every variable in a symptom severity/frequency dataframe
#'
#' This function plots a ridge plot for every variable the dataframe.
#'
#' @importFrom magrittr %>%
#' @param df Dataframe representing symptom severity/frequency (PRONA format)
#' @param color Color of the ridge plot (Default: Green)
#' @param x_label X-axis label (Default: "severity distribution")
#' @return A ggplot2 object representing the ridge plot of symptom severity/frequency
#' @export

plot_density <- function(df, color = "#4caf4a", x_label = "severity\ndistribution") {
    check_data_format(df) # Check if df is formatted properly
    df <- dplyr::select(df, -ID) # Get rid of ID column

    # Format data as needed to make a ridge plot
    df_long <- df %>%
        data.frame() %>%
        tidyr::pivot_longer(tidyr::everything(), names_to = "symptom", values_to="severity") %>%
        dplyr::rowwise() %>%
        dplyr::mutate(symptom = gsub("\\.", "-", symptom))

    # Plot ridge plot
    ridge_plot <- ggplot2::ggplot(df_long, ggplot2::aes(x=severity, y=symptom)) +
        ggridges::geom_density_ridges(fill = color) +
        ggplot2::xlab(x_label) +
        ggplot2::ylab("") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 12)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::coord_cartesian(clip = "off")

    return(ridge_plot)
}