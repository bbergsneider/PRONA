#' Construct a Gaussian Graphical Model (GGM)
#'
#' This function takes in a dataframe of symptom severity/frequency data
#' (in PRONA format) at creates a Gaussian Graphical Model (GGM) for it
#' using the EGA function from the R package EGAnet.
#'
#' @param df Dataframe of symptom severity/frequency data
#' @param normal Whether to consider all variables normally distrubted.
#' If false, conducts a non-paranormal transformation (Default: FALSE)
#' @return The output of the EGA function from EGAnet
#' @export

construct_ggm <- function(df, normal = FALSE) {
    check_data_format(df) # Check if df is formatted properly
    df <- dplyr::select(df, -ID) # Remove ID column

    # Calculate correlation matrix
    cor_matrix <- qgraph::cor_auto(df, detectOrdinal = TRUE, npn.SKEPTIC = !normal, forcePD = TRUE)

    # Construct an EBIC-GLASSO regularized GGM using the EGA function from EGAnet
    # EGA also identifies symptom clusters using the walktrap algorithm
    ega <- EGAnet::EGA(cor_matrix, n = nrow(df), plot.EGA = FALSE)

    return(ega)
}


#' Plot a GGM
#'
#' This function takes the output from construct_ggm and plots it
#'
#' @param ega Output from construct_ggm
#' @param colors A list of colors to color the clusters with (optional)
#' @param legend.names A list of strings to label the clusters with (optional)
#' @return A plot of the GGM
#' @export

plot_ggm <- function(ega, colors = NULL, legend.names = NULL) {
    if (is.null(colors)) {
        if (is.null(legend.names)) {
            ega_plot <- plot(ega, plot.args = list(node.size = 8, label.size = 3.5))
        } else {
            ega_plot <- plot(ega, plot.args = list(node.size = 8, label.size = 3.5, legend.names = legend.names))
        }
    } else {
        if (is.null(legend.names)) {
            ega_plot <- plot(ega, plot.args = list(node.size = 8, label.size = 3.5, color.palette = colors))
        } else {
            ega_plot <- plot(ega, plot.args = list(node.size = 8, label.size = 3.5, color.palette = colors, legend.names = legend.names))
        }
    }
    return(ega_plot)
}



#' Returns a dataframe of the weights of each edge in a GGM
#'
#' This function takes the output from construct_ggm and returns
#' a dataframe of the weights on each edge in the GGM
#'
#' @param ega Output from construct_ggm
#' @return A dataframe of the weights of each edge in the GGM
#' @export

get_ggm_weights <- function(ega) {
    m <- ega$network
    df <- data.frame(row = rownames(m)[row(m)[upper.tri(m)]],
           col = colnames(m)[col(m)[upper.tri(m)]],
           weight = m[upper.tri(m)])
    return(df)
}