#' Perform BootEGA
#'
#' This function analyzes symptom cluster stability using the boot.ega function
#' from EGAnet
#'
#' @param df Dataframe of symptom severity/frequency data
#' @param normal Whether to consider all variables normally distrubted.
#' If false, conducts a non-paranormal transformation (Default: FALSE)
#' @param iter Number of permutations for bootstrap analysis (Default: 10,000)
#' @return The output of the BootEGA function from EGAnet
#' @export

run_bootEGA <- function(df, normal = FALSE, iter = 10000) {
    check_data_format(df) # Check if df is formatted properly
    df <- dplyr::select(df, -ID) # Remove ID column
    
    # Calculate correlation matrix
    cor_matrix <- qgraph::cor_auto(df, detectOrdinal = FALSE, npn.SKEPTIC = !normal, forcePD = TRUE)
    
    # Analyze stability of dimensions using parametric bootstrapping with 10,000 iterations
    boot.ega <- EGAnet::bootEGA(cor_matrix, n = nrow(df), iter = iter, type = 'parametric')
    return(boot.ega)
}



#' Plot median network structure from bootEGA
#'
#' This function plots the median determined by BootEGA stability analysis
#'
#' @param boot.ega Output of run_bootEGA function
#' @return A plot of the median GGM
#' @export

plot_bootEGA <- function(boot.ega) {
    median_boot_ega <- plot(boot.ega, plot.args = list(node.size = 8, label.size = 4))
    return(median_boot_ega)
}


#' Calculate item and dimension stability
#'
#' This function calculates item and dimension stability
#' using the EGAnet function dimensionStability
#'
#' @param boot.ega Output of run_bootEGA function
#' @return Output of EGAnet dimensionStability function
#' @export

calculate_dimStab <- function(boot.ega) {
    dimStab <- EGAnet::dimensionStability(boot.ega)
    return(dimStab)
}


#' Plot item stability
#'
#' This function plots item stability using the output of the
#' calculate_dimStab function
#'
#' @param dimStab Output of calculate_dimStab function
#' @param colors A list of colors to color the dimension stability
#' plot with (optional)
#' @return A graph of item stability values
#' @export

plot_item_stability <- function(dimStab, colors = NULL) {
    if (!is.null(colors)) {
        stab_plot <- dimStab$item.stability$plot +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values=colors) +
            ggplot2::ylab("replication proportion") + 
            ggplot2::xlab("") + 
            ggplot2::theme(axis.text.y = ggplot2::element_text(size=11),
                    axis.text.x = ggplot2::element_text(size=11),
                    strip.text = ggplot2::element_text(size=11), 
                    axis.title.y = ggplot2::element_text(size=11),
                    legend.position = "none")
    } else {
        stab_plot <- dimStab$item.stability$plot +
            ggplot2::theme_minimal() +
            ggplot2::ylab("replication proportion") + 
            ggplot2::xlab("") + 
            ggplot2::theme(axis.text.y = ggplot2::element_text(size=11),
                    axis.text.x = ggplot2::element_text(size=11),
                    strip.text = ggplot2::element_text(size=11), 
                    axis.title.y = ggplot2::element_text(size=11),
                    legend.position = "none")
    }

    return(stab_plot)
}