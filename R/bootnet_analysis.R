#' Network estimation function that assumes all variables are normally distributed
#'
#' This function is used by bootnet functions to specify how
#' the GGM should be calculated
#'
#' @param df Dataframe of symptom severity/frequency data
#' @return The output of the EBICglasso function from qgraph, which is
#' a partial correlation matrix
#' @export

network_estimation_fun_normal <- function(df) {
  check_data_format(df) # Check if df is formatted properly
  df <- dplyr::select(df, -ID) # Remove ID column

  # Calculate correlation matrix
  cor_matrix <- qgraph::cor_auto(df, detectOrdinal=FALSE, forcePD=TRUE, npn.SKEPTIC = FALSE)
  # Perform EBIC glasso
  EBIC_matrix <- qgraph::EBICglasso(cor_matrix, nrow(df), 0.5)

  return(EBIC_matrix)
}

#' Network estimation function that assumes all variables are non-normally distributed
#'
#' This function is used by bootnet functions to specify how
#' the GGM should be calculated
#'
#' @param df Dataframe of symptom severity/frequency data
#' @return The output of the EBICglasso function from qgraph, which is
#' a partial correlation matrix
#' @export

network_estimation_fun_non_normal <- function(df) {
  check_data_format(df) # Check if df is formatted properly
  df <- dplyr::select(df, -ID) # Remove ID column

  # Calculate correlation matrix
  cor_matrix <- qgraph::cor_auto(df, detectOrdinal=FALSE, forcePD=TRUE, npn.SKEPTIC = TRUE)
  # Perform EBIC glasso
  EBIC_matrix <- qgraph::EBICglasso(cor_matrix, nrow(df), 0.5)

  return(EBIC_matrix)
}


#' Nonparametric bootstrapping
#'
#' This function performs nonparametric bootstrapping to analyze edge weight
#' accuracy.
#'
#' @param df Dataframe of symptom severity/frequency data
#' @param normal Whether to consider all variables normally distrubted.
#' If false, conducts a non-paranormal transformation (Default: FALSE)
#' @param nBoots Number of bootstrapping iterations to perform (Default: 2500)
#' @param statistics Vector indicating which statistics to store. See bootnet
#' for options. (Default: c('edge','strength','closeness','betweenness')
#' @param communities If you are running bootnet on bridge centrality measures, use
#' this parameter to set community labels (set this to ega$wc, where ega is the
#' output from the construct_ggm function)
#' @return The output of the bootnet function when run with nonparametric
#' bootstrapping
#' @export

nonparam_boot <- function(df, normal = FALSE, nBoots = 2500, statistics = c('edge','strength','closeness','betweenness'), communities = NULL) {
    if (normal) {
        nonparam.boot <-bootnet::bootnet(df, nBoots=nBoots, statistics=statistics, fun='network_estimation_fun_normal', communities=communities)
    } else {
        nonparam.boot <-bootnet::bootnet(df, nBoots=nBoots, statistics=statistics, fun='network_estimation_fun_non_normal', communities=communities)
    }
    return(nonparam.boot)
}


#' Plot nonparametric bootstrapping
#'
#' This function performs plots the results of nonparametric bootstrapping
#' to analyze edge weight accuracy
#'
#' @param nonparam.boot Output of nonparam_boot function
#' @return A graph of edge weight accuracy
#' @export

plot_nonparam_boot <- function(nonparam.boot) {
    p <- plot(nonparam.boot, labels = FALSE, order = "sample", sampleColor=red, bootColor=blue, meanColor=blue, meanlwd=0.4, bootAlpha=0.2, bootlwd=0.4, areaAlpha=0.5)
    pp <- p +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", axis.text.y = ggplot2::element_blank(), 
                                  legend.text = ggplot2::element_text(size=12),
                                  axis.text.x = ggplot2::element_text(size=12),
                                  strip.text = ggplot2::element_text(size=12), 
                                  panel.grid.major.y = ggplot2::element_blank(),
                                  panel.grid.minor.y = ggplot2::element_blank())
    return(pp)
}


#' Get nonparam boot summary table
#'
#' This function returns a summary table of the object that
#' is outputted from the nonparam_boot function 
#'
#' @param nonparam.boot Output of nonparam_boot function
#' @return A dataframe summarizing the nonparametric bootstrapping
#' results
#' @export

summarize_nonparam_boot <- function(nonparam.boot) {
    summary_table <- summary(nonparam.boot, statistics = "edge")
    return(summary_table)
}


#' Casedrop bootstrapping
#'
#' This function performs casedrop bootstrapping to analyze edge weight
#' and centrality measure stability
#'
#' @param df Dataframe of symptom severity/frequency data
#' @param normal Whether to consider all variables normally distrubted.
#' If false, conducts a non-paranormal transformation (Default: FALSE)
#' @param nBoots Number of bootstrapping iterations to perform (Default: 2500)
#' @param statistics Vector indicating which statistics to store. See bootnet
#' for options. (Default: c('edge','strength','closeness','betweenness')
#' @param communities If you are running bootnet on bridge centrality measures, use
#' this parameter to set community labels (set this to ega$wc, where ega is the
#' output from the construct_ggm function)
#' @return The output of the bootnet function when run with casedrop
#' bootstrapping
#' @export

casedrop_boot <- function(df, normal = FALSE, nBoots = 2500, statistics = c('edge','strength','closeness','betweenness'), communities = NULL) {
    if (normal) {
        casedrop.boot <-bootnet::bootnet(df, nBoots=nBoots, type="case", statistics=statistics, fun='network_estimation_fun_normal', communities=communities)
    } else {
        casedrop.boot <-bootnet::bootnet(df, nBoots=nBoots, type="case", statistics=statistics, fun='network_estimation_fun_non_normal', communities=communities)
    }
    return(casedrop.boot)
}


#' Plot casedrop bootstrapping results
#'
#' This function performs plots the results of casedrop bootstrapping
#' to analyze edge weight and centrality measure stability
#'
#' @param casedrop.boot Output of casedrop_boot function
#' @param type String indicating which measurements to plot. Options are
#' "edge" (which plots only casedrop results for edge weight), "centralities"
#' (which plots the casedrop results for strength, closeness, and betweenness),
#' and "bridge centralities" (which plots the casedrop results for bridge
#' strength, bridge closeness, and bridge betweenness) (Default: "centralities")
#' @return A graph of edge weight/centrality measure stability
#' @export

plot_casedrop_boot <- function(casedrop.boot, type = 'centralities') {
    # Set color palette for plot
    orange="#ff7f00"
    red="#e4211c"
    blue="#387db8"

    edge_colors <- c("edge"=red)
    centrality_colors <- c("strength"=red, "closeness"=blue, "betweenness"=orange)
    bridge_centrality_colors <- c("bridgeStrength"=red, "bridgeCloseness"=blue, "bridgeBetweenness"=orange)
    
    # Set color and statistics parameters for plotting depending on type of plot
    if (type == 'edge'){
        colors_use <- edge_colors
        statistics = c('edge')
    } else if (type == 'centralities') {
        colors_use <- centrality_colors
        statistics = c('strength','closeness','betweenness')
    } else if (type == 'bridge centralities') {
        colors_use <- bridge_centrality_colors
        statistics = c('bridgeStrength', 'bridgeCloseness', 'bridgeBetweenness')
    } else {
        stop("statistics parameter must be 'edge', 'centralities', or 'bridge centralities'")
    }

    # Plot casedrop bootstrapping results
    p <- plot(casedrop.boot, statistics=statistics) +
        ggplot2::theme_minimal() +
        ggplot2::scale_color_manual(values=colors_use, labels=names(colors_use)) +
        ggplot2::scale_fill_manual(values=colors_use, labels=names(colors_use)) +
        ggplot2::xlab("% sampled cases") +
        ggplot2::ylab("correlation with original") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size=11), axis.text.x = ggplot2::element_text(size=10), legend.position = "top") +
        ggplot2::guides(fill = ggplot2::guide_legend(title=""), color = ggplot2::guide_legend(title=""))
    
    return(p)
}


#' Perform correlation stability analysis
#'
#' This function uses the corStability function from bootnet
#' to perform correlation stability analysis
#'
#' @param casedrop.boot Output of casedrop_boot function
#' @return Output of bootnet corStability function
#' @export

cor_stability_analysis <- function(casedrop.boot) {
    return(bootnet::corStability(casedrop.boot))
}


#' Nodedrop bootstrapping
#'
#' This function performs nodedrop bootstrapping to analyze edge weight
#' and centrality measure stability
#'
#' @param df Dataframe of symptom severity/frequency data
#' @param normal Whether to consider all variables normally distrubted.
#' If false, conducts a non-paranormal transformation (Default: FALSE)
#' @param nBoots Number of bootstrapping iterations to perform (Default: 2500)
#' @param statistics Vector indicating which statistics to store. See bootnet
#' for options. (Default: c('edge','strength','closeness','betweenness')
#' @return The output of the bootnet function when run with nodedrop
#' bootstrapping
#' @export

nodedrop_boot <- function(df, normal = FALSE, nBoots = 2500, statistics = c('edge','strength','closeness','betweenness')) {
    if (normal) {
        nodedrop.boot <-bootnet::bootnet(df, nBoots=nBoots, type="node", statistics=statistics, fun='network_estimation_fun_normal')
    } else {
        nodedrop.boot <-bootnet::bootnet(df, nBoots=nBoots, type="node", statistics=statistics, fun='network_estimation_fun_non_normal')
    }
    return(nodedrop.boot)
}


#' Plot nodedrop bootstrapping results
#'
#' This function performs plots the results of nodedrop bootstrapping
#' to analyze edge weight and centrality measurement stability
#'
#' @param nodedrop.boot Output of casedrop_boot function
#' @param type String indicating which measurements to plot. Options are
#' "edge" (which plots only casedrop results for edge weight) or "centralities"
#' (which plots the casedrop results for strength, closeness, and betweenness)
#' @return A graph of edge weight/centrality measure stability
#' @export

plot_nodedrop_boot <- function(nodedrop.boot, type = 'centralities') {
    # Set color palette for plot
    orange="#ff7f00"
    red="#e4211c"
    blue="#387db8"

    edge_colors <- c("edge"=red)
    centrality_colors <- c("strength"=red, "closeness"=blue, "betweenness"=orange)
    
    # Set color and statistics parameters for plotting depending on type of plot
    if (type == 'edge'){
        colors_use <- edge_colors
        statistics = c('edge')
    } else if (type == 'centralities') {
        colors_use <- centrality_colors
        statistics = c('strength','closeness','betweenness')
    } else {
        stop("statistics parameter must be either 'edge' or 'centralities'")
    }

    # Plot nodedrop bootstrapping results
    p <- plot(nodedrop.boot, statistics=statistics) +
        ggplot2::theme_minimal() +
        ggplot2::scale_color_manual(values=colors_use, labels=names(colors_use)) +
        ggplot2::scale_fill_manual(values=colors_use, labels=names(colors_use)) +
        ggplot2::xlab("% sampled nodes") +
        ggplot2::ylab("correlation with original") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size=11), axis.text.x = ggplot2::element_text(size=10), legend.position = "top") +
        ggplot2::guides(fill = ggplot2::guide_legend(title=""), color = ggplot2::guide_legend(title=""))
    
    return(p)
}