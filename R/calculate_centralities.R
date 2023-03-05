#' Calculate the centralities of a GGM
#'
#' This function takes in the output of construct_ggm and
#' returns a dataframe containing the betweenness, closeness,
#' strength, and expected influence for each node in the network.
#'
#' @param ega Output of construct_ggm
#' @return A dataframe containing centrality measurements for each node
#' @export

calculate_centralities <- function(ega) {
    centrality_stats <- qgraph::centrality_auto(ega$network)
    centrality_table <- centrality_stats$node.centrality
    return(centrality_table)
}


#' Calculate the bridge centralities of a GGM
#'
#' This function takes in the output of construct_ggm and
#' returns a dataframe containing the bridge betweenness, bridge closeness,
#' and bridge strength for each node in the network.
#'
#' @param ega Output of construct_ggm
#' @return A dataframe containing bridge centrality measurements for each node
#' @export

calculate_bridge_centralities <- function(ega) {
    bridge_centralities <- networktools::bridge(ega$network, communities = ega$wc)
    bridge_centrality_table <- as.data.frame(list(bridge_centralities$'Bridge Strength', bridge_centralities$'Bridge Betweenness', bridge_centralities$'Bridge Closeness'), col.names = list('Bridge Strength','Bridge Betweenness','Bridge Closeness'))
    return(bridge_centrality_table)
}


#' Plot the centralities of a GGM
#'
#' This function takes in the output of construct_ggm and
#' plots the strength, closeness, and betweenness centrality
#' of the network.
#'
#' @importFrom magrittr %>%
#' @param ega Output of construct_ggm
#' @return Plot of strength, closeness, and betweenness centrality
#' @export

plot_centralities <- function(ega) {
    # Get centralities from the ega input
    centralities <- qgraph::centralityTable(ega$network) %>%
        dplyr::filter(measure %in% c("Strength", "Closeness", "Betweenness")) %>%
        dplyr::select(node, measure, value)

    return(make_centrality_graph(centralities, bridge = FALSE))
}


#' Plot the bridge centralities of a GGM
#'
#' This function takes in the output of construct_ggm and
#' plots the bridge strength, bridge closeness, and bridge
#' betweenness centrality of the network.
#'
#' @importFrom magrittr %>%
#' @param ega Output of construct_ggm
#' @return Plot of bridge strength, closeness, and betweenness centrality
#' @export

plot_bridge_centralities <- function(ega) {
    # Get table of bridge centralities
    bridge_centrality_table <- calculate_bridge_centralities(ega)

    # Calculate bridge centrality z-scores
    bridge_centrality_scaled <- apply(bridge_centrality_table, 2, scale) %>%
        as.data.frame()
    bridge_centrality_scaled$node <- rownames(bridge_centrality_table)
    bridge_centrality_scaled <- bridge_centrality_scaled %>%
        tidyr::pivot_longer(!node, names_to = "measure", values_to = "value")
    bridge_centrality_scaled <- bridge_centrality_scaled %>%
        dplyr::rowwise() %>%
        dplyr::mutate(measure = gsub("\\.", "\n", measure))

    return(make_centrality_graph(bridge_centrality_scaled, bridge = TRUE))
}


#' Plot centrality graph given a dataframe of centrality measurements
#'
#' This is an internal function that is used by plot_centralities and
#' plot_bridge centralities to make a plot of centrality measurements
#'
#' @param centralities A dataframe of centrality measurements
#' @param bridge Boolean. Whether to graph normal centrality (FALSE)
#' or bridge centrality (TRUE)
#' @return A plot of the centrality measurements

make_centrality_graph <- function(centralities, bridge) {
    # Order the variables by strength centrality
    if (!bridge) {
        filter_var <- "Strength"
    } else {
        filter_var <- "Bridge\nStrength"
    }
    orderbystrength <- centralities %>%
        dplyr::filter(measure == filter_var) %>%
        dplyr::arrange(value) %>%
        dplyr::select(node)

    # Order the entire centralities dataframe by strength centrality
    if (!bridge) {
        labels <- c("Strength", "Closeness", "Betweenness")
    } else {
        labels <- c("Bridge\nStrength", "Bridge\nCloseness", "Bridge\nBetweenness")
    }
    centralities <- centralities %>%
        dplyr::group_by(measure) %>%
        dplyr::mutate(node = factor(node, levels = orderbystrength$node), measure = factor(measure, levels = labels)) %>%
        dplyr::arrange(node)

    # Create a ggplot of the centralities
    if (!bridge) {
        centrality_colors <- c("Strength"="#e4211c", "Closeness"="#387db8", "Betweenness"="#ff7f00")
    } else {
        centrality_colors <- c("Bridge\nStrength"="#e4211c", "Bridge\nCloseness"="#387db8", "Bridge\nBetweenness"="#ff7f00")
    }
    p <- ggplot2::ggplot(centralities %>% dplyr::group_by(measure), ggplot2::aes(x=value, y=node, group=measure, color=measure)) + 
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_path(size = 0.7) +
        ggplot2::scale_color_manual(values = centrality_colors) +
        ggplot2::facet_wrap(~measure, ncol = 3) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::theme(axis.text.y=ggplot2::element_text(size=12), axis.text.x=ggplot2::element_text(size=12), strip.text=ggplot2::element_text(size=12), legend.position = "none")

    return(p)
}