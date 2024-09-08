#' Plot heatmap of symptom severities in different communities
#'
#' This function takes in the output of the get_communities function
#' and plots a heatmap of the symptom severities in each community
#'
#' @param data The output of the get_communities function, which is
#' a dataframe of symptom severities with a column representing which
#' community each patient belongs to
#' @param cluster_rows Boolean representing whether or not to cluster
#' the rows via heirarchical clustering (Default: TRUE)
#' @param network An optional parameter that lets you pass in a network
#' object output by construct_ggm. If a network is passed into this function,
#' the heatmap will automatically order rows (symptoms) by the symptom
#' clusters identified in the network. This overrides the cluster_rows
#' parameter. (Default: NULL)
#' @param row_label_size The font size of the row labels (Default: 10)
#' @return A heatmap of the symptom severities in each community
#' @export

plot_community_heatmap <- function(data, cluster_rows = TRUE, network = NULL, row_label_size = 10) {
  # Check if data is formatted properly
  check_data_format_communities(data)
 
  m <- data.matrix(t(data[,3:ncol(data)]))
 
  # Combine all the data to find the global minimum and maximum values
  min_value <- min(m, na.rm = TRUE)
  max_value <- max(m, na.rm = TRUE)
  col_fun = circlize::colorRamp2(c(min_value, max_value), c('white', 'red'))
 
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  communities <- paste("PC", levels(as.factor(unique(data$community))))
  ncommunities <- length(communities)
  community_colors <- gg_color_hue(ncommunities+1)[2:(ncommunities+1)]
  names(community_colors) <- communities
  column_ha = ComplexHeatmap::HeatmapAnnotation(`PC` = paste("PC", data$community), col=list(`PC`=community_colors), show_legend = FALSE)
 
  if (is.null(network)){
    heatmap_plot <- ComplexHeatmap::Heatmap(m, col=col_fun,
                                            column_split=paste("PC", data$community),
                                            cluster_column_slices = FALSE,
                                            name="rating", show_column_dend = FALSE,
                                            top_annotation = column_ha,
                                            row_names_gp = grid::gpar(fontsize=row_label_size),
                                            cluster_rows = cluster_rows)
  } else {
    row_community <- network$wc
    row_communities <- paste("SC", levels(as.factor(unique(network$wc))))
    row_community_colors <- gg_color_hue(length(row_communities))
    names(row_community_colors) <- row_communities
    symptom_order <- rownames(m)[order(row_community)]
    row_ha = ComplexHeatmap::rowAnnotation(`SC` = paste("SC", row_community[order(row_community)]), col=list(`SC`=row_community_colors), show_legend = FALSE)
   
    heatmap_plot <- ComplexHeatmap::Heatmap(m[symptom_order, ], col=col_fun, column_split=paste("PC", data$community),
                               cluster_column_slices = FALSE, name="rating", show_column_dend = FALSE,
                               top_annotation = column_ha, cluster_rows=FALSE, right_annotation = row_ha, row_names_gp = grid::gpar(fontsize=row_label_size))
   
  }
  return(ComplexHeatmap::draw(heatmap_plot, heatmap_legend_side = "left"))
}


#' Get a summary of each symptom in each community
#'
#' Function that takes the output of get_communities and
#' returns a dataframe with the mean, median, upper bound,
#' and lower bound of each symptom in each community
#'
#' @param data output of get_communities
#' @return dataframe with summary statistics for each symptom
#' in each community
#' @export
#'
get_community_summary <- function(data) {
    # Check if data is formatted properly
    check_data_format_communities(data)
    
    # Bind all community symptom data frames together in long format
    all_comms_long_symptom_data <- reshape2::melt(data %>% dplyr::select(-ID), id='community')

    # Make dataframe summarizing mean, median, and 95% CI of each variable in each community
    severity_summary_df <- all_comms_long_symptom_data %>%
        dplyr::group_by(community, variable) %>%
        dplyr::summarize(
        mean = mean(value, na.rm = TRUE),
        median = median(value, na.rm = TRUE),
        upper = mean + 1.96 * sd(value, na.rm = TRUE)/sqrt(sum(!is.na(value))),
        lower = mean - 1.96 * sd(value, na.rm = TRUE)/sqrt(sum(!is.na(value)))
        )

    return(severity_summary_df)
}


#' Plot a line plot of the severities for the communities you specify
#'
#' This function plots a line plot of avergae severity (with 95% CIs) for the
#' communities you specify
#'
#' @param data Output of get_communities
#' @param communities Vector specifying which communities to plot. If you
#' include 0 in the vector, it will plot the severity for all patients
#' combined. (Default: c(0))
#' @return A line plot of symptom severities in each community
#' @export
#'
community_line_plot <- function(data, communities = c(0)) {
    # Check if data is formatted properly
    check_data_format_communities(data)

    # Get summary dataframe for data
    summary_df <- get_community_summary(data)

    # Order symptoms by mean expression in all patients
    col_order <- order_cols_by_mean(data %>% dplyr::select(-ID,-community))

    # If 0 is included in communities, add community 0 to severity_summary_df
    # that represents summary stats for entire dataset
    if (0 %in% communities) {
        all_comms_long_symptom_data <- reshape2::melt(data %>% dplyr::select(-ID), id='community')
        combined_severity_summary_df <- all_comms_long_symptom_data %>%
            dplyr::group_by(variable) %>%
            dplyr::summarize(
            mean = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            upper = mean + 1.96 * sd(value, na.rm = TRUE)/sqrt(sum(!is.na(value))),
            lower = mean - 1.96 * sd(value, na.rm = TRUE)/sqrt(sum(!is.na(value)))
            )
        combined_severity_summary_df <- cbind(community=0,combined_severity_summary_df)
        summary_df <- rbind(combined_severity_summary_df, summary_df)
    }

    # Make plot
    p <- ggplot2::ggplot(summary_df[summary_df$community %in% communities,], ggplot2::aes(x = factor(variable, level=col_order), y = mean, colour = factor(community), group = factor(community))) + 
        ggplot2::geom_line(size=1) +
        ggplot2::geom_point(ggplot2::aes(shape = factor(community)), size = 3) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2, size = 1) +
        ggplot2::theme_classic() +
        ggplot2::scale_color_discrete(labels = communities) +
        ggplot2::scale_shape_discrete(labels = communities) +
        ggplot2::theme(legend.position = "bottom", text = ggplot2::element_text(size = 8)) +
        ggplot2::labs(x = "Symptom", y = "Severity") +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 60,vjust=1,hjust=1)) +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 8), legend.text = ggplot2::element_text(size = 8))
    
    return(p)
}


#' Order columns of a dataframe from highest to lowest mean value
#'
#' This function aims to order the columns of a symptom severity dataframe
#' from highest to lowest mean severity. It is a helper function for
#' community_line_plot
#'
#' @param data symptom severity data (without ID and community columns)
#' @return List of column names, ordered by mean severity
#' @export
#'
order_cols_by_mean <- function(df){
  # calculate mean of each column
  mean_cols <- colMeans(df, na.rm = TRUE)
  
  # sort column means in descending order
  mean_cols <- sort(mean_cols, decreasing = FALSE)
  
  # match column means to column names
  col_order <- names(mean_cols)
  
  # return column order
  return(col_order)
}