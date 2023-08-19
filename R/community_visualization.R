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
#' @param symptom_order A vector of symptom names indicating the order
#' to plot symptoms on the y-axis. This order will only show up if
#' cluster_rows = FALSE. The default order if cluster_rows = FALSE is
#' the order of column names in the dataframe (minus the ID and
#' community columns)
#' @return A heatmap of the symptom severities in each community
#' @export

plot_community_heatmap <- function(data, cluster_rows = TRUE, symptom_order = NULL) {
    # Check if data is formatted properly
    check_data_format_communities(data)

    # If symptom_order is NULL, set it to the order of column names in the
    # dataframe by default
    if (is.null(symptom_order)) {
        symptom_order <- colnames(data)[c(3:length(colnames(data)))]
    }

    # Splitting the data by community and then selecting the necessary columns
    community_symptom_data_list <- split(data, data$community)
    community_symptom_data_list <- lapply(community_symptom_data_list, function(x) dplyr::select(x, -ID, -community))

    # Combine all the data to find the global minimum and maximum values
    all_data_combined <- do.call(rbind, community_symptom_data_list)
    min_value <- min(all_data_combined, na.rm = TRUE)
    max_value <- max(all_data_combined, na.rm = TRUE)

    # Create the color function based on the minimum and maximum values
    col_fun = circlize::colorRamp2(c(min_value, max_value), c('white', 'red'))

    # Create a heatmap for each community
    heatmaps <- lapply(seq_along(community_symptom_data_list), function(i) {
        create_community_heatmap(community_symptom_data_list, community_symptom_data_list[[i]], i, col_fun, cluster_rows, symptom_order)
    })

    # Combine the heatmaps
    heatmap_plot <- Reduce("+", heatmaps)
    return(heatmap_plot)
}


#' Function for plotting heatmap for a single community
#'
#' Function for plotting a heatmap for a single community.
#' Helper function for plot_community_heatmap.
#'
#' @param community_symptom_data_list list of symptom data for each community
#' @param data symptom severity for single community
#' @param cluster_number number of community of interest
#' @param col_fun palette to use for heatmap
#' @param cluster_rows Boolean representing whether or not to cluster
#' the rows via heirarchical clustering
#' @param symptom_order A vector of symptom names indicating the order
#' to plot symptoms on the y-axis.
#' @return A heatmap of the symptom severities in one community
#' @export
#'
create_community_heatmap <- function(community_symptom_data_list, data, cluster_number, col_fun, cluster_rows, symptom_order) {
  # Order rows by symptom_order
  data_ordered <- order_rows_by_vector(data, order_vector = symptom_order)

  # Set up heatmap arguments
  heatmap_args <- list(
    data.matrix(t(data_ordered)),
    col = col_fun,
    name = paste("Cluster", cluster_number),
    column_title = paste('PC', cluster_number),
    show_heatmap_legend = (cluster_number == 1),
    cluster_rows = cluster_rows,
    column_names_gp = grid::gpar(fontsize = 0)
  )

  # Check if it's the last community, and if so, add the row_names_gp parameter
  if (cluster_number == length(community_symptom_data_list)) {
    heatmap_args$row_names_gp = grid::gpar(fontsize = 8)
  }

  do.call(ComplexHeatmap::Heatmap, heatmap_args)
}


#' Function for ordering rows of a dataframe
#'
#' Function for order rows of a dataframe based on an input vector.
#' This is a helper function create_community_heatmap
#'
#' @param data dataframe of symptom severity
#' @param order_vector vector of how you want symptoms ordered
#' @return dataframe with rows ordered as specified
#' @export
#'
order_rows_by_vector <- function(data, order_vector) {
    data_ordered = data[, order_vector]
    return(data_ordered)
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
        ggplot2::theme(legend.position = "bottom", text = ggplot2::element_text(size = 6)) +
        ggplot2::labs(x = "Symptom", y = "Severity") +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 60,vjust=1,hjust=1)) +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 6), legend.text = ggplot2::element_text(size = 6))
    
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