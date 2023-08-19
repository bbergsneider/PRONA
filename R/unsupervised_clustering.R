#' Calculate the concordance matrix for a vector
#'
#' This function takes a vector (vec) and computes a concordance
#' matrix using the formula (vec*vec^T)/100. The result of this
#' function is a square matrix that represents the concordance
#' between elements of the vector. This is a helper function for
#' the vec2concords function
#'
#' @param vec input vector
#' @export

vec2concord <- function(vec){
  (vec %*% t(vec))/100
}


#' Calculate a concordance matrix for each row in a matrix
#'
#' This function takes a matrix (mat) where each row represents
#' a vector. It then applies the vec2concord function to each row
#' in the matrix, resulting in a list of concordance matrices for
#' each row in mat. This is a helper function for the similarity_matrix
#' function
#'
#' @param mat input matrix
#' @export

vecs2concords <- function(mat) {
  # aggregate concordance matrices in a list
  lapply(1:nrow(mat), function(i) vec2concord(mat[i,]))
}


#' Compute a similarity matrix based on an input matrix
#'
#' This function computes a similarity matrix based on an input
#' matrix (mat) using either Euclidean disease or the Adjusted
#' Rand Index (ARI) as the similarity measure. See Henry et al.,
#' Plos One, 2018 for more details. This is a helper function
#' for the get_communities function.
#'
#' @param mat The input matrix for which the similarity matrix
#' will be calculated
#' @param concords By default, this argument takes the value returned by
#' the vec2concords(mat) function, which calculates the concordance
#' matrices for each row of mat
#' @param simil_measure The similarity measure to use, either 'Euclidean'
#' or 'ARI' (Adjusted Rand Index) (Default: ARI)
#' @export

similarity_matrix <- function(mat, concords=vecs2concords(mat), simil_measure='ARI') {
  
  # Compute similarities
  sims <- NULL

  if (simil_measure == 'Euclidean') {
    a.conc <- t(apply(mat, 1, function(x) (x%*%t(x))[upper.tri(x%*%t(x), diag=FALSE)])) # taking just upper triangle w/o diagonal
    dists <- as.matrix(dist(a.conc))
    sims <- sqrt(1/(1+dists))
  } else {
    sims <- matrix(0, nrow=length(concords), ncol=length(concords))

    # vectorize concords
    vConcords <- lapply(1:length(concords), function(i) as.vector(concords[[i]][upper.tri(concords[[i]], diag = F)]))

    # compute similarity matrix by calculating pairwise ARI
    # (upper triangular part - not including diagonal)
    for (i in 1:(nrow(sims)-1)) {
      for (j in (i+1):ncol(sims)) {
        sims[i,j] <- mclust::adjustedRandIndex(vConcords[[i]], vConcords[[j]])
      }
    }

    # reflect upper triangle across diagonal (cause symmetric)
    sims <- sims + t(sims)
    
    # diagonal (all 1s)
    for (i in 1:nrow(mat)) {
      sims[i,i] <- 1
    }
  }

  # set negative entries to 0
  sims[sims < 0] <- 0

  sims
}


#' Detect patient communities based on similarity of symptom patterns
#'
#' This function takes in PRO symtom data and detects communities
#' of patients with shared symptom patterns based on concordance network
#' clustering. See Henry et al., Plos One, 2018 for more details on
#' methodology.
#'
#' @param data Dataframe of symptom severities (formatted as required by PRONA)
#' @param nrows Number of rows in the dataframe (Defaults to the number of rows
#' in data)
#' @param ncols Number of columns in the dataframe (Defaults to the number of
#' columns in data)
#' @param detectAlgo The community detection algorithm to use. Options include
#' FG (fast greedy), IM (infomap), LP (label propagation), LE (leading eigen),
#' LV (louvain), or WT (walktrap). (Default: WT)
#' @param simil_measure The similarity measure to use, either 'Euclidean'
#' or 'ARI' (Adjusted Rand Index) (Default: ARI)
#' @param simplify_graphs Boolean that indicates whether to simplify the graph
#' by removing multi-edges and loops (Default: TRUE)
#' @return A new dataframe that includes community grouping information along
#' with original symptom severity data
#' @export

get_communities <- function(data, nrows=nrow(data), ncols=ncol(data), detectAlgo='WT', simil_measure='ARI', simplify_graphs=TRUE) {
   # Check if data is formatted properly
  check_data_format(data)
  output_data <- data
  data <- dplyr::select(data, -ID)
  
  # convert dataframe to matrix
  mat <- as.matrix(data)

  # if matrix is empty throw error
  if (nrow(mat) == 0){
    stop('Matrix was empty.')
  }
  
  # ensure reproducibility
  set.seed(42)
  
  # fix row names
  rownames(mat) <- NULL
  
  # compute between patient similarity (generate similarity matrix)
  concords <- vecs2concords(mat)
  sims <- similarity_matrix(mat, concords, simil_measure=simil_measure)
  
  # convert sims to a graph to prepare for random walk
  g <- igraph::graph_from_adjacency_matrix(sims, mode='undirected', weighted=TRUE, diag=FALSE)
  
  # remove multi-edges and loops
  if (simplify_graphs) {
    g <- igraph::simplify(g)
  }
  
  # run community detection algorithm to get communities (the default is Random Walk (walktrap))
  comms <- NULL
  
  if (detectAlgo == 'FG'){
    comms <- igraph::cluster_fast_greedy(g)
  } else if (detectAlgo == 'IM'){
    comms <- igraph::cluster_infomap(g)
  } else if (detectAlgo == 'LP'){
    comms <- igraph::cluster_label_prop(g)
  } else if (detectAlgo == 'LE'){
    comms <- igraph::cluster_leading_eigen(g)
  } else if (detectAlgo == 'LV'){
    comms <- igraph::cluster_louvain(g)
  } else {
    comms <- igraph::cluster_walktrap(g)
  }
  
  # convert comms to regular r list
  comms <- lapply(igraph::groups(comms), function(x) as.numeric(x))
  
  # add community IDs to symptom severity output dataframe
  output_data <- output_data %>% tibble::add_column(community=0, .after='ID')
  for (i in c(1:length(comms))) {
    output_data[comms[[i]],]['community'] = i
  }

  return(output_data)
}


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
    dplyr::group_by(community,variable) %>%
    dplyr::summarize(mean=mean(value), median=median(value), upper = mean + 1.96 * sd(value)/sqrt(length(value)), lower = mean - 1.96 * sd(value)/sqrt(length(value)))

    return(severity_summary_df)
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