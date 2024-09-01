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
#' This function takes in PRO symptom data and detects communities
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


#' Subsample patients and re-run clustering
#'
#' This function takes subsamples of the patient dataset and re-runs concordance
#' network clustering to assess the stability of patient clusters. Specifically,
#' it takes random subsamples of 100%, 99%, 95%, 90%, 80%, 70%, 60%, 50%, 40%,
#' 30%, 20%, and 10% of the patient dataset. It re-runs clustering and outputs
#' a dataframe containing the community membership of each patient in each of
#' the subsampling runs.
#'
#' @param data Dataframe of symptom severities (formatted as required by PRONA)
#' @param detectAlgo The community detection algorithm to use. Options include
#' FG (fast greedy), IM (infomap), LP (label propagation), LE (leading eigen),
#' LV (louvain), or WT (walktrap). (Default: WT)
#' @param simil_measure The similarity measure to use, either 'Euclidean'
#' or 'ARI' (Adjusted Rand Index) (Default: ARI)
#' @param simplify_graphs Boolean that indicates whether to simplify the graph
#' by removing multi-edges and loops (Default: TRUE)
#' @param sampling_rates Vector containing the percentages to subset the dataset
#' by. Default is c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1).
#' If you change this, the vector must start with 1 (meaning the first subsample
#' is 100% of the dataset) for downstream analysis to work.
#' @return A dataframe containing the community membership of each patient in
#' each of the subsampling runs.
#' @export

cluster_subsamples <- function(data, detectAlgo='WT', simil_measure='ARI', simplify_graphs=TRUE, sampling_rates=c(1, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)) {
    # Initialize 'subsample_communities_df' dataframe
    subsample_communities_df <- dplyr::select(data, ID)

    # Loop over each sampling rate
    for(i in seq_along(sampling_rates)) {
      set.seed(42)
      data_sample <- data[sample(1:nrow(data), size = nrow(data)*sampling_rates[i]), ] 
      clusters <- get_communities(data_sample, detectAlgo, simil_measure, simplify_graphs)[, c('ID','community')]
      
      # Rename the 'community' column to reflect the sampling rate
      colnames(clusters)[which(names(clusters) == "community")] <- paste0("community_", sampling_rates[i])
      
      # Merge the clusters with 'subsample_communities_df'
      subsample_communities_df <- merge(subsample_communities_df, clusters, by = "ID", all.x = TRUE)
    }

    return(subsample_communities_df)
}


#' Compute stability of a single cluster across multiple subsampling runs
#'
#' This is a helper function for compute_cluster_stabilities. It computes the
#' stability of a single cluster across different runs (iterations) of
#' concordance network-based clustering.
#'
#' @param cluster_IDs The IDs of the patients in the original cluster of
#' interest
#' @param cluster_df A dataframe of the cluster assignments across different
#' subsampling iterations for only the patients in the original cluster of
#' interest
#' @param cluster_cols The names of the columns for each clustering
#' iteration
#' @return A vector of the stability of the cluster of iterest across each
#' of the subsampling iterations
#' @export

compute_cluster_stability <- function(cluster_IDs, cluster_df, cluster_cols) {
  counts = sapply(cluster_cols, function(col) {
    same_cluster = cluster_df[cluster_df$ID %in% cluster_IDs, col]
    freq_clusters = table(same_cluster)
    return(max(freq_clusters[names(freq_clusters) != ""]))
  })
  return(counts / length(cluster_IDs))
}

#' Compute stabilities of multiple clusters across multiple subsampling runs
#'
#' This function computes the stabilities of the clusters in the original
#' clustering run across multiple subsampling runs. "Stability" is how often
#' data points that belong to the same cluster in the original data are still
#' clustered together in the subsamples (this is also known as Prediction
#' Strength). It is calculated by taking all the patients that belong to
#' one of the individual clusters, finding which communities those patients
#' have been grouped into in the new clustering iteration, and dividing the size
#' of the largest community in the new iteration by the size of original
#' community. For example, say clustering on a dataset of 1000 patients yields
#' three communities, one of 500 patients (community 1), one of 400 (community
#' 2), and one of 100 (community 3). We rerun clustering on a subsample of the
#' data, and out of the 500 patients that originally belonged to community 1,
#' 400 of them still belong to community 1,  50 belong to community 2, and 50
#' were randomly removed during subsampling. The "stability" of the
#' original community 1 for this subsampling iteration is 400/500 = 0.8
#' (aka 80%).
#'
#' @param subsample_communities_df The output of the cluster_subsamples function
#' @return A dataframe containing the stability of each cluster in each
#' subsampling iteration.
#' @export

compute_cluster_stabilities <- function(subsample_communities_df) {
  # Get the columns of the different subsampling runs
  cluster_cols = grep('community_', names(subsample_communities_df), value=TRUE)

  # Loop over each cluster in the original clustering
  original_clusters = unique(subsample_communities_df$community_1)

  # Calculate the stability of each cluster in each subsampling run
  stabilities = sapply(original_clusters, function(cluster) {
    cluster_IDs = subsample_communities_df$ID[subsample_communities_df$community_1 == cluster]
    stability = compute_cluster_stability(cluster_IDs, subsample_communities_df, cluster_cols)
  })

  # Convert stabilities to a dataframe
  stabilities_df <- as.data.frame(stabilities)
  colnames(stabilities_df) <- paste("Cluster", seq_len(ncol(stabilities_df)), "Stability")
  stabilities_df$SamplingRate <- gsub("community_", "", rownames(stabilities_df))
  stabilities_df <- dplyr::select(stabilities_df, SamplingRate, dplyr::everything())
  rownames(stabilities_df) <- NULL

  return(stabilities_df)
}


#' Plot cluster stabilities
#'
#' This function plots the results of compute_cluster_stabilities.
#'
#' @param stabilities_df The output of compute_cluster_stabilities
#' @return A ggplot object of a plot of the cluster stabilities
#' @export

plot_cluster_stabilities <- function(stabilities_df) {
  # Change dataframe to long format
  stability_df <- stabilities_df %>%
    dplyr::mutate(SamplingRate = as.numeric(SamplingRate)) %>%
    tidyr::gather("Cluster", "Stability", -SamplingRate) %>%
    dplyr::mutate(dplyr::across(c(SamplingRate, Stability), ~ifelse(is.finite(.), ., NA))) %>%
    na.omit()

  # Plot line graph of stabilities for each cluster
  p <- ggplot2::ggplot(stability_df, ggplot2::aes(x=SamplingRate, y=Stability, color=Cluster)) + 
       ggplot2::geom_point(size = 2) +
       ggplot2::geom_smooth(method = 'lm', linetype = 1, se = TRUE, size = 0.5) +
       ggplot2::theme_minimal() +
       ggplot2::xlab("Sampling Rate") +
       ggplot2::ylab("Stability") +
       ggplot2::theme(axis.text.y=ggplot2::element_text(size=10), 
                      axis.text.x=ggplot2::element_text(size=10),
                      legend.position = "none",
                      strip.text = ggplot2::element_text(size = 10)) +
       ggplot2::facet_wrap(~Cluster, ncol = length(unique(stability_df$Cluster)))

  return(p)
}