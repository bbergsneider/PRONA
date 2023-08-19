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