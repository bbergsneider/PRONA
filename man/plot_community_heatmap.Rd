% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/community_visualization.R
\name{plot_community_heatmap}
\alias{plot_community_heatmap}
\title{Plot heatmap of symptom severities in different communities}
\usage{
plot_community_heatmap(
  data,
  cluster_rows = TRUE,
  network = NULL,
  row_label_size = 10
)
}
\arguments{
\item{data}{The output of the get_communities function, which is
a dataframe of symptom severities with a column representing which
community each patient belongs to}

\item{cluster_rows}{Boolean representing whether or not to cluster
the rows via heirarchical clustering (Default: TRUE)}

\item{network}{An optional parameter that lets you pass in a network
object output by construct_ggm. If a network is passed into this function,
the heatmap will automatically order rows (symptoms) by the symptom
clusters identified in the network. This overrides the cluster_rows
parameter. (Default: NULL)}

\item{row_label_size}{The font size of the row labels (Default: 10)}
}
\value{
A heatmap of the symptom severities in each community
}
\description{
This function takes in the output of the get_communities function
and plots a heatmap of the symptom severities in each community
}
