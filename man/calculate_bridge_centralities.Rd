% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_centralities.R
\name{calculate_bridge_centralities}
\alias{calculate_bridge_centralities}
\title{Calculate the bridge centralities of a GGM}
\usage{
calculate_bridge_centralities(ega)
}
\arguments{
\item{ega}{Output of construct_ggm}
}
\value{
A dataframe containing bridge centrality measurements for each node
}
\description{
This function takes in the output of construct_ggm and
returns a dataframe containing the bridge betweenness, bridge closeness,
and bridge strength for each node in the network.
}
