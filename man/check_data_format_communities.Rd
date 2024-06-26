% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_data_format.R
\name{check_data_format_communities}
\alias{check_data_format_communities}
\title{Check if a dataframe is formatted properly for community severity
analysis}
\usage{
check_data_format_communities(df)
}
\arguments{
\item{df}{Datframe representing output of get_communities}
}
\description{
This function checks if a dataframe is formatted properly for community
severity analysis. This is the format output by the get_communities
function.
}
\details{
Requirements include:
\enumerate{
\item First column must be named 'ID'
\item Second column must be named 'community'
\item All remaining columns must be numerical values
}
}
