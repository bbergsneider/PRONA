% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootnet_analysis.R
\name{network_estimation_fun_normal_gamma0.25}
\alias{network_estimation_fun_normal_gamma0.25}
\title{Network estimation function that assumes all variables are normally distributed
and uses gamma = 0.25}
\usage{
network_estimation_fun_normal_gamma0.25(df)
}
\arguments{
\item{df}{Dataframe of symptom severity/frequency data}
}
\value{
The output of the EBICglasso function from qgraph, which is
a partial correlation matrix
}
\description{
This function is used by bootnet functions to specify how
the GGM should be calculated
}
