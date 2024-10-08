% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/construct_ggm.R
\name{construct_ggm}
\alias{construct_ggm}
\title{Construct a Gaussian Graphical Model (GGM)}
\usage{
construct_ggm(df, normal = FALSE, ...)
}
\arguments{
\item{df}{Dataframe of symptom severity/frequency data}

\item{normal}{Whether to consider all variables normally distrubted.
If false, conducts a non-paranormal transformation (Default: FALSE)}
}
\value{
The output of the EGA function from EGAnet
}
\description{
This function takes in a dataframe of symptom severity/frequency data
(in PRONA format) at creates a Gaussian Graphical Model (GGM) for it
using the EGA function from the R package EGAnet. Any unspecified parameters
are passed into the EGA function in EGAnet.
}
