% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_severity_frequency_occurrence.R
\name{plot_frequency}
\alias{plot_frequency}
\title{Plot bar chart of frequency}
\usage{
plot_frequency(df, color = "#4caf4a", y_label = "frequency")
}
\arguments{
\item{df}{Dataframe of symptom frequencies (formatted as required by PRONA)}

\item{color}{Color of the barplot (default green)}

\item{y_label}{Y-axis label}
}
\value{
A ggplot2 object of the frequency barplot
}
\description{
Same function as plot_severity, just with different default y-axis label
It is provided to make it clear to users that PRONA can be used to analyze
both severity and frequency data.
}
