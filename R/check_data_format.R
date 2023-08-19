#' Check if a dataframe is formatted according to PRONA specifications
#'
#' This function checks if a dataframe is formatted properly for PRONA
#' and throws an error if not
#'
#' Requirements include:
#' 1. First column must be named 'ID'
#' 2. All remaining columns must be numerical values
#' 3. There must not be any NA values in the dataframe
#'
#' @param df Datframe representing symptom severity/frequency
#' @export

check_data_format <- function(df) {
    if (colnames(df)[1] != "ID") { # Check if first column is named 'ID'
        stop("Data formatted improperly: Make sure first column is named 'ID'")
    } else if (!all(sapply(dplyr::select(df, -ID), is.numeric))) { # Check if all remaining columns are numeric values
        stop("Data formatted improperly: Not all symptom data is numeric")
    } else if (any(is.na(df))) { # Check if there are any NA values in the dataframe
        stop("Data formatted improperly: Make sure there are no NA values in the dataframe")
    }
}


#' Check if a dataframe is formatted properly for community severity
#' analysis
#'
#' This function checks if a dataframe is formatted properly for community
#' severity analysis. This is the format output by the get_communities
#' function.
#'
#' Requirements include:
#' 1. First column must be named 'ID'
#' 2. Second column must be named 'community'
#' 3. All remaining columns must be numerical values
#'
#' @param df Datframe representing output of get_communities
#' @export

check_data_format_communities <- function(df) {
    if (colnames(df)[1] != "ID") { # Check if first column is named 'ID'
        stop("Data formatted improperly: Make sure first column is named 'ID'")
    } else if (colnames(df)[2] != "community") { # Check if second column is named 'community'
        stop("Data formatted improperly: Make sure second column is named 'community'")
    } else if (!all(sapply(dplyr::select(df, -ID), is.numeric))) { # Check if all remaining columns are numeric values
        stop("Data formatted improperly: Not all symptom data is numeric")
    } else if (any(is.na(df))) { # Check if there are any NA values in the dataframe
        stop("Data formatted improperly: Make sure there are no NA values in the dataframe")
    }
}