#' Check if a dataframe is formatted according to PRONA specifications
#'
#' This function checks if a dataframe is formatted properly for PRONA
#' and throws an error if not
#'
#' Requirements include:
#' 1. First column must be named 'ID'
#' 2. All remaining columns must be numerical values
#'
#' @param df Datframe representing symptom severity/frequency
#' @export

check_data_format <- function(df) {
    if (colnames(df)[1] != "ID") { # Check if first column is named 'ID'
        stop("Data formatted improperly: Make sure first column is named 'ID'")
    } else if (!all(sapply(dplyr::select(df, -ID), is.numeric))) { # Check if all remaining columns are numeric values
        stop("Data formatted improperly: Not all symptom data is numeric")
    }
}