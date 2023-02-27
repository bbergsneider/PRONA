#' Format PTSD data
#'
#' This function takes the qs.csv file available at
#' https://datashare.nida.nih.gov/study/nida-ctn-0015
#' and converts it to the format necessary for PRONA.
#' This function is provided as a reference for those
#' interested in directly downloading and working with
#' the raw data; however, running this function is not
#' necessary as the reformatted data is already provided
#' on the GitHub and within the installed package.
#'
#' These steps are adapted from the supplmentary materials
#' of Epskamp, Borsboom, & Fried (2018)
#' DOI: 10.3758/s13428-017-0862-1
#'
#' To download the raw data, go to:
#' https://datashare.nida.nih.gov/study/nida-ctn-0015
#' and click on 'CTN-0015 Data Files'. Fill out the
#' required data sharing agreement and navigate to the
#' file titled 'qs.csv'.
#'
#' @importFrom magrittr %>%
#' @param FullData Datframe representing qs.csv
#' @return A dataframe of the qs.csv data reformatted to match the PRONA format
#' @export
format_ptsd_data <- function(FullData){
    Data <- FullData %>%
        dplyr::filter(EPOCH == "BASELINE",
        grepl("^PSSR\\d+A$",QSTESTCD)) %>%
        dplyr::select(USUBJID,QSTEST,QSORRES) %>%
        tidyr::spread(QSTEST, QSORRES) %>%
        dplyr::rename("ID" = "USUBJID") %>%
        dplyr::mutate_each(dplyr::funs(replace(.,.=="NOT ANSWERED",NA))) %>%
        dplyr::mutate_each(dplyr::funs(replace(.,.=="NOT AT ALL",0))) %>%
        dplyr::mutate_each(dplyr::funs(replace(.,.=="ONCE A WEEK",1))) %>%
        dplyr::mutate_each(dplyr::funs(replace(.,.=="2-4 TIMES PER WEEK/HALF THE TIME",2))) %>%
        dplyr::mutate_each(dplyr::funs(replace(.,.=="5 OR MORE TIMES PER WEEK/ALMOST ALWAYS",3)))

    return(Data)
}