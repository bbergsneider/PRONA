#' Perform Unique Variable Analysis (UVA)
#'
#' This function performs Unique Variable Analysis (UVA) on a
#' dataframe of symptom severity/frequency data (in PRONA format).
#' The function allows the user to manually choose which variables
#' to consolidate
#'
#' It saves the output from UVA in a .RDS file, and the new dataframe
#' with consolidated latent variables in a .csv file
#'
#' For more information on UVA, see Christensen, Garrido, & Colino
#' (PsyArXiv, doi: 10.31234/osf.io/4kra2)
#'
#' @param df Dataframe representing symptom severity/frequency
#' @param reduce Boolean. Should redundancy reduction be performed?
#' (default: FALSE)
#' @param scale The scale on which symptom severity/frequency is measured.
#' For example if symptoms are measured on a 0-10 scale, scale = 10.
#' (default: 1) (only applicable if reduce = TRUE)
#' @param auto Boolean. Should redundancy reduction be automated?
#' Defaults to FALSE for manual selection.
#' @param label.latent 	Boolean. Should latent variables be manually labelled?
#' Defaults to TRUE. Set to FALSE for arbitrary labelling (i.e., "LV_").
#' @param output_dir Directory in which the .RDS and .csv outputs of UVA should
#' be saved (only applicable if reduce = TRUE) (defaults to current directory)
#' @return The list of the outputs from the UVA function as implemented
#' in the R package EGAnet
#' @export

perform_uva <- function(df, reduce = FALSE, scale = 1, auto = FALSE, label.latent = FALSE,output_dir = "") {
    check_data_format(df) # Check if df is formatted properly
    ids <- df["ID"] # Save ID column
    df <- dplyr::select(df, -ID) # Get rid of ID column

    # Conduct UVA
    uva_output <- EGAnet::UVA(data = df, reduce = reduce, auto = auto, label.latent = label.latent, adhoc = FALSE)

    if (reduce) {
        # Save reduced data in a new dataframe
        reduced_df <- uva_output$reduced$data
        # Re-add ID column to reduced data
        reduced_df <- cbind(ids, reduced_df)

        # Latent variable reduction yields variable values that are not
        # on the same scale as the original data. In order to keep latent
        # variable values consistent with other variables, latent variable
        # values must be normalized to be on the same scale
        for (lv in rownames(uva_output$reduced$merged)) {
            reduced_df[,lv] <- (reduced_df[,lv] - min(reduced_df[,lv])) / (max(reduced_df[,lv]) - min(reduced_df[,lv])) * scale
        }

        # Save uva_output and reduced_df in output files
        if (output_dir == "") {
            saveRDS(uva_output, "uvaResults.RDS")
            write.csv(reduced_df, "reduced_data.csv", row.names = FALSE)
        } else {
            saveRDS(uva_output, paste(output_dir, "uvaResults.RDS", sep = "/"))
            write.csv(reduced_df, paste(output_dir, "reduced_data.csv", sep = "/"), row.names = FALSE)
        }
    }

    # Return uva_output
    return(uva_output)
}




#' Plot weighted topological overlap (wTO) of Unique Variable Analysis (UVA)
#'
#' Plot a bar chart of the weighted topological overlaps (wTOs) of the
#' perform_uva output
#'
#' @param uva_output List output by perform_uva
#' @param color Color of the barplot (Default: orange)
#' @return A bar chart of the wTOs
#' @export

plot_wTO <- function(uva_output, color = "#ff7f00") {
    wTO_df <- data.frame(uva_output$redundancy$descriptives$centralTendency)
    wTO_df <- cbind(pair = rownames(wTO_df), wTO_df)
    rownames(wTO_df) <- 1:nrow(wTO_df)

    wTO_plot <- ggplot2::ggplot(data=head(wTO_df,20), ggplot2::aes(x=reorder(pair, wTO),y=wTO)) +
        ggplot2::geom_bar(stat="identity", color='black', fill=color, alpha=0.75, width=0.8) +
        ggplot2::ylab('wTO\n') +
        ggplot2::xlab('') +
        ggplot2::ylim(0,0.5) +
        ggplot2::geom_text(ggplot2::aes(label=sprintf("%.3f",wTO)), hjust=-0.3, size=4) +
        ggplot2::theme_minimal() +
        ggplot2::coord_flip() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 12))

    return(wTO_plot)
}