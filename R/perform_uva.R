#' Calculate wTO
#'
#' This function calculates the wTO of symptoms in a dataframe based on the
#' Unique Variable Analysis (UVA) method described in Christensen, Garrido,
#' & Colino (Multivariate Behavioral Research, 2023)
#'
#' @param df Dataframe representing symptom severity/frequency
#' @return The list of the outputs from the UVA function as implemented
#' in the R package EGAnet
#' @export
#' 
calculate_wTO <- function(df) {
    check_data_format(df) # Check if df is formatted properly
    df <- dplyr::select(df, -ID) # Get rid of ID column

    # Conduct UVA
    uva_output <- EGAnet::UVA(data = df, reduce = FALSE, reduce.method = "latent")

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
#' 
plot_wTO <- function(uva_output, color = "#ff7f00") {
    wTO_df <- uva_output$wto$pairwise
    wTO_df$pair <- paste(wTO_df$node_i, wTO_df$node_j, sep = "-")

    wTO_plot <- ggplot2::ggplot(data=head(wTO_df,20), ggplot2::aes(x=reorder(pair, wto),y=wto)) +
        ggplot2::geom_bar(stat="identity", color='black', fill=color, alpha=0.75, width=0.8) +
        ggplot2::ylab('wTO\n') +
        ggplot2::xlab('') +
        ggplot2::ylim(0,0.5) +
        ggplot2::geom_text(ggplot2::aes(label=sprintf("%.3f",wto)), hjust=-0.3, size=4) +
        ggplot2::theme_minimal() +
        ggplot2::coord_flip() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 12))

    return(wTO_plot)
}


#' Perform Unique Variable Analysis (UVA) consolidation
#'
#' This function performs Unique Variable Analysis (UVA) consolidation
#' on a symptom severity dataframe, given a certain wTO cutoff.
#'
#' For more information on UVA, see Christensen, Garrido, & Colino
#' (Multivariate Behavioral Research, 2023)
#'
#' @param df Dataframe representing symptom severity/frequency
#' @param scale The scale on which symptom severity/frequency is measured.
#' For example if symptoms are measured on a 0-10 scale, scale = 10.
#' (default: 1)
#' @param cut.off Cut-off used to determine when pairwise wto values are
#' considered locally dependent (or redundant). Must be values between 0 and 1.
#' Defaults to 0.25.
#' @param reduce.method Method to reduce redundancies. Available options:
#' "latent", "mean", "remove", "sum". See Christensen, Garrido, & Colino
#' for more details (Default: latent)
#' @param new.names Vector of new names to give to consolidated variables.
#' Variable pairs will be renamed in descending order of wTO. If this vector
#' is not given, new variable pairs will be renamed "CV1", "CV2", etc...
#' Moreover, if reduce.method = "remove", this vector will not be used.
#' (Default: NULL)
#' @param output_dir Directory in which the .RDS and .csv outputs of UVA should
#' be saved (defaults to current directory)
#' @return The list of the outputs from the UVA function as implemented
#' in the R package EGAnet
#' @export

perform_uva_consolidation <- function(df, scale = 1, cut.off = 0.25, reduce.method = "latent", new.names = NULL, output_dir = "") {
    check_data_format(df) # Check if df is formatted properly
    ids <- df["ID"] # Save ID column
    df <- dplyr::select(df, -ID) # Get rid of ID column

    # Conduct UVA
    uva_output <- EGAnet::UVA(data = df, reduce.method = reduce.method, reduce = TRUE, cut.off = cut.off)

    # Save reduced data in a new dataframe
    reduced_df <- uva_output[['reduced_data']]
    # Re-add ID column to reduced data
    reduced_df <- cbind(ids, reduced_df)

    # If latent varaible reduction has been performed:
    # Latent variable reduction yields variable values that are not
    # on the same scale as the original data. In order to keep latent
    # variable values consistent with other variables, latent variable
    # values must be normalized to be on the same scale
    if (reduce.method == "latent") {
        for (lv in c(1:length(uva_output$redundant))) {
            lv <- paste("CV", lv, sep = "")
            reduced_df[,lv] <- (reduced_df[,lv] - min(reduced_df[,lv])) / (max(reduced_df[,lv]) - min(reduced_df[,lv])) * scale
        }
    }

    # Add new variable names if specified
    if (!is.null(new.names) && reduce.method != "remove") {
        if (length(new.names) != length(uva_output$redundant)) {
            stop("Length of new.names does not match how many variable pairs were consolidated")
        }
        for (index in c(1:length(uva_output$redundant))) {
            lv <- paste("CV", index, sep = "")
            colnames(reduced_df)[colnames(reduced_df) == lv] <- new.names[index]
        }
    }

    # Save uva_output and reduced_df in output files
    if (output_dir == "") {
        saveRDS(uva_output, "uvaResults.RDS")
        write.csv(reduced_df, "reduced_data.csv", row.names = FALSE)
    } else {
        saveRDS(uva_output, paste(output_dir, "uvaResults.RDS", sep = "/"))
        write.csv(reduced_df, paste(output_dir, "reduced_data.csv", sep = "/"), row.names = FALSE)
    }

    # Return uva_output
    return(uva_output)
}