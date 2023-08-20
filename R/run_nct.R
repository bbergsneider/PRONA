#' Run Network Comparison Test (NCT)
#'
#' This function runs a network comparison test between two dataframes
#' of symptom severity data
#'
#' @param df1 First dataframe of symptom severity/frequency data
#' @param df2 Second dataframe of symptom severity/frequency data
#' @param normal Boolean. Whether to consider all variables normally distrubted.
#' If false, conducts a non-paranormal transformation (Default: FALSE)
#' @param it Number of bootstrapping iterations to run (Default: 2500)
#' @param p.adjust.methods Character. Can be one of "holm", "hochberg",
#' "hommel", "bonferroni", "BH", "BY", "fdr", or "none". To control
#' (or not) for testing of multiple edges. Defaults to "none".
#' @param test.edges Boolean. Whether to test differences in individual
#' edge weights. (Default: TRUE)
#' @param test.centrality Boolean. Whether to test differences in centrality
#' measures (Default: TRUE)
#' @param centrality Vector of which centrality measures to test (Default:
#' c('closeness','betweenness','strength','expectedInfluence'))
#' @return The output of the NCT function from the NetworkComparisonTest
#' package
#' @export
#' 
run_NCT <- function(df1, df2, normal=FALSE, it = 2500, p.adjust.methods = "none", test.edges=TRUE, test.centrality=TRUE, centrality=c('closeness','betweenness','strength','expectedInfluence')) {
    # Check if data is formatted properly
    check_data_format(df1)
    check_data_format(df2)
    df1 <- dplyr::select(df1, -ID)
    df2 <- dplyr::select(df2, -ID)
    
    if (normal) {
        estimator = 'network_estimation_fun_normal'
    } else {
        estimator = 'network_estimation_fun_non_normal'
    }

    res <- NetworkComparisonTest::NCT(df1, df2, it=it, estimator = estimator, test.edges = test.edges, edges = 'all', p.adjust.methods = p.adjust.methods, test.centrality = test.centrality, centrality = centrality, nodes = 'all')
    return(res)
}


#' Plot Network Comparison Test (NCT)
#'
#' This function plots the results of the run_NCT function
#'
#' @param res The output of the run_NCT function
#' @param what Defines what has to be plotted: results pertaining
#' to test on invariance of global strength ("strength"), network
#' structure ("network"), edge strength ("edge"), or specific
#' centrality measure ("centrality")
#' @return Plot of NCT results
#' @export
#' 
plot_NCT <- function(res, what) {
    return(plot(res, what = what))
}


#' Gets p-values from running NCT
#'
#' This function gets the p-values from running NCT
#'
#' @param res The output of the run_NCT function
#' @param what Defines what measure to get p-values for.
#' Options: "edges" or "centralities"
#' @return Dataframe of p-values for either individual edges
#' or centrality measures
#' @export
#' 
get_NCT_pvalues <- function(res, what) {
    if (what == "edges") {
        return(res$einv.pvals)
    } else if (what == "centralities") {
        return(res$diffcen.pval)
    } else {
        stop("what parameter must equal either 'edges' or 'centralities'")
    }
}