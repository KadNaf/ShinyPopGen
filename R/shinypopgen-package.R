#' @keywords internal
#'
#' @import shiny
#' @import bslib
#' @import shinydashboard
#' @import ggplot2
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
#' @importFrom waiter Waiter useWaiter
#' @importFrom htmltools tags HTML tagList withTags
#' @importFrom leaflet leafletOutput renderLeaflet leafletProxy
#' @importFrom DT DTOutput renderDT
#' @importFrom plotly plotlyOutput renderPlotly ggplotly
#' @importFrom utils combn
#' @importFrom stats aggregate median rbinom sd
#' @importFrom utils head write.csv write.table
#' @importFrom graphics plot.new text
#' @importFrom waiter spin_3 transparent
#' @useDynLib shinypopgen, .registration = TRUE
"_PACKAGE"

# Suppress R CMD check NOTEs for symbols that cannot be resolved statically.
utils::globalVariables(c(
  # C++ functions (RcppExports)
  "boot_indiv_wc84_fst",
  "boot_indiv_wc84_fst_parallel",
  "observed_wc84_fst",
  "batch_permute_wc84_stats",
  "batch_permute_fit_global",
  "batch_permute_wc_fis",
  "boot_indiv_wc_fis",
  "boot_popblock_wc_fis",
  "fis_wc_cpp",
  "calculate_observed_fis",
  "compute_ld_pvalues",
  "wc84_components_fst",
  # base R functions R CMD check cannot trace
  "setNames",
  # legacy parallel helpers (parallel / doParallel)
  "makeCluster", "stopCluster", "clusterExport",
  "registerDoParallel", "foreach", "%dopar%",
  # data-masking column names used in dplyr / ggplot2 aes()
  "pop", "Allele", "Frequency", "Population", "Marker",
  "Missing_Proportion", "Na",
  # ggplot2 / waiter column references in server_general_stats
  "Hs", "ID", "Observed_FIS", "Observed_FIT", "P_value",
  "Significant", "CI_L", "CI_U",
  "Observed_FST", "Observed_HT", "Observed_HS",
  # dplyr helpers
  "filter_all", "median",
  # htmltools bare tag functions inside withTags() context
  "thead", "tr", "th", "table",
  # server_null_alleles column names
  "Null_Freq", "Ho", "He", "Locus",
  "Null_Freq_Display", "Ho_Display", "He_Display",
  "N_Samples", "N_Missing", "Impact_Level", "Coding_Method"
))
