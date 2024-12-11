#' Launch the zAMPExplorer Application
#'
#' Launches the Shiny application by calling the `app_ui` and `app_server` functions.
#'
#' @param ... Arguments to pass to golem options. See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
#' @importFrom ggplot2 ggplot aes geom_bar
#' @importFrom utils combn write.csv write.table
#' @importFrom shinydashboard dashboardBody dashboardHeader dashboardPage dashboardSidebar menuItem tabItems tabItem
#' @importFrom DT datatable DTOutput
#' @importFrom plotly plotlyOutput renderPlotly
#' @importFrom InteractiveComplexHeatmap InteractiveComplexHeatmapOutput
#' @importFrom MicrobiotaProcess ggrarecurve
#' @importFrom RColorBrewer brewer.pal
#' @importFrom phyloseq phyloseq distance otu_table sample_data rank_names sample_variables
#' @importFrom phyloseq tax_table prune_samples sample_sums sample_names sample_data
#' @importFrom dplyr filter mutate select
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom vegan capscale adonis vegdist
#' @importFrom writexl write_xlsx
#' @importFrom fs path file_exists
#' @importFrom microViz comp_barplot ord_plot tax_agg tax_sort dist_calc
#' @importFrom microbiome transform aggregate_taxa core core_members
#' @importFrom shinyFiles shinyFilesButton shinyDirButton
#' @importFrom shinyWidgets pickerInput
#'
if (!requireNamespace("microViz", quietly = TRUE)) {
  stop("The 'microViz' package is not installed. Please install it before running zAMPExplorer:\n",
       "install.packages('microViz') or remotes::install_github('david-barnett/microViz')")
}

zAMPExplorer_app <- function(
    onStart = NULL,
    options = list(),
    enableBookmarking = NULL,
    uiPattern = "/",
    ...
) {
  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}

