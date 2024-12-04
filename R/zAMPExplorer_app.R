#' Launch the zAMPExplorer Shiny App
#'
#' This function launches the Shiny app for zAMPExplorer. It assumes that the app files
#' (ui.R and server.R) are located in the `inst/shiny` directory.
#' 
#' @export
#'
#' @import shiny
#' @importFrom shinydashboard dashboardBody dashboardHeader dashboardPage dashboardSidebar menuItem tabItems tabItem
#' @importFrom DT datatable DTOutput
#' @importFrom plotly plotlyOutput renderPlotly
#' @importFrom InteractiveComplexHeatmap InteractiveComplexHeatmapOutput
#' @importFrom MicrobiotaProcess ggrarecurve
#' @importFrom RColorBrewer brewer.pal
#' @importFrom phyloseq phyloseq distance otu_table sample_data rank_names sample_variables
#' @importFrom phyloseq tax_table prune_samples sample_sums sample_names
#' @importFrom dplyr filter mutate select
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom vegan capscale adonis vegdist
#' @importFrom writexl write_xlsx
#' @importFrom fs path file_exists
#' @importFrom microViz comp_barplot ord_plot tax_agg tax_sort dist_calc
#' @importFrom microbiome transform aggregate_taxa core core_members
#' @importFrom shinyFiles shinyFilesButton shinyDirButton
#' @importFrom shinyWidgets pickerInput
zAMPExplorer <- function() {
  # Check if the app directory exists
  app_dir <- system.file("shiny", package = "zAMPExplorer")
  if (app_dir == "") {
    stop("Shiny app directory not found. Please ensure the package is correctly installed.")
  }
  # Launch the Shiny app
  shiny::runApp(app_dir)
}
