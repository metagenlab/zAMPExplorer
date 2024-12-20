#' Launch the zAMPExplorer Application
#'
#' zAMPExplorer provides an interactive interface for analyzing and visualizing microbiome sequencing data,
#' particularly the output of the zAMP pipeline (phyloseq objects). It includes functionalities such as:
#' compositional barplots, heatmaps, diversity metrics, beta-diversity analyses, and more.
#'
#' @section Features:
#' \itemize{
#'   \item Read QC, Taxa overviews.
#'   \item Compositional barplots, relative abundance heatmaps.
#'   \item Alpha and beta diversity visualizations.
#'   \item Differential abundance testing using MaAsLin2.
#'   \item Community typing using Dirichlet Multinomial Mixtures.
#'   \item Interactive RDA ordination plots.
#' }
#'
#' @seealso Useful links:
#' \itemize{
#'   \item Documentation: \url{https://github.com/metagenlab/zAMPExplorer}
#'   \item zAMP pipeline: \url{https://zamp.readthedocs.io/en/latest/}
#' }
#' @return A `shiny.appobj` object representing the running Shiny application.
#'
#' @param ... Arguments to pass to golem options. See \code{\link[golem]{get_golem_options}} for more details.
#' @inheritParams shiny::shinyApp
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
#' @importFrom shinydashboard dashboardBody dashboardHeader dashboardPage dashboardSidebar menuItem tabItems tabItem
#' @importFrom fs path file_exists
#' @examples
#' if (interactive()) {
#'   zAMPExplorer_app()
#' }
#' @export

zAMPExplorer_app <- function(
    onStart = NULL,
    options = list(),
    enableBookmarking = NULL,
    uiPattern = "/",
    ...) {
  # Check if microViz is installed
  if (!requireNamespace("microViz", quietly = TRUE)) {
    stop(
      "The 'microViz' package is required for zAMPExplorer. Please install it first using one of the following methods:\n",
      "- From R-universe: \n",
      "  install.packages('microViz', repos = c(davidbarnett = 'https://david-barnett.r-universe.dev', getOption('repos')))\n",
      "- Or from GitHub: \n",
      "  remotes::install_github('david-barnett/microViz')\n"
    )
  }



  # Launch the Shiny app
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
