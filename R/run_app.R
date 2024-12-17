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
#' @importFrom shinydashboard dashboardBody dashboardHeader dashboardPage dashboardSidebar menuItem tabItems tabItem
#' @importFrom fs path file_exists



zAMPExplorer_app <- function(
    onStart = NULL,
    options = list(),
    enableBookmarking = NULL,
    uiPattern = "/",
    ...
) {
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
