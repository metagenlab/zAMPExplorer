#' Launch the zAMPExplorer Shiny App
#'
#' This function launches the Shiny app for zAMPExplorer. It assumes that the app files
#' (ui.R and server.R) are located in the `inst/shiny` directory.
#' @export
zAMPExplorer <- function() {
  # Check if the app directory exists
  app_dir <- system.file("shiny", package = "zAMPExplorer")
  if (app_dir == "") {
    stop("Shiny app directory not found. Please ensure the package is correctly installed.")
  }
  # Launch the Shiny app
  shiny::runApp(app_dir)
}
