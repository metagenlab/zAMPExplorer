#' @title Launch the Shiny App
#' @description This function launches the Shiny app.
#' @export
zAMPExplorer_app <- function() {
  appDir <- system.file("shiny", package = "zAMPExplorer")
  if (appDir == "") {
    stop("Could not find the Shiny app directory. Try re-installing the `zAMPExplorer` package.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
