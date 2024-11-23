#' @title Launch the Shiny App
#' @description This function launches the Shiny app and installs any missing dependencies.
#' @export
zAMPExplorer <- function() {
  library(shiny)
  shinyApp(ui = ui, server = server)
}
