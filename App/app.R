
Sys.setenv(CHROMOTE_CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")

source("ui.R")
source("server.R")

# Create the Shiny app
shinyApp(ui = ui, server = server)


