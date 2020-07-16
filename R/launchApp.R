#' launches the shinyAppDemo app
#'
#' @export launchApp
#'
#' @return shiny application object
#'
#'
#' @import shiny
#'


# wrapper for shiny::shinyApp()
launchApp <- function() {
  library(shiny)
  shinyApp(ui = shinyAppUI, server = shinyAppServer)
}
