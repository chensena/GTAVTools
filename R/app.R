# R/app.R
#' Launch the Shiny app
#'
#' This function starts the Shiny app.
#'
#' @export
library(DT)


GTAVTools <- function() {
  shiny::shinyApp(ui = ui, server = server)
}
