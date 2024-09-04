# R/app.R
#' Launch the Shiny app
#'
#' This function starts the Shiny app.
#'
#' @export



GTAVTools <- function() {
  shiny::shinyApp(ui = ui, server = server)
}
