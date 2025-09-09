#' Launch the gsEstimators Shiny App
#'
#' This function runs the Shiny application bundled in the package.
#'
#' @export
run_app <- function(...) {
  shiny::shinyApp(ui = ui, server = server)
}
