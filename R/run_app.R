#' Run the ShinyPopGen application
#'
#' Launches the app inside RStudio's own Viewer pane when running from
#' RStudio (self-contained, works fully offline, no separate browser
#' window) — falls back to the system's default browser when RStudio isn't
#' available (plain R console, Rscript).
#'
#' @param ... Arguments passed to \code{\link[shiny]{runApp}} (e.g. `port`,
#'   `host`). Pass `launch.browser = TRUE` yourself to force the system
#'   browser even inside RStudio.
#' @return Invisibly, the result of \code{\link[shiny]{runApp}}.
#' @examples
#' \dontrun{
#'   run_app()
#' }
#' @export
run_app <- function(...) {
  options(shiny.maxRequestSize = 500 * 1024^2)  # 500 MB upload limit

  dots <- list(...)
  if (is.null(dots$launch.browser)) {
    in_rstudio <- requireNamespace("rstudioapi", quietly = TRUE) &&
      tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE)
    dots$launch.browser <- if (in_rstudio) rstudioapi::viewer else TRUE
  }

  app <- shiny::shinyApp(ui = app_ui, server = app_server)
  do.call(shiny::runApp, c(list(appDir = app), dots))
}
