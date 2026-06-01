#' Application Server
#'
#' Handles both the splash screen and the main application.
#' When the user clicks "Launch" on the splash, the UI is swapped
#' to the full \code{app_ui()} via \code{shiny::renderUI}.
#'
#' @param input,output,session Shiny server arguments.
#' @importFrom shiny reactiveValues observe renderUI
#' @export
app_server <- function(input, output, session) {

  # ── Shared reactive state ────────────────────────────────────────────────
  rv <- shiny::reactiveValues(
    raw            = NULL,
    data           = NULL,
    formatted_data = NULL,
    launched       = FALSE      # FALSE = splash, TRUE = main app
  )

  # ── Swap to main app when Launch is clicked ──────────────────────────────
  shiny::observeEvent(input$btn_launch, {
    rv$launched <- TRUE
  })

  # ── Render the correct UI layer ──────────────────────────────────────────
  output$main_ui <- shiny::renderUI({
    if (!rv$launched) {
      splash_ui()
    } else {
      app_ui()
    }
  })

  # ── Wire up analysis modules (only active after launch) ──────────────────
  # We use observe + req so that module servers only run once the data UI
  # is actually present in the DOM.
  shiny::observe({
    shiny::req(rv$launched)

    server_import_data("import",          rv)
    server_allele_frequencies("allele",   rv)
    server_general_stats("general_stats", rv)
    server_LD("ld",                       rv)
    server_null_alleles("null_alleles",   rv)
    server_isolation_by_distance("ibd",   rv)
  })
}


#' Run the ShinyPopGen application
#'
#' Starts with the splash screen; clicking "Launch" transitions to the
#' full analysis UI without a full page reload.
#'
#' @param ... Arguments passed to \code{\link[shiny]{shinyApp}}.
#' @return A \code{shiny.appobj} (invisibly).
#' @examples
#' \dontrun{
#'   run_app()
#' }
#' @export
run_app <- function(...) {
  options(shiny.maxRequestSize = 500 * 1024^2)   # 500 MB upload limit

  # Minimal bootstrap wrapper: a single uiOutput that the server swaps
  # between splash_ui() and app_ui().
  root_ui <- function() {
    shiny::fluidPage(
      # No chrome of its own — the child UIs bring their own <head> content.
      shiny::uiOutput("main_ui")
    )
  }

  shiny::shinyApp(
    ui      = root_ui,
    server  = app_server,
    options = list(...)
  )
}
