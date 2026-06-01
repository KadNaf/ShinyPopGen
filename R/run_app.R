#' Application Server
#'
#' @param input,output,session Shiny server arguments.
#' @importFrom shiny reactiveValues observeEvent observe req renderUI
#' @export
app_server <- function(input, output, session) {

  rv <- shiny::reactiveValues(
    raw            = NULL,
    data           = NULL,
    formatted_data = NULL,
    launched       = FALSE
  )

  # Swap to main app when Launch is clicked
  shiny::observeEvent(input$btn_launch, {
    rv$launched <- TRUE
  })

  # Render the correct UI layer
  output$main_ui <- shiny::renderUI({
    if (!rv$launched) {
      splash_ui()
    } else {
      app_ui()
    }
  })

  # Wire analysis modules — only after launch
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
#' @param ... Arguments passed to \code{\link[shiny]{shinyApp}}.
#' @return A \code{shiny.appobj} (invisibly).
#' @examples
#' \dontrun{
#'   run_app()
#' }
#' @export
run_app <- function(...) {
  options(shiny.maxRequestSize = 500 * 1024^2)

  # Minimal root page: a single uiOutput that the server swaps.
  # No CSS, no theme — completely neutral so splash and app_ui
  # each own their own styles 100%.
  root_ui <- function() {
    shiny::fluidPage(
      shiny::uiOutput("main_ui")
    )
  }

  shiny::shinyApp(
    ui      = root_ui,
    server  = app_server,
    options = list(...)
  )
}