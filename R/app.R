# app.R — lets RStudio's "Run App" button (and plain `shiny::runApp()`)
# launch SPG directly from the project folder, without needing to run
# `devtools::load_all()` first or know this is an R-package-structured app.
#
# Runs fully offline, inside RStudio's own Viewer pane (no external browser,
# no internet connection required) — set shiny.launch.browser = FALSE so it
# never tries to open a system/web browser window.

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
} else {
  # Fallback if pkgload/devtools are not installed: source the package's R/
  # files directly (same effect for a Shiny app with no compiled exports
  # beyond what's already built in src/).
  for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
    tryCatch(source(f, local = FALSE), error = function(e) NULL)
  }
}

options(
  shiny.launch.browser = FALSE,  # stay inside RStudio's Viewer pane
  shiny.maxRequestSize = 500 * 1024^2
)

shiny::shinyApp(ui = app_ui, server = app_server)
