# ui_import_data.R
import_data_ui <- function(id) {
  ns <- NS(id)

  box_title_style <- "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;"

  fluidPage(

    # CSS: equal-height top row boxes
    tags$style(HTML("
      .import-top-row > .row {
        display: flex !important;
        flex-wrap: wrap;
        align-items: stretch;
      }
      .import-top-row > .row > [class*='col-'] {
        display: flex !important;
        flex-direction: column;
      }
      .import-top-row > .row > [class*='col-'] > .box {
        flex: 1 1 auto !important;
      }
    ")),

    module_banner("upload", "Import Data",
      "Upload CSV/TXT \u00b7 Map columns \u00b7 Assign populations and markers",
      "#6B64EF"),

    # ── Row 1: two config boxes side by side, same height ───────────────────
    tags$div(
      class = "import-top-row",
      fluidRow(

      column(
        width = 6,
        box(
          width = 12,
          title = div(style = box_title_style, icon("magic"), " Upload your data, automatic formatting"),
          solidHeader = TRUE,

          fileInput(ns("file1"), "Choose CSV File", multiple = FALSE,
                    accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
          checkboxInput(ns("header"), "Header", TRUE),
          radioButtons(ns("sep"), "Separator",
                       choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                       selected = "\t"),
          br(),
          actionButton(ns("load_user_data"),    "Load Data",         icon = icon("rocket"),   class = "btn-action-primary"),
          actionButton(ns("load_default_data"), "Load Default Data", icon = icon("database"), class = "btn-action-secondary")
        )
      ),

      column(
        width = 6,
        box(
          width = 12,
          title = div(style = box_title_style, icon("sliders-h"), " Manual column assignment and formatting"),
          solidHeader = TRUE,
          footer = "* mandatory fields",

          selectizeInput(ns("pop_data"),       "Population name*",    choices = NULL, options = list(placeholder = "select")),
          selectizeInput(ns("latitude_data"),  "Latitude",            choices = NULL, options = list(placeholder = "select")),
          selectizeInput(ns("longitude_data"), "Longitude",           choices = NULL, options = list(placeholder = "select")),
          textInput(ns("metadata_ranges"),
                    "Metadata columns (indices) \u2014 format: 1-4 or 5:10 or 1-3,5",
                    value = ""),
          textInput(ns("col_ranges_data"),
                    "Marker locus columns* \u2014 format: 1-4"),
          textInput(ns("missing_code"), "Code for missing data", value = 0),
          br(),
          actionButton(ns("run_assign"), "Assign metadata", icon = icon("rocket"), class = "btn-action-primary")
        )
      )
    )),  # closes fluidRow + tags$div(.import-top-row)

    # ── Row 2: Imported data preview ─────────────────────────────────────────
    fluidRow(
      box(
        width = 12,
        title = div(style = box_title_style, icon("table"), " Imported data preview"),
        solidHeader = TRUE,

        div(style = "overflow-x: auto; width: 100%;",
            DT::DTOutput(ns("formatted_table"))),
        br(),
        downloadButton(ns("download_csv_transformed"), ".csv", class = "btn-download-primary"),
        downloadButton(ns("download_txt_transformed"), ".txt", class = "btn-download-secondary")
      )
    ),

    # ── Row 3: Map (user-resizable, taller default) ───────────────────────────
    fluidRow(
      box(
        width = 12,
        title = div(style = box_title_style, icon("map-marked-alt"), " Map"),
        solidHeader = TRUE,

        # resize: vertical needs overflow: auto to show the drag handle
        div(
          id    = ns("map_resize_container"),
          class = "spg-map-resize",
          style = paste0(
            "resize: vertical; overflow: auto; ",
            "height: 650px; min-height: 300px; max-height: 90vh; ",
            "border: 1px solid #D9D0D3; border-radius: 4px;"
          ),
          leafletOutput(ns("map"), height = "100%")
        ),
        # JS: ResizeObserver invalidates the leaflet map when container is dragged
        tags$script(HTML(sprintf("
          (function() {
            var containerId = '%s';
            function invalidateMap(container) {
              var widgets = container.querySelectorAll('.leaflet.html-widget');
              widgets.forEach(function(el) {
                if (window.HTMLWidgets) {
                  var b = HTMLWidgets.find('#' + el.id);
                  if (b && typeof b.getMap === 'function') {
                    var m = b.getMap();
                    if (m) { setTimeout(function() { m.invalidateSize({animate:false}); }, 20); }
                  }
                }
              });
            }
            function init() {
              var container = document.getElementById(containerId);
              if (!container) { setTimeout(init, 300); return; }
              if (!window.ResizeObserver) return;
              var ro = new ResizeObserver(function(entries) {
                entries.forEach(function(e) { invalidateMap(e.target); });
              });
              ro.observe(container);
            }
            if (document.readyState === 'loading') {
              document.addEventListener('DOMContentLoaded', init);
            } else { setTimeout(init, 300); }
          })();
        ", ns("map_resize_container")))),
        br(),
        downloadButton(ns("download_map"), ".png", class = "btn-download-primary")
      )
    ),

    # ── Row 4: Imported object summary (adaptive height) ─────────────────────
    fluidRow(
      box(
        width = 12,
        title = div(style = box_title_style, icon("info-circle"), " Imported object summary"),
        solidHeader = TRUE,

        div(
          style = "overflow: auto; max-height: 60vh; background: #F8F8F8; padding: 10px; border-radius: 4px;",
          verbatimTextOutput(ns("formatted_summary"))
        )
      )
    )
  )
}
