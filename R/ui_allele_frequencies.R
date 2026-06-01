ui_allele_frequencies <- function(id) {
  ns <- NS(id)

  custom_css <- tags$style(HTML("
    .af-vbox {
      background: white; border: 0.5px solid rgba(0,0,0,.1);
      border-radius: 10px; padding: .75rem 1rem;
      display: flex; align-items: center; gap: 10px; margin-bottom: 4px;
    }
    .af-vbox-icon {
      width: 36px; height: 36px; border-radius: 8px;
      display: flex; align-items: center; justify-content: center;
      font-size: 15px; flex-shrink: 0;
    }
    .af-vbox-label { font-size: 11px; color: #6b7280; margin-bottom: 1px; }
    .af-vbox-val   { font-size: 22px; font-weight: 500; line-height: 1.1; }
    .nav-pills > li > a { font-size: 12px; padding: 4px 12px; }
    .box-body .nav-pills { margin-bottom: 10px; }
  "))

  fluidPage(
    custom_css,

    # ── Value boxes ────────────────────────────────────────────────────────
    fluidRow(
      style = "margin-bottom:10px;",
      column(3, tags$div(class = "af-vbox",
        tags$div(class = "af-vbox-icon",
          style = "background:#E6F1FB;color:#185FA5;", icon("users")),
        tags$div(tags$div(class = "af-vbox-label", "Individuals"),
          uiOutput(ns("vb_individuals")))
      )),
      column(3, tags$div(class = "af-vbox",
        tags$div(class = "af-vbox-icon",
          style = "background:#EAF3DE;color:#3B6D11;", icon("map-marker-alt")),
        tags$div(tags$div(class = "af-vbox-label", "Populations"),
          uiOutput(ns("vb_populations")))
      )),
      column(3, tags$div(class = "af-vbox",
        tags$div(class = "af-vbox-icon",
          style = "background:#EEEDFE;color:#534AB7;", icon("dna")),
        tags$div(tags$div(class = "af-vbox-label", "Markers"),
          uiOutput(ns("vb_markers")))
      )),
      column(3, tags$div(class = "af-vbox",
        tags$div(class = "af-vbox-icon",
          style = "background:#FAEEDA;color:#854F0B;", icon("exclamation-triangle")),
        tags$div(tags$div(class = "af-vbox-label", "Missing data"),
          uiOutput(ns("vb_missing")))
      ))
    ),

    # ── Fstat-style view ───────────────────────────────────────────────────
    br(),

    box(title = tagList(icon("sliders-h"), " Parameters"),
      status = "primary", solidHeader = FALSE, width = NULL,
      fluidRow(
        column(4,
          selectInput(ns("fstat_population"), "Population:",
            choices = c("All populations" = "all"), multiple = FALSE)),
        column(4,
          selectizeInput(ns("fstat_marker"), "Marker:",
            choices = NULL, multiple = FALSE,
            options = list(placeholder = "All markers"))),
        column(2,
          tags$div(style = "margin-top:25px;",
            actionButton(ns("update_fstat"),
              label = tagList(icon("play"), tags$strong(" Compute")),
              class = "btn-primary"))),
        column(2,
          tags$div(style = "margin-top:25px; display:flex; gap:4px;",
            downloadButton(ns("download_fstat_csv"), "CSV",
              class = "btn btn-default btn-sm"),
            downloadButton(ns("download_fstat_txt"), "TXT",
              class = "btn btn-default btn-sm")))
      )
    ),

    box(
      title = tagList(
        icon("th"),
        " Allele frequencies by locus \u2014 populations in columns",
        tags$span(
          style = paste0("font-size:11px;color:#6b7280;",
                         "margin-left:8px;font-weight:400;"),
          "(N genotyped \u00b7 N missing \u00b7 allele frequencies incl. 0 \u00b7",
          " Na / Ne / He / Ho / Fis)"
        )
      ),
      status = "info", solidHeader = FALSE, width = NULL,
      DT::DTOutput(ns("fstat_table"))
    )
  )
}