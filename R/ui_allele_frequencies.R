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
    .af-export-row {
      display: flex; align-items: center; gap: 6px;
      padding-top: .6rem; border-top: 0.5px solid rgba(0,0,0,.07); margin-top: .6rem;
    }
    .af-export-label { font-size: 11px; color: #9ca3af; }
    .nav-pills > li > a { font-size: 12px; padding: 4px 12px; }
    .box-body .nav-pills { margin-bottom: 10px; }
  "))

  export_row <- function(csv_id, txt_id) {
    tags$div(class = "af-export-row",
      tags$span(class = "af-export-label", "Export:"),
      downloadButton(ns(csv_id), "CSV",
        class = "btn btn-default btn-xs",
        style = "padding:2px 10px; font-size:11px;"),
      downloadButton(ns(txt_id), "TXT",
        class = "btn btn-default btn-xs",
        style = "padding:2px 10px; font-size:11px;")
    )
  }

  fluidPage(
    custom_css,

    module_banner("chart-pie", "Allele Frequencies",
      "Per-population allele tables, frequencies and missing data overview",
      "#2CBF9F"),

    # ── Value boxes ────────────────────────────────────────────────────────
    fluidRow(
      style = "margin-bottom:10px;",
      column(3, tags$div(class = "af-vbox",
        tags$div(class = "af-vbox-icon",
          style = "background:#D6EEF4;color:#3B9AB2;", icon("users")),         # Zissou1 teal
        tags$div(tags$div(class = "af-vbox-label", "Individuals"),
          uiOutput(ns("vb_individuals")))
      )),
      column(3, tags$div(class = "af-vbox",
        tags$div(class = "af-vbox-icon",
          style = "background:#D8EAE4;color:#74A089;", icon("map-marker-alt")), # Royal2 sage
        tags$div(tags$div(class = "af-vbox-label", "Populations"),
          uiOutput(ns("vb_populations")))
      )),
      column(3, tags$div(class = "af-vbox",
        tags$div(class = "af-vbox-icon",
          style = "background:#EAE5EF;color:#9986A5;", icon("dna")),            # IsleofDogs1 mauve
        tags$div(tags$div(class = "af-vbox-label", "Markers"),
          uiOutput(ns("vb_markers")))
      )),
      column(3, tags$div(class = "af-vbox",
        tags$div(class = "af-vbox-icon",
          style = "background:#FBF0D0;color:#E1AF00;", icon("exclamation-triangle")), # Zissou1 amber
        tags$div(tags$div(class = "af-vbox-label", "Missing data"),
          uiOutput(ns("vb_missing")))
      ))
    ),

    # ── Main tabsetPanel ───────────────────────────────────────────────────
    tabsetPanel(
      id   = ns("main_tabs"),
      type = "tabs",

      # ══════════════════════════════════════════════════════════════════ #
      # TAB 1 — Data summary                                              #
      # ══════════════════════════════════════════════════════════════════ #
      tabPanel(
        title = tagList(icon("table"), " Data summary"),
        value = "summary",
        br(),

        fluidRow(
          column(12,
            actionButton(ns("generate_summary"),
              label = tagList(icon("calculator"), tags$strong("Generate data summary")),
              class = "btn-action-primary")
          )
        ),

        br(),

        box(
          title       = tagList(icon("map-marker-alt"), " Population summary"),
          status = "primary",
          solidHeader = TRUE,
          width       = NULL,
          collapsible = TRUE,
          div(style = "overflow-x: auto;",
            tableOutput(ns("populations_Summary_table")))
        ),

        br(),

        box(
          title = tagList(icon("database"),
            " Sample sizes, genotyped and missing data"),
          status = "primary", solidHeader = TRUE, width = NULL,

          tabsetPanel(type = "pills",
            tabPanel("By population \u00d7 marker", br(),
              DT::DTOutput(ns("data_summary_by_pop_locus")),
              export_row("download_summary_pop_locus",
                         "download_summary_pop_locus_txt")),
            tabPanel("By marker \u2014 mean", br(),
              DT::DTOutput(ns("data_summary_by_locus_mean")),
              export_row("download_summary_locus_mean",
                         "download_summary_locus_mean_txt")),
            tabPanel("By marker \u2014 selected pop.", br(),
              fluidRow(column(4,
                selectInput(ns("selected_population_subsamples"), "Population:",
                  choices = c("All populations" = "all"), selected = "all")
              )),
              DT::DTOutput(ns("data_summary_by_Subsamples_locus_sum")),
              export_row("download_summary_Subsamples_locus_sum",
                         "download_summary_Subsamples_locus_sum_txt")),
            tabPanel("By marker \u2014 sum", br(),
              DT::DTOutput(ns("data_summary_by_locus_sum")),
              export_row("download_summary_locus_sum",
                         "download_summary_locus_sum_txt")),
            tabPanel("By population \u2014 sum", br(),
              DT::DTOutput(ns("data_summary_by_pop_sum")),
              export_row("download_summary_pop_sum",
                         "download_summary_pop_sum_txt")),
            tabPanel("Global", br(),
              DT::DTOutput(ns("data_summary_global")),
              export_row("download_summary_global",
                         "download_summary_global_txt"))
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════ #
      # TAB 2 — Allele frequencies (existing multi-subtab view)           #
      # ══════════════════════════════════════════════════════════════════ #
      tabPanel(
        title = tagList(icon("dna"), " Allele frequencies"),
        value = "frequencies",
        br(),

        box(title = tagList(icon("sliders-h"), " Analysis parameters"),
          status = "primary", solidHeader = TRUE, width = NULL,
          fluidRow(
            column(4,
              selectInput(ns("selected_population"), "Population:",
                choices = c("All populations" = "all"), multiple = FALSE)),
            column(4,
              selectizeInput(ns("selected_marker"), "Marker:",
                choices = NULL, multiple = FALSE,
                options = list(placeholder = "All markers"))),
            column(2,
              tags$div(style = "margin-top:25px;",
                actionButton(ns("update_analysis"),
                  label = tagList(icon("play"), tags$strong(" Compute")),
                  class = "btn-action-primary")))
          )
        ),

        box(title = tagList(icon("chart-bar"), " Results"),
          status = "primary", solidHeader = TRUE, width = NULL,

          tabsetPanel(type = "pills",
            tabPanel("By population", br(),
              fluidRow(
                column(9, DT::DTOutput(ns("allele_freq_by_pop"))),
                column(3,
                  tags$div(
                    style = "background:#f9fafb;border-radius:8px;padding:.7rem .85rem;",
                    tags$div(style = "font-weight:500;font-size:12px;margin-bottom:6px;",
                      "Quick statistics"),
                    verbatimTextOutput(ns("quick_stats_pop"), placeholder = TRUE)
                  ),
                  br(),
                  downloadButton(ns("download_freq_pop_csv"), "CSV",
                    class = "btn btn-default btn-sm",
                    style = "width:100%;margin-bottom:4px;"),
                  downloadButton(ns("download_freq_pop_txt"), "TXT",
                    class = "btn btn-default btn-sm",
                    style = "width:100%;")
                )
              )
            ),
            tabPanel("By marker (global)", br(),
              fluidRow(
                column(9, DT::DTOutput(ns("allele_freq_global"))),
                column(3,
                  tags$div(
                    style = "background:#f9fafb;border-radius:8px;padding:.7rem .85rem;",
                    tags$div(style = "font-weight:500;font-size:12px;margin-bottom:6px;",
                      "Global statistics"),
                    verbatimTextOutput(ns("quick_stats_global"), placeholder = TRUE)
                  ),
                  br(),
                  downloadButton(ns("download_freq_global_csv"), "CSV",
                    class = "btn btn-default btn-sm",
                    style = "width:100%;margin-bottom:4px;"),
                  downloadButton(ns("download_freq_global_txt"), "TXT",
                    class = "btn btn-default btn-sm",
                    style = "width:100%;")
                )
              )
            ),
            tabPanel("Plots", br(),
              plotOutput(ns("allele_freq_plot"),     height = "400px"), br(),
              plotOutput(ns("missing_data_plot"),    height = "400px"), br(),
              plotOutput(ns("allele_richness_plot"), height = "400px"), br(),
              downloadButton(ns("download_allele_plots"), ".png",
                class = "btn-download-primary")
            )
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════ #
      # TAB 3 — Fstat-style view  (NEW)                                   #
      # One table per locus: N, missing, all allele freqs (incl. 0),      #
      # Na/Ne/He/Ho/Fis — populations as columns + Global column          #
      # ══════════════════════════════════════════════════════════════════ #
      tabPanel(
        title = tagList(icon("th"), " Frequencies by locus"),
        value = "fstat_view",
        br(),

        box(title = tagList(icon("sliders-h"), " Parameters"),
          status = "primary", solidHeader = TRUE, width = NULL,
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
                  class = "btn-action-primary"))),
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
          status = "primary", solidHeader = TRUE, width = NULL,
          DT::DTOutput(ns("fstat_table"))
        )
      ),

      # ══════════════════════════════════════════════════════════════════ #
      # TAB 4 — Marker details                                            #
      # ══════════════════════════════════════════════════════════════════ #
      tabPanel(
        title = tagList(icon("microscope"), " Marker details"),
        value = "markers",
        br(),
        box(title = tagList(icon("list-alt"),
            " Per-marker summary (allele count, missing, observed alleles)"),
          status = "primary", solidHeader = TRUE, width = NULL,
          DT::DTOutput(ns("marker_details")))
      ),

      # ══════════════════════════════════════════════════════════════════ #
      # TAB 5 — Complete matrix                                           #
      # ══════════════════════════════════════════════════════════════════ #
      tabPanel(
        title = tagList(icon("th-large"), " Complete matrix"),
        value = "matrix",
        br(),

        tags$div(
          style = paste0("background:#D9D0D3;border-left:3px solid #8D8680;",
                         "border-radius:6px;padding:.55rem .9rem;",
                         "font-size:12px;color:#39312F;margin-bottom:1rem;",
                         "display:flex;align-items:flex-start;gap:7px;"),
          icon("info-circle"),
          tags$span(
            "All alleles for each marker across all populations,",
            " including absent alleles (frequency\u00a0=\u00a00).",
            " Same number of rows per marker for every population."
          )
        ),

        box(
          title = tagList(icon("table"),
            " Complete allele frequency matrix"),
          status = "primary", solidHeader = TRUE, width = NULL,
          fluidRow(column(4,
            actionButton(ns("generate_complete_matrix"),
              label = tagList(icon("cogs"), tags$strong(" Generate complete matrix")),
              class = "btn-action-primary")
          )),
          br(),
          DT::DTOutput(ns("complete_freq_matrix")),
          export_row("download_complete_matrix_csv",
                     "download_complete_matrix_txt")
        )
      )
    )
  )
}