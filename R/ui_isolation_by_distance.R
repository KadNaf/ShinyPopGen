# ui_isolation_by_distance.R
# Isolation by Distance (Rousset 1997) + generic Mantel test.
#
# This module is a CONTINUATION of the "Null alleles" module — it reuses the
# pairwise FST/FST-ENA/DCSE/DCSE-INA (+ bootstrap CI) already computed there
# (rv$null_alleles_results) instead of recomputing them, so the tabs that
# used to duplicate that module (Null allele frequencies, Fst/Fst-ENA,
# DCSE/DCSE-INA, Per-locus x pair) have been removed from here.
#
# Workflow: go to the "Null alleles" module first, choose your per-locus
# coding, click "Compute + Bootstrap + Export" — THEN come here.

isolation_by_distance_UI <- function(id) {
  ns <- NS(id)

  fluidPage(
    tags$head(gs_head()),

    module_banner(
      "atom",
      "Isolation by Distance \u00b7 Mantel Test",
      "Rousset (1997) IBD regression \u00b7 Mantel (1967) permutation test \u2014 built on the pairwise FST-ENA / DCSE-INA already computed in the Null Alleles module",
      "#0c4a6e"
    ),

    tags$div(
      class = "spg-method-note", style = "border-left-color:#0c4a6e;",
      HTML(paste0(
        "This module reuses the pairwise F<sub>ST</sub> / F<sub>ST</sub>-ENA / D<sub>CSE</sub> / ",
        "D<sub>CSE</sub>-INA (+ bootstrap CI) already computed in the <b>Null Alleles</b> module — ",
        "nothing is recomputed here. Go there first, choose your per-locus coding, and click ",
        "<b>\"Compute + Bootstrap + Export\"</b>, then come back to this module.",
        "<br>Geographic distance (D<sub>geo</sub>) is the <b>Vincenty ellipsoidal geodesic distance</b> ",
        "(WGS84), in metres, from each population's GPS centroid (mean Latitude/Longitude of its individuals)."
      ))
    ),

    fluidRow(
      box(
        width = 12, solidHeader = TRUE, status = "primary",
        title = div(style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
                    icon("chart-bar"), " Summary (from Null Alleles module)"),
        fluidRow(
          column(3, valueBoxOutput(ns("box_nloci"),  width = NULL)),
          column(3, valueBoxOutput(ns("box_npops"),  width = NULL)),
          column(3, valueBoxOutput(ns("box_fstena"), width = NULL)),
          column(3, valueBoxOutput(ns("box_nboot"),  width = NULL))
        ),
        uiOutput(ns("ui_run_status"))
      )
    ),

    tabsetPanel(
      id = ns("ibd_tabs"), type = "tabs",

      # ══════════════════════════════════════════════════════════════════
      # TAB 1 — Isolation by Distance (now the FIRST tab)
      # ══════════════════════════════════════════════════════════════════
      tabPanel(title = tagList(icon("chart-line"), " Isolation by Distance"), value = "tab_ibd",
        br(),
        fluidRow(
          box(width = 3, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("sliders-h"), " Parameters"),
            radioButtons(ns("ibd_model"), "Habitat model:",
              choices = c("2D (F_R ~ ln(D_geo))" = "2D",
                          "1D (F_R ~ D_geo)"      = "1D"),
              selected = "2D"),
            radioButtons(ns("ibd_metric"), "Genetic distance metric:",
              choices = c("F_R (raw FST)"      = "raw",
                          "F_R (FST-ENA)"      = "ena"),
              selected = "ena"),
            tags$hr(),
            tags$p(style="color:#777;font-size:11px;",
              "Requires GPS (Latitude/Longitude) set at import for at least 2 populations. ",
              "Population centroid is the mean GPS of its individuals; distance is the ",
              "Vincenty ellipsoidal geodesic distance (WGS84), in metres."),
            actionButton(ns("run_ibd"), "Run IBD Regression",
                         icon = icon("rocket"), class = "btn-action-primary btn-block",
                         style = "font-weight:bold;")
          ),
          box(width = 9, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("chart-line"), " Regression summary (slope / b / Nb / Nem)"),
            DT::DTOutput(ns("dt_ibd_reg")),
            tags$br(),
            uiOutput(ns("ui_ibd_interpretation"))
          )
        ),

        fluidRow(
          box(width = 12, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("table"), " Full pairwise table (Farm1, Farm2, D_geo, FST-FreeNA, F_R, D_CSE-INA, D_CSE)"),
            DT::DTOutput(ns("dt_ibd_table")),
            tags$br(),
            downloadButton(ns("dl_ibd_csv"), "Download full table (CSV)", class = "btn-action-secondary btn-sm")
          )
        ),

        fluidRow(
          box(width = 12, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("chart-area"), " IBD plot"),
            plotly::plotlyOutput(ns("ibd_plot"), height = "440px")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # TAB 2 — Mantel test
      # ══════════════════════════════════════════════════════════════════
      tabPanel(title = tagList(icon("project-diagram"), " Mantel Test"), value = "tab_mantel",
        br(),

        tags$div(
          class = "spg-method-note", style = "border-left-color:#7A5DC7;",
          HTML(paste0(
            "Generic Mantel permutation test between any two pairwise distances (genetic, ",
            "geographic, temporal, ecological or categorical), using a table of pairs in rows / ",
            "distances in columns (RT, Manly 2018; Fstat 2.9.4 convention). Supports ",
            "<b>rectangular matrices</b>: pairs can be excluded (e.g. to keep contemporaneous ",
            "pairs only) without dropping every pair involving the corresponding sub-samples. ",
            "Permutation is by <b>joint row/column relabelling</b> of one matrix, which stays ",
            "valid when either matrix is incomplete. Statistic: Pearson's r (Fstat convention) ",
            "or the slope of the Rousset (1997) regression (Genepop convention for IBD). ",
            "One-sided p-value = (b+1)/(m+1), b = number of permuted statistics \u2265 observed."
          ))
        ),

        fluidRow(
          box(width = 4, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("database"), " Data source"),
            radioButtons(ns("mt_source"), NULL,
              choices = c(
                "Internal pairwise table (from Null Alleles module)" = "internal",
                "Upload external column file"                        = "upload"
              ), selected = "internal"),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'upload'", ns("mt_source")),
              fileInput(ns("mt_file"), "File (Pop1, Pop2, dist1, dist2, ...):",
                        accept = c(".csv", ".txt", ".tsv")),
              radioButtons(ns("mt_sep"), "Separator:",
                choices = c("Comma"=",", "Tab"="\t", "Semicolon"=";"),
                selected = ",", inline = TRUE),
              checkboxInput(ns("mt_header"), "File has header row", value = TRUE)
            ),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'internal'", ns("mt_source")),
              checkboxInput(ns("mt_use_extra"),
                "Merge an extra distance file (temporal / ecological / categorical)",
                value = FALSE),
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns("mt_use_extra")),
                fileInput(ns("mt_extra_file"), "Extra file (first 2 cols = Pop1, Pop2 IDs):",
                          accept = c(".csv", ".txt", ".tsv")),
                radioButtons(ns("mt_extra_sep"), "Separator:",
                  choices = c("Comma"=",", "Tab"="\t", "Semicolon"=";"),
                  selected = ",", inline = TRUE),
                checkboxInput(ns("mt_extra_header"), "File has header row", value = TRUE)
              )
            ),

            tags$hr(),
            tags$div(style="font-size:12px;color:#555;margin-bottom:6px;", "Column assignment:"),
            uiOutput(ns("mt_col_pop1_ui")),
            uiOutput(ns("mt_col_pop2_ui")),
            uiOutput(ns("mt_col_x_ui")),
            uiOutput(ns("mt_col_y_ui"))
          ),

          box(width = 8, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("sliders-h"), " Mantel parameters & results"),
            fluidRow(
              column(4,
                radioButtons(ns("mt_stat"), "Statistic:",
                  choices = c("Pearson r" = "r", "Regression slope (Rousset)" = "b"),
                  selected = "r"),
                checkboxInput(ns("mt_log_x"), "ln(transform) X", value = FALSE)
              ),
              column(4,
                numericInput(ns("mt_n_perm"), "Permutations:",
                             value = 10000, min = 99, max = 200000, step = 1000),
                tags$p(style="color:#777;font-size:11px;", "Advised \u2265 1000.")
              ),
              column(4,
                textInput(ns("mt_exclude"), "Exclude pairs ('ID1-ID2', comma-sep):", value = ""),
                actionButton(ns("run_mantel"), "Run Mantel Test",
                             icon = icon("random"), class = "btn-action-primary btn-block",
                             style = "font-weight:bold;")
              )
            ),
            tags$hr(),
            fluidRow(
              column(3, valueBoxOutput(ns("box_m_stat"), width = NULL)),
              column(3, valueBoxOutput(ns("box_m_pval"), width = NULL)),
              column(3, valueBoxOutput(ns("box_m_n"),    width = NULL)),
              column(3, valueBoxOutput(ns("box_m_r2"),   width = NULL))
            ),
            uiOutput(ns("ui_mantel_summary"))
          )
        ),

        fluidRow(
          box(width = 6, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("chart-line"), " Scatter plot"),
            plotly::plotlyOutput(ns("mantel_scatter"), height = "360px")
          ),
          box(width = 6, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("chart-area"), " Permutation distribution"),
            plotly::plotlyOutput(ns("mantel_hist"), height = "360px")
          )
        ),

        fluidRow(
          box(width = 12, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("table"), " Data used in the last Mantel run"),
            DT::DTOutput(ns("dt_mantel_data")),
            tags$br(),
            downloadButton(ns("dl_mantel_csv"), "Download data used", class = "btn-action-secondary btn-sm")
          )
        )
      )
    )
  )
}
