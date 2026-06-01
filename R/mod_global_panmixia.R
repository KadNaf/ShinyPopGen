# mod_global_panmixia.R
# Tab: Global panmixia
# Overall HWE across all populations - FIT (Weir & Cockerham): bootstrap CI + permutation test.
# Golem module UI - server: server_general_stats("general_stats", rv)

mod_global_panmixia_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    tags$head(gs_head()),

    module_banner("globe", "Global Panmixia \u2014 FIT",
      "Overall HWE across all populations \u00b7 Weir & Cockerham (1984) \u00b7 Block bootstrap + permutation p-value",
      "#E1AF00"),
    tags$div(class = "spg-method-note", style = "border-left-color:#E1AF00;",
      HTML(paste0(
        "Global panmixia: the <em>entire</em> dataset is at HWE \u2014 all populations mate as a single unit. ",
        "<br><br>",
        "<b>H<sub>0</sub>:</b> FIT = 0 globally (no departure from HWE across all samples). &nbsp;",
        "<b>Bootstrap:</b> populations are the resampling unit (block bootstrap); percentile CI per locus. &nbsp;",
        "<b>Permutation:</b> alleles shuffled across all individuals (ignoring populations); two-sided |FIT| test."
      ))
    ),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("globe-americas"),
                    "FIT: CI & p-value parameters"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            h4(icon("sliders"), "Parameters"),
            numericInput(ns("n_perm_fit"),    "Number of Permutations:",        value = 5000, min = 100, max = 20000, step = 100),
            numericInput(ns("n_boot_fit"),    "Number of Bootstrap Replicates:", value = 5000, min = 100, max = 20000, step = 100),
            numericInput(ns("conf_level_fit"),"Confidence Level:",               value = 0.95, min = 0.80, max = 0.99, step = 0.01),
            actionButton(ns("Run_FIT_Analysis"), "Run FIT Analysis",
                         icon = icon("rocket"),
                         class = "btn-action-primary btn-block", style = "font-weight: bold;")
          ),
          column(9,
            h4(icon("chart-line"), "FIT Analysis Summary",
               style = "font-weight: 600; color: #2c3e50; margin-bottom: 15px;"),
            fluidRow(
              column(3,
                valueBoxOutput(ns("global_fit_box"),       width = NULL),
                valueBoxOutput(ns("fit_ci_width_box"),     width = NULL)
              ),
              column(3,
                valueBoxOutput(ns("global_fit_pvalue_box"),width = NULL),
                valueBoxOutput(ns("fit_power_box"),        width = NULL)
              ),
              column(3,
                valueBoxOutput(ns("significant_loci_fit_box"),  width = NULL),
                valueBoxOutput(ns("fit_convergence_box"),       width = NULL)
              ),
              column(3,
                valueBoxOutput(ns("analysis_time_fit_box"),width = NULL),
                valueBoxOutput(ns("fit_quality_box"),      width = NULL)
              )
            ),
            fluidRow(
              column(12,
                h5("Analysis Progress", style = "margin-top: 15px; font-weight: 600;"),
                shinyWidgets::progressBar(id = ns("fit_progress"), value = 0,
                                          title = "Overall Progress")
              )
            )
          )
        )
      )
    ),

    h2("FIT \u2014 Bootstrap CI and permutation results", class = "section-title"),
    tags$p(HTML(paste0(
      "Population-block bootstrap confidence intervals (populations are the resampling unit) and ",
      "permutation p-values for the global HWE test across <em>all</em> populations combined. ",
      "<br>Bootstrap CI and p-values are given per locus; the Overall row uses the ratio-of-sums FIT across all loci. ",
      "<br>A CI excluding zero or a small p-value indicates a significant global departure from HWE."
    )), style = "font-size: 16px; line-height: 1.5; color: #2c3e50;"),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("table"),
                    "Results"),
        solidHeader = TRUE, status = "primary",
        tabsetPanel(
          tabPanel("P-value and confidence intervals",
            h4(icon("info-circle"), "FIT estimates with bootstrap CI and permutation p-values"),
            p("FIT estimates per locus. Bootstrap CI: population-block resampling (populations resampled with replacement).
              Permutation p-values: global allele shuffle, two-sided |FIT| test, consistent with the FIS permutation test."),
            DTOutput(ns("fit_results_table")), br(),
            fluidRow(
              column(6, downloadButton(ns("download_fit_table"),     ".csv", class = "btn-download-primary btn-block")),
              column(6, downloadButton(ns("download_fit_table_txt"), ".txt", class = "btn-download-secondary btn-block"))
            )
          ),
          tabPanel("Visualization",
            h4(icon("chart-line"), "FIT estimates by locus"),
            plotOutput(ns("fit_plot"), height = "400px"), br(),
            downloadButton(ns("download_fit_plot"), ".png", class = "btn-download-primary")
          ),
        ),
        style = "padding: 10px;"
      )
    )
  )
}
