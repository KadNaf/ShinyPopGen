# mod_subdivision.R
# Tab: Subdivision
# Population differentiation - FST (Weir & Cockerham): bootstrap CI + permutation test.
# Note: HS and HT are also computed here but displayed in the Genetic diversities tab.
# Golem module UI - server: server_general_stats("general_stats", rv)

mod_subdivision_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    tags$head(gs_head()),

    module_banner("sitemap", "Population Subdivision \u2014 FST",
      "Population differentiation \u00b7 Weir & Cockerham (1984) \u00b7 Block bootstrap CI + permutation p-value",
      "#B40F20"),
    tags$div(class = "spg-method-note", style = "border-left-color:#B40F20;",
      HTML(paste0(
        "Population subdivision: allele-frequency differences among populations (FST > 0). ",
        "HS and HT are also computed here \u2014 see the <b>Diversities</b> tab for locus bootstrap. ",
        "<br><br>",
        "<b>H<sub>0</sub>:</b> FST = 0 (no differentiation). &nbsp;",
        "<b>Bootstrap:</b> population-block resampling; percentile CI. &nbsp;",
        "<b>Permutation:</b> genotypes randomly reassigned among populations; one-sided test."
      ))
    ),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("sitemap"),
                    "FST: CI & p-value parameters"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            h4(icon("sliders"), "Parameters"),
            numericInput(ns("n_perm_fst"),    "Number of Permutations:",        value = 5000, min = 100, max = 20000, step = 100),
            numericInput(ns("n_boot_fst"),    "Number of Bootstrap Replicates:", value = 5000, min = 100, max = 20000, step = 100),
            numericInput(ns("conf_level_fst"),"Confidence Level:",               value = 0.95, min = 0.80, max = 0.99, step = 0.01),
            actionButton(ns("run_FST_Analysis"), "Run FST Analysis",
                         icon = icon("rocket"),
                         class = "btn-action-primary btn-block", style = "font-weight: bold;"),
            tags$small(
              style = "color: #666; margin-top: 6px; display: block;",
              icon("info-circle"),
              "Also computes HS, HT and locus bootstrap (see Genetic diversities tab)."
            )
          ),
          column(9,
            h4(icon("chart-line"), "FST Analysis Summary",
               style = "font-weight: 600; color: #2c3e50; margin-bottom: 15px;"),
            fluidRow(
              column(3,
                valueBoxOutput(ns("global_fst_box"),       width = NULL),
                valueBoxOutput(ns("fst_ci_width_box"),     width = NULL)
              ),
              column(3,
                valueBoxOutput(ns("global_fst_pvalue_box"),width = NULL),
                valueBoxOutput(ns("fst_power_box"),        width = NULL)
              ),
              column(3,
                valueBoxOutput(ns("significant_loci_fst_box"),  width = NULL),
                valueBoxOutput(ns("fst_convergence_box"),       width = NULL)
              ),
              column(3,
                valueBoxOutput(ns("analysis_time_fst_box"),width = NULL),
                valueBoxOutput(ns("fst_quality_box"),      width = NULL)
              )
            ),
            fluidRow(
              column(12,
                h5("Analysis Progress", style = "margin-top: 15px; font-weight: 600;"),
                shinyWidgets::progressBar(id = ns("fst_progress"), value = 0,
                                          title = "Overall Progress")
              )
            )
          )
        )
      )
    ),

    h2("FST \u2014 Bootstrap CI and permutation results", class = "section-title"),
    tags$p(HTML(paste0(
      "FST per locus with population-block bootstrap confidence intervals. ",
      "Permutation p-values derived from shuffling population labels (one-sided test, FST &ge; observed). ",
      "<br>HS and HT computed in the same run are reported in the ",
      "<b>Genetic Diversities</b> tab."
    )), style = "font-size: 16px; line-height: 1.5; color: #2c3e50;"),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("table"),
                    "Results"),
        solidHeader = TRUE, status = "primary",
        tabsetPanel(
          tabPanel("FST results",
            h4(icon("info-circle"), "Bootstrap confidence intervals"),
            p("FST per locus with population-block bootstrap CI and
              permutation p-values (population labels shuffled, one-sided test)."),
            DTOutput(ns("fst_results_table")), br(),
            fluidRow(
              column(6, downloadButton(ns("download_fst_table"),     ".csv", class = "btn-download-primary btn-block")),
              column(6, downloadButton(ns("download_fst_table_txt"), ".txt", class = "btn-download-secondary btn-block"))
            )
          ),
          tabPanel("Visualization",
            h4(icon("chart-line"), "FST estimates by locus"),
            plotOutput(ns("fst_plot"), height = "400px"), br(),
            downloadButton(ns("download_fst_plot"), ".png", class = "btn-download-primary")
          ),
        ),
        style = "padding: 10px;"
      )
    )
  )
}
