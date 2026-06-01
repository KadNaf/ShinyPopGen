# mod_local_panmixia.R
# Tab: Local panmixia
# Within-population HWE - FIS (Weir & Cockerham): bootstrap CI + permutation test.
# Golem module UI - server: server_general_stats("general_stats", rv)

mod_local_panmixia_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    tags$head(gs_head()),

    module_banner("flask", "Local Panmixia \u2014 FIS",
      "Within-population HWE \u00b7 Weir & Cockerham (1984) \u00b7 Bootstrap CI + permutation p-value",
      "#9986A5"),
    tags$div(class = "spg-method-note", style = "border-left-color:#9986A5;",
      HTML(paste0(
        "Local panmixia means that each sub-population is at ",
        "Hardy-Weinberg equilibrium (HWE) \u2014 individuals mate randomly ",
        "<em>within</em> their population. ",
        "<br><br>",
        "<b>H<sub>0</sub>:</b> FIS = 0 within each population (no departure from HWE). &nbsp;",
        "<b>Bootstrap:</b> individuals resampled with replacement within populations; percentile CI. &nbsp;",
        "<b>Permutation:</b> alleles reshuffled within each population; two-sided |FIS| test."
      ))
    ),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("rocket"),
                    "FIS: CI & p-value parameters"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            h4(icon("sliders"), "Parameters"),
            numericInput(ns("n_perm"),    "Number of Permutations:",        value = 5000, min = 100, max = 20000),
            numericInput(ns("n_boot"),    "Number of Bootstrap Replicates:", value = 5000, min = 100, max = 20000),
            numericInput(ns("conf_level"),"Confidence Level:",               value = 0.95, min = 0.80, max = 0.99, step = 0.01),
            selectInput(ns("analysis_level"), "Analysis Level:",
                        choices = c("By Locus", "By Population"), selected = "By Locus"),
            actionButton(ns("Run_FIS_Analysis"), "Run FIS Analysis",
                         icon = icon("rocket"),
                         class = "btn-action-primary btn-block", style = "font-weight: bold;")
          ),
          column(9,
            h4(icon("chart-line"), "Analysis Summary",
               style = "font-weight: 600; color: #2c3e50; margin-bottom: 15px;"),
            fluidRow(
              column(3, valueBoxOutput(ns("global_fis_box"),        width = NULL)),
              column(3, valueBoxOutput(ns("global_pvalue_box"),     width = NULL)),
              column(3, valueBoxOutput(ns("significant_loci_box"),  width = NULL)),
              column(3, valueBoxOutput(ns("analysis_time_box"),     width = NULL))
            ),
            fluidRow(
              column(12,
                h5("Analysis Progress", style = "margin-top: 15px; font-weight: 600;"),
                shinyWidgets::progressBar(id = ns("fis_progress"), value = 0,
                                          title = "Overall Progress")
              )
            )
          )
        )
      )
    ),

    h2("FIS \u2014 Bootstrap CI and permutation results", class = "section-title"),
    tags$p(HTML(paste0(
      "Bootstrap confidence intervals derived from resampling individuals within populations. ",
      "Permutation p-values from allele shuffling within populations (two-sided |FIS| test). ",
      "<br>A CI excluding zero indicates a significant departure from HWE at that locus / population."
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
            h4(icon("info-circle"),
               "FIS estimates with bootstrap CI and permutation p-values"),
            p("Bootstrap CI derived from resampling individuals within populations.
              P-values from allele permutation within populations (two-sided |FIS| test).
              CI excluding zero indicates significant departure from HWE."),
            DTOutput(ns("fis_results_table")), br(),
            fluidRow(
              column(6, downloadButton(ns("download_fis_table"),     ".csv", class = "btn-download-primary btn-block")),
              column(6, downloadButton(ns("download_fis_table_txt"), ".txt", class = "btn-download-secondary btn-block"))
            )
          ),
          tabPanel("Visualization",
            h4(icon("chart-line"), "Bootstrap-based FIS inference"),
            plotOutput(ns("fis_plot"), height = "400px"), br(),
            downloadButton(ns("download_fis_plot"), ".png", class = "btn-download-primary")
          ),
        ),
        style = "padding: 10px;"
      )
    ),

    h2("FIS \u2014 By Locus \u00d7 Population", class = "section-title"),
    tags$p(HTML(paste0(
      "WC84 FIS and permutation p-values for every locus \u00d7 population combination. ",
      "Permutation only (no bootstrap CI). Run independently of the main analysis above."
    )), style = "font-size: 16px; line-height: 1.5; color: #2c3e50;"),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("table"),
                    "FIS (WC84) per locus \u00d7 population"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            h4(icon("sliders"), "Parameters"),
            numericInput(ns("fis_lp_n_perm"), "Number of Permutations:",
                         value = 5000, min = 100, max = 20000)
          ),
          column(9,
            br(),
            actionButton(ns("run_fis_locus_pop"),
                         label = tagList(icon("calculator"), tags$strong("Compute")),
                         class = "btn-action-primary")
          )
        ),
        br(), br(),
        h5("Observed FIS (WC84)"),
        DTOutput(ns("fis_locus_pop_obs")), br(),
        h5("Permutation p-values (two-sided)"),
        DTOutput(ns("fis_locus_pop_pval")), br(),
        fluidRow(
          column(6, downloadButton(ns("download_fis_locus_pop"),     ".csv", class = "btn-download-primary btn-block")),
          column(6, downloadButton(ns("download_fis_locus_pop_txt"), ".txt", class = "btn-download-secondary btn-block"))
        ),
        style = "padding: 10px;"
      )
    )
  )
}
