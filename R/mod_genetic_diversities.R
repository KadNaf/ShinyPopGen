# mod_genetic_diversities.R
# Tab: Genetic diversities
# HS, HT (within- and total-gene diversity) + locus bootstrap for all multilocus estimators.
# Results are populated by running FST analysis from the Subdivision tab.
# Golem module UI - server: server_general_stats("general_stats", rv)

mod_genetic_diversities_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    tags$head(gs_head()),

    module_banner("chart-line", "Genetic Diversities \u2014 HS, HT",
      "Within- and total gene diversity \u00b7 Locus bootstrap for all multilocus WC84 estimators",
      "#78B7C5"),
    tags$div(class = "spg-method-note", style = "border-left-color:#78B7C5;",
      HTML(paste0(
        "<b>HS</b> (within-population gene diversity) and <b>HT</b> (total gene diversity) ",
        "from Weir &amp; Cockerham (1984), reported per locus and as multilocus estimates. ",
        "<br><br>",
        "<b>Confidence intervals are computed by three resampling schemes:</b>",
        "<ul style='margin:4px 0 0 16px;'>",
        "<li><b>Individuals</b> (HS per locus and per population): individuals resampled with replacement within each population.</li>",
        "<li><b>Populations</b> (HS and HT per locus): populations resampled with replacement.</li>",
        "<li><b>Loci</b> (overall HS and HT only): loci resampled with replacement.</li>",
        "</ul>"
      ))
    ),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("chart-line"),
                    "Genetic Diversity Analysis parameters"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            h4(icon("sliders"), "Parameters"),
            numericInput(ns("n_perm_fst_div"),    "Number of Permutations:",        value = 5000, min = 100, max = 20000, step = 100),
            numericInput(ns("n_boot_fst_div"),    "Number of Bootstrap Replicates:", value = 5000, min = 100, max = 20000, step = 100),
            actionButton(ns("run_FST_Analysis_div"), "Run Diversity Analysis",
                         icon = icon("rocket"),
                         class = "btn-action-primary btn-block", style = "font-weight: bold;"),
            tags$small(
              style = "color: #666; margin-top: 6px; display: block;",
              icon("info-circle"),
              "Also populates FST/FIT results in the Population subdivision tab."
            )
          ),
          column(9,
            h4(icon("chart-line"), "Diversity Analysis Summary",
               style = "font-weight: 600; color: #2c3e50; margin-bottom: 15px;"),
            fluidRow(
              column(3, valueBoxOutput(ns("global_fst_div_box"),    width = NULL)),
              column(3, valueBoxOutput(ns("global_hs_box"),         width = NULL)),
              column(3, valueBoxOutput(ns("global_ht_box"),         width = NULL)),
              column(3, valueBoxOutput(ns("analysis_time_div_box"), width = NULL))
            ),
            fluidRow(
              column(12,
                h5("Analysis Progress", style = "margin-top: 15px; font-weight: 600;"),
                shinyWidgets::progressBar(id = ns("fst_progress_div"), value = 0,
                                          title = "Overall Progress")
              )
            )
          )
        )
      )
    ),

    h2("HS, HT \u2014 per-locus results and locus bootstrap", class = "section-title"),
    tags$p(HTML(paste0(
      "<b>HS</b> (within-population gene diversity) and <b>HT</b> (total gene diversity) per locus. ",
      "The <em>locus bootstrap</em> resamples L loci with replacement (B replicates) to obtain ",
      "SE and percentile CI for the multilocus estimators: ",
      "FST, FIT, FIS, HS and HT. ",
      "<br>Results are produced by the run in the <b>Subdivision</b> tab or by the independent run above."
    )), style = "font-size: 16px; line-height: 1.5; color: #2c3e50;"),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("dna"), "Results"),
        solidHeader = TRUE, status = "primary",
        tabsetPanel(
          tabPanel("HS per locus \u2014 individuals",
            br(),
            h4(icon("info-circle"), "HS per locus"),
            p("Confidence interval obtained by resampling individuals with replacement within each population."),
            DTOutput(ns("hs_indiv_table")), br()
          ),
          tabPanel("HS per locus \u2014 populations",
            br(),
            h4(icon("info-circle"), "HS per locus"),
            p("Confidence interval obtained by resampling populations with replacement."),
            DTOutput(ns("hs_pop_table")), br(),
            fluidRow(
              column(6, downloadButton(ns("download_hs_table"),     ".csv", class = "btn-download-primary btn-block")),
              column(6, downloadButton(ns("download_hs_table_txt"), ".txt", class = "btn-download-secondary btn-block"))
            )
          ),
          tabPanel("HS per population",
            br(),
            h4(icon("users"), "HS per population"),
            p("HS per population, averaged across loci. Confidence interval obtained by resampling individuals with replacement within each population."),
            DTOutput(ns("hs_per_pop_table")), br()
          ),
          tabPanel("Overall HS \u2014 loci",
            br(),
            h4(icon("retweet"), "Overall HS"),
            p("Confidence interval obtained by resampling loci with replacement. Applies to the overall multilocus HS only."),
            DTOutput(ns("hs_locus_table")), br()
          ),
          tabPanel("HS visualization",
            h4(icon("chart-line"), "HS per locus \u2014 CI from resampling populations"),
            plotOutput(ns("hs_plot"), height = "400px"), br(),
            downloadButton(ns("download_hs_plot"), ".png", class = "btn-download-primary")
          ),
          tabPanel("HT results",
            br(),
            h4(icon("info-circle"), "Total gene diversity (HT)"),
            p("HT per locus. Confidence interval obtained by resampling populations with replacement. The Overall row also shows CI from resampling loci with replacement."),
            DTOutput(ns("ht_results_table")), br(),
            fluidRow(
              column(6, downloadButton(ns("download_ht_table"),     ".csv", class = "btn-download-primary btn-block")),
              column(6, downloadButton(ns("download_ht_table_txt"), ".txt", class = "btn-download-secondary btn-block"))
            )
          ),
          tabPanel("HT visualization",
            h4(icon("chart-line"), "HT estimates by locus"),
            plotOutput(ns("ht_plot"), height = "400px"), br(),
            downloadButton(ns("download_ht_plot"), ".png", class = "btn-download-primary")
          ),
          tabPanel("Locus bootstrap",
            h4(icon("retweet"), "Multilocus estimators \u2014 locus bootstrap"),
            p("Bootstrap SE and percentile CI for FST, FIT, FIS, HS, HT
              computed by resampling L loci with replacement (B replicates)."),
            DTOutput(ns("locus_boot_table"))
          )
        ),
        style = "padding: 10px;"
      )
    )
  )
}
