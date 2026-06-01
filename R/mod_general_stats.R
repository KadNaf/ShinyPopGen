# mod_general_stats.R
# Tab: General stats
# Basic observed statistics + per-allele F-statistics.
# Golem module UI — server: server_general_stats("general_stats", rv)

mod_general_stats_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    useWaiter(),
    tags$head(gs_head()),

    module_banner("table", "General Statistics",
      "Na \u00b7 Ne \u00b7 Ho \u00b7 He \u00b7 F-statistics per allele (Weir & Cockerham 1984)",
      "#3B9AB2"),

    fluidRow(
      box(
        width = 3,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("chart-bar"), "Statistics selection"),
        solidHeader = TRUE, status = "primary",
        h4("Select Genetic Indices"),
        h5("Diversity Measures:"),
        checkboxInput(ns("ho_checkbox"),  "Ho (Observed Heterozygosity)", TRUE),
        checkboxInput(ns("hs_checkbox"),  "Hs (Expected Heterozygosity within populations)", TRUE),
        checkboxInput(ns("ht_checkbox"),  "Ht (Total expected heterozygosity)", TRUE),
        h5("F-statistics (Weir & Cockerham):"),
        checkboxInput(ns("fit_wc_checkbox"), "FIT (Weir & Cockerham estimator)", TRUE),
        checkboxInput(ns("fis_wc_checkbox"), "FIS (Weir & Cockerham estimator)", TRUE),
        checkboxInput(ns("fst_wc_checkbox"), "FST (Weir & Cockerham estimator)", TRUE),
        h5("Advanced Statistics:"),
        checkboxInput(ns("fst_max_checkbox"),  "Fst-max (Maximum differentiation, Meirmans)", FALSE),
        checkboxInput(ns("fst_prim_checkbox"), "Fst' (Meirmans) (Empirical standardisation)", FALSE),
        checkboxInput(ns("GST_checkbox"),      "GST (Nei's genetic differentiation)", FALSE),
        checkboxInput(ns("GST_sec_checkbox"),  "GST'' (Hedrick's correction)", FALSE),
        tags$hr(),
        actionButton(ns("run_basic_stats"),
                     label = tagList(icon("calculator"), tags$strong("Compute Statistics")),
                     class = "btn-action-primary btn-block")
      ),
      box(
        width = 9,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("table"), "General statistics"),
        solidHeader = TRUE, status = "primary",
        tabsetPanel(
          tabPanel("Statistic estimates",
            DTOutput(ns("basic_stats_table")), br(),
            fluidRow(
              column(6, downloadButton(ns("download_basic_stats"),     ".csv", class = "btn-download-primary btn-block")),
              column(6, downloadButton(ns("download_basic_stats_txt"), ".txt", class = "btn-download-secondary btn-block"))
            )
          ),
          tabPanel("Gene Diversity by Population",
            h5("Expected heterozygosity (Hs) per locus and population",
               style = "margin-top: 10px;"),
            DTOutput(ns("gene_diversity_table")), br(),
            fluidRow(
              column(6, downloadButton(ns("download_gene_diversity"),     ".csv", class = "btn-download-primary btn-block")),
              column(6, downloadButton(ns("download_gene_diversity_txt"), ".txt", class = "btn-download-secondary btn-block"))
            )
          ),
          tabPanel("By Population",
            h5("All populations — Ho, Hs, Fis (WC) averaged over loci",
               style = "margin-top: 10px;"),
            tableOutput(ns("overall_by_pop")), br(),
            fluidRow(
              column(6, downloadButton(ns("download_overall_by_pop"),     ".csv", class = "btn-download-primary btn-block")),
              column(6, downloadButton(ns("download_overall_by_pop_txt"), ".txt", class = "btn-download-secondary btn-block"))
            ),
            tags$hr(),
            h5("Per-locus detail for selected population"),
            selectInput(ns("selected_pop_overall"), "Select Population:", choices = NULL),
            DTOutput(ns("basic_stats_by_pop_selected")), br(),
            fluidRow(
              column(6, downloadButton(ns("download_pop_stats"),     ".csv", class = "btn-download-primary btn-block")),
              column(6, downloadButton(ns("download_pop_stats_txt"), ".txt", class = "btn-download-secondary btn-block"))
            )
          )
        ),
        style = "overflow-y: auto; max-height: 600px; padding: 10px;"
      )
    ),

    h2("F-statistics per allele (Weir & Cockerham)", class = "section-title"),
    tags$p(HTML(paste0(
      "For each allele at each locus, WC84 variance components ",
      "(a = between-pop, b = between-indiv, c = within-indiv) are computed ",
      "and the three F-statistics derived: ",
      "<b>FIS</b> = b/(b+c), <b>FST</b> = a/(a+b+c), <b>FIT</b> = (a+b)/(a+b+c). ",
      "<br>High <b>FIS</b> for a specific allele may indicate amplification dropout. ",
      "Outlier <b>FST</b> may signal selection or local adaptation."
    )), style = "font-size: 16px; line-height: 1.5; color: #2c3e50;"),

    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("dna"), "F-statistics per allele (Weir & Cockerham)"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            actionButton(ns("compute_allele_fstats"),
              label = tagList(icon("calculator"),
                              tags$strong("Compute F-statistics per allele")),
              class = "btn-action-primary btn-block", style = "font-weight: bold;")
          ),
          column(9,
            h5("Results table"),
            p("FIS, FST and FIT are reported per allele and per locus from WC84
              variance components. Run independently of the bootstrap/permutation
              analyses.")
          )
        ),
        br(),
        DTOutput(ns("fis_allele_table")), br(),
        fluidRow(
          column(6, downloadButton(ns("download_fis_allele_table"),     ".csv", class = "btn-download-primary btn-block")),
          column(6, downloadButton(ns("download_fis_allele_table_txt"), ".txt", class = "btn-download-secondary btn-block"))
        ),
        style = "padding: 10px;"
      )
    )
  )
}
