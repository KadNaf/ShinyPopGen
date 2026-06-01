# ui_LD.R  (MODULE UI)

linkage_desequilibrium_UI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    useWaiter(),
    tags$head(
      tags$style(HTML("
        .ld-heatmap-container {
          overflow-x: auto;
          overflow-y: auto;
          max-height: 600px;
        }
      "))
    ),
    
    module_banner("link", "Linkage Disequilibrium",
      "Pairwise LD among all loci \u00b7 Permutation p-values",
      "#EBCC2A"),
    
    fluidRow(
      box(
        width = 4,
        title = div(
          style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
          icon("cogs"), "Linkage disequilibrium parameters"
        ),
        solidHeader = TRUE,
        status = "primary",
        
        h4(icon("sliders"), "Parameters"),
        checkboxInput(ns("include_missing"), "Include Missing Data", value = TRUE),
        numericInput(ns("n_iterations"), "Number of Permutations:",
                     value = 10000, min = 1000, max = 100000, step = 1000),
        
        tags$hr(),
        actionButton(ns("run_LD"), "Run LD Analysis",
                     icon = icon("rocket"),
                     class = "btn-action-primary btn-block")
      ),
      
      box(
        width = 8,
        title = div(
          style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
          icon("chart-line"), "Analysis Summary"
        ),
        solidHeader = TRUE,
        status = "primary",
        
        fluidRow(
          column(3, valueBoxOutput(ns("total_pairs_box"), width = NULL)),
          column(3, valueBoxOutput(ns("significant_pairs_box"), width = NULL)),
          column(3, valueBoxOutput(ns("mean_pvalue_box"), width = NULL)),
          column(3, valueBoxOutput(ns("analysis_time_ld_box"), width = NULL))
        ),
        
        br(),
        
        fluidRow(
          column(
            12,
            h5("Analysis Progress", style = "margin-top: 15px; font-weight: 600;"),
            shinyWidgets::progressBar(id = ns("LD_progress"), value = 0, title = "Overall Progress")
          )
        ),
        
        br(),
        
        div(
          style = paste(
            "margin-top:12px; padding:10px 14px;",
            "background:#D9D0D3; border-left:4px solid #8D8680;",
            "border-radius:4px; font-size:13px; line-height:1.6;"
          ),
          h5(icon("calculator"), "Additive property of the G-statistic",
             style = "margin-top:0; color:#39312F; font-weight:600;"),
          tags$p(
            "The G-statistic is ", tags$strong("additive across subsamples"),
            " (populations): the global G over all populations equals the sum ",
            "of the per-population G values for the same locus pair."
          ),
          tags$p(HTML(
            "<b>G<sub>global</sub> = G<sub>pop1</sub> + G<sub>pop2</sub> + \u2026 + G<sub>popk</sub></b>"
          )),
          tags$p(
            "The permutation null distribution of G\u2090\u2097\u2092\u2071\u2090\u2097 is built ",
            "by summing the permuted G values ", tags$em("from the same replicate"),
            " across populations. This yields a valid global p-value that accounts ",
            "for population structure, without requiring independence between ",
            "subsamples."
          ),
          tags$p(
            style = "margin-bottom:0; color:#555;",
            "Note: p-values are Monte Carlo estimates ",
            HTML("p = (n<sub>\u2265obs</sub> + 1) / B"),
            ". At least ", tags$strong("B = 1\u202f000 permutations"),
            " are recommended; B \u2265 10\u202f000 for publication."
          )
        )
      )
    ),
    
    # ===== SECTION 2: RESULTS TABLE =====
    h2("Linkage Disequilibrium Results", class = "section-title"),
    
    fluidRow(
      box(
        width = 12,
        title = div(
          style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
          icon("table"), "LD P-values Table"
        ),
        solidHeader = TRUE,
        status = "primary",
        
        fluidRow(
          column(
            3,
            h4(icon("cog"), "Table Options"),
            selectInput(ns("table_view"), "Display Mode:",
                        choices = c("All Pairs" = "all",
                                    "Significant Only (p < 0.05)" = "sig_05",
                                    "Highly Significant (p < 0.01)" = "sig_01",
                                    "Very Highly Significant (p < 0.001)" = "sig_001")),
            selectInput(ns("sort_by_ld"), "Sort By:",
                        choices = c("Locus Pair" = "pair",
                                    "P-value (Ascending)" = "pval_asc",
                                    "P-value (Descending)" = "pval_desc")),
            numericInput(ns("decimal_places"), "Decimal Places:",
                         value = 5, min = 2, max = 10, step = 1),
            checkboxInput(ns("highlight_sig"), "Highlight Significant", TRUE),
            downloadButton(ns("download_LD_csv"), ".csv", class = "btn-download-primary"),
            downloadButton(ns("download_LD_txt"), ".txt", class = "btn-download-secondary")
          ),
          column(
            9,
            DTOutput(ns("summary_output")),
            br(),
            verbatimTextOutput(ns("table_summary_stats"))
          )
        )
      )
    )
  )
}
