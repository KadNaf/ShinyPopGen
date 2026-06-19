# ui_isolation_by_distance.R
# Two tabs:
#   1. Pairwise Distances (FreeNA: FST WC84 + Cavalli-Sforza + ENA/INA + bootstrap)
#   2. Mantel Test (on matrices or rectangular/long format)

isolation_by_distance_UI <- function(id) {
  ns <- NS(id)

  fluidPage(
    tags$head(gs_head()),

    module_banner(
      "map-marker-alt",
      "Isolation by Distance",
      "FreeNA ENA/INA correction \u00b7 Bootstrap on loci \u00b7 Mantel test",
      "#2CBF9F"
    ),

    tags$div(
      style = paste(
        "display:flex; align-items:center; gap:12px;",
        "background:#FFF8E1; border:2px solid #E1AF00;",
        "border-radius:6px; padding:10px 16px; margin-bottom:16px;"
      ),
      tags$span(style = "font-size:1.8em; line-height:1;", "\U0001f6a7"),
      tags$div(
        tags$strong(style = "color:#7B5800;", "Module under construction"),
        tags$span(style = "color:#7B5800; margin-left:8px; font-size:0.9em;",
          "Results are functional but the module is still being validated. Use with caution.")
      )
    ),

    tabsetPanel(
      id = ns("main_tabs"),
      type = "tabs",

      # ══════════════════════════════════════════════════════════════════════
      # TAB 1: Pairwise Distances (FreeNA)
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = "Pairwise Distances",
        icon = icon("project-diagram"),
        tags$div(
          class = "spg-method-note", style = "border-left-color:#2CBF9F;",
          HTML(paste0(
            "Computes pairwise F<sub>ST</sub> (Weir & Cockerham 1984) and Cavalli-Sforza & Edwards (1967) distances ",
            "between populations. Null allele frequencies are estimated with the EM algorithm of ",
            "Dempster, Laird & Rubin (1977). The <b>ENA correction</b> (Chapuis & Estoup 2007) is applied ",
            "to F<sub>ST</sub>; the <b>INA correction</b> to Cavalli-Sforza distances. ",
            "Bootstrap resampling over <b>loci</b> provides 95% confidence intervals. ",
            "<b>Requires genotype data with population assignments (imported into DuckDB).</b><br><br>",
            "<b>References:</b> Chapuis & Estoup 2007 <em>Mol Ecol Notes</em> 7:221\u2013231 \u00b7 ",
            "Dempster et al. 1977 <em>J R Stat Soc B</em> 39:1\u201338 \u00b7 ",
            "Weir & Cockerham 1984 <em>Evolution</em> 38:1358\u20131370"
          ))
        ),

        # Configuration
        fluidRow(
          box(
            width = 12,
            title = div(
              style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
              icon("sliders-h"), " FreeNA parameters"
            ),
            solidHeader = TRUE, status = "primary",
            fluidRow(
              column(3,
                numericInput(ns("n_boot_loci"), "Bootstrap replicates (loci):",
                             value = 1000, min = 100, max = 10000, step = 100),
                checkboxInput(ns("calc_fst"),  "Compute FST (WC84)", value = TRUE),
                checkboxInput(ns("ena_corr"),  "Apply ENA correction (FST)", value = TRUE),
                checkboxInput(ns("calc_cs"),   "Compute Cavalli-Sforza distance", value = TRUE),
                checkboxInput(ns("ina_corr"),  "Apply INA correction (CS)", value = TRUE),
                tags$hr(),
                actionButton(
                  ns("run_freena"), "Run FreeNA Analysis",
                  icon  = icon("calculator"),
                  class = "btn-action-primary btn-block",
                  style = "font-weight:bold;"
                )
              ),
              column(9,
                h4(icon("chart-bar"), "Summary",
                   style = "font-weight:600; color:#2c3e50; margin-bottom:15px;"),
                fluidRow(
                  column(3, valueBoxOutput(ns("box_ninds"),    width = NULL)),
                  column(3, valueBoxOutput(ns("box_npops"),    width = NULL)),
                  column(3, valueBoxOutput(ns("box_nloci"),    width = NULL)),
                  column(3, valueBoxOutput(ns("box_npairs"),   width = NULL))
                ),
                tags$h5("Null allele frequencies (EM Dempster)",
                        style = "font-weight:600; margin-top:14px; color:#2c3e50;"),
                DT::DTOutput(ns("rd_table"))
              )
            )
          )
        ),

        # Pairwise matrices
        fluidRow(
          column(6,
            box(
              width = 12,
              title = div(
                style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
                icon("th"), " Pairwise FST matrix"
              ),
              solidHeader = FALSE,
              radioButtons(ns("fst_matrix_choice"), "Which FST:",
                           choices = c("Uncorrected" = "fst", "ENA-corrected" = "fst_ena"),
                           selected = "fst_ena", inline = TRUE),
              DT::DTOutput(ns("fst_matrix_table"))
            )
          ),
          column(6,
            box(
              width = 12,
              title = div(
                style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
                icon("th"), " Cavalli-Sforza distance matrix"
              ),
              solidHeader = FALSE,
              radioButtons(ns("cs_matrix_choice"), "Which CS:",
                           choices = c("Uncorrected" = "cs", "INA-corrected" = "cs_ina"),
                           selected = "cs_ina", inline = TRUE),
              DT::DTOutput(ns("cs_matrix_table"))
            )
          )
        ),

        # Detailed pairwise table with CI
        fluidRow(
          box(
            width = 12,
            title = div(
              style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
              icon("list"), " Detailed pairwise results with 95% bootstrap CI"
            ),
            solidHeader = FALSE,
            DT::DTOutput(ns("pairwise_detail_table")),
            tags$br(),
            fluidRow(
              column(6, downloadButton(ns("dl_pairwise_csv"), "Download pairwise table (CSV)",
                                       class = "btn-action-secondary btn-sm btn-block")),
              column(6, downloadButton(ns("dl_freena_log"),   "Download FreeNA log (TXT)",
                                       class = "btn-action-secondary btn-sm btn-block"))
            )
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════════
      # TAB 2: Mantel Test
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(
        title = "Mantel Test",
        icon = icon("chart-line"),
        tags$div(
          class = "spg-method-note", style = "border-left-color:#2CBF9F;",
          HTML(paste0(
            "Mantel test (Mantel 1967) assesses the correlation between two distance matrices. ",
            "Significance is tested by random permutation of rows/columns of one matrix. ",
            "Supports <b>square matrices</b> or <b>rectangular (column-wise)</b> format as in Fstat/RT: ",
            "first row = population pairs, subsequent rows = distance values (e.g., from different loci or measures). ",
            "For IBD, one matrix is genetic (computed in Tab 1 or uploaded), the other is geographic ",
            "(from GPS coordinates or uploaded).<br><br>",
            "<b>Reference:</b> Mantel N. 1967. <em>Math Geol</em> 15:65\u201375."
          ))
        ),

        fluidRow(
          box(
            width = 12,
            title = div(
              style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
              icon("sliders-h"), " Mantel test parameters"
            ),
            solidHeader = TRUE, status = "primary",
            fluidRow(
              column(3,
                # --- Matrix 1 (Genetic) ---
                h5(tags$strong("Matrix 1 (genetic)")),
                radioButtons(
                  ns("mat1_source"),
                  "Source:",
                  choices = c(
                    "From Tab 1 (FST)"     = "tab1_fst",
                    "From Tab 1 (CS dist)" = "tab1_cs",
                    "Upload file"          = "upload1"
                  ),
                  selected = "tab1_fst"
                ),
                conditionalPanel(
                  condition = "input.mat1_source == 'tab1_fst'", ns = ns,
                  radioButtons(ns("mat1_fst_choice"), NULL,
                               choices = c("Uncorrected" = "fst", "ENA-corrected" = "fst_ena"),
                               selected = "fst_ena", inline = TRUE)
                ),
                conditionalPanel(
                  condition = "input.mat1_source == 'tab1_cs'", ns = ns,
                  radioButtons(ns("mat1_cs_choice"), NULL,
                               choices = c("Uncorrected" = "cs", "INA-corrected" = "cs_ina"),
                               selected = "cs_ina", inline = TRUE)
                ),
                conditionalPanel(
                  condition = "input.mat1_source == 'upload1'", ns = ns,
                  fileInput(ns("file_mat1"), "Distance file (CSV):",
                            accept = c(".csv", ".txt", ".tab")),
                  radioButtons(ns("mat1_format"), "Format:",
                               choices = c("Square matrix" = "square",
                                          "Rectangular (column-wise)" = "rectangular"),
                               selected = "square")
                ),

                tags$hr(),

                # --- Matrix 2 (Geographic) ---
                h5(tags$strong("Matrix 2 (geographic)")),
                radioButtons(
                  ns("mat2_source"),
                  "Source:",
                  choices = c(
                    "From GPS coordinates" = "gps",
                    "Upload file"          = "upload2"
                  ),
                  selected = "gps"
                ),
                conditionalPanel(
                  condition = "input.mat2_source == 'upload2'", ns = ns,
                  fileInput(ns("file_mat2"), "Distance file (CSV):",
                            accept = c(".csv", ".txt", ".tab")),
                  radioButtons(ns("mat2_format"), "Format:",
                               choices = c("Square matrix" = "square",
                                          "Rectangular (column-wise)" = "rectangular"),
                               selected = "square")
                ),
                conditionalPanel(
                  condition = "input.mat2_source == 'gps'", ns = ns,
                  checkboxInput(ns("use_log_dist"), "Use ln(distance)", value = TRUE)
                ),

                tags$hr(),
                numericInput(ns("n_perm_mantel"), "Permutations:",
                             value = 9999, min = 99, max = 99999, step = 1000),
                selectInput(ns("mantel_method"), "Correlation:",
                            choices = c("Pearson" = "pearson",
                                       "Spearman" = "spearman"),
                            selected = "pearson"),
                tags$hr(),
                actionButton(
                  ns("run_mantel"), "Run Mantel Test",
                  icon  = icon("play"),
                  class = "btn-action-primary btn-block",
                  style = "font-weight:bold;"
                )
              ),

              column(9,
                h4(icon("chart-bar"), "Mantel test results",
                   style = "font-weight:600; color:#2c3e50; margin-bottom:15px;"),
                fluidRow(
                  column(4, valueBoxOutput(ns("box_mantel_r"), width = NULL)),
                  column(4, valueBoxOutput(ns("box_mantel_p"), width = NULL)),
                  column(4, valueBoxOutput(ns("box_mantel_n"), width = NULL))
                ),
                verbatimTextOutput(ns("mantel_summary"))
              )
            )
          )
        ),

        # Mantel plot + Interpretation
        fluidRow(
          box(
            width = 8,
            title = div(
              style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
              icon("chart-scatter"), " Mantel scatter plot"
            ),
            solidHeader = FALSE,
            plotly::plotlyOutput(ns("mantel_plot"), height = "500px")
          ),
          box(
            width = 4,
            title = div(
              style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
              icon("info-circle"), " Interpretation"
            ),
            solidHeader = FALSE,
            tags$div(
              style = "font-size:13px; line-height:1.8;",
              tags$p(tags$strong("Mantel r:"), " Pearson/Spearman correlation between corresponding ",
                     "entries of the two distance matrices."),
              tags$p(tags$strong("P-value:"), " Proportion of permuted correlations \u2265 observed r (one-sided)."),
              tags$p(tags$strong("Significance:"), " p < 0.05 \u2192 significant correlation (e.g. IBD)."),
              tags$hr(),
              tags$p(tags$strong("File formats accepted:")),
              tags$ul(
                style = "font-size:12px; padding-left:16px; line-height:1.9;",
                tags$li(tags$em("Square matrix:"), " N\u00d7N symmetric matrix with population labels as row/col names"),
                tags$li(tags$em("Rectangular:"), " First row = pair labels (e.g. 'Pop1-Pop2'), columns = distances. ",
                        "Values are averaged over rows to obtain a single distance per pair.")
              ),
              tags$hr(),
              tags$p(style = "color:#777; font-size:12px;",
                "Mantel N. 1967. Math Geol 15:65-75.")
            )
          )
        )
      )
    )
  )
}