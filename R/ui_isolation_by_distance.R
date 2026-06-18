# ui_isolation_by_distance.R
# Tab: Isolation by Distance
# Two tabs: Pairwise Genetic Distances + Mantel Test

isolation_by_distance_UI <- function(id) {
  ns <- NS(id)

  fluidPage(
    tags$head(gs_head()),

    module_banner(
      "map-marker-alt",
      "Isolation by Distance",
      "Pairwise genetic distances · Bootstrap CI · Mantel test",
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

      # ── TAB 1: Pairwise Genetic Distances ─────────────────────────────────
      tabPanel(
        title = "Pairwise Distances",
        icon = icon("project-diagram"),
        tags$div(
          class = "spg-method-note", style = "border-left-color:#2CBF9F;",
          HTML(paste0(
            "Computes pairwise genetic distances (Cavalli-Sforza & Edwards 1967) ",
            "and F<sub>ST</sub> (Weir & Cockerham 1984) between all population pairs. ",
            "Bootstrap resampling over <b>loci</b> provides 95% confidence intervals. ",
            "Both uncorrected and ENA-corrected (Chapuis & Estoup 2007) estimates are computed. ",
            "<b>Requires genotype data imported with population assignments.</b>"
          ))
        ),

        # Configuration
        fluidRow(
          box(
            width = 12,
            title = div(
              style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
              icon("sliders-h"), " Pairwise distances parameters"
            ),
            solidHeader = TRUE, status = "primary",
            fluidRow(
              column(3,
                numericInput(ns("n_boot_dist"), "Bootstrap replicates (loci):",
                             value = 1000, min = 100, max = 10000, step = 100),
                checkboxInput(ns("calc_fst"), "Compute FST matrix", value = TRUE),
                checkboxInput(ns("calc_cs"), "Compute Cavalli-Sforza distance", value = TRUE),
                checkboxInput(ns("ena_correction"), "Apply ENA correction (null alleles)", value = TRUE),
                tags$hr(),
                actionButton(
                  ns("run_dist"), "Compute Pairwise Distances",
                  icon  = icon("calculator"),
                  class = "btn-action-primary btn-block",
                  style = "font-weight:bold;"
                )
              ),
              column(9,
                h4(icon("table"), "Results summary",
                   style = "font-weight:600; color:#2c3e50; margin-bottom:15px;"),
                fluidRow(
                  column(3, valueBoxOutput(ns("box_npops_dist"), width = NULL)),
                  column(3, valueBoxOutput(ns("box_nloci_dist"), width = NULL)),
                  column(3, valueBoxOutput(ns("box_npairs_dist"), width = NULL)),
                  column(3, valueBoxOutput(ns("box_boot_dist"), width = NULL))
                ),
                tags$h5("Pairwise distances matrix",
                        style = "font-weight:600; margin-top:14px; color:#2c3e50;"),
                DT::DTOutput(ns("dist_matrix_table"))
              )
            )
          )
        ),

        # Detailed pairwise table
        fluidRow(
          box(
            width = 12,
            title = div(
              style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
              icon("list"), " Detailed pairwise results"
            ),
            solidHeader = FALSE,
            DT::DTOutput(ns("pairwise_detail_table")),
            tags$br(),
            downloadButton(ns("dl_pairwise_csv"), "Download CSV",
                           class = "btn-action-secondary btn-sm")
          )
        )
      ),

      # ── TAB 2: Mantel Test ─────────────────────────────────────────────────
      tabPanel(
        title = "Mantel Test",
        icon = icon("chart-line"),
        tags$div(
          class = "spg-method-note", style = "border-left-color:#2CBF9F;",
          HTML(paste0(
            "Mantel test (Mantel 1967) assesses the correlation between two distance matrices ",
            "(e.g., genetic vs geographic distances). Significance is tested by permutation of rows/columns. ",
            "The test can be run on <b>rectangular format</b> (column-wise data as in Fstat/RT) ",
            "or on full distance matrices. ",
            "<b>Requires either computed pairwise distances + GPS coordinates, or uploaded distance files.</b>"
          ))
        ),

        # Configuration
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
                radioButtons(
                  ns("mantel_data_source"),
                  "Data source:",
                  choices = c(
                    "Use computed pairwise distances" = "computed",
                    "Upload distance files" = "upload"
                  ),
                  selected = "computed"
                ),
                conditionalPanel(
                  condition = "input.mantel_data_source == 'upload'",
                  ns = ns,
                  fileInput(ns("file_gen_dist"), "Genetic distance file (CSV):",
                            accept = c(".csv", ".txt")),
                  fileInput(ns("file_geo_dist"), "Geographic distance file (CSV):",
                            accept = c(".csv", ".txt"))
                ),
                numericInput(ns("n_perm_mantel"), "Permutations:",
                             value = 9999, min = 99, max = 99999, step = 1000),
                selectInput(ns("mantel_method"), "Correlation method:",
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
                tags$h5("Test details",
                        style = "font-weight:600; margin-top:14px; color:#2c3e50;"),
                verbatimTextOutput(ns("mantel_summary"))
              )
            )
          )
        ),

        # Mantel plot
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
              tags$p(tags$strong("Mantel statistic (r):")),
              tags$ul(
                style = "font-size:12px; padding-left:16px; line-height:1.9;",
                tags$li("Correlation between two distance matrices"),
                tags$li("Ranges from -1 to 1"),
                tags$li("Positive r: distances increase together")
              ),
              tags$p(tags$strong("P-value:"), " Probability of observing such correlation by chance (permutation test)"),
              tags$p(tags$strong("Significance:"), " p < 0.05 indicates significant isolation by distance"),
              tags$hr(),
              tags$p(style = "color:#777; font-size:12px;",
                "Mantel N. 1967. Math. Geol. 15:65-75.")
            )
          )
        )
      )
    )
  )
}