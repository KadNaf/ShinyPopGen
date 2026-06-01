# module/ui_null_alleles.R

null_alleles_UI <- function(id) {
  ns <- NS(id)

  tagList(
    module_banner("circle-notch", "Null Alleles",
      "Frequency estimation by locus \u00d7 population \u00b7 FreeNA EM algorithm (Chapuis & Estoup 2007)",
      "#8D8680"),

    # ── Method banner ──────────────────────────────────────────────────────
    fluidRow(
      column(12,
        div(style = "border-left: 4px solid #8D8680; background:#D9D0D3; border-radius:4px; padding:12px 16px; margin-bottom: 18px;",
            icon("info-circle"), strong(" Method: "),
            "EM algorithm (Dempster, Laird & Rubin 1977) as implemented in ",
            strong("FreeNA"), " (Chapuis & Estoup 2007).",
            br(),
            "Frequencies are estimated ", strong("per locus and per sub-population"),
            ", then a weighted mean is calculated:",
            tags$em(HTML(" p̄<sub>null</sub> = &Sigma;(p<sub>null-i</sub> &times; N<sub>i</sub>
              &times; Ĥ<sub>S-i</sub>) / &Sigma;(N<sub>i</sub> &times; Ĥ<sub>S-i</sub>)")),
            "where N", tags$sub("i"), " = sample size of sub-population i and Ĥ",
            tags$sub("S-i"), " = its genetic diversity at the considered locus."
        )
      )
    ),

    fluidRow(
      # ── Left panel: configuration ─────────────────────────────────────
      box(
        width = 4,
        title = div(style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
                    icon("cogs"), " Null alleles parameters"),
        solidHeader = TRUE, status = "primary",

        h4(icon("code"), " Missing Data Coding"),
        div(style = "background:#D9D0D3; border-left:4px solid #8D8680; border-radius:4px; font-size:0.88em; padding:8px;",
            icon("exclamation-triangle"),
            HTML("<br><b>FreeNA convention:</b><br>
                 &bull; <code>999999/999999</code> &rarr; <b>null homozygote</b>
                 (counted in EM algorithm)<br>
                 &bull; <code>0/0</code> &rarr; <b>missing genotype</b>
                 (excluded from analysis)<br>
                 &bull; Null <em>heterozygotes</em> are <b>forbidden</b> in FreeNA<br><br>
                 If your data uses different codes, modify them below.")
        ),

        br(),
        h5(icon("sliders-h"), " Default recoding (all loci)"),
        selectInput(
          ns("default_missing_recode"),
          label   = "Treat missing values as:",
          choices = c("0/0 — missing genotype (excluded)" = "0",
                      "999999/999999 — null homozygote (included)" = "999999"),
          selected = "0"
        ),

        hr(),
        h5(icon("list"), " Per-locus overrides"),
        p(style = "font-size:0.88em; color:#555;",
          'Choose "Use default" to keep the global choice.',
          "Otherwise, select the specific recoding for this locus."),
        uiOutput(ns("locus_coding_ui")),

        hr(),
        actionButton(ns("run_null_alleles"),
                     label = "Estimate null allele frequencies",
                     icon  = icon("calculator"),
                     class = "btn-action-primary btn-block",
                     style = "font-weight:bold;")
      ),

      # ── Right panel: results ───────────────────────────────────────────
      box(
        width = 8,
        title = div(style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
                    icon("table"), " Results"),
        solidHeader = TRUE, status = "primary",

        # Tabs for navigation
        tabsetPanel(
          id = ns("results_tabs"),
          type = "tabs",
          
          # Tab 1: Detailed results by locus × sub-population
          tabPanel(
            title = tagList(icon("table"), "Detail by sub-population"),
            br(),
            p(style = "font-size:0.88em; color:#555;",
              "Each row corresponds to a (locus, sub-population) pair.",
              "Columns indicate:"),
            tags$ul(style = "font-size:0.85em; color:#555;",
              tags$li(tags$strong("p̂_null"), " : estimated null allele frequency"),
              tags$li(tags$strong("Ĥe"), " : observed genetic diversity (He)"),
              tags$li(tags$strong("N valid"), " : individuals with non-missing genotype"),
              tags$li(tags$strong("N null hom."), " : null homozygote individuals (999999/999999)"),
              tags$li(tags$strong("Iterations"), " : number of EM iterations"),
              tags$li(tags$strong("Converged"), " : ✓ if algorithm converged")
            ),
            br(),
            DT::DTOutput(ns("null_allele_detail_table"))
          ),
          
          # Tab 2: Weighted means (summary)
          tabPanel(
            title = tagList(icon("chart-line"), "Weighted means"),
            br(),
            p(style = "font-size:0.88em; color:#555;",
              "Inter-population weighted means by N", tags$sub("i"), " × Ĥ", tags$sub("S-i"),
              " as in FreeNA."),
            DT::DTOutput(ns("null_allele_summary_table"))
          ),
          
          # Tab 3: Visualization
          tabPanel(
            title = tagList(icon("chart-bar"), "Visualization"),
            br(),
            fluidRow(
              column(6,
                plotOutput(ns("null_allele_dist_plot"), height = "400px")
              ),
              column(6,
                plotOutput(ns("null_allele_heatmap"), height = "400px")
              )
            ),
            hr(),
            fluidRow(
              column(12,
                plotOutput(ns("null_allele_convergence_plot"), height = "300px")
              )
            ),
            br(),
            downloadButton(ns("download_null_alleles_png"), ".png",
                           class = "btn-download-primary")
          )
        ),

        br(),
        fluidRow(
          column(4, downloadButton(ns("download_null_alleles_detail"),  "Detail (.csv)",  class = "btn-download-primary btn-block")),
          column(4, downloadButton(ns("download_null_alleles_summary"), "Summary (.csv)", class = "btn-download-primary btn-block")),
          column(4, downloadButton(ns("download_null_alleles_txt"),     "Report (.txt)",  class = "btn-download-secondary btn-block"))
        ),

        hr(),
        h5(icon("align-left"), " Statistical report"),
        verbatimTextOutput(ns("null_allele_summary"))
      )
    )
  )
}