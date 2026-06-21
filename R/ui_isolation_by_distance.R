# ui_null_alleles.R
# Null allele frequency estimation (EM) + Fst/Fst-ENA (Weir 1996/Genepop) +
# DCSE/DCSE-INA (Cavalli-Sforza & Edwards 1967), faithfully translated from
# the Pascal reference FreeNA_optm2R.pas (engine_freena.R).
#
# Bootstrap: OVER LOCI ONLY (as in the Pascal source) — one shared resampled
# loci sequence per replicate, applied simultaneously to every statistic.
# Minimum requirements (faithful to Pascal): nperm >= 100, nloc > 4.
#
# Tabs 5-6 (Isolation by Distance, Mantel test) consume the pairwise
# FST/FST-ENA/DCSE/DCSE-INA (+ bootstrap CI) computed in tabs 1-4 — nothing
# is recomputed, following the SPG-V1 module specification:
#   - Rousset (1997) regression: FR ~ Dgeo (model 1) or FR ~ ln(Dgeo) (model 2),
#     FR = FST/(1-FST) or FST-ENA/(1-FST-ENA); 3 lines fitted (point estimate,
#     lower CI, upper CI) — IBD signal if all 3 slopes are positive.
#   - Generic Mantel test (RT / Fstat 2.9.4 convention): Pearson r or Rousset
#     slope, joint row/column permutation (valid on rectangular matrices),
#     p = (b+1)/(m+1), one-sided positive/negative, % variance explained.

isolation_by_distance_UI <- function(id) {
  ns <- NS(id)

  fluidPage(
    tags$head(gs_head()),

    module_banner(
      "atom",
      "Null Allele Estimation \u00b7 Fst-ENA \u00b7 DCSE-INA \u00b7 Isolation by Distance \u00b7 Mantel",
      "EM algorithm (Dempster, Laird & Rubin 1977) \u00b7 FreeNA (Chapuis & Estoup 2007) \u00b7 Weir (1996) \u00b7 Cavalli-Sforza & Edwards (1967) \u00b7 Rousset (1997)",
      "#0c4a6e"
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
          "Engine validated against an independent WC84 implementation and unit tests; ",
          "UI/DB integration still being checked end-to-end.")
      )
    ),

    tags$div(
      class = "spg-method-note", style = "border-left-color:#0c4a6e;",
      HTML(paste0(
        "Faithful R translation of the Pascal reference program ",
        "<b>FreeNA_optm2R.pas</b>. Null allele frequencies (r<sub>d</sub>) are estimated ",
        "per locus &times; population with the EM algorithm; F<sub>ST</sub> (raw, Weir 1996 / ",
        "Genepop method) and F<sub>ST</sub>-ENA (null-allele corrected, Chapuis &amp; Estoup 2007) ",
        "use the same <b>double n<sub>c</sub>-weighting</b> scheme as the Pascal source when ",
        "combining loci into a multilocus estimate. D<sub>CSE</sub> (Cavalli-Sforza &amp; Edwards ",
        "1967 chord distance) and D<sub>CSE</sub>-INA (null allele appended as an extra, ",
        "non-renormalised category) follow the same logic. ",
        "<b>95% CIs are obtained by bootstrap over loci only</b> (one shared resampled-loci ",
        "sequence per replicate, applied to every statistic simultaneously) \u2014 exactly as in ",
        "the Pascal source. Bootstrap requires at least 100 replicates and more than 4 loci."
      ))
    ),

    fluidRow(
      box(
        width = 3,
        title = div(style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
                    icon("sliders-h"), " Parameters"),
        solidHeader = TRUE, status = "primary",

        textInput(ns("null_code"), "Null allele code:", value = "999999",
                  placeholder = "e.g. 999999 (3-digit) or 99 (2-digit)"),
        tags$p(style = "color:#777; font-size:11px;",
          "Genepop convention: missing = \"0\"/\"00\"/\"000\"/\"0000\"/\"000000\" ",
          "(handled automatically); null homozygote = the code above, repeated for ",
          "both alleles (e.g. \"999999/999999\")."),

        tags$hr(),

        numericInput(ns("n_boot"), "Bootstrap replicates (loci):",
                     value = 2000, min = 0, max = 20000, step = 100),
        sliderInput(ns("conf_level"), "Confidence level (%):",
                    min = 80, max = 99, value = 95, step = 1),
        tags$p(style = "color:#777; font-size:11px;",
          "Set replicates to 0 to skip bootstrapping. Per Pascal source: bootstrap is ",
          "skipped if replicates < 100 or number of loci \u2264 4."),

        tags$hr(),
        actionButton(
          ns("run_all"), "Compute Everything",
          icon  = icon("play"),
          class = "btn-action-primary btn-block",
          style = "font-weight:bold;"
        ),
        uiOutput(ns("ui_run_status"))
      ),

      box(
        width = 9,
        title = div(style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
                    icon("chart-bar"), " Summary"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3, valueBoxOutput(ns("box_nloci"),  width = NULL)),
          column(3, valueBoxOutput(ns("box_npops"),  width = NULL)),
          column(3, valueBoxOutput(ns("box_avgrd"),  width = NULL)),
          column(3, valueBoxOutput(ns("box_fstena"), width = NULL))
        )
      )
    ),

    tabsetPanel(
      id = ns("na_tabs"), type = "tabs",

      # ── TAB 1: Null allele frequencies ──────────────────────────────────
      tabPanel(title = tagList(icon("dna"), " Null allele frequencies"), value = "tab_rd",
        br(),
        box(width = 12, solidHeader = TRUE, status = "primary",
            title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                        icon("list"), " r", tags$sub("d"), " per locus \u00d7 population (EM algorithm)"),
          DT::DTOutput(ns("dt_rd")),
          tags$br(),
          downloadButton(ns("dl_rd_csv"), "Download CSV", class = "btn-action-secondary btn-sm")
        )
      ),

      # ── TAB 2: Fst / Fst-ENA ─────────────────────────────────────────────
      tabPanel(title = tagList(icon("chart-bar"), " Fst / Fst-ENA"), value = "tab_fst",
        br(),
        box(width = 12, solidHeader = TRUE, status = "primary",
            title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                        icon("globe"), " Global multilocus F", tags$sub("ST")),
          uiOutput(ns("ui_global_fst"))
        ),
        box(width = 12, solidHeader = FALSE,
            title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                        icon("th"), " Pairwise F", tags$sub("ST"), " (lower triangle)"),
          radioButtons(ns("fst_display"), "Display:",
            choices = c("Raw" = "raw", "ENA-corrected" = "ena", "Both" = "both"),
            selected = "both", inline = TRUE),
          uiOutput(ns("ui_fst_matrix")),
          tags$br(),
          DT::DTOutput(ns("dt_fst_pair_ci")),
          downloadButton(ns("dl_fst_pair_csv"), "Download pairwise Fst (CSV)",
                         class = "btn-action-secondary btn-sm")
        )
      ),

      # ── TAB 3: DCSE / DCSE-INA ───────────────────────────────────────────
      tabPanel(title = tagList(icon("ruler-combined"), " DCSE / DCSE-INA"), value = "tab_dc",
        br(),
        box(width = 12, solidHeader = TRUE, status = "primary",
            title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                        icon("th"), " Pairwise D", tags$sub("CSE"), " (lower triangle)"),
          radioButtons(ns("dc_display"), "Display:",
            choices = c("Raw" = "raw", "INA-corrected" = "ina", "Both" = "both"),
            selected = "both", inline = TRUE),
          uiOutput(ns("ui_dc_matrix")),
          tags$br(),
          DT::DTOutput(ns("dt_dc_pair_ci")),
          downloadButton(ns("dl_dc_pair_csv"), "Download pairwise DCSE (CSV)",
                         class = "btn-action-secondary btn-sm")
        )
      ),

      # ── TAB 4: Per-locus x pair ──────────────────────────────────────────
      tabPanel(title = tagList(icon("table"), " Per-locus \u00d7 pair"), value = "tab_lp",
        br(),
        box(width = 12, solidHeader = TRUE, status = "primary",
            title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                        icon("list"), " Per-locus breakdown (outlier detection)"),
          fluidRow(
            column(4, selectInput(ns("fl_locus"), "Locus:", choices = c("All loci" = "all"))),
            column(4, selectInput(ns("fl_pop1"), "Population 1:", choices = c("All pairs" = "all"))),
            column(4, selectInput(ns("fl_pop2"), "Population 2:", choices = c("All pairs" = "all")))
          ),
          DT::DTOutput(ns("dt_locus_pair")),
          downloadButton(ns("dl_locus_pair_csv"), "Download CSV", class = "btn-action-secondary btn-sm")
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # TAB 5: Isolation by Distance (Rousset 1997)
      # ══════════════════════════════════════════════════════════════════
      tabPanel(title = tagList(icon("map-marker-alt"), " Isolation by Distance"), value = "tab_ibd",
        br(),

        tags$div(
          class = "spg-method-note", style = "border-left-color:#2CBF9F;",
          HTML(paste0(
            "Rousset (1997) regression: model 1 (1D) F<sub>R</sub> \u223c D<sub>geo</sub>, or ",
            "model 2 (2D) F<sub>R</sub> \u223c ln(D<sub>geo</sub>), where F<sub>R</sub> = F<sub>ST</sub>/(1\u2212F<sub>ST</sub>) ",
            "or F<sub>ST</sub>-ENA/(1\u2212F<sub>ST</sub>-ENA) \u2014 computed here directly from the F<sub>ST</sub> / ",
            "F<sub>ST</sub>-ENA values and their bootstrap CIs already calculated in the previous tabs ",
            "(nothing is recomputed). Three regression lines are fitted: the point estimate (F<sub>R</sub>), ",
            "and its lower/upper confidence bounds (F<sub>R</sub>-l, F<sub>R</sub>-u). ",
            "<b>If all three slopes are positive, this supports isolation by distance.</b> ",
            "If the lower-bound slope is negative while the others are positive, this may indicate ",
            "low power of the per-pair bootstrap rather than a true absence of IBD \u2014 a Mantel test ",
            "(next tab), ideally using D<sub>CSE</sub>, can help confirm the conclusion (S\u00e9r\u00e9 et al. 2017)."
          ))
        ),

        fluidRow(
          box(width = 3, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("sliders-h"), " Parameters"),
            radioButtons(ns("ibd_model"), "Model:",
              choices = c("Model 1 (1D): FR ~ Dgeo"        = "1D",
                          "Model 2 (2D): FR ~ ln(Dgeo)"     = "2D"),
              selected = "2D"),
            radioButtons(ns("ibd_metric"), "Genetic distance:",
              choices = c("FR (raw Fst)"      = "raw",
                          "FR-ENA (corrected)" = "ena"),
              selected = "ena"),
            tags$hr(),
            tags$p(style="color:#777;font-size:11px;",
              "Requires GPS (Latitude/Longitude) set at import for at least 2 populations. ",
              "Population centroid is the mean GPS of its individuals; distance is the ",
              "great-circle (Haversine) distance."),
            actionButton(ns("run_ibd"), "Run IBD Regression",
                         icon = icon("rocket"), class = "btn-action-primary btn-block",
                         style = "font-weight:bold;")
          ),
          box(width = 9, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("chart-line"), " Regression results"),
            DT::DTOutput(ns("dt_ibd_reg")),
            tags$br(),
            uiOutput(ns("ui_ibd_interpretation"))
          )
        ),

        fluidRow(
          box(width = 7, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("chart-area"), " IBD plot"),
            plotly::plotlyOutput(ns("ibd_plot"), height = "440px")
          ),
          box(width = 5, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("table"), " Pair table used"),
            DT::DTOutput(ns("dt_ibd_table")),
            tags$br(),
            downloadButton(ns("dl_ibd_csv"), "Download CSV", class = "btn-action-secondary btn-sm")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # TAB 6: Mantel test
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
                "Internal pairwise table (this module)" = "internal",
                "Upload external column file"            = "upload"
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