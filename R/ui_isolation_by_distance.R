# ui_isolation_by_distance.R
# Isolation by Distance (Rousset 1997) + generic Mantel test.
#
# This module is a CONTINUATION of the "Null alleles" module — it reuses the
# pairwise FST/FST-ENA/DCSE/DCSE-INA (+ bootstrap CI) already computed there
# (rv$null_alleles_results) instead of recomputing them, so the tabs that
# used to duplicate that module (Null allele frequencies, Fst/Fst-ENA,
# DCSE/DCSE-INA, Per-locus x pair) have been removed from here.
#
# Workflow: go to the "Null alleles" module first, choose your per-locus
# coding, click "Compute + Bootstrap + Export" — THEN come here.

isolation_by_distance_UI <- function(id) {
  ns <- NS(id)

  fluidPage(
    tags$head(gs_head()),

    module_banner(
      "atom",
      "Isolation by Distance \u00b7 Mantel Test",
      "Rousset (1997) IBD regression \u00b7 Mantel (1967) permutation test \u2014 built on the pairwise FST-ENA / DCSE-INA already computed in the Null Alleles module",
      "#0c4a6e"
    ),

    tags$div(
      class = "spg-method-note", style = "border-left-color:#0c4a6e;",
      HTML(paste0(
        "This module reuses the pairwise F<sub>ST</sub> / F<sub>ST</sub>-ENA / D<sub>CSE</sub> / ",
        "D<sub>CSE</sub>-INA (+ bootstrap CI) already computed in the <b>Null Alleles</b> module — ",
        "nothing is recomputed here. Go there first, choose your per-locus coding, and click ",
        "<b>\"Compute + Bootstrap + Export\"</b>, then come back to this module.",
        "<br>Geographic distance (D<sub>geo</sub>) is the <b>Vincenty ellipsoidal geodesic distance</b> ",
        "(WGS84), in metres, from each population's GPS centroid (mean Latitude/Longitude of its individuals) ",
        "\u2014 or, alternatively, distances loaded from an <b>external Pop1/Pop2/Distance file</b> ",
        "(e.g. the subsample-pairs template exported by the <b>Subdivision</b> module and edited by hand)."
      ))
    ),

    fluidRow(
      box(
        width = 12, solidHeader = TRUE, status = "primary",
        title = div(style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
                    icon("chart-bar"), " Summary (from Null Alleles module)"),
        fluidRow(
          column(3, valueBoxOutput(ns("box_nloci"),  width = NULL)),
          column(3, valueBoxOutput(ns("box_npops"),  width = NULL)),
          column(3, valueBoxOutput(ns("box_fstena"), width = NULL)),
          column(3, valueBoxOutput(ns("box_nboot"),  width = NULL))
        ),
        uiOutput(ns("ui_run_status"))
      )
    ),

    tabsetPanel(
      id = ns("ibd_tabs"), type = "tabs",

      # ══════════════════════════════════════════════════════════════════
      # TAB 1 — Isolation by Distance (now the FIRST tab)
      # ══════════════════════════════════════════════════════════════════
      tabPanel(title = tagList(icon("chart-line"), " Isolation by Distance"), value = "tab_ibd",
        br(),
        fluidRow(
          box(width = 3, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("sliders-h"), " Parameters"),
            radioButtons(ns("ibd_source"), "Data source:",
              choices = c("Null Alleles module (this session)" = "internal",
                          "Re-load an exported pairwise file"   = "external"),
              selected = "internal"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'external'", ns("ibd_source")),
              tags$p(style = "color:#777;font-size:11px;",
                "Load a pairwise file previously exported by this app (e.g. the Null Alleles ",
                "module's \u201cpairwise_long_format\u201d file) \u2014 or that same file freely edited by ",
                "hand (rows removed, values corrected). This lets the IBD module run entirely on ",
                "its own, without needing the Null Alleles module in this session."),
              fileInput(ns("ibd_ext_file"), "Pairwise file (Pop1, Pop2, FST_ENA, DCSE_INA, Dgeo_m\u2026):",
                        accept = c(".csv", ".txt", ".tsv")),
              radioButtons(ns("ibd_ext_sep"), "Separator:",
                choices = c("Tab" = "\t", "Comma" = ",", "Semicolon" = ";"),
                selected = "\t", inline = TRUE),
              checkboxInput(ns("ibd_ext_header"), "File has header row", value = TRUE)
            ),
            tags$hr(),
            radioButtons(ns("ibd_model"), "Habitat model:",
              choices = c("2D (F_R ~ ln(D_geo))" = "2D",
                          "1D (F_R ~ D_geo)"      = "1D"),
              selected = "2D"),
            radioButtons(ns("ibd_metric"), "Genetic distance metric:",
              choices = c("F_R (raw FST)"      = "raw",
                          "F_R (FST-ENA)"      = "ena"),
              selected = "ena"),
            tags$hr(),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'internal'", ns("ibd_source")),
              radioButtons(ns("ibd_dgeo_source"), "Distance (D_geo) source:",
                choices = c("GPS centroid (auto, Vincenty)"       = "gps",
                            "External pairs/distances file"       = "external"),
                selected = "gps"),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'external'", ns("ibd_dgeo_source")),
                tags$p(style = "color:#777;font-size:11px;",
                  "Load a file with one row per subsample pair (columns: Pop1, Pop2, Distance). ",
                  "Use the template generated by the ", tags$b("Subdivision"), " module: open it in a ",
                  "spreadsheet or text editor, delete any pairs you don't want, and fill in the ",
                  "Distance column with your own values \u2014 geographic or otherwise."),
                fileInput(ns("ibd_dgeo_file"), "Pairs/distances file (Pop1, Pop2, Distance):",
                          accept = c(".csv", ".txt", ".tsv")),
                radioButtons(ns("ibd_dgeo_sep"), "Separator:",
                  choices = c("Comma" = ",", "Tab" = "\t", "Semicolon" = ";"),
                  selected = ",", inline = TRUE),
                checkboxInput(ns("ibd_dgeo_header"), "File has header row", value = TRUE),
                tags$p(style = "color:#777;font-size:11px;",
                  "Only pairs present in the file are used; pairs you deleted from the file are excluded from the analysis.")
              ),
              tags$p(style="color:#777;font-size:11px;",
                "GPS mode requires Latitude/Longitude set at import for at least 2 populations. ",
                "Population centroid is the mean GPS of its individuals; distance is the ",
                "Vincenty ellipsoidal geodesic distance (WGS84), in metres.")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'external'", ns("ibd_source")),
              tags$p(style = "color:#777;font-size:11px;",
                "D_geo is read directly from the uploaded file's Dgeo_m / lnDgeo columns when present.")
            ),
            tags$hr(),
            tags$div(style="font-size:12px;color:#555;margin-bottom:4px;", "Output file (regression + full table):"),
            fluidRow(
              column(7, textInput(ns("ibd_out_root"), "Root:", value = "",
                                   placeholder = "auto-filled from imported file")),
              column(5, textInput(ns("ibd_out_suffix"), "Suffix:", value = ""))
            ),
            tags$p(style="color:#777;font-size:11px;",
              "File name = ", tags$code("<root>-IBD-<suffix>.txt"), " \u2014 edit either field freely."),
            actionButton(ns("run_ibd"), "Run IBD Regression",
                         icon = icon("rocket"), class = "btn-action-primary btn-block",
                         style = "font-weight:bold;")
          ),
          box(width = 9, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("chart-line"), " Regression summary (slope / b / Nb / Nem)"),
            DT::DTOutput(ns("dt_ibd_reg")),
            tags$br(),
            uiOutput(ns("ui_ibd_interpretation"))
          )
        ),

        fluidRow(
          box(width = 12, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("table"), " Full pairwise table (Farm1, Farm2, D_geo, FST-FreeNA, F_R, D_CSE-INA, D_CSE)"),
            DT::DTOutput(ns("dt_ibd_table")),
            tags$br(),
            downloadButton(ns("dl_ibd_txt"), "Download full table + regression summary (.txt)", class = "btn-action-secondary btn-sm")
          )
        ),

        fluidRow(
          box(width = 12, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("chart-area"), " IBD plot"),
            plotly::plotlyOutput(ns("ibd_plot"), height = "440px")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # TAB 2 — Mantel test
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
            "valid when either matrix is incomplete. Statistic: Pearson's r or Spearman's rho ",
            "(Fstat convention) or the slope of the Rousset (1997) regression (Genepop convention ",
            "for IBD). ",
            "One-sided p-value = (b+1)/(m+1), b = number of permuted statistics \u2265 observed."
          ))
        ),

        fluidRow(
          box(width = 4, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("database"), " Data source"),
            radioButtons(ns("mt_source"), NULL,
              choices = c(
                "Internal pairwise table (from Null Alleles module)" = "internal",
                "Upload external column file"                        = "upload"
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
                  choices = c("Pearson r" = "r", "Spearman rho" = "spearman",
                              "Regression slope (Rousset)" = "b"),
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
      ),

      # ══════════════════════════════════════════════════════════════════
      # TAB 3 — Partial Mantel test (multiple matrices)
      # ══════════════════════════════════════════════════════════════════
      tabPanel(title = tagList(icon("layer-group"), " Partial Mantel"), value = "tab_partial_mantel",
        br(),

        tags$div(
          class = "spg-method-note", style = "border-left-color:#B45309;",
          HTML(paste0(
            "Partial Mantel test for <b>more than 2\u20133 matrices at once</b> (up to 10), via ",
            "multiple regression on distance matrices (MRM; Legendre, Lapointe &amp; Casgrain 1994) ",
            "\u2014 the standard generalisation of the classic partial Mantel test. Base R packages ",
            "(ade4, vegan, ecodist) cap partial Mantel at 2\u20133 matrices and Pearson/Spearman/Kendall ",
            "only; this tab removes both limits and works on the <b>same rectangular/incomplete data</b> ",
            "as the Mantel Test tab (uses the same data source and Pop1/Pop2 columns \u2014 set those there first).",
            "<br><b>Caveat:</b> like the classic partial Mantel test, this MRM generalisation can have ",
            "inflated type I error when the matrices being partialled out (e.g. geographic distance) are ",
            "themselves spatially autocorrelated (Guillot &amp; Rousset 2013; Crabot et al. 2019, ",
            "<i>Methods Ecol Evol</i> 10:532\u2013540). A plain Mantel test per predictor (previous tab) ",
            "and, in the future, the Procrustes association metric (Lisboa et al. 2014, ",
            "<i>PLoS ONE</i> 9(6):e101238) as a less-controversial alternative, are worth cross-checking ",
            "any result against. Borcard &amp; Legendre (2012, <i>Ecology</i> 93:1473\u20131481) found the ",
            "plain Mantel test itself has acceptable power for most ecological uses."
          ))
        ),

        fluidRow(
          box(width = 4, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("database"), " Variables (from the Mantel Test tab's data source)"),
            uiOutput(ns("pm_col_y_ui")),
            uiOutput(ns("pm_col_x_ui")),
            tags$hr(),
            checkboxInput(ns("pm_standardize"), "Standardize (z-score) all variables before fitting",
                          value = TRUE),
            numericInput(ns("pm_n_perm"), "Permutations:", value = 999, min = 99, max = 20000, step = 100),
            tags$p(style="color:#777;font-size:11px;",
              "Permutation is by joint row/column relabelling of the response (Y) matrix, keeping all ",
              "predictor matrices fixed \u2014 the standard MRM permutation scheme, valid on incomplete data."),
            actionButton(ns("run_partial_mantel"), "Run Partial Mantel (MRM)",
                         icon = icon("random"), class = "btn-action-primary btn-block",
                         style = "font-weight:bold;")
          ),
          box(width = 8, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("sliders-h"), " Results"),
            DT::DTOutput(ns("dt_partial_mantel")),
            uiOutput(ns("ui_partial_mantel_summary")),
            tags$br(),
            downloadButton(ns("dl_partial_mantel_csv"), "Download results (CSV)", class = "btn-action-secondary btn-sm")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # TAB 4 — Procrustes analysis & PAM (Lisboa et al. 2014)
      # ══════════════════════════════════════════════════════════════════
      tabPanel(title = tagList(icon("shapes"), " Procrustes / PAM"), value = "tab_pam",
        br(),

        tags$div(
          class = "spg-method-note", style = "border-left-color:#B45309;",
          HTML(paste0(
            "Procrustes analysis and the <b>Procrustean Association Metric (PAM)</b> \u2014 ",
            "Lisboa, Peres-Neto, Chaer, Jesus, Mitchell, Chapman &amp; Berbara (2014, <i>PLoS ONE</i> ",
            "9(6):e101238) \u2014 a less-controversial, more powerful alternative to (partial) Mantel tests. ",
            "Instead of comparing pairwise distance matrices directly, each of the two chosen matrices ",
            "is first ordinated (PCoA) into a configuration of populations in k-dimensional space; the two ",
            "configurations are then superimposed (Procrustes rotation) and compared with a global ",
            "significance test (<b>PROTEST</b>) <i>and</i> a <b>per-population residual (PAM)</b> \u2014 ",
            "showing which specific populations drive an overall (mis)match, something the Mantel tests ",
            "above cannot do. The PAM vector can itself be used as a response variable in a follow-up ",
            "regression or ANOVA (e.g. against herd management, host breed, or other covariates), as ",
            "illustrated in the paper.",
            "<br><b>Note:</b> unlike the Mantel tests, Procrustes/PAM needs a <b>complete</b> square ",
            "matrix \u2014 populations with any missing pairwise value are dropped automatically (reported ",
            "below), rather than silently tolerated as in the rectangular-safe Mantel tests."
          ))
        ),

        fluidRow(
          box(width = 4, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("database"), " Matrices (from the Mantel Test tab's data source)"),
            uiOutput(ns("pam_col_x_ui")),
            uiOutput(ns("pam_col_y_ui")),
            tags$hr(),
            numericInput(ns("pam_k"), "Number of PCoA axes (k):", value = 2, min = 1, max = 10, step = 1),
            numericInput(ns("pam_n_perm"), "PROTEST permutations:", value = 999, min = 99, max = 20000, step = 100),
            tags$p(style="color:#777;font-size:11px;",
              "Both configurations are ordinated with the same number of axes (Lisboa et al. 2014, Fig. 2b) ",
              "and rescaled to unit sum of squares before fitting, so m\u00b2 and r are comparable across pairs ",
              "of variables."),
            actionButton(ns("run_pam"), "Run Procrustes + PROTEST",
                         icon = icon("rocket"), class = "btn-action-primary btn-block",
                         style = "font-weight:bold;")
          ),
          box(width = 8, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("chart-bar"), " Global fit (PROTEST)"),
            fluidRow(
              column(4, valueBoxOutput(ns("box_pam_r"),  width = NULL)),
              column(4, valueBoxOutput(ns("box_pam_m2"), width = NULL)),
              column(4, valueBoxOutput(ns("box_pam_p"),  width = NULL))
            ),
            valueBoxOutput(ns("box_pam_n"), width = 12)
          )
        ),

        fluidRow(
          box(width = 6, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("table"), " PAM per population"),
            DT::DTOutput(ns("dt_pam")),
            uiOutput(ns("ui_pam_interpretation")),
            tags$br(),
            downloadButton(ns("dl_pam_csv"), "Download PAM (CSV)", class = "btn-action-secondary btn-sm")
          ),
          box(width = 6, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("chart-bar"), " PAM per population (bar chart)"),
            plotly::plotlyOutput(ns("pam_plot"), height = "420px")
          )
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # TAB 5 — Correlogram (Borcard & Legendre 2012) & MSR-Mantel (Crabot
      # et al. 2019)
      # ══════════════════════════════════════════════════════════════════
      tabPanel(title = tagList(icon("wave-square"), " Correlogram / MSR-Mantel"), value = "tab_correlogram",
        br(),

        tags$div(
          class = "spg-method-note", style = "border-left-color:#0c4a6e;",
          HTML(paste0(
            "Two complementary approaches to spatial autocorrelation, both reusing the Mantel Test ",
            "tab's data source and Pop1/Pop2 columns \u2014 set those there first.<br>",
            "<b>Mantel correlogram</b> (Borcard &amp; Legendre 2014, <i>Ecology</i> 93:1473\u20131481; ",
            "Sokal 1986; Oden &amp; Sokal 1986): splits distance into classes and tests, class by class, ",
            "whether the response matrix is more similar within than among that class \u2014 shows ",
            "<b>at what spatial scale</b> the pattern occurs, which a single Mantel r cannot. Borcard &amp; ",
            "Legendre found it has acceptable power for most ecological uses.<br>",
            "<b>MSR-Mantel</b> (Wagner &amp; Dray 2015; Crabot, Clappe, Dray &amp; Datry 2019, ",
            "<i>Methods Ecol Evol</i> 10:532\u2013540): corrects the classic Mantel test's inflated type I ",
            "error when <b>both</b> matrices are spatially autocorrelated, by comparing the observed ",
            "statistic to a null distribution built from spatially-constrained random replicates (rather ",
            "than fully random permutations) of one matrix."
          ))
        ),

        fluidRow(
          box(width = 12, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("chart-line"), " Mantel correlogram"),
            fluidRow(
              column(3, uiOutput(ns("corr_col_y_ui"))),
              column(3, uiOutput(ns("corr_col_d_ui"))),
              column(2, numericInput(ns("corr_n_classes"), "Distance classes:", value = 7, min = 2, max = 30, step = 1)),
              column(2, numericInput(ns("corr_n_perm"), "Permutations:", value = 999, min = 99, max = 20000, step = 100)),
              column(2, numericInput(ns("corr_alpha"), "Alpha:", value = 0.05, min = 0.001, max = 0.5, step = 0.01))
            ),
            actionButton(ns("run_correlogram"), "Run Mantel Correlogram",
                         icon = icon("random"), class = "btn-action-primary",
                         style = "font-weight:bold;"),
            tags$hr(),
            fluidRow(
              column(6, DT::DTOutput(ns("dt_correlogram"))),
              column(6, plotly::plotlyOutput(ns("correlogram_plot"), height = "380px"))
            ),
            tags$br(),
            downloadButton(ns("dl_correlogram_csv"), "Download correlogram (CSV)", class = "btn-action-secondary btn-sm")
          )
        ),

        fluidRow(
          box(width = 4, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("map-marked-alt"), " MSR-Mantel parameters"),
            uiOutput(ns("msr_col_x_ui")),
            uiOutput(ns("msr_col_y_ui")),
            uiOutput(ns("msr_col_w_ui")),
            numericInput(ns("msr_k"), "PCoA axes for Matrix X (k):", value = 2, min = 1, max = 10, step = 1),
            numericInput(ns("msr_n_perm"), "MSR replicates:", value = 99, min = 49, max = 999, step = 50),
            tags$p(style="color:#777;font-size:11px;",
              "W is a row-standardized inverse-distance matrix built from the chosen distance column ",
              "(a simple default \u2014 see the module's info panel for caveats vs. graph-based W)."),
            actionButton(ns("run_msr"), "Run MSR-Mantel",
                         icon = icon("rocket"), class = "btn-action-primary btn-block",
                         style = "font-weight:bold;")
          ),
          box(width = 8, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("chart-bar"), " MSR-Mantel results"),
            fluidRow(
              column(4, valueBoxOutput(ns("box_msr_robs"), width = NULL)),
              column(4, valueBoxOutput(ns("box_msr_corrected"), width = NULL)),
              column(4, valueBoxOutput(ns("box_msr_p"), width = NULL))
            ),
            valueBoxOutput(ns("box_msr_n"), width = 12),
            uiOutput(ns("ui_msr_summary")),
            plotly::plotlyOutput(ns("msr_hist"), height = "320px")
          )
        )
      )
    )
  )
}
