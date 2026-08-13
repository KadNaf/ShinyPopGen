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
                    icon("chart-bar"), " Summary "),
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
              checkboxInput(ns("ibd_ext_header"), "File has header row", value = TRUE),
              uiOutput(ns("ibd_ext_file_status"))
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
                uiOutput(ns("ibd_dgeo_file_status")),
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
            "distances in columns. Supports ",
            "<b>rectangular matrices</b>: pairs can be excluded (e.g. to keep contemporaneous ",
            "pairs only) without dropping every pair involving the corresponding sub-samples. ",
            "Permutation is by <b>joint row/column relabelling</b> of one matrix, which stays ",
            "valid when either matrix is incomplete. Statistic: Pearson's correlation, Spearman's ",
            "rank correlation, or the slope of an ordinary least-squares regression (Rousset 1997, ",
            "for isolation-by-distance). ",
            "One-sided p-value = (b+1)/(m+1) (b = number of permuted statistics at least as extreme ",
            "as observed, m = number of permutations) \u2014 a bias-corrected proportion that avoids a ",
            "p-value of exactly 0 or 1 from a finite number of replicates (Davison &amp; Hinkley 1997). ",
            "A companion tab uses the same test with a plain (uncorrected) proportion instead, kept ",
            "separate so the two formulas are never mixed up.",
            "<br>Two computation engines are available: a <b>C++</b> engine running the permutation loop ",
            "natively for speed, and a fully portable <b>R</b> engine giving the same statistics.",
            "<br><b>Tip:</b> for isolation-by-distance specifically, the regression-slope statistic is usually ",
            "run on linearised genetic distance (e.g. F<sub>ST</sub>/(1-F<sub>ST</sub>)) against ",
            "<b>ln(geographic distance)</b> (2D habitat model) \u2014 make sure you selected the slope statistic, ",
            "picked a linearised distance as Y, and ticked \"ln(transform) X\" (or chose an already-logged ",
            "distance column) if that is the model you intend to test."
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
              checkboxInput(ns("mt_header"), "File has header row", value = TRUE),
              uiOutput(ns("mt_file_status"))
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
                checkboxInput(ns("mt_extra_header"), "File has header row", value = TRUE),
                uiOutput(ns("mt_extra_file_status"))
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
                checkboxInput(ns("mt_log_x"), "ln(transform) X", value = FALSE),
                uiOutput(ns("mt_double_log_warning"))
              ),
              column(4,
                numericInput(ns("mt_n_perm"), "Permutations:",
                             value = 10000, min = 99, max = 200000, step = 1000),
                tags$p(style="color:#777;font-size:11px;", "Advised \u2265 1000."),
                radioButtons(ns("mt_engine"), "Engine:",
                  choices = c("C++ \u2014 native, faster" = "cpp",
                              "R \u2014 portable" = "r"),
                  selected = "cpp"),
                tags$p(style="color:#777;font-size:11px;",
                  icon("lock"), " p-value formula: ", tags$strong("(b+1)/(m+1) (bias-corrected proportion)")),
                tags$div(style = "display:none;",
                  numericInput(ns("mt_seed"), "Random seed:", value = 67144630, min = 1, max = 2147483647, step = 1)
                )
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
                          icon("table"), " Result summary"),
            DT::DTOutput(ns("dt_mantel_summary"))
          ),
          box(width = 6, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("table"), " Null distribution quantiles (permutations)"),
            DT::DTOutput(ns("dt_mantel_quantiles")),
            tags$p(style="color:#777;font-size:11px;margin-top:6px;",
              "The observed statistic can be compared directly against these percentile thresholds of the ",
              "permuted null distribution.")
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
      # TAB 2b — Mantel Test (b/m): dedicated tab for the plain-proportion
      # p-value formula, kept separate from the generic Mantel Test tab
      # above (which uses the bias-corrected formula) so the two never get
      # mixed up / accidentally compared under the wrong setting.
      # ══════════════════════════════════════════════════════════════════
      tabPanel(title = tagList(icon("divide"), " Mantel Test (b/m)"), value = "tab_mantel_genepop",
        br(),

        tags$div(
          class = "spg-method-note", style = "border-left-color:#0c4a6e;",
          HTML(paste0(
            "Same Mantel permutation test as the previous tab, using a different one-sided p-value formula: ",
            "a <b>plain proportion</b>, <code>p = b/m</code> (b = number of permuted statistics at least as ",
            "extreme as observed, m = number of permutations) \u2014 no bias correction is applied here, unlike ",
            "the previous tab's <code>(b+1)/(m+1)</code>. Kept in its own tab so the two formulas are never ",
            "mixed up when comparing results; the test statistic and permutation scheme (joint row/column ",
            "relabelling) are otherwise identical.<br>",
            "Two computation engines are available: a <b>C++</b> engine running the permutation loop natively ",
            "for speed, and a fully portable <b>R</b> engine giving the same statistics via R's own permutation ",
            "mechanism \u2014 pick whichever suits your workflow."
          ))
        ),

        fluidRow(
          box(width = 4, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("database"), " Data source"),
            radioButtons(ns("gf_source"), NULL,
              choices = c("Null Alleles module (this session)" = "internal",
                          "Upload a file"                       = "upload"),
              selected = "internal"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'upload'", ns("gf_source")),
              fileInput(ns("gf_file"), "Browse\u2026 (pairwise file: Pop1, Pop2, distances)",
                        accept = c(".csv", ".txt", ".tsv")),
              radioButtons(ns("gf_sep"), "Separator:",
                choices = c("Tab" = "\t", "Comma" = ",", "Semicolon" = ";"),
                selected = "\t", inline = TRUE),
              checkboxInput(ns("gf_header"), "File has header row", value = TRUE)
            ),
            uiOutput(ns("gf_file_status")),
            tags$hr(),
            uiOutput(ns("gf_col_pop1_ui")),
            uiOutput(ns("gf_col_pop2_ui"))
          ),
          box(width = 8, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("sliders-h"), " Mantel parameters"),
            fluidRow(
              column(4,
                uiOutput(ns("gf_col_x_ui")),
                uiOutput(ns("gf_col_y_ui"))
              ),
              column(4,
                radioButtons(ns("gf_stat"), "Statistic:",
                  choices = c("Pearson r" = "r", "Spearman rho" = "spearman",
                              "Regression slope (Rousset)" = "b"),
                  selected = "r"),
                checkboxInput(ns("gf_log_x"), "ln(transform) X", value = FALSE),
                uiOutput(ns("gf_double_log_warning"))
              ),
              column(4,
                radioButtons(ns("gf_engine"), "Engine:",
                  choices = c("C++ \u2014 native, faster" = "cpp",
                              "R \u2014 portable" = "r"),
                  selected = "cpp"),
                numericInput(ns("gf_n_perm"), "Permutations:",
                             value = 10000, min = 99, max = 200000, step = 1000),
                tags$p(style="color:#777;font-size:11px;",
                  icon("lock"), " p-value formula: ", tags$strong("b/m (plain proportion, no correction)")),
                tags$div(style = "display:none;",
                  numericInput(ns("gf_seed"), "Random seed:", value = 67144630, min = 1, max = 2147483647, step = 1)
                ),
                actionButton(ns("run_gf_mantel"), "Run Mantel Test",
                             icon = icon("random"), class = "btn-action-primary btn-block",
                             style = "font-weight:bold;")
              )
            ),
            tags$hr(),
            fluidRow(
              column(6, DT::DTOutput(ns("dt_gf_summary"))),
              column(6, DT::DTOutput(ns("dt_gf_quantiles")))
            )
          )
        ),

        fluidRow(
          box(width = 12, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("table"), " Data used in the last run"),
            DT::DTOutput(ns("dt_gf_data")),
            tags$br(),
            downloadButton(ns("dl_gf_csv"), "Download data used", class = "btn-action-secondary btn-sm")
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
            "Two ways to test for an association between one matrix and another while controlling for a ",
            "third (or more) \u2014 pick whichever fits how many matrices you need to control for.<br>",
            "<b>MRM</b> (multiple regression on distance matrices; Legendre, Lapointe &amp; Casgrain 1994) fits ",
            "a joint multiple regression of the response on <b>up to 10 predictor matrices at once</b>, testing ",
            "each predictor's own regression coefficient by permuting the response (Y) matrix.<br>",
            "<b>Classic partial Mantel</b> uses the Yule (1907) partial-correlation formula ",
            "<code>(rxy \u2212 rxz\u00b7ryz) / \u221a(1\u2212rxz\u00b2) / \u221a(1\u2212ryz\u00b2)</code> between the variable of interest ",
            "(X) and the response (Y) while controlling for exactly <b>one</b> matrix (Z), tested by permuting X's ",
            "row/column labels while Y and Z stay fixed. Limited to a single control matrix, but reports a partial ",
            "correlation coefficient rather than a regression slope \u2014 the two methods answer a related but ",
            "genuinely different statistical question, so their numbers are not meant to match each other.<br>",
            "Both support the same <b>rectangular/incomplete data</b> as the Mantel Test tab, and both are ",
            "available with a native <b>C++</b> engine (faster) or a portable <b>R</b> engine (same statistics).",
            "<br><b>Caveat (applies to both methods):</b> partial Mantel tests can have inflated type I error when ",
            "the matrix being partialled out (e.g. geographic distance) is itself spatially autocorrelated (Guillot ",
            "&amp; Rousset 2013; Crabot et al. 2019, <i>Methods Ecol Evol</i> 10:532\u2013540). A plain Mantel test per ",
            "predictor (previous tabs) and, in the future, the Procrustes association metric (Lisboa et al. 2014, ",
            "<i>PLoS ONE</i> 9(6):e101238) as a less-controversial alternative, are worth cross-checking any result ",
            "against. Borcard &amp; Legendre (2012, <i>Ecology</i> 93:1473\u20131481) found the plain Mantel test itself ",
            "has acceptable power for most ecological uses."
          ))
        ),

        fluidRow(
          box(width = 4, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("database"), " Data source"),
            radioButtons(ns("pm_source"), NULL,
              choices = c("Null Alleles module (this session)" = "internal",
                          "Upload a file"                       = "upload"),
              selected = "internal"),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'upload'", ns("pm_source")),
              fileInput(ns("pm_file"), "Browse\u2026 (pairwise file: Pop1, Pop2, distances)",
                        accept = c(".csv", ".txt", ".tsv")),
              radioButtons(ns("pm_sep"), "Separator:",
                choices = c("Tab" = "\t", "Comma" = ",", "Semicolon" = ";"),
                selected = "\t", inline = TRUE),
              checkboxInput(ns("pm_header"), "File has header row", value = TRUE)
            ),
            uiOutput(ns("pm_file_status")),
            tags$hr(),
            uiOutput(ns("pm_col_pop1_ui")),
            uiOutput(ns("pm_col_pop2_ui")),
            tags$hr(),
            radioButtons(ns("pm_method"), "Method:",
              choices = c("MRM \u2014 up to 10 predictors at once" = "mrm",
                          "Classic partial Mantel \u2014 1 control matrix" = "classic"),
              selected = "mrm"),
            tags$hr(),
            uiOutput(ns("pm_col_y_ui")),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'mrm'", ns("pm_method")),
              uiOutput(ns("pm_col_x_ui")),
              checkboxInput(ns("pm_standardize"), "Standardize (z-score) variables before fitting",
                            value = FALSE),
              tags$p(style="color:#777;font-size:11px;",
                "Off by default so coefficients are reported in the original (raw) units of each matrix. Turn ",
                "on only if you want to compare the relative strength of predictors measured in different units.")
            ),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'classic'", ns("pm_method")),
              uiOutput(ns("pm_col_x1_ui")),
              uiOutput(ns("pm_col_z_ui")),
              radioButtons(ns("pm_classic_stat"), "Correlation method:",
                choices = c("Pearson" = "pearson", "Spearman" = "spearman"), selected = "pearson", inline = TRUE)
            ),
            numericInput(ns("pm_n_perm"), "Permutations:", value = 999, min = 99, max = 20000, step = 100),
            radioButtons(ns("pm_engine"), "Engine:",
              choices = c("C++ \u2014 native, faster" = "cpp",
                          "R \u2014 portable" = "r"),
              selected = "cpp"),
            tags$p(style="color:#777;font-size:11px;",
              "Applies to whichever method is selected above. MRM permutes the response (Y) matrix; classic ",
              "mode permutes the variable-of-interest (X) matrix. Both use joint row/column relabelling, ",
              "valid on incomplete data."),
            tags$div(style = "display:none;",
              numericInput(ns("pm_seed"), "Random seed:", value = 67144630, min = 1, max = 2147483647, step = 1)
            ),
            actionButton(ns("run_partial_mantel"), "Run Partial Mantel",
                         icon = icon("random"), class = "btn-action-primary btn-block",
                         style = "font-weight:bold;")
          ),
          box(width = 8, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("sliders-h"), " Results"),
            DT::DTOutput(ns("dt_partial_mantel")),
            uiOutput(ns("ui_partial_mantel_summary")),
            tags$br(),
            downloadButton(ns("dl_partial_mantel_csv"), "Download results (CSV)", class = "btn-action-secondary btn-sm"),
            tags$hr(),
            tags$div(style="font-weight:600;color:#333a43;margin-bottom:6px;",
                     "Correlation matrix among Y, X and Z (collinearity check)"),
            DT::DTOutput(ns("dt_pm_corr")),
            tags$p(style="color:#777;font-size:11px;margin-top:6px;",
              "High correlations (|r| > ~0.7) between two matrices make coefficients involving them ",
              "unstable/hard to interpret \u2014 check this table before trusting a single predictor's ",
              "coefficient in the results above.")
          )
        )
      )
    )
  )
}
