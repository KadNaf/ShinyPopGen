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
            "distances in columns (RT, Manly 2018; Fstat 2.9.4 convention). Supports ",
            "<b>rectangular matrices</b>: pairs can be excluded (e.g. to keep contemporaneous ",
            "pairs only) without dropping every pair involving the corresponding sub-samples. ",
            "Permutation is by <b>joint row/column relabelling</b> of one matrix, which stays ",
            "valid when either matrix is incomplete. Statistic: Pearson's r or Spearman's rho ",
            "(Fstat convention) or the slope of the Rousset (1997) regression (Genepop convention ",
            "for IBD). ",
            "One-sided p-value = (b+1)/(m+1), b = number of permuted statistics \u2265 observed.",
            "<br><b>If your numbers don't match Fstat:</b> Fstat's own Isolation-by-Distance Mantel test ",
            "uses the regression <b>slope b</b> (not Pearson r) of <b>F_R = FST/(1-FST)</b> against ",
            "<b>ln(distance)</b> (2D habitat) or raw distance (1D) \u2014 make sure you selected ",
            "\"Regression slope (Rousset)\", picked <code>FR</code>/<code>FR_raw</code> as Y, and ticked ",
            "\"ln(transform) X\" if comparing to a 2D Fstat run. Also check the Exclude-pairs field is empty ",
            "if Fstat used every pair."
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
                          icon("table"), " Result summary"),
            DT::DTOutput(ns("dt_mantel_summary"))
          ),
          box(width = 6, solidHeader = FALSE,
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("table"), " Null distribution quantiles (permutations)"),
            DT::DTOutput(ns("dt_mantel_quantiles")),
            tags$p(style="color:#777;font-size:11px;margin-top:6px;",
              "Same style of output as Fstat's permutation table: the observed statistic can be ",
              "compared directly against these percentile thresholds of the permuted null distribution.")
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
            "plain Mantel test itself has acceptable power for most ecological uses.",
            "<br><b>If your numbers don't match Fstat:</b> this is expected to some degree \u2014 Fstat's own ",
            "partial Mantel does not use the same algorithm as this tab (MRM, a joint multiple regression ",
            "across all matrices at once). Coefficients are reported in <b>raw units by default</b> ",
            "(untick \"Standardize\" is the default) to stay comparable to Fstat's own output; check the ",
            "correlation matrix below \u2014 if two predictors are highly correlated, no algorithm's individual ",
            "coefficients for them will be directly comparable across methods."
          ))
        ),

        fluidRow(
          box(width = 4, solidHeader = TRUE, status = "primary",
              title = div(style="background:#FFFFFF;padding:10px;color:#333a43;font-weight:600;",
                          icon("database"), " Variables (from the Mantel Test tab's data source)"),
            uiOutput(ns("pm_col_y_ui")),
            uiOutput(ns("pm_col_x_ui")),
            tags$hr(),
            checkboxInput(ns("pm_standardize"), "Standardize (z-score) variables before fitting",
                          value = FALSE),
            tags$p(style="color:#777;font-size:11px;",
              "Off by default so coefficients are in the same raw units Fstat reports. Turn on only if ",
              "you want to compare the relative strength of predictors measured in different units."),
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
            downloadButton(ns("dl_partial_mantel_csv"), "Download results (CSV)", class = "btn-action-secondary btn-sm"),
            tags$hr(),
            tags$div(style="font-weight:600;color:#333a43;margin-bottom:6px;",
                     "Correlation matrix among Y and all predictors (collinearity check)"),
            DT::DTOutput(ns("dt_pm_corr")),
            tags$p(style="color:#777;font-size:11px;margin-top:6px;",
              "High correlations (|r| > ~0.7) between two predictors make their individual coefficients ",
              "unstable/hard to interpret \u2014 check this table before trusting a single predictor's ",
              "coefficient in the results above.")
          )
        )
      )
    )
  )
}
