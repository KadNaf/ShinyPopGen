# ui_isolation_by_distance.R
# Isolation by Distance & Mantel Test — SPG-V1 specification
#
# Three tabs:
#   1. Pairwise Distances — load external dataset OR compute Dgeo from GPS
#   2. IBD Regression (Rousset 1997) — 1D/2D models, FR/FR-ENA, 3 regression lines
#   3. Mantel Test — rectangular matrices, Pearson or Rousset slope, p=(b+1)/(m+1)

isolation_by_distance_UI <- function(id) {
  ns <- NS(id)

  custom_css <- tags$style(HTML("
    @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap');

    .ibd-module * { font-family: 'IBM Plex Sans', sans-serif; }

    /* ── Header ─────────────────────────────────────────────────────── */
    .ibd-header {
      background: linear-gradient(135deg, #0f172a 0%, #1e293b 55%, #0c4a6e 100%);
      border-radius: 10px; padding: 1.2rem 1.6rem; margin-bottom: 1rem;
      position: relative; overflow: hidden;
    }
    .ibd-header::before {
      content: ''; position: absolute; inset: 0;
      background: repeating-linear-gradient(
        -45deg, transparent, transparent 28px,
        rgba(255,255,255,.018) 28px, rgba(255,255,255,.018) 29px);
    }
    .ibd-header-title { font-size:1.05rem; font-weight:600; color:#f1f5f9; letter-spacing:.01em; margin-bottom:.2rem; }
    .ibd-header-sub   { font-size:.75rem; color:#94a3b8; font-family:'IBM Plex Mono',monospace; }
    .ibd-badges { display:flex; gap:6px; margin-top:.5rem; flex-wrap:wrap; }
    .ibd-badge  { display:inline-block; border-radius:20px; padding:2px 10px; font-size:.67rem; font-family:'IBM Plex Mono',monospace; }
    .ibd-badge-blue   { background:rgba(56,189,248,.15);  border:1px solid rgba(56,189,248,.3);  color:#38bdf8; }
    .ibd-badge-green  { background:rgba(74,222,128,.12);  border:1px solid rgba(74,222,128,.3);  color:#4ade80; }
    .ibd-badge-amber  { background:rgba(251,191,36,.12);  border:1px solid rgba(251,191,36,.3);  color:#fbbf24; }
    .ibd-badge-teal   { background:rgba(20,184,166,.15);  border:1px solid rgba(20,184,166,.3);  color:#2dd4bf; }
    .ibd-badge-purple { background:rgba(168,85,247,.15);  border:1px solid rgba(168,85,247,.3);  color:#a855f7; }
    .ibd-badge-rose   { background:rgba(244,63,94,.15);   border:1px solid rgba(244,63,94,.3);   color:#fb7185; }

    /* ── Value boxes ─────────────────────────────────────────────────── */
    .ibd-vbox-row { display:flex; gap:9px; margin-bottom:1rem; flex-wrap:wrap; }
    .ibd-vbox { flex:1; min-width:110px; background:#fff; border:1px solid #e2e8f0; border-radius:9px; padding:.6rem .85rem; display:flex; align-items:center; gap:9px; }
    .ibd-vbox-icon  { width:30px; height:30px; border-radius:7px; display:flex; align-items:center; justify-content:center; font-size:12px; flex-shrink:0; }
    .ibd-vbox-label { font-size:10px; color:#94a3b8; text-transform:uppercase; letter-spacing:.06em; margin-bottom:1px; }
    .ibd-vbox-val   { font-size:18px; font-weight:600; color:#0f172a; line-height:1.1; font-family:'IBM Plex Mono',monospace; }

    /* ── Panels ──────────────────────────────────────────────────────── */
    .ibd-panel { background:#fff; border:1px solid #e2e8f0; border-radius:9px; margin-bottom:.85rem; overflow:hidden; }
    .ibd-panel-head { background:#f8fafc; border-bottom:1px solid #e2e8f0; padding:.55rem .9rem; }
    .ibd-panel-title { font-size:12px; font-weight:600; color:#1e293b; display:flex; align-items:center; gap:6px; flex-wrap:wrap; }
    .ibd-panel-body { padding:.85rem; }

    /* ── Info strips ─────────────────────────────────────────────────── */
    .ibd-info    { background:#eff6ff; border:1px solid #bfdbfe; border-radius:7px; padding:.45rem .8rem; font-size:11.5px; color:#1d4ed8; margin-bottom:.85rem; line-height:1.65; }
    .ibd-warn    { background:#fffbeb; border:1px solid #fcd34d; border-radius:7px; padding:.45rem .8rem; font-size:11.5px; color:#92400e; margin-bottom:.85rem; line-height:1.65; }
    .ibd-success { background:#f0fdf4; border:1px solid #86efac; border-radius:7px; padding:.45rem .8rem; font-size:11.5px; color:#166534; margin-bottom:.85rem; line-height:1.65; }
    .ibd-result  { background:#faf5ff; border:1px solid #d8b4fe; border-radius:8px; padding:.65rem 1rem; font-size:11.5px; color:#3b0764; font-family:'IBM Plex Mono',monospace; line-height:1.9; margin-top:.75rem; }
    .ibd-result strong { color:#6d28d9; }

    /* ── Matrix table ────────────────────────────────────────────────── */
    .ibd-matrix-wrap { overflow-x:auto; margin-top:.5rem; }
    .ibd-matrix { border-collapse:collapse; font-size:11px; font-family:'IBM Plex Mono',monospace; width:100%; }
    .ibd-matrix th { background:#f8fafc; color:#475569; font-weight:600; padding:4px 9px; border:1px solid #e2e8f0; font-size:10.5px; white-space:nowrap; }
    .ibd-matrix td { padding:4px 9px; border:1px solid #e2e8f0; color:#1e293b; text-align:right; white-space:nowrap; font-size:11px; }
    .ibd-matrix tr:nth-child(even) td { background:#f8fafc; }
    .ibd-matrix .diag  { background:#f1f5f9 !important; color:#94a3b8; text-align:center; }
    .ibd-matrix .upper { color:#cbd5e1; text-align:center; }
    .ibd-matrix .lbl   { font-weight:700; color:#0f172a; text-align:left; white-space:nowrap; }

    /* ── Buttons ─────────────────────────────────────────────────────── */
    .ibd-btn-run {
      background:linear-gradient(135deg,#0369a1,#0c4a6e) !important;
      border:none !important; color:#fff !important; border-radius:7px !important;
      font-weight:600 !important; font-size:13px !important; padding:7px 22px !important;
      box-shadow:0 2px 8px rgba(3,105,161,.3) !important;
    }
    .ibd-btn-run:hover { opacity:.9; }
    .ibd-btn-mantel {
      background:linear-gradient(135deg,#7c3aed,#4c1d95) !important;
      border:none !important; color:#fff !important; border-radius:7px !important;
      font-weight:600 !important; font-size:13px !important; padding:7px 22px !important;
      box-shadow:0 2px 8px rgba(124,58,237,.3) !important;
    }
    .ibd-btn-mantel:hover { opacity:.9; }

    /* ── Download row ────────────────────────────────────────────────── */
    .ibd-dl-row { display:flex; gap:6px; flex-wrap:wrap; margin-top:.5rem; }
    .ibd-dl-row .btn { font-size:11px; padding:3px 12px; }

    /* ── DT tweaks ───────────────────────────────────────────────────── */
    .ibd-module .dataTables_wrapper { font-size:12px; }
    .ibd-module table.dataTable thead th {
      background:#f8fafc !important; color:#475569 !important;
      font-family:'IBM Plex Mono',monospace !important;
      font-size:10.5px !important; font-weight:600 !important;
    }
    .ibd-module table.dataTable tbody td {
      font-family:'IBM Plex Mono',monospace !important;
      font-size:11px !important; color:#1e293b !important;
    }
    .ibd-module .nav-tabs > li > a { font-size:12px; font-weight:500; color:#475569; padding:5px 13px; }
    .ibd-module .nav-tabs > li.active > a { color:#0f172a; font-weight:600; }
  "))

  tags$div(class="ibd-module", custom_css,

    # ── Header ─────────────────────────────────────────────────────────────
    tags$div(class="ibd-header",
      tags$div(class="ibd-header-title",
        icon("map-marked-alt"), " Isolation by Distance & Mantel Test"),
      tags$div(class="ibd-header-sub",
        "Rousset (1997) regression \u00b7 geosphere Haversine \u00b7 Mantel (1967) permutation test \u00b7 ",
        "Rectangular matrices (RT / Fstat format)"),
      tags$div(class="ibd-badges",
        tags$span(class="ibd-badge ibd-badge-blue",   "FR = FST/(1\u2212FST)"),
        tags$span(class="ibd-badge ibd-badge-teal",   "FR-ENA (null alleles)"),
        tags$span(class="ibd-badge ibd-badge-green",  "1D / 2D models"),
        tags$span(class="ibd-badge ibd-badge-amber",  "CI threshold tuning"),
        tags$span(class="ibd-badge ibd-badge-purple", "Mantel \u2014 Pearson / Rousset b"),
        tags$span(class="ibd-badge ibd-badge-rose",   "p = (b+1)/(m+1)")
      )
    ),

    # ════════════════════════════════════════════════════════════════════════
    # 3 TABS
    # ════════════════════════════════════════════════════════════════════════
    tabsetPanel(id = ns("ibd_tabs"), type = "tabs",

      # ══════════════════════════════════════════════════════════════════════
      # TAB 1: PAIRWISE DATASET
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(title = tagList(icon("table"), " Pairwise Dataset"),
               value = "tab_data", br(),

        tags$div(class="ibd-info",
          icon("info-circle"), " ",
          tags$strong("Pairwise subsample distances."), " Load an external dataset with headers, ",
          "or compute geographic distances (Dgeo) from GPS coordinates (decimal degrees) using the ",
          tags$em("geosphere"), " package (Hijmans et al. 2019). ",
          "Other distance types (temporal, ecological, categorical) can be loaded as additional columns. ",
          "The user can then filter pairs (e.g. keep only contemporaneous pairs for IBD testing)."
        ),

        # ── Data source ─────────────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("folder-open"), " Data source")),
          tags$div(class="ibd-panel-body",
            fluidRow(
              column(6,
                radioButtons(ns("data_source"), "Source of pairwise data:",
                  choices = c(
                    "Load external dataset (CSV/TXT)" = "external",
                    "Compute Dgeo from imported GPS"  = "gps"
                  ),
                  selected = "external"),
                conditionalPanel(
                  condition = "input.data_source == 'external'", ns = ns,
                  fileInput(ns("file_ext"), "Pairwise distance file:",
                            accept = c(".csv", ".txt", ".tab")),
                  tags$div(class="ibd-info", style="margin-top:.5rem;margin-bottom:0;",
                    icon("info-circle"), " ",
                    "Expected columns: ", tags$strong("Subsample1, Subsample2"),
                    " (mandatory), then any combination of: ",
                    tags$code("Dgeo"), ", ", tags$code("lnDgeo"), ", ",
                    tags$code("FST"), ", ", tags$code("FST_ENA"), ", ",
                    tags$code("FR"), ", ", tags$code("FR_ENA"), ", ",
                    tags$code("FR_l"), ", ", tags$code("FR_u"), ", ",
                    tags$code("DCSE"), ", ", tags$code("DCSE_INA"), ", etc."
                  )
                ),
                conditionalPanel(
                  condition = "input.data_source == 'gps'", ns = ns,
                  tags$div(class="ibd-info", style="margin-bottom:0;",
                    icon("info-circle"), " ",
                    "Dgeo will be computed from the GPS coordinates (Latitude/Longitude in decimal degrees) ",
                    "of each population using ", tags$code("geosphere::distHaversine"), "."
                  )
                )
              ),
              column(6,
                tags$div(class="ibd-vbox-row",
                  tags$div(class="ibd-vbox",
                    tags$div(class="ibd-vbox-icon",style="background:#e0f2fe;color:#0369a1;",icon("project-diagram")),
                    tags$div(tags$div(class="ibd-vbox-label","Pairs"),
                             tags$div(class="ibd-vbox-val",uiOutput(ns("vb_pairs"))))),
                  tags$div(class="ibd-vbox",
                    tags$div(class="ibd-vbox-icon",style="background:#dcfce7;color:#166534;",icon("map-marker-alt")),
                    tags$div(tags$div(class="ibd-vbox-label","Subsamples"),
                             tags$div(class="ibd-vbox-val",uiOutput(ns("vb_pops"))))),
                  tags$div(class="ibd-vbox",
                    tags$div(class="ibd-vbox-icon",style="background:#fef9c3;color:#854d0e;",icon("columns")),
                    tags$div(tags$div(class="ibd-vbox-label","Columns"),
                             tags$div(class="ibd-vbox-val",uiOutput(ns("vb_cols")))))
                )
              )
            )
          )
        ),

        # ── Pair filter ─────────────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("filter"), " Filter subsample pairs (optional)")),
          tags$div(class="ibd-panel-body",
            tags$div(class="ibd-warn",
              icon("exclamation-triangle"), " ",
              "Isolation by geographic distance must be tested between ", tags$strong("contemporaneous pairs only"),
              " to avoid confusing temporal effects on genetic distances. ",
              "Use this filter to exclude non-contemporaneous pairs, or pairs involving specific subsamples."
            ),
            fluidRow(
              column(6,
                selectizeInput(ns("filter_pop1"), "Include pairs where Subsample 1 is in:",
                               choices = NULL, multiple = TRUE,
                               options = list(placeholder = "Leave empty to include all")),
                selectizeInput(ns("filter_pop2"), "Include pairs where Subsample 2 is in:",
                               choices = NULL, multiple = TRUE,
                               options = list(placeholder = "Leave empty to include all"))
              ),
              column(6,
                numericInput(ns("dgeo_min"), "Min Dgeo (km):", value = NA, min = 0),
                numericInput(ns("dgeo_max"), "Max Dgeo (km):", value = NA, min = 0)
              )
            ),
            tags$hr(),
            actionButton(ns("apply_filter"), "Apply filter",
                         icon = icon("filter"), class = "ibd-btn-run btn")
          )
        ),

        # ── Data table ──────────────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("list"), " Pairwise dataset (filtered)")),
          tags$div(class="ibd-panel-body",
            DT::DTOutput(ns("dt_pairwise"))),
          tags$div(class="ibd-panel-body",
            tags$div(class="ibd-dl-row",
              downloadButton(ns("dl_pairwise_csv"), "Download CSV", class = "btn btn-default btn-sm"),
              downloadButton(ns("dl_pairwise_txt"), "Download TXT", class = "btn btn-default btn-sm"))))
      ),

      # ══════════════════════════════════════════════════════════════════════
      # TAB 2: IBD REGRESSION (Rousset 1997)
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(title = tagList(icon("chart-line"), " IBD Regression (Rousset 1997)"),
               value = "tab_ibd", br(),

        tags$div(class="ibd-info",
          icon("info-circle"), " ",
          tags$strong("Rousset's (1997) regression for isolation by distance."), " ",
          "Two models: ", tags$strong("Model 1D"), " (FR ~ Dgeo) for one-dimensional IBD; ",
          tags$strong("Model 2D"), " (FR ~ ln(Dgeo)) for two-dimensional IBD at migration-mutation-drift equilibrium. ",
          "Genetic distance is ", tags$strong("FR = FST/(1\u2212FST)"), " or ",
          tags$strong("FR-ENA = FST-ENA/(1\u2212FST-ENA)"), " if null alleles are present. ",
          "Three regression lines are fitted: through FR, FR-l (lower CI) and FR-u (upper CI). ",
          tags$br(), tags$br(),
          tags$strong("Original feature:"), " the CI level of bootstraps can be tuned to obtain the exact threshold ",
          "for significance of the slope. If the slope and its CI are all positive \u2192 IBD confirmed. ",
          "Reference: Rousset F. 1997. Genetics 145:1219\u20131228."
        ),

        # ── Configuration ───────────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("sliders-h"), " IBD regression configuration")),
          tags$div(class="ibd-panel-body",
            fluidRow(
              column(4,
                radioButtons(ns("ibd_model"), "IBD model:",
                  choices = c(
                    "Model 1D \u2014 FR ~ Dgeo"       = "1D",
                    "Model 2D \u2014 FR ~ ln(Dgeo)"   = "2D"
                  ),
                  selected = "2D"),
                radioButtons(ns("ibd_fr"), "Genetic distance:",
                  choices = c(
                    "FR     = FST / (1 \u2212 FST)"         = "FR",
                    "FR-ENA = FST-ENA / (1 \u2212 FST-ENA)" = "FR_ENA"
                  ),
                  selected = "FR")
              ),
              column(4,
                selectInput(ns("ibd_col_pop1"), "Column: Subsample 1",
                            choices = NULL),
                selectInput(ns("ibd_col_pop2"), "Column: Subsample 2",
                            choices = NULL),
                selectInput(ns("ibd_col_dgeo"), "Column: geographic distance",
                            choices = NULL),
                selectInput(ns("ibd_col_fr"),   "Column: genetic distance (FR or FR-ENA)",
                            choices = NULL)
              ),
              column(4,
                selectInput(ns("ibd_col_frl"), "Column: FR lower CI (FR-l)",
                            choices = NULL),
                selectInput(ns("ibd_col_fru"), "Column: FR upper CI (FR-u)",
                            choices = NULL),
                tags$hr(),
                actionButton(ns("run_ibd"), "Run IBD Regression",
                             icon = icon("chart-line"),
                             class = "ibd-btn-run btn", width = "100%")
              )
            )
          )
        ),

        # ── CI threshold tuning (ORIGINAL FEATURE) ──────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("sliders"), " CI threshold tuning \u2014 find the exact significance level")),
          tags$div(class="ibd-panel-body",
            tags$div(class="ibd-warn",
              icon("lightbulb"), " ",
              tags$strong("Original feature of SPG-V1."), " Since the level of confidence interval of bootstraps ",
              "can be modulated, one can compute the exact percentage needed to get the threshold value for ",
              "significance. Move the slider to find the CI level at which the lower regression line (FR-l) ",
              "has a slope that just crosses zero \u2014 this gives the exact p-value of the IBD test."
            ),
            fluidRow(
              column(6,
                sliderInput(ns("ci_threshold"), "Confidence level for threshold finding:",
                            min = 0.50, max = 0.999, value = 0.95, step = 0.001)
              ),
              column(6,
                uiOutput(ns("ui_ibd_threshold_result"))
              )
            )
          )
        ),

        # ── Regression results ──────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("chart-bar"), " Regression parameters")),
          tags$div(class="ibd-panel-body",
            fluidRow(
              column(3, tags$div(class="ibd-vbox",
                tags$div(class="ibd-vbox-icon",style="background:#dcfce7;color:#166534;",icon("check-circle")),
                tags$div(tags$div(class="ibd-vbox-label","IBD status"),
                         tags$div(class="ibd-vbox-val",uiOutput(ns("vb_ibd_status")))))),
              column(3, tags$div(class="ibd-vbox",
                tags$div(class="ibd-vbox-icon",style="background:#e0f2fe;color:#0369a1;",icon("sort-amount-up")),
                tags$div(tags$div(class="ibd-vbox-label","Slope b"),
                         tags$div(class="ibd-vbox-val",uiOutput(ns("vb_ibd_b")))))),
              column(3, tags$div(class="ibd-vbox",
                tags$div(class="ibd-vbox-icon",style="background:#fef9c3;color:#854d0e;",icon("users")),
                tags$div(tags$div(class="ibd-vbox-label","Nb = 1/b"),
                         tags$div(class="ibd-vbox-val",uiOutput(ns("vb_ibd_nb")))))),
              column(3, tags$div(class="ibd-vbox",
                tags$div(class="ibd-vbox-icon",style="background:#f3e8ff;color:#7e22ce;",icon("exchange-alt")),
                tags$div(tags$div(class="ibd-vbox-label","Nem = 1/(2\u03c0b)"),
                         tags$div(class="ibd-vbox-val",uiOutput(ns("vb_ibd_nem"))))))
            ),
            DT::DTOutput(ns("dt_ibd_reg")),
            uiOutput(ns("ui_ibd_interpretation"))
          )
        ),

        # ── IBD plot ────────────────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("chart-area"), " IBD regression plot")),
          tags$div(class="ibd-panel-body",
            plotly::plotlyOutput(ns("ibd_plot"), height = "520px"))),

        # ── Download ────────────────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("file-download"), " Download IBD regression results")),
          tags$div(class="ibd-panel-body",
            tags$div(class="ibd-dl-row",
              downloadButton(ns("dl_ibd_csv"), "Regression table (CSV)", class = "btn btn-default btn-sm"))))
      ),

      # ══════════════════════════════════════════════════════════════════════
      # TAB 3: MANTEL TEST
      # ══════════════════════════════════════════════════════════════════════
      tabPanel(title = tagList(icon("project-diagram"), " Mantel Test"),
               value = "tab_mantel", br(),

        tags$div(class="ibd-info",
          icon("info-circle"), " ",
          tags$strong("Mantel test (Mantel 1967)."), " Assesses the correlation between two distance matrices ",
          "measured between any pair of subsamples. Distances can be genetic (FST, FST-ENA, FR, FR-ENA, DCSE, DCSE-INA), ",
          "geographic, temporal, ecological or categorical (0/1). ",
          tags$br(), tags$br(),
          tags$strong("Rectangular matrices"), " (as in RT / Fstat 2.9.4) are supported: ",
          "this is convenient to exclude some pairs without excluding all pairs involving a subsample. ",
          tags$br(), tags$br(),
          tags$strong("Statistic used:"), " Pearson's correlation coefficient (as in Fstat 2.9.4) ",
          "OR the slope of Rousset's regression (as in Genepop). ",
          "Cells of one matrix are randomized m times; the one-sided p-value is computed as ",
          tags$code("p = (b+1)/(m+1)"), " where b = number of permutations with statistic \u2265 observed. ",
          "References: Mantel 1967; Manly 2018 (RT); Séré et al. 2017."
        ),

        # ── Configuration ───────────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("sliders-h"), " Mantel test configuration")),
          tags$div(class="ibd-panel-body",
            fluidRow(
              # ── Matrix 1 ──────────────────────────────────────────────────
              column(4,
                tags$div(class="ibd-panel", style="border-color:#e9d5ff;",
                  tags$div(class="ibd-panel-head", style="background:#faf5ff;",
                    tags$div(class="ibd-panel-title", style="color:#6d28d9;",
                      icon("dna"), " Matrix 1 (X)")),
                  tags$div(class="ibd-panel-body",
                    radioButtons(ns("m1_source"), "Source:",
                      choices = c(
                        "From pairwise dataset (column)" = "col",
                        "Upload file"                   = "upload"
                      ),
                      selected = "col"),
                    conditionalPanel(
                      condition = "input.m1_source == 'col'", ns = ns,
                      selectInput(ns("m1_col"), "Column of pairwise dataset:",
                                  choices = NULL)
                    ),
                    conditionalPanel(
                      condition = "input.m1_source == 'upload'", ns = ns,
                      fileInput(ns("m1_file"), "Distance file (CSV/TXT):",
                                accept = c(".csv", ".txt", ".tab")),
                      radioButtons(ns("m1_format"), "Format:",
                        choices = c("Square matrix" = "square",
                                    "Rectangular (column-wise)" = "rectangular"),
                        selected = "rectangular")
                    )
                  )
                )
              ),

              # ── Matrix 2 ──────────────────────────────────────────────────
              column(4,
                tags$div(class="ibd-panel", style="border-color:#99f6e4;",
                  tags$div(class="ibd-panel-head", style="background:#f0fdfa;",
                    tags$div(class="ibd-panel-title", style="color:#0d9488;",
                      icon("globe"), " Matrix 2 (Y)")),
                  tags$div(class="ibd-panel-body",
                    radioButtons(ns("m2_source"), "Source:",
                      choices = c(
                        "From pairwise dataset (column)" = "col",
                        "Computed Dgeo from GPS"         = "gps_km",
                        "Computed ln(Dgeo) from GPS"     = "gps_ln",
                        "Upload file"                    = "upload"
                      ),
                      selected = "gps_km"),
                    conditionalPanel(
                      condition = "input.m2_source == 'col'", ns = ns,
                      selectInput(ns("m2_col"), "Column of pairwise dataset:",
                                  choices = NULL)
                    ),
                    conditionalPanel(
                      condition = "input.m2_source == 'upload'", ns = ns,
                      fileInput(ns("m2_file"), "Distance file (CSV/TXT):",
                                accept = c(".csv", ".txt", ".tab")),
                      radioButtons(ns("m2_format"), "Format:",
                        choices = c("Square matrix" = "square",
                                    "Rectangular (column-wise)" = "rectangular"),
                        selected = "rectangular")
                    )
                  )
                )
              ),

              # ── Test parameters ───────────────────────────────────────────
              column(4,
                tags$div(class="ibd-panel", style="border-color:#fcd34d;",
                  tags$div(class="ibd-panel-head", style="background:#fffbeb;",
                    tags$div(class="ibd-panel-title", style="color:#92400e;",
                      icon("cog"), " Test parameters")),
                  tags$div(class="ibd-panel-body",
                    radioButtons(ns("mantel_stat"), "Statistic:",
                      choices = c(
                        "Pearson r (Fstat 2.9.4)"       = "pearson",
                        "Rousset slope b (Genepop)"     = "rousset"
                      ),
                      selected = "pearson"),
                    numericInput(ns("n_perm"), "Number of permutations (m):",
                                 value = 10000, min = 1000, max = 999999, step = 1000),
                    radioButtons(ns("mantel_side"), "Alternative hypothesis:",
                      choices = c(
                        "Positive correlation (one-sided)" = "greater",
                        "Negative correlation (one-sided)" = "less",
                        "Two-sided"                        = "two.sided"
                      ),
                      selected = "greater"),
                    tags$hr(),
                    actionButton(ns("run_mantel"), "Run Mantel Test",
                                 icon = icon("play"),
                                 class = "ibd-btn-mantel btn", width = "100%")
                  )
                )
              )
            )
          )
        ),

        # ── Mantel results ──────────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("chart-bar"), " Mantel test results")),
          tags$div(class="ibd-panel-body",
            fluidRow(
              column(2, tags$div(class="ibd-vbox",
                tags$div(class="ibd-vbox-icon",style="background:#f3e8ff;color:#7e22ce;",icon("chart-line")),
                tags$div(tags$div(class="ibd-vbox-label","Statistic"),
                         tags$div(class="ibd-vbox-val",uiOutput(ns("vb_mantel_stat")))))),
              column(2, tags$div(class="ibd-vbox",
                tags$div(class="ibd-vbox-icon",style="background:#dcfce7;color:#166534;",icon("check-circle")),
                tags$div(tags$div(class="ibd-vbox-label","p-value"),
                         tags$div(class="ibd-vbox-val",uiOutput(ns("vb_mantel_p")))))),
              column(2, tags$div(class="ibd-vbox",
                tags$div(class="ibd-vbox-icon",style="background:#fef9c3;color:#854d0e;",icon("hashtag")),
                tags$div(tags$div(class="ibd-vbox-label","Pairs (n)"),
                         tags$div(class="ibd-vbox-val",uiOutput(ns("vb_mantel_n")))))),
              column(2, tags$div(class="ibd-vbox",
                tags$div(class="ibd-vbox-icon",style="background:#e0f2fe;color:#0369a1;",icon("percent")),
                tags$div(tags$div(class="ibd-vbox-label","R\u00b2 (%)"),
                         tags$div(class="ibd-vbox-val",uiOutput(ns("vb_mantel_r2")))))),
              column(2, tags$div(class="ibd-vbox",
                tags$div(class="ibd-vbox-icon",style="background:#fce7f3;color:#9d174d;",icon("sort-amount-up")),
                tags$div(tags$div(class="ibd-vbox-label","b \u2265 b_obs"),
                         tags$div(class="ibd-vbox-val",uiOutput(ns("vb_mantel_b")))))),
              column(2, tags$div(class="ibd-vbox",
                tags$div(class="ibd-vbox-icon",style="background:#ccfbf1;color:#0d9488;",icon("exchange-alt")),
                tags$div(tags$div(class="ibd-vbox-label","Pops aligned"),
                         tags$div(class="ibd-vbox-val",uiOutput(ns("vb_mantel_pops"))))))
            ),
            uiOutput(ns("ui_mantel_result")),
            br(),
            tags$div(class="ibd-dl-row",
              downloadButton(ns("dl_mantel_csv"), "Download results (CSV)",
                             class = "btn btn-default btn-sm"))
          )
        ),

        # ── Mantel scatter plot ─────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("chart-scatter"), " Mantel scatter plot")),
          tags$div(class="ibd-panel-body",
            plotly::plotlyOutput(ns("mantel_plot"), height = "500px")))
      )

    ) # end tabsetPanel
  )   # end tags$div.ibd-module
}