# ui_isolation_by_distance.R
# Isolation by Distance module
# Rousset (1997) regression: FR = FST/(1-FST) vs geographic distance
# Mantel test on rectangular matrices

isolation_by_distance_UI <- function(id) {
  ns <- NS(id)

  custom_css <- tags$style(HTML("
    @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap');

    .ibd-module * { font-family: 'IBM Plex Sans', sans-serif; }

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
    .ibd-badge-teal   { background:rgba(20,184,166,.15);  border:1px solid rgba(20,184,166,.3);  color:#2dd4bf; }
    .ibd-badge-blue   { background:rgba(56,189,248,.15);  border:1px solid rgba(56,189,248,.3);  color:#38bdf8; }
    .ibd-badge-amber  { background:rgba(251,191,36,.12);  border:1px solid rgba(251,191,36,.3);  color:#fbbf24; }
    .ibd-badge-purple { background:rgba(192,132,252,.12); border:1px solid rgba(192,132,252,.3); color:#c084fc; }
    .ibd-badge-pink   { background:rgba(244,114,182,.12); border:1px solid rgba(244,114,182,.3); color:#f472b6; }

    .ibd-vbox-row { display:flex; gap:9px; margin-bottom:1rem; flex-wrap:wrap; }
    .ibd-vbox { flex:1; min-width:110px; background:#fff; border:1px solid #e2e8f0; border-radius:9px; padding:.6rem .85rem; display:flex; align-items:center; gap:9px; }
    .ibd-vbox-icon  { width:30px; height:30px; border-radius:7px; display:flex; align-items:center; justify-content:center; font-size:12px; flex-shrink:0; }
    .ibd-vbox-label { font-size:10px; color:#94a3b8; text-transform:uppercase; letter-spacing:.06em; margin-bottom:1px; }
    .ibd-vbox-val   { font-size:18px; font-weight:600; color:#0f172a; line-height:1.1; font-family:'IBM Plex Mono',monospace; }

    .ibd-panel { background:#fff; border:1px solid #e2e8f0; border-radius:9px; margin-bottom:.85rem; overflow:hidden; }
    .ibd-panel-head { background:#f8fafc; border-bottom:1px solid #e2e8f0; padding:.55rem .9rem; }
    .ibd-panel-title { font-size:12px; font-weight:600; color:#1e293b; display:flex; align-items:center; gap:6px; flex-wrap:wrap; }
    .ibd-panel-body { padding:.85rem; }

    .ibd-info { background:#eff6ff; border:1px solid #bfdbfe; border-radius:7px; padding:.45rem .8rem; font-size:11.5px; color:#1d4ed8; margin-bottom:.85rem; line-height:1.65; }
    .ibd-warn { background:#fffbeb; border:1px solid #fcd34d; border-radius:7px; padding:.45rem .8rem; font-size:11.5px; color:#92400e; margin-bottom:.85rem; line-height:1.65; }

    .ibd-btn-run {
      background:linear-gradient(135deg,#0369a1,#0c4a6e) !important;
      border:none !important; color:#fff !important; border-radius:7px !important;
      font-weight:600 !important; font-size:13px !important; padding:7px 22px !important;
      box-shadow:0 2px 8px rgba(3,105,161,.3) !important;
    }
    .ibd-btn-run:hover { opacity:.9; }

    .ibd-panel-boot { background:#faf5ff; border:1px solid #e9d5ff; border-radius:9px; margin-bottom:.85rem; overflow:hidden; }
    .ibd-panel-boot-head { background:#f3e8ff; border-bottom:1px solid #e9d5ff; padding:.55rem .9rem; }
    .ibd-panel-boot-title { font-size:12px; font-weight:600; color:#4c1d95; display:flex; align-items:center; gap:6px; }

    .ibd-boot-result {
      background:#faf5ff; border:1px solid #d8b4fe; border-radius:8px;
      padding:.65rem 1rem; font-size:11.5px; color:#3b0764;
      font-family:'IBM Plex Mono',monospace; line-height:1.9;
      margin-top:.75rem;
    }
    .ibd-boot-result strong { color:#6d28d9; }

    .ibd-dl-row { display:flex; gap:6px; flex-wrap:wrap; margin-top:.5rem; }
    .ibd-dl-row .btn { font-size:11px; padding:3px 12px; }

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
        icon("map-marker-alt"), " Isolation by Distance \u00b7 Rousset (1997)"),
      tags$div(class="ibd-header-sub",
        "FR = F\u209b\u209c/(1\u2212F\u209b\u209c) vs geographic distance \u00b7 Three regression lines ",
        "\u00b7 Mantel test on rectangular matrices"),
      tags$div(class="ibd-badges",
        tags$span(class="ibd-badge ibd-badge-teal",   "FR-ENA \u2014 FreeNA correction"),
        tags$span(class="ibd-badge ibd-badge-blue",   "Mantel test \u2014 permutation"),
        tags$span(class="ibd-badge ibd-badge-amber",  "Model 1D \u2014 distance km"),
        tags$span(class="ibd-badge ibd-badge-purple", "Model 2D \u2014 ln(distance)"),
        tags$span(class="ibd-badge ibd-badge-pink",   "N\u2093 = 1/b")
      )
    ),

    # ── Value boxes ─────────────────────────────────────────────────────────
    tags$div(class="ibd-vbox-row",
      tags$div(class="ibd-vbox",
        tags$div(class="ibd-vbox-icon",style="background:#e0f2fe;color:#0369a1;",icon("dna")),
        tags$div(tags$div(class="ibd-vbox-label","Loci"),
                 tags$div(class="ibd-vbox-val",uiOutput(ns("vb_loci"))))),
      tags$div(class="ibd-vbox",
        tags$div(class="ibd-vbox-icon",style="background:#dcfce7;color:#166534;",icon("map-marker-alt")),
        tags$div(tags$div(class="ibd-vbox-label","Populations"),
                 tags$div(class="ibd-vbox-val",uiOutput(ns("vb_pops"))))),
      tags$div(class="ibd-vbox",
        tags$div(class="ibd-vbox-icon",style="background:#f3e8ff;color:#7e22ce;",icon("project-diagram")),
        tags$div(tags$div(class="ibd-vbox-label","Pairs"),
                 tags$div(class="ibd-vbox-val",uiOutput(ns("vb_pairs"))))),
      tags$div(class="ibd-vbox",
        tags$div(class="ibd-vbox-icon",style="background:#fef9c3;color:#854d0e;",icon("chart-line")),
        tags$div(tags$div(class="ibd-vbox-label","Mantel r"),
                 tags$div(class="ibd-vbox-val",uiOutput(ns("vb_mantel_r"))))),
      tags$div(class="ibd-vbox",
        tags$div(class="ibd-vbox-icon",style="background:#ccfbf1;color:#0d9488;",icon("ruler")),
        tags$div(tags$div(class="ibd-vbox-label","N\u2093"),
                 tags$div(class="ibd-vbox-val",uiOutput(ns("vb_nb")))))
    ),

    # ════════════════════════════════════════════════════════════════════════
    # TABS
    # ════════════════════════════════════════════════════════════════════════
    tabsetPanel(id = ns("ibd_tabs"), type = "tabs",

      # ── TAB 1: Pairwise IBD Analysis ─────────────────────────────────────
      tabPanel(
        title = tagList(icon("project-diagram"), " Pairwise IBD"),
        value = "tab_pairwise", br(),

        # ── Configuration ─────────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("sliders-h"), " IBD parameters")),
          tags$div(class="ibd-panel-body",
            fluidRow(
              column(3,
                radioButtons(
                  ns("model"),
                  "Rousset's model:",
                  choices = c(
                    "2D \u2014 ln(distance km)" = "2D",
                    "1D \u2014 distance km"     = "1D"
                  ),
                  selected = "2D"
                ),
                numericInput(ns("n_boot_pw"), "Bootstrap per pair (CI):",
                             value = 500, min = 100, max = 5000, step = 100),
                numericInput(ns("n_boot_loci"), "Bootstrap over loci:",
                             value = 1000, min = 100, max = 10000, step = 100),
                selectInput(ns("boot_ci_level"),
                  label = "Confidence interval level:",
                  choices = c(
                    "99.99% (alpha = 0.0001)" = "0.0001",
                    "99.9%  (alpha = 0.001)"  = "0.001",
                    "99%    (alpha = 0.01)"   = "0.01",
                    "95%    (alpha = 0.05)"   = "0.05",
                    "90%    (alpha = 0.10)"   = "0.10"
                  ),
                  selected = "0.05"),
                tags$hr(),
                actionButton(
                  ns("run_ibd"), "Run IBD Analysis",
                  icon  = icon("rocket"),
                  class = "ibd-btn-run btn",
                  style = "font-weight:bold; width:100%;"
                ),
                br(), br(),
                uiOutput(ns("ui_run_status"))
              ),
              column(9,
                tags$h5("Regression parameters (Rousset 1997)",
                        style = "font-weight:600; margin-top:0; color:#2c3e50;"),
                DT::DTOutput(ns("reg_table"))
              )
            )
          )
        ),

        # ── IBD plot ───────────────────────────────────────────────────────
        fluidRow(
          column(8,
            tags$div(class="ibd-panel",
              tags$div(class="ibd-panel-head",
                tags$div(class="ibd-panel-title",
                  icon("chart-line"), " IBD plot \u2014 three regression lines")),
              tags$div(class="ibd-panel-body",
                plotly::plotlyOutput(ns("ibd_plot"), height = "460px")
              )
            )
          ),
          column(4,
            tags$div(class="ibd-panel",
              tags$div(class="ibd-panel-head",
                tags$div(class="ibd-panel-title",
                  icon("info-circle"), " Interpretation")),
              tags$div(class="ibd-panel-body", style="font-size:12px; line-height:1.8;",
                tags$p(tags$strong("Three regression lines"), style="margin-bottom:4px;"),
                tags$ul(style="padding-left:16px; line-height:1.9;",
                  tags$li(tags$span(style="color:#333a43; font-weight:600;", "Average"), 
                          " \u2014 through point estimates (FR)"),
                  tags$li(tags$span(style="color:#B40F20; font-weight:600;", "Upper CI (FR-s)"), 
                          " \u2014 through upper 95% CI bounds"),
                  tags$li(tags$span(style="color:#3B9AB2; font-weight:600;", "Lower CI (FR-i)"), 
                          " \u2014 through lower 95% CI bounds")
                ),
                tags$p(tags$strong("Slope b"), 
                       " \u2014 in 2D model: b = 1/N",tags$sub("b")),
                tags$p(tags$strong("N",tags$sub("b")," = 1/b"), 
                       " \u2014 neighbourhood size"),
                tags$p(tags$strong("N",tags$sub("em")," = 1/(2\u03c0b)"), 
                       " \u2014 effective migrants per generation"),
                tags$hr(style="margin:6px 0;"),
                tags$p(style="color:#777; font-size:10.5px;",
                  "Rousset (1997) Genetics 145:1219-1228.")
              )
            )
          )
        ),

        # ── Pairwise tables ────────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("table"), " Pairwise distances with bootstrap confidence intervals")),
          tags$div(class="ibd-panel-body",
            DT::DTOutput(ns("pairwise_table")),
            br(),
            uiOutput(ns("ui_dl_pairwise"))
          )
        )
      ),

      # ── TAB 2: Mantel Test ──────────────────────────────────────────────
      tabPanel(
        title = tagList(icon("chart-bar"), " Mantel Test"),
        value = "tab_mantel", br(),

        tags$div(class="ibd-info",
          icon("info-circle"), " ",
          tags$strong("Mantel test on rectangular matrices"),
          " \u2014 permutation test for correlation between any two distance matrices.",
          " Can test genetic vs geographic, temporal, ecological or categorical distances.",
          tags$br(),
          "Works on rectangular matrices (as in RT or Fstat 2.9.4) where some pairs are excluded."
        ),

        # ── Data input ─────────────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("table"), " Distance data input")),
          tags$div(class="ibd-panel-body",
            fluidRow(
              column(6,
                tags$div(class="ibd-warn",
                  icon("info-circle"), " ",
                  tags$strong("Rectangular matrix format"),
                  tags$br(),
                  "Two first columns: subsample pairs (Pop1, Pop2)",
                  tags$br(),
                  "Subsequent columns: any distance measures",
                  tags$br(),
                  "Rows can be filtered to exclude specific pairs"
                )
              ),
              column(6,
                fileInput(ns("dist_file"), "Upload distance matrix (.csv, .tsv or .txt)",
                          accept = c(".csv", ".tsv", ".txt", ".tab")),
                tags$div(style="margin-top:10px;",
                  actionButton(ns("load_dist"), "Load data",
                               icon = icon("upload"),
                               class = "btn btn-default btn-sm")
                ),
                br(),
                uiOutput(ns("ui_dist_loaded"))
              )
            )
          )
        ),

        # ── Mantel test parameters ────────────────────────────────────────
        tags$div(class="ibd-panel-boot",
          tags$div(class="ibd-panel-boot-head",
            tags$div(class="ibd-panel-boot-title",
              icon("random"), " Mantel test parameters")),
          tags$div(class="ibd-panel-body",
            fluidRow(
              column(4,
                numericInput(ns("n_perm_mantel"), "Number of permutations:",
                             value = 10000, min = 100, max = 99999, step = 1000)
              ),
              column(4,
                selectInput(ns("mantel_x_col"), "Distance X (independent):",
                            choices = c("Select column..." = ""), selected = "")
              ),
              column(4,
                selectInput(ns("mantel_y_col"), "Distance Y (dependent):",
                            choices = c("Select column..." = ""), selected = "")
              )
            ),
            fluidRow(
              column(4,
                radioButtons(ns("mantel_stat"),
                  label = "Test statistic:",
                  choices = c(
                    "Pearson correlation (r)" = "pearson",
                    "Regression slope (b)" = "slope"
                  ),
                  selected = "pearson"
                )
              ),
              column(4,
                radioButtons(ns("mantel_alternative"),
                  label = "Alternative hypothesis:",
                  choices = c(
                    "Positive correlation (r > 0)" = "greater",
                    "Negative correlation (r < 0)" = "less",
                    "Two-sided (r != 0)" = "two.sided"
                  ),
                  selected = "greater"
                )
              ),
              column(4,
                tags$div(style="margin-top:25px;",
                  actionButton(
                    ns("run_mantel"), "Run Mantel Test",
                    icon = icon("play"),
                    class = "ibd-btn-run btn",
                    style = "font-weight:bold;"
                  )
                )
              )
            ),
            uiOutput(ns("ui_mantel_status"))
          )
        ),

        # ── Mantel results ─────────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("chart-line"), " Mantel test results")),
          tags$div(class="ibd-panel-body",
            fluidRow(
              column(6,
                tags$div(class="ibd-boot-result", style="font-size:14px;",
                  tags$p(tags$strong("Observed statistic:"), 
                         style="font-size:16px;", id = ns("mantel_stat_text")),
                  tags$p(tags$strong("p-value:"), 
                         style="font-size:16px;", id = ns("mantel_p_text")),
                  tags$p(tags$strong("Number of permutations:"), 
                         style="font-size:16px;", id = ns("mantel_n_text")),
                  tags$p(tags$strong("Number of pairs:"), 
                         style="font-size:16px;", id = ns("mantel_pairs_text"))
                )
              ),
              column(6,
                tags$div(class="ibd-boot-result", style="font-size:14px;",
                  tags$p(tags$strong("Regression slope (b):"), 
                         style="font-size:16px;", id = ns("mantel_slope_text")),
                  tags$p(tags$strong("Intercept:"), 
                         style="font-size:16px;", id = ns("mantel_intercept_text")),
                  tags$p(tags$strong("R\u00B2:"), 
                         style="font-size:16px;", id = ns("mantel_r2_text")),
                  tags$p(tags$strong("N\u2093 = 1/b:"), 
                         style="font-size:16px;", id = ns("mantel_nb_text"))
                )
              )
            ),
            br(),
            plotly::plotlyOutput(ns("mantel_plot"), height = "400px")
          )
        ),

        # ── Mantel permutation histogram ──────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("chart-bar"), " Permutation distribution")),
          tags$div(class="ibd-panel-body",
            plotly::plotlyOutput(ns("mantel_hist"), height = "350px"),
            br(),
            uiOutput(ns("ui_dl_mantel"))
          )
        )
      )

    ) # end tabsetPanel
  )   # end tags$div.ibd-module
}