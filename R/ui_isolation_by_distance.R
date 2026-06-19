# ui_isolation_by_distance.R
# Tab: Isolation by Distance
# Rousset (1997) linearised FST/(1-FST) vs geographic distance, Mantel test.
# Three regression lines: average, upper 95% CI, lower 95% CI.
# Two tabs: IBD Analysis and Bootstrap Confidence Intervals

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
    .ibd-badge-teal   { background:rgba(20,184,166,.15);  border:1px solid rgba(20,184,166,.3);  color:#2dd4bf; }
    .ibd-badge-blue   { background:rgba(56,189,248,.15);  border:1px solid rgba(56,189,248,.3);  color:#38bdf8; }
    .ibd-badge-amber  { background:rgba(251,191,36,.12);  border:1px solid rgba(251,191,36,.3);  color:#fbbf24; }
    .ibd-badge-purple { background:rgba(192,132,252,.12); border:1px solid rgba(192,132,252,.3); color:#c084fc; }

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
    .ibd-info { background:#eff6ff; border:1px solid #bfdbfe; border-radius:7px; padding:.45rem .8rem; font-size:11.5px; color:#1d4ed8; margin-bottom:.85rem; line-height:1.65; }
    .ibd-warn { background:#fffbeb; border:1px solid #fcd34d; border-radius:7px; padding:.45rem .8rem; font-size:11.5px; color:#92400e; margin-bottom:.85rem; line-height:1.65; }

    /* ── Buttons ─────────────────────────────────────────────────────── */
    .ibd-btn-run {
      background:linear-gradient(135deg,#0369a1,#0c4a6e) !important;
      border:none !important; color:#fff !important; border-radius:7px !important;
      font-weight:600 !important; font-size:13px !important; padding:7px 22px !important;
      box-shadow:0 2px 8px rgba(3,105,161,.3) !important;
    }
    .ibd-btn-run:hover { opacity:.9; }

    /* ── Bootstrap panel ─────────────────────────────────────────────── */
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

    /* ── Matrix table ────────────────────────────────────────────────── */
    .ibd-matrix-wrap { overflow-x:auto; margin-top:.5rem; }
    .ibd-matrix { border-collapse:collapse; font-size:11px; font-family:'IBM Plex Mono',monospace; width:100%; }
    .ibd-matrix th { background:#f8fafc; color:#475569; font-weight:600; padding:4px 9px; border:1px solid #e2e8f0; font-size:10.5px; white-space:nowrap; }
    .ibd-matrix td { padding:4px 9px; border:1px solid #e2e8f0; color:#1e293b; text-align:right; white-space:nowrap; font-size:11px; }
    .ibd-matrix tr:nth-child(even) td { background:#f8fafc; }
    .ibd-matrix .diag  { background:#f1f5f9 !important; color:#94a3b8; text-align:center; }
    .ibd-matrix .upper { color:#cbd5e1; text-align:center; }
    .ibd-matrix .lbl   { font-weight:700; color:#0f172a; text-align:left; white-space:nowrap; }

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

  # ── Shared download row ──────────────────────────────────────────────────
  dlrow <- function(...) tags$div(class="ibd-dl-row", ...)

  tags$div(class="ibd-module", custom_css,

    # ── Header ─────────────────────────────────────────────────────────────
    tags$div(class="ibd-header",
      tags$div(class="ibd-header-title",
        icon("map-marker-alt"), " Isolation by Distance \u00b7 Rousset (1997)"),
      tags$div(class="ibd-header-sub",
        "Linearised F\u209b\u209c/(1\u2212F\u209b\u209c) vs geographic distance \u00b7 Mantel test ",
        "\u00b7 Three regression lines (average, upper CI, lower CI)"),
      tags$div(class="ibd-badges",
        tags$span(class="ibd-badge ibd-badge-teal",   "F\u209b\u209c-ENA \u2014 FreeNA correction"),
        tags$span(class="ibd-badge ibd-badge-blue",   "Mantel test \u2014 permutation"),
        tags$span(class="ibd-badge ibd-badge-amber",  "1D \u2014 distance km"),
        tags$span(class="ibd-badge ibd-badge-purple", "2D \u2014 ln(distance)"),
        tags$span(class="ibd-badge ibd-badge-teal",   "N\u2093 = 1/b")
      )
    ),

    # ── Value boxes ─────────────────────────────────────────────────────────
    tags$div(class="ibd-vbox-row",
      tags$div(class="ibd-vbox",
        tags$div(class="ibd-vbox-icon",style="background:#e0f2fe;color:#0369a1;",icon("map-marker-alt")),
        tags$div(tags$div(class="ibd-vbox-label","Populations"),
                 tags$div(class="ibd-vbox-val",uiOutput(ns("vb_pops"))))),
      tags$div(class="ibd-vbox",
        tags$div(class="ibd-vbox-icon",style="background:#dcfce7;color:#166534;",icon("project-diagram")),
        tags$div(tags$div(class="ibd-vbox-label","Pairs"),
                 tags$div(class="ibd-vbox-val",uiOutput(ns("vb_pairs"))))),
      tags$div(class="ibd-vbox",
        tags$div(class="ibd-vbox-icon",style="background:#f3e8ff;color:#7e22ce;",icon("chart-line")),
        tags$div(tags$div(class="ibd-vbox-label","Mantel r"),
                 tags$div(class="ibd-vbox-val",uiOutput(ns("vb_mantel_r"))))),
      tags$div(class="ibd-vbox",
        tags$div(class="ibd-vbox-icon",style="background:#fef9c3;color:#854d0e;",icon("check-circle")),
        tags$div(tags$div(class="ibd-vbox-label","p-value"),
                 tags$div(class="ibd-vbox-val",uiOutput(ns("vb_pval"))))),
      tags$div(class="ibd-vbox",
        tags$div(class="ibd-vbox-icon",style="background:#ccfbf1;color:#0d9488;",icon("ruler")),
        tags$div(tags$div(class="ibd-vbox-label","N\u2093 (avg)"),
                 tags$div(class="ibd-vbox-val",uiOutput(ns("vb_nb")))))
    ),

    # ════════════════════════════════════════════════════════════════════════
    # TABS
    # ════════════════════════════════════════════════════════════════════════
    tabsetPanel(id = ns("ibd_tabs"), type = "tabs",

      # ── TAB 1: IBD Analysis ──────────────────────────────────────────────
      tabPanel(
        title = tagList(icon("chart-line"), " IBD Analysis"),
        value = "tab_ibd", br(),

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
                  "Habitat model:",
                  choices = c(
                    "2D \u2014 ln(distance km)" = "2D",
                    "1D \u2014 distance km"     = "1D"
                  ),
                  selected = "2D"
                ),
                numericInput(ns("n_boot_pw"), "Bootstrap per pair (CI):",
                             value = 500, min = 100, max = 5000, step = 100),
                numericInput(ns("n_perm"), "Mantel permutations:",
                             value = 9999, min = 99, max = 99999, step = 1000),
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
                tags$h5("Regression parameters",
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
                          " \u2014 through point estimates"),
                  tags$li(tags$span(style="color:#B40F20; font-weight:600;", "Upper CI (ls)"), 
                          " \u2014 through upper 95% CI bounds"),
                  tags$li(tags$span(style="color:#3B9AB2; font-weight:600;", "Lower CI (li)"), 
                          " \u2014 through lower 95% CI bounds")
                ),
                tags$p(tags$strong("Slope b"), 
                       " \u2014 in 2D model: b = 1/N",tags$sub("b")),
                tags$p(tags$strong("N",tags$sub("b")," = 1/b"), 
                       " \u2014 neighbourhood size"),
                tags$p(tags$strong("N",tags$sub("em")," = 1/(2\u03c0b)"), 
                       " \u2014 effective migrants"),
                tags$hr(style="margin:6px 0;"),
                tags$p(style="color:#777; font-size:10.5px;",
                  "Rousset (1997) Genetics 145:1219. de Mee\u00fbs (2006).")
              )
            )
          )
        ),

        # ── Pairwise tables ────────────────────────────────────────────────
        fluidRow(
          column(7,
            tags$div(class="ibd-panel",
              tags$div(class="ibd-panel-head",
                tags$div(class="ibd-panel-title",
                  icon("table"), " Pairwise F\u209b\u209c-ENA & linearised values")),
              tags$div(class="ibd-panel-body",
                DT::DTOutput(ns("fst_table")),
                br(),
                uiOutput(ns("ui_dl_fst"))
              )
            )
          ),
          column(5,
            tags$div(class="ibd-panel",
              tags$div(class="ibd-panel-head",
                tags$div(class="ibd-panel-title",
                  icon("ruler"), " Pairwise distances (km)")),
              tags$div(class="ibd-panel-body",
                DT::DTOutput(ns("dist_table")),
                br(),
                uiOutput(ns("ui_dl_dist"))
              )
            )
          )
        )
      ),

      # ── TAB 2: Bootstrap Confidence Intervals ──────────────────────────
      tabPanel(
        title = tagList(icon("braces"), " Bootstrap CIs"),
        value = "tab_boot", br(),

        tags$div(class="ibd-info",
          icon("info-circle"), " ",
          tags$strong("Bootstrap over loci"), 
          " \u2014 resample loci with replacement (FreeNA approach) to compute 95% CI ",
          "for pairwise F<sub>ST</sub>-ENA and linearised values. ",
          "Chapuis & Estoup (2007) / FreeNA method."
        ),

        # ── Bootstrap parameters ──────────────────────────────────────────
        tags$div(class="ibd-panel-boot",
          tags$div(class="ibd-panel-boot-head",
            tags$div(class="ibd-panel-boot-title",
              icon("random"), " Bootstrap parameters")),
          tags$div(class="ibd-panel-body",
            fluidRow(
              column(4,
                numericInput(ns("n_boot_loci"), "Number of bootstrap replicates:",
                             value = 1000, min = 100, max = 10000, step = 100)
              ),
              column(4,
                selectInput(ns("boot_ci_level"),
                  label = "Confidence interval level:",
                  choices = c(
                    "99.99% (alpha = 0.0001)" = "0.0001",
                    "99.9%  (alpha = 0.001)"  = "0.001",
                    "99%    (alpha = 0.01)"   = "0.01",
                    "95%    (alpha = 0.05)"   = "0.05",
                    "90%    (alpha = 0.10)"   = "0.10"
                  ),
                  selected = "0.05")
              ),
              column(4,
                tags$div(style="margin-top:25px;",
                  actionButton(
                    ns("run_boot"), "Run Bootstrap",
                    icon = icon("play"),
                    class = "ibd-btn-run btn",
                    style = "font-weight:bold;"
                  )
                )
              )
            ),
            uiOutput(ns("ui_boot_status"))
          )
        ),

        # ── Bootstrap summary boxes ──────────────────────────────────────
        tags$div(class="ibd-vbox-row",
          tags$div(class="ibd-vbox",
            tags$div(class="ibd-vbox-icon",style="background:#e0f2fe;color:#0369a1;",icon("dna")),
            tags$div(tags$div(class="ibd-vbox-label","Loci"),
                     tags$div(class="ibd-vbox-val",uiOutput(ns("boot_n_loci"))))),
          tags$div(class="ibd-vbox",
            tags$div(class="ibd-vbox-icon",style="background:#f3e8ff;color:#7e22ce;",icon("repeat")),
            tags$div(tags$div(class="ibd-vbox-label","Replicates"),
                     tags$div(class="ibd-vbox-val",uiOutput(ns("boot_n_reps"))))),
          tags$div(class="ibd-vbox",
            tags$div(class="ibd-vbox-icon",style="background:#dcfce7;color:#166534;",icon("check-circle")),
            tags$div(tags$div(class="ibd-vbox-label","Valid pairs"),
                     tags$div(class="ibd-vbox-val",uiOutput(ns("boot_n_valid")))))
        ),

        # ── Bootstrap table ──────────────────────────────────────────────
        tags$div(class="ibd-panel",
          tags$div(class="ibd-panel-head",
            tags$div(class="ibd-panel-title",
              icon("table"), " Bootstrap CI results")),
          tags$div(class="ibd-panel-body",
            DT::DTOutput(ns("boot_table")),
            br(),
            uiOutput(ns("ui_dl_boot"))
          )
        ),

        # ── Bootstrap plots ──────────────────────────────────────────────
        fluidRow(
          column(6,
            tags$div(class="ibd-panel",
              tags$div(class="ibd-panel-head",
                tags$div(class="ibd-panel-title",
                  icon("chart-line"), " FST-ENA CI plot")),
              tags$div(class="ibd-panel-body",
                plotly::plotlyOutput(ns("boot_fst_plot"), height = "400px")
              )
            )
          ),
          column(6,
            tags$div(class="ibd-panel",
              tags$div(class="ibd-panel-head",
                tags$div(class="ibd-panel-title",
                  icon("chart-line"), " FR (linearised) CI plot")),
              tags$div(class="ibd-panel-body",
                plotly::plotlyOutput(ns("boot_fr_plot"), height = "400px")
              )
            )
          )
        )
      )

    ) # end tabsetPanel
  )   # end tags$div.ibd-module
}