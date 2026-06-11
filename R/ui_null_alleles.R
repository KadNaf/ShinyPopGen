# module/ui_null_alleles.R
# Null allele frequency estimation (EM), FST-ENA, DCSE-INA

null_alleles_UI <- function(id) {
  ns <- NS(id)

  custom_css <- tags$style(HTML("
    @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap');

    .na-module * { font-family: 'IBM Plex Sans', sans-serif; }

    .na-header {
      background: linear-gradient(135deg, #0f172a 0%, #1e293b 55%, #0c4a6e 100%);
      border-radius: 10px; padding: 1.2rem 1.6rem; margin-bottom: 1rem;
      position: relative; overflow: hidden;
    }
    .na-header::before {
      content: ''; position: absolute; inset: 0;
      background: repeating-linear-gradient(
        -45deg, transparent, transparent 28px,
        rgba(255,255,255,.018) 28px, rgba(255,255,255,.018) 29px);
    }
    .na-header-title { font-size:1.05rem; font-weight:600; color:#f1f5f9; letter-spacing:.01em; margin-bottom:.2rem; }
    .na-header-sub   { font-size:.75rem; color:#94a3b8; font-family:'IBM Plex Mono',monospace; }
    .na-badges { display:flex; gap:6px; margin-top:.5rem; flex-wrap:wrap; }
    .na-badge  { display:inline-block; border-radius:20px; padding:2px 10px; font-size:.67rem; font-family:'IBM Plex Mono',monospace; }
    .na-badge-blue   { background:rgba(56,189,248,.15);  border:1px solid rgba(56,189,248,.3);  color:#38bdf8; }
    .na-badge-green  { background:rgba(74,222,128,.12);  border:1px solid rgba(74,222,128,.3);  color:#4ade80; }
    .na-badge-amber  { background:rgba(251,191,36,.12);  border:1px solid rgba(251,191,36,.3);  color:#fbbf24; }
    .na-badge-teal   { background:rgba(20,184,166,.15);  border:1px solid rgba(20,184,166,.3);  color:#2dd4bf; }

    .na-vbox-row { display:flex; gap:9px; margin-bottom:1rem; flex-wrap:wrap; }
    .na-vbox { flex:1; min-width:110px; background:#fff; border:1px solid #e2e8f0; border-radius:9px; padding:.6rem .85rem; display:flex; align-items:center; gap:9px; }
    .na-vbox-icon  { width:30px; height:30px; border-radius:7px; display:flex; align-items:center; justify-content:center; font-size:12px; flex-shrink:0; }
    .na-vbox-label { font-size:10px; color:#94a3b8; text-transform:uppercase; letter-spacing:.06em; margin-bottom:1px; }
    .na-vbox-val   { font-size:18px; font-weight:600; color:#0f172a; line-height:1.1; font-family:'IBM Plex Mono',monospace; }

    .na-panel { background:#fff; border:1px solid #e2e8f0; border-radius:9px; margin-bottom:.85rem; overflow:hidden; }
    .na-panel-head { background:#f8fafc; border-bottom:1px solid #e2e8f0; padding:.55rem .9rem; }
    .na-panel-title { font-size:12px; font-weight:600; color:#1e293b; display:flex; align-items:center; gap:6px; flex-wrap:wrap; }
    .na-panel-body { padding:.85rem; }

    .na-info { background:#eff6ff; border:1px solid #bfdbfe; border-radius:7px; padding:.45rem .8rem; font-size:11.5px; color:#1d4ed8; margin-bottom:.85rem; line-height:1.65; }
    .na-warn { background:#fffbeb; border:1px solid #fcd34d; border-radius:7px; padding:.45rem .8rem; font-size:11.5px; color:#92400e; margin-bottom:.85rem; line-height:1.65; }

    .na-locus-grid { display:flex; flex-wrap:wrap; gap:8px; margin-top:.5rem; }
    .na-locus-item {
      background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px;
      padding:.45rem .7rem; min-width:160px; flex:1;
    }
    .na-locus-item .control-label { display:none; }
    .na-locus-name {
      font-size:11px; font-weight:700; color:#1e293b;
      font-family:'IBM Plex Mono',monospace; margin-bottom:3px;
    }
    .na-locus-item .radio { margin:2px 0; }
    .na-locus-item .radio label { font-size:11px; color:#475569; }

    .na-btn-run {
      background:linear-gradient(135deg,#0369a1,#0c4a6e) !important;
      border:none !important; color:#fff !important; border-radius:7px !important;
      font-weight:600 !important; font-size:13px !important; padding:7px 22px !important;
      box-shadow:0 2px 8px rgba(3,105,161,.3) !important;
    }
    .na-btn-run:hover { opacity:.9; }

    .na-boot-result {
      background:#faf5ff; border:1px solid #d8b4fe; border-radius:8px;
      padding:.65rem 1rem; font-size:11.5px; color:#3b0764;
      font-family:'IBM Plex Mono',monospace; line-height:1.9;
      margin-top:.75rem;
    }
    .na-boot-result strong { color:#6d28d9; }

    .na-matrix-wrap { overflow-x:auto; margin-top:.5rem; }
    .na-matrix { border-collapse:collapse; font-size:11px; font-family:'IBM Plex Mono',monospace; width:100%; }
    .na-matrix th { background:#f8fafc; color:#475569; font-weight:600; padding:4px 9px; border:1px solid #e2e8f0; font-size:10.5px; white-space:nowrap; }
    .na-matrix td { padding:4px 9px; border:1px solid #e2e8f0; color:#1e293b; text-align:right; white-space:nowrap; font-size:11px; }
    .na-matrix tr:nth-child(even) td { background:#f8fafc; }
    .na-matrix .diag  { background:#f1f5f9 !important; color:#94a3b8; text-align:center; }
    .na-matrix .upper { color:#cbd5e1; text-align:center; }
    .na-matrix .lbl   { font-weight:700; color:#0f172a; text-align:left; white-space:nowrap; }

    .na-dl-row { display:flex; gap:6px; flex-wrap:wrap; margin-top:.5rem; }
    .na-dl-row .btn { font-size:11px; padding:3px 12px; }

    .na-module .dataTables_wrapper { font-size:12px; }
    .na-module table.dataTable thead th {
      background:#f8fafc !important; color:#475569 !important;
      font-family:'IBM Plex Mono',monospace !important;
      font-size:10.5px !important; font-weight:600 !important;
    }
    .na-module table.dataTable tbody td {
      font-family:'IBM Plex Mono',monospace !important;
      font-size:11px !important; color:#1e293b !important;
    }
    .na-module .nav-tabs > li > a { font-size:12px; font-weight:500; color:#475569; padding:5px 13px; }
    .na-module .nav-tabs > li.active > a { color:#0f172a; font-weight:600; }
  "))

  tags$div(class="na-module", custom_css,

    tags$div(class="na-header",
      tags$div(class="na-header-title",
        icon("atom"), " Null Allele Estimation \u00b7 FST-ENA \u00b7 DCSE-INA"),
      tags$div(class="na-header-sub",
        "EM algorithm \u00b7 Dempster, Laird & Rubin (1977) \u00b7 FreeNA \u2014 Chapuis & Estoup (2007)",
        " \u00b7 Weir (1996) \u00b7 Cavalli-Sforza & Edwards (1967)"),
      tags$div(class="na-badges",
        tags$span(class="na-badge na-badge-blue",  "EM \u2014 null allele frequency"),
        tags$span(class="na-badge na-badge-teal",  "ENA \u2014 FST corrected"),
        tags$span(class="na-badge na-badge-green", "INA \u2014 DCSE corrected"),
        tags$span(class="na-badge na-badge-amber", "Bootstrap CI \u2014 loci & sub-samples")
      )
    ),

    tags$div(class="na-vbox-row",
      tags$div(class="na-vbox",
        tags$div(class="na-vbox-icon",style="background:#e0f2fe;color:#0369a1;",icon("dna")),
        tags$div(tags$div(class="na-vbox-label","Loci"),
                 tags$div(class="na-vbox-val",uiOutput(ns("vb_loci"))))),
      tags$div(class="na-vbox",
        tags$div(class="na-vbox-icon",style="background:#dcfce7;color:#166534;",icon("map-marker-alt")),
        tags$div(tags$div(class="na-vbox-label","Populations"),
                 tags$div(class="na-vbox-val",uiOutput(ns("vb_pops"))))),
      tags$div(class="na-vbox",
        tags$div(class="na-vbox-icon",style="background:#f3e8ff;color:#7e22ce;",icon("users")),
        tags$div(tags$div(class="na-vbox-label","Individuals"),
                 tags$div(class="na-vbox-val",uiOutput(ns("vb_n"))))),
      tags$div(class="na-vbox",
        tags$div(class="na-vbox-icon",style="background:#fef9c3;color:#854d0e;",icon("percentage")),
        tags$div(tags$div(class="na-vbox-label","Avg p_nulls"),
                 tags$div(class="na-vbox-val",uiOutput(ns("vb_avg_null"))))),
      tags$div(class="na-vbox",
        tags$div(class="na-vbox-icon",style="background:#fce7f3;color:#9d174d;",icon("exclamation-triangle")),
        tags$div(tags$div(class="na-vbox-label","Max p_nulls"),
                 tags$div(class="na-vbox-val",uiOutput(ns("vb_max_null"))))),
      tags$div(class="na-vbox",
        tags$div(class="na-vbox-icon",style="background:#ccfbf1;color:#0d9488;",icon("chart-bar")),
        tags$div(tags$div(class="na-vbox-label","Global FST-ENA"),
                 tags$div(class="na-vbox-val",uiOutput(ns("vb_fst_ena")))))
    ),

    tags$div(class="na-panel",
      tags$div(class="na-panel-head",
        tags$div(class="na-panel-title",
          icon("sliders-h"), " Setup \u2014 3 parameters to configure")),
      tags$div(class="na-panel-body",

        tags$div(class="na-warn",
          icon("exclamation-triangle"), " ",
          tags$strong("(1) Missing genotype coding per locus"),
          tags$br(),
          tags$span(style="font-size:11px;",
            tags$strong("000000"), " \u2014 missing coded as absent / PCR failure (recommended default).",
            tags$br(),
            tags$strong("999999"), " \u2014 missing coded as null homozygote."
          )
        ),
        uiOutput(ns("locus_coding_ui")),

        tags$hr(style="margin:1rem 0;"),

        tags$strong("(2) Bootstrap parameters", style="font-size:12px; color:#1e293b;"),
        tags$br(), tags$br(),
        fluidRow(
          column(4,
            numericInput(ns("nboot"),
              label = "Number of replicates:",
              value = 5000, min = 100, max = 99999, step = 1000)),
          column(4,
            selectInput(ns("ci_level"),
              label = "Confidence interval level:",
              choices = c(
                "99.99% (alpha = 0.0001)" = "0.0001",
                "99.9%  (alpha = 0.001)"  = "0.001",
                "99%    (alpha = 0.01)"   = "0.01",
                "95%    (alpha = 0.05)"   = "0.05",
                "90%    (alpha = 0.10)"   = "0.10"
              ),
              selected = "0.05")),
          column(4,
            tags$div(style="margin-top:25px;font-size:11px;color:#64748b;",
              icon("info-circle"), " Bootstrap over loci: vectorised, fast.",
              tags$br(),
              "Bootstrap over sub-samples: re-runs EM per replicate."
            ))
        ),

        tags$hr(style="margin:1rem 0;"),

        tags$strong("(3) Run all computations + generate output files",
                    style="font-size:12px; color:#1e293b;"),
        tags$br(), tags$br(),
        fluidRow(
          column(4,
            actionButton(ns("run_all"),
              label = tagList(icon("play"), tags$strong("  Compute + Bootstrap + Export")),
              class = "na-btn-run btn",
              width = "100%"))
        ),
        br(),
        uiOutput(ns("ui_run_status"))
      )
    ),

    tags$div(class="na-panel",
      tags$div(class="na-panel-head",
        tags$div(class="na-panel-title",
          icon("file-download"), " Output files \u2014 automatically generated after computation")),
      tags$div(class="na-panel-body",
        tags$div(class="na-info",
          icon("info-circle"), " ",
          "All four files are generated automatically when you click Compute above."
        ),
        fluidRow(
          column(3,
            tags$div(class="na-panel", style="border-color:#bfdbfe;",
              tags$div(class="na-panel-head", style="background:#eff6ff;",
                tags$div(class="na-panel-title", style="color:#1d4ed8;",
                  icon("file-alt"), " File 1 \u2014 Null allele frequencies")),
              tags$div(class="na-panel-body", style="font-size:11px;color:#334155;",
                "p_nulls per locus \u00d7 population",
                tags$br(), "Global weighted mean per locus",
                tags$br(), "Locus coding reminder",
                tags$br(), br(),
                uiOutput(ns("ui_dl_file1"))
              )
            )
          ),
          column(3,
            tags$div(class="na-panel", style="border-color:#99f6e4;",
              tags$div(class="na-panel-head", style="background:#f0fdfa;",
                tags$div(class="na-panel-title", style="color:#0d9488;",
                  icon("chart-bar"), " File 2 \u2014 Global FST & FST-ENA")),
              tags$div(class="na-panel-body", style="font-size:11px;color:#334155;",
                "Per locus + multilocus FST / FST-ENA",
                tags$br(), "CI from bootstrap over loci",
                tags$br(), "CI from bootstrap over sub-samples",
                uiOutput(ns("ui_dl_file2"))
              )
            )
          ),
          column(3,
            tags$div(class="na-panel", style="border-color:#e9d5ff;",
              tags$div(class="na-panel-head", style="background:#faf5ff;",
                tags$div(class="na-panel-title", style="color:#7c3aed;",
                  icon("table"), " File 3 \u2014 Pairwise long format")),
              tags$div(class="na-panel-body", style="font-size:11px;color:#334155;",
                "FST, FST-ENA, DCSE, DCSE-INA",
                tags$br(), "Per pair of sub-samples",
                tags$br(), "CI from bootstrap over loci",
                uiOutput(ns("ui_dl_file3"))
              )
            )
          ),
          column(3,
            tags$div(class="na-panel", style="border-color:#fcd34d;",
              tags$div(class="na-panel-head", style="background:#fffbeb;",
                tags$div(class="na-panel-title", style="color:#92400e;",
                  icon("th"), " File 4 \u2014 Per-locus half-matrices")),
              tags$div(class="na-panel-body", style="font-size:11px;color:#334155;",
                "FST, FST-ENA, DCSE, DCSE-INA",
                tags$br(), "Half-matrix per locus",
                tags$br(), uiOutput(ns("ui_dl_file4"))
              )
            )
          )
        )
      )
    ),

    tabsetPanel(id = ns("na_tabs"), type = "tabs",

      tabPanel(title = tagList(icon("dna"), " Null allele frequencies"),
               value = "tab_na", br(),
        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              icon("list"), " p_nulls per locus \u00d7 population (EM algorithm)")),
          tags$div(class="na-panel-body",
            DT::DTOutput(ns("dt_t1")))),
        tags$br(),
        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              icon("globe"), " Global summary per locus (N-weighted mean)")),
          tags$div(class="na-panel-body",
            DT::DTOutput(ns("dt_t2"))))
      ),

      tabPanel(title = tagList(icon("chart-bar"), " FST / FST-ENA"),
               value = "tab_fst", br(),

        tags$div(class="na-info",
          icon("info-circle"), " ",
          tags$strong("Global multilocus FST"), " \u2014 Weir (1996) / Genepop method. ",
          tags$strong("FST-ENA"), ": EM-corrected frequencies, Excluding Null Alleles."
        ),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              icon("list"), " Per-locus FST and FST-ENA")),
          tags$div(class="na-panel-body",
            DT::DTOutput(ns("dt_fst_global")))),
        br(),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              icon("random"), " Bootstrap CI \u2014 Global FST and FST-ENA")),
          tags$div(class="na-panel-body",
            uiOutput(ns("ui_boot_global_fst")))),
        br(),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              icon("exchange-alt"), " Pairwise FST and FST-ENA")),
          tags$div(class="na-panel-body",
            fluidRow(
              column(5,
                radioButtons(ns("fst_pair_display"), "Display:",
                  choices = c("Raw FST" = "raw", "FST-ENA" = "ena", "Both" = "both"),
                  selected = "both", inline = TRUE))),
            uiOutput(ns("ui_fst_pair_matrix")))),
        br(),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              icon("random"), " Bootstrap CI \u2014 Pairwise FST-ENA")),
          tags$div(class="na-panel-body",
            uiOutput(ns("ui_boot_pair_fst"))))
      ),

      tabPanel(title = tagList(icon("ruler-combined"), " DCSE / DCSE-INA"),
               value = "tab_dc", br(),

        tags$div(class="na-info",
          icon("info-circle"), " ",
          tags$strong("Cavalli-Sforza & Edwards (1967) chord distance."),
          " DCSE-INA includes the null allele as an extra state."
        ),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              icon("th"), " Pairwise DCSE and DCSE-INA")),
          tags$div(class="na-panel-body",
            fluidRow(
              column(5,
                radioButtons(ns("dc_display"), "Display:",
                  choices = c("Raw DCSE" = "raw", "DCSE-INA" = "ina", "Both" = "both"),
                  selected = "both", inline = TRUE))),
            uiOutput(ns("ui_dc_matrix")))),
        br(),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              icon("random"), " Bootstrap CI \u2014 Pairwise DCSE-INA")),
          tags$div(class="na-panel-body",
            uiOutput(ns("ui_boot_pair_dc"))))
      ),

      tabPanel(title = tagList(icon("table"), " Per-locus \u00d7 pair"),
               value = "tab_locus_pair", br(),

        tags$div(class="na-info",
          icon("info-circle"), " ",
          "FST, FST-ENA, DCSE and DCSE-INA for each locus \u00d7 pair of populations."
        ),

        fluidRow(
          column(3, selectInput(ns("fl_locus"), "Locus:",
            choices = c("All loci" = "all"), selected = "all")),
          column(3, selectInput(ns("fl_pop1"), "Population 1:",
            choices = c("All pairs" = "all"), selected = "all")),
          column(3, selectInput(ns("fl_pop2"), "Population 2:",
            choices = c("All pairs" = "all"), selected = "all"))
        ),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              icon("list"), " FST and FST-ENA per locus \u00d7 pair")),
          tags$div(class="na-panel-body",
            DT::DTOutput(ns("dt_fst_locus")))),
        br(),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              icon("list"), " DCSE and DCSE-INA per locus \u00d7 pair")),
          tags$div(class="na-panel-body",
            DT::DTOutput(ns("dt_dc_locus"))))
      )

    )
  )
}