# module/ui_null_alleles.R
# Null allele frequency estimation (EM), FST-ENA, DCSE-INA
# Simplified UI per supervisor feedback:
#   - Radio buttons for missing genotype coding, default = 000000
#   - Single bootstrap panel: n replicates + CI level choice
#   - 4 automatic output files
#
# References:
#   Dempster, Laird & Rubin (1977)  — EM algorithm
#   Chapuis & Estoup (2007)         — FreeNA: ENA and INA corrections
#   Weir (1996)                     — FST following Genepop method
#   Cavalli-Sforza & Edwards (1967) — Chord genetic distance (DCSE)

null_alleles_UI <- function(id) {
  ns <- NS(id)

  custom_css <- tags$style(HTML("
    

    .na-module * { font-family: 'Segoe UI','IBM Plex Sans',Arial,sans-serif; }

    /* ── Header ─────────────────────────────────────────────────────── */
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
    .na-header-sub   { font-size:.75rem; color:#94a3b8; font-family:'Consolas','IBM Plex Mono',monospace; }
    .na-badges { display:flex; gap:6px; margin-top:.5rem; flex-wrap:wrap; }
    .na-badge  { display:inline-block; border-radius:20px; padding:2px 10px; font-size:.67rem; font-family:'Consolas','IBM Plex Mono',monospace; }
    .na-badge-blue   { background:rgba(56,189,248,.15);  border:1px solid rgba(56,189,248,.3);  color:#38bdf8; }
    .na-badge-green  { background:rgba(74,222,128,.12);  border:1px solid rgba(74,222,128,.3);  color:#4ade80; }
    .na-badge-amber  { background:rgba(251,191,36,.12);  border:1px solid rgba(251,191,36,.3);  color:#fbbf24; }
    .na-badge-teal   { background:rgba(20,184,166,.15);  border:1px solid rgba(20,184,166,.3);  color:#2dd4bf; }

    /* ── Value boxes ─────────────────────────────────────────────────── */
    .na-vbox-row { display:flex; gap:9px; margin-bottom:1rem; flex-wrap:wrap; }
    .na-vbox { flex:1; min-width:110px; background:#fff; border:1px solid #e2e8f0; border-radius:9px; padding:.6rem .85rem; display:flex; align-items:center; gap:9px; }
    .na-vbox-icon  { width:30px; height:30px; border-radius:7px; display:flex; align-items:center; justify-content:center; font-size:12px; flex-shrink:0; }
    .na-vbox-label { font-size:10px; color:#94a3b8; text-transform:uppercase; letter-spacing:.06em; margin-bottom:1px; }
    .na-vbox-val   { font-size:18px; font-weight:600; color:#0f172a; line-height:1.1; font-family:'Consolas','IBM Plex Mono',monospace; }

    /* ── Panels ──────────────────────────────────────────────────────── */
    .na-panel { background:#fff; border:1px solid #e2e8f0; border-radius:9px; margin-bottom:.85rem; overflow:hidden; }
    .na-panel-head { background:#f8fafc; border-bottom:1px solid #e2e8f0; padding:.55rem .9rem; }
    .na-panel-title { font-size:12px; font-weight:600; color:#1e293b; display:flex; align-items:center; gap:6px; flex-wrap:wrap; }
    .na-panel-body { padding:.85rem; }

    /* ── Bootstrap panel ─────────────────────────────────────────────── */
    .na-panel-boot { background:#faf5ff; border:1px solid #e9d5ff; border-radius:9px; margin-bottom:.85rem; overflow:hidden; }
    .na-panel-boot-head { background:#f3e8ff; border-bottom:1px solid #e9d5ff; padding:.55rem .9rem; }
    .na-panel-boot-title { font-size:12px; font-weight:600; color:#4c1d95; display:flex; align-items:center; gap:6px; }

    /* ── Info strips ─────────────────────────────────────────────────── */
    .na-info { background:#eff6ff; border:1px solid #bfdbfe; border-radius:7px; padding:.45rem .8rem; font-size:11.5px; color:#1d4ed8; margin-bottom:.85rem; line-height:1.65; }
    .na-warn { background:#fffbeb; border:1px solid #fcd34d; border-radius:7px; padding:.45rem .8rem; font-size:11.5px; color:#92400e; margin-bottom:.85rem; line-height:1.65; }

    /* ── Locus coding grid — radio buttons ───────────────────────────── */
    .na-locus-grid { display:flex; flex-wrap:wrap; gap:8px; margin-top:.5rem; }
    .na-locus-item {
      background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px;
      padding:.5rem .8rem; min-width:180px; flex:1;
    }
    .na-locus-item .control-label { display:none; } /* hide redundant label */
    .na-locus-name {
      font-size:22px; font-weight:800; color:#0f172a;
      font-family:'Consolas','IBM Plex Mono',monospace; margin-bottom:5px;
    }
    .na-locus-item .radio { margin:2px 0; }
    .na-locus-item .radio label { font-size:13px; color:#334155; }

    /* ── Buttons ─────────────────────────────────────────────────────── */
    .na-btn-run {
      background:linear-gradient(135deg,#0369a1,#0c4a6e) !important;
      border:none !important; color:#fff !important; border-radius:7px !important;
      font-weight:600 !important; font-size:13px !important; padding:7px 22px !important;
      box-shadow:0 2px 8px rgba(3,105,161,.3) !important;
    }
    .na-btn-run:hover { opacity:.9; }
    .na-btn-boot {
      background:linear-gradient(135deg,#7c3aed,#4c1d95) !important;
      border:none !important; color:#fff !important; border-radius:7px !important;
      font-weight:600 !important; font-size:13px !important; padding:7px 22px !important;
      box-shadow:0 2px 8px rgba(124,58,237,.3) !important;
    }
    .na-btn-boot:hover { opacity:.9; }

    /* ── Bootstrap result ────────────────────────────────────────────── */
    .na-boot-result {
      background:#faf5ff; border:1px solid #d8b4fe; border-radius:8px;
      padding:.65rem 1rem; font-size:11.5px; color:#3b0764;
      font-family:'Consolas','IBM Plex Mono',monospace; line-height:1.9;
      margin-top:.75rem;
    }
    .na-boot-result strong { color:#6d28d9; }

    /* ── Matrix table ────────────────────────────────────────────────── */
    .na-matrix-wrap { overflow-x:auto; margin-top:.5rem; }
    .na-matrix { border-collapse:collapse; font-size:11px; font-family:'Consolas','IBM Plex Mono',monospace; width:100%; }
    .na-matrix th { background:#f8fafc; color:#475569; font-weight:600; padding:4px 9px; border:1px solid #e2e8f0; font-size:10.5px; white-space:nowrap; }
    .na-matrix td { padding:4px 9px; border:1px solid #e2e8f0; color:#1e293b; text-align:right; white-space:nowrap; font-size:11px; }
    .na-matrix tr:nth-child(even) td { background:#f8fafc; }
    .na-matrix .diag  { background:#f1f5f9 !important; color:#94a3b8; text-align:center; }
    .na-matrix .upper { color:#cbd5e1; text-align:center; }
    .na-matrix .lbl   { font-weight:700; color:#0f172a; text-align:left; white-space:nowrap; }

    /* ── Download row ────────────────────────────────────────────────── */
    .na-dl-row { display:flex; gap:6px; flex-wrap:wrap; margin-top:.5rem; }
    .na-dl-row .btn { font-size:11px; padding:3px 12px; }

    /* ── DT tweaks ───────────────────────────────────────────────────── */
    .na-module .dataTables_wrapper { font-size:12px; }
    .na-module table.dataTable thead th {
      background:#f8fafc !important; color:#475569 !important;
      font-family:'Consolas','IBM Plex Mono',monospace !important;
      font-size:10.5px !important; font-weight:600 !important;
    }
    .na-module table.dataTable tbody td {
      font-family:'Consolas','IBM Plex Mono',monospace !important;
      font-size:11px !important; color:#1e293b !important;
    }
    .na-module .nav-tabs > li > a { font-size:12px; font-weight:500; color:#475569; padding:5px 13px; }
    .na-module .nav-tabs > li.active > a { color:#0f172a; font-weight:600; }
  "))

  # ── Shared download row ──────────────────────────────────────────────────
  dlrow <- function(...) tags$div(class="na-dl-row", ...)

  tags$div(class="na-module", custom_css,

    # ── Header ─────────────────────────────────────────────────────────────
    tags$div(class="na-header",
      tags$div(class="na-header-title",
        " Null Allele Estimation \u00b7 FST-ENA \u00b7 DCSE-INA"),
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

    # ── Value boxes ─────────────────────────────────────────────────────────
    tags$div(class="na-vbox-row",
      tags$div(class="na-vbox",
        tags$div(class="na-vbox-icon",style="background:#e0f2fe;color:#0369a1;"),
        tags$div(tags$div(class="na-vbox-label","Loci"),
                 tags$div(class="na-vbox-val",uiOutput(ns("vb_loci"))))),
      tags$div(class="na-vbox",
        tags$div(class="na-vbox-icon",style="background:#dcfce7;color:#166534;"),
        tags$div(tags$div(class="na-vbox-label","Populations"),
                 tags$div(class="na-vbox-val",uiOutput(ns("vb_pops"))))),
      tags$div(class="na-vbox",
        tags$div(class="na-vbox-icon",style="background:#f3e8ff;color:#7e22ce;"),
        tags$div(tags$div(class="na-vbox-label","Individuals"),
                 tags$div(class="na-vbox-val",uiOutput(ns("vb_n"))))),
      tags$div(class="na-vbox",
        tags$div(class="na-vbox-icon",style="background:#fef9c3;color:#854d0e;"),
        tags$div(tags$div(class="na-vbox-label","Avg p_nulls"),
                 tags$div(class="na-vbox-val",uiOutput(ns("vb_avg_null"))))),
      tags$div(class="na-vbox",
        tags$div(class="na-vbox-icon",style="background:#fce7f3;color:#9d174d;"),
        tags$div(tags$div(class="na-vbox-label","Max p_nulls"),
                 tags$div(class="na-vbox-val",uiOutput(ns("vb_max_null"))))),
      tags$div(class="na-vbox",
        tags$div(class="na-vbox-icon",style="background:#ccfbf1;color:#0d9488;"),
        tags$div(tags$div(class="na-vbox-label","Global FST-ENA"),
                 tags$div(class="na-vbox-val",uiOutput(ns("vb_fst_ena")))))
    ),

    # ════════════════════════════════════════════════════════════════════════
    # SETUP PANEL
    # ════════════════════════════════════════════════════════════════════════
    tags$div(class="na-panel",
      tags$div(class="na-panel-head",
        tags$div(class="na-panel-title","Setup")),
      tags$div(class="na-panel-body",

        # ── (1) Missing genotype coding per locus ────────────────────────────
        tags$div(class="na-warn",
          tags$p(style="margin:.25rem 0;",
            "Please choose how to code missing data for each locus:", tags$br(),
            tags$strong("0"), " = true missing data (ignored by the algorithm);", tags$br(),
            tags$strong("999999"), " = homozygote for allele 999 (code for all null alleles)"),
          tags$p(style="margin:.5rem 0 0;font-weight:600;color:#92400e;",
            "Please make sure you do not already have any allele coded as 999.")
        ),
        uiOutput(ns("locus_coding_ui")),

        tags$hr(style="margin:1rem 0;"),

        # ── (2) Bootstrap parameters ────────────────────────────────────────
        tags$strong("Bootstrap parameters", style="font-size:12px; color:#1e293b;"),
        tags$br(), tags$br(),
        fluidRow(
          column(4,
            numericInput(ns("nboot"),
              label = "Number of replicates (bootstrap over loci):",
              value = 5000, min = 100, max = 99999, step = 1000)),
          column(4,
            numericInput(ns("nboot_subs"),
              label = "Number of replicates (bootstrap over sub-samples):",
              value = 5000, min = 100, max = 99999, step = 1000)),
          column(4,
            numericInput(ns("alpha"),
              label = "Confidence interval — alpha:",
              value = 0.05, min = 0.0001, max = 0.5, step = 0.01))
        ),

        tags$hr(style="margin:1rem 0;"),

        # ── (3) Output files ─────────────────────────────────────────────────
        tags$strong("Output files", style="font-size:12px; color:#1e293b;"),
        tags$br(), tags$br(),
        fluidRow(
          column(6,
            tags$div(style="display:flex; align-items:flex-end; gap:8px;",
              tags$div(style="flex:1;",
                textInput(ns("out_dir_display"), "Please choose a folder for output files:",
                          value = "", placeholder = "(no folder chosen \u2014 files download to your browser instead)")),
              shinyFiles::shinyDirButton(ns("out_dir_browse"), "Browse", "Choose output folder",
                                          class = "btn-action-secondary", style="margin-bottom:15px;"))),
          column(6,
            textInput(ns("out_root"), "Please choose a root for the name of output files:",
                      value = "", placeholder = "e.g. BoophilusAdultsDataCattle"))
        ),
        tags$p(style="color:#777;font-size:11px;",
          "File names = root + description (e.g. ", tags$code("<root>null_allele_frequencies.txt"),
          "). No date is added (already shown by the file explorer) \u2014 if you re-run with a ",
          "different missing-data coding and want to keep both results, add your own suffix below."),
        textInput(ns("out_suffix"), "Optional suffix to distinguish this run (e.g. \"1\"):", value = ""),
        tags$p(style="color:#777;font-size:11px;",
          "Files are saved as tab-delimited ", tags$strong(".txt"), " (not .csv)."),

        tags$hr(style="margin:1rem 0;"),

        # ── (3) Run computation ─────────────────────────────────────────────
        tags$strong("(3) Run all computations + generate output files",
                    style="font-size:12px; color:#1e293b;"),
        tags$br(), tags$br(),
        fluidRow(
          column(4,
            actionButton(ns("run_all"),
              label = tagList(tags$strong("  Compute + Bootstrap + Export")),
              class = "na-btn-run btn",
              width = "100%"))
        ),
        br(),
        uiOutput(ns("ui_run_status"))
      )
    ),

    # ════════════════════════════════════════════════════════════════════════
    # OUTPUT FILES PANEL
    # ════════════════════════════════════════════════════════════════════════
    tags$div(class="na-panel",
      tags$div(class="na-panel-head",
        tags$div(class="na-panel-title", " Output files")),
      tags$div(class="na-panel-body",
        fluidRow(
          # File 1
          column(2,
            tags$div(class="na-panel", style="border-color:#bfdbfe;",
              tags$div(class="na-panel-head", style="background:#eff6ff;",
                tags$div(class="na-panel-title", style="color:#1d4ed8;font-size:11px;",
                  uiOutput(ns("ui_filename_1"), inline=TRUE))),
              tags$div(class="na-panel-body", style="font-size:11px;color:#334155;",
                "p_nulls per locus \u00d7 subsample",
                tags$br(), "Global weighted mean per locus",
                tags$br(), br(),
                uiOutput(ns("ui_dl_file1"))
              )
            )
          ),
          # File 2
          column(2,
            tags$div(class="na-panel", style="border-color:#99f6e4;",
              tags$div(class="na-panel-head", style="background:#f0fdfa;",
                tags$div(class="na-panel-title", style="color:#0d9488;font-size:11px;",
                  uiOutput(ns("ui_filename_2"), inline=TRUE))),
              tags$div(class="na-panel-body", style="font-size:11px;color:#334155;",
                "Per locus + multilocus FST / FST-ENA",
                tags$br(), "CI over loci and over sub-samples",
                uiOutput(ns("ui_dl_file2"))
              )
            )
          ),
          # File 3
          column(2,
            tags$div(class="na-panel", style="border-color:#e9d5ff;",
              tags$div(class="na-panel-head", style="background:#faf5ff;",
                tags$div(class="na-panel-title", style="color:#7c3aed;font-size:11px;",
                  uiOutput(ns("ui_filename_3"), inline=TRUE))),
              tags$div(class="na-panel-body", style="font-size:11px;color:#334155;",
                "FST, FST-ENA, DCSE, DCSE-INA",
                tags$br(), "Per pair of sub-samples, all loci combined",
                uiOutput(ns("ui_dl_file3"))
              )
            )
          ),
          # File 4
          column(2,
            tags$div(class="na-panel", style="border-color:#fcd34d;",
              tags$div(class="na-panel-head", style="background:#fffbeb;",
                tags$div(class="na-panel-title", style="color:#92400e;font-size:11px;",
                  uiOutput(ns("ui_filename_4"), inline=TRUE))),
              tags$div(class="na-panel-body", style="font-size:11px;color:#334155;",
                "FST, FST-ENA, DCSE, DCSE-INA",
                tags$br(), "Half-matrix per locus, per pair",
                uiOutput(ns("ui_dl_file4"))
              )
            )
          ),
          # File 5
          column(4,
            tags$div(class="na-panel", style="border-color:#fca5a5;",
              tags$div(class="na-panel-head", style="background:#fef2f2;",
                tags$div(class="na-panel-title", style="color:#991b1b;font-size:11px;",
                  uiOutput(ns("ui_filename_5"), inline=TRUE))),
              tags$div(class="na-panel-body", style="font-size:11px;color:#334155;",
                "All bootstrap replicate values (over loci and over sub-samples)",
                tags$br(),
                uiOutput(ns("ui_dl_file5"))
                # plotly::plotlyOutput(ns("boot_dist_plot"), height="220px")
              )
            )
          )
        )
      )
    ),

    # ════════════════════════════════════════════════════════════════════════
    # RESULTS TABS — for visual inspection
    # ════════════════════════════════════════════════════════════════════════
    tabsetPanel(id = ns("na_tabs"), type = "tabs",

      # ── TAB 1: Null allele frequencies ────────────────────────────────── #
      tabPanel(title = tagList(" Null allele frequencies"),
               value = "tab_na", br(),
        tags$div(class="na-info",
          "Reproduces FreeNA's own null-allele-frequency report: the EM algorithm ",
          "(Dempster, Laird & Rubin 1977) estimated per locus \u00d7 population below, ",
          "and the N-weighted per-locus summary (Av(p_nulls), Av(N_exp_blanks), ",
          "f(expBlanks), one-sided binomial test p-value, and chosen blank coding) further down."
        ),
        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              " p_nulls per locus \u00d7 population (EM algorithm)")),
          tags$div(class="na-panel-body",
            DT::DTOutput(ns("dt_t1")))),
        tags$br(),
        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              " Per-locus summary (N-weighted mean, FreeNA report format)")),
          tags$div(class="na-panel-body",
            DT::DTOutput(ns("dt_t2"))))
      ),

      # ── TAB 2: FST & FST-ENA ──────────────────────────────────────────── #
      tabPanel(title = tagList(" FST / FST-ENA"),
               value = "tab_fst", br(),

        tags$div(class="na-info",
          tags$strong("Global multilocus FST"), " \u2014 Weir (1996) / Genepop method. ",
          tags$strong("FST-ENA"), ": EM-corrected frequencies, Excluding Null Alleles \u2014 Chapuis & Estoup (2007).",
          tags$br(),
          "Bootstrap CI over loci (resample loci with replacement) and over sub-samples ",
          "(resample individuals within each population with replacement)."
        ),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              " Per-locus FST and FST-ENA")),
          tags$div(class="na-panel-body",
            DT::DTOutput(ns("dt_fst_global")))),
        br(),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              " Bootstrap CI \u2014 Global FST and FST-ENA")),
          tags$div(class="na-panel-body",
            uiOutput(ns("ui_boot_global_fst")))),
        br(),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              " Pairwise FST and FST-ENA \u2014 lower triangle matrix")),
          tags$div(class="na-panel-body",
            fluidRow(
              column(5,
                radioButtons(ns("fst_pair_display"), "Display:",
                  choices = c(
                    "Raw FST (uncorrected)" = "raw",
                    "FST-ENA (corrected)"   = "ena",
                    "Both side by side"     = "both"),
                  selected = "both", inline = TRUE))),
            uiOutput(ns("ui_fst_pair_matrix")))),
        br(),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              " Bootstrap CI \u2014 Pairwise FST-ENA (over loci)")),
          tags$div(class="na-panel-body",
            uiOutput(ns("ui_boot_pair_fst"))))
      ),

      # ── TAB 3: DCSE / DCSE-INA ────────────────────────────────────────── #
      tabPanel(title = tagList(" DCSE / DCSE-INA"),
               value = "tab_dc", br(),

        tags$div(class="na-info",
          tags$strong("Cavalli-Sforza & Edwards (1967) chord distance."),
          " DCSE-INA includes the null allele as an extra state \u2014 Chapuis & Estoup (2007).",
          tags$br(),
          "DCSE(i,j) = (2/\u03c0)\u00d7\u221a[2\u00d7(1\u2212\u03a3\u221a(p_ik\u00d7p_jk))]  ",
          "INA: corrdgenefreq + null allele appended (freq = rd[locus, pop])."
        ),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              " Pairwise DCSE and DCSE-INA \u2014 lower triangle matrix")),
          tags$div(class="na-panel-body",
            fluidRow(
              column(5,
                radioButtons(ns("dc_display"), "Display:",
                  choices = c(
                    "Raw DCSE (uncorrected)" = "raw",
                    "DCSE-INA (corrected)"   = "ina",
                    "Both side by side"      = "both"),
                  selected = "both", inline = TRUE))),
            uiOutput(ns("ui_dc_matrix")))),
        br(),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              " Bootstrap CI \u2014 Pairwise DCSE-INA (over loci)")),
          tags$div(class="na-panel-body",
            uiOutput(ns("ui_boot_pair_dc"))))
      ),

      # ── TAB 4: Per-locus x pair ───────────────────────────────────────── #
      tabPanel(title = tagList(" Per-locus \u00d7 pair"),
               value = "tab_locus_pair", br(),

        tags$div(class="na-info",
          "FST, FST-ENA, DCSE and DCSE-INA for each locus \u00d7 pair of populations.",
          " Useful for detecting outlier loci."
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
              " FST and FST-ENA per locus \u00d7 pair")),
          tags$div(class="na-panel-body",
            DT::DTOutput(ns("dt_fst_locus")))),
        br(),

        tags$div(class="na-panel",
          tags$div(class="na-panel-head",
            tags$div(class="na-panel-title",
              " DCSE and DCSE-INA per locus \u00d7 pair")),
          tags$div(class="na-panel-body",
            DT::DTOutput(ns("dt_dc_locus"))))
      )

    ) # end tabsetPanel
  )   # end tags$div.na-module
}
