# module/ui_null_alleles.R
# Null allele frequency estimation (EM), FST-ENA, DCSE-INA
# Simplified UI per supervisor feedback:
#   - Radio buttons for missing genotype coding, default = 000000
#   - Single bootstrap panel: n replicates + CI level choice
#   - 4 automatic output files
#
# UI rebuilt to match the shared ShinyPopGen design system used by every
# other module (module_banner, shinydashboard box(), valueBox, spg-method-note,
# btn-action-primary/secondary, btn-download-primary/secondary, tabsetPanel).
# Server output IDs are unchanged so server_null_alleles.R keeps working as-is.
#
# References:
#   Dempster, Laird & Rubin (1977)  — EM algorithm
#   Chapuis & Estoup (2007)         — FreeNA: ENA and INA corrections
#   Weir & Cockerham (1984)         — FST unbiased moment estimator
#   Cavalli-Sforza & Edwards (1967) — Chord genetic distance (DCSE)

null_alleles_UI <- function(id) {
  ns <- NS(id)

  # ── Supplemental CSS — only the bits with no shared-system equivalent:
  #    per-locus coding grid, pairwise matrices, bootstrap-result callouts.
  #    Recolored to the app's own palette (no more standalone dark/neon theme).
  supplemental_css <- tags$style(HTML("
    .na-info {
      background:#f5f7fa; border-left:4px solid #8D8680; border-radius:3px;
      padding:10px 14px; font-size:13px; line-height:1.65; color:#2c3e50; margin-bottom:14px;
    }
    .na-warn {
      background:#fff8e6; border-left:4px solid #E1AF00; border-radius:3px;
      padding:10px 14px; font-size:13px; line-height:1.65; color:#5c4400; margin-bottom:14px;
    }
    .na-locus-grid { display:flex; flex-wrap:wrap; gap:8px; margin-top:8px; }
    .na-locus-item {
      background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px;
      padding:.5rem .8rem; min-width:150px; flex:1;
    }
    .na-locus-item .control-label { display:none; } /* hide redundant label */
    .na-locus-name {
      font-size:20px; font-weight:700; color:#333a43;
      font-family:'Consolas','IBM Plex Mono',monospace; margin-bottom:4px;
    }
    .na-locus-item .radio { margin:2px 0; }
    .na-locus-item .radio label { font-size:13px; color:#333a43; }

    .na-boot-result {
      background:#f5f7fa; border:1px solid #8ea1b9; border-radius:6px;
      padding:.65rem 1rem; font-size:12.5px; color:#333a43;
      font-family:'Consolas','IBM Plex Mono',monospace; line-height:1.9; margin-top:.5rem;
    }
    .na-boot-result strong { color:#0c4a6e; }

    .na-matrix-wrap { overflow-x:auto; margin-top:.5rem; }
    .na-matrix { border-collapse:collapse; font-size:11px; font-family:'Consolas','IBM Plex Mono',monospace; width:100%; }
    .na-matrix th { background:#f8fafc; color:#475569; font-weight:600; padding:4px 9px; border:1px solid #e2e8f0; font-size:10.5px; white-space:nowrap; }
    .na-matrix td { padding:4px 9px; border:1px solid #e2e8f0; color:#1e293b; text-align:right; white-space:nowrap; font-size:11px; }
    .na-matrix tr:nth-child(even) td { background:#f8fafc; }
    .na-matrix .diag  { background:#f1f5f9 !important; color:#94a3b8; text-align:center; }
    .na-matrix .upper { color:#cbd5e1; text-align:center; }
    .na-matrix .lbl   { font-weight:700; color:#333a43; text-align:left; white-space:nowrap; }

    .na-dl-row { display:flex; gap:6px; flex-wrap:wrap; margin-top:.5rem; }
    .na-dl-row .btn { font-size:11px; padding:3px 12px; }

    .na-filecard { margin-bottom: 14px; }
    .na-filecard .fname { margin-top:8px; font-size:11px; word-break:break-all; color:#666; }
  "))

  box_title_style <- "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;"

  fluidPage(
    tags$head(gs_head()),
    supplemental_css,

    module_banner("circle-notch", "Null Allele Estimation \u00b7 FST-ENA \u00b7 DCSE-INA",
      "EM algorithm \u00b7 Dempster, Laird & Rubin (1977) \u00b7 FreeNA \u2014 Chapuis & Estoup (2007) \u00b7 Weir & Cockerham (1984) \u00b7 Cavalli-Sforza & Edwards (1967)",
      "#8D8680"),

    tags$div(class = "spg-method-note", style = "border-left-color:#8D8680;",
      HTML(paste0(
        "<b>EM algorithm</b> (Dempster, Laird &amp; Rubin 1977) estimates the null allele frequency ",
        "at each locus \u00d7 population, following the <b>FreeNA</b> approach (Chapuis &amp; Estoup 2007). ",
        "<br><br>",
        "<b>F<sub>ST</sub>-ENA</b>: multilocus and pairwise F<sub>ST</sub> (Weir &amp; Cockerham 1984), ",
        "Excluding Null Alleles \u2014 corrected using the EM-estimated frequencies. ",
        "<b>D<sub>CSE</sub>-INA</b>: Cavalli-Sforza &amp; Edwards (1967) chord distance, Including the ",
        "null allele as an extra allelic state.",
        "<br><br>",
        "<b>Bootstrap confidence intervals</b> are computed by two resampling schemes:",
        "<ul style='margin:4px 0 0 16px;'>",
        "<li><b>Loci</b>, resampled with replacement across the whole locus set (multilocus estimates only).</li>",
        "<li><b>Sub-samples</b> (populations), resampled as whole blocks with replacement (multilocus and per-locus).</li>",
        "</ul>"
      ))
    ),

    # ════════════════════════════════════════════════════════════════════
    # SUMMARY — value boxes
    # ════════════════════════════════════════════════════════════════════
    # fluidRow(
    #   box(
    #     width = 12, solidHeader = TRUE, status = "primary",
    #     title = div(style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
    #                 icon("chart-bar"), " Summary "),
    #     fluidRow(
    #       column(2, valueBoxOutput(ns("vb_loci"),  width = NULL)),
    #       column(2, valueBoxOutput(ns("vb_pops"),  width = NULL)),
    #       column(2, valueBoxOutput(ns("vb_n"),     width = NULL)),
    #       column(2, valueBoxOutput(ns("vb_avg_null"), width = NULL)),
    #       column(2, valueBoxOutput(ns("vb_max_null"), width = NULL)),
    #       column(2, valueBoxOutput(ns("vb_fst_ena"),  width = NULL))
    #     )
    #   )
    # ),

    # ════════════════════════════════════════════════════════════════════
    # SETUP
    # ════════════════════════════════════════════════════════════════════
    fluidRow(
      box(
        width = 12,
        title = div(style = box_title_style, icon("sliders"), "Setup"),
        solidHeader = TRUE, status = "primary",

        h4(icon("code-branch"), "(1) Missing genotype coding per locus"),
        tags$div(class = "na-warn",
          tags$p(style = "margin:.25rem 0;",
            "Please choose how to code missing data for each locus:", tags$br(),
            tags$strong("0"), " = true missing data (ignored by the algorithm);", tags$br(),
            tags$strong("999999"), " = homozygote for allele 999 (code for all null alleles)"),
          tags$p(style = "margin:.5rem 0 0;font-weight:600;",
            "Please make sure you do not already have any allele coded as 999.")
        ),
        uiOutput(ns("locus_coding_ui")),

        tags$hr(),

        h4(icon("dice"), "(2) Bootstrap parameters"),
        fluidRow(
          column(3,
            numericInput(ns("nboot"),
              label = "Number of replicates (bootstrap over loci):",
              value = 5000, min = 100, max = 99999, step = 1000)),
          column(3,
            numericInput(ns("nboot_subs"),
              label = "Number of replicates (bootstrap over sub-samples):",
              value = 5000, min = 100, max = 99999, step = 1000)),
          column(3,
            numericInput(ns("alpha"),
              label = "Confidence interval \u2014 alpha:",
              value = 0.05, min = 0.0001, max = 0.5, step = 0.01)),
          column(3,
            numericInput(ns("boot_seed"),
              label = "Random seed (reproducibility):",
              value = 12345, min = 1, max = 2147483647, step = 1))
        ),
        tags$p(style = "color:#666;font-size:12px;",
          icon("info-circle"), " ",
          "Bootstrap resampling is random: point estimates (FST, FST-ENA, DCSE\u2026) never change, ",
          "but confidence interval bounds will shift slightly from run to run unless the seed is kept ",
          "the same. Re-run with the same seed, same data and same number of replicates to reproduce ",
          "the exact same confidence intervals \u2014 the seed used is recorded in every exported file."),

        tags$hr(),

        h4(icon("save"), "(3) Output files"),
        fluidRow(
          column(6,
            tags$div(style = "display:flex; align-items:flex-end; gap:8px;",
              tags$div(style = "flex:1;",
                textInput(ns("out_dir_display"), "Save files to this folder (optional):",
                          value = "", placeholder = "(not set \u2014 use the .txt buttons below instead)")),
              shinyFiles::shinyDirButton(ns("out_dir_browse"), "Browse", "Choose output folder",
                                          class = "btn-action-secondary", style = "margin-bottom:15px;"))),
          column(6,
            textInput(ns("out_root"), "Root for the name of output files:",
                      value = "", placeholder = "auto-filled from the imported data file name"))
        ),
        tags$p(style = "color:#666;font-size:12px;",
          icon("info-circle"), " ", tags$strong("Browse\u2026"), " opens a folder picker for ",
          tags$strong("this computer"), " (the one running this app). Pick a folder and every file will ",
          "be saved there automatically each time you click Compute \u2014 no need to click the .txt buttons ",
          "one by one. Leave it empty to just use the .txt download buttons below each result instead."),
        tags$p(style = "color:#666;font-size:12px;",
          "The root is proposed automatically from the name of the data file you imported ",
          "and you can freely edit or extend it \u2014 e.g. add your own notes such as which loci ",
          "were recoded to 999999. File names = root + description (e.g. ",
          tags$code("<root>null_allele_frequencies.txt"), "). No date is added (already shown by your ",
          "computer's file browser) \u2014 if you re-run with a different missing-data coding and want to ",
          "keep both results, use the suffix field to tell them apart."),
        textInput(ns("out_suffix"), "Optional suffix to distinguish this run (e.g. \"1\"):", value = ""),
        tags$p(style = "color:#666;font-size:12px;",
          "Files are saved as tab-delimited ", tags$strong(".txt"), " (not .csv)."),

        tags$hr(),

        h4(icon("rocket"), "(4) Run all computations + generate output files"),
        fluidRow(
          column(4,
            actionButton(ns("run_all"),
              label = tagList(icon("rocket"), tags$strong(" Compute + Bootstrap + Export")),
              class = "btn-action-primary btn-block",
              style = "font-weight: bold;"))
        ),
        br(),
        uiOutput(ns("ui_run_status"))
      )
    ),

    # ════════════════════════════════════════════════════════════════════
    # OUTPUT FILES — one card per exported file
    # ════════════════════════════════════════════════════════════════════
    h2("Output files", class = "section-title"),
    fluidRow(
      box(
        width = 12,
        title = div(style = box_title_style, icon("file-export"),
                    "5 files generated by \u201cCompute + Bootstrap + Export\u201d"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(2, tags$div(class = "spg-module-card na-filecard",
            tags$div(class = "card-icon", icon("table")),
            h5("p_nulls / locus"),
            p("Per locus \u00d7 subsample, plus global weighted mean per locus."),
            tags$div(class = "fname", uiOutput(ns("ui_filename_1"), inline = TRUE)),
            uiOutput(ns("ui_dl_file1"))
          )),
          column(2, tags$div(class = "spg-module-card na-filecard",
            tags$div(class = "card-icon", icon("chart-bar")),
            h5("FST / FST-ENA"),
            p("Per locus + multilocus, CI over loci and over sub-samples."),
            tags$div(class = "fname", uiOutput(ns("ui_filename_2"), inline = TRUE)),
            uiOutput(ns("ui_dl_file2"))
          )),
          column(2, tags$div(class = "spg-module-card na-filecard",
            tags$div(class = "card-icon", icon("route")),
            h5("Pairwise (long format)"),
            p("FST, FST-ENA, DCSE, DCSE-INA per pair of sub-samples, all loci combined."),
            tags$div(class = "fname", uiOutput(ns("ui_filename_3"), inline = TRUE)),
            uiOutput(ns("ui_dl_file3"))
          )),
          column(2, tags$div(class = "spg-module-card na-filecard",
            tags$div(class = "card-icon", icon("th")),
            h5("Per-locus half-matrices"),
            p("FST, FST-ENA, DCSE, DCSE-INA \u2014 per locus, per pair."),
            tags$div(class = "fname", uiOutput(ns("ui_filename_4"), inline = TRUE)),
            uiOutput(ns("ui_dl_file4"))
          )),
          column(4, tags$div(class = "spg-module-card na-filecard",
            tags$div(class = "card-icon", icon("dice")),
            h5("Bootstrap distributions"),
            p("All bootstrap replicate values (over loci and over sub-samples)."),
            tags$div(class = "fname", uiOutput(ns("ui_filename_5"), inline = TRUE)),
            uiOutput(ns("ui_dl_file5"))
          ))
        )
      )
    ),

    # ════════════════════════════════════════════════════════════════════
    # RESULTS TABS — for visual inspection
    # ════════════════════════════════════════════════════════════════════
    h2("Results", class = "section-title"),
    fluidRow(
      box(
        width = 12,
        title = div(style = box_title_style, icon("table"), "Detailed results"),
        solidHeader = TRUE, status = "primary",

        tabsetPanel(id = ns("na_tabs"), type = "tabs",

          # ── TAB 1: Null allele frequencies ────────────────────────────────── #
          tabPanel(title = tagList(icon("percent"), " Null allele frequencies"),
                   value = "tab_na", br(),
            tags$div(class = "na-info",
              "Reproduces FreeNA's own null-allele-frequency report: the EM algorithm ",
              "(Dempster, Laird & Rubin 1977) estimated per locus \u00d7 population below, ",
              "and the N-weighted per-locus summary (Av(p_nulls), Av(N_exp_blanks), ",
              "f(expBlanks), one-sided binomial test p-value, and chosen blank coding) further down."
            ),
            h4(icon("info-circle"), "p_nulls per locus \u00d7 population (EM algorithm)"),
            DT::DTOutput(ns("dt_t1")), br(),
            h4(icon("info-circle"), "Per-locus summary (N-weighted mean, FreeNA report format)"),
            DT::DTOutput(ns("dt_t2"))
          ),

          # ── TAB 2: FST & FST-ENA ──────────────────────────────────────────── #
          tabPanel(title = tagList(icon("chart-bar"), " FST / FST-ENA"),
                   value = "tab_fst", br(),
            tags$div(class = "na-info",
              tags$strong("Global multilocus FST"), " \u2014 Weir & Cockerham (1984) unbiased moment estimator. ",
              tags$strong("FST-ENA"), ": EM-corrected frequencies, Excluding Null Alleles \u2014 Chapuis & Estoup (2007).",
              tags$br(),
              "Bootstrap CI over loci (resample loci with replacement, multilocus estimates only) and over ",
              "sub-samples (resample populations as whole blocks with replacement, available both for the ",
              "multilocus estimate and per locus \u2014 see the per-locus table below)."
            ),
            h4(icon("table"), "Per-locus FST and FST-ENA"),
            DT::DTOutput(ns("dt_fst_global")), br(),

            h4(icon("chart-area"), "Bootstrap CI \u2014 Global FST and FST-ENA"),
            uiOutput(ns("ui_boot_global_fst")), br(),

            h4(icon("th"), "Pairwise FST and FST-ENA \u2014 lower triangle matrix"),
            fluidRow(
              column(5,
                radioButtons(ns("fst_pair_display"), "Display:",
                  choices = c(
                    "Raw FST (uncorrected)" = "raw",
                    "FST-ENA (corrected)"   = "ena",
                    "Both side by side"     = "both"),
                  selected = "both", inline = TRUE))),
            uiOutput(ns("ui_fst_pair_matrix")), br(),

            h4(icon("chart-area"), "Bootstrap CI \u2014 Pairwise FST-ENA (over loci)"),
            uiOutput(ns("ui_boot_pair_fst"))
          ),

          # ── TAB 3: DCSE / DCSE-INA ────────────────────────────────────────── #
          tabPanel(title = tagList(icon("route"), " DCSE / DCSE-INA"),
                   value = "tab_dc", br(),
            tags$div(class = "na-info",
              tags$strong("Cavalli-Sforza & Edwards (1967) chord distance."),
              " DCSE-INA includes the null allele as an extra state \u2014 Chapuis & Estoup (2007).",
              tags$br(),
              "DCSE(i,j) = (2/\u03c0)\u00d7\u221a[2\u00d7(1\u2212\u03a3\u221a(p_ik\u00d7p_jk))]  ",
              "INA: corrdgenefreq + null allele appended (freq = rd[locus, pop])."
            ),
            h4(icon("th"), "Pairwise DCSE and DCSE-INA \u2014 lower triangle matrix"),
            fluidRow(
              column(5,
                radioButtons(ns("dc_display"), "Display:",
                  choices = c(
                    "Raw DCSE (uncorrected)" = "raw",
                    "DCSE-INA (corrected)"   = "ina",
                    "Both side by side"      = "both"),
                  selected = "both", inline = TRUE))),
            uiOutput(ns("ui_dc_matrix")), br(),

            h4(icon("chart-area"), "Bootstrap CI \u2014 Pairwise DCSE-INA (over loci)"),
            uiOutput(ns("ui_boot_pair_dc"))
          ),

          # ── TAB 4: Per-locus x pair ───────────────────────────────────────── #
          tabPanel(title = tagList(icon("border-all"), " Per-locus \u00d7 pair"),
                   value = "tab_locus_pair", br(),
            tags$div(class = "na-info",
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
            h4(icon("table"), "FST and FST-ENA per locus \u00d7 pair"),
            DT::DTOutput(ns("dt_fst_locus")), br(),
            h4(icon("table"), "DCSE and DCSE-INA per locus \u00d7 pair"),
            DT::DTOutput(ns("dt_dc_locus"))
          )
        ),
        style = "padding: 10px;"
      )
    )
  )
}