# ui_null_alleles.R
# Null allele frequency estimation (EM) + Fst/Fst-ENA (Weir 1996/Genepop) +
# DCSE/DCSE-INA (Cavalli-Sforza & Edwards 1967), faithfully translated from
# the Pascal reference FreeNA_optm2R.pas (engine_freena.R).
#
# Bootstrap: OVER LOCI ONLY (as in the Pascal source) — one shared resampled
# loci sequence per replicate, applied simultaneously to every statistic.
# Minimum requirements (faithful to Pascal): nperm >= 100, nloc > 4.

isolation_by_distance_UI <- function(id) {
  ns <- NS(id)

  fluidPage(
    tags$head(gs_head()),

    module_banner(
      "atom",
      "Null Allele Estimation \u00b7 Fst-ENA \u00b7 DCSE-INA",
      "EM algorithm (Dempster, Laird & Rubin 1977) \u00b7 FreeNA (Chapuis & Estoup 2007) \u00b7 Weir (1996) \u00b7 Cavalli-Sforza & Edwards (1967)",
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
      )
    )
  )
}