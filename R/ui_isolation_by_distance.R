# ui_isolation_by_distance.R
# Tab: Isolation by Distance
# Rousset (1997) linearised FST/(1-FST) vs geographic distance, Mantel test.
# Three regression lines: average, upper 95% CI, lower 95% CI.

isolation_by_distance_UI <- function(id) {
  ns <- NS(id)

  fluidPage(
    tags$head(gs_head()),

    module_banner(
      "map-marker-alt",
      "Isolation by Distance",
      "Rousset (1997) linearisation \u00b7 FST\u2044(1\u2212FST) vs geographic distance \u00b7 Mantel test",
      "#2CBF9F"
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
          "Results are functional but the module is still being validated. Use with caution.")
      )
    ),

    tags$div(
      class = "spg-method-note", style = "border-left-color:#2CBF9F;",
      HTML(paste0(
        "Tests whether genetic differentiation increases with geographic distance (IBD). ",
        "Pairwise F<sub>ST</sub> (WC84) is linearised as F<sub>ST</sub>\u2044(1\u2212F<sub>ST</sub>) ",
        "and regressed against geographic distance (1D) or ln(distance) (2D \u2014 Rousset 1997). ",
        "Three regression lines are fitted: through the mean (Average), ",
        "the upper 95% CI bound (ls) and the lower 95% CI bound (li). ",
        "The slope <b>b</b> estimates 1/N<sub>b</sub> (neighbourhood size). ",
        "Significance is assessed with a Mantel permutation test. ",
        "<b>Requires latitude/longitude columns set during data import.</b><br><br>",
        "<b>References:</b> Rousset F. 1997. <em>Genetics</em> 145:1219\u20131228. ",
        "| de Mee\u00fbs T <em>et al.</em> 2006. <em>Infect Genet Evol</em>."
      ))
    ),

    # ── Configuration ─────────────────────────────────────────────────────────
    fluidRow(
      box(
        width = 12,
        title = div(
          style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
          icon("sliders-h"), " IBD parameters"
        ),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            radioButtons(
              ns("model"),
              "Habitat model",
              choices = c(
                "2D \u2014 ln(distance km)" = "2D",
                "1D \u2014 distance km"     = "1D"
              ),
              selected = "2D"
            ),
            numericInput(ns("n_boot_pw"),  "Bootstrap per pair (CI):",
                         value = 500, min = 100, max = 5000, step = 100),
            numericInput(ns("n_perm"),     "Mantel permutations:",
                         value = 9999, min = 99, max = 99999, step = 1000),
            tags$hr(),
            actionButton(
              ns("run"), "Run IBD Analysis",
              icon  = icon("rocket"),
              class = "btn-action-primary btn-block",
              style = "font-weight:bold;"
            )
          ),
          column(9,
            h4(icon("chart-line"), "Results summary",
               style = "font-weight:600; color:#2c3e50; margin-bottom:15px;"),
            fluidRow(
              column(3, valueBoxOutput(ns("box_npops"),    width = NULL)),
              column(3, valueBoxOutput(ns("box_npairs"),   width = NULL)),
              column(3, valueBoxOutput(ns("box_mantel_r"), width = NULL)),
              column(3, valueBoxOutput(ns("box_pval"),     width = NULL))
            ),
            # Regression summary table (b, Nb, Nem for 3 lines)
            tags$h5("Regression parameters",
                    style = "font-weight:600; margin-top:14px; color:#2c3e50;"),
            DT::DTOutput(ns("reg_table"))
          )
        )
      )
    ),

    # ── IBD plot ───────────────────────────────────────────────────────────────
    fluidRow(
      box(
        width = 8,
        title = div(
          style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
          icon("chart-line"), " IBD plot"
        ),
        solidHeader = FALSE,
        plotly::plotlyOutput(ns("ibd_plot"), height = "460px")
      ),
      box(
        width = 4,
        title = div(
          style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
          icon("info-circle"), " Interpretation"
        ),
        solidHeader = FALSE,
        tags$div(
          style = "font-size:13px; line-height:1.8;",
          tags$p(tags$strong("Three regression lines")),
          tags$ul(
            style = "font-size:12px; padding-left:16px; line-height:1.9;",
            tags$li(tags$span(style="color:#333a43; font-weight:600;", "Average"), " \u2014 regression through point estimates of F",tags$sub("ST"),"\u2044(1\u2212F",tags$sub("ST"),")"),
            tags$li(tags$span(style="color:#B40F20; font-weight:600;", "Upper CI (ls)"), " \u2014 regression through upper 95% CI bounds"),
            tags$li(tags$span(style="color:#3B9AB2; font-weight:600;", "Lower CI (li)"), " \u2014 regression through lower 95% CI bounds")
          ),
          tags$p(tags$strong("Slope b"), " \u2014 in the 2D model: b = 1/N",tags$sub("b")," where N",tags$sub("b")," is the neighbourhood size."),
          tags$p(tags$strong("N",tags$sub("b")," = 1/b"), " \u2014 number of individuals in the dispersal neighbourhood."),
          tags$p(tags$strong("N",tags$sub("em")," = 1/(2\u03c0b)"), " \u2014 effective number of migrants per generation."),
          tags$hr(),
          tags$p(style = "color:#777; font-size:12px;",
            "Rousset (1997) Genetics 145:1219. de Mee\u00fbs (2006) Infect Genet Evol.")
        )
      )
    ),

    # ── Pairwise FST table ─────────────────────────────────────────────────────
    fluidRow(
      box(
        width = 7,
        title = div(
          style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
          icon("table"), " Pairwise F\u209b\u209c & linearised values"
        ),
        solidHeader = FALSE,
        DT::DTOutput(ns("fst_table")),
        tags$br(),
        downloadButton(ns("dl_fst"), "Download table",
                       class = "btn-action-secondary btn-sm")
      ),
      box(
        width = 5,
        title = div(
          style = "background:#FFFFFF; padding:10px; color:#333a43; font-weight:600;",
          icon("ruler"), " Pairwise distances (km)"
        ),
        solidHeader = FALSE,
        DT::DTOutput(ns("dist_table")),
        tags$br(),
        downloadButton(ns("dl_dist"), "Download distances",
                       class = "btn-action-secondary btn-sm")
      )
    )
  )
}
