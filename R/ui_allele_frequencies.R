ui_allele_frequencies <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    tags$head(gs_head()),
    
    module_banner("table", "Allele Frequencies",
      "Population-specific allele frequencies · N genotyped · N missing · Na / Ne / He / Ho / Fis",
      "#78B7C5"),
    
    tags$div(class = "spg-method-note", style = "border-left-color:#78B7C5;",
      HTML(paste0(
        "<b>Allele frequency analysis</b> following the Fstat format. ",
        "Displays allele frequencies for each marker across populations, ",
        "with sample sizes (N genotyped, N missing) and diversity indices ",
        "per locus-population combination.",
        "<br><br>",
        "<b>Diversity indices reported:</b>",
        "<ul style='margin:4px 0 0 16px;'>",
        "<li><b>Na</b>: Number of alleles</li>",
        "<li><b>Ne</b>: Effective number of alleles</li>",
        "<li><b>He</b>: Expected heterozygosity (gene diversity)</li>",
        "<li><b>Ho</b>: Observed heterozygosity</li>",
        "<li><b>Fis</b>: Inbreeding coefficient (per locus-population)</li>",
        "</ul>",
        "Allele frequencies include zeros for missing alleles."
      ))
    ),
    
    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("sliders-h"),
                    "Allele Frequency Analysis Parameters"),
        solidHeader = TRUE, status = "primary",
        fluidRow(
          column(3,
            h4(icon("filter"), "Selection"),
            selectInput(ns("fstat_population"), "Population:",
              choices = c("All populations" = "all"), multiple = FALSE),
            selectizeInput(ns("fstat_marker"), "Marker:",
              choices = NULL, multiple = FALSE,
              options = list(placeholder = "Select a marker")),
            br(),
            actionButton(ns("update_fstat"),
              label = tagList(icon("play"), " Compute Frequencies"),
              class = "btn-action-primary btn-block", style = "font-weight: bold;")
          ),
          column(9,
            h4(icon("chart-bar"), "Dataset Summary",
               style = "font-weight: 600; color: #2c3e50; margin-bottom: 15px;"),
            fluidRow(
              column(3,
                div(class = "af-vbox",
                  div(class = "af-vbox-icon", style = "background:#E6F1FB;color:#185FA5;", icon("users")),
                  div(
                    div(class = "af-vbox-label", "Individuals"),
                    uiOutput(ns("vb_individuals"))
                  )
                )
              ),
              column(3,
                div(class = "af-vbox",
                  div(class = "af-vbox-icon", style = "background:#EAF3DE;color:#3B6D11;", icon("map-marker-alt")),
                  div(
                    div(class = "af-vbox-label", "Populations"),
                    uiOutput(ns("vb_populations"))
                  )
                )
              ),
              column(3,
                div(class = "af-vbox",
                  div(class = "af-vbox-icon", style = "background:#EEEDFE;color:#534AB7;", icon("dna")),
                  div(
                    div(class = "af-vbox-label", "Markers"),
                    uiOutput(ns("vb_markers"))
                  )
                )
              ),
              column(3,
                div(class = "af-vbox",
                  div(class = "af-vbox-icon", style = "background:#FAEEDA;color:#854F0B;", icon("exclamation-triangle")),
                  div(
                    div(class = "af-vbox-label", "Missing data"),
                    uiOutput(ns("vb_missing"))
                  )
                )
              )
            ),
            br(),
            fluidRow(
              column(12,
                div(style = "display: flex; justify-content: flex-end; gap: 8px;",
                  downloadButton(ns("download_fstat_csv"), ".csv", class = "btn-download-primary"),
                  downloadButton(ns("download_fstat_txt"), ".txt", class = "btn-download-secondary")
                )
              )
            )
          )
        )
      )
    ),
    
    h2("Allele frequencies \u2014 populations in columns", class = "section-title"),
    tags$p(HTML(paste0(
      "<b>Allele frequencies by locus</b> with populations as columns. ",
      "Each cell shows: allele (frequency) with N genotyped · N missing. ",
      "The last columns report Na (number of alleles), Ne (effective number of alleles), ",
      "He (expected heterozygosity), Ho (observed heterozygosity), and Fis (inbreeding coefficient)."
    )), style = "font-size: 16px; line-height: 1.5; color: #2c3e50;"),
    
    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("th"), "Allele Frequency Table"),
        solidHeader = TRUE, status = "primary",
        DT::DTOutput(ns("fstat_table")),
        style = "padding: 10px;"
      )
    )
  )
}