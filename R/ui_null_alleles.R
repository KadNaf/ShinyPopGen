# module/ui_null_alleles.R
# Null allele frequency estimation (EM), FST-ENA correction, and DCSE-INA genetic distance
# Bootstrap 95% CI: over loci (global + per-locus) and over individuals
#
# References:
#   Dempster, Laird & Rubin (1977)  — EM algorithm
#   Chapuis & Estoup (2007)         — FreeNA: ENA and INA corrections for null alleles
#   Weir (1996)                     — FST following Genepop method
#   Cavalli-Sforza & Edwards (1967) — Chord genetic distance (DCSE)

null_alleles_UI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    tags$head(gs_head()),
    
    module_banner("atom", "Null Allele Estimation, FST-ENA Correction & DCSE-INA Genetic Distance",
      "EM algorithm · FreeNA — Chapuis & Estoup (2007) · Weir (1996) · Cavalli-Sforza & Edwards (1967)",
      "#78B7C5"),
    
    tags$div(class = "spg-method-note", style = "border-left-color:#78B7C5;",
      HTML(paste0(
        "<b>Null allele frequency estimation</b> using the EM algorithm (Dempster, Laird & Rubin, 1977). ",
        "FST correction following the ENA method (Excluding Null Alleles) and DCSE genetic distance ",
        "following the INA method (Including Null Alleles) — Chapuis & Estoup (2007).",
        "<br><br>",
        "<b>Missing genotype coding:</b>",
        "<ul style='margin:4px 0 0 16px;'>",
        "<li><b>999999</b> — missing coded as null homozygote → higher p_nulls (inferred from excess homozygosity)</li>",
        "<li><b>000000</b> — missing coded as absent / PCR failure → lower p_nulls (no null allele signal from missing data)</li>",
        "</ul>",
        "<b>Bootstrap confidence intervals (95% CI):</b> ",
        "locus bootstrap (resample loci with replacement) and individual bootstrap (resample individuals within populations)."
      ))
    ),
    
    # ── Value boxes ───────────────────────────────────────────────────────────
    fluidRow(
      column(12,
        fluidRow(
          column(2,
            div(class = "af-vbox",
              div(class = "af-vbox-icon", style = "background:#E6F1FB;color:#185FA5;", icon("dna")),
              div(
                div(class = "af-vbox-label", "Loci"),
                div(class = "af-vbox-val", uiOutput(ns("vb_loci")))
              )
            )
          ),
          column(2,
            div(class = "af-vbox",
              div(class = "af-vbox-icon", style = "background:#EAF3DE;color:#3B6D11;", icon("map-marker-alt")),
              div(
                div(class = "af-vbox-label", "Populations"),
                div(class = "af-vbox-val", uiOutput(ns("vb_pops")))
              )
            )
          ),
          column(2,
            div(class = "af-vbox",
              div(class = "af-vbox-icon", style = "background:#EEEDFE;color:#534AB7;", icon("users")),
              div(
                div(class = "af-vbox-label", "Individuals"),
                div(class = "af-vbox-val", uiOutput(ns("vb_n")))
              )
            )
          ),
          column(2,
            div(class = "af-vbox",
              div(class = "af-vbox-icon", style = "background:#FAEEDA;color:#854F0B;", icon("percentage")),
              div(
                div(class = "af-vbox-label", "Avg p_nulls"),
                div(class = "af-vbox-val", uiOutput(ns("vb_avg_null")))
              )
            )
          ),
          column(2,
            div(class = "af-vbox",
              div(class = "af-vbox-icon", style = "background:#FCE7F3;color:#9D174D;", icon("exclamation-triangle")),
              div(
                div(class = "af-vbox-label", "Max p_nulls"),
                div(class = "af-vbox-val", uiOutput(ns("vb_max_null")))
              )
            )
          ),
          column(2,
            div(class = "af-vbox",
              div(class = "af-vbox-icon", style = "background:#CCFBF1;color:#0D9488;", icon("chart-bar")),
              div(
                div(class = "af-vbox-label", "Global FST-ENA"),
                div(class = "af-vbox-val", uiOutput(ns("vb_fst_ena")))
              )
            )
          )
        )
      )
    ),
    
    # ── Per-locus treatment selector ──────────────────────────────────────────
    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("cog"), "Missing Genotype Coding — Per Locus"),
        solidHeader = TRUE, status = "primary",
        p("Select the missing genotype coding used in the original Genepop file for each locus:",
          style = "margin-bottom: 10px; color: #666;"),
        tags$div(class = "alert alert-warning", style = "padding: 10px; border-radius: 5px; background-color: #fffbeb; border-left: 4px solid #fcd34d;",
          icon("exclamation-triangle"),
          tags$strong("Coding guide:"),
          tags$br(),
          tags$strong("999999"), " — missing coded as null homozygote → higher p_nulls (inferred from excess homozygosity).",
          tags$br(),
          tags$strong("000000"), " — missing coded as absent / PCR failure → lower p_nulls (no null allele signal from missing data)."
        ),
        uiOutput(ns("locus_treatment_ui"))
      )
    ),
    
    # ── tabsetPanel ──────────────────────────────────────────────────────────
    h2("Null allele analysis results", class = "section-title"),
    
    fluidRow(
      box(
        width = 12,
        title = div(style = "background-color: #FFFFFF; padding: 10px; color: #333a43; font-weight: 600;",
                    icon("dna"), "Results"),
        solidHeader = TRUE, status = "primary",
        
        tabsetPanel(id = ns("na_tabs"), type = "tabs",
          
          # ══ TAB 1 ═════════════════════════════════════════════════════════════ #
          tabPanel(title = tagList(icon("table"), "Per locus × population"),
            br(),
            p("Null allele frequency estimated by EM per locus × population.",
              tags$strong("p_nulls"), ": estimated null allele frequency.",
              tags$strong("N"), ": total individuals in population.",
              tags$strong("N_exp_blanks = N × p_nulls²"), ": expected null homozygote count.",
              tags$strong("p_nulls×N"), ": expected null allele copies.",
              style = "font-size: 14px; line-height: 1.6; color: #2c3e50; margin-bottom: 15px;"),
            
            fluidRow(
              column(4, selectInput(ns("t1_locus"), "Locus:", choices = c("All loci" = "all"), selected = "all")),
              column(4, selectInput(ns("t1_pop"), "Population:", choices = c("All populations" = "all"), selected = "all")),
              column(4, actionButton(ns("run_t1"), label = tagList(icon("play"), " Compute"), 
                                     class = "btn-action-primary", style = "margin-top: 25px; width: 100%;"))
            ),
            br(),
            DT::DTOutput(ns("dt_t1")),
            br(),
            fluidRow(
              column(6, downloadButton(ns("dl_t1_csv"), ".csv", class = "btn-download-primary")),
              column(6, downloadButton(ns("dl_t1_txt"), ".txt", class = "btn-download-secondary"))
            )
          ),
          
          # ══ TAB 2 ═════════════════════════════════════════════════════════════ #
          tabPanel(title = tagList(icon("globe"), "Global summary per locus"),
            br(),
            p("Global summary across all populations per locus.",
              tags$strong("Av(N_exp_blanks)"), " = Σ(Nᵢ × pᵢ²): total expected null homozygotes.",
              tags$strong("Av(p_nulls)"), " = Σ(Nᵢ × pᵢ) / N_tot: N-weighted mean.",
              tags$strong("f(expBlanks) = Av(N_exp_blanks) / N_tot"), ".",
              style = "font-size: 14px; line-height: 1.6; color: #2c3e50; margin-bottom: 15px;"),
            
            fluidRow(
              column(4, selectInput(ns("t2_locus"), "Locus:", choices = c("All loci" = "all"), selected = "all")),
              column(4, actionButton(ns("run_t2"), label = tagList(icon("play"), " Compute"),
                                     class = "btn-action-primary", style = "margin-top: 25px; width: 100%;"))
            ),
            br(),
            DT::DTOutput(ns("dt_t2")),
            br(),
            fluidRow(
              column(6, downloadButton(ns("dl_t2_csv"), ".csv", class = "btn-download-primary")),
              column(6, downloadButton(ns("dl_t2_txt"), ".txt", class = "btn-download-secondary"))
            )
          ),
          
          # ══ TAB 3 — Global FST (ENA) + Bootstrap ══════════════════════════════ #
          tabPanel(title = tagList(icon("chart-bar"), "Global FST (ENA)"),
            br(),
            tags$div(class = "alert alert-info", style = "padding: 12px; border-radius: 5px; background-color: #eff6ff; border-left: 4px solid #0d9488;",
              icon("info-circle"),
              tags$strong("Global multilocus FST"), " — Weir (1996) following Genepop's method.",
              tags$br(),
              tags$strong("Raw FST"), ": computed from observed allele frequencies (null homozygotes excluded from denominator).",
              tags$br(),
              tags$strong("FST-ENA"), ": computed from EM-corrected allele frequencies (Excluding Null Alleles) — Chapuis & Estoup (2007)."
            ),
            
            tags$div(class = "spg-method-note", style = "border-left-color: #0d9488; margin-bottom: 15px;",
              HTML(paste0(
                "<b>Weir (1996):</b><br>",
                "FST = S1/S3   S1 = Σ[s²P×nc]   S3 = Σ[(s²P+s²I+s²G)×nc]<br>",
                "nc = (N_tot−Σni²/N_tot)/(r−1)   ENA: nA = corrdgenefreq×2ni ; AA_corr = AA×p/(p+2r)"
              ))
            ),
            
            fluidRow(
              column(3, actionButton(ns("run_fst_global"), label = tagList(icon("play"), " Compute"),
                                     class = "btn-action-primary", style = "width: 100%;"))
            ),
            br(),
            DT::DTOutput(ns("dt_fst_global")),
            br(),
            fluidRow(
              column(6, downloadButton(ns("dl_fst_global_csv"), ".csv", class = "btn-download-primary")),
              column(6, downloadButton(ns("dl_fst_global_txt"), ".txt", class = "btn-download-secondary"))
            ),
            
            hr(),
            h4(icon("random"), "Bootstrap confidence intervals (95% CI) — Locus bootstrap", style = "color: #4c1d95;"),
            p("Resample loci with replacement — 95% CI on the multilocus statistic + per-locus observed values.",
              "5000 replicates · percentile method · vectorised (fast, a few seconds).",
              style = "font-size: 13px; color: #4c1d95; margin-bottom: 15px;"),
            
            fluidRow(
              column(4, numericInput(ns("run_boot_fst_global_nboot"), label = "Replicates:", 
                                     value = 5000, min = 999, max = 9999, step = 1000)),
              column(4, selectInput(ns("run_boot_fst_global_type"), label = "Bootstrap type:",
                                    choices = c("Bootstrap over loci" = "loci"), selected = "loci")),
              column(4, actionButton(ns("run_boot_fst_global"), label = tagList(icon("random"), " Run Bootstrap"),
                                     class = "btn-download-primary", style = "margin-top: 25px; width: 100%; background: #7c3aed;"))
            ),
            uiOutput(ns("ui_boot_fst_global")),
            fluidRow(
              column(6, downloadButton(ns("dl_boot_fst_global_csv"), ".csv", class = "btn-download-primary")),
              column(6, downloadButton(ns("dl_boot_fst_global_txt"), ".txt", class = "btn-download-secondary"))
            )
          ),
          
          # ══ TAB 4 — Pairwise FST (ENA) + Bootstrap ════════════════════════════ #
          tabPanel(title = tagList(icon("exchange-alt"), "Pairwise FST (ENA)"),
            br(),
            p("Pairwise FST — Weir (1996) for each pair of populations.",
              "Lower triangle: raw FST (uncorrected) and FST-ENA (ENA-corrected).",
              tags$strong("NA"), ": insufficient sample size for the pair.",
              style = "font-size: 14px; line-height: 1.6; color: #2c3e50; margin-bottom: 15px;"),
            
            fluidRow(
              column(5, radioButtons(ns("fst_pair_type"), "Display:",
                choices = c("Raw FST (uncorrected)" = "raw", 
                           "FST-ENA (corrected)" = "ena", 
                           "Both side by side" = "both"),
                selected = "both", inline = FALSE)),
              column(3, actionButton(ns("run_fst_pair"), label = tagList(icon("play"), " Compute"),
                                     class = "btn-action-primary", style = "margin-top: 25px; width: 100%;"))
            ),
            br(),
            uiOutput(ns("ui_fst_pair_matrix")),
            br(),
            fluidRow(
              column(6, downloadButton(ns("dl_fst_pair_csv"), ".csv", class = "btn-download-primary")),
              column(6, downloadButton(ns("dl_fst_pair_txt"), ".txt", class = "btn-download-secondary"))
            ),
            br(),
            DT::DTOutput(ns("dt_fst_pair")),
            br(),
            fluidRow(
              column(6, downloadButton(ns("dl_fst_pair_long_csv"), ".csv", class = "btn-download-primary")),
              column(6, downloadButton(ns("dl_fst_pair_long_txt"), ".txt", class = "btn-download-secondary"))
            ),
            
            hr(),
            h4(icon("random"), "Bootstrap confidence intervals (95% CI)", style = "color: #4c1d95;"),
            p("Bootstrap over loci and/or individuals — 5000 replicates · percentile method.",
              style = "font-size: 13px; color: #4c1d95; margin-bottom: 15px;"),
            
            fluidRow(
              column(4, numericInput(ns("run_boot_fst_pair_nboot"), label = "Replicates:", 
                                     value = 5000, min = 999, max = 9999, step = 1000)),
              column(4, selectInput(ns("run_boot_fst_pair_type"), label = "Bootstrap type:",
                                    choices = c("Bootstrap over loci" = "loci",
                                               "Bootstrap over individuals" = "indiv",
                                               "Both (loci + individuals)" = "both_boot"),
                                    selected = "both_boot")),
              column(4, actionButton(ns("run_boot_fst_pair"), label = tagList(icon("random"), " Run Bootstrap"),
                                     class = "btn-download-primary", style = "margin-top: 25px; width: 100%; background: #7c3aed;"))
            ),
            uiOutput(ns("ui_boot_fst_pair")),
            fluidRow(
              column(6, downloadButton(ns("dl_boot_fst_pair_csv"), ".csv", class = "btn-download-primary")),
              column(6, downloadButton(ns("dl_boot_fst_pair_txt"), ".txt", class = "btn-download-secondary"))
            )
          ),
          
          # ══ TAB 5 — DCSE distance (INA) + Bootstrap ═══════════════════════════ #
          tabPanel(title = tagList(icon("ruler-combined"), "DCSE distance (INA)"),
            br(),
            p("Cavalli-Sforza & Edwards (1967) chord genetic distance — pairwise DCSE.",
              tags$strong("Raw DCSE"), ": computed from observed allele frequencies (null allele excluded).",
              tags$strong("DCSE-INA"), ": null allele included as an extra allelic state (Including Null Alleles) — Chapuis & Estoup (2007).",
              style = "font-size: 14px; line-height: 1.6; color: #2c3e50; margin-bottom: 15px;"),
            
            tags$div(class = "spg-method-note", style = "border-left-color: #0d9488; margin-bottom: 15px;",
              HTML(paste0(
                "<b>Cavalli-Sforza & Edwards (1967):</b><br>",
                "DCSE(i,j) = (2/π)×√[2×(1−Σ_k√(p_ik×p_jk))]   averaged over valid loci (CSprod≤1)<br>",
                "<b>INA:</b> corrdgenefreq + null allele appended (freq = rd[locus, pop])"
              ))
            ),
            
            fluidRow(
              column(5, radioButtons(ns("dc_type"), "Display:",
                choices = c("Raw DCSE (uncorrected)" = "raw", 
                           "DCSE-INA (corrected)" = "ina", 
                           "Both side by side" = "both"),
                selected = "both", inline = FALSE)),
              column(3, actionButton(ns("run_dc"), label = tagList(icon("play"), " Compute"),
                                     class = "btn-action-primary", style = "margin-top: 25px; width: 100%;"))
            ),
            br(),
            uiOutput(ns("ui_dc_matrix")),
            br(),
            fluidRow(
              column(6, downloadButton(ns("dl_dc_csv"), ".csv", class = "btn-download-primary")),
              column(6, downloadButton(ns("dl_dc_txt"), ".txt", class = "btn-download-secondary"))
            ),
            br(),
            DT::DTOutput(ns("dt_dc")),
            br(),
            fluidRow(
              column(6, downloadButton(ns("dl_dc_long_csv"), ".csv", class = "btn-download-primary")),
              column(6, downloadButton(ns("dl_dc_long_txt"), ".txt", class = "btn-download-secondary"))
            ),
            
            hr(),
            h4(icon("random"), "Bootstrap confidence intervals (95% CI)", style = "color: #4c1d95;"),
            p("Bootstrap over loci and/or individuals — 5000 replicates · percentile method.",
              style = "font-size: 13px; color: #4c1d95; margin-bottom: 15px;"),
            
            fluidRow(
              column(4, numericInput(ns("run_boot_dc_nboot"), label = "Replicates:", 
                                     value = 5000, min = 999, max = 9999, step = 1000)),
              column(4, selectInput(ns("run_boot_dc_type"), label = "Bootstrap type:",
                                    choices = c("Bootstrap over loci" = "loci",
                                               "Bootstrap over individuals" = "indiv",
                                               "Both (loci + individuals)" = "both_boot"),
                                    selected = "both_boot")),
              column(4, actionButton(ns("run_boot_dc"), label = tagList(icon("random"), " Run Bootstrap"),
                                     class = "btn-download-primary", style = "margin-top: 25px; width: 100%; background: #7c3aed;"))
            ),
            uiOutput(ns("ui_boot_dc")),
            fluidRow(
              column(6, downloadButton(ns("dl_boot_dc_csv"), ".csv", class = "btn-download-primary")),
              column(6, downloadButton(ns("dl_boot_dc_txt"), ".txt", class = "btn-download-secondary"))
            )
          ),
          
          # ══ TAB 6 — FST per locus x pair ══════════════════════════════════════ #
          tabPanel(title = tagList(icon("table"), "FST per locus × pair"),
            br(),
            p("Per-locus FST for each pair of populations.",
              "Useful for detecting outlier loci and comparing raw versus ENA-corrected estimates locus by locus.",
              style = "font-size: 14px; line-height: 1.6; color: #2c3e50; margin-bottom: 15px;"),
            
            fluidRow(
              column(3, selectInput(ns("fl_locus"), "Locus:", choices = c("All loci" = "all"), selected = "all")),
              column(3, selectInput(ns("fl_pop1"), "Population 1:", choices = c("All pairs" = "all"), selected = "all")),
              column(3, selectInput(ns("fl_pop2"), "Population 2:", choices = c("All pairs" = "all"), selected = "all")),
              column(3, actionButton(ns("run_fst_locus"), label = tagList(icon("play"), " Compute"),
                                     class = "btn-action-primary", style = "margin-top: 25px; width: 100%;"))
            ),
            br(),
            DT::DTOutput(ns("dt_fst_locus")),
            br(),
            fluidRow(
              column(6, downloadButton(ns("dl_fst_locus_csv"), ".csv", class = "btn-download-primary")),
              column(6, downloadButton(ns("dl_fst_locus_txt"), ".txt", class = "btn-download-secondary"))
            )
          )
        ),
        style = "padding: 10px;"
      )
    )
  )
}