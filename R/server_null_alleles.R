# module/server_null_alleles.R
#
# Null allele frequency estimation — EM algorithm (Dempster et al. 1977)
# STRICT IMPLEMENTATION of FreeNA (Chapuis & Estoup 2007)
# Faithful translation of Pascal source code (rDempster_per_locus)
#
# FreeNA coding convention:
#   - 0/0           → missing genotype (absentgeno, excluded)
#   - 999999/999999 → null homozygote (nnullhomo, included)
#   - Null heterozygotes → FORBIDDEN (fatal error)
#   - Other formats → normal alleles
#
# Inter-population weighted mean:
#   p̄_null = Σ(p_null_i × N_i × Ĥe_i) / Σ(N_i × Ĥe_i)
# ---------------------------------------------------------------------------

server_null_alleles <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns
    `%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

    # ════════════════════════════════════════════════════════════════════════
    # 1. DATA READING
    # ════════════════════════════════════════════════════════════════════════
    formatted_data <- reactive({
      tick <- rv$db_tick
      req(rv$con)
      con     <- rv$con
      tbl_raw <- rv$tbl_raw %||% "raw"
      message("[NULL] formatted_data refresh, db_tick = ", tick)

      ok_raw <- tryCatch(DBI::dbExistsTable(con, tbl_raw), error = function(e) FALSE)
      ok_par <- tryCatch(DBI::dbExistsTable(con, "params"), error = function(e) FALSE)
      if (!isTRUE(ok_raw) || !isTRUE(ok_par)) return(NULL)

      marker_json <- tryCatch(
        DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='marker_cols_raw'")$value[1],
        error = function(e) NA_character_
      )
      marker_cols_raw <- if (!is.na(marker_json) && nzchar(marker_json)) {
        tryCatch(jsonlite::fromJSON(marker_json), error = function(e) character(0))
      } else character(0)
      if (!length(marker_cols_raw)) stop("No marker_cols_raw found in DuckDB params.")

      genotype_format <- tryCatch(
        DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='genotype_format'")$value[1],
        error = function(e) NA_character_
      )
      if (is.na(genotype_format) || !nzchar(genotype_format))
        genotype_format <- if (any(grepl("(_1|\\.[0-9]+)$", marker_cols_raw))) "paired" else "string"

      pop_col <- tryCatch(
        DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='pop_col'")$value[1],
        error = function(e) NA_character_
      )
      if (is.na(pop_col) || !nzchar(pop_col)) pop_col <- NULL

      keep     <- unique(marker_cols_raw)
      keep_sql <- paste(vapply(keep, function(x) as.character(DBI::dbQuoteIdentifier(con, x)), character(1)), collapse = ", ")
      pop_sql  <- if (!is.null(pop_col)) paste0(", ", as.character(DBI::dbQuoteIdentifier(con, pop_col))) else ""
      raw_sql  <- sprintf("SELECT rowid AS individual%s, %s FROM %s;",
                           pop_sql, keep_sql, as.character(DBI::dbQuoteIdentifier(con, tbl_raw)))
      raw_df   <- as.data.frame(DBI::dbGetQuery(con, raw_sql), stringsAsFactors = FALSE)
      if (!nrow(raw_df)) stop("No genotype rows returned from DuckDB raw table.")

      pop_vector <- if (!is.null(pop_col) && pop_col %in% names(raw_df)) {
        as.character(raw_df[[pop_col]])
      } else rep("All", nrow(raw_df))

      pick_b <- function(locus, nms) {
        cands <- c(paste0(locus, "_1"), paste0(locus, "_2"), paste0(locus, ".", 1:9), paste0(locus, "_", 1:9))
        hit   <- cands[cands %in% nms]; if (length(hit)) hit[1] else NA_character_
      }

      if (identical(genotype_format, "paired")) {
        nms  <- names(raw_df)
        loci <- unique(sub("(_1|_2|\\.[0-9]+)$", "", marker_cols_raw))
        hap_df <- data.frame(row.names = seq_len(nrow(raw_df)))
        for (locus in loci) {
          b <- pick_b(locus, nms)
          if (!locus %in% nms || is.na(b) || !b %in% nms) next
          a_val <- as.character(raw_df[[locus]]); a_val[is.na(a_val) | trimws(a_val) == ""] <- "0"
          b_val <- as.character(raw_df[[b]]);     b_val[is.na(b_val) | trimws(b_val) == ""] <- "0"
          already <- grepl("/", a_val, fixed = TRUE) | grepl("-", a_val, fixed = TRUE)
          hap_df[[locus]] <- ifelse(already, a_val, paste0(a_val, "/", b_val))
        }
      } else {
        hap_df <- as.data.frame(raw_df[, marker_cols_raw, drop = FALSE], stringsAsFactors = FALSE)
        for (j in seq_along(hap_df)) {
          x <- as.character(hap_df[[j]]); x[is.na(x) | trimws(x) == ""] <- "0/0"; hap_df[[j]] <- x
        }
      }
      if (!ncol(hap_df)) stop("No locus columns could be reconstructed.")
      message("[NULL] haplotype dim: ", nrow(hap_df), " x ", ncol(hap_df))
      list(haplotype = hap_df, pop_vector = pop_vector)
    })

    # ════════════════════════════════════════════════════════════════════════
    # 2. REACTIVE STATE
    # ════════════════════════════════════════════════════════════════════════
    null_allele_results   <- reactiveVal(NULL)
    locus_coding_settings <- reactiveValues()

    observeEvent(rv$db_tick, ignoreInit = TRUE, {
      null_allele_results(NULL)
      for (nm in names(reactiveValuesToList(locus_coding_settings)))
        locus_coding_settings[[nm]] <- NULL
    })

    # ════════════════════════════════════════════════════════════════════════
    # 3. DYNAMIC UI — per-locus overrides
    # ════════════════════════════════════════════════════════════════════════
    output$locus_coding_ui <- renderUI({
      fd <- formatted_data()
      if (is.null(fd) || !ncol(fd$haplotype))
        return(div(style = "color:#6c757d; padding:8px;", "Please import a dataset first."))
      loci <- colnames(fd$haplotype)
      for (locus in loci)
        if (is.null(locus_coding_settings[[locus]])) locus_coding_settings[[locus]] <- "default"
      div(
        style = "max-height:280px; overflow-y:auto; border:1px solid #ddd; padding:10px; border-radius:5px;",
        lapply(loci, function(locus) {
          selectInput(
            inputId  = ns(paste0("coding_", locus)),
            label    = locus,
            choices  = c("Use default" = "default",
                        "0/0 = missing" = "0",
                        "999999/999999 = null homozygote" = "999999"),
            selected = locus_coding_settings[[locus]] %||% "default",
            width    = "100%"
          )
        })
      )
    })

    # ════════════════════════════════════════════════════════════════════════
    # 4. EM ALGORITHM — STRICT FreeNA IMPLEMENTATION
    # ════════════════════════════════════════════════════════════════════════
    
    em_freena <- function(genotypes, null_code = "999999", miss_code = "0",
                           tol = 1e-6, max_iter = 10000L) {
      
      # Clean and parse genotypes
      g <- as.character(genotypes)
      g <- g[!is.na(g) & nzchar(trimws(g))]
      
      parse2 <- function(gg) {
        # Handle / or - separators
        if (grepl("/", gg, fixed = TRUE)) {
          al <- strsplit(gg, "/", fixed = TRUE)[[1]]
        } else if (grepl("-", gg, fixed = TRUE)) {
          al <- strsplit(gg, "-", fixed = TRUE)[[1]]
        } else {
          al <- c(gg, gg)  # Homozygote genotype without separator
        }
        if (length(al) == 2) al else c(al[1], al[1])
      }
      
      # Genotype classification (strictly as in FreeNA)
      n_missing <- 0L      # absentgeno
      n_nullhom <- 0L      # nnullhomo (null homozygote individuals)
      allele_list <- character(0)
      genos_obs <- list()
      
      for (gg in g) {
        al <- parse2(gg)
        a1 <- trimws(al[1])
        a2 <- trimws(al[2])
        
        # Check missing codes
        both_miss <- (a1 %in% c("0", "", miss_code)) && (a2 %in% c("0", "", miss_code))
        both_null <- (a1 == null_code) && (a2 == null_code)
        one_null <- (a1 == null_code) || (a2 == null_code)
        
        if (both_miss) {
          n_missing <- n_missing + 1L
          next
        }
        if (both_null) {
          n_nullhom <- n_nullhom + 1L
          next
        }
        if (one_null) {
          # FreeNA: fatal error on null heterozygote
          warning("Null heterozygote ignored (forbidden in FreeNA): ", gg)
          next
        }
        
        # Valid genotype
        allele_list <- c(allele_list, a1, a2)
        genos_obs[[length(genos_obs) + 1]] <- c(a1, a2)
      }
      
      # N = efpop - absentgeno (individuals with any genotype)
      N <- n_nullhom + length(genos_obs)
      
      if (N == 0) {
        return(list(p_null = NA_real_, He = NA_real_, N = 0L,
                    N_missing = n_missing, N_nullhom = 0L, 
                    cpt = 0L, converged = FALSE))
      }
      
      # Only null homozygotes
      if (length(genos_obs) == 0 && n_nullhom > 0) {
        return(list(p_null = 1.0, He = 0.0, N = N,
                    N_missing = n_missing, N_nullhom = n_nullhom, 
                    cpt = 0L, converged = TRUE))
      }
      
      # Unique alleles
      alleles_uniq <- sort(unique(allele_list))
      A <- length(alleles_uniq)
      
      # Observed allele frequencies (genefreq in FreeNA)
      N_obs <- length(genos_obs)  # = N - n_nullhom
      al_table <- table(factor(allele_list, levels = alleles_uniq))
      genefreq <- as.numeric(al_table) / (2 * N_obs)
      names(genefreq) <- alleles_uniq
      
      # Count homozygotes and heterozygotes per allele
      n_hom <- setNames(integer(A), alleles_uniq)
      n_het <- setNames(integer(A), alleles_uniq)
      
      for (pair in genos_obs) {
        a1 <- pair[1]
        a2 <- pair[2]
        if (a1 == a2) {
          n_hom[a1] <- n_hom[a1] + 1L
        } else {
          n_het[a1] <- n_het[a1] + 1L
          n_het[a2] <- n_het[a2] + 1L
        }
      }

      # Initialize rd (p_null)
      if (n_nullhom > 0) {
        rd <- sqrt(n_nullhom / N)
      } else {
        rd <- sqrt(1 / (N + 1))
      }

      # Initialize cq[k] (corrdgenefreq)
      cq <- setNames(numeric(A), alleles_uniq)
      
      for (k in seq_len(A)) {
        a <- alleles_uniq[k]
        if (genefreq[k] <= 0) next
        
        ii <- n_hom[a]  # homozygotes a/a
        jj <- n_het[a]  # heterozygotes a/x with x != a
        
        if (n_nullhom > 0) {
          # FreeNA: cq[k] = 1 - sqrt((n_nullhom + N - ii - jj) / N)
          cq[k] <- 1 - sqrt((n_nullhom + N - ii - jj) / N)
        } else {
          # FreeNA: cq[k] = 1 - sqrt((1 + N - ii - jj) / (N + 1))
          cq[k] <- 1 - sqrt((1 + N - ii - jj) / (N + 1))
        }

        # Clamp to [0,1]
        cq[k] <- min(max(cq[k], 1e-10), 1 - 1e-10)
      }

      # EM loop
      cpt <- 0L
      old_rd <- rd
      old_cq <- cq
      converged <- FALSE

      repeat {
        rdi <- 0.0

        # Update corrected allele frequencies (cq)
        for (k in seq_len(A)) {
          a <- alleles_uniq[k]
          if (genefreq[k] <= 0) next

          ii <- n_hom[a]
          jj <- n_het[a]
          cq_old <- cq[k]

          # FreeNA update formula:
          # cq_new = (cq_old + rd)/(cq_old + 2*rd) * (ii/N) + jj/(2*N)
          denom <- cq_old + 2 * rd
          if (denom > 0) {
            cq[k] <- ((cq_old + rd) / denom) * (ii / N) + jj / (2 * N)
          } else {
            cq[k] <- cq_old
          }

          # Contribution to rdi
          # rdi += rd/(cq_old + 2*rd) * (ii/N)
          if (denom > 0) {
            rdi <- rdi + (rd / denom) * (ii / N)
          }

          # Clamp
          cq[k] <- min(max(cq[k], 1e-10), 1 - 1e-10)
        }

        # Update rd
        # FreeNA: rd = rdi + n_nullhom/N
        rd_new <- rdi + n_nullhom / N
        rd_new <- min(max(rd_new, 1e-10), 1 - 1e-10)

        # Check convergence
        rd_change <- abs(rd_new - rd)
        cq_change <- max(abs(cq - old_cq), na.rm = TRUE)

        rd <- rd_new
        cpt <- cpt + 1L

        if (rd_change < tol && cq_change < tol) {
          converged <- TRUE
          break
        }

        if (cpt >= max_iter) {
          warning("EM did not converge after ", max_iter, " iterations")
          break
        }

        old_rd <- rd
        old_cq <- cq
      }

      # He (observed genetic diversity)
      He <- 1 - sum(genefreq^2)

      list(p_null = rd, He = He, N = N, 
           N_missing = n_missing, N_nullhom = n_nullhom, 
           cpt = cpt, converged = converged)
    }

    # ════════════════════════════════════════════════════════════════════════
    # 5. ESTIMATION PER LOCUS AND SUB-POPULATION
    # ════════════════════════════════════════════════════════════════════════
    estimate_all <- function(fd, default_recode, locus_specific) {

      hap <- fd$haplotype
      pop_vector <- fd$pop_vector
      loci <- colnames(hap)
      pops <- sort(unique(pop_vector))

      get_codes <- function(locus) {
        choice <- locus_specific[[locus]] %||% "default"
        if (choice == "default") choice <- default_recode
        list(null_code = if (choice == "999999") "999999" else "0",
             miss_code = "0")
      }

      results_list <- list()

      # Debug output
      message("\n=== NULL ALLELE ESTIMATION ===")
      message("Loci: ", paste(loci, collapse = ", "))
      message("Populations: ", paste(pops, collapse = ", "))
      message("Default recode: ", default_recode)

      for (locus in loci) {
        codes <- get_codes(locus)
        message("\n--- Locus: ", locus, " ---")
        message("  null_code = ", codes$null_code, ", miss_code = ", codes$miss_code)

        pop_results <- list()

        for (pop in pops) {
          g_pop <- as.character(hap[pop_vector == pop, locus])
          n_total <- length(g_pop)
          n_non_na <- sum(!is.na(g_pop) & nzchar(trimws(g_pop)))

          message("  Population: ", pop, " (n=", n_total, ", non-NA=", n_non_na, ")")

          # Preview first genotypes
          if (n_non_na > 0) {
            preview <- head(g_pop[!is.na(g_pop) & nzchar(trimws(g_pop))], 5)
            message("    Preview: ", paste(preview, collapse = ", "))
          }

          em <- tryCatch(
            em_freena(g_pop, null_code = codes$null_code, miss_code = codes$miss_code),
            error = function(e) {
              warning("[NULL] EM error locus=", locus, " pop=", pop, ": ", e$message)
              list(p_null = NA_real_, He = NA_real_, N = 0L,
                   N_missing = NA_integer_, N_nullhom = 0L, 
                   cpt = 0L, converged = FALSE)
            }
          )

          message("    p_null = ", round(em$p_null, 6), 
                  ", He = ", round(em$He, 6),
                  ", N = ", em$N,
                  ", nullhom = ", em$N_nullhom,
                  ", missing = ", em$N_missing,
                  ", converged = ", em$converged)

          pop_results[[pop]] <- em

          results_list[[length(results_list) + 1]] <- data.frame(
            Locus = locus,
            Population = pop,
            p_null = em$p_null,
            He = em$He,
            N_total = length(g_pop),
            N_valid = em$N,
            N_missing = em$N_missing,
            N_nullhom = em$N_nullhom,
            Iterations = em$cpt,
            Converged = em$converged,
            stringsAsFactors = FALSE
          )
        }

        # Calculate inter-population weighted mean
        w_p_num <- 0
        w_p_den <- 0
        w_h_num <- 0
        w_h_den <- 0

        for (pop in pops) {
          em <- pop_results[[pop]]
          if (!is.na(em$p_null) && !is.na(em$He) && em$N > 0) {
            w <- em$N * em$He
            w_p_num <- w_p_num + em$p_null * w
            w_p_den <- w_p_den + w
          }
          if (!is.na(em$He) && em$N > 0) {
            w_h_num <- w_h_num + em$He * em$N
            w_h_den <- w_h_den + em$N
          }
        }

        results_list[[length(results_list) + 1]] <- data.frame(
          Locus = locus,
          Population = "WEIGHTED MEAN",
          p_null = if (w_p_den > 0) w_p_num / w_p_den else NA_real_,
          He = if (w_h_den > 0) w_h_num / w_h_den else NA_real_,
          N_total = sum(sapply(pop_results, function(x) x$N + x$N_missing), na.rm = TRUE),
          N_valid = sum(sapply(pop_results, function(x) x$N), na.rm = TRUE),
          N_missing = sum(sapply(pop_results, function(x) x$N_missing), na.rm = TRUE),
          N_nullhom = sum(sapply(pop_results, function(x) x$N_nullhom), na.rm = TRUE),
          Iterations = NA_integer_,
          Converged = NA,
          stringsAsFactors = FALSE
        )
      }

      message("\n=== ESTIMATION COMPLETE ===\n")
      do.call(rbind, results_list)
    }

    # ════════════════════════════════════════════════════════════════════════
    # 6. OBSERVER — "Estimate" button
    # ════════════════════════════════════════════════════════════════════════
    observeEvent(input$run_null_alleles, ignoreInit = TRUE, {
      tryCatch({
        fd <- formatted_data()
        if (is.null(fd)) {
          showNotification("Please import a dataset first.", type = "warning", duration = 6)
          return()
        }
        req(fd$haplotype)

        showNotification("EM algorithm running (FreeNA)...", type = "message", duration = 3)

        loci <- colnames(fd$haplotype)
        locus_specific <- list()
        for (locus in loci) {
          val <- input[[paste0("coding_", locus)]]
          if (!is.null(val) && val != "default") locus_specific[[locus]] <- val
        }

        results <- estimate_all(fd = fd,
                                default_recode = input$default_missing_recode,
                                locus_specific = locus_specific)

        null_allele_results(results)

        # Report
        mean_rows <- results[results$Population == "WEIGHTED MEAN", ]
        n_high <- sum(mean_rows$p_null > 0.20, na.rm = TRUE)
        n_converged <- sum(results$Converged == TRUE, na.rm = TRUE)

        showNotification(
          paste0("Complete — ", nrow(mean_rows), " loci | ",
                 "Converged: ", n_converged, " | ",
                 "p_null > 0.20: ", n_high),
          type = "message", duration = 6
        )
      }, error = function(e) {
        message("[NULL] ERROR: ", e$message)
        showNotification(paste("Error:", e$message), type = "error", duration = 10)
      })
    })

        # ════════════════════════════════════════════════════════════════════════
    # 7. DETAIL TABLE (per locus × sub-population)
    # ════════════════════════════════════════════════════════════════════════
    output$null_allele_detail_table <- DT::renderDT({
      req(null_allele_results())
      df <- null_allele_results()
      
      # Filter to keep only sub-population rows (not weighted mean)
      detail_df <- df[df$Population != "WEIGHTED MEAN", ]
      
      if(nrow(detail_df) == 0) {
        return(DT::datatable(data.frame(Message = "No detailed results available")))
      }
      
      # Formatting
      detail_df$p_null_fmt <- ifelse(is.na(detail_df$p_null), "—", sprintf("%.6f", detail_df$p_null))
      detail_df$He_fmt <- ifelse(is.na(detail_df$He), "—", sprintf("%.6f", detail_df$He))
      
      # Create Converged_fmt column for display
      detail_df$Converged_fmt <- ifelse(is.na(detail_df$Converged), "—", 
                                        ifelse(detail_df$Converged, "✓", "✗"))
      
      # Select columns for display
      display_df <- detail_df[, c("Locus", "Population", "p_null_fmt", "He_fmt",
                                   "N_total", "N_valid", "N_missing", "N_nullhom",
                                   "Iterations", "Converged_fmt")]
      
      # Rename the column for display purposes
      colnames(display_df)[colnames(display_df) == "Converged_fmt"] <- "Converged"
      
      DT::datatable(
        display_df,
        rownames = FALSE,
        colnames = c("Locus", "Sub-population", "p̂_null (EM)", "Ĥe",
                     "N total", "N valid", "N missing", "N null hom.",
                     "Iterations", "Converged"),
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv"),
          order = list(list(0, "asc"), list(1, "asc"))
        )
      ) %>%
        DT::formatStyle("Converged",
          backgroundColor = DT::styleEqual(c("✓", "✗"), c("#d4edda", "#f8d7da"))
        )
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # 8. SUMMARY TABLE (weighted means)
    # ════════════════════════════════════════════════════════════════════════
    output$null_allele_summary_table <- DT::renderDT({
      req(null_allele_results())
      df <- null_allele_results()
      
      # Filter to keep only weighted mean rows
      summary_df <- df[df$Population == "WEIGHTED MEAN", ]
      
      if(nrow(summary_df) == 0) {
        return(DT::datatable(data.frame(Message = "No summary results available")))
      }
      
      # Impact classification
      summary_df$Impact <- dplyr::case_when(
        is.na(summary_df$p_null) ~ "Undetermined",
        summary_df$p_null < 0.05 ~ "Negligible",
        summary_df$p_null < 0.10 ~ "Weak",
        summary_df$p_null < 0.20 ~ "Moderate",
        TRUE ~ "High"
      )
      
      summary_df$p_null_fmt <- ifelse(is.na(summary_df$p_null), "—", sprintf("%.6f", summary_df$p_null))
      summary_df$He_fmt <- ifelse(is.na(summary_df$He), "—", sprintf("%.6f", summary_df$He))
      
      display_df <- summary_df[, c("Locus", "p_null_fmt", "He_fmt", 
                                    "N_valid", "N_nullhom", "Impact")]
      
      DT::datatable(
        display_df,
        rownames = FALSE,
        colnames = c("Locus", "p̄_null (weighted)", "Ĥe mean",
                     "N total", "N null hom.", "Impact"),
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv"),
          order = list(list(1, "desc"))
        )
      ) %>%
        DT::formatStyle("Impact",
          backgroundColor = DT::styleEqual(
            c("Negligible", "Weak", "Moderate", "High", "Undetermined"),
            c("#d4edda", "#fff3cd", "#ffe5b4", "#f8d7da", "#e2e3e5")
          )
        )
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # 9. PLOTS — CORRECTED FOR SINGLE POPULATION CASE
    # ════════════════════════════════════════════════════════════════════════
    
    # p_null distribution
    output$null_allele_dist_plot <- renderPlot({
      req(null_allele_results())
      df <- null_allele_results()
      
      detail_df <- df[df$Population != "WEIGHTED MEAN" & !is.na(df$p_null), ]
      
      if(nrow(detail_df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data available", cex = 1.5)
        return()
      }
      
      n_pops <- length(unique(detail_df$Population))
      n_loci <- length(unique(detail_df$Locus))
      
      if(n_pops == 1) {
        # Single population - bar plot by locus
        p <- ggplot(detail_df, aes(x = Locus, y = p_null, fill = Locus)) +
          geom_bar(stat = "identity") +
          geom_hline(yintercept = c(0.05, 0.10, 0.20), linetype = "dashed", color = "red", alpha = 0.7) +
          geom_text(aes(label = sprintf("%.3f", p_null), y = p_null + 0.02), 
                    size = 3, angle = 45, hjust = 0) +
          scale_fill_viridis_d() +
          theme_minimal() +
          labs(title = "Null allele frequencies by locus",
               subtitle = paste("Population:", unique(detail_df$Population)[1]),
               x = "Locus", y = "p̂_null") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none")
        print(p)
      } else if(n_loci == 1) {
        # Single locus - distribution across populations
        p <- ggplot(detail_df, aes(x = Population, y = p_null, fill = Population)) +
          geom_bar(stat = "identity") +
          geom_hline(yintercept = c(0.05, 0.10, 0.20), linetype = "dashed", color = "red", alpha = 0.7) +
          geom_text(aes(label = sprintf("%.3f", p_null), y = p_null + 0.02), 
                    size = 3, angle = 45, hjust = 0) +
          scale_fill_viridis_d() +
          theme_minimal() +
          labs(title = paste("Null allele frequencies for locus:", unique(detail_df$Locus)[1]),
               x = "Population", y = "p̂_null") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none")
        print(p)
      } else {
        # Multiple populations and loci - faceted histogram
        p <- ggplot(detail_df, aes(x = p_null, fill = Population)) +
          geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
          geom_vline(xintercept = c(0.05, 0.10, 0.20), linetype = "dashed", color = "red", alpha = 0.7) +
          facet_wrap(~Locus, scales = "free_y") +
          theme_minimal() +
          labs(title = "Distribution of null allele frequencies by locus",
               x = "p̂_null", y = "Number of sub-populations") +
          theme(legend.position = "bottom",
                strip.text = element_text(face = "bold"))
        print(p)
      }
    })
    
    # Heatmap of p_null by locus × population — CORRECTED
    output$null_allele_heatmap <- renderPlot({
      req(null_allele_results())
      df <- null_allele_results()
      
      detail_df <- df[df$Population != "WEIGHTED MEAN" & !is.na(df$p_null), 
                      c("Locus", "Population", "p_null")]
      
      if(nrow(detail_df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data available", cex = 1.5)
        return()
      }
      
      # Check dimensions for heatmap
      n_pops <- length(unique(detail_df$Population))
      n_loci <- length(unique(detail_df$Locus))
      
      if(n_pops < 2 || n_loci < 2) {
        # Fallback to bar plot instead of heatmap
        if(n_pops == 1 && n_loci >= 1) {
          p <- ggplot(detail_df, aes(x = Locus, y = p_null, fill = Locus)) +
            geom_bar(stat = "identity") +
            geom_hline(yintercept = c(0.05, 0.10, 0.20), linetype = "dashed", color = "red", alpha = 0.7) +
            geom_text(aes(label = sprintf("%.3f", p_null), y = p_null + 0.02), 
                      size = 3, angle = 45, hjust = 0) +
            scale_fill_viridis_d() +
            theme_minimal() +
            labs(title = "Null allele frequencies by locus",
                 subtitle = paste("Population:", unique(detail_df$Population)[1]),
                 x = "Locus", y = "p̂_null") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none")
          print(p)
        } else if(n_loci < 2) {
          plot.new()
          text(0.5, 0.5, "Need at least 2 loci for heatmap visualization", cex = 1.2)
        }
        return()
      }
      
      # Create matrix for heatmap
      heatmap_data <- tidyr::pivot_wider(detail_df, 
                                          names_from = "Population", 
                                          values_from = "p_null",
                                          values_fill = NA)
      
      mat <- as.matrix(heatmap_data[, -1])
      rownames(mat) <- heatmap_data$Locus
      
      # Double-check matrix dimensions
      if(nrow(mat) < 2 || ncol(mat) < 2) {
        plot.new()
        text(0.5, 0.5, "Insufficient data for heatmap (need at least 2 loci and 2 populations)", 
             cex = 1.2)
        return()
      }
      
      # Use base R heatmap if gplots not available
      if(!requireNamespace("gplots", quietly = TRUE)) {
        heatmap(mat,
                main = "Null allele frequencies by locus and sub-population",
                xlab = "Sub-population", ylab = "Locus",
                col = grDevices::colorRampPalette(c("white", "yellow", "orange", "red"))(100),
                margins = c(8, 8))
        return()
      }
      
      # Use gplots::heatmap.2
      gplots::heatmap.2(mat,
                main = "Null allele frequencies by locus and sub-population",
                xlab = "Sub-population", ylab = "Locus",
                col = grDevices::colorRampPalette(c("white", "yellow", "orange", "red"))(100),
                margins = c(8, 8),
                key = TRUE,
                key.title = NA,
                key.xlab = "p̂_null",
                trace = "none",
                density.info = "none",
                cellnote = round(mat, 3),
                notecex = 0.7,
                notecol = "black")
    })
    
    # Convergence plot
    output$null_allele_convergence_plot <- renderPlot({
      req(null_allele_results())
      df <- null_allele_results()
      
      detail_df <- df[df$Population != "WEIGHTED MEAN", ]
      
      if(nrow(detail_df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data available", cex = 1.5)
        return()
      }
      
      # Summary by locus
      conv_summary <- detail_df %>%
        dplyr::group_by(Locus) %>%
        dplyr::summarise(
          n_pops = dplyr::n(),
          n_converged = sum(Converged == TRUE, na.rm = TRUE),
          pct_converged = 100 * n_converged / n_pops
        )
      
      if(nrow(conv_summary) == 0) {
        plot.new()
        text(0.5, 0.5, "No convergence data available", cex = 1.5)
        return()
      }
      
      p <- ggplot(conv_summary, aes(x = reorder(Locus, -pct_converged), y = pct_converged)) +
        geom_bar(stat = "identity", aes(fill = pct_converged == 100)) +
        geom_hline(yintercept = 100, linetype = "dashed", color = "red", alpha = 0.7) +
        scale_fill_manual(values = c("TRUE" = "#28a745", "FALSE" = "#ffc107"),
                          name = "100% converged") +
        theme_minimal() +
        labs(title = "EM algorithm convergence rate by locus",
             x = "Locus", y = "Convergence (%)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Add value labels if not too many loci
      if(nrow(conv_summary) <= 20) {
        p <- p + geom_text(aes(label = sprintf("%.0f%%", pct_converged), 
                               y = pct_converged + 2), 
                           size = 3)
      }
      
      print(p)
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # 10. STATISTICAL REPORT
    # ════════════════════════════════════════════════════════════════════════
    output$null_allele_summary <- renderPrint({
      req(null_allele_results())
      df <- null_allele_results()
      mean_rows <- df[df$Population == "WEIGHTED MEAN", ]
      pop_details <- unique(df$Population[df$Population != "WEIGHTED MEAN"])
      
      cat("=== NULL ALLELE FREQUENCIES ===\n")
      cat("Method: EM algorithm (Dempster, Laird & Rubin 1977)\n")
      cat("         FreeNA implementation (Chapuis & Estoup 2007)\n\n")
      
      cat("FreeNA convention:\n")
      cat("  • 999999/999999 = null homozygote (counted in EM)\n")
      cat("  • 0/0           = missing genotype (excluded)\n")
      cat("  • Null heterozygotes → forbidden\n\n")
      
      cat("Results by sub-population:\n")
      for(pop in pop_details) {
        pop_df <- df[df$Population == pop, ]
        cat(sprintf("  %s: %d loci analyzed\n", pop, nrow(pop_df)))
      }
      cat("\n")
      
      cat("--- Inter-population weighted means ---\n")
      cat("Formula: p̄_null = Σ(p_null_i × N_i × Ĥe_i) / Σ(N_i × Ĥe_i)\n\n")
      
      cat(sprintf("  Mean    : %.6f\n", mean(mean_rows$p_null, na.rm = TRUE)))
      cat(sprintf("  Median  : %.6f\n", median(mean_rows$p_null, na.rm = TRUE)))
      cat(sprintf("  SD      : %.6f\n", sd(mean_rows$p_null, na.rm = TRUE)))
      cat(sprintf("  Min     : %.6f\n", min(mean_rows$p_null, na.rm = TRUE)))
      cat(sprintf("  Max     : %.6f\n", max(mean_rows$p_null, na.rm = TRUE)))
      cat("\n")
      
      cat("--- Impact classification (p̂_null) ---\n")
      cat(sprintf("  Negligible (< 0.05)  : %d loci\n", 
                  sum(mean_rows$p_null < 0.05, na.rm = TRUE)))
      cat(sprintf("  Weak   (0.05–0.10)   : %d loci\n", 
                  sum(mean_rows$p_null >= 0.05 & mean_rows$p_null < 0.10, na.rm = TRUE)))
      cat(sprintf("  Moderate (0.10–0.20) : %d loci\n", 
                  sum(mean_rows$p_null >= 0.10 & mean_rows$p_null < 0.20, na.rm = TRUE)))
      cat(sprintf("  High    (≥ 0.20)     : %d loci\n", 
                  sum(mean_rows$p_null >= 0.20, na.rm = TRUE)))
      cat("\n")
      
      cat("--- EM algorithm convergence ---\n")
      conv_df <- df[df$Population != "WEIGHTED MEAN" & !is.na(df$Converged), ]
      cat(sprintf("  Converged : %d / %d (%.1f%%)\n",
                  sum(conv_df$Converged), nrow(conv_df),
                  if(nrow(conv_df) > 0) 100 * sum(conv_df$Converged) / nrow(conv_df) else 0))
      cat("\n")
      
      cat("--- Top 10 loci by p̂_null ---\n")
      top10 <- mean_rows[order(-mean_rows$p_null, na.last = TRUE), ]
      top10 <- head(top10, 10)
      if(nrow(top10) > 0) {
        print(top10[, c("Locus", "p_null", "He", "N_valid", "N_nullhom")], 
              row.names = FALSE, digits = 6)
      }
      cat("\n")
      
      cat("--- Recommendations ---\n")
      n_high <- sum(mean_rows$p_null > 0.20, na.rm = TRUE)
      if (n_high > 0) {
        cat("⚠️  WARNING:", n_high, "loci have p̂_null > 0.20.\n")
        cat("   → Use FST-ENA correction (Chapuis & Estoup 2007)\n")
        cat("   → Or exclude these loci from further analyses\n")
      } else {
        cat("✓ No locus with p̂_null > 0.20 detected.\n")
      }
      cat("\n")
      
      cat("--- References ---\n")
      cat("  Dempster AP, Laird NM, Rubin DB (1977). J R Stat Soc B 39:1-38.\n")
      cat("  Chapuis MP, Estoup A (2007). Mol Biol Evol 24:621-631.\n")
      cat("  FreeNA : http://www1.montpellier.inra.fr/URLB/softwares/freena/\n")
    })
    
    # ════════════════════════════════════════════════════════════════════════
    # 11. DOWNLOAD HANDLERS
    # ════════════════════════════════════════════════════════════════════════
    
    output$download_null_alleles_detail <- downloadHandler(
      filename = function() paste0("null_alleles_detail_", Sys.Date(), ".csv"),
      content = function(file) {
        req(null_allele_results())
        df <- null_allele_results()
        detail_df <- df[df$Population != "WEIGHTED MEAN", ]
        write.csv(detail_df, file, row.names = FALSE)
      }
    )
    
    output$download_null_alleles_summary <- downloadHandler(
      filename = function() paste0("null_alleles_summary_", Sys.Date(), ".csv"),
      content = function(file) {
        req(null_allele_results())
        df <- null_allele_results()
        summary_df <- df[df$Population == "WEIGHTED MEAN", ]
        write.csv(summary_df, file, row.names = FALSE)
      }
    )
    
    output$download_null_alleles_png <- downloadHandler(
      filename = function() paste0("null_alleles_plots_", Sys.Date(), ".png"),
      content = function(file) {
        req(null_allele_results())
        p1 <- local({
          df <- null_allele_results()
          mean_df <- df[df$Population == "WEIGHTED MEAN", ]
          req(nrow(mean_df) > 0)
          zpal <- colorRampPalette(c("#3B9AB2","#78B7C5","#EBCC2A","#E1AF00","#F21A00"))
          ggplot2::ggplot(mean_df, ggplot2::aes(x = p_null)) +
            ggplot2::geom_histogram(bins = 20, fill = "#3B9AB2", color = "white") +
            ggplot2::labs(title = "Distribution of null allele frequencies",
                          x = "p_null", y = "Count") +
            ggplot2::theme_minimal()
        })
        p2 <- local({
          df <- null_allele_results()
          pop_df <- df[df$Population != "WEIGHTED MEAN", ]
          req(nrow(pop_df) > 0)
          ggplot2::ggplot(pop_df, ggplot2::aes(x = Locus, y = Population, fill = p_null)) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::scale_fill_gradientn(
              colors = c("#3B9AB2","#78B7C5","#EBCC2A","#E1AF00","#F21A00"),
              name = "p_null") +
            ggplot2::labs(title = "Null allele frequency heatmap") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        })
        png(file, width = 1200, height = 1000, res = 120)
        gridExtra::grid.arrange(p1, p2, ncol = 1)
        dev.off()
      }
    )

    output$download_null_alleles_txt <- downloadHandler(
      filename = function() paste0("null_alleles_report_", Sys.Date(), ".txt"),
      content = function(file) {
        req(null_allele_results())
        sink(file)
        cat("NULL ALLELE FREQUENCY REPORT\n")
        cat("============================\n\n")
        cat("Generated: ", Sys.Date(), "\n\n")
        
        df <- null_allele_results()
        mean_rows <- df[df$Population == "WEIGHTED MEAN", ]
        
        cat("Summary statistics:\n")
        cat(sprintf("  Number of loci: %d\n", nrow(mean_rows)))
        cat(sprintf("  Mean p_null: %.6f\n", mean(mean_rows$p_null, na.rm = TRUE)))
        cat(sprintf("  Median p_null: %.6f\n", median(mean_rows$p_null, na.rm = TRUE)))
        cat(sprintf("  Loci with p_null > 0.20: %d\n", 
                    sum(mean_rows$p_null > 0.20, na.rm = TRUE)))
        sink()
      }
    )
    
  })
}