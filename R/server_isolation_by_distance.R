# server_null_alleles.R
# Calls engine_freena.R (faithful translation of FreeNA_optm2R.pas) for all
# genetic-distance statistics. Tabs 5-6 (Isolation by Distance, Mantel test)
# reuse the pairwise FST/FST-ENA/DCSE/DCSE-INA + bootstrap CI computed here —
# nothing is recomputed — following the SPG-V1 module specification.
#
# IMPORTANT: source engine_freena.R before this file (e.g. in app.R / global.R
# / help.R):
#   source("help.R")                 # contains engine_freena.R functions
#   source("ui_null_alleles.R")
#   source("server_null_alleles.R")

server_isolation_by_distance <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    `%||%` <- function(a, b) if (!is.null(a) && !(length(a) == 1 && is.na(a))) a else b

    # ── DB plumbing (same conventions as other modules) ────────────────────
    db_tick    <- reactive({ rv$db_tick })
    con_r      <- reactive({ shiny::req(rv$con); rv$con })
    tbl_meta_r <- reactive({ rv$tbl_meta %||% "meta" })

    db_ready <- reactive({
      db_tick(); con <- con_r()
      shiny::req(isTRUE(rv$db_ready))
      shiny::validate(shiny::need(DBI::dbExistsTable(con, tbl_meta_r()),
                                   "DuckDB meta table missing."))
      TRUE
    })

    # ── Raw genotypes -> hap_df (individuals x loci) + pop_vector ──────────
    raw_genos_r <- reactive({
      db_ready()
      con     <- con_r()
      tbl_raw <- rv$tbl_raw %||% "raw"
      shiny::validate(shiny::need(
        DBI::dbExistsTable(con, tbl_raw),
        "Raw genotype table not found. Please re-import the dataset."))

      ok_par <- tryCatch(DBI::dbExistsTable(con, "params"), error = function(e) FALSE)
      shiny::validate(shiny::need(ok_par, "params table not found."))

      marker_json <- tryCatch(
        DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='marker_cols_raw'")$value[1L],
        error = function(e) NA_character_)
      marker_cols_raw <- if (!is.na(marker_json) && nzchar(marker_json))
        tryCatch(jsonlite::fromJSON(marker_json), error = function(e) character(0))
      else character(0)
      shiny::validate(shiny::need(length(marker_cols_raw) > 0L, "No marker_cols_raw in params."))

      geno_fmt <- tryCatch(
        DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='genotype_format'")$value[1L],
        error = function(e) NA_character_)
      if (is.na(geno_fmt) || !nzchar(geno_fmt))
        geno_fmt <- if (any(grepl("(_1|\\.[0-9]+)$", marker_cols_raw))) "paired" else "string"

      keep     <- unique(marker_cols_raw)
      keep_sql <- paste(vapply(keep, function(x)
        as.character(DBI::dbQuoteIdentifier(con, x)), character(1L)), collapse = ", ")
      raw_df <- as.data.frame(
        DBI::dbGetQuery(con, sprintf(
          "SELECT rowid AS individual, %s FROM %s", keep_sql,
          as.character(DBI::dbQuoteIdentifier(con, tbl_raw)))),
        stringsAsFactors = FALSE)
      shiny::validate(shiny::need(nrow(raw_df) > 0L, "No rows in raw table."))

      meta_pop <- DBI::dbGetQuery(con, sprintf(
        "SELECT individual, Population FROM %s WHERE Population IS NOT NULL",
        sql_ident(con, tbl_meta_r())))
      raw_df$Population <- meta_pop$Population[match(raw_df$individual, meta_pop$individual)]
      pop_vector <- as.character(raw_df$Population)

      pick_b <- function(locus, nms) {
        cands <- c(paste0(locus, "_1"), paste0(locus, "_2"),
                   paste0(locus, ".", 1:9), paste0(locus, "_", 1:9))
        hit <- cands[cands %in% nms]; if (length(hit)) hit[1L] else NA_character_
      }

      if (identical(geno_fmt, "paired")) {
        nms  <- names(raw_df)
        loci <- unique(sub("(_1|_2|\\.[0-9]+)$", "", marker_cols_raw))
        hap_df <- data.frame(row.names = seq_len(nrow(raw_df)))
        for (locus in loci) {
          b <- pick_b(locus, nms)
          if (!locus %in% nms || is.na(b) || !b %in% nms) next
          a_v <- as.character(raw_df[[locus]])
          b_v <- as.character(raw_df[[b]])
          a_v[is.na(a_v) | trimws(a_v) == ""] <- "0"
          b_v[is.na(b_v) | trimws(b_v) == ""] <- "0"
          already <- grepl("/", a_v, fixed = TRUE) | grepl("-", a_v, fixed = TRUE)
          hap_df[[locus]] <- ifelse(already, a_v, paste0(a_v, "/", b_v))
        }
      } else {
        hap_df <- as.data.frame(raw_df[, marker_cols_raw, drop = FALSE],
                                 stringsAsFactors = FALSE)
        for (j in seq_along(hap_df)) {
          x <- as.character(hap_df[[j]])
          x[is.na(x) | trimws(x) == ""] <- "0/0"
          hap_df[[j]] <- x
        }
      }
      shiny::validate(shiny::need(ncol(hap_df) > 0L,
        "No locus columns could be reconstructed."))
      list(hap_df = hap_df, pop_vector = pop_vector, n_loci = ncol(hap_df))
    })

    # ── Population GPS centroids (optional, for IBD tab) ──────────────────
    coords_r <- reactive({
      db_ready()
      con  <- con_r()
      cols <- tryCatch(DBI::dbGetQuery(con, sprintf(
        "SELECT column_name FROM information_schema.columns WHERE table_name = '%s'",
        tbl_meta_r()))$column_name, error = function(e) character(0))
      if (!all(c("Latitude", "Longitude") %in% cols)) return(NULL)
      df <- tryCatch(DBI::dbGetQuery(con, sprintf(
        "SELECT Population,
                AVG(CAST(Latitude  AS DOUBLE)) AS Latitude,
                AVG(CAST(Longitude AS DOUBLE)) AS Longitude
         FROM %s
         WHERE Population IS NOT NULL
           AND Latitude IS NOT NULL AND Longitude IS NOT NULL
         GROUP BY Population ORDER BY Population",
        sql_ident(con, tbl_meta_r()))), error = function(e) NULL)
      if (is.null(df) || nrow(df) < 2L) return(NULL)
      df
    })

    # ── Populate per-locus / per-pop filter selectors ───────────────────────
    observe({
      rg <- tryCatch(raw_genos_r(), error = function(e) NULL)
      if (is.null(rg)) return(invisible(NULL))
      loci <- colnames(rg$hap_df)
      pops <- sort(unique(rg$pop_vector))
      updateSelectInput(session, "fl_locus", choices = c("All loci" = "all", stats::setNames(loci, loci)))
      updateSelectInput(session, "fl_pop1",  choices = c("All pairs" = "all", stats::setNames(pops, pops)))
      updateSelectInput(session, "fl_pop2",  choices = c("All pairs" = "all", stats::setNames(pops, pops)))
    })

    # ══════════════════════════════════════════════════════════════════════
    # MAIN COMPUTATION — engine_freena.R does all the work
    # ══════════════════════════════════════════════════════════════════════
    results_r <- eventReactive(input$run_all, {
      shiny::req(db_ready())
      rg <- raw_genos_r()
      shiny::validate(shiny::need(
        length(unique(rg$pop_vector)) >= 2L,
        "At least 2 populations are required."))

      null_code <- trimws(input$null_code)
      conf      <- input$conf_level / 100
      n_boot    <- as.integer(input$n_boot)

      withProgress(message = "Computing EM, Fst, Fst-ENA, DCSE, DCSE-INA\u2026", value = 0.1, {
        res <- .fr_compute_all(rg$hap_df, rg$pop_vector, null_code = null_code)
        setProgress(0.6, detail = "Bootstrap over loci\u2026")

        # Faithful to Pascal: skip bootstrap if nperm<100 or nloc<=4
        boot <- NULL
        boot_skipped_reason <- NULL
        if (n_boot >= 100L && res$nloc > 4L) {
          boot <- .fr_bootstrap_loci(res, n_boot = n_boot, conf = conf)
        } else if (n_boot > 0L) {
          boot_skipped_reason <- if (n_boot < 100L)
            "Bootstrap skipped: replicates must be >= 100 (Pascal source requirement)."
          else
            "Bootstrap skipped: more than 4 loci are required (Pascal source requirement)."
        }
        setProgress(1.0)
      })

      list(res = res, boot = boot, boot_skipped_reason = boot_skipped_reason,
           null_code = null_code, conf = conf, n_boot = n_boot)
    })

    # ── Status banner ────────────────────────────────────────────────────
    output$ui_run_status <- renderUI({
      r <- tryCatch(results_r(), error = function(e) NULL)
      if (is.null(r)) return(NULL)
      msgs <- list(
        tags$div(style = "margin-top:10px; padding:8px 10px; background:#eff6ff; border:1px solid #bfdbfe; border-radius:6px; font-size:11.5px; color:#1d4ed8;",
          icon("check-circle"), sprintf(" Done \u2014 %d loci, %d populations.", r$res$nloc, r$res$npop))
      )
      if (!is.null(r$boot_skipped_reason)) {
        msgs[[length(msgs) + 1L]] <- tags$div(
          style = "margin-top:6px; padding:8px 10px; background:#fffbeb; border:1px solid #fcd34d; border-radius:6px; font-size:11.5px; color:#92400e;",
          icon("exclamation-triangle"), " ", r$boot_skipped_reason)
      }
      tagList(msgs)
    })

    # ── Value boxes ──────────────────────────────────────────────────────
    output$box_nloci <- renderValueBox({
      valueBox(results_r()$res$nloc, "Loci", icon = icon("dna"), color = "navy")
    })
    output$box_npops <- renderValueBox({
      valueBox(results_r()$res$npop, "Populations", icon = icon("users"), color = "teal")
    })
    output$box_avgrd <- renderValueBox({
      r <- results_r()$res
      rds <- vapply(r$loci, function(lo) {
        vals <- vapply(r$pops, function(p) {
          e <- r$em_cache[[lo]][[p]]
          if (isTRUE(e$notappl)) NA_real_ else e$rd
        }, numeric(1))
        mean(vals, na.rm = TRUE)
      }, numeric(1))
      v <- round(mean(rds, na.rm = TRUE), 4)
      col <- if (is.na(v)) "navy" else if (v > 0.2) "red" else if (v > 0.1) "yellow" else "green"
      valueBox(if (is.na(v)) "NA" else v, HTML("Avg r<sub>d</sub>"),
               icon = icon("percentage"), color = col)
    })
    output$box_fstena <- renderValueBox({
      v <- round(results_r()$res$fst_global_ena, 4)
      col <- if (is.na(v)) "navy" else if (v > 0.15) "red" else if (v > 0.05) "yellow" else "green"
      valueBox(if (is.na(v)) "NA" else v, HTML("Global F<sub>ST</sub>-ENA"),
               icon = icon("chart-bar"), color = col)
    })

    # ══════════════════════════════════════════════════════════════════════
    # TAB 1 — null allele frequencies (rd) table
    # ══════════════════════════════════════════════════════════════════════
    rd_table_r <- reactive({
      r <- results_r()$res
      rows <- list(); k <- 1L
      for (lo in r$loci) {
        for (p in r$pops) {
          e <- r$em_cache[[lo]][[p]]
          pr <- r$parsed_cache[[lo]][[p]]
          rows[[k]] <- data.frame(
            Locus = lo, Population = p,
            rd        = if (isTRUE(e$notappl)) NA_real_ else round(e$rd, 6),
            N_total   = r$efpop[[p]],
            N_absent  = pr$n_absent + pr$n_nullhet,
            N_nullhomo= pr$n_nullhomo,
            N_valid   = pr$n_valid,
            stringsAsFactors = FALSE)
          k <- k + 1L
        }
      }
      do.call(rbind, rows)
    })

    output$dt_rd <- DT::renderDT({
      df <- rd_table_r()
      DT::datatable(df, rownames = FALSE,
        options = list(pageLength = 20, scrollX = TRUE, dom = "lftip"),
        class = "compact hover stripe") %>%
        DT::formatStyle("rd", backgroundColor = DT::styleInterval(
          c(0.05, 0.10, 0.20, 0.30),
          c("#f0fdf4", "#dcfce7", "#fefce8", "#fff7ed", "#fef2f2")))
    })
    output$dl_rd_csv <- downloadHandler(
      filename = function() paste0("null_allele_frequencies_", Sys.Date(), ".csv"),
      content  = function(file) write.csv(rd_table_r(), file, row.names = FALSE)
    )

    # ══════════════════════════════════════════════════════════════════════
    # TAB 2 — Fst / Fst-ENA
    # ══════════════════════════════════════════════════════════════════════

    output$ui_global_fst <- renderUI({
      rr <- results_r(); r <- rr$res; boot <- rr$boot
      ci_pct <- paste0(round(rr$conf * 100, 1), "%")
      ci_raw <- if (!is.null(boot)) sprintf("[ %.6f , %.6f ]", boot$global_raw_ci[1], boot$global_raw_ci[2]) else "NA (bootstrap not run)"
      ci_ena <- if (!is.null(boot)) sprintf("[ %.6f , %.6f ]", boot$global_ena_ci[1], boot$global_ena_ci[2]) else "NA (bootstrap not run)"
      tags$div(style = "font-family:monospace; font-size:13px; line-height:2;",
        tags$div(tags$strong("Raw F"), tags$sub("ST"), sprintf(": %.6f", r$fst_global_raw),
                 tags$br(), sprintf("%s CI (bootstrap over loci): %s", ci_pct, ci_raw)),
        tags$br(),
        tags$div(tags$strong("F"), tags$sub("ST"), "-ENA", sprintf(": %.6f", r$fst_global_ena),
                 tags$br(), sprintf("%s CI (bootstrap over loci): %s", ci_pct, ci_ena))
      )
    })

    # Continuous-gradient colored half-matrix
    .render_mat_html <- function(mat, digits = 4, thr = c(0.05, 0.15, 0.25),
                                  clrs = c("#f0fdf4", "#dcfce7", "#fefce8", "#fef2f2")) {
      labs <- rownames(mat); n <- length(labs)
      cell <- function(i, j) {
        if (i == j) return('<td style="background:#f1f5f9;color:#94a3b8;text-align:center;padding:4px 9px;">\u2014</td>')
        if (i < j)  return('<td style="color:#cbd5e1;text-align:center;padding:4px 9px;">\u00b7</td>')
        v <- mat[i, j]
        if (!is.finite(v)) return('<td style="color:#94a3b8;text-align:center;padding:4px 9px;">NA</td>')
        bg <- clrs[findInterval(v, thr) + 1L]
        sprintf('<td style="background:%s;text-align:right;padding:4px 9px;">%s</td>', bg, round(v, digits))
      }
      thead <- paste0('<tr><th></th>', paste(sprintf('<th style="padding:4px 9px;">%s</th>', labs[-n]), collapse = ""), '</tr>')
      tbody <- paste(sapply(seq_len(n), function(i) {
        if (i == 1L) return("")
        paste0('<tr><td style="font-weight:700;white-space:nowrap;padding:4px 9px;">', labs[i], '</td>',
               paste(sapply(seq_len(n), function(j) cell(i, j)), collapse = ""), '</tr>')
      }), collapse = "")
      HTML(sprintf('<div style="overflow-x:auto;"><table style="border-collapse:collapse;font-size:11px;width:100%%;"><thead>%s</thead><tbody>%s</tbody></table></div>', thead, tbody))
    }

    .build_pair_matrix <- function(r, value_col) {
      pops <- r$pops; n <- length(pops)
      m <- matrix(NA_real_, n, n, dimnames = list(pops, pops))
      for (k in seq_len(nrow(r$pair_df))) {
        i <- which(pops == r$pair_df$Pop1[k]); j <- which(pops == r$pair_df$Pop2[k])
        m[j, i] <- r$pair_df[[value_col]][k]
      }
      m
    }

    output$ui_fst_matrix <- renderUI({
      r <- results_r()$res
      typ <- input$fst_display %||% "both"
      m_raw <- .build_pair_matrix(r, "FST_raw")
      m_ena <- .build_pair_matrix(r, "FST_ENA")
      if (identical(typ, "both"))
        tags$div(tags$p(tags$strong("Raw")), .render_mat_html(m_raw), tags$br(),
                 tags$p(tags$strong("ENA-corrected")), .render_mat_html(m_ena))
      else if (identical(typ, "raw")) .render_mat_html(m_raw)
      else .render_mat_html(m_ena)
    })

    fst_pair_ci_r <- reactive({
      rr <- results_r(); r <- rr$res; boot <- rr$boot
      df <- r$pair_df[, c("Pop1", "Pop2", "FST_raw", "FST_ENA")]
      if (!is.null(boot)) {
        df$CI_lo_raw <- round(boot$pair_raw_ci[, 1], 6); df$CI_hi_raw <- round(boot$pair_raw_ci[, 2], 6)
        df$CI_lo_ena <- round(boot$pair_ena_ci[, 1], 6); df$CI_hi_ena <- round(boot$pair_ena_ci[, 2], 6)
      }
      df$FST_raw <- round(df$FST_raw, 6); df$FST_ENA <- round(df$FST_ENA, 6)
      df
    })
    output$dt_fst_pair_ci <- DT::renderDT({
      DT::datatable(fst_pair_ci_r(), rownames = FALSE,
        options = list(pageLength = 25, scrollX = TRUE, dom = "lftip"),
        class = "compact hover stripe")
    })
    output$dl_fst_pair_csv <- downloadHandler(
      filename = function() paste0("pairwise_FST_FST-ENA_", Sys.Date(), ".csv"),
      content  = function(file) write.csv(fst_pair_ci_r(), file, row.names = FALSE)
    )

    # ══════════════════════════════════════════════════════════════════════
    # TAB 3 — DCSE / DCSE-INA
    # ══════════════════════════════════════════════════════════════════════
    output$ui_dc_matrix <- renderUI({
      r <- results_r()$res
      typ <- input$dc_display %||% "both"
      thr <- c(0.1, 0.25, 0.4); clrs <- c("#eff6ff", "#dbeafe", "#fef9c3", "#fef2f2")
      m_raw <- .build_pair_matrix(r, "DCSE_raw")
      m_ina <- .build_pair_matrix(r, "DCSE_INA")
      if (identical(typ, "both"))
        tags$div(tags$p(tags$strong("Raw")), .render_mat_html(m_raw, thr = thr, clrs = clrs), tags$br(),
                 tags$p(tags$strong("INA-corrected")), .render_mat_html(m_ina, thr = thr, clrs = clrs))
      else if (identical(typ, "raw")) .render_mat_html(m_raw, thr = thr, clrs = clrs)
      else .render_mat_html(m_ina, thr = thr, clrs = clrs)
    })

    dc_pair_ci_r <- reactive({
      rr <- results_r(); r <- rr$res; boot <- rr$boot
      df <- r$pair_df[, c("Pop1", "Pop2", "DCSE_raw", "DCSE_INA")]
      if (!is.null(boot)) {
        df$CI_lo_raw <- round(boot$dc_raw_ci[, 1], 6); df$CI_hi_raw <- round(boot$dc_raw_ci[, 2], 6)
        df$CI_lo_ina <- round(boot$dc_ena_ci[, 1], 6); df$CI_hi_ina <- round(boot$dc_ena_ci[, 2], 6)
      }
      df$DCSE_raw <- round(df$DCSE_raw, 6); df$DCSE_INA <- round(df$DCSE_INA, 6)
      df
    })
    output$dt_dc_pair_ci <- DT::renderDT({
      DT::datatable(dc_pair_ci_r(), rownames = FALSE,
        options = list(pageLength = 25, scrollX = TRUE, dom = "lftip"),
        class = "compact hover stripe")
    })
    output$dl_dc_pair_csv <- downloadHandler(
      filename = function() paste0("pairwise_DCSE_DCSE-INA_", Sys.Date(), ".csv"),
      content  = function(file) write.csv(dc_pair_ci_r(), file, row.names = FALSE)
    )

    # ══════════════════════════════════════════════════════════════════════
    # TAB 4 — Per-locus x pair
    # ══════════════════════════════════════════════════════════════════════
    locus_pair_table_r <- reactive({
      r <- results_r()$res
      npairs <- nrow(r$pair_df)
      rows <- vector("list", npairs * r$nloc); k <- 1L
      for (pi in seq_len(npairs)) {
        for (li in seq_len(r$nloc)) {
          a1 <- r$w_s1_raw_pair[pi, li]; a3 <- r$w_s3_raw_pair[pi, li]
          a1e<- r$w_s1_ena_pair[pi, li]; a3e<- r$w_s3_ena_pair[pi, li]
          rows[[k]] <- data.frame(
            Locus = r$loci[li], Pop1 = r$pair_df$Pop1[pi], Pop2 = r$pair_df$Pop2[pi],
            FST_raw  = if (a3 != 0)  round(a1 / a3, 6)   else NA_real_,
            FST_ENA  = if (a3e != 0) round(a1e / a3e, 6) else NA_real_,
            DCSE_raw = round(r$dc_raw_pair[pi, li], 6),
            DCSE_INA = round(r$dc_ena_pair[pi, li], 6),
            stringsAsFactors = FALSE)
          k <- k + 1L
        }
      }
      do.call(rbind, rows)
    })

    output$dt_locus_pair <- DT::renderDT({
      df <- locus_pair_table_r()
      sl <- input$fl_locus %||% "all"; sp1 <- input$fl_pop1 %||% "all"; sp2 <- input$fl_pop2 %||% "all"
      if (!identical(sl, "all"))  df <- df[df$Locus == sl, , drop = FALSE]
      if (!identical(sp1, "all")) df <- df[df$Pop1 == sp1 | df$Pop2 == sp1, , drop = FALSE]
      if (!identical(sp2, "all")) df <- df[df$Pop2 == sp2 | df$Pop1 == sp2, , drop = FALSE]
      shiny::validate(shiny::need(nrow(df) > 0L, "No rows for selected filters."))
      DT::datatable(df, rownames = FALSE,
        options = list(pageLength = 25, scrollX = TRUE, dom = "lftip"),
        class = "compact hover stripe")
    })
    output$dl_locus_pair_csv <- downloadHandler(
      filename = function() paste0("per_locus_pair_", Sys.Date(), ".csv"),
      content  = function(file) write.csv(locus_pair_table_r(), file, row.names = FALSE)
    )

    # ══════════════════════════════════════════════════════════════════════
    # SHARED HELPER — full pairwise table (FST, FST-ENA, DCSE, DCSE-INA, FR,
    # FR-ENA + bootstrap CI bounds, linearised) + Dgeo/lnDgeo if GPS present.
    # Used by BOTH the Isolation by Distance tab and the Mantel tab.
    # ══════════════════════════════════════════════════════════════════════

    .linearise <- function(x) { x <- pmin(pmax(x, 0), 0.9999); x / (1 - x) }

    .haversine_km <- function(lat1, lon1, lat2, lon2) {
      R <- 6371.0
      dlat <- (lat2 - lat1) * pi / 180; dlon <- (lon2 - lon1) * pi / 180
      a <- sin(dlat / 2)^2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2)^2
      2 * R * asin(sqrt(a))
    }

    full_pair_table_r <- reactive({
      rr <- results_r(); r <- rr$res; boot <- rr$boot
      df <- r$pair_df[, c("Pop1", "Pop2", "FST_raw", "FST_ENA", "DCSE_raw", "DCSE_INA")]

      if (!is.null(boot)) {
        df$FST_raw_lo <- boot$pair_raw_ci[, 1]; df$FST_raw_hi <- boot$pair_raw_ci[, 2]
        df$FST_ENA_lo <- boot$pair_ena_ci[, 1]; df$FST_ENA_hi <- boot$pair_ena_ci[, 2]
      } else {
        df$FST_raw_lo <- NA_real_; df$FST_raw_hi <- NA_real_
        df$FST_ENA_lo <- NA_real_; df$FST_ENA_hi <- NA_real_
      }

      # Rousset's FR = FST/(1-FST); CI bounds linearised too (monotonic transform)
      df$FR     <- .linearise(df$FST_raw)
      df$FR_lo  <- .linearise(df$FST_raw_lo)
      df$FR_hi  <- .linearise(df$FST_raw_hi)
      df$FR_ENA    <- .linearise(df$FST_ENA)
      df$FR_ENA_lo <- .linearise(df$FST_ENA_lo)
      df$FR_ENA_hi <- .linearise(df$FST_ENA_hi)

      # Geographic distance, if GPS available
      coords <- tryCatch(coords_r(), error = function(e) NULL)
      if (!is.null(coords)) {
        get_d <- function(p1, p2) {
          c1 <- coords[coords$Population == p1, ]; c2 <- coords[coords$Population == p2, ]
          if (nrow(c1) >= 1L && nrow(c2) >= 1L)
            .haversine_km(c1$Latitude[1L], c1$Longitude[1L], c2$Latitude[1L], c2$Longitude[1L])
          else NA_real_
        }
        df$Dgeo_km <- mapply(get_d, df$Pop1, df$Pop2)
        df$lnDgeo  <- ifelse(df$Dgeo_km > 0, log(df$Dgeo_km), NA_real_)
      } else {
        df$Dgeo_km <- NA_real_; df$lnDgeo <- NA_real_
      }

      df
    })

    # ══════════════════════════════════════════════════════════════════════
    # TAB 5 — Isolation by Distance (Rousset 1997)
    # ══════════════════════════════════════════════════════════════════════

    .fit_line <- function(y, x) {
      ok <- is.finite(y) & is.finite(x)
      if (sum(ok) < 3L) return(list(slope = NA_real_, intercept = NA_real_, r2 = NA_real_))
      m <- lm(y[ok] ~ x[ok])
      list(slope = unname(coef(m)[2L]), intercept = unname(coef(m)[1L]), r2 = summary(m)$r.squared)
    }

    ibd_results_r <- eventReactive(input$run_ibd, {
      df <- full_pair_table_r()
      shiny::validate(shiny::need(
        any(is.finite(df$Dgeo_km)),
        "No geographic distances available. Set Latitude/Longitude at import for at least 2 populations."))

      use_log <- identical(input$ibd_model, "2D")
      x <- if (use_log) df$lnDgeo else df$Dgeo_km
      x_label <- if (use_log) "ln(Dgeo)" else "Dgeo (km)"

      if (identical(input$ibd_metric, "ena")) {
        y_avg <- df$FR_ENA; y_lo <- df$FR_ENA_lo; y_hi <- df$FR_ENA_hi
        y_label <- "FR-ENA"
      } else {
        y_avg <- df$FR; y_lo <- df$FR_lo; y_hi <- df$FR_hi
        y_label <- "FR"
      }

      reg_avg <- .fit_line(y_avg, x)
      reg_lo  <- .fit_line(y_lo,  x)
      reg_hi  <- .fit_line(y_hi,  x)

      list(df = df, x = x, y_avg = y_avg, y_lo = y_lo, y_hi = y_hi,
           x_label = x_label, y_label = y_label,
           reg_avg = reg_avg, reg_lo = reg_lo, reg_hi = reg_hi,
           use_log = use_log, metric = input$ibd_metric)
    })

    output$dt_ibd_reg <- DT::renderDT({
      r <- ibd_results_r()
      fmt <- function(x) if (is.na(x)) "NA" else formatC(x, format = "f", digits = 6)
      tab <- data.frame(
        Line      = c(paste0(r$y_label, " (point estimate)"),
                       paste0(r$y_label, "-l (lower CI)"),
                       paste0(r$y_label, "-u (upper CI)")),
        Slope     = c(fmt(r$reg_avg$slope), fmt(r$reg_lo$slope), fmt(r$reg_hi$slope)),
        Intercept = c(fmt(r$reg_avg$intercept), fmt(r$reg_lo$intercept), fmt(r$reg_hi$intercept)),
        R2        = c(fmt(r$reg_avg$r2), fmt(r$reg_lo$r2), fmt(r$reg_hi$r2)),
        stringsAsFactors = FALSE
      )
      DT::datatable(tab, rownames = FALSE,
        options = list(dom = "t", pageLength = 3, ordering = FALSE),
        class = "compact stripe")
    })

    output$ui_ibd_interpretation <- renderUI({
      r <- ibd_results_r()
      slopes <- c(r$reg_avg$slope, r$reg_lo$slope, r$reg_hi$slope)
      if (any(is.na(slopes))) {
        return(tags$div(class = "spg-method-note", style = "border-left-color:#999;",
          "Could not fit all three regression lines (insufficient valid pairs)."))
      }
      all_pos <- all(slopes > 0)
      lo_neg  <- r$reg_lo$slope < 0 && r$reg_avg$slope > 0 && r$reg_hi$slope > 0
      if (all_pos) {
        tags$div(style = "padding:10px; background:#dcfce7; border:1px solid #86efac; border-radius:6px; color:#166534; font-size:13px;",
          icon("check-circle"), tags$strong(" All three slopes are positive: "),
          "this supports isolation by distance.")
      } else if (lo_neg) {
        tags$div(style = "padding:10px; background:#fffbeb; border:1px solid #fcd34d; border-radius:6px; color:#92400e; font-size:13px;",
          icon("exclamation-triangle"), tags$strong(" Lower-bound slope is negative: "),
          "this may indicate low power of the per-pair bootstrap rather than a true absence of IBD. ",
          "Consider running the Mantel test (next tab), ideally with DCSE, to confirm.")
      } else {
        tags$div(style = "padding:10px; background:#fef2f2; border:1px solid #fca5a5; border-radius:6px; color:#991b1b; font-size:13px;",
          icon("times-circle"), tags$strong(" No consistent positive trend: "),
          "no clear evidence of isolation by distance with this dataset/model.")
      }
    })

    output$ibd_plot <- plotly::renderPlotly({
      r <- ibd_results_r()
      shiny::req(any(is.finite(r$x)))
      x_seq <- seq(min(r$x, na.rm = TRUE), max(r$x, na.rm = TRUE), length.out = 100)
      mkline <- function(reg) if (is.na(reg$slope)) NULL else
        data.frame(x = x_seq, y = reg$intercept + reg$slope * x_seq)
      l_avg <- mkline(r$reg_avg); l_lo <- mkline(r$reg_lo); l_hi <- mkline(r$reg_hi)

      p <- plotly::plot_ly() %>%
        plotly::add_segments(x = ~r$x, xend = ~r$x, y = ~r$y_lo, yend = ~r$y_hi,
          line = list(color = "rgba(100,100,100,0.35)", width = 1),
          showlegend = FALSE, hoverinfo = "none") %>%
        plotly::add_markers(
          x = r$x, y = r$y_avg,
          text = paste0(r$df$Pop1, " \u2013 ", r$df$Pop2),
          hoverinfo = "text",
          marker = list(color = "#2CBF9F", size = 7, opacity = 0.85),
          name = "Pairs")

      if (!is.null(l_avg)) p <- p %>% plotly::add_lines(data = l_avg, x = ~x, y = ~y,
        line = list(color = "#333a43", width = 2), name = sprintf("Avg b=%.4f", r$reg_avg$slope))
      if (!is.null(l_lo)) p <- p %>% plotly::add_lines(data = l_lo, x = ~x, y = ~y,
        line = list(color = "#3B9AB2", width = 1.5, dash = "dot"), name = sprintf("Lower CI b=%.4f", r$reg_lo$slope))
      if (!is.null(l_hi)) p <- p %>% plotly::add_lines(data = l_hi, x = ~x, y = ~y,
        line = list(color = "#B40F20", width = 1.5, dash = "dash"), name = sprintf("Upper CI b=%.4f", r$reg_hi$slope))

      p %>% plotly::layout(
        xaxis = list(title = r$x_label), yaxis = list(title = r$y_label),
        legend = list(x = 0.02, y = 0.98, bgcolor = "rgba(255,255,255,0.8)"),
        margin = list(t = 30))
    })

    output$dt_ibd_table <- DT::renderDT({
      r <- ibd_results_r()
      df <- r$df[, c("Pop1", "Pop2", "Dgeo_km", "lnDgeo", "FR", "FR_lo", "FR_hi", "FR_ENA", "FR_ENA_lo", "FR_ENA_hi")]
      num_cols <- setdiff(names(df), c("Pop1", "Pop2"))
      df[num_cols] <- lapply(df[num_cols], round, 5)
      DT::datatable(df, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 10, dom = "lrtip"),
        class = "compact stripe hover")
    })
    output$dl_ibd_csv <- downloadHandler(
      filename = function() paste0("IBD_pairwise_", Sys.Date(), ".csv"),
      content  = function(file) write.csv(ibd_results_r()$df, file, row.names = FALSE)
    )

    # ══════════════════════════════════════════════════════════════════════
    # TAB 6 — Mantel test (joint row/column permutation; rectangular-safe)
    # ══════════════════════════════════════════════════════════════════════

    .mt_build_square <- function(df, id1, id2, value_col, all_labels) {
      n <- length(all_labels)
      m <- matrix(NA_real_, n, n, dimnames = list(all_labels, all_labels))
      for (k in seq_len(nrow(df))) {
        i <- as.character(df[[id1]][k]); j <- as.character(df[[id2]][k]); v <- df[[value_col]][k]
        if (i %in% all_labels && j %in% all_labels && is.finite(v)) { m[i, j] <- v; m[j, i] <- v }
      }
      m
    }

    .mt_mantel_matrix <- function(mat1, mat2, n_perm = 9999L, stat = "r") {
      common <- intersect(rownames(mat1), rownames(mat2))
      if (length(common) < 3L)
        return(list(stat_obs = NA_real_, p_pos = NA_real_, p_neg = NA_real_, n_pairs = 0L,
                    slope = NA_real_, intercept = NA_real_, r2 = NA_real_,
                    x = numeric(0), y = numeric(0), common = common, perm_stats = numeric(0)))
      m1 <- mat1[common, common, drop = FALSE]; m2 <- mat2[common, common, drop = FALSE]
      n  <- length(common)
      lower_idx <- which(lower.tri(matrix(TRUE, n, n)))
      x_all <- m1[lower_idx]; y_all <- m2[lower_idx]
      stat_fn <- function(xx, yy) {
        ok <- is.finite(xx) & is.finite(yy)
        if (sum(ok) < 3L) return(NA_real_)
        if (stat == "b") unname(coef(lm(yy[ok] ~ xx[ok]))[2L]) else suppressWarnings(cor(xx[ok], yy[ok]))
      }
      ok_obs   <- is.finite(x_all) & is.finite(y_all)
      stat_obs <- stat_fn(x_all, y_all)
      perm_stats <- vapply(seq_len(n_perm), function(.b) {
        perm <- sample.int(n); m2p <- m2[perm, perm, drop = FALSE]
        stat_fn(x_all, m2p[lower_idx])
      }, numeric(1L))
      perm_fin <- perm_stats[is.finite(perm_stats)]
      p_pos <- if (length(perm_fin) > 0L && is.finite(stat_obs)) mean(perm_fin >= stat_obs) else NA_real_
      lm0 <- tryCatch(lm(y_all[ok_obs] ~ x_all[ok_obs]), error = function(e) NULL)
      list(stat_obs = stat_obs, p_pos = p_pos, p_neg = 1 - p_pos, n_pairs = sum(ok_obs),
           slope = if (!is.null(lm0)) unname(coef(lm0)[2L]) else NA_real_,
           intercept = if (!is.null(lm0)) unname(coef(lm0)[1L]) else NA_real_,
           r2 = if (!is.null(lm0)) summary(lm0)$r.squared else NA_real_,
           x = x_all[ok_obs], y = y_all[ok_obs], common = common, perm_stats = perm_fin)
    }

    .mt_read_file <- function(fileinfo, sep, header) {
      df <- tryCatch(read.table(fileinfo$datapath, header = header, sep = sep,
                                stringsAsFactors = FALSE, check.names = FALSE,
                                fill = TRUE, quote = "\""),
                     error = function(e) NULL)
      shiny::validate(shiny::need(!is.null(df) && nrow(df) >= 3L,
        "Could not parse the file. Check separator / header settings."))
      df
    }

    mt_base_df_r <- reactive({
      if (input$mt_source == "internal") {
        df <- full_pair_table_r()
        if (isTRUE(input$mt_use_extra)) {
          shiny::req(input$mt_extra_file)
          extra <- .mt_read_file(input$mt_extra_file, input$mt_extra_sep, input$mt_extra_header)
          shiny::validate(shiny::need(ncol(extra) >= 3L,
            "Extra file must have 2 ID columns + at least 1 distance column."))
          id_cols  <- names(extra)[1:2]
          val_cols <- setdiff(names(extra), id_cols)
          extra_keep <- extra[, val_cols, drop = FALSE]
          key <- function(a, b) { a<-as.character(a); b<-as.character(b); ifelse(a<=b, paste(a,b,sep="__"), paste(b,a,sep="__")) }
          extra_keep$.key <- key(extra[[1L]], extra[[2L]])
          extra_keep <- extra_keep[!duplicated(extra_keep$.key), , drop = FALSE]
          df$.key <- key(df$Pop1, df$Pop2)
          df <- merge(df, extra_keep, by = ".key", all.x = TRUE, sort = FALSE)
          df$.key <- NULL
        }
        df
      } else {
        shiny::req(input$mt_file)
        .mt_read_file(input$mt_file, input$mt_sep, input$mt_header)
      }
    })

    .guess_col <- function(cols, patterns, fallback) {
      for (pat in patterns) { hit <- grep(pat, cols, value = TRUE, ignore.case = TRUE); if (length(hit)) return(hit[1L]) }
      fallback
    }

    output$mt_col_pop1_ui <- renderUI({
      cols <- tryCatch(names(mt_base_df_r()), error = function(e) character(0))
      selectInput(session$ns("mt_col_pop1"), "Population 1 column:", choices = cols,
                  selected = .guess_col(cols, c("^Pop1$"), cols[1]))
    })
    output$mt_col_pop2_ui <- renderUI({
      cols <- tryCatch(names(mt_base_df_r()), error = function(e) character(0))
      selectInput(session$ns("mt_col_pop2"), "Population 2 column:", choices = cols,
                  selected = .guess_col(cols, c("^Pop2$"), cols[min(2L, length(cols))]))
    })
    output$mt_col_x_ui <- renderUI({
      df <- tryCatch(mt_base_df_r(), error = function(e) NULL)
      cols <- if (is.null(df)) character(0) else names(df)[sapply(df, is.numeric)]
      selectInput(session$ns("mt_col_x"), "X column:", choices = cols,
                  selected = .guess_col(cols, c("Dgeo", "lnDgeo"), if (length(cols)) cols[1] else NULL))
    })
    output$mt_col_y_ui <- renderUI({
      df <- tryCatch(mt_base_df_r(), error = function(e) NULL)
      cols <- if (is.null(df)) character(0) else names(df)[sapply(df, is.numeric)]
      selectInput(session$ns("mt_col_y"), "Y column:", choices = cols,
                  selected = .guess_col(cols, c("^FR_ENA$", "^FR$", "FST_ENA", "DCSE"),
                                        if (length(cols) >= 2L) cols[2] else NULL))
    })

    mantel_result_r <- eventReactive(input$run_mantel, {
      df <- mt_base_df_r()
      shiny::req(input$mt_col_pop1, input$mt_col_pop2, input$mt_col_x, input$mt_col_y)
      p1c <- input$mt_col_pop1; p2c <- input$mt_col_pop2; xcol <- input$mt_col_x; ycol <- input$mt_col_y

      shiny::validate(
        shiny::need(all(c(p1c, p2c, xcol, ycol) %in% names(df)), "Selected columns not found."),
        shiny::need(p1c != p2c, "Population 1 and 2 must differ."),
        shiny::need(xcol != ycol, "X and Y must differ.")
      )

      if (nzchar(trimws(input$mt_exclude %||% ""))) {
        excl <- trimws(strsplit(input$mt_exclude, ",")[[1L]]); excl <- excl[nzchar(excl)]
        if (length(excl)) {
          key <- function(a,b){a<-as.character(a);b<-as.character(b);ifelse(a<=b,paste(a,b,sep="__"),paste(b,a,sep="__"))}
          key_df <- key(df[[p1c]], df[[p2c]])
          key_excl <- vapply(excl, function(s) {
            ids <- trimws(strsplit(s, "-")[[1L]]); if (length(ids) == 2L) key(ids[1], ids[2]) else NA_character_
          }, character(1L))
          df <- df[!(key_df %in% key_excl), , drop = FALSE]
        }
      }

      x <- suppressWarnings(as.numeric(df[[xcol]]))
      y <- suppressWarnings(as.numeric(df[[ycol]]))
      if (isTRUE(input$mt_log_x)) x <- ifelse(x > 0, log(x), NA_real_)

      all_labels <- sort(unique(c(as.character(df[[p1c]]), as.character(df[[p2c]]))))
      tmp <- data.frame(P1 = as.character(df[[p1c]]), P2 = as.character(df[[p2c]]), X = x, Y = y)
      m_x <- .mt_build_square(tmp, "P1", "P2", "X", all_labels)
      m_y <- .mt_build_square(tmp, "P1", "P2", "Y", all_labels)

      n_perm <- as.integer(input$mt_n_perm); stat <- input$mt_stat
      withProgress(message = "Running Mantel test\u2026", value = 0.2, {
        res <- .mt_mantel_matrix(m_y, m_x, n_perm = n_perm, stat = stat)
        setProgress(1.0)
      })
      res$x_label <- paste0(xcol, if (isTRUE(input$mt_log_x)) " (ln)" else "")
      res$y_label <- ycol
      res$stat_label <- if (stat == "b") "Slope b" else "Pearson r"
      res
    })

    output$box_m_stat <- renderValueBox({
      r <- mantel_result_r()
      valueBox(round(r$stat_obs, 4), HTML(paste0(r$stat_label, "<br>(observed)")),
               icon = icon("chart-line"), color = "purple")
    })
    output$box_m_pval <- renderValueBox({
      r <- mantel_result_r(); pv <- r$p_pos
      col <- if (is.na(pv)) "yellow" else if (pv < 0.05) "green" else if (pv < 0.10) "yellow" else "red"
      valueBox(if (is.na(pv)) "NA" else formatC(pv, format = "f", digits = 4),
               HTML("p-value<br>(one-sided)"), icon = icon("check-circle"), color = col)
    })
    output$box_m_n <- renderValueBox({
      valueBox(mantel_result_r()$n_pairs, "Pairs used", icon = icon("project-diagram"), color = "blue")
    })
    output$box_m_r2 <- renderValueBox({
      r2 <- mantel_result_r()$r2
      valueBox(if (is.na(r2)) "NA" else paste0(round(r2 * 100, 1), "%"),
               HTML("Variance<br>explained (R\u00b2)"), icon = icon("percentage"), color = "teal")
    })

    output$ui_mantel_summary <- renderUI({
      r <- mantel_result_r()
      tags$div(style = "margin-top:8px; font-family:monospace; font-size:12px; color:#555;",
        sprintf("Slope = %.6f, Intercept = %.6f", r$slope, r$intercept), tags$br(),
        sprintf("One-sided p (negative association) = %s",
                if (is.na(r$p_neg)) "NA" else formatC(r$p_neg, format = "f", digits = 4)), tags$br(),
        sprintf("Common populations: %d \u2014 %s", length(r$common), paste(r$common, collapse = ", "))
      )
    })

    output$mantel_scatter <- plotly::renderPlotly({
      r <- mantel_result_r()
      shiny::req(length(r$x) > 0L)
      x_s <- seq(min(r$x), max(r$x), length.out = 100); y_s <- r$intercept + r$slope * x_s
      plotly::plot_ly() %>%
        plotly::add_markers(x = r$x, y = r$y, marker = list(color = "#7A5DC7", size = 8, opacity = 0.8), name = "Pairs") %>%
        plotly::add_lines(x = x_s, y = y_s, line = list(color = "#B40F20", width = 2),
          name = sprintf("OLS: b=%.4f, R\u00b2=%.4f", r$slope, r$r2)) %>%
        plotly::layout(xaxis = list(title = r$x_label), yaxis = list(title = r$y_label),
          title = list(text = sprintf("%s=%.4f, p=%.4f", r$stat_label, r$stat_obs, r$p_pos), font = list(size = 12)),
          legend = list(x = 0.02, y = 0.98, bgcolor = "rgba(255,255,255,0.8)"), margin = list(t = 40))
    })

    output$mantel_hist <- plotly::renderPlotly({
      r <- mantel_result_r()
      shiny::req(length(r$perm_stats) > 0L)
      plotly::plot_ly() %>%
        plotly::add_histogram(x = r$perm_stats, nbinsx = 60,
          marker = list(color = "rgba(122,93,199,0.55)", line = list(color = "rgba(122,93,199,1)", width = 0.4))) %>%
        plotly::layout(
          shapes = list(list(type = "line", x0 = r$stat_obs, x1 = r$stat_obs, y0 = 0, y1 = 1, yref = "paper",
                              line = list(color = "#B40F20", width = 2, dash = "dash"))),
          xaxis = list(title = r$stat_label), yaxis = list(title = "Count"),
          title = list(text = sprintf("n = %d permutations", length(r$perm_stats)), font = list(size = 11)),
          margin = list(t = 40), showlegend = FALSE)
    })

    output$dt_mantel_data <- DT::renderDT({
      r <- mantel_result_r()
      df <- data.frame(X = round(r$x, 6), Y = round(r$y, 6))
      DT::datatable(df, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 10, dom = "lrtip"),
        class = "compact stripe hover")
    })
    output$dl_mantel_csv <- downloadHandler(
      filename = function() paste0("mantel_data_", Sys.Date(), ".csv"),
      content  = function(file) {
        r <- mantel_result_r()
        write.csv(data.frame(X = r$x, Y = r$y), file, row.names = FALSE)
      }
    )

  })
}