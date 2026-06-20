# server_null_alleles.R
# Calls engine_freena.R (faithful translation of FreeNA_optm2R.pas) for all
# statistics. This file only handles: DB plumbing (genotype retrieval),
# wiring engine outputs to the UI, and downloads.
#
# IMPORTANT: source engine_freena.R before this file (e.g. in app.R / global.R):
#   source("engine_freena.R")
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
    # (Same reconstruction logic as the other modules in this app.)
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

  })
}