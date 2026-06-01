# server_LD.R 


server_LD <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    
    `%||%` <- function(a, b) if (!is.null(a)) a else b
    
    db_tick <- reactive({ rv$db_tick })
    con_r   <- reactive({ req(rv$con); rv$con })
    
    tbl_meta_r <- reactive({ rv$tbl_meta %||% "meta" })
    
    tbl_hf_r <- reactive({
      con <- con_r()
      if (exists("duck_tbl_exists", mode = "function", inherits = TRUE) &&
          exists(".duckdb_get_params", mode = "function", inherits = TRUE) &&
          duck_tbl_exists(con, "params")) {
        p <- .duckdb_get_params(con)
        th <- p$tbl_hf %||% "hf"
        return(as.character(th))
      }
      "hf"
    })
    
    db_ready <- reactive({
      db_tick()
      con <- con_r()
      shiny::req(isTRUE(rv$db_ready))
      
      shiny::validate(
        shiny::need(DBI::dbExistsTable(con, tbl_meta_r()), "DuckDB meta table missing."),
        shiny::need(DBI::dbExistsTable(con, tbl_hf_r()),   "DuckDB hf table missing.")
      )
      TRUE
    })
    
    base_r <- reactive({
      db_ready()
      
      b <- rv$base_ld %||% rv$base %||% rv$base_r %||% rv$genotype_base
      b <- suppressWarnings(as.integer(b))
      if (length(b) == 1L && is.finite(b) && b > 1L) return(as.integer(b))
      
      con <- con_r()
      if (DBI::dbExistsTable(con, "params") &&
          exists(".duckdb_get_params", mode = "function", inherits = TRUE)) {
        p <- .duckdb_get_params(con)
        b <- suppressWarnings(as.integer(p$base %||% p$base_scalar_full %||% p$base_scalar_preview))
        if (length(b) == 1L && is.finite(b) && b > 1L) return(as.integer(b))
      }
      
      1000L
    })
    
    ld_timing <- reactiveVal(NULL)
    
    ld_data <- reactive({
      db_ready()
      
      con      <- con_r()
      hf_tbl   <- tbl_hf_r()
      meta_tbl <- tbl_meta_r()
      base     <- base_r()
      
      shiny::validate(shiny::need(!is.null(con),      "LD: DuckDB connection is NULL."))
      shiny::validate(shiny::need(!is.null(hf_tbl),   "LD: hf DuckDB table name is NULL."))
      shiny::validate(shiny::need(!is.null(meta_tbl), "LD: meta DuckDB table name is NULL."))
      
      hf_info <- DBI::dbGetQuery(
        con,
        sprintf("PRAGMA table_info(%s)", DBI::dbQuoteIdentifier(con, hf_tbl))
      )
      meta_info <- DBI::dbGetQuery(
        con,
        sprintf("PRAGMA table_info(%s)", DBI::dbQuoteIdentifier(con, meta_tbl))
      )
      
      hf_cols <- hf_info$name
      meta_cols <- meta_info$name
      
      if (all(c("individual", "locus", "g") %in% hf_cols)) {
        hf_ind_col   <- "individual"
        hf_locus_col <- "locus"
        hf_gt_col    <- "g"
      } else if (all(c("indiv_id", "locus_id", "gt") %in% hf_cols)) {
        hf_ind_col   <- "indiv_id"
        hf_locus_col <- "locus_id"
        hf_gt_col    <- "gt"
      } else {
        shiny::validate(
          shiny::need(
            FALSE,
            "LD: hf must contain either (individual,locus,g) or (indiv_id,locus_id,gt)."
          )
        )
      }
      
      if ("individual" %in% meta_cols) {
        meta_ind_col <- "individual"
      } else if ("indiv_id" %in% meta_cols) {
        meta_ind_col <- "indiv_id"
      } else {
        shiny::validate(shiny::need(FALSE, "LD: no individual column found in meta table."))
        }
      
      pop_candidates <- c("Population", "population", "pop", "pop_code")
      pop_col <- pop_candidates[pop_candidates %in% meta_cols][1]
      shiny::validate(shiny::need(!is.na(pop_col), "LD: no population column found in meta table."))
      
      hf_tbl_q   <- as.character(DBI::dbQuoteIdentifier(con, hf_tbl))
      meta_tbl_q <- as.character(DBI::dbQuoteIdentifier(con, meta_tbl))
      hf_ind_q   <- as.character(DBI::dbQuoteIdentifier(con, hf_ind_col))
      hf_locus_q <- as.character(DBI::dbQuoteIdentifier(con, hf_locus_col))
      hf_gt_q    <- as.character(DBI::dbQuoteIdentifier(con, hf_gt_col))
      meta_ind_q <- as.character(DBI::dbQuoteIdentifier(con, meta_ind_col))
      pop_q      <- as.character(DBI::dbQuoteIdentifier(con, pop_col))
      
      loci_sql <- sprintf("
        SELECT DISTINCT CAST(%s AS VARCHAR) AS locus
        FROM %s
        WHERE %s IS NOT NULL
        ORDER BY locus
      ", hf_locus_q, hf_tbl_q, hf_locus_q)
      
      loci_df <- DBI::dbGetQuery(con, loci_sql)
      loci <- as.character(loci_df$locus)
      shiny::validate(shiny::need(length(loci) > 1, "LD: need at least 2 loci."))
      
      base <- as.integer(base)

      case_exprs <- vapply(loci, function(loc) {
        loc_q   <- as.character(DBI::dbQuoteString(con, loc))
        alias_q <- as.character(DBI::dbQuoteIdentifier(con, loc))
        sprintf("MAX(CASE WHEN locus = %s THEN %s END) AS %s",
                loc_q, hf_gt_q, alias_q)
      }, character(1))

      wide_sql <- sprintf("
        WITH long AS (
          SELECT
            CAST(h.%s AS VARCHAR) AS individual,
            CAST(h.%s AS VARCHAR) AS locus,
            CAST(m.%s AS VARCHAR) AS Population,
            h.%s                  AS gt
          FROM %s h
          LEFT JOIN %s m
            ON CAST(h.%s AS VARCHAR) = CAST(m.%s AS VARCHAR)
          WHERE h.%s IS NOT NULL AND h.%s > 0
        )
        SELECT
          Population,
          individual,
          %s
        FROM long
        GROUP BY Population, individual
        ORDER BY Population, individual
      ",
                          hf_ind_q,
                          hf_locus_q,
                          pop_q,
                          hf_gt_q,
                          hf_tbl_q,
                          meta_tbl_q,
                          hf_ind_q,
                          meta_ind_q,
                          hf_gt_q, hf_gt_q,
                          paste(case_exprs, collapse = ",\n          ")
      )
      
      out <- DBI::dbGetQuery(con, wide_sql)
      
      shiny::validate(shiny::need(nrow(out) > 1, "LD: not enough individuals after reshaping."))
      shiny::validate(shiny::need(ncol(out) > 3, "LD: need at least 2 loci after reshaping."))
      
      fixed_cols <- c("Population", "individual")
      locus_cols <- setdiff(names(out), fixed_cols)
      out <- out[, c(fixed_cols, sort(locus_cols)), drop = FALSE]
      
      out
    })
    
    loci_names <- reactive({
      df <- ld_data()
      loci <- setdiff(names(df), c("Population", "individual", "Individual"))
      shiny::validate(shiny::need(length(loci) > 1, "Need at least 2 loci to compute LD."))
      loci
    })
    
    # Warn when fewer than 1000 permutations are requested
    shiny::observeEvent(input$run_LD, {
      np <- suppressWarnings(as.integer(input$n_iterations))
      if (!is.na(np) && np < 1000L) {
        shinyalert::shinyalert(
          title = "Low number of permutations",
          text  = paste0(
            "You have set B = ", np, " permutations.\n\n",
            "The Monte Carlo p-value formula p = (n\u2265obs + 1) / B gives ",
            "a slight overestimation when B is small, which may mislead ",
            "naive users and produce unreliable significance calls.\n\n",
            "A minimum of 1000 permutations is strongly recommended; ",
            "10\u202f000 or more for publication-quality results."
          ),
          type = "warning",
          confirmButtonText = "Continue anyway",
          showCancelButton  = FALSE
        )
      }
    }, ignoreInit = TRUE)

    ld_results_reactive <- eventReactive(input$run_LD, {
      df <- ld_data()
      loci <- loci_names()

      shiny::validate(shiny::need(nrow(df) > 1, "Not enough individuals to compute LD."))
      
      start_time <- Sys.time()
      shinyWidgets::updateProgressBar(session, "LD_progress", value = 10)
      
      geno_mat <- as.matrix(df[, loci, drop = FALSE])
      storage.mode(geno_mat) <- "integer"

      nbperms <- as.integer(input$n_iterations)
      if (is.na(nbperms) || nbperms < 1L) nbperms <- 10000L

      set.seed(1)
      res_cpp <- tryCatch({
        ld_pvalues_cpp(Population = df$Population,
                       geno_mat   = geno_mat,
                       base       = base_r(),
                       nbperms    = nbperms)
      }, error = function(e) {
        showNotification(paste("Error in LD C++ computation:", e$message), type = "error")
        NULL
      })
      if (is.null(res_cpp)) return(NULL)
      
      pv <- data.frame(
        Pair = rownames(res_cpp),
        as.data.frame(res_cpp, check.names = FALSE),
        row.names = NULL,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
      shinyWidgets::updateProgressBar(session, "LD_progress", value = 100)
      ld_timing(round(difftime(Sys.time(), start_time, units = "secs"), 1))
      
      pv
    }, ignoreInit = TRUE)
    
    # ---- Table formatting / filtering ----
    summary_table_reactive <- reactive({
      req(ld_results_reactive())
      pv <- ld_results_reactive()
      if (is.null(pv) || nrow(pv) == 0) return(data.frame(Message = "No data available"))

      dec_places <- input$decimal_places %||% 4
      pv_formatted <- pv

      for (j in 2:ncol(pv_formatted)) {
        nums <- pv_formatted[[j]]
        pv_formatted[[j]] <- ifelse(is.na(nums), "NA",
                                    sprintf(paste0("%.", dec_places, "f"), nums))
      }

      # compute once, reuse for both filtering and sorting
      min_pvals <- min_pvals_by_pair(pv)

      if (!is.null(input$table_view) && input$table_view != "all") {
        threshold <- switch(input$table_view,
                            "sig_05"  = 0.05,
                            "sig_01"  = 0.01,
                            "sig_001" = 0.001,
                            1.0)
        keep <- min_pvals < threshold
        pv_formatted <- pv_formatted[keep, , drop = FALSE]
        min_pvals    <- min_pvals[keep]
      }

      if (!is.null(input$sort_by_ld) && input$sort_by_ld %in% c("pval_asc", "pval_desc")) {
        ord <- order(min_pvals, decreasing = (input$sort_by_ld == "pval_desc"))
        pv_formatted <- pv_formatted[ord, , drop = FALSE]
      }

      pv_formatted
    })
    
    # ---- Value boxes ----
    output$total_pairs_box <- renderValueBox({
      pv <- ld_results_reactive()
      total <- if (is.null(pv)) 0 else nrow(pv)
      
      valueBox(
        value = total,
        subtitle = HTML("<small>Total locus pairs<br>tested for LD</small>"),
        color = "blue",
        icon = icon("link")
      )
    })
    
    output$significant_pairs_box <- renderValueBox({
      pv <- ld_results_reactive()
      if (is.null(pv) || nrow(pv) == 0) {
        return(valueBox("0 (0%)", HTML("<small>Significant pairs<br>No data</small>"), color = "green", icon = icon("exclamation-triangle")))
      }
      
      min_pvals <- min_pvals_by_pair(pv)
      alpha <- input$alpha_level %||% 0.05
      sig_count <- sum(min_pvals < alpha, na.rm = TRUE)
      pct <- round(100 * sig_count / length(min_pvals), 1)
      
      valueBox(
        value = paste0(sig_count, " (", pct, "%)"),
        subtitle = HTML(paste0("<small>Significant pairs<br>p < ", alpha, "</small>")),
        color = if (sig_count > 0) "yellow" else "green",
        icon = icon("exclamation-triangle")
      )
    })
    
    output$mean_pvalue_box <- renderValueBox({
      pv <- ld_results_reactive()
      if (is.null(pv) || nrow(pv) == 0) {
        return(valueBox("N/A", HTML("<small>Mean p-value<br>No data</small>"), color = "red", icon = icon("calculator")))
      }
      
      all_pvals <- suppressWarnings(as.numeric(unlist(pv[, -1, drop = FALSE])))
      all_pvals <- all_pvals[!is.na(all_pvals) & is.finite(all_pvals)]
      mean_p <- if (length(all_pvals) == 0) NA_real_ else mean(all_pvals)
      
      p_display <- ifelse(is.na(mean_p), "N/A",
                          ifelse(mean_p < 0.001, "< 0.001", format(round(mean_p, 4), nsmall = 4)))
      
      valueBox(
        value = p_display,
        subtitle = HTML("<small>Mean p-value<br>across all tests</small>"),
        color = if (is.na(mean_p)) "red" else if (mean_p < 0.05) "red" else "aqua",
        icon = icon("calculator")
      )
    })
    
    output$analysis_time_ld_box <- renderValueBox({
      req(input$run_LD > 0, !is.null(ld_timing()))
      time_sec <- ld_timing()
      time_display <- if (time_sec < 60) paste0(time_sec, " s") else paste0(round(time_sec / 60, 1), " min")
      valueBox(time_display, HTML("<small>Analysis time<br>LD computation</small>"), color = "aqua", icon = icon("clock"))
    })
    
    
    # -----------------------------#
    # helpers (put inside moduleServer)
    # -----------------------------#
    min_pvals_by_pair <- function(pv) {
      if (is.null(pv) || nrow(pv) == 0) return(numeric(0))
      mat <- as.matrix(pv[, -1, drop = FALSE])
      storage.mode(mat) <- "double"
      v <- apply(mat, 1, min, na.rm = TRUE)
      v[!is.finite(v)] <- 1.0
      v
    }
    
    safe_empty_plot <- function(msg) {
      plot.new()
      text(0.5, 0.5, msg, cex = 1.2)
    }
    
    # -----------------------------#
    # summary table
    # -----------------------------#
    output$summary_output <- DT::renderDT({
      df <- summary_table_reactive()
      
      if (is.null(df) || nrow(df) == 0) {
        return(DT::datatable(data.frame(Message = "No data available. Click 'Run LD Analysis' to start.")))
      }
      
      dt <- DT::datatable(
        df,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = c("copy"),
          pageLength = 25,
          scrollX = TRUE,
          scrollY = "400px"
        ),
        rownames = FALSE
      )
      
      # keep highlighting only if you still have that checkbox in UI
      if (isTRUE(input$highlight_sig) && ncol(df) > 1) {
        for (j in 2:ncol(df)) {
          dt <- dt %>%
            DT::formatStyle(
              colnames(df)[j],
              backgroundColor = DT::styleInterval(
                c(0.001, 0.01, 0.05),
                c("#ffcdd2", "#ffecb3", "#fff9c4", "#ffffff")
              )
            )
        }
      }
      
      dt
    })
    
    # -----------------------------#
    # downloads
    # -----------------------------#
    output$download_LD_csv <- downloadHandler(
      filename = function() paste0("LD_results_", Sys.Date(), ".csv"),
      content  = function(file) {
        df <- summary_table_reactive()
        if (!is.null(df) && nrow(df) > 0) write.csv(df, file, row.names = FALSE, quote = FALSE)
      }
    )
    
    output$download_LD_txt <- downloadHandler(
      filename = function() paste0("LD_results_", Sys.Date(), ".txt"),
      content  = function(file) {
        df <- summary_table_reactive()
        if (!is.null(df) && nrow(df) > 0) write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
  })
}