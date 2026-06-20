# server_isolation_by_distance.R
# Isolation by Distance & Mantel Test — SPG-V1 specification
#
# Tab 1: Pairwise dataset (load external OR compute Dgeo from GPS)
# Tab 2: IBD Regression (Rousset 1997) — 1D/2D, FR/FR-ENA, 3 regression lines, CI threshold tuning
# Tab 3: Mantel Test — rectangular matrices, Pearson r OR Rousset slope b, p=(b+1)/(m+1)

server_isolation_by_distance <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    # ── Helpers ────────────────────────────────────────────────────────────────
    `%||%` <- function(a, b) if (!is.null(a)) a else b
    sql_id <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))

    # ── Haversine via geosphere (Hijmans et al. 2019) ──────────────────────────
    haversine_km <- function(lat1, lon1, lat2, lon2) {
      if (requireNamespace("geosphere", quietly = TRUE)) {
        geosphere::distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000
      } else {
        # Fallback
        R <- 6371.0
        dlat <- (lat2 - lat1) * pi / 180
        dlon <- (lon2 - lon1) * pi / 180
        a <- sin(dlat / 2)^2 +
             cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2)^2
        2 * R * asin(sqrt(a))
      }
    }

    geo_dist_matrix <- function(coords) {
      pops <- coords$Population
      n <- nrow(coords)
      mat <- matrix(0.0, n, n, dimnames = list(pops, pops))
      for (i in seq_len(n - 1L)) {
        for (j in (i + 1L):n) {
          d <- haversine_km(coords$Latitude[i], coords$Longitude[i],
                            coords$Latitude[j], coords$Longitude[j])
          mat[i, j] <- mat[j, i] <- d
        }
      }
      mat
    }

    # ── Parse distance file (square or rectangular column-wise) ────────────────
    parse_dist_file <- function(file_path, format = "rectangular") {
      first_lines <- readLines(file_path, n = 5, warn = FALSE)
      sep <- if (any(grepl("\t", first_lines))) "\t" else ","
      raw <- tryCatch(
        read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE, sep = sep),
        error = function(e) read.table(file_path, header = TRUE, sep = sep,
                                       check.names = FALSE, stringsAsFactors = FALSE))

      if (format == "square") {
        rn <- as.character(raw[[1L]])
        mat <- as.matrix(raw[, -1L, drop = FALSE])
        storage.mode(mat) <- "double"
        rownames(mat) <- rn; colnames(mat) <- rn
        return(mat)
      }

      # Rectangular: column headers = pair labels
      pair_names <- names(raw)
      split_pair <- function(s) {
        for (pat in c("-", "_", " vs ", " - ", " _ ")) {
          parts <- strsplit(s, pat, fixed = TRUE)[[1L]]
          if (length(parts) == 2L) return(trimws(parts))
        }
        NULL
      }

      pops <- character(0); pair_list <- list()
      for (pn in pair_names) {
        parts <- split_pair(pn)
        if (!is.null(parts)) { pair_list[[pn]] <- parts; pops <- c(pops, parts) }
      }
      pops <- unique(pops)

      avg_vals <- sapply(pair_names, function(cn) {
        v <- suppressWarnings(as.numeric(raw[[cn]]))
        mean(v, na.rm = TRUE)
      })

      n <- length(pops)
      mat <- matrix(NA_real_, n, n, dimnames = list(pops, pops))
      for (pn in pair_names) {
        parts <- pair_list[[pn]]
        if (!is.null(parts)) {
          p1 <- parts[1L]; p2 <- parts[2L]
          if (p1 %in% pops && p2 %in% pops) mat[p1, p2] <- mat[p2, p1] <- avg_vals[[pn]]
        }
      }
      diag(mat) <- 0
      mat
    }

    # ── Mantel test — p = (b+1)/(m+1) ─────────────────────────────────────────
    # stat = "pearson" or "rousset" (slope of regression)
    mantel_test <- function(mat1, mat2, n_perm = 10000L,
                            stat = "pearson", side = "greater") {
      ok <- !is.na(mat1) & !is.na(mat2) & lower.tri(mat1)
      v1 <- mat1[ok]; v2 <- mat2[ok]
      n <- length(v1)
      if (n < 3L) {
        return(list(stat = NA_real_, p_value = NA_real_, n = n,
                    b_ge = 0L, v1 = v1, v2 = v2, r2 = NA_real_,
                    message = "Not enough data points (need >= 3 pairs)"))
      }

      # Observed statistic
      compute_stat <- function(x, y) {
        if (stat == "pearson") cor(x, y, use = "complete.obs")
        else {
          fit <- suppressWarnings(lm(y ~ x))
          unname(coef(fit)[2L])
        }
      }
      stat_obs <- compute_stat(v1, v2)

      # R² (only for Pearson)
      r2 <- if (stat == "pearson") stat_obs^2 else {
        fit <- suppressWarnings(lm(v2 ~ v1))
        summary(fit)$r.squared
      }

      # Permutations — permute rows/columns of mat2 jointly
      n_mat <- nrow(mat2)
      perm_stats <- numeric(n_perm)
      for (k in seq_len(n_perm)) {
        idx <- sample(n_mat)
        perm_mat2 <- mat2[idx, idx]
        v2_perm <- perm_mat2[ok]
        perm_stats[k] <- compute_stat(v1, v2_perm)
      }
      perm_stats <- perm_stats[is.finite(perm_stats)]

      # b = number of permutations with stat >= stat_obs (one-sided)
      b_ge <- sum(perm_stats >= stat_obs)
      m <- length(perm_stats)

      # p-value: (b+1)/(m+1)
      p_value <- switch(side,
        "greater"   = (b_ge + 1) / (m + 1),
        "less"      = (m - b_ge + 1) / (m + 1),
        "two.sided" = 2 * min((b_ge + 1) / (m + 1), (m - b_ge + 1) / (m + 1))
      )

      list(stat = stat_obs, p_value = p_value, n = n,
           b_ge = b_ge, m = m, r2 = r2,
           v1 = v1, v2 = v2, perm_stats = perm_stats,
           stat_name = stat, side = side,
           message = "OK")
    }

    # ── DB plumbing ────────────────────────────────────────────────────────────
    db_tick    <- reactive({ rv$db_tick })
    con_r      <- reactive({ req(rv$con); rv$con })
    tbl_meta_r <- reactive({ rv$tbl_meta %||% "meta" })

    db_ready <- reactive({
      db_tick(); con <- con_r()
      shiny::req(isTRUE(rv$db_ready))
      shiny::validate(
        shiny::need(DBI::dbExistsTable(con, tbl_meta_r()), "DuckDB meta table missing.")
      )
      TRUE
    })

    meta_schema_r <- reactive({
      db_ready(); con <- con_r()
      info <- DBI::dbGetQuery(con,
        sprintf("PRAGMA table_info(%s)", DBI::dbQuoteIdentifier(con, tbl_meta_r())))
      cols    <- info$name
      ind_col <- if ("individual" %in% cols) "individual"
                 else if ("indiv_id" %in% cols) "indiv_id"
                 else shiny::validate(shiny::need(FALSE, "No individual column found in meta."))
      pop_col <- c("Population","population","pop","pop_code")[
        c("Population","population","pop","pop_code") %in% cols][1]
      shiny::validate(shiny::need(!is.na(pop_col), "No population column found in meta."))
      list(ind_col = ind_col, pop_col = pop_col,
           has_lat = "Latitude" %in% cols, has_lon = "Longitude" %in% cols)
    })

    # ── GPS centroids ──────────────────────────────────────────────────────────
    gps_coords_r <- reactive({
      db_ready(); con <- con_r(); ms <- meta_schema_r()
      shiny::validate(shiny::need(
        ms$has_lat && ms$has_lon,
        "No Latitude/Longitude columns in meta table. Please re-import with GPS data (decimal degrees)."))
      df <- DBI::dbGetQuery(con, sprintf(
        "SELECT CAST(%s AS VARCHAR) AS Population,
                AVG(CAST(Latitude  AS DOUBLE)) AS Latitude,
                AVG(CAST(Longitude AS DOUBLE)) AS Longitude
         FROM %s
         WHERE %s IS NOT NULL
           AND Latitude IS NOT NULL AND Longitude IS NOT NULL
         GROUP BY Population ORDER BY Population",
        sql_id(con, ms$pop_col),
        sql_id(con, tbl_meta_r()),
        sql_id(con, ms$pop_col)))
      df
    })

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 1 — PAIRWISE DATASET
    # ══════════════════════════════════════════════════════════════════════════

    # Raw pairwise data (before filtering)
    raw_pairwise_r <- reactive({
      src <- input$data_source
      if (src == "external") {
        shiny::req(input$file_ext)
        first_lines <- readLines(input$file_ext$datapath, n = 5, warn = FALSE)
        sep <- if (any(grepl("\t", first_lines))) "\t" else ","
        df <- tryCatch(
          read.csv(input$file_ext$datapath, check.names = FALSE,
                   stringsAsFactors = FALSE, sep = sep),
          error = function(e) read.table(input$file_ext$datapath, header = TRUE,
                                         sep = sep, check.names = FALSE,
                                         stringsAsFactors = FALSE))
        shiny::validate(shiny::need(nrow(df) >= 1, "Empty dataset."))
        df
      } else {
        # Compute from GPS
        coords <- gps_coords_r()
        shiny::validate(shiny::need(nrow(coords) >= 2,
          "At least 2 populations with GPS required."))
        geo_mat <- geo_dist_matrix(coords)
        pops <- coords$Population
        n <- length(pops)
        rows <- list()
        for (i in seq_len(n - 1L)) {
          for (j in (i + 1L):n) {
            d <- geo_mat[pops[i], pops[j]]
            rows[[length(rows) + 1L]] <- data.frame(
              Subsample1 = pops[i], Subsample2 = pops[j],
              Dgeo = d, lnDgeo = log(d),
              stringsAsFactors = FALSE)
          }
        }
        do.call(rbind, rows)
      }
    })

    # Update column selectors in Tab 2 and Tab 3
    observe({
      df <- raw_pairwise_r()
      cols <- names(df)
      updateSelectInput(session, "ibd_col_pop1", choices = cols,
                        selected = if ("Subsample1" %in% cols) "Subsample1" else cols[1])
      updateSelectInput(session, "ibd_col_pop2", choices = cols,
                        selected = if ("Subsample2" %in% cols) "Subsample2" else cols[2])
      updateSelectInput(session, "ibd_col_dgeo", choices = cols,
                        selected = if ("lnDgeo" %in% cols) "lnDgeo"
                                   else if ("Dgeo" %in% cols) "Dgeo" else cols[3])
      updateSelectInput(session, "ibd_col_fr", choices = cols,
                        selected = if ("FR" %in% cols) "FR"
                                   else if ("FR_ENA" %in% cols) "FR_ENA" else cols[4])
      updateSelectInput(session, "ibd_col_frl", choices = cols,
                        selected = if ("FR_l" %in% cols) "FR_l" else "")
      updateSelectInput(session, "ibd_col_fru", choices = cols,
                        selected = if ("FR_u" %in% cols) "FR_u" else "")
      updateSelectInput(session, "m1_col", choices = cols)
      updateSelectInput(session, "m2_col", choices = cols)

      # Update filter selectors
      p1_col <- input$ibd_col_pop1 %||% "Subsample1"
      p2_col <- input$ibd_col_pop2 %||% "Subsample2"
      all_pops <- sort(unique(c(df[[p1_col]], df[[p2_col]])))
      updateSelectizeInput(session, "filter_pop1", choices = all_pops, server = TRUE)
      updateSelectizeInput(session, "filter_pop2", choices = all_pops, server = TRUE)
    })

    # Filtered pairwise data
    filtered_pairwise_r <- eventReactive(input$apply_filter, {
      df <- raw_pairwise_r()
      p1_col <- input$ibd_col_pop1 %||% "Subsample1"
      p2_col <- input$ibd_col_pop2 %||% "Subsample2"

      if (!is.null(input$filter_pop1) && length(input$filter_pop1) > 0) {
        df <- df[df[[p1_col]] %in% input$filter_pop1, , drop = FALSE]
      }
      if (!is.null(input$filter_pop2) && length(input$filter_pop2) > 0) {
        df <- df[df[[p2_col]] %in% input$filter_pop2, , drop = FALSE]
      }
      if ("Dgeo" %in% names(df)) {
        if (!is.na(input$dgeo_min)) df <- df[df$Dgeo >= input$dgeo_min, , drop = FALSE]
        if (!is.na(input$dgeo_max)) df <- df[df$Dgeo <= input$dgeo_max, , drop = FALSE]
      }
      df
    }, ignoreInit = FALSE)

    # Initialize filtered = raw on first load
    observe({
      raw_pairwise_r()
      # Trigger apply_filter silently at start
    })

    # Value boxes
    output$vb_pairs <- renderUI({
      tryCatch({
        df <- tryCatch(filtered_pairwise_r(), error = function(e) raw_pairwise_r())
        tags$span(nrow(df))
      }, error = function(e) tags$span("\u2014"))
    })
    output$vb_pops <- renderUI({
      tryCatch({
        df <- tryCatch(filtered_pairwise_r(), error = function(e) raw_pairwise_r())
        p1 <- input$ibd_col_pop1 %||% "Subsample1"
        p2 <- input$ibd_col_pop2 %||% "Subsample2"
        tags$span(length(unique(c(df[[p1]], df[[p2]]))))
      }, error = function(e) tags$span("\u2014"))
    })
    output$vb_cols <- renderUI({
      tryCatch({
        df <- tryCatch(filtered_pairwise_r(), error = function(e) raw_pairwise_r())
        tags$span(ncol(df))
      }, error = function(e) tags$span("\u2014"))
    })

    # Pairwise table
    output$dt_pairwise <- DT::renderDT({
      df <- tryCatch(filtered_pairwise_r(), error = function(e) raw_pairwise_r())
      shiny::validate(shiny::need(nrow(df) > 0, "No data."))
      DT::datatable(df, rownames = FALSE,
        options = list(pageLength = 25, scrollX = TRUE, dom = "lftip"),
        class = "compact hover stripe") |>
        DT::formatRound(columns = sapply(df, is.numeric) |> which() |> names(), digits = 6)
    }, server = TRUE)

    # Downloads
    output$dl_pairwise_csv <- downloadHandler(
      filename = function() paste0("pairwise_dataset_", Sys.Date(), ".csv"),
      content = function(file) {
        df <- tryCatch(filtered_pairwise_r(), error = function(e) raw_pairwise_r())
        write.csv(df, file, row.names = FALSE)
      }
    )
    output$dl_pairwise_txt <- downloadHandler(
      filename = function() paste0("pairwise_dataset_", Sys.Date(), ".txt"),
      content = function(file) {
        df <- tryCatch(filtered_pairwise_r(), error = function(e) raw_pairwise_r())
        write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 2 — IBD REGRESSION (Rousset 1997)
    # ══════════════════════════════════════════════════════════════════════════

    ibd_results_r <- eventReactive(input$run_ibd, {
      df <- tryCatch(filtered_pairwise_r(), error = function(e) raw_pairwise_r())
      shiny::validate(shiny::need(nrow(df) >= 3, "Need at least 3 pairs."))

      p1_col  <- input$ibd_col_pop1
      p2_col  <- input$ibd_col_pop2
      dgeo_col <- input$ibd_col_dgeo
      fr_col  <- input$ibd_col_fr
      frl_col <- input$ibd_col_frl
      fru_col <- input$ibd_col_fru

      shiny::validate(shiny::need(
        all(c(p1_col, p2_col, dgeo_col, fr_col) %in% names(df)),
        "Some selected columns are missing from the dataset."))

      model <- input$ibd_model
      x_raw <- df[[dgeo_col]]
      x <- if (model == "2D") log(x_raw) else x_raw
      y <- df[[fr_col]]

      # Regression on FR
      ok <- is.finite(x) & is.finite(y)
      fit_avg <- lm(y[ok] ~ x[ok])
      b_avg <- unname(coef(fit_avg)[2L])
      int_avg <- unname(coef(fit_avg)[1L])

      # Regression on FR-l and FR-u if available
      fit_l <- fit_u <- NULL
      b_l <- b_u <- int_l <- int_u <- NA_real_
      if (!is.null(frl_col) && nzchar(frl_col) && frl_col %in% names(df)) {
        yl <- df[[frl_col]]
        ok_l <- is.finite(x) & is.finite(yl)
        if (sum(ok_l) >= 3) {
          fit_l <- lm(yl[ok_l] ~ x[ok_l])
          b_l <- unname(coef(fit_l)[2L]); int_l <- unname(coef(fit_l)[1L])
        }
      }
      if (!is.null(fru_col) && nzchar(fru_col) && fru_col %in% names(df)) {
        yu <- df[[fru_col]]
        ok_u <- is.finite(x) & is.finite(yu)
        if (sum(ok_u) >= 3) {
          fit_u <- lm(yu[ok_u] ~ x[ok_u])
          b_u <- unname(coef(fit_u)[2L]); int_u <- unname(coef(fit_u)[1L])
        }
      }

      # Nb = 1/b, Nem = 1/(2*pi*b) — only for 2D model
      Nb  <- if (is.finite(b_avg) && b_avg > 0) 1 / b_avg else NA_real_
      Nem <- if (is.finite(b_avg) && b_avg > 0) 1 / (2 * pi * b_avg) else NA_real_

      # IBD status
      status <- if (is.finite(b_avg) && is.finite(b_l) && b_avg > 0 && b_l > 0) "IBD confirmed"
                else if (is.finite(b_avg) && b_avg > 0) "IBD (no lower CI)"
                else "No IBD"

      list(
        df = df, model = model, x = x, y = y, x_raw = x_raw,
        p1_col = p1_col, p2_col = p2_col,
        fr_col = fr_col, frl_col = frl_col, fru_col = fru_col,
        fit_avg = fit_avg, fit_l = fit_l, fit_u = fit_u,
        b_avg = b_avg, int_avg = int_avg,
        b_l = b_l, int_l = int_l,
        b_u = b_u, int_u = int_u,
        Nb = Nb, Nem = Nem, status = status
      )
    })

    # CI threshold tuning (ORIGINAL FEATURE)
    output$ui_ibd_threshold_result <- renderUI({
      r <- tryCatch(ibd_results_r(), error = function(e) NULL)
      if (is.null(r)) return(tags$p("Run IBD regression first.", style = "color:#94a3b8;"))

      ci_level <- input$ci_threshold
      alpha <- 1 - ci_level

      # The FR-l and FR-u columns already contain the CI bounds at some level.
      # We recompute the slope of the lower line at the requested CI level
      # by linear interpolation between the observed slope and the bounds.
      # This is an approximation — the exact method would require raw bootstrap replicates.
      b_avg <- r$b_avg
      b_l   <- r$b_l
      b_u   <- r$b_u

      if (!is.finite(b_avg) || !is.finite(b_l)) {
        return(tags$div(class = "ibd-warn",
          "Cannot compute threshold: FR-l column is missing or insufficient data."))
      }

      # Find the CI level where the lower slope crosses zero
      # Linear interpolation: slope_l(ci) = b_avg - k * (1-ci)
      # At ci0 (the CI used in data), slope_l = b_l
      # Assume slope_l(ci) = b_avg - (b_avg - b_l) * (1-ci) / (1-ci0)
      # For ci0 = 0.95 (default):
      ci0 <- 0.95
      k <- (b_avg - b_l) / (1 - ci0)
      # slope_l(ci) = b_avg - k * (1-ci) = 0  =>  ci* = 1 - b_avg/k
      if (k > 0) {
        ci_star <- 1 - b_avg / k
        if (is.finite(ci_star) && ci_star > 0 && ci_star < 1) {
          p_approx <- 1 - ci_star
          tags$div(class = "ibd-result",
            tags$strong(sprintf("At CI = %.1f%% (current slider):", ci_level * 100)),
            tags$br(),
            sprintf("Approx. lower slope at this CI: %.6f",
                    b_avg - k * (1 - ci_level)),
            tags$br(), tags$br(),
            tags$strong(sprintf("Exact threshold CI: %.3f%%", ci_star * 100)),
            tags$br(),
            sprintf("Approx. p-value of IBD slope: %.4f", p_approx),
            tags$br(),
            if (ci_level >= ci_star)
              tags$span(style = "color:#166534;font-weight:600;",
                        "\u2714 At this CI level, the lower slope is positive \u2192 IBD significant")
            else
              tags$span(style = "color:#9d174d;font-weight:600;",
                        "\u2718 At this CI level, the lower slope is negative \u2192 IBD not significant")
          )
        } else {
          tags$div(class = "ibd-warn", "Threshold CI out of range [0,1].")
        }
      } else {
        tags$div(class = "ibd-warn", "Cannot compute threshold (k \u2264 0).")
      }
    })

    # Value boxes
    output$vb_ibd_status <- renderUI({
      tryCatch({
        r <- ibd_results_r()
        col <- switch(r$status,
          "IBD confirmed"      = "#166534",
          "IBD (no lower CI)"  = "#854d0e",
          "#9d174d")
        tags$span(style = paste0("color:", col, ";font-size:13px;"), r$status)
      }, error = function(e) tags$span("\u2014"))
    })
    output$vb_ibd_b <- renderUI({
      tryCatch({
        r <- ibd_results_r()
        tags$span(formatC(r$b_avg, format = "e", digits = 3))
      }, error = function(e) tags$span("\u2014"))
    })
    output$vb_ibd_nb <- renderUI({
      tryCatch({
        r <- ibd_results_r()
        if (is.na(r$Nb)) tags$span("\u2014") else tags$span(round(r$Nb, 2))
      }, error = function(e) tags$span("\u2014"))
    })
    output$vb_ibd_nem <- renderUI({
      tryCatch({
        r <- ibd_results_r()
        if (is.na(r$Nem)) tags$span("\u2014") else tags$span(round(r$Nem, 3))
      }, error = function(e) tags$span("\u2014"))
    })

    # Regression table
    output$dt_ibd_reg <- DT::renderDT({
      r <- ibd_results_r()
      df <- data.frame(
        Line = c("Average (FR)", "Lower CI (FR-l)", "Upper CI (FR-u)"),
        Intercept = c(r$int_avg, r$int_l, r$int_u),
        Slope_b   = c(r$b_avg, r$b_l, r$b_u),
        stringsAsFactors = FALSE
      )
      if (r$model == "2D") {
        df$Nb  <- c(r$Nb, if (is.finite(r$b_l) && r$b_l > 0) 1/r$b_l else NA,
                          if (is.finite(r$b_u) && r$b_u > 0) 1/r$b_u else NA)
        df$Nem <- c(r$Nem, if (is.finite(r$b_l) && r$b_l > 0) 1/(2*pi*r$b_l) else NA,
                           if (is.finite(r$b_u) && r$b_u > 0) 1/(2*pi*r$b_u) else NA)
      }
      DT::datatable(df, rownames = FALSE,
        options = list(dom = "t", pageLength = 3, ordering = FALSE),
        class = "compact stripe") |>
        DT::formatStyle("Line", target = "row",
          backgroundColor = DT::styleEqual(
            c("Average (FR)", "Lower CI (FR-l)", "Upper CI (FR-u)"),
            c("#f5f5f5", "#fff0f0", "#f0f8ff")))
    })

    # Interpretation
    output$ui_ibd_interpretation <- renderUI({
      r <- ibd_results_r()
      model_txt <- if (r$model == "2D")
        "FR ~ ln(Dgeo) (Rousset 1997, 2D model for migration-mutation-drift equilibrium)"
      else
        "FR ~ Dgeo (1D model)"
      tags$div(class = "ibd-result",
        tags$strong("Model: "), model_txt,
        tags$br(),
        tags$strong("Genetic distance: "), r$fr_col,
        tags$br(),
        tags$strong("Slope b (average): "), formatC(r$b_avg, format = "e", digits = 4),
        tags$br(),
        if (!is.na(r$Nb))
          HTML(paste0(tags$strong("Nb = 1/b: "), round(r$Nb, 2),
                      " (neighbourhood size)"))
        else "",
        tags$br(),
        if (!is.na(r$Nem))
          HTML(paste0(tags$strong("Nem = 1/(2\u03c0b): "), round(r$Nem, 3),
                      " (effective migrants per generation)"))
        else "",
        tags$br(), tags$br(),
        if (r$status == "IBD confirmed")
          tags$span(style = "color:#166534;font-weight:600;",
                    "\u2714 The slope and its CI are all positive \u2192 population under isolation by distance.")
        else if (r$status == "IBD (no lower CI)")
          tags$span(style = "color:#854d0e;font-weight:600;",
                    "\u26a0 Average slope is positive but lower CI is missing \u2014 interpret with caution.")
        else
          tags$span(style = "color:#9d174d;font-weight:600;",
                    "\u2718 No evidence for isolation by distance (slope \u2264 0). ",
                    "Consider a Mantel test with DCSE for more power (Séré et al. 2017).")
      )
    })

    # IBD plot
    output$ibd_plot <- plotly::renderPlotly({
      r <- ibd_results_r()
      x <- r$x; y <- r$y
      x_label <- if (r$model == "2D") "ln(geographic distance)" else "Geographic distance"

      x_seq <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 100)
      line_avg <- data.frame(x = x_seq, y = r$int_avg + r$b_avg * x_seq)

      p <- plotly::plot_ly() |>
        plotly::add_markers(
          x = x, y = y,
          marker = list(color = "#2CBF9F", size = 7, opacity = 0.85),
          name = "Pairs",
          hoverinfo = "text",
          text = ~paste0(r$df[[r$p1_col]], " \u2013 ", r$df[[r$p2_col]],
                        "<br>", r$fr_col, ": ", round(y, 5),
                        "<br>", r$ibd_col_dgeo %||% "Dgeo", ": ", round(r$x_raw, 2))
        ) |>
        plotly::add_lines(
          data = line_avg, x = ~x, y = ~y,
          line = list(color = "#333a43", width = 2.5),
          name = sprintf("Average  b = %.4e", r$b_avg)
        )

      if (is.finite(r$b_l)) {
        line_l <- data.frame(x = x_seq, y = r$int_l + r$b_l * x_seq)
        p <- p |> plotly::add_lines(
          data = line_l, x = ~x, y = ~y,
          line = list(color = "#B40F20", width = 1.5, dash = "dash"),
          name = sprintf("Lower CI  b = %.4e", r$b_l))
      }
      if (is.finite(r$b_u)) {
        line_u <- data.frame(x = x_seq, y = r$int_u + r$b_u * x_seq)
        p <- p |> plotly::add_lines(
          data = line_u, x = ~x, y = ~y,
          line = list(color = "#3B9AB2", width = 1.5, dash = "dot"),
          name = sprintf("Upper CI  b = %.4e", r$b_u))
      }

      p |> plotly::layout(
        title = list(text = sprintf("IBD regression \u2014 %s \u2014 %s",
                                    if (r$model == "2D") "Model 2D" else "Model 1D",
                                    r$fr_col),
                     font = list(size = 14)),
        xaxis = list(title = x_label),
        yaxis = list(title = r$fr_col),
        legend = list(x = 0.02, y = 0.98, bgcolor = "rgba(255,255,255,0.8)")
      )
    })

    # Download IBD regression
    output$dl_ibd_csv <- downloadHandler(
      filename = function() paste0("IBD_regression_", Sys.Date(), ".csv"),
      content = function(file) {
        r <- ibd_results_r()
        df <- data.frame(
          Line = c("Average (FR)", "Lower CI (FR-l)", "Upper CI (FR-u)"),
          Intercept = c(r$int_avg, r$int_l, r$int_u),
          Slope_b   = c(r$b_avg, r$b_l, r$b_u),
          Nb  = c(r$Nb, if (is.finite(r$b_l) && r$b_l > 0) 1/r$b_l else NA,
                        if (is.finite(r$b_u) && r$b_u > 0) 1/r$b_u else NA),
          Nem = c(r$Nem, if (is.finite(r$b_l) && r$b_l > 0) 1/(2*pi*r$b_l) else NA,
                        if (is.finite(r$b_u) && r$b_u > 0) 1/(2*pi*r$b_u) else NA),
          Model = r$model,
          Genetic_distance = r$fr_col,
          stringsAsFactors = FALSE)
        write.csv(df, file, row.names = FALSE)
      }
    )

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 3 — MANTEL TEST
    # ══════════════════════════════════════════════════════════════════════════

    # Build matrix from a column of pairwise dataset
    build_matrix_from_col <- function(df, p1_col, p2_col, val_col) {
      pops <- sort(unique(c(df[[p1_col]], df[[p2_col]])))
      n <- length(pops)
      mat <- matrix(NA_real_, n, n, dimnames = list(pops, pops))
      diag(mat) <- 0
      for (i in seq_len(nrow(df))) {
        p1 <- df[[p1_col]][i]; p2 <- df[[p2_col]][i]
        v <- suppressWarnings(as.numeric(df[[val_col]][i]))
        if (p1 %in% pops && p2 %in% pops && is.finite(v)) {
          mat[p1, p2] <- mat[p2, p1] <- v
        }
      }
      mat
    }

    # Matrix 1
    mat1_r <- reactive({
      if (input$m1_source == "col") {
        shiny::req(input$m1_col)
        df <- tryCatch(filtered_pairwise_r(), error = function(e) raw_pairwise_r())
        p1_col <- input$ibd_col_pop1 %||% "Subsample1"
        p2_col <- input$ibd_col_pop2 %||% "Subsample2"
        build_matrix_from_col(df, p1_col, p2_col, input$m1_col)
      } else {
        shiny::req(input$m1_file)
        parse_dist_file(input$m1_file$datapath, input$m1_format)
      }
    })

    # Matrix 2
    mat2_r <- reactive({
      src <- input$m2_source
      if (src == "col") {
        shiny::req(input$m2_col)
        df <- tryCatch(filtered_pairwise_r(), error = function(e) raw_pairwise_r())
        p1_col <- input$ibd_col_pop1 %||% "Subsample1"
        p2_col <- input$ibd_col_pop2 %||% "Subsample2"
        build_matrix_from_col(df, p1_col, p2_col, input$m2_col)
      } else if (src %in% c("gps_km", "gps_ln")) {
        coords <- gps_coords_r()
        mat <- geo_dist_matrix(coords)
        if (src == "gps_ln") {
          mat[mat > 0] <- log(mat[mat > 0])
          diag(mat) <- 0
        }
        mat
      } else {
        shiny::req(input$m2_file)
        parse_dist_file(input$m2_file$datapath, input$m2_format)
      }
    })

    # Run Mantel
    mantel_results_r <- eventReactive(input$run_mantel, {
      m1 <- mat1_r(); m2 <- mat2_r()
      common <- intersect(rownames(m1), rownames(m2))
      shiny::validate(shiny::need(
        length(common) >= 3L,
        sprintf("Need at least 3 populations in BOTH matrices. Found: %d.",
                length(common))))
      m1 <- m1[common, common, drop = FALSE]
      m2 <- m2[common, common, drop = FALSE]

      n_perm <- as.integer(input$n_perm %||% 10000L)
      stat   <- input$mantel_stat %||% "pearson"
      side   <- input$mantel_side %||% "greater"

      withProgress(message = "Running Mantel test...", value = 0.5, {
        res <- mantel_test(m1, m2, n_perm = n_perm, stat = stat, side = side)
        incProgress(0.5, detail = "Done")
      })

      res$mat1 <- m1; res$mat2 <- m2
      res$n_perm <- n_perm; res$pops <- common
      res
    })

    # Value boxes
    output$vb_mantel_stat <- renderUI({
      tryCatch({
        r <- mantel_results_r()
        tags$span(formatC(r$stat, digits = 4, format = "f"))
      }, error = function(e) tags$span("\u2014"))
    })
    output$vb_mantel_p <- renderUI({
      tryCatch({
        r <- mantel_results_r()
        sig <- if (r$p_value < 0.001) "***" else if (r$p_value < 0.01) "**"
               else if (r$p_value < 0.05) "*" else "ns"
        col <- if (r$p_value < 0.05) "#166534" else "#9d174d"
        tags$span(style = paste0("color:", col, ";"),
                  sprintf("%.4f %s", r$p_value, sig))
      }, error = function(e) tags$span("\u2014"))
    })
    output$vb_mantel_n <- renderUI({
      tryCatch(tags$span(mantel_results_r()$n), error = function(e) tags$span("\u2014"))
    })
    output$vb_mantel_r2 <- renderUI({
      tryCatch({
        r <- mantel_results_r()
        if (is.na(r$r2)) tags$span("\u2014")
        else tags$span(sprintf("%.2f%%", r$r2 * 100))
      }, error = function(e) tags$span("\u2014"))
    })
    output$vb_mantel_b <- renderUI({
      tryCatch({
        r <- mantel_results_r()
        tags$span(sprintf("%d / %d", r$b_ge, r$m))
      }, error = function(e) tags$span("\u2014"))
    })
    output$vb_mantel_pops <- renderUI({
      tryCatch(tags$span(length(mantel_results_r()$pops)),
               error = function(e) tags$span("\u2014"))
    })

    # Mantel result text
    output$ui_mantel_result <- renderUI({
      r <- tryCatch(mantel_results_r(), error = function(e) NULL)
      if (is.null(r))
        return(tags$p("Configure matrices and click 'Run Mantel Test'.",
                      style = "color:#94a3b8;"))
      stat_label <- if (r$stat_name == "pearson") "Pearson r" else "Rousset slope b"
      side_label <- switch(r$side,
        "greater"   = "positive correlation (one-sided)",
        "less"      = "negative correlation (one-sided)",
        "two.sided" = "two-sided")
      tags$div(class = "ibd-result",
        tags$strong("Mantel Test Summary"),
        tags$br(),
        sprintf("Statistic:          %s = %.6f", stat_label, r$stat),
        tags$br(),
        sprintf("Permutations (m):   %d", r$n_perm),
        tags$br(),
        sprintf("b (\u2265 observed):    %d", r$b_ge),
        tags$br(),
        sprintf("p-value formula:    (b+1)/(m+1) = (%d+1)/(%d+1)", r$b_ge, r$m),
        tags$br(),
        sprintf("p-value:            %.6f", r$p_value),
        tags$br(),
        sprintf("R\u00b2 (variance):     %.4f%%", r$r2 * 100),
        tags$br(),
        sprintf("Alternative:        %s", side_label),
        tags$br(), tags$br(),
        if (r$p_value < 0.05)
          tags$span(style = "color:#166534;font-weight:600;",
                    sprintf("\u2714 Significant association between the two distance matrices (p < 0.05)."))
        else
          tags$span(style = "color:#9d174d;font-weight:600;",
                    "\u2718 No significant association (p \u2265 0.05).")
      )
    })

    # Download Mantel results
    output$dl_mantel_csv <- downloadHandler(
      filename = function() paste0("mantel_test_", Sys.Date(), ".csv"),
      content = function(file) {
        r <- mantel_results_r()
        hdr <- c(
          "# Mantel Test Results (SPG-V1)",
          paste0("# Statistic: ", if (r$stat_name == "pearson") "Pearson r" else "Rousset slope b"),
          paste0("# Permutations (m): ", r$n_perm),
          paste0("# b (>= observed): ", r$b_ge),
          paste0("# p-value formula: (b+1)/(m+1)"),
          paste0("# p-value: ", r$p_value),
          paste0("# R²: ", r$r2),
          paste0("# Pairs (n): ", r$n),
          paste0("# Populations: ", length(r$pops)),
          "#"
        )
        writeLines(hdr, con = file)
        df <- data.frame(Distance1 = r$v1, Distance2 = r$v2)
        write.table(df, file = file, sep = ",", row.names = FALSE,
                    quote = FALSE, append = TRUE, col.names = TRUE)
      }
    )

    # Mantel scatter plot
    output$mantel_plot <- plotly::renderPlotly({
      r <- tryCatch(mantel_results_r(), error = function(e) NULL)
      if (is.null(r) || length(r$v1) < 3L) {
        return(plotly::plot_ly() |>
          plotly::layout(title = "Run Mantel test to see the scatter plot"))
      }
      df <- data.frame(x = r$v1, y = r$v2)
      fit <- lm(y ~ x, data = df)
      stat_label <- if (r$stat_name == "pearson") "Pearson r" else "Rousset slope b"

      plotly::plot_ly() |>
        plotly::add_markers(
          data = df, x = ~x, y = ~y,
          marker = list(color = "#2CBF9F", size = 7, opacity = 0.85),
          name = "Pairs", hoverinfo = "text",
          text = ~paste0("Matrix 1: ", round(x, 4),
                        "<br>Matrix 2: ", round(y, 4))) |>
        plotly::add_lines(
          data = data.frame(x = df$x, y = fitted(fit)),
          x = ~x, y = ~y,
          line = list(color = "#B40F20", width = 2),
          name = sprintf("Regression (R\u00b2 = %.3f)", r$r2)) |>
        plotly::layout(
          title = list(
            text = sprintf("Mantel test \u2014 %s = %.4f, p = %.4f",
                           stat_label, r$stat, r$p_value),
            font = list(size = 14)),
          xaxis = list(title = "Matrix 1 (X)"),
          yaxis = list(title = "Matrix 2 (Y)"),
          legend = list(x = 0.02, y = 0.98, bgcolor = "rgba(255,255,255,0.8)")
        )
    })

  }) # end moduleServer
}