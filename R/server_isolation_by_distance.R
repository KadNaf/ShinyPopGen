# server_isolation_by_distance.R
# Isolation by Distance (Rousset 1997) + generic Mantel test.
#
# This module is a CONTINUATION of the "Null alleles" module: it no longer
# recomputes the null-allele EM / FST / FST-ENA / DCSE / DCSE-INA statistics
# (that used to be tabs 1-4 here, duplicating the Null Alleles module) — it
# reuses the results already computed there, shared through `rv`:
#
#   rv$null_alleles_results   (set by server_null_alleles.R after each
#                              "Compute + Bootstrap + Export" run)
#
# So: go to the Null Alleles module, click Compute, THEN come here.
#
# Geographic distance (D_geo) is the Vincenty ellipsoidal geodesic distance
# (WGS84), in metres — NOT a simple spherical Haversine approximation and
# NOT a planar UTM distance; both were checked against a FreeNA-style
# reference and only Vincenty reproduced it to the metre.
#
# References:
#   Rousset (1997)  — Isolation by distance regression: FR = FST/(1-FST)
#                      regressed on ln(geographic distance) (2D habitat) or
#                      on raw distance (1D habitat); Nb = 1/slope,
#                      Nem = Nb/(2*pi).
#   Mantel (1967) / RT (Manly) / Fstat 2.9.4 convention — permutation test
#      by joint row/column relabelling of one distance matrix; one-sided
#      p = (b+1)/(m+1), b = number of permuted statistics >= observed.

server_isolation_by_distance <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    `%||%` <- function(a, b) if (!is.null(a) && !(length(a) == 1 && is.na(a))) a else b

    # ── DB plumbing (same conventions as other modules) ──────────────────────
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

    # ── Null Alleles module results (shared via rv — nothing recomputed) ────
    na_results_r <- reactive({
      r <- rv$null_alleles_results
      shiny::validate(shiny::need(!is.null(r),
        paste0("No results yet. Go to the \"Null alleles\" module, choose your ",
               "per-locus coding, and click \"Compute + Bootstrap + Export\" first — ",
               "this module reuses those results directly (nothing is recomputed here).")))
      r
    })

    output$box_nloci <- renderValueBox({
      valueBox(length(na_results_r()$markers), "Loci", icon = icon("dna"), color = "navy")
    })
    output$box_npops <- renderValueBox({
      valueBox(length(na_results_r()$pops), "Populations", icon = icon("users"), color = "teal")
    })
    output$box_fstena <- renderValueBox({
      v <- round(na_results_r()$fst_global$global_ena, 4)
      col <- if (is.na(v)) "navy" else if (v > 0.15) "red" else if (v > 0.05) "yellow" else "green"
      valueBox(if (is.na(v)) "NA" else v, HTML("Global F<sub>ST</sub>-ENA"),
               icon = icon("chart-bar"), color = col)
    })
    output$box_nboot <- renderValueBox({
      r <- na_results_r()
      valueBox(r$nboot, "Bootstrap replicates (loci)", icon = icon("dice"), color = "purple")
    })

    output$ui_run_status <- renderUI({
      r <- tryCatch(na_results_r(), error = function(e) NULL)
      if (is.null(r)) return(NULL)
      tags$div(class = "na-info",
        icon("check-circle"), " ",
        sprintf("Using Null Alleles results: %d loci \u00b7 %d populations \u00b7 %d pairwise combinations.",
                length(r$markers), length(r$pops), nrow(r$fst_pair$long)))
    })

    # ── Population GPS centroids (needed for D_geo; IBD-specific) ───────────
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

    # ══════════════════════════════════════════════════════════════════════
    #  GEOGRAPHIC DISTANCE — Vincenty ellipsoidal geodesic (WGS84), in metres
    #  (validated against a FreeNA-style isolation-by-distance reference file:
    #  matches to the metre; a spherical Haversine approximation was off by
    #  10-20 m and a planar UTM distance was off by ~20 m over ~66 km)
    # ══════════════════════════════════════════════════════════════════════
    .vincenty_m <- function(lat1, lon1, lat2, lon2) {
      a <- 6378137.0; f <- 1/298.257223563; b <- (1 - f) * a
      L  <- (lon2 - lon1) * pi / 180
      U1 <- atan((1 - f) * tan(lat1 * pi / 180))
      U2 <- atan((1 - f) * tan(lat2 * pi / 180))
      sinU1 <- sin(U1); cosU1 <- cos(U1); sinU2 <- sin(U2); cosU2 <- cos(U2)
      lam <- L
      for (i in seq_len(200L)) {
        sinLam <- sin(lam); cosLam <- cos(lam)
        sinSigma <- sqrt((cosU2*sinLam)^2 + (cosU1*sinU2 - sinU1*cosU2*cosLam)^2)
        if (sinSigma == 0) return(0)
        cosSigma <- sinU1*sinU2 + cosU1*cosU2*cosLam
        sigma <- atan2(sinSigma, cosSigma)
        sinAlpha <- cosU1*cosU2*sinLam/sinSigma
        cosSqAlpha <- 1 - sinAlpha^2
        cos2SigmaM <- if (cosSqAlpha != 0) cosSigma - 2*sinU1*sinU2/cosSqAlpha else 0
        C <- f/16*cosSqAlpha*(4 + f*(4 - 3*cosSqAlpha))
        lamPrev <- lam
        lam <- L + (1 - C)*f*sinAlpha*(sigma + C*sinSigma*(cos2SigmaM + C*cosSigma*(-1 + 2*cos2SigmaM^2)))
        if (abs(lam - lamPrev) < 1e-12) break
      }
      uSq <- cosSqAlpha*(a^2 - b^2)/b^2
      A <- 1 + uSq/16384*(4096 + uSq*(-768 + uSq*(320 - 175*uSq)))
      B <- uSq/1024*(256 + uSq*(-128 + uSq*(74 - 47*uSq)))
      deltaSigma <- B*sinSigma*(cos2SigmaM + B/4*(cosSigma*(-1 + 2*cos2SigmaM^2) -
                    B/6*cos2SigmaM*(-3 + 4*sinSigma^2)*(-3 + 4*cos2SigmaM^2)))
      b * A * (sigma - deltaSigma)
    }

    # ══════════════════════════════════════════════════════════════════════
    #  FULL PAIRWISE TABLE — matches the reference layout exactly:
    #  Pop1, Pop2, D_geo, FST-FreeNA(+CI), ln(D_geo), F_R(+CI), D_CSE-INA, D_CSE
    #  Sourced ENTIRELY from rv$null_alleles_results — nothing recomputed here
    #  except D_geo/ln(D_geo) (which the Null Alleles module doesn't compute).
    # ══════════════════════════════════════════════════════════════════════
    .linearise <- function(x) { x <- pmin(pmax(x, 0), 0.9999); x / (1 - x) }

    full_pair_table_r <- reactive({
      na <- na_results_r()
      fst_long <- na$fst_pair$long                     # Pop1,Pop2,FST_raw,FST_ENA,Delta_FST
      dc_long  <- na$dc_pair$long                       # Pop1,Pop2,DCSE_raw,DCSE_INA,Delta_DCSE
      bf       <- na$boot_pair_fst                      # Pop1,Pop2,FST_ENA_obs,CI_lo_loci,...,FST_raw_obs,CI_lo_raw,CI_hi_raw
      bd       <- na$boot_pair_dc

      df <- merge(fst_long, dc_long[, c("Pop1","Pop2","DCSE_raw","DCSE_INA")],
                  by = c("Pop1","Pop2"), sort = FALSE)

      if (!is.null(bf) && nrow(bf) > 0L) {
        bf2 <- bf[, c("Pop1","Pop2","CI_lo_raw","CI_hi_raw","CI_lo_loci","CI_hi_loci")]
        names(bf2)[3:6] <- c("FST_raw_lo","FST_raw_hi","FST_ENA_lo","FST_ENA_hi")
        df <- merge(df, bf2, by = c("Pop1","Pop2"), sort = FALSE)
      } else {
        df$FST_raw_lo <- NA_real_; df$FST_raw_hi <- NA_real_
        df$FST_ENA_lo <- NA_real_; df$FST_ENA_hi <- NA_real_
      }
      if (!is.null(bd) && nrow(bd) > 0L) {
        bd2 <- bd[, c("Pop1","Pop2","CI_lo_raw","CI_hi_raw","CI_lo_loci","CI_hi_loci")]
        names(bd2)[3:6] <- c("DCSE_raw_lo","DCSE_raw_hi","DCSE_INA_lo","DCSE_INA_hi")
        df <- merge(df, bd2, by = c("Pop1","Pop2"), sort = FALSE)
      } else {
        df$DCSE_raw_lo <- NA_real_; df$DCSE_raw_hi <- NA_real_
        df$DCSE_INA_lo <- NA_real_; df$DCSE_INA_hi <- NA_real_
      }

      # Rousset's FR = FST/(1-FST) — reference "F_R" is based on FST-ENA
      # ("FST-FreeNA"); FR based on raw FST is kept too (model choice = "raw").
      df$FR        <- .linearise(df$FST_ENA)
      df$FR_lo     <- .linearise(df$FST_ENA_lo)
      df$FR_hi     <- .linearise(df$FST_ENA_hi)
      df$FR_raw    <- .linearise(df$FST_raw)
      df$FR_raw_lo <- .linearise(df$FST_raw_lo)
      df$FR_raw_hi <- .linearise(df$FST_raw_hi)

      # Geographic distance (Vincenty, metres), if GPS available
      coords <- tryCatch(coords_r(), error = function(e) NULL)
      if (!is.null(coords)) {
        get_d <- function(p1, p2) {
          c1 <- coords[coords$Population == p1, ]; c2 <- coords[coords$Population == p2, ]
          if (nrow(c1) >= 1L && nrow(c2) >= 1L)
            .vincenty_m(c1$Latitude[1L], c1$Longitude[1L], c2$Latitude[1L], c2$Longitude[1L])
          else NA_real_
        }
        df$Dgeo_m <- mapply(get_d, df$Pop1, df$Pop2)
        df$lnDgeo <- ifelse(df$Dgeo_m > 0, log(df$Dgeo_m), NA_real_)
      } else {
        df$Dgeo_m <- NA_real_; df$lnDgeo <- NA_real_
      }

      df
    })

    # ══════════════════════════════════════════════════════════════════════
    #  TAB 1 — Isolation by Distance (Rousset 1997)  [now the FIRST tab]
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
        any(is.finite(df$Dgeo_m)),
        "No geographic distances available. Set Latitude/Longitude at import for at least 2 populations."))

      use_log <- identical(input$ibd_model, "2D")
      x <- if (use_log) df$lnDgeo else df$Dgeo_m
      x_label <- if (use_log) "ln(D_geo)" else "D_geo (m)"

      if (identical(input$ibd_metric, "raw")) {
        y_avg <- df$FR_raw; y_lo <- df$FR_raw_lo; y_hi <- df$FR_raw_hi
        y_label <- "F_R (raw FST)"
      } else {
        y_avg <- df$FR; y_lo <- df$FR_lo; y_hi <- df$FR_hi
        y_label <- "F_R (FST-ENA)"
      }

      reg_avg <- .fit_line(y_avg, x)
      reg_lo  <- .fit_line(y_lo,  x)
      reg_hi  <- .fit_line(y_hi,  x)

      # Nb = 1/slope, Nem = Nb/(2*pi) — using the DISPLAYED (4-decimal-rounded)
      # slope, matching the reference tool's own convention.
      nbnem <- function(reg) {
        b <- round(reg$slope, 4)
        if (is.na(b) || b == 0) return(c(b = b, Nb = NA_real_, Nem = NA_real_))
        Nb <- 1 / b
        c(b = b, Nb = Nb, Nem = Nb / (2 * pi))
      }
      summ <- rbind(
        c(Line = "Average",  nbnem(reg_avg)),
        c(Line = "95%CI-i",  nbnem(reg_lo)),
        c(Line = "95%CI-s",  nbnem(reg_hi))
      )

      list(df = df, x = x, y_avg = y_avg, y_lo = y_lo, y_hi = y_hi,
           x_label = x_label, y_label = y_label,
           reg_avg = reg_avg, reg_lo = reg_lo, reg_hi = reg_hi,
           summary = summ, use_log = use_log, metric = input$ibd_metric)
    })

    # Full reference-style table: Pop1, Pop2, D_geo, FST-FreeNA(+CI),
    # ln(D_geo), F_R(+CI), D_CSE-INA, D_CSE — displayed completely, in the
    # same column order as the reference tool.
    output$dt_ibd_table <- DT::renderDT({
      r <- ibd_results_r()
      d <- r$df
      out <- data.frame(
        Farm1          = d$Pop1,
        Farm2          = d$Pop2,
        D_geo          = round(d$Dgeo_m, 4),
        `FST-FreeNA`   = round(d$FST_ENA, 6),
        `FST-FreeNA-i` = round(d$FST_ENA_lo, 6),
        `FST-FreeNA-s` = round(d$FST_ENA_hi, 6),
        `ln(D_geo)`    = round(d$lnDgeo, 6),
        F_R            = round(d$FR, 6),
        `F_R-i`        = round(d$FR_lo, 6),
        `F_R-s`        = round(d$FR_hi, 6),
        `D_CSE-INA`    = round(d$DCSE_INA, 6),
        D_CSE          = round(d$DCSE_raw, 6),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      DT::datatable(out, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 28, dom = "lrtip"),
        class = "compact stripe hover") |>
        DT::formatRound(c("D_geo","FST-FreeNA","FST-FreeNA-i","FST-FreeNA-s",
                           "ln(D_geo)","F_R","F_R-i","F_R-s","D_CSE-INA","D_CSE"), 6)
    })
    output$dl_ibd_csv <- downloadHandler(
      filename = function() paste0("IBD_pairwise_", Sys.Date(), ".csv"),
      content  = function(file) write.csv(ibd_results_r()$df, file, row.names = FALSE)
    )

    # Regression summary: slope (b) / Nb / Nem for the 3 fitted lines
    output$dt_ibd_reg <- DT::renderDT({
      r <- ibd_results_r()
      s <- as.data.frame(r$summary, stringsAsFactors = FALSE)
      s$b   <- round(as.numeric(s$b), 4)
      s$Nb  <- round(as.numeric(s$Nb), 4)
      s$Nem <- round(as.numeric(s$Nem), 4)
      names(s) <- c("slope", "b", "Nb", "Nem")
      DT::datatable(s, rownames = FALSE,
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

    # ══════════════════════════════════════════════════════════════════════
    #  TAB 2 — Mantel test (joint row/column permutation; rectangular-safe)
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

    # Generic Mantel permutation test: joint row/column relabelling of one
    # matrix (valid on rectangular/incomplete matrices too), Pearson r or
    # Rousset regression slope as the statistic.
    # p-value = (b+1)/(m+1)  [b = permuted statistics >= observed, m = total
    # valid permutations] — the standard correction that avoids ever reporting
    # p = 0 (Davison & Hinkley 1997; also the Fstat/RT convention this module
    # documents). Previously the code computed a plain proportion with no
    # +1/+1 correction — fixed here to match the method actually documented.
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
      m_valid  <- length(perm_fin)
      if (m_valid > 0L && is.finite(stat_obs)) {
        b_pos <- sum(perm_fin >= stat_obs)
        b_neg <- sum(perm_fin <= stat_obs)
        p_pos <- (b_pos + 1) / (m_valid + 1)
        p_neg <- (b_neg + 1) / (m_valid + 1)
      } else {
        p_pos <- NA_real_; p_neg <- NA_real_
      }
      lm0 <- tryCatch(lm(y_all[ok_obs] ~ x_all[ok_obs]), error = function(e) NULL)
      list(stat_obs = stat_obs, p_pos = p_pos, p_neg = p_neg, n_pairs = sum(ok_obs),
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
                  selected = .guess_col(cols, c("^FR$", "^FR_raw$", "FST_ENA", "DCSE_INA"),
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
               HTML("p-value<br>(one-sided, (b+1)/(m+1))"), icon = icon("check-circle"), color = col)
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
        sprintf("One-sided p (positive association, IBD) = %s",
                if (is.na(r$p_pos)) "NA" else formatC(r$p_pos, format = "f", digits = 4)), tags$br(),
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
