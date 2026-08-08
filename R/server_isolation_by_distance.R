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
      r <- tryCatch(na_results_r(), error = function(e) NULL)
      valueBox(if (is.null(r)) "\u2014" else length(r$markers), "Loci", icon = icon("dna"), color = "navy")
    })
    output$box_npops <- renderValueBox({
      r <- tryCatch(na_results_r(), error = function(e) NULL)
      valueBox(if (is.null(r)) "\u2014" else length(r$pops), "Populations", icon = icon("users"), color = "teal")
    })
    output$box_fstena <- renderValueBox({
      r <- tryCatch(na_results_r(), error = function(e) NULL)
      v <- if (is.null(r)) NA_real_ else round(r$fst_global$global_ena, 4)
      col <- if (is.na(v)) "navy" else if (v > 0.15) "red" else if (v > 0.05) "yellow" else "green"
      valueBox(if (is.na(v)) "\u2014" else v, HTML("Global F<sub>ST</sub>-ENA"),
               icon = icon("chart-bar"), color = col)
    })
    output$box_nboot <- renderValueBox({
      r <- tryCatch(na_results_r(), error = function(e) NULL)
      valueBox(if (is.null(r)) "\u2014" else r$nboot, "Bootstrap replicates (loci)", icon = icon("dice"), color = "purple")
    })

    output$ui_run_status <- renderUI({
      if (isTRUE(identical(input$ibd_source, "external"))) {
        return(tags$div(class = "na-info", icon("file-import"), " ",
          "Using a re-loaded external pairwise file \u2014 the Null Alleles module is not needed for this run."))
      }
      r <- tryCatch(na_results_r(), error = function(e) NULL)
      if (is.null(r)) return(NULL)
      tags$div(class = "na-info",
        icon("check-circle"), " ",
        sprintf("Using Null Alleles results: %d loci \u00b7 %d populations \u00b7 %d pairwise combinations.",
                length(r$markers), length(r$pops), nrow(r$fst_pair$long)))
    })

    # ══════════════════════════════════════════════════════════════════════
    #  OUTPUT FILE NAMING — same convention as the Null Alleles module: root
    #  auto-proposed from the imported data file's name (editable, never
    #  silently overwritten once the user has typed their own), + optional
    #  suffix. This module is meant to be usable on its own (e.g. re-loading
    #  a previously exported/edited pairwise file — see ibd_source/mt_source
    #  "external" options below), so the root also falls back gracefully
    #  when no data file was ever imported in this session.
    # ══════════════════════════════════════════════════════════════════════
    last_auto_root_ibd <- reactiveVal("")
    observeEvent(rv$dataset_filename, {
      fn <- rv$dataset_filename
      if (is.null(fn) || !nzchar(trimws(fn))) return(invisible(NULL))
      root_guess <- tools::file_path_sans_ext(basename(trimws(fn)))
      cur <- trimws(input$ibd_out_root %||% "")
      if (!nzchar(cur) || identical(cur, last_auto_root_ibd())) {
        updateTextInput(session, "ibd_out_root", value = root_guess)
        last_auto_root_ibd(root_guess)
      }
    }, ignoreInit = FALSE, ignoreNULL = TRUE)

    ibd_out_root_r   <- reactive({ r <- trimws(input$ibd_out_root %||% ""); if (nzchar(r)) r else "SPG_" })
    ibd_out_suffix_r <- reactive({ trimws(input$ibd_out_suffix %||% "") })
    ibd_out_filename <- function(desc) {
      suf <- ibd_out_suffix_r()
      paste0(ibd_out_root_r(), "-", desc, if (nzchar(suf)) paste0("-", suf) else "", ".txt")
    }

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
      if (isTRUE(identical(input$ibd_source, "external"))) {
        full_pair_table_external_r()
      } else {
        full_pair_table_internal_r()
      }
    })

    # ── EXTERNAL SOURCE: re-load a previously exported (and freely edited)
    #    pairwise file — e.g. the "pairwise_long_format" file exported by the
    #    Null Alleles module, or that same file hand-edited by the operator
    #    (rows removed, values corrected, extra columns added). This makes
    #    the IBD module fully standalone: it never has to touch the Null
    #    Alleles module in the same session. Column names matching the
    #    Null Alleles export (FST_raw, FST_ENA, DCSE_raw, DCSE_INA, their
    #    _lo/_hi CI bounds, FR/FR_raw + CI, Dgeo_m, lnDgeo) are used
    #    directly when present; anything missing is derived when possible
    #    (FR from FST via Rousset's linearisation, lnDgeo from Dgeo_m) or
    #    left NA otherwise (e.g. GPS-based Dgeo_m is only available in
    #    "internal" mode, since it depends on this session's imported data).
    full_pair_table_external_r <- reactive({
      shiny::req(input$ibd_ext_file)
      ext <- .mt_read_file(input$ibd_ext_file, input$ibd_ext_sep, input$ibd_ext_header)
      nm <- names(ext)
      pop1_col <- .guess_col(nm, c("^Pop1$", "^Farm1$", "^ID1$"), nm[1])
      pop2_col <- .guess_col(nm, c("^Pop2$", "^Farm2$", "^ID2$"), nm[2])
      shiny::validate(shiny::need(!identical(pop1_col, pop2_col),
        "Could not identify two distinct Pop1/Pop2 columns in the uploaded file."))

      df <- data.frame(Pop1 = as.character(ext[[pop1_col]]),
                        Pop2 = as.character(ext[[pop2_col]]),
                        stringsAsFactors = FALSE)

      num_col <- function(pats) {
        c <- .guess_col(nm, pats, NA_character_)
        if (is.na(c) || !(c %in% nm)) rep(NA_real_, nrow(ext))
        else suppressWarnings(as.numeric(ext[[c]]))
      }
      df$FST_raw     <- num_col(c("^FST_raw$"))
      df$FST_ENA     <- num_col(c("^FST_ENA$"))
      df$FST_raw_lo  <- num_col(c("^FST_raw_CI_lo_loci$", "^FST_raw_lo$"))
      df$FST_raw_hi  <- num_col(c("^FST_raw_CI_hi_loci$", "^FST_raw_hi$"))
      df$FST_ENA_lo  <- num_col(c("^FST_ENA_CI_lo_loci$", "^FST_ENA_lo$"))
      df$FST_ENA_hi  <- num_col(c("^FST_ENA_CI_hi_loci$", "^FST_ENA_hi$"))
      df$DCSE_raw    <- num_col(c("^DCSE_raw$"))
      df$DCSE_INA    <- num_col(c("^DCSE_INA$"))
      df$DCSE_raw_lo <- num_col(c("^DCSE_raw_CI_lo_loci$", "^DCSE_raw_lo$"))
      df$DCSE_raw_hi <- num_col(c("^DCSE_raw_CI_hi_loci$", "^DCSE_raw_hi$"))
      df$DCSE_INA_lo <- num_col(c("^DCSE_INA_CI_lo_loci$", "^DCSE_INA_lo$"))
      df$DCSE_INA_hi <- num_col(c("^DCSE_INA_CI_hi_loci$", "^DCSE_INA_hi$"))

      # FR (Rousset's linearised FST): use the file's own FR columns if
      # present, else derive them from FST_raw/FST_ENA (+ CI).
      fr_col <- function(pats, fallback_from) {
        c <- .guess_col(nm, pats, NA_character_)
        if (!is.na(c) && c %in% nm) suppressWarnings(as.numeric(ext[[c]]))
        else .linearise(fallback_from)
      }
      df$FR        <- fr_col(c("^FR$"),        df$FST_ENA)
      df$FR_lo     <- fr_col(c("^FR_lo$"),     df$FST_ENA_lo)
      df$FR_hi     <- fr_col(c("^FR_hi$"),     df$FST_ENA_hi)
      df$FR_raw    <- fr_col(c("^FR_raw$"),    df$FST_raw)
      df$FR_raw_lo <- fr_col(c("^FR_raw_lo$"), df$FST_raw_lo)
      df$FR_raw_hi <- fr_col(c("^FR_raw_hi$"), df$FST_raw_hi)

      df$Dgeo_m <- num_col(c("^Dgeo_m$", "^D_geo$", "^Distance$"))
      lnd <- num_col(c("^lnDgeo$", "^ln\\(D_geo\\)$"))
      df$lnDgeo <- ifelse(is.finite(lnd), lnd,
                           ifelse(is.finite(df$Dgeo_m) & df$Dgeo_m > 0, log(df$Dgeo_m), NA_real_))
      df
    })

    # ── INTERNAL SOURCE (default): built entirely from the Null Alleles
    #    module's results shared via rv$null_alleles_results.
    #    Matches the reference layout exactly:
    #    Pop1, Pop2, D_geo, FST-FreeNA(+CI), ln(D_geo), F_R(+CI), D_CSE-INA, D_CSE
    full_pair_table_internal_r <- reactive({
      na <- na_results_r()
      fst_long <- na$fst_pair$long                     # Pop1,Pop2,FST_raw,FST_ENA
      dc_long  <- na$dc_pair$long                       # Pop1,Pop2,DCSE_raw,DCSE_INA
      bf       <- na$boot_pair_fst                      # Pop1,Pop2,FST_ENA_obs,FST_ENA_CI_lo_loci,...,FST_raw_obs,FST_raw_CI_lo_loci,FST_raw_CI_hi_loci
      bd       <- na$boot_pair_dc

      df <- merge(fst_long, dc_long[, c("Pop1","Pop2","DCSE_raw","DCSE_INA")],
                  by = c("Pop1","Pop2"), sort = FALSE)

      if (!is.null(bf) && nrow(bf) > 0L) {
        bf2 <- bf[, c("Pop1","Pop2","FST_raw_CI_lo_loci","FST_raw_CI_hi_loci",
                      "FST_ENA_CI_lo_loci","FST_ENA_CI_hi_loci")]
        names(bf2)[3:6] <- c("FST_raw_lo","FST_raw_hi","FST_ENA_lo","FST_ENA_hi")
        df <- merge(df, bf2, by = c("Pop1","Pop2"), sort = FALSE)
      } else {
        df$FST_raw_lo <- NA_real_; df$FST_raw_hi <- NA_real_
        df$FST_ENA_lo <- NA_real_; df$FST_ENA_hi <- NA_real_
      }
      if (!is.null(bd) && nrow(bd) > 0L) {
        bd2 <- bd[, c("Pop1","Pop2","DCSE_raw_CI_lo_loci","DCSE_raw_CI_hi_loci",
                      "DCSE_INA_CI_lo_loci","DCSE_INA_CI_hi_loci")]
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

      # Geographic / pairwise distance (D_geo): either the Vincenty distance
      # from GPS centroids, or an external Pop1/Pop2/Distance file (e.g. the
      # subsample-pairs template exported from the Subdivision module, edited
      # by the operator to keep/exclude pairs and fill in distances of any
      # kind — not necessarily geographic).
      use_external_dgeo <- isTRUE(identical(input$ibd_dgeo_source, "external"))

      if (use_external_dgeo) {
        shiny::req(input$ibd_dgeo_file)
        ext <- .mt_read_file(input$ibd_dgeo_file, input$ibd_dgeo_sep, input$ibd_dgeo_header)
        shiny::validate(shiny::need(ncol(ext) >= 3L,
          "External distance file must have at least 3 columns: Pop1, Pop2, Distance."))
        nm <- names(ext)
        pop1_col <- .guess_col(nm, c("^Pop1$", "^Farm1$", "^ID1$"), nm[1])
        pop2_col <- .guess_col(nm, c("^Pop2$", "^Farm2$", "^ID2$"), nm[2])
        dist_col <- .guess_col(nm, c("^Distance$", "^Dgeo", "^Dist$"), nm[3])

        key <- function(a, b) { a <- as.character(a); b <- as.character(b)
                                 ifelse(a <= b, paste(a, b, sep = "__"), paste(b, a, sep = "__")) }

        ext2 <- data.frame(
          .key     = key(ext[[pop1_col]], ext[[pop2_col]]),
          Dgeo_ext = suppressWarnings(as.numeric(ext[[dist_col]])),
          stringsAsFactors = FALSE
        )
        ext2 <- ext2[!duplicated(ext2$.key), , drop = FALSE]

        df$.key <- key(df$Pop1, df$Pop2)
        # Inner join: ONLY pairs present in the external file are kept — pairs
        # the operator deleted from the file are excluded from the analysis.
        df <- merge(df, ext2, by = ".key", sort = FALSE)
        df$.key <- NULL
        df$Dgeo_m <- df$Dgeo_ext
        df$Dgeo_ext <- NULL
        df$lnDgeo <- ifelse(is.finite(df$Dgeo_m) & df$Dgeo_m > 0, log(df$Dgeo_m), NA_real_)
      } else {
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
        "No distances available. Either set Latitude/Longitude at import for at least 2 populations (GPS mode), or upload a Pop1/Pop2/Distance file (external file mode)."))

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
    output$dl_ibd_txt <- downloadHandler(
      filename = function() ibd_out_filename("regression"),
      content  = function(file) {
        r <- ibd_results_r()
        s <- as.data.frame(r$summary, stringsAsFactors = FALSE)
        hdr <- c(
          "Isolation by Distance \u2014 Rousset (1997) regression",
          sprintf("Habitat model: %s", if (r$use_log) "2D (F_R ~ ln(D_geo))" else "1D (F_R ~ D_geo)"),
          sprintf("Genetic distance metric: %s", if (identical(r$metric, "raw")) "F_R (raw FST)" else "F_R (FST-ENA)"),
          sprintf("Data source: %s", if (isTRUE(identical(input$ibd_source, "external"))) "external re-loaded pairwise file" else "Null Alleles module (this session)"),
          sprintf("Slope (b) / Nb / Nem \u2014 average: b=%.6f Nb=%.6f Nem=%.6f",
                   r$reg_avg$slope, 1/r$reg_avg$slope, (1/r$reg_avg$slope)/(2*pi)),
          ""
        )
        con <- file(file, open = "w", encoding = "UTF-8"); on.exit(close(con))
        writeLines(hdr, con = con, useBytes = TRUE)
        writeLines("Regression summary (slope / b / Nb / Nem for average and CI bounds):", con = con)
        write.table(s, file = con, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
        writeLines("", con = con)
        writeLines("Full pairwise table:", con = con)
        write.table(r$df, file = con, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
      }
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
        if (stat == "b") unname(coef(lm(yy[ok] ~ xx[ok]))[2L])
        else if (stat == "spearman") suppressWarnings(cor(xx[ok], yy[ok], method = "spearman"))
        else suppressWarnings(cor(xx[ok], yy[ok]))
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
      res$stat_label <- switch(stat, b = "Slope b", spearman = "Spearman rho", "Pearson r")
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

    # ══════════════════════════════════════════════════════════════════════
    #  TAB 3 — Partial Mantel test for MORE THAN 2-3 matrices, i.e. multiple
    #  regression on distance matrices (MRM; Legendre, Lapointe & Casgrain
    #  1994; Lichstein 2007), the standard generalisation of the classic
    #  2-3-matrix partial Mantel test to up to 10 predictor matrices at once
    #  (Fstat 2.9.4 convention). Base R packages (ade4, vegan, ecodist) cap
    #  partial Mantel at 2-3 matrices and Pearson/Spearman/Kendall only.
    #  Reuses the SAME data source, column pickers and rectangular-matrix
    #  handling (joint row/column relabelling, valid on incomplete pairwise
    #  data) as the simple Mantel test above (Tab 2).
    #
    #  CAVEAT the module's info panel also states: Guillot & Rousset (2013)
    #  and Crabot et al. (2019, Methods Ecol Evol 10:532-540) showed that
    #  classic partial Mantel tests have inflated type I error when the
    #  matrices being partialled out (e.g. geographic distance) are
    #  themselves spatially autocorrelated — a limitation this MRM
    #  generalisation inherits, since it relies on the same joint-relabelling
    #  permutation scheme. Borcard & Legendre (2012, Ecology 93:1473-1481)
    #  found the plain (non-partial) Mantel test has acceptable power for
    #  most ecological applications, so a simple Mantel per predictor (Tab 2)
    #  remains a reasonable cross-check. Lisboa et al. (2014, PLoS ONE
    #  9(6):e101238) proposed the Procrustes association metric (PAM) as a
    #  less-controversial, more powerful alternative to partial Mantel;
    #  PAM is NOT implemented here but is flagged as a natural next step.
    # ══════════════════════════════════════════════════════════════════════

    output$pm_col_y_ui <- renderUI({
      df <- tryCatch(mt_base_df_r(), error = function(e) NULL)
      cols <- if (is.null(df)) character(0) else names(df)[sapply(df, is.numeric)]
      selectInput(session$ns("pm_col_y"), "Response (Y):", choices = cols,
                  selected = .guess_col(cols, c("^FR$", "^FST_ENA$", "^DCSE_INA$"),
                                        if (length(cols)) cols[1] else NULL))
    })
    output$pm_col_x_ui <- renderUI({
      df <- tryCatch(mt_base_df_r(), error = function(e) NULL)
      cols <- if (is.null(df)) character(0) else names(df)[sapply(df, is.numeric)]
      cols <- setdiff(cols, input$pm_col_y %||% "")
      selectInput(session$ns("pm_col_x"), "Predictors (X1\u2026X10) \u2014 pick up to 10:",
                  choices = cols, selected = NULL, multiple = TRUE)
    })

    .pm_std <- function(v) {
      s <- stats::sd(v, na.rm = TRUE)
      if (!is.finite(s) || s == 0) v - mean(v, na.rm = TRUE) else (v - mean(v, na.rm = TRUE)) / s
    }

    partial_mantel_result_r <- eventReactive(input$run_partial_mantel, {
      df <- mt_base_df_r()
      p1c <- input$mt_col_pop1; p2c <- input$mt_col_pop2
      ycol <- input$pm_col_y; xcols <- input$pm_col_x

      shiny::validate(
        shiny::need(length(xcols) >= 1L, "Choose at least one predictor matrix."),
        shiny::need(length(xcols) <= 10L, "Up to 10 predictor matrices are supported (Fstat convention)."),
        shiny::need(!(ycol %in% xcols), "Y cannot also be used as a predictor."),
        shiny::need(all(c(p1c, p2c, ycol, xcols) %in% names(df)), "Selected columns not found.")
      )

      all_labels <- sort(unique(c(as.character(df[[p1c]]), as.character(df[[p2c]]))))
      build <- function(valcol) {
        tmp <- data.frame(P1 = as.character(df[[p1c]]), P2 = as.character(df[[p2c]]),
                           V = suppressWarnings(as.numeric(df[[valcol]])))
        .mt_build_square(tmp, "P1", "P2", "V", all_labels)
      }
      mat_y  <- build(ycol)
      mats_x <- stats::setNames(lapply(xcols, build), xcols)

      common <- Reduce(intersect, c(list(rownames(mat_y)), lapply(mats_x, rownames)))
      shiny::validate(shiny::need(length(common) >= 4L,
        "Not enough sub-samples shared by Y and all predictor matrices."))

      n <- length(common)
      lower_idx <- which(lower.tri(matrix(TRUE, n, n)))
      mat_y_c <- mat_y[common, common, drop = FALSE]
      y_all   <- mat_y_c[lower_idx]
      x_all   <- lapply(mats_x, function(m) m[common, common, drop = FALSE][lower_idx])

      compute_ok <- function(yv) { ok <- is.finite(yv); for (xv in x_all) ok <- ok & is.finite(xv); ok }
      ok0 <- compute_ok(y_all)
      shiny::validate(shiny::need(sum(ok0) >= (length(xcols) + 3L),
        "Too few dyads with complete data across Y and all predictors for this many predictors."))

      n_perm   <- as.integer(input$pm_n_perm)
      use_std  <- isTRUE(input$pm_standardize)

      fit_once <- function(yv) {
        okk <- compute_ok(yv)
        if (sum(okk) < (length(xcols) + 3L)) return(NULL)
        d <- as.data.frame(x_all)[okk, , drop = FALSE]; names(d) <- xcols
        d$Y <- yv[okk]
        if (use_std) { d$Y <- .pm_std(d$Y); for (nmc in xcols) d[[nmc]] <- .pm_std(d[[nmc]]) }
        m <- tryCatch(lm(Y ~ ., data = d), error = function(e) NULL)
        if (is.null(m)) return(NULL)
        list(coef = coef(m)[-1L], r2 = summary(m)$r.squared)
      }

      obs <- fit_once(y_all)
      shiny::validate(shiny::need(!is.null(obs), "Could not fit the model (check for collinear predictors)."))

      withProgress(message = "Running partial Mantel (MRM)\u2026", value = 0.1, {
        perm_coef <- matrix(NA_real_, n_perm, length(xcols), dimnames = list(NULL, xcols))
        perm_r2   <- numeric(n_perm)
        report_every <- max(1L, round(n_perm / 20))
        for (b in seq_len(n_perm)) {
          perm <- sample.int(n)
          y_p <- mat_y_c[perm, perm, drop = FALSE][lower_idx]
          fp  <- fit_once(y_p)
          if (!is.null(fp)) { perm_coef[b, ] <- fp$coef[xcols]; perm_r2[b] <- fp$r2 }
          if (b %% report_every == 0L) setProgress(0.1 + 0.85 * b / n_perm)
        }
      })

      tbl <- do.call(rbind, lapply(xcols, function(nmc) {
        pc <- perm_coef[, nmc]; pc <- pc[is.finite(pc)]
        b_obs   <- unname(obs$coef[nmc])
        m_valid <- length(pc)
        p_pos <- if (m_valid > 0L) (sum(pc >= b_obs) + 1) / (m_valid + 1) else NA_real_
        p_neg <- if (m_valid > 0L) (sum(pc <= b_obs) + 1) / (m_valid + 1) else NA_real_
        data.frame(Predictor = nmc, Coefficient = b_obs, p_positive = p_pos, p_negative = p_neg,
                   stringsAsFactors = FALSE)
      }))

      r2v  <- perm_r2[is.finite(perm_r2)]
      p_r2 <- if (length(r2v) > 0L) (sum(r2v >= obs$r2) + 1) / (length(r2v) + 1) else NA_real_

      list(table = tbl, r2 = obs$r2, p_r2 = p_r2, n_dyads = sum(ok0), n_pops = n,
           standardized = use_std, y_label = ycol, x_labels = xcols)
    })

    output$dt_partial_mantel <- DT::renderDT({
      r <- partial_mantel_result_r()
      d <- r$table
      names(d) <- c("Predictor", if (r$standardized) "Standardized coef." else "Coefficient",
                     "p (positive assoc.)", "p (negative assoc.)")
      DT::datatable(d, rownames = FALSE,
        options = list(dom = "t", pageLength = 10, ordering = FALSE),
        class = "compact stripe hover") |>
        DT::formatRound(names(d)[2:4], 4)
    })

    output$ui_partial_mantel_summary <- renderUI({
      r <- partial_mantel_result_r()
      tags$div(style = "margin-top:8px; font-size:13px; color:#333;",
        sprintf("Full model: R\u00b2 = %.4f, p = %s (permutation test on R\u00b2)",
                r$r2, if (is.na(r$p_r2)) "NA" else formatC(r$p_r2, format = "f", digits = 4)), tags$br(),
        sprintf("Dyads used: %d \u2014 sub-samples in common: %d", r$n_dyads, r$n_pops), tags$br(),
        sprintf("Response: %s \u2014 Predictors: %s", r$y_label, paste(r$x_labels, collapse = ", "))
      )
    })

    output$dl_partial_mantel_csv <- downloadHandler(
      filename = function() paste0("partial_mantel_", Sys.Date(), ".csv"),
      content  = function(file) write.csv(partial_mantel_result_r()$table, file, row.names = FALSE)
    )

    # ══════════════════════════════════════════════════════════════════════
    #  TAB 4 — Procrustes analysis & the Procrustean Association Metric (PAM)
    #  Lisboa, Peres-Neto, Chaer, Jesus, Mitchell, Chapman & Berbara (2014,
    #  PLoS ONE 9(6):e101238) — the "less controversial" alternative to
    #  partial Mantel flagged in the previous pass. Unlike the Mantel tests
    #  above (which compare pairwise DISTANCE matrices directly), Procrustes
    #  compares the CONFIGURATIONS of populations in two ordination spaces
    #  and, crucially, yields a per-population residual (PAM) rather than a
    #  single global statistic — so it can show WHICH populations drive an
    #  overall (mis)match between e.g. genetic and geographic structure.
    #
    #  Implementation (base R only, no new package dependency):
    #  1. PCoA (classical scaling, cmdscale() with the Cailliez correction
    #     for negative eigenvalues) turns each of the two chosen pairwise
    #     matrices into an n_pops x k coordinate configuration.
    #  2. Procrustes superimposition (translation + optimal rotation +
    #     isotropic scaling) of the Y configuration onto X, via SVD — the
    #     standard algorithm (Gower 1971; Peres-Neto & Jackson 2001).
    #     Both configurations are first scaled to unit sum of squares so
    #     the fit statistic m2 (Gower's statistic) lies in [0, 1] and the
    #     trace correlation r = sqrt(1 - m2) is comparable across pairs of
    #     variables.
    #  3. PROTEST (Jackson 1995): permutation test of r by randomly
    #     reshuffling which Y-configuration row corresponds to which X-row,
    #     one-sided p = (n_perm_r_obs_or_higher + 1) / (n_perm + 1).
    #  4. PAM: per-population Euclidean residual between each X point and
    #     its fitted (rotated/scaled) Y point — the vector highlighted by
    #     Lisboa et al. as usable in further regression/ANOVA analyses.
    #
    #  Like the simple Mantel test, this works on the SAME data source as
    #  Tabs 2-3 (Pop1/Pop2 long-format table), but — unlike Mantel/partial
    #  Mantel, which tolerate incomplete/rectangular pairwise data — PCoA
    #  and Procrustes both need COMPLETE square matrices, so populations
    #  with any missing pairwise value against the common set are dropped
    #  (reported explicitly to the user rather than silently).
    # ══════════════════════════════════════════════════════════════════════

    output$pam_col_x_ui <- renderUI({
      df <- tryCatch(mt_base_df_r(), error = function(e) NULL)
      cols <- if (is.null(df)) character(0) else names(df)[sapply(df, is.numeric)]
      selectInput(session$ns("pam_col_x"), "Matrix X:", choices = cols,
                  selected = .guess_col(cols, c("^FST_ENA$", "^FR$"), if (length(cols)) cols[1] else NULL))
    })
    output$pam_col_y_ui <- renderUI({
      df <- tryCatch(mt_base_df_r(), error = function(e) NULL)
      cols <- if (is.null(df)) character(0) else names(df)[sapply(df, is.numeric)]
      cols <- setdiff(cols, input$pam_col_x %||% "")
      selectInput(session$ns("pam_col_y"), "Matrix Y:", choices = cols,
                  selected = .guess_col(cols, c("^Dgeo_m$", "^lnDgeo$", "^DCSE_INA$"),
                                        if (length(cols)) cols[1] else NULL))
    })

    # Complete-case pruning: iteratively drop the population with the most
    # missing pairwise values (against the current common set, across BOTH
    # matrices) until none remain — needed because PCoA/Procrustes require
    # complete square matrices, unlike the Mantel tests above.
    .pam_prune_complete <- function(matX, matY) {
      common <- intersect(rownames(matX), rownames(matY))
      repeat {
        if (length(common) < 4L) break
        mx <- matX[common, common, drop = FALSE]; diag(mx) <- 0
        my <- matY[common, common, drop = FALSE]; diag(my) <- 0
        miss <- rowSums(is.na(mx)) + rowSums(is.na(my))
        if (all(miss == 0)) break
        drop_i <- which.max(miss)
        common <- common[-drop_i]
      }
      list(common = common)
    }

    .pam_pcoa <- function(mat, labels, k) {
      m <- mat[labels, labels, drop = FALSE]; diag(m) <- 0
      m <- (m + t(m)) / 2  # enforce exact symmetry (numerical safety)
      k_eff <- min(k, length(labels) - 1L)
      pc <- tryCatch(stats::cmdscale(stats::as.dist(m), k = k_eff, eig = TRUE, add = TRUE),
                     error = function(e) NULL)
      if (is.null(pc)) return(NULL)
      pts <- pc$points
      if (ncol(pts) < k) pts <- cbind(pts, matrix(0, nrow(pts), k - ncol(pts)))
      rownames(pts) <- labels
      pts
    }

    # Procrustes rotation of Yc onto Xc (both pre-centred, pre-scaled to
    # unit sum of squares) via SVD — Gower (1971) / Peres-Neto & Jackson (2001).
    .procrustes_fit <- function(Xc, Yc) {
      s <- svd(crossprod(Yc, Xc))
      Rmat <- s$u %*% t(s$v)
      d_sum <- sum(s$d)
      Yfit <- Yc %*% Rmat
      list(Yfit = Yfit, r = d_sum, m2 = 1 - d_sum^2)
    }

    .procrustes_prepare <- function(pts) {
      c0 <- scale(pts, center = TRUE, scale = FALSE)
      ss <- sum(c0^2)
      if (!is.finite(ss) || ss <= 0) return(c0)
      c0 / sqrt(ss)
    }

    pam_result_r <- eventReactive(input$run_pam, {
      df <- mt_base_df_r()
      p1c <- input$mt_col_pop1; p2c <- input$mt_col_pop2
      xcol <- input$pam_col_x; ycol <- input$pam_col_y
      shiny::validate(
        shiny::need(all(c(p1c, p2c, xcol, ycol) %in% names(df)), "Selected columns not found."),
        shiny::need(xcol != ycol, "Matrix X and Matrix Y must differ.")
      )

      all_labels <- sort(unique(c(as.character(df[[p1c]]), as.character(df[[p2c]]))))
      build <- function(valcol) {
        tmp <- data.frame(P1 = as.character(df[[p1c]]), P2 = as.character(df[[p2c]]),
                           V = suppressWarnings(as.numeric(df[[valcol]])))
        .mt_build_square(tmp, "P1", "P2", "V", all_labels)
      }
      matX <- build(xcol); matY <- build(ycol)

      pr <- .pam_prune_complete(matX, matY)
      common <- pr$common
      shiny::validate(shiny::need(length(common) >= 4L,
        "Not enough populations with complete pairwise data in both matrices (need >= 4)."))
      n_dropped <- length(all_labels) - length(common)

      k <- max(1L, min(as.integer(input$pam_k), length(common) - 1L))
      Xpts <- .pam_pcoa(matX, common, k); Ypts <- .pam_pcoa(matY, common, k)
      shiny::validate(shiny::need(!is.null(Xpts) && !is.null(Ypts),
        "PCoA failed \u2014 check that both matrices contain valid numeric distances."))

      Xc <- .procrustes_prepare(Xpts); Yc <- .procrustes_prepare(Ypts)
      fit <- .procrustes_fit(Xc, Yc)

      pam_vec <- sqrt(rowSums((Xc - fit$Yfit)^2))
      names(pam_vec) <- common

      n_perm <- as.integer(input$pam_n_perm)
      withProgress(message = "Running PROTEST permutation test\u2026", value = 0.1, {
        r_perm <- vapply(seq_len(n_perm), function(b) {
          perm <- sample.int(length(common))
          .procrustes_fit(Xc, Yc[perm, , drop = FALSE])$r
        }, numeric(1L))
        setProgress(1)
      })
      r_perm <- r_perm[is.finite(r_perm)]
      p_val <- if (length(r_perm) > 0L) (sum(r_perm >= fit$r) + 1) / (length(r_perm) + 1) else NA_real_

      list(pam = pam_vec, m2 = fit$m2, r = fit$r, p_value = p_val,
           n_pops = length(common), n_dropped = n_dropped, k = k,
           common = common, x_label = xcol, y_label = ycol, r_perm = r_perm)
    })

    output$box_pam_r <- renderValueBox({
      r <- pam_result_r()
      valueBox(round(r$r, 4), HTML("Trace correlation<br>r = sqrt(1-m\u00b2)"),
               icon = icon("link"), color = "purple")
    })
    output$box_pam_m2 <- renderValueBox({
      r <- pam_result_r()
      valueBox(round(r$m2, 4), HTML("Gower's statistic<br>m\u00b2 (lower = better fit)"),
               icon = icon("compress-arrows-alt"), color = "blue")
    })
    output$box_pam_p <- renderValueBox({
      r <- pam_result_r(); pv <- r$p_value
      col <- if (is.na(pv)) "yellow" else if (pv < 0.05) "green" else if (pv < 0.10) "yellow" else "red"
      valueBox(if (is.na(pv)) "NA" else formatC(pv, format = "f", digits = 4),
               HTML("PROTEST p-value<br>(one-sided)"), icon = icon("check-circle"), color = col)
    })
    output$box_pam_n <- renderValueBox({
      r <- pam_result_r()
      valueBox(r$n_pops, sprintf("Populations used (%d dropped: incomplete data)", r$n_dropped),
               icon = icon("users"), color = "teal")
    })

    output$dt_pam <- DT::renderDT({
      r <- pam_result_r()
      d <- data.frame(Population = names(r$pam), PAM = round(unname(r$pam), 6), stringsAsFactors = FALSE)
      d <- d[order(-d$PAM), , drop = FALSE]
      DT::datatable(d, rownames = FALSE,
        options = list(dom = "t", pageLength = nrow(d), ordering = FALSE),
        class = "compact stripe hover") |>
        DT::formatStyle("PAM", background = DT::styleColorBar(range(d$PAM), "#fca5a5"))
    })

    output$pam_plot <- plotly::renderPlotly({
      r <- pam_result_r()
      d <- data.frame(Population = names(r$pam), PAM = unname(r$pam))
      d <- d[order(-d$PAM), , drop = FALSE]
      plotly::plot_ly(d, x = ~reorder(Population, -PAM), y = ~PAM, type = "bar",
                       marker = list(color = "#B45309")) %>%
        plotly::layout(xaxis = list(title = "", tickangle = -45),
                       yaxis = list(title = sprintf("PAM (%s vs %s, %d axes)", r$x_label, r$y_label, r$k)),
                       margin = list(b = 90))
    })

    output$ui_pam_interpretation <- renderUI({
      r <- pam_result_r()
      tags$div(style = "margin-top:8px; font-size:13px; color:#333;",
        sprintf("Higher PAM = that population's relative position differs more between the %s and %s configurations.",
                r$x_label, r$y_label), tags$br(),
        "The PAM values above are a plain numeric vector (one row per population) \u2014 copy them ",
        "into a spreadsheet to regress against covariates (e.g. herd management, host breed) or run an ",
        "ANOVA by group, exactly as illustrated in Lisboa et al. (2014)."
      )
    })

    output$dl_pam_csv <- downloadHandler(
      filename = function() paste0("PAM_", Sys.Date(), ".csv"),
      content  = function(file) {
        r <- pam_result_r()
        write.csv(data.frame(Population = names(r$pam), PAM = unname(r$pam)), file, row.names = FALSE)
      }
    )

    # ══════════════════════════════════════════════════════════════════════
    #  TAB 5a — Mantel correlogram (Borcard & Legendre 2012, Ecology 93:1473-1481;
    #  Sokal 1986; Oden & Sokal 1986). Splits the chosen distance variable
    #  (e.g. D_geo) into equal-width classes; for each class, tests whether
    #  the response matrix (e.g. FST-ENA) is more similar within that class
    #  than among classes, using a Mantel test against a 0/1 "same class"
    #  model matrix. One-tailed in the direction of the observed sign (the
    #  paper's own convention, since the first classes are usually expected
    #  a priori to show positive spatial correlation and later ones negative);
    #  no correction for multiple testing across classes, matching the paper.
    #  Reuses the exact same .mt_mantel_matrix() permutation machinery (and
    #  therefore the same rectangular/incomplete-matrix tolerance) as the
    #  simple Mantel test above.
    # ══════════════════════════════════════════════════════════════════════

    output$corr_col_y_ui <- renderUI({
      df <- tryCatch(mt_base_df_r(), error = function(e) NULL)
      cols <- if (is.null(df)) character(0) else names(df)[sapply(df, is.numeric)]
      selectInput(session$ns("corr_col_y"), "Response matrix (Y):", choices = cols,
                  selected = .guess_col(cols, c("^FST_ENA$", "^FR$"), if (length(cols)) cols[1] else NULL))
    })
    output$corr_col_d_ui <- renderUI({
      df <- tryCatch(mt_base_df_r(), error = function(e) NULL)
      cols <- if (is.null(df)) character(0) else names(df)[sapply(df, is.numeric)]
      selectInput(session$ns("corr_col_d"), "Distance classes from:", choices = cols,
                  selected = .guess_col(cols, c("^Dgeo_m$", "^lnDgeo$"), if (length(cols)) cols[1] else NULL))
    })

    correlogram_result_r <- eventReactive(input$run_correlogram, {
      df <- mt_base_df_r()
      p1c <- input$mt_col_pop1; p2c <- input$mt_col_pop2
      ycol <- input$corr_col_y; dcol <- input$corr_col_d
      shiny::validate(
        shiny::need(all(c(p1c, p2c, ycol, dcol) %in% names(df)), "Selected columns not found."),
        shiny::need(ycol != dcol, "Response and distance columns must differ.")
      )

      all_labels <- sort(unique(c(as.character(df[[p1c]]), as.character(df[[p2c]]))))
      build <- function(valcol) {
        tmp <- data.frame(P1 = as.character(df[[p1c]]), P2 = as.character(df[[p2c]]),
                           V = suppressWarnings(as.numeric(df[[valcol]])))
        .mt_build_square(tmp, "P1", "P2", "V", all_labels)
      }
      matY <- build(ycol); matD <- build(dcol)

      n_classes <- max(2L, as.integer(input$corr_n_classes))
      d_vals <- matD[lower.tri(matD)]
      d_vals <- d_vals[is.finite(d_vals)]
      shiny::validate(shiny::need(length(d_vals) >= n_classes,
        "Not enough valid pairwise distances to build that many classes."))
      brks <- seq(min(d_vals), max(d_vals), length.out = n_classes + 1L)
      brks[1] <- brks[1] - 1e-9; brks[length(brks)] <- brks[length(brks)] + 1e-9

      n_perm <- as.integer(input$corr_n_perm); alpha <- as.numeric(input$corr_alpha)

      withProgress(message = "Computing Mantel correlogram\u2026", value = 0, {
        rows <- vector("list", n_classes)
        for (k in seq_len(n_classes)) {
          in_class <- matD > brks[k] & matD <= brks[k + 1L]
          in_class[is.na(in_class)] <- FALSE
          modelD <- matrix(0, nrow(matD), ncol(matD), dimnames = dimnames(matD))
          modelD[in_class] <- 1
          diag(modelD) <- 0
          res <- .mt_mantel_matrix(matY, modelD, n_perm = n_perm, stat = "r")
          p_one <- if (is.na(res$stat_obs)) NA_real_ else if (res$stat_obs >= 0) res$p_pos else res$p_neg
          rows[[k]] <- data.frame(
            Class     = k,
            D_lo      = round(brks[k] + 1e-9, 6),
            D_hi      = round(brks[k + 1L] - if (k == n_classes) 1e-9 else 0, 6),
            D_mean    = round(mean(d_vals[d_vals > brks[k] & d_vals <= brks[k + 1L]]), 6),
            n_pairs   = sum(in_class) / 2,
            r_Mantel  = round(res$stat_obs, 6),
            p_value   = round(p_one, 6),
            significant = !is.na(p_one) && p_one < alpha,
            stringsAsFactors = FALSE
          )
          setProgress(k / n_classes)
        }
      })
      tbl <- do.call(rbind, rows)
      list(table = tbl, y_label = ycol, d_label = dcol, alpha = alpha)
    })

    output$dt_correlogram <- DT::renderDT({
      r <- correlogram_result_r()
      d <- r$table; names(d)[names(d) == "significant"] <- sprintf("p < %.3g", r$alpha)
      DT::datatable(d, rownames = FALSE,
        options = list(dom = "t", pageLength = nrow(d), ordering = FALSE),
        class = "compact stripe hover")
    })

    output$correlogram_plot <- plotly::renderPlotly({
      r <- correlogram_result_r(); d <- r$table
      cols <- ifelse(d$significant, "#B45309", "rgba(180,83,9,0.25)")
      plotly::plot_ly(d, x = ~D_mean, y = ~r_Mantel, type = "scatter", mode = "lines+markers",
                       marker = list(size = 12, color = cols, line = list(color = "#B45309", width = 1.5)),
                       line = list(color = "#B45309", width = 1)) %>%
        plotly::layout(
          shapes = list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 0, y1 = 0,
                              line = list(color = "grey", width = 1, dash = "dot"))),
          xaxis = list(title = r$d_label), yaxis = list(title = sprintf("Mantel r (%s)", r$y_label)),
          title = list(text = "Filled = significant, open = not significant", font = list(size = 11)))
    })

    output$dl_correlogram_csv <- downloadHandler(
      filename = function() paste0("mantel_correlogram_", Sys.Date(), ".csv"),
      content  = function(file) write.csv(correlogram_result_r()$table, file, row.names = FALSE)
    )

    # ══════════════════════════════════════════════════════════════════════
    #  TAB 5b — Spatially-constrained (MSR) Mantel test
    #  Wagner & Dray (2015, Methods Ecol Evol 6:1169-1178); Crabot, Clappe,
    #  Dray & Datry (2019, Methods Ecol Evol 10:532-540). Addresses the
    #  inflation of the classic (and partial) Mantel test's type I error
    #  rate when BOTH matrices are spatially autocorrelated (Guillot &
    #  Rousset 2013) — the caveat flagged for the partial-Mantel tab.
    #
    #  Procedure (Crabot et al. 2019, section 2.3-ii, "data obtained
    #  directly as distances"): (1) PCoA of DX -> X_PCO configuration;
    #  (2) build a spatial weighting matrix W and its Moran's Eigenvector
    #  Maps (MEM); (3) generate nMSR random replicates of X_PCO that
    #  preserve its spatial autocorrelation structure (Moran spectral
    #  randomization) -> Euclidean distances D_X-MSR; (4) Mantel statistic
    #  between D_X-MSR and the FIXED DY, nMSR times, to build the null
    #  distribution and the corrected statistic r*_obs = robs - E(rM-MSR).
    #
    #  IMPLEMENTATION NOTES / SIMPLIFICATIONS (base R only, no spdep/
    #  adespatial dependency):
    #  - W is built here as a simple ROW-STANDARDIZED INVERSE-DISTANCE
    #    matrix from the chosen distance column (not a Gabriel/Delaunay/
    #    minimum-spanning-tree graph as in the paper's own comparisons) —
    #    a common, defensible default, but a different W choice can change
    #    results (the paper itself shows MSR Mantel is sensitive to W
    #    misspecification), so treat this as one reasonable option, not
    #    the only correct one.
    #  - MSR replicates are generated per PCoA axis by keeping each MEM
    #    coefficient's magnitude fixed (|a_i| = |cor(v_i, x)|) and applying
    #    independent random signs. This satisfies the paper's defining
    #    constraints exactly (identical Moran's I / multiscale structure,
    #    Eqs 5-7) but is a coarser (discrete, 2^(n-1) possible replicates)
    #    randomization than adespatial's continuous-rotation algorithm —
    #    adequate for a permutation test's null distribution, but not a
    #    byte-for-byte reimplementation of the original method.
    # ══════════════════════════════════════════════════════════════════════

    output$msr_col_x_ui <- renderUI({
      df <- tryCatch(mt_base_df_r(), error = function(e) NULL)
      cols <- if (is.null(df)) character(0) else names(df)[sapply(df, is.numeric)]
      selectInput(session$ns("msr_col_x"), "Matrix X (will be spatially randomized):", choices = cols,
                  selected = .guess_col(cols, c("^FST_ENA$", "^FR$"), if (length(cols)) cols[1] else NULL))
    })
    output$msr_col_y_ui <- renderUI({
      df <- tryCatch(mt_base_df_r(), error = function(e) NULL)
      cols <- if (is.null(df)) character(0) else names(df)[sapply(df, is.numeric)]
      cols <- setdiff(cols, input$msr_col_x %||% "")
      selectInput(session$ns("msr_col_y"), "Matrix Y (fixed):", choices = cols,
                  selected = .guess_col(cols, c("^Dgeo_m$", "^lnDgeo$", "^DCSE_INA$"),
                                        if (length(cols)) cols[1] else NULL))
    })
    output$msr_col_w_ui <- renderUI({
      df <- tryCatch(mt_base_df_r(), error = function(e) NULL)
      cols <- if (is.null(df)) character(0) else names(df)[sapply(df, is.numeric)]
      selectInput(session$ns("msr_col_w"), "Distance used to build W (spatial weights):", choices = cols,
                  selected = .guess_col(cols, c("^Dgeo_m$"), if (length(cols)) cols[1] else NULL))
    })

    # Moran's Eigenvector Maps from a row-standardized inverse-distance W
    .mem_from_dist <- function(matD, labels) {
      d <- matD[labels, labels, drop = FALSE]; diag(d) <- 0
      W <- ifelse(is.finite(d) & d > 0, 1 / d, 0)
      rs <- rowSums(W); rs[rs == 0] <- 1
      W <- W / rs                                    # row-standardized
      n <- length(labels)
      Wc <- (diag(n) - matrix(1/n, n, n)) %*% W %*% (diag(n) - matrix(1/n, n, n))
      Wc <- (Wc + t(Wc)) / 2                          # numerical symmetry
      e <- eigen(Wc, symmetric = TRUE)
      keep <- which(abs(e$values) > 1e-8)
      keep <- keep[seq_len(min(length(keep), n - 1L))]
      list(V = e$vectors[, keep, drop = FALSE], lambda = e$values[keep])
    }

    .msr_replicate <- function(x, V) {
      xc <- x - mean(x)
      ss <- sum(xc^2)
      if (!is.finite(ss) || ss == 0) return(x)
      corr <- as.numeric(crossprod(V, xc)) / sqrt(colSums(V^2) * ss)   # cor(v_i, x)
      signs <- sample(c(-1, 1), length(corr), replace = TRUE)
      a <- corr * signs
      scl <- sqrt(ss) / sqrt(sum(a^2))                # rescale so xMSR has the same variance as x
      mean(x) + as.numeric(V %*% (a * scl))
    }

    .mantel_stat_only <- function(mat1, mat2) {
      common <- intersect(rownames(mat1), rownames(mat2))
      if (length(common) < 3L) return(NA_real_)
      m1 <- mat1[common, common, drop = FALSE]; m2 <- mat2[common, common, drop = FALSE]
      idx <- which(lower.tri(m1)); xx <- m1[idx]; yy <- m2[idx]
      ok <- is.finite(xx) & is.finite(yy)
      if (sum(ok) < 3L) return(NA_real_)
      suppressWarnings(cor(xx[ok], yy[ok]))
    }

    msr_mantel_result_r <- eventReactive(input$run_msr, {
      df <- mt_base_df_r()
      p1c <- input$mt_col_pop1; p2c <- input$mt_col_pop2
      xcol <- input$msr_col_x; ycol <- input$msr_col_y; wcol <- input$msr_col_w
      shiny::validate(
        shiny::need(all(c(p1c, p2c, xcol, ycol, wcol) %in% names(df)), "Selected columns not found."),
        shiny::need(xcol != ycol, "Matrix X and Matrix Y must differ.")
      )

      all_labels <- sort(unique(c(as.character(df[[p1c]]), as.character(df[[p2c]]))))
      build <- function(valcol) {
        tmp <- data.frame(P1 = as.character(df[[p1c]]), P2 = as.character(df[[p2c]]),
                           V = suppressWarnings(as.numeric(df[[valcol]])))
        .mt_build_square(tmp, "P1", "P2", "V", all_labels)
      }
      matX <- build(xcol); matY <- build(ycol); matW <- build(wcol)

      pr <- .pam_prune_complete(matX, matY)
      common <- intersect(pr$common, .pam_prune_complete(matX, matW)$common)
      shiny::validate(shiny::need(length(common) >= 5L,
        "Not enough populations with complete data across X, Y and the spatial-weights distance (need >= 5)."))
      n_dropped <- length(all_labels) - length(common)

      k <- max(1L, min(as.integer(input$msr_k), length(common) - 2L))
      Xpts <- .pam_pcoa(matX, common, k)
      shiny::validate(shiny::need(!is.null(Xpts), "PCoA on Matrix X failed."))

      mem <- .mem_from_dist(matW, common)
      shiny::validate(shiny::need(ncol(mem$V) >= 1L, "Could not build spatial eigenvectors (MEM) from the chosen distance."))

      matY_c <- matY[common, common, drop = FALSE]; diag(matY_c) <- 0
      Dx_obs <- as.matrix(stats::dist(Xpts))
      r_obs  <- .mantel_stat_only(Dx_obs, matY_c)

      n_msr <- as.integer(input$msr_n_perm)
      withProgress(message = "Running MSR-Mantel\u2026", value = 0, {
        r_msr <- vapply(seq_len(n_msr), function(b) {
          Xmsr <- apply(Xpts, 2, .msr_replicate, V = mem$V)
          rownames(Xmsr) <- common
          Dmsr <- as.matrix(stats::dist(Xmsr))
          .mantel_stat_only(Dmsr, matY_c)
        }, numeric(1L))
        setProgress(1)
      })
      r_msr <- r_msr[is.finite(r_msr)]
      e_rmsr <- if (length(r_msr) > 0L) mean(r_msr) else NA_real_
      r_corrected <- r_obs - e_rmsr
      p_val <- if (length(r_msr) > 0L) (sum(r_msr >= r_obs) + 1) / (length(r_msr) + 1) else NA_real_

      list(r_obs = r_obs, e_rmsr = e_rmsr, r_corrected = r_corrected, p_value = p_val,
           r_msr = r_msr, n_pops = length(common), n_dropped = n_dropped, k = k,
           x_label = xcol, y_label = ycol, w_label = wcol)
    })

    output$box_msr_robs <- renderValueBox({
      r <- msr_mantel_result_r()
      valueBox(round(r$r_obs, 4), "Observed Mantel r", icon = icon("chart-line"), color = "purple")
    })
    output$box_msr_corrected <- renderValueBox({
      r <- msr_mantel_result_r()
      valueBox(round(r$r_corrected, 4), HTML("Corrected r*<br>(r_obs \u2212 E[r_MSR])"),
               icon = icon("map-marker-alt"), color = "blue")
    })
    output$box_msr_p <- renderValueBox({
      r <- msr_mantel_result_r(); pv <- r$p_value
      col <- if (is.na(pv)) "yellow" else if (pv < 0.05) "green" else if (pv < 0.10) "yellow" else "red"
      valueBox(if (is.na(pv)) "NA" else formatC(pv, format = "f", digits = 4),
               "MSR p-value (one-sided)", icon = icon("check-circle"), color = col)
    })
    output$box_msr_n <- renderValueBox({
      r <- msr_mantel_result_r()
      valueBox(r$n_pops, sprintf("Populations used (%d dropped: incomplete data)", r$n_dropped),
               icon = icon("users"), color = "teal")
    })

    output$msr_hist <- plotly::renderPlotly({
      r <- msr_mantel_result_r()
      shiny::req(length(r$r_msr) > 0L)
      plotly::plot_ly() %>%
        plotly::add_histogram(x = r$r_msr, nbinsx = 40,
          marker = list(color = "rgba(180,83,9,0.55)", line = list(color = "rgba(180,83,9,1)", width = 0.4))) %>%
        plotly::layout(
          shapes = list(
            list(type = "line", x0 = r$r_obs, x1 = r$r_obs, y0 = 0, y1 = 1, yref = "paper",
                 line = list(color = "#B40F20", width = 2, dash = "dash")),
            list(type = "line", x0 = r$e_rmsr, x1 = r$e_rmsr, y0 = 0, y1 = 1, yref = "paper",
                 line = list(color = "#333a43", width = 2, dash = "dot"))),
          xaxis = list(title = "Mantel r (MSR replicates)"), yaxis = list(title = "Count"),
          title = list(text = "Red dashed = observed r \u00b7 Grey dotted = E(r_MSR)", font = list(size = 11)),
          margin = list(t = 40), showlegend = FALSE)
    })

    output$ui_msr_summary <- renderUI({
      r <- msr_mantel_result_r()
      tags$div(style = "margin-top:8px; font-size:13px; color:#333;",
        sprintf("X = %s (spatially randomized, k=%d PCoA axes) \u2014 Y = %s (fixed) \u2014 W built from %s.",
                r$x_label, r$k, r$y_label, r$w_label), tags$br(),
        "If r_obs sits well outside the MSR null distribution (histogram), the association between X and Y ",
        "is unlikely to be a spurious artefact of both matrices sharing spatial structure."
      )
    })

  })
}
