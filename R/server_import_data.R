

server_import_data <- function(id, rv) {
  
  shiny::moduleServer(id, function(input, output, session) {
    rv$preview_raw  <- NULL
    rv$preview_meta <- NULL
    uploaded_file <- shiny::reactiveVal(NULL)
    rv$populationsLL_grouped <- shiny::reactiveVal(NULL)
    
    # ---------------------------------------------------------#
    # Helpers ####
    # ---------------------------------------------------------#
    show_map_connection_message <- function(
    text = "No internet connection detected. Leaflet could not download the map tiles."
    ) {
      leaflet::leafletProxy(session$ns("map")) %>%
        leaflet::clearControls() %>%
        leaflet::addControl(
          html = htmltools::HTML(
            paste0(
              "<div style='",
              "background: rgba(255,255,255,0.96);",
              "padding: 10px 12px;",
              "border-left: 4px solid #dc3545;",
              "border-radius: 6px;",
              "box-shadow: 0 1px 6px rgba(0,0,0,0.25);",
              "font-size: 13px;",
              "line-height: 1.35;",
              "max-width: 280px;",
              "'>",
              htmltools::htmlEscape(text),
              "</div>"
            )
          ),
          position = "topright",
          className = "leaflet-offline-warning"
        )
    }
    
    clear_map_connection_message <- function() {
      leaflet::leafletProxy(session$ns("map")) %>%
        leaflet::clearControls()
    }
    render_map_from_db <- function() {
      shiny::req(rv$con, rv$tbl_meta)
      
      ok_tbl <- tryCatch(DBI::dbExistsTable(rv$con, rv$tbl_meta), error = function(e) FALSE)
      if (!isTRUE(ok_tbl)) {
        rv$populationsLL_grouped(NULL)
        leaflet::leafletProxy(session$ns("map")) %>%
          leaflet::clearMarkers() %>%
          leaflet::clearShapes()
        return(invisible(NULL))
      }
      
      cols <- DBI::dbListFields(rv$con, rv$tbl_meta)
      if (!all(c("Population","Latitude","Longitude") %in% cols)) {
        rv$populationsLL_grouped(NULL)
        leaflet::leafletProxy("map") %>%
          leaflet::clearMarkers() %>%
          leaflet::clearShapes()
        return(invisible(NULL))
      }
      
      q <- sprintf("
    SELECT Population, Latitude, Longitude, COUNT(*) AS Population_size
    FROM %s
    WHERE Latitude IS NOT NULL AND Longitude IS NOT NULL
    GROUP BY Population, Latitude, Longitude
  ", .sql_ident(rv$tbl_meta))
      
      populationsLL_grouped <- .time_it("MAP query grouped points", {
        DBI::dbGetQuery(rv$con, q)
      })
      
      if (!nrow(populationsLL_grouped)) {
        rv$populationsLL_grouped(NULL)
        leaflet::leafletProxy("map") %>%
          leaflet::clearMarkers() %>%
          leaflet::clearShapes()
        return(invisible(NULL))
      }
      
      rv$populationsLL_grouped(populationsLL_grouped)
      
      .time_it("MAP leafletProxy update", {
        
        # ---- dynamic radius scaling in PIXELS 
        n <- populationsLL_grouped$Population_size
        n_sqrt <- sqrt(n)
        
        r_min <- 4
        r_max <- 18
        
        if (length(unique(n_sqrt)) > 1) {
          r_scaled <- r_min + (n_sqrt - min(n_sqrt)) /
            (max(n_sqrt) - min(n_sqrt)) * (r_max - r_min)
        } else {
          r_scaled <- rep((r_min + r_max) / 2, length(n_sqrt))
        }
        
        populationsLL_grouped$radius_px <- r_scaled
        
        proxy <- leaflet::leafletProxy(session$ns("map"), data = populationsLL_grouped) %>%
          leaflet::clearMarkers() %>%
          leaflet::clearShapes() %>%
          leaflet::addCircleMarkers(
            lng = ~Longitude, lat = ~Latitude,
            radius = ~radius_px,
            stroke = FALSE,
            fillOpacity = 0.7,
            popup = ~paste0("Location: ", htmltools::htmlEscape(Population), "<br>Population size: ", Population_size)
          )

        lng_rng <- range(populationsLL_grouped$Longitude, na.rm = TRUE)
        lat_rng <- range(populationsLL_grouped$Latitude,  na.rm = TRUE)
        
        if (all(is.finite(c(lng_rng, lat_rng)))) {
          proxy %>% leaflet::fitBounds(lng_rng[1], lat_rng[1], lng_rng[2], lat_rng[2])
        }
      })
      invisible(populationsLL_grouped)
    }
    
    # ---------------------------------------------------------#
    # Core formatting pipeline (shared by AUTO + MANUAL) ####
    # ---------------------------------------------------------#
    do_assign_and_format <- function(colnames_all,
                                     pop_data,
                                     col_ranges_data,
                                     metadata_ranges,
                                     missing_code,
                                     ploidy,
                                     latitude_data,
                                     longitude_data,
                                     selected_levels = character(0),
                                     make_map = TRUE,
                                     con = NULL,
                                     tbl_raw = NULL) {
      
      # -------------------------#
      # Missing handling in duckDB
      # -------------------------#
      miss_info   <- normalize_missing_code(missing_code)
      missing_gt  <- miss_info$miss_gt
      missing_set <- miss_info$miss_set
      
      # -------------------------#
      # Required inputs
      # -------------------------#
      if (is.null(pop_data) || length(pop_data) != 1L || is.na(pop_data) || !nzchar(pop_data)) {
        shinyalert::shinyalert("Error", "You need to select populations.", type = "error")
        return(invisible(FALSE))
      }
      
      if (is.null(col_ranges_data) || length(col_ranges_data) != 1L || is.na(col_ranges_data) || !nzchar(col_ranges_data)) {
        shinyalert::shinyalert("Error", "You need to select a marker range.", type = "error")
        return(invisible(FALSE))
      }
      if (is.null(colnames_all) || !length(colnames_all)) {
        shinyalert::shinyalert("Error", "No column names available (header not read).", type = "error")
        return(invisible(FALSE))
      }
      
      n_all <- length(colnames_all)
      
      # -------------------------#
      # Marker range validation
      # -------------------------#
      rvv <- unlist(strsplit(col_ranges_data, "[:-]"))
      rvv <- suppressWarnings(as.integer(rvv))
      
      if (length(rvv) != 2 || any(is.na(rvv)) || rvv[1] < 1 || rvv[2] > n_all) {
        shinyalert::shinyalert("Error", "Try again, your range is out of bounds.", type = "error")
        return(invisible(FALSE))
      }
      
      marker_idx <- seq.int(min(rvv), max(rvv))
      marker_raw_names <- colnames_all[marker_idx]
      
      # -------------------------#
      # GPS detection
      # -------------------------#
      has_gps <- !is.null(latitude_data)  && nzchar(latitude_data) &&
        !is.null(longitude_data) && nzchar(longitude_data) &&
        latitude_data  %in% colnames_all &&
        longitude_data %in% colnames_all &&
        latitude_data != pop_data &&
        longitude_data != pop_data
      
      if (has_gps && (latitude_data %in% marker_raw_names || longitude_data %in% marker_raw_names)) {
        shinyalert::shinyalert("Error", "Latitude/Longitude cannot be inside the marker range.", type = "error")
        return(invisible(FALSE))
      }
      
      # -------------------------#
      # Metadata ranges
      # -------------------------#
      meta_idx   <- parse_col_index_ranges(metadata_ranges, n_max = n_all)
      meta_names <- if (length(meta_idx)) colnames_all[meta_idx] else character(0)
      
      # forbid overlaps
      meta_names <- setdiff(meta_names, c(pop_data, latitude_data, longitude_data, marker_raw_names))
      
      # -------------------------#
      # Build keep columns (names)
      # -------------------------#
      keep <- c(pop_data)
      if (has_gps) keep <- c(keep, latitude_data, longitude_data)
      if (length(meta_names)) keep <- c(keep, meta_names)
      
      # all markers in preview (rows limited to preview_n by SQL LIMIT)
      marker_preview_names <- marker_raw_names
      keep <- unique(c(keep, marker_preview_names))
      keep <- intersect(keep, colnames_all)
      
      if (!length(keep)) {
        shinyalert::shinyalert("Error", "No columns selected (check separator/header).", type = "error")
        return(invisible(FALSE))
      }

      # -------------------------#
      # Read only required columns
      # -------------------------#
      if (is.null(con) || is.null(tbl_raw) || !nzchar(tbl_raw)) {
        shinyalert::shinyalert("Error", "Database not initialised.", type = "error")
        return(invisible(FALSE))
      }
      
      ok_raw <- tryCatch(DBI::dbExistsTable(con, tbl_raw), error = function(e) FALSE)
      if (!isTRUE(ok_raw)) {
        shinyalert::shinyalert("Error", "Raw table not found in DuckDB (import step missing).", type = "error")
        return(invisible(FALSE))
      }
      
      preview_n <- 100L
      
      # quote identifiers safely (columns + table)
      keep_sql <- paste(
        vapply(keep, function(x) as.character(DBI::dbQuoteIdentifier(con, x)), character(1)),
        collapse = ", "
      )
      tbl_sql <- as.character(DBI::dbQuoteIdentifier(con, tbl_raw))
      
      # IMPORTANT: bring rowid as stable individual id (DuckDB)
      sql_full <- sprintf("SELECT rowid AS individual, %s FROM %s;", keep_sql, tbl_sql)
      
      new_df <- .time_it("dbGetQuery selected cols -> R", {
        tryCatch(DBI::dbGetQuery(con, sql_full), error = function(e) NULL)
      })
      if (is.null(new_df) || !nrow(new_df)) {
        shinyalert::shinyalert("Error", "No rows returned from DuckDB for selected columns.", type = "error")
        return(invisible(FALSE))
      }
      
      new_df <- as.data.frame(new_df, stringsAsFactors = FALSE)
      
      
      
      # -------------------------#
      # Standardize key column names
      # -------------------------#
      names(new_df)[names(new_df) == pop_data] <- "Population"
      
      if (has_gps) {
        names(new_df)[names(new_df) == latitude_data]  <- "Latitude"
        names(new_df)[names(new_df) == longitude_data] <- "Longitude"
        new_df$Latitude  <- suppressWarnings(as.numeric(new_df$Latitude))
        new_df$Longitude <- suppressWarnings(as.numeric(new_df$Longitude))
      }
      
      # -------------------------#
      # Marker frame
      # -------------------------#
      marker_raw_names_ok <- marker_raw_names
      if (!length(marker_raw_names_ok)) {
        shinyalert::shinyalert("Error", "No marker columns found (check marker range).", type = "error")
        return(invisible(FALSE))
      }
      
      # All markers in preview (rows already limited to 100 by LIMIT clause)
      marker_preview_names <- marker_raw_names_ok
      marker_preview_ok <- intersect(marker_preview_names, names(new_df))
      if (!length(marker_preview_ok)) {
        shinyalert::shinyalert("Error", "Marker preview columns not found in selected columns (check marker range).", type = "error")
        return(invisible(FALSE))
      }
      
      marker_df <- new_df[, marker_preview_ok, drop = FALSE]
      
      # -------------------------#
      # PREVIEW ONLY (NO DECODING)
      # -------------------------#
      .preview_start <- Sys.time()
      
      meta_part <- .time_it("preview: extract meta part (R)", {
        new_df[, setdiff(names(new_df), marker_preview_ok), drop = FALSE]
      })
      
      marker_preview <- .time_it("preview: extract marker subset (R)", {
        new_df[, marker_preview_ok, drop = FALSE]
      })
      
      new_df_prev <- .time_it("preview: cbind meta + marker subset (R)", {
        out <- cbind(meta_part, marker_preview, stringsAsFactors = FALSE)
        
        # ---- merge paired allele columns into "a1/a2" for preview 
        # Example: B12 + B12_1  ->  B12 = "195/197"  (and remove B12_1)
        pick_b <- function(locus, nms) {
          cands <- c(paste0(locus, "_1"), paste0(locus, ".", 1:9))
          hit <- cands[cands %in% nms]
          if (length(hit)) hit[1] else NA_character_
        }
        
        nms <- names(out)
        loci <- unique(sub("(_1|\\.[0-9]+)$", "", nms))
        
        for (locus in loci) {
          a <- locus
          b <- pick_b(locus, nms)
          
          if (!a %in% nms) next
          if (is.na(b) || !b %in% nms) next
          
          out[[locus]] <- paste0(out[[a]], "/", out[[b]])
          
          # drop the allele columns (keep only merged locus)
          a_val <- as.character(out[[a]])
          b_val <- as.character(out[[b]])
          
          already_diploid <- grepl("/", a_val, fixed = TRUE) | grepl("-", a_val, fixed = TRUE)
          
          out[[a]] <- ifelse(
            already_diploid,
            a_val,                             
            paste0(a_val, "/", b_val)         
          )
          
          out[[b]] <- NULL 
          
          # refresh names after deletion
          nms <- names(out)
        }
        
        out
      })
      
      .time_it("preview: dbWriteTable formatted_preview", {
        DBI::dbWriteTable(con, rv$tbl_formatted_preview, new_df_prev, overwrite = TRUE)
      })
      
      base <- .time_it("base: infer from marker strings (R)", {
        infer_base_from_marker_strings(marker_df, sep = "/", default_base = 1000L)
      })
      haplotype_length <- suppressWarnings(as.integer(round(log10(as.numeric(base)))))
      if (!is.finite(haplotype_length) || haplotype_length < 1L) haplotype_length <- 3L
      
      .time_it("params: store downstream config", {
        .duckdb_set_params(con, list(
          is_preview       = 1L,
          preview_n        = as.integer(preview_n),
          pop_col_raw      = pop_data,
          ploidy           = as.integer(ploidy),
          missing_code_raw = as.character(missing_code),
          missing_gt       = 0L,
          missing_set      = paste(sort(unique(missing_set)), collapse = ","),
          sep              = "/",
          base             = as.integer(base),
          haplotype_length   = as.integer(haplotype_length),
          width_scalar_full  = as.integer(haplotype_length),
          base_scalar_full   = as.integer(base),
          marker_cols      = marker_raw_names_ok,   
          tbl_raw          = tbl_raw,
          tbl_meta         = rv$tbl_meta %||% "meta",
          tbl_formatted_preview = rv$tbl_formatted_preview %||% "formatted_preview",
          tbl_hf           = rv$tbl_hf %||% "hf",
          marker_cols_raw  = marker_raw_names_ok,                         # <-- CHANGE
          locus_cols       = .get_locus_cols_from_marker_cols(marker_raw_names_ok),
          genotype_format  = if (any(grepl("(_1|\\.[0-9]+)$", marker_raw_names_ok))) "paired" else "string"
        ))
      })
      
      .time_it("preview: formatted_preview from DB (LIMIT 100)", {
        rv$preview_raw <- DBI::dbGetQuery(
          con,
          sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_formatted_preview))
        )
      })
      
      rm(marker_df, meta_part, marker_preview, new_df_prev, new_df)
      gc(FALSE)
      
      invisible(TRUE)
    }
    
    
    # ---------------------------------------------------------#
    ## Session-scoped DuckDB ####
    # ---------------------------------------------------------#
    session_dir <- file.path(tempdir(), paste0("pgacmdr_", session$token))
    dir.create(session_dir, showWarnings = FALSE, recursive = TRUE)
    
    rv$tbl_raw               <- "raw"
    rv$tbl_meta              <- "meta"
    rv$tbl_formatted         <- "formatted"         # full (option B, later)
    rv$tbl_formatted_preview <- "formatted_preview" # UI-only
    rv$tbl_hf                <- "hf"
    rv$db_tick <- 0L
    
    session$userData$db_path <- file.path(session_dir, "pgacmdr.duckdb")
    
    con <- DBI::dbConnect(
      duckdb::duckdb(),
      dbdir = session$userData$db_path,
      read_only = FALSE
    )
    
    .duckdb_tune_for_big_import(con)
    
    rv$con <- con
    rv$db_ready <- TRUE

    # base map render 
    output$map <- leaflet::renderLeaflet({
      m <- leaflet::leaflet() %>%
        leaflet::addTiles()
      
      htmlwidgets::onRender(
        m,
        "
    function(el, x) {
      var map = this;
      var tileLayer = null;

      map.eachLayer(function(layer) {
        if (!tileLayer && layer instanceof L.TileLayer) {
          tileLayer = layer;
        }
      });

      if (!tileLayer) return;

      var loaded = false;
      var errorReported = false;

      function report(status) {
        if (HTMLWidgets.shinyMode) {
          Shiny.setInputValue(
            el.id + '_tile_status',
            { status: status, nonce: Date.now() },
            { priority: 'event' }
          );
        }
      }

      // Immediate browser offline check
      if (!navigator.onLine) {
        report('offline');
        return;
      }

      // Success
      tileLayer.on('load', function() {
        loaded = true;
        errorReported = false;
        report('loaded');
      });

      // Tile failure
      tileLayer.on('tileerror', function() {
        if (!errorReported) {
          errorReported = true;
          report('tileerror');
        }
      });

      // 1-second timeout test
      setTimeout(function() {
        if (!loaded) {
          report('timeout');
        }
      }, 1000);

      // React to network changes
      window.addEventListener('offline', function() {
        report('offline');
      });

      window.addEventListener('online', function() {
        loaded = false;
        errorReported = false;
        tileLayer.redraw();

        setTimeout(function() {
          if (!loaded) {
            report('timeout');
          }
        }, 1000);
      });
    }
    "
      )
    })
    shiny::observeEvent(input$map_tile_status, {
      status <- input$map_tile_status$status %||% ""
      
      if (status == "loaded") {
        clear_map_connection_message()
      } else if (status %in% c("offline", "tileerror", "timeout")) {
        show_map_connection_message(
          "Map unavailable: no internet connection or tile server unreachable."
        )
      }
    }, ignoreInit = TRUE)
    output$download_csv_transformed <- shiny::downloadHandler(
      filename = function() sprintf("ShinyPopGen_%s.csv", Sys.Date()),
      content = function(file) {
        shiny::req(rv$con)
        
        tbl <- if (DBI::dbExistsTable(rv$con, rv$tbl_formatted)) {
          rv$tbl_formatted
        } else if (DBI::dbExistsTable(rv$con, rv$tbl_formatted_preview)) {
          rv$tbl_formatted_preview
        } else {
          rv$tbl_raw
        }
        
        DBI::dbExecute(
          rv$con,
          sprintf(
            "COPY (SELECT * FROM %s) TO %s (HEADER, DELIMITER ',');",
            .sql_ident(tbl),
            DBI::dbQuoteString(rv$con, file)
          )
        )
      }
    )
    
    output$download_txt_transformed <- shiny::downloadHandler(
      filename = function() sprintf("ShinyPopGen_%s.txt", Sys.Date()),
      content = function(file) {
        shiny::req(rv$con)
        
        tbl <- if (DBI::dbExistsTable(rv$con, rv$tbl_formatted)) {
          rv$tbl_formatted
        } else if (DBI::dbExistsTable(rv$con, rv$tbl_formatted_preview)) {
          rv$tbl_formatted_preview
        } else {
          rv$tbl_raw
        }
        
        DBI::dbExecute(
          rv$con,
          sprintf(
            "COPY (SELECT * FROM %s) TO %s (HEADER, DELIMITER '\t');",
            .sql_ident(tbl),
            DBI::dbQuoteString(rv$con, file)
          )
        )
      }
    )
    
    output$download_map <- shiny::downloadHandler(
      filename = function() sprintf("ShinyPopGen_map_%s.png", Sys.Date()),
      content = function(file) {
        pops <- shiny::isolate(rv$populationsLL_grouped())
        shiny::req(!is.null(pops), nrow(pops) > 0)
        n      <- pops$Population_size
        n_sqrt <- sqrt(n)
        r_min  <- 4; r_max <- 18
        pops$radius_px <- if (length(unique(n_sqrt)) > 1)
          r_min + (n_sqrt - min(n_sqrt)) / (max(n_sqrt) - min(n_sqrt)) * (r_max - r_min)
        else rep((r_min + r_max) / 2, length(n_sqrt))
        m <- leaflet::leaflet(pops) %>%
          leaflet::addTiles() %>%
          leaflet::addCircleMarkers(
            lng = ~Longitude, lat = ~Latitude,
            radius = ~radius_px,
            stroke = FALSE, fillOpacity = 0.7,
            popup = ~paste0("Location: ", htmltools::htmlEscape(Population),
                            "<br>Population size: ", Population_size)
          ) %>%
          leaflet::fitBounds(
            min(pops$Longitude, na.rm = TRUE), min(pops$Latitude,  na.rm = TRUE),
            max(pops$Longitude, na.rm = TRUE), max(pops$Latitude,  na.rm = TRUE)
          )
        tmp_html <- tempfile(fileext = ".html")
        htmlwidgets::saveWidget(m, tmp_html, selfcontained = TRUE)
        webshot2::webshot(tmp_html, file = file, vwidth = 900, vheight = 600, delay = 1.5)
        unlink(tmp_html)
      }
    )

    session$onSessionEnded(function() {
      try(DBI::dbDisconnect(rv$con, shutdown = TRUE), silent = TRUE)
      try(duckdb::duckdb_shutdown(), silent = TRUE)
      try(unlink(session$userData$db_path, force = TRUE), silent = TRUE)
      if (dir.exists(session_dir)) unlink(session_dir, recursive = TRUE, force = TRUE)
    })
    
    # ---------------------------------------------------------#
    # Load_default_data ####
    # ---------------------------------------------------------#
    
    shiny::observeEvent(input$load_default_data, {
      
      reset_downstream_state(rv)
      
      default_path <- system.file("extdata", "default_dataset.csv", package = "shinypopgen")
      if (!nzchar(default_path) || !file.exists(default_path)) {
        shinyalert::shinyalert("Error", "Default dataset not found in package (inst/extdata/default_dataset.csv).", type = "error")
        return()
      }

      rv$file_path <- default_path
      rv$sep <- "\t"
      rv$header <- input$header
      
      # ------------------------------#
      # State 1: import RAW
      # ------------------------------#
      ok_db <- tryCatch(
        .duckdb_import_raw(
          con       = rv$con,
          tbl_raw   = rv$tbl_raw,
          file_path = rv$file_path,
          sep       = rv$sep,
          header    = rv$header
        ),
        error = function(e) {
          shinyalert::shinyalert("DuckDB import failed", conditionMessage(e), type = "error")
          FALSE
        }
      )
      if (!isTRUE(ok_db)) return()
      
      # ------------------------------#
      # State 2: header + raw preview from DuckDB
      # ------------------------------#
      rv$colnames_all <- tryCatch(DBI::dbListFields(rv$con, rv$tbl_raw), error = function(e) character(0))
      if (!length(rv$colnames_all)) {
        shinyalert::shinyalert("Header read failed", "No columns detected in DuckDB raw table.", type = "error")
        return()
      }
      
      update_metadata_choices(session, rv$colnames_all)
      
      rv$preview_raw <- tryCatch(
        DBI::dbGetQuery(rv$con, sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_raw))),
        error = function(e) NULL
      )
      if (is.null(rv$preview_raw) || !nrow(rv$preview_raw)) {
        shinyalert::shinyalert("No data rows", "Raw table has 0 rows after import.", type = "warning")
        return()
      }
      
      # ------------------------------#
      # State 3: auto-detect on raw preview
      # ------------------------------#
      miss_auto <- normalize_missing_code(input$missing_code)
      det <- detect_columns_auto(rv$preview_raw, missing_info = miss_auto)
      rv$det <- det
      
      # populate the manual panel selections
      populate_manual_from_detection(session, det, rv$colnames_all)
      
      # meta_ranges used BOTH for meta table + for auto-format call
      meta_idx <- match(det$metadata_cols, rv$colnames_all)
      meta_idx <- meta_idx[!is.na(meta_idx)]
      meta_ranges <- .compress_idx_ranges(meta_idx)
      
      # If population missing -> stop (manual UI is ready)
      if (is.null(det$population) || is.na(det$population) || !nzchar(det$population)) {
        shinyalert::shinyalert(
          "Auto-detect incomplete",
          "Population column not detected. Either change the separator or use the manual mode.",
          type = "warning"
        )
        return()
      }
      
      # ------------------------------#
      # State 3b: AUTO-FORMAT 
      # ------------------------------#
      if (!is.null(det$population) && nzchar(det$population) &&
          !is.null(det$marker_range) && nzchar(det$marker_range)) {
        
        ok_fmt <- do_assign_and_format(
          colnames_all    = rv$colnames_all,
          pop_data        = det$population,
          col_ranges_data = det$marker_range,
          metadata_ranges = meta_ranges,
          missing_code    = input$missing_code,
          ploidy          = as.numeric(input$ploidy),
          latitude_data   = det$latitude %|||% "",
          longitude_data  = det$longitude %|||% "",
          selected_levels = character(0),
          make_map        = FALSE,
          con             = rv$con,
          tbl_raw         = rv$tbl_raw
        )
        
        if (isTRUE(ok_fmt)) {
          
          # full DB: formatted is already written by do_assign_and_format()
          ok_fmt_tbl <- tryCatch(DBI::dbExistsTable(rv$con, rv$tbl_formatted_preview), error = function(e) FALSE)
          if (!isTRUE(ok_fmt_tbl)) {
            shinyalert::shinyalert("Error", "Formatted table not found in DuckDB.", type = "error")
            return()
          }
        }
      }
      
      # ------------------------------#
      # State 4: build META + preview + map (DB meta used for map)
      # ------------------------------#
      tryCatch(
        .time_it("META build", {
          .duckdb_build_meta(
            con = rv$con,
            tbl_raw = rv$tbl_raw,
            tbl_meta = rv$tbl_meta,
            colnames_all = rv$colnames_all,
            pop_data = det$population,
            latitude_data = det$latitude,
            longitude_data = det$longitude,
            metadata_ranges = meta_ranges
          )
        })
        ,
        error = function(e) {
          shinyalert::shinyalert("DuckDB meta build failed", conditionMessage(e), type = "error")
          NULL
        }
      )
      
      # ---- ensure meta.individual exists and matches raw.rowid 
      DBI::dbExecute(
        rv$con,
        sprintf(
          "
    ALTER TABLE %s ADD COLUMN IF NOT EXISTS individual BIGINT;

    UPDATE %s m
    SET individual = r.rowid
    FROM %s r
    WHERE m.rowid = r.rowid;
    ",
          .sql_ident(rv$tbl_meta),
          .sql_ident(rv$tbl_meta),
          .sql_ident(rv$tbl_raw)
        )
      )
      
    # ---- build HF (chunked) AFTER meta exists
      .build_hf_from_params <- function(con, rv, missing_code_raw, batch_size = 10000L) {
        shiny::req(con)
        
        # marker_cols_raw (raw columns, may include suffixes like _1 / .1)
        marker_json <- tryCatch(
          DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='marker_cols_raw'")$value[1],
          error = function(e) NA_character_
        )
        
        marker_cols_raw <- if (!is.na(marker_json) && nzchar(marker_json)) {
          tryCatch(jsonlite::fromJSON(marker_json), error = function(e) character(0))
        } else character(0)
        
        if (!length(marker_cols_raw)) stop("No marker_cols_raw found in DuckDB params.")
        
        # base
        base <- tryCatch(
          DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='base'")$value[1],
          error = function(e) NA_character_
        )
        base <- suppressWarnings(as.integer(base))
        if (!is.finite(base) || base <= 0L) stop("Invalid base in DuckDB params.")
        
        # genotype format
        fmt <- .duckdb_get_param(con, "genotype_format", default = "auto")
        if (identical(fmt, "auto")) {
          fmt <- if (any(grepl("(_1|\\.[0-9]+)$", marker_cols_raw))) "paired" else "string"
        }
        
        .time_it("HF build (chunked)", {
          .duckdb_build_hf_from_raw_chunked(
            con              = con,
            tbl_raw          = rv$tbl_raw,
            tbl_meta         = rv$tbl_meta,
            tbl_hf           = rv$tbl_hf,
            marker_cols      = marker_cols_raw,   # <-- FIX
            genotype_format  = fmt,
            base             = base,
            missing_code_raw = missing_code_raw,
            missing_gt       = 0L,
            batch_size       = as.integer(batch_size),
            attach_pop_code  = TRUE
          )
        })
        
        invisible(TRUE)
      }
      
      tryCatch(
        {
          .build_hf_from_params(rv$con, rv, missing_code_raw = input$missing_code, batch_size = 10000L)
          
          rv$preview_meta <- tryCatch(
            DBI::dbGetQuery(rv$con, sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_meta))),
            error = function(e) NULL
          )
          
          render_map_from_db()
          rv$db_tick <- rv$db_tick + 1L
        },
        error = function(e) {
          shinyalert::shinyalert("HF build failed", conditionMessage(e), type = "error")
        }
      )
    
      
    })
    
    
    # ---------------------------------------------------------#
    # Upload user data ####
    # ---------------------------------------------------------#
    shiny::observeEvent(input$load_user_data, {
      shiny::req(input$file1)
      shiny::req(input$file1$datapath)
      
      uploaded_file(input$file1$datapath[1])
      # OR directly: rv$file_path <- input$file1$datapath[1]
      
      reset_downstream_state(rv)
      
      rv$file_path <- uploaded_file()
      rv$sep <- input$sep
      rv$header <- input$header
      # ---------------------------------------------------------#
      # State 1: import RAW
      # ---------------------------------------------------------#
      
      ok_db <- tryCatch(
        .duckdb_import_raw(
          con       = rv$con,
          tbl_raw   = rv$tbl_raw,
          file_path = rv$file_path,
          sep       = rv$sep,
          header    = rv$header
        ),
        error = function(e) {
          shinyalert::shinyalert("DuckDB import failed", conditionMessage(e), type = "error")
          FALSE
        }
      )
      if (!isTRUE(ok_db)) return()
      # ---------------------------------------------------------#
      # State 2: header + raw preview
      # ---------------------------------------------------------#
      
      rv$colnames_all <- tryCatch(DBI::dbListFields(rv$con, rv$tbl_raw), error = function(e) character(0))
      if (!length(rv$colnames_all)) {
        shinyalert::shinyalert("Header read failed", "No columns detected in DuckDB raw table.", type = "error")
        return()
      }
      
      update_metadata_choices(session, rv$colnames_all)
      
      rv$preview_raw <- tryCatch(
        DBI::dbGetQuery(rv$con, sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_raw))),
        error = function(e) NULL
      )
      if (is.null(rv$preview_raw) || !nrow(rv$preview_raw)) {
        shinyalert::shinyalert("No data rows", "Raw table has 0 rows after import.", type = "warning")
        return()
      }
      # ---------------------------------------------------------#
      # State 3: auto-detect on raw preview
      # ---------------------------------------------------------#
      
      miss_auto <- normalize_missing_code(input$missing_code)
      det <- detect_columns_auto(rv$preview_raw, missing_info = miss_auto)
      rv$det <- det
      
      populate_manual_from_detection(session, det, rv$colnames_all)
      
      # keep meta_ranges for build_meta (same as before)
      meta_idx <- match(det$metadata_cols, rv$colnames_all)
      meta_idx <- meta_idx[!is.na(meta_idx)]
      meta_ranges <- .compress_idx_ranges(meta_idx)
      
      # State 4: if population missing -> stop (manual UI ready)
      if (is.null(det$population) || is.na(det$population) || !nzchar(det$population)) {
        shinyalert::shinyalert(
          "Auto-detect incomplete",
          "Population column not detected. Either change the separator or use the manual mode.",
          type = "warning"
        )
        return()
      }
      
      # ---------------------------------------------------------#
      # State 5: AUTO format
      # ---------------------------------------------------------#
      if (!is.null(det$marker_range) && nzchar(det$marker_range)) {
        
        ok_fmt <- do_assign_and_format(
          colnames_all    = rv$colnames_all,
          pop_data        = det$population,
          col_ranges_data = det$marker_range,
          metadata_ranges = meta_ranges,
          missing_code    = input$missing_code,
          ploidy          = as.numeric(input$ploidy),
          latitude_data   =  det$latitude  %|||% "",
          longitude_data  = det$longitude %|||% "",
          selected_levels = character(0),
          make_map        = FALSE,
          con             = rv$con,
          tbl_raw         = rv$tbl_raw
        )
        
        if (!isTRUE(ok_fmt)) {
          shinyalert::shinyalert("Auto-format failed",
                                 "Auto-detection ran but formatting did not complete. Use the manual panel and click 'Assign metadata'.",
                                 type = "warning")
          return()
        }
        
        ok_fmt_tbl <- tryCatch(DBI::dbExistsTable(rv$con, rv$tbl_formatted_preview), error = function(e) FALSE)
        if (!isTRUE(ok_fmt_tbl)) {
          shinyalert::shinyalert("Auto-format failed",
                                 "Formatting reported success, but the DuckDB formatted table is missing.",
                                 type = "error")
          return()
        }
        
        rv$preview_raw <- tryCatch(
          DBI::dbGetQuery(rv$con, sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_formatted_preview))),
          error = function(e) NULL
        )
      }
      
      # ---- UI preview
      rv$preview_raw <- tryCatch(
        DBI::dbGetQuery(
          rv$con,
          sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_formatted_preview))
        ),
        error = function(e) NULL
      )
      
      if (is.null(rv$preview_raw) || !nrow(rv$preview_raw)) {
        shinyalert::shinyalert(
          "Formatted preview empty",
          "Formatted table exists but returned 0 rows (unexpected).",
          type = "warning"
        )
      }
      
      # ------------------------------#
      # State 6: build META (Population + GPS + metadata only)
      # ------------------------------#
      tryCatch(
        .time_it("META build", {
          .duckdb_build_meta(
            con = rv$con,
            tbl_raw = rv$tbl_raw,
            tbl_meta = rv$tbl_meta,
            colnames_all = rv$colnames_all,
            pop_data = det$population,
            latitude_data = det$latitude,
            longitude_data = det$longitude,
            metadata_ranges = meta_ranges
          )
        })
        ,
        error = function(e) {
          shinyalert::shinyalert("DuckDB meta build failed", conditionMessage(e), type = "error")
          NULL
        }
      )
      
      # ---- ensure meta.individual exists and matches raw.rowid
      DBI::dbExecute(
        rv$con,
        sprintf(
          "
  ALTER TABLE %s ADD COLUMN IF NOT EXISTS individual BIGINT;

  UPDATE %s m
  SET individual = r.rowid
  FROM %s r
  WHERE m.rowid = r.rowid;
",
          .sql_ident(rv$tbl_meta),
          .sql_ident(rv$tbl_meta),
          .sql_ident(rv$tbl_raw)
        )
      )

      # ---- build HF (chunked) AFTER meta exists
      .build_hf_from_params <- function(con, rv, missing_code_raw, batch_size = 10000L) {
        shiny::req(con)
        
        # marker_cols_raw (raw columns, may include suffixes like _1 / .1)
        marker_json <- tryCatch(
          DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='marker_cols_raw'")$value[1],
          error = function(e) NA_character_
        )
        
        marker_cols_raw <- if (!is.na(marker_json) && nzchar(marker_json)) {
          tryCatch(jsonlite::fromJSON(marker_json), error = function(e) character(0))
        } else character(0)
        
        if (!length(marker_cols_raw)) stop("No marker_cols_raw found in DuckDB params.")
        
        # base
        base <- tryCatch(
          DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='base'")$value[1],
          error = function(e) NA_character_
        )
        base <- suppressWarnings(as.integer(base))
        if (!is.finite(base) || base <= 0L) stop("Invalid base in DuckDB params.")
        
        # genotype format
        fmt <- .duckdb_get_param(con, "genotype_format", default = "auto")
        if (identical(fmt, "auto")) {
          fmt <- if (any(grepl("(_1|\\.[0-9]+)$", marker_cols_raw))) "paired" else "string"
        }
        
        .time_it("HF build (chunked)", {
          .duckdb_build_hf_from_raw_chunked(
            con              = con,
            tbl_raw          = rv$tbl_raw,
            tbl_meta         = rv$tbl_meta,
            tbl_hf           = rv$tbl_hf,
            marker_cols      = marker_cols_raw,   # <-- FIX
            genotype_format  = fmt,
            base             = base,
            missing_code_raw = missing_code_raw,
            missing_gt       = 0L,
            batch_size       = as.integer(batch_size),
            attach_pop_code  = TRUE
          )
        })
        
        invisible(TRUE)
      }
      
      tryCatch(
        {
          .build_hf_from_params(rv$con, rv, missing_code_raw = input$missing_code, batch_size = 10000L)
          
          rv$preview_meta <- tryCatch(
            DBI::dbGetQuery(rv$con, sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_meta))),
            error = function(e) NULL
          )
          
          render_map_from_db()
          rv$db_tick <- rv$db_tick + 1L
        },
        error = function(e) {
          shinyalert::shinyalert("HF build failed", conditionMessage(e), type = "error")
        }
      )
    
      
    })
      
    # ---------------------------------------------------------#
    # MANUAL: run_assign  ####
    # ---------------------------------------------------------#
    shiny::observeEvent(input$run_assign, {
      shiny::req(rv$file_path, rv$colnames_all)
      
      ok <- do_assign_and_format(
        colnames_all    = rv$colnames_all,
        pop_data        = input$pop_data,
        col_ranges_data = input$col_ranges_data,
        metadata_ranges = input$metadata_ranges,
        missing_code    = input$missing_code,
        ploidy          = as.numeric(input$ploidy),
        latitude_data   = input$latitude_data,
        longitude_data  = input$longitude_data,
        make_map        = FALSE,
        con             = rv$con,
        tbl_raw         = rv$tbl_raw
      )
      if (!isTRUE(ok)) return()
      
      rv$preview_raw <- DBI::dbGetQuery(
        rv$con,
        sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_formatted_preview))
      )
      
      # ------------------------------#
      # State: build META 
      # ------------------------------#
      meta_ranges <- input$metadata_ranges %||% ""
      
      tryCatch(
        .time_it("META build", {
          .duckdb_build_meta(
            con             = rv$con,
            tbl_raw         = rv$tbl_raw,
            tbl_meta        = rv$tbl_meta,
            colnames_all    = rv$colnames_all,
            pop_data        = input$pop_data,
            latitude_data   = input$latitude_data,
            longitude_data  = input$longitude_data,
            metadata_ranges = meta_ranges
          )
        }),
        error = function(e) {
          shinyalert::shinyalert("DuckDB meta build failed", conditionMessage(e), type = "error")
          NULL
        }
      )
      # ---- ensure meta.individual exists and matches raw.rowid
      DBI::dbExecute(
        rv$con,
        sprintf(
          "
  ALTER TABLE %s ADD COLUMN IF NOT EXISTS individual BIGINT;

  UPDATE %s m
  SET individual = r.rowid
  FROM %s r
  WHERE m.rowid = r.rowid;
",
          .sql_ident(rv$tbl_meta),
          .sql_ident(rv$tbl_meta),
          .sql_ident(rv$tbl_raw)
        )
      )
      
      
      # ---- build HF (chunked) AFTER meta exists
      .build_hf_from_params <- function(con, rv, missing_code_raw, batch_size = 10000L) {
        shiny::req(con)
        
        # marker_cols_raw (raw columns, may include suffixes like _1 / .1)
        marker_json <- tryCatch(
          DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='marker_cols_raw'")$value[1],
          error = function(e) NA_character_
        )
        
        marker_cols_raw <- if (!is.na(marker_json) && nzchar(marker_json)) {
          tryCatch(jsonlite::fromJSON(marker_json), error = function(e) character(0))
        } else character(0)
        
        if (!length(marker_cols_raw)) stop("No marker_cols_raw found in DuckDB params.")
        
        # base
        base <- tryCatch(
          DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='base'")$value[1],
          error = function(e) NA_character_
        )
        base <- suppressWarnings(as.integer(base))
        if (!is.finite(base) || base <= 0L) stop("Invalid base in DuckDB params.")
        
        # genotype format
        fmt <- .duckdb_get_param(con, "genotype_format", default = "auto")
        if (identical(fmt, "auto")) {
          fmt <- if (any(grepl("(_1|\\.[0-9]+)$", marker_cols_raw))) "paired" else "string"
        }
        
        .time_it("HF build (chunked)", {
          .duckdb_build_hf_from_raw_chunked(
            con              = con,
            tbl_raw          = rv$tbl_raw,
            tbl_meta         = rv$tbl_meta,
            tbl_hf           = rv$tbl_hf,
            marker_cols      = marker_cols_raw,   # <-- FIX
            genotype_format  = fmt,
            base             = base,
            missing_code_raw = missing_code_raw,
            missing_gt       = 0L,
            batch_size       = as.integer(batch_size),
            attach_pop_code  = TRUE
          )
        })
        
        invisible(TRUE)
      }
      
      
      tryCatch(
        {
          .build_hf_from_params(rv$con, rv, missing_code_raw = input$missing_code, batch_size = 10000L)
          
          rv$preview_meta <- tryCatch(
            DBI::dbGetQuery(rv$con, sprintf("SELECT * FROM %s LIMIT 100;", .sql_ident(rv$tbl_meta))),
            error = function(e) NULL
          )
          
          render_map_from_db()
          rv$db_tick <- rv$db_tick + 1L
        },
        error = function(e) {
          shinyalert::shinyalert("HF build failed", conditionMessage(e), type = "error")
        }
      )
      
    })
    
    # ---------------------------------------------------------#
    # Ouput formatted data table ####
    # ---------------------------------------------------------#
    output$formatted_table <- DT::renderDT({
      df <- rv$preview_raw %||% rv$preview_meta
      if (is.null(df)) df <- rv$preview_meta
      shiny::req(df)
      
      df <- df[, !names(df) %in% "individual", drop = FALSE]
      col_nums <- seq_along(names(df))
      sketch <- htmltools::withTags(
        table(
          class = "display nowrap",
          thead(
            tr(lapply(names(df), function(nm)
              th(nm, style = "text-align:center; white-space:nowrap; min-width:90px;")
            )),
            tr(lapply(col_nums, function(i)
              th(i, style = "font-weight:normal; color:#888; font-size:11px; text-align:center; min-width:90px;")
            ))
          )
        )
      )

      DT::datatable(
        df,
        container = sketch,
        rownames = FALSE,
        escape = FALSE,
        class = "nowrap",
        options = list(
          pageLength = 10,
          dom = 't<"bottom"lip>',
          scrollX = TRUE,
          autoWidth = FALSE,
          columnDefs = list(list(className = "dt-center", targets = "_all")),
          orderCellsTop = TRUE
        )
      ) %>%
        DT::formatStyle(columns = names(df), fontSize = "14px")
    })
    
    # ---------------------------------------------------------#
    # Output formatted object summary ####
    # ---------------------------------------------------------#
    output$formatted_summary <- shiny::renderPrint({
      
      tick <- rv$db_tick  # force reactivity
      
      cat("Database import summary\n")
      cat("=======================\n\n")
      
      # ---- guards
      if (is.null(rv$con)) {
        cat("No database connection.\n")
        return()
      }
      
      tbls <- tryCatch(DBI::dbListTables(rv$con), error = function(e) character(0))
      cat("Tables in DuckDB:", if (length(tbls)) paste(tbls, collapse = ", ") else "(none)", "\n\n")
      
      # ---- HF status (Design 1)
      cat("HF table:", ifelse(rv$tbl_hf %in% tbls, "YES", "NO"), "\n")
      if (rv$tbl_hf %in% tbls) {
        n_hf <- DBI::dbGetQuery(rv$con, sprintf("SELECT COUNT(*) AS n FROM %s;", .sql_ident(rv$tbl_hf)))$n
        cat("HF rows:", n_hf, "\n")
      }
      cat("\n")
      
      if (is.null(rv$tbl_raw) || !nzchar(rv$tbl_raw) || !(rv$tbl_raw %in% tbls)) {
        cat("Raw table not available yet.\n")
        return()
      }
      
      # quote identifiers safely
      raw_sql  <- as.character(DBI::dbQuoteIdentifier(rv$con, rv$tbl_raw))
      meta_sql <- if (!is.null(rv$tbl_meta) && nzchar(rv$tbl_meta) && rv$tbl_meta %in% tbls)
        as.character(DBI::dbQuoteIdentifier(rv$con, rv$tbl_meta)) else NA_character_
      
      # ---- core counts
      n_raw <- DBI::dbGetQuery(rv$con, paste0("SELECT COUNT(*) AS n FROM ", raw_sql))$n
      cols_raw <- DBI::dbListFields(rv$con, rv$tbl_raw)
      
      cat("Raw table:", rv$tbl_raw, "\n")
      cat("Individuals:", n_raw, "\n")
      cat("Variables:", length(cols_raw), "\n\n")
      
      # DB-first: prefer locus_cols (collapsed loci), fallback to marker_cols_raw (raw cols) then collapse
      marker_names  <- character(0)
      marker_source <- "none"
      
      locus_json <- tryCatch(
        DBI::dbGetQuery(rv$con, "SELECT value FROM params WHERE key='locus_cols'")$value[1],
        error = function(e) NA_character_
      )
      
      if (!is.na(locus_json) && nzchar(locus_json)) {
        marker_names  <- tryCatch(jsonlite::fromJSON(locus_json), error = function(e) character(0))
        marker_source <- "duckdb params (locus_cols)"
      } else {
        raw_json <- tryCatch(
          DBI::dbGetQuery(rv$con, "SELECT value FROM params WHERE key='marker_cols_raw'")$value[1],
          error = function(e) NA_character_
        )
        
        raw_cols <- if (!is.na(raw_json) && nzchar(raw_json)) {
          tryCatch(jsonlite::fromJSON(raw_json), error = function(e) character(0))
        } else character(0)
        
        if (length(raw_cols)) {
          marker_names  <- unique(sub("(_1|\\.[0-9]+)$", "", as.character(raw_cols)))
          marker_source <- "duckdb params (marker_cols_raw \u2192 collapsed)"
        } else if (!is.null(rv$det$locus_cols) && length(rv$det$locus_cols)) {
          # optional fallback if you ever add it to detect_columns_auto()
          marker_names  <- rv$det$locus_cols
          marker_source <- "auto-detect (rv$det$locus_cols)"
        } else if (!is.null(rv$det$marker_cols) && length(rv$det$marker_cols)) {
          # last resort: collapse what auto-detect found
          marker_names  <- unique(sub("(_1|\\.[0-9]+)$", "", as.character(rv$det$marker_cols)))
          marker_source <- "auto-detect (rv$det$marker_cols \u2192 collapsed)"
        }
      }
      
      
      has_markers <- length(marker_names) > 0
      cat("Has markers:", ifelse(has_markers, "TRUE", "FALSE"), "\n")
      cat("Marker set source:", marker_source, "\n")
      cat("Markers (n):", length(marker_names), "\n")
      
      if (length(marker_names)) {
        show_n <- min(10, length(marker_names))
        cat("Marker names (first", show_n, "):\n")
        cat(" - ", paste(marker_names[seq_len(show_n)], collapse = ", "), "\n", sep = "")
        if (length(marker_names) > 10) cat(" - ...\n")
      }
      
      # warn if formatted still contains suffixes like _1 / .1
      if (length(marker_names) && any(grepl("(\\.|_)\\d+$", marker_names))) {
        cat("WARNING: marker_cols contain allele suffixes (e.g. _1 / .1). Collapsing may not have occurred.\n")
      }
      
      cat("\n")
      
      # ---- metadata columns (in addition to Population)
      cat("Meta table created:", ifelse(!is.na(meta_sql), "YES", "NO"), "\n")
      if (!is.na(meta_sql)) {
        cols_meta <- DBI::dbListFields(rv$con, rv$tbl_meta)
        meta_extra <- setdiff(cols_meta, c("Population", "Latitude", "Longitude"))
        
        show_n <- min(10, length(meta_extra))
        cat("Metadata columns (excluding Population) (first", show_n, "):\n")
        if (show_n == 0) {
          cat(" - (none)\n")
        } else {
          cat(" - ", paste(meta_extra[seq_len(show_n)], collapse = ", "), "\n", sep = "")
          if (length(meta_extra) > 10) cat(" - ...\n")
        }
        cat("\n")
      }
      
      # ---- independent locations (from meta, based on unique lat/lon pairs)
      if (!is.na(meta_sql) && all(c("Latitude","Longitude") %in% DBI::dbListFields(rv$con, rv$tbl_meta))) {
        
        loc_df <- DBI::dbGetQuery(
          rv$con,
          paste0(
            "SELECT Latitude, Longitude, COUNT(*) AS n_indiv ",
            "FROM ", meta_sql, " ",
            "WHERE Latitude IS NOT NULL AND Longitude IS NOT NULL ",
            "GROUP BY Latitude, Longitude ",
            "ORDER BY n_indiv DESC"
          )
        )
        
        n_loc <- nrow(loc_df)
        cat("Independent locations (unique Lat/Lon):", n_loc, "\n")
        
        if (n_loc > 0) {
          show_n <- min(10, n_loc)
          cat("Locations (first", show_n, "):\n")
          for (i in seq_len(show_n)) {
            cat(sprintf(" - (%.6f, %.6f) [n=%d]\n",
                        loc_df$Latitude[i], loc_df$Longitude[i], loc_df$n_indiv[i]))
          }
          if (n_loc > 10) cat(" - ...\n")
        }
        
      } else {
        cat("Independent locations: NA (Latitude/Longitude not available in meta)\n")
      }
    })
    
    
  })
}
