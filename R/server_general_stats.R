# server_general_stats.R
# ==============================================================================#
# server_general_stats.R  (DB-first refactor contract)
#
# Source of truth:
#   DuckDB tables only: params, meta, hf  (optionally raw for diagnostics)
#
# Genotypes:
#   Stored in DuckDB hf.gt as packed int: gt = a*base + b ; missing = 0
#   This module MUST NOT decode gt to "a/b" strings in R for computation.
#   (Decoding is allowed only for display/debug, never for stats pipelines.)
#
# Idempotency:
#   Every computation is a pure function of (hf, meta, params, user inputs).
#   No mutation, no hidden state. Re-running yields same result.
#
# Caching:
#   Heavy computations MUST be cached by a stable cache key:
#     key = list(db_tick, tbl_hf, tbl_meta, params_hash, user_inputs_hash)
#   Invalidation happens ONLY when db_tick changes or inputs change.
#
# Outputs:
#   All results tables are standardized:
#     - ID column (Locus or Population label) + numeric columns
#     - "Overall" row last
# ==============================================================================#


## =========================================================#
# Helpers ####
## =========================================================#
hs_by_pop_locus_from_mat <- function(mat, base) {
  stopifnot(is.matrix(mat), ncol(mat) >= 2L, base > 1L)
  
  pop_codes  <- as.integer(mat[, 1])
  pop_levels <- attr(mat, "pop_levels")
  loci       <- colnames(mat)[-1L]
  
  pops <- sort(unique(pop_codes[is.finite(pop_codes) & pop_codes > 0L]))
  out  <- vector("list", length(loci) * length(pops))
  ii   <- 1L
  
  for (j in seq_along(loci)) {
    g_all <- as.integer(mat[, j + 1L])
    
    for (pp in pops) {
      idx <- which(pop_codes == pp)
      g   <- g_all[idx]
      
      ok_gt <- is.finite(g) & g > 0L
      if (!any(ok_gt)) next
      
      g  <- g[ok_gt]
      a1 <- g %/% base
      a2 <- g %%  base
      
      ok <- a1 > 0L & a2 > 0L
      if (!any(ok)) next
      
      a1 <- a1[ok]
      a2 <- a2[ok]
      n  <- length(a1)   # number of diploid genotypes
      
      hs <- if (n <= 1L) {
        NA_real_
      } else {
        cnt <- table(c(a1, a2))
        p   <- as.numeric(cnt) / (2 * n)
        (2 * n / (2 * n - 1)) * (1 - sum(p^2))
      }
      
      pop_name <- if (!is.null(pop_levels) &&
                      pp >= 1L && pp <= length(pop_levels)) {
        as.character(pop_levels[[pp]])
      } else {
        as.character(pp)
      }
      
      out[[ii]] <- data.frame(
        Locus      = loci[j],
        Population = pop_name,
        n_mat      = n,
        Hs_mat     = hs,
        stringsAsFactors = FALSE
      )
      ii <- ii + 1L
    }
  }
  
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0L) {
    return(data.frame(
      Locus = character(),
      Population = character(),
      n_mat = integer(),
      Hs_mat = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  out <- do.call(rbind, out)
  out[order(out$Locus, out$Population), , drop = FALSE]
}

duck_get_base <- function(con) {
  p <-  .duckdb_get_params(con)
  
  # 1) if base already stored, trust it
  base <- suppressWarnings(as.integer(p$base %|||% p$base_scalar_preview %|||% p$base_scalar_full))
  if (length(base) == 1L && is.finite(base) && base > 0L) {
    return(as.integer(base))
  }
  
  # 2) otherwise compute from haplotype length / width
  hl <- suppressWarnings(as.integer(
    p$haplotype_length %|||%
      p$width_scalar_preview %|||%
      p$width_scalar_full %|||%
      p$length_haplotype
  ))
  
  if (length(hl) != 1L || is.na(hl) || hl <= 0L) {
    stop("params must define base (positive integer) or haplotype_length/width_scalar_* to compute base=10^L")
  }
  
  as.integer(10L ^ hl)
}

duck_get_haplotype_length <- function(con, default = 3L) {
  p <- .duckdb_get_params(con)
  
  # 1) preferred explicit key
  hl <- suppressWarnings(as.integer(p[["haplotype_length"]]))
  if (is.finite(hl) && hl >= 1L) return(hl)
  
  # 2) backward-compatible alias
  hl <- suppressWarnings(as.integer(p[["width_scalar_full"]]))
  if (is.finite(hl) && hl >= 1L) return(hl)
  
  # 3) derive from base if needed
  b <- suppressWarnings(as.numeric(p[["base"]]))
  if (is.finite(b) && b > 1) {
    hl <- as.integer(round(log10(b)))
    if (is.finite(hl) && hl >= 1L) return(hl)
  }
  
  stop("Invalid haplotype_length in params")
}


duck_get_loci <- function(con, tbl_hf = "hf") {
  DBI::dbGetQuery(con, sprintf(
    "SELECT DISTINCT locus_id AS Locus FROM %s ORDER BY 1",
    sql_ident(con, tbl_hf)
  ))$Locus
}
gs_dt_options <- function(pageLength = 10L) {
  list(
    dom = paste0(
      "<'row'<'col-sm-6'l><'col-sm-6'f>>",   
      "<'row'<'col-sm-12'B>>",               
      "<'row'<'col-sm-12'tr>>",
      "<'row'<'col-sm-5'i><'col-sm-7'p>>"
    ),
    buttons = list(
      list(
        extend = "copy",
        text = "Copy",
        title = NULL,
        exportOptions = list(columns = ":visible")
      )
    ),
    pageLength = pageLength,
    scrollX = TRUE,
    autoWidth = FALSE
  )
}

## =========================================================#
# Server_general_stats ####
## =========================================================#

server_general_stats <- function(id, rv) {
  
  moduleServer(id, function(input, output, session) {
    
    
    # ------------------------------------------------------------------#
    # Parallel controls (server-side defaults)
    # - If you already have UI inputs for threads/seed, replace these.
    # ------------------------------------------------------------------#
    .n_threads <- reactive({
      if (!is.null(input$n_threads) && is.finite(input$n_threads)) {
        return(max(1L, as.integer(input$n_threads)))
        
      }
      # sensible default: use all available cores (or leave one free)
      nc <- NA_integer_
      if (requireNamespace("parallel", quietly = TRUE)) {
        nc <- suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)))
      }
      if (!is.finite(nc) || length(nc) != 1L || nc < 1L) nc <- 1L
      max(1L, nc - 1L)
    })
    
    .seed <- reactive({
      # prefer an existing input if you have it
      if (!is.null(input$seed) && is.finite(input$seed)) {
        return(as.numeric(input$seed))
      }
      1L
    })
    
    ## =========================================================#
    ## DB truth layer (ONLY sources of truth)
    ## =========================================================#
    
    db_tick <- reactive({ rv$db_tick })
    con_r   <- reactive({ shiny::req(rv$con); rv$con })
    
    tbl_meta_r <- reactive({ rv$tbl_meta %|||% "meta" })
    
    # Prefer params key tbl_hf if you store it; fallback to "hf"
    tbl_hf_r <- reactive({
      con <- con_r()
      if (duck_tbl_exists(con, "params")) {
        p <-  .duckdb_get_params(con)
        th <- p$tbl_hf %|||% "hf"
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
    
    n_pop_db_r <- reactive({
      db_ready()
      con <- con_r()
      tbl <- tbl_meta_r()
      DBI::dbGetQuery(con, sprintf(
        "SELECT COUNT(DISTINCT Population) AS n FROM %s WHERE Population IS NOT NULL",
        sql_ident(con, tbl)
      ))$n[[1]]
    })
    
    params_r <- reactive({
      db_ready()
      .duckdb_get_params(con_r())
    })
    
    base_r <- reactive({
      p <- params_r()
      
      base <- suppressWarnings(as.integer(p$base %|||% p$base_scalar_full %|||% p$base_scalar_preview))
      if (length(base) == 1L && is.finite(base) && base > 1L) return(as.integer(base))
      
      hl <- suppressWarnings(as.integer(p$haplotype_length %|||% p$length_haplotype %|||% p$width_scalar_full))
      shiny::validate(need(length(hl) == 1L && is.finite(hl) && hl > 0L,
                    "params must define base OR haplotype_length/length_haplotype to compute base = 10^L"))
      as.integer(10L ^ hl)
    })
    
    meta_r <- reactive({
      db_ready()
      con <- con_r()
      
      DBI::dbGetQuery(con, sprintf("
    SELECT individual, Population
    FROM %s
    WHERE Population IS NOT NULL
    ORDER BY individual
  ", sql_ident(con, tbl_meta_r())))
    })
    
    hf_mat_r <- reactive({
      db_ready()
      con  <- con_r()
      meta <- meta_r()
      tbl_hf <- tbl_hf_r()
      
      shiny::validate(need(nrow(meta) > 0, "meta table is empty"))
      
      pop_levels <- sort(unique(meta$Population))
      pop_code   <- match(meta$Population, pop_levels)
      
      loci <- DBI::dbGetQuery(con, sprintf(
        "SELECT DISTINCT locus_id FROM %s ORDER BY 1",
        sql_ident(con, tbl_hf)
      ))$locus_id
      shiny::validate(need(length(loci) > 0, "hf table has no loci"))
      
      hf <- DBI::dbGetQuery(con, sprintf("
    SELECT indiv_id AS individual, locus_id, gt AS g
    FROM %s
  ", sql_ident(con, tbl_hf)))
      
      hf$ind_idx   <- match(hf$individual, meta$individual)
      hf$locus_idx <- match(hf$locus_id, loci)
      
      
      N <- nrow(meta); L <- length(loci)
      mat <- matrix(0L, nrow = N, ncol = L + 1L)
      mat[, 1] <- as.integer(pop_code)
      
      ok <- !is.na(hf$ind_idx) & !is.na(hf$locus_idx)
      mat[cbind(hf$ind_idx[ok], hf$locus_idx[ok] + 1L)] <- as.integer(hf$g[ok])
      
      colnames(mat) <- c("pop", as.character(loci))
      attr(mat, "pop_levels") <- pop_levels
      mat
    })
    
    
    ## =========================================================#
    ## Containers
    ## =========================================================#
    
    # Basic stats
    result_stats_reactive <- reactiveVal(NULL)
    plot_output <- reactiveValues()
    result_stats_download <- reactiveVal(NULL)
    result_stats_numeric_reactive <- reactiveVal(NULL)
    result_stats_display_rv <- reactiveVal(NULL)
    
    # Bootstrap FIS Analysis
    fis_boot_results <- reactiveVal(NULL)
    fis_boot_timing <- reactiveVal(NULL)
    perm_results <- reactiveVal(NULL)
    fis_allele_results <- reactiveVal(NULL)
    
    # For "Overall" stats selection
    result_stats_select_reactive <- reactiveVal(NULL)
    
    ## =========================================================#
    ## Convenience reactives for metadata ####
    ## =========================================================#
    
    n_pop <- reactive({
      db_ready()
      m <- meta_r()
      length(unique(m$Population))
    })
    
    pops <- reactive({
      db_ready()
      sort(unique(meta_r()$Population))
    })
    
    n_indv <- reactive({
      db_ready()
      nrow(meta_r())
    })
    
    nloci <- reactive({
      db_ready()
      ncol(hf_mat_r()) - 1L
    })
    
    loci_names <- reactive({
      db_ready()
      colnames(hf_mat_r())[-1L]
    })
    
    hs_by_pop_wide_r <- reactive({
      db_ready()
      con  <- con_r()
      base <- base_r()
      
      long <- duck_hs_by_pop_locus_long(
        con       = con,
        tbl_hf    = tbl_hf_r(),
        tbl_meta  = tbl_meta_r(),
        base      = base,
        missing_code = 0L
      )
      
      shiny::validate(need(nrow(long) > 0, "No Hs results available (check hf/meta tables)."))
      
      # Wide matrix: rows = Locus, cols = Population
      wide <- tidyr::pivot_wider(
        long,
        names_from  = Population,
        values_from = Hs
      )
      
      wide <- as.data.frame(wide, stringsAsFactors = FALSE)
      wide
    })

    ## =========================================================#
    ## Basic stats cache (DB-first)
    ## =========================================================#
    
    basic_cache <- reactiveVal(list(key = NULL, value = NULL))
    
    .get_basic_stats_cached <- function() {
      db_ready()
      
      mat  <- hf_mat_r()
      base <- base_r()
      k    <- n_pop_db_r()
      ui_hash <- .hash_key(list(
        missing_code = 0L
      ))
      key <- .hash_key(list(
        db_tick = db_tick(),
        tbl_hf  = tbl_hf_r(),
        tbl_meta= tbl_meta_r(),
        params_hash = .hash_key(params_r()),
        user_inputs_hash = ui_hash
      ))
      
      
      cur <- basic_cache()
      if (isTRUE(identical(cur$key, key)) && !is.null(cur$value)) {
        return(cur$value)
      }
      
      res <- .compute_basic_stats(mat = mat, base = base, k = k)
      basic_cache(list(key = key, value = res))
      res
    }
    
    ## =========================================================#
    ## Observer Basic stats ####
    ## =========================================================#
    
    observeEvent(input$run_basic_stats, {
      db_ready()
      
      result_stats <- .get_basic_stats_cached()
      
      # keep if any other outputs depend on the full table
      result_stats_reactive(result_stats)
      
      selected_stats <- c(
        "Ho"                 = isTRUE(input$ho_checkbox),
        "Hs"                 = isTRUE(input$hs_checkbox),
        "Ht"                 = isTRUE(input$ht_checkbox),
        "Fit (W&C)"          = isTRUE(input$fit_wc_checkbox),
        "Fis (W&C)"          = isTRUE(input$fis_wc_checkbox),
        "Fst (W&C)"          = isTRUE(input$fst_wc_checkbox),
        "Fst-max (Meirmans)"  = isTRUE(input$fst_max_checkbox),
        "Fst' (Meirmans)"     = isTRUE(input$fst_prim_checkbox),
        # "Fst' (Hedrick)"      = isTRUE(input$fst_prim_hedrick_checkbox),
        "GST"                = isTRUE(input$GST_checkbox),
        "GST''"              = isTRUE(input$GST_sec_checkbox)
      )
      
      # Only keep columns that exist in result_stats
      keep <- names(selected_stats)[selected_stats]
      
      keep <- intersect(keep, colnames(result_stats))
      
      if (length(keep) > 0) {
        
        result_stats_select <- result_stats[, c("ID", keep), drop = FALSE]
        
        # guarantee Overall last
        if ("Overall" %in% result_stats_select$ID) {
          result_stats_select <- rbind(
            result_stats_select[result_stats_select$ID != "Overall", , drop = FALSE],
            result_stats_select[result_stats_select$ID == "Overall", , drop = FALSE]
          )
        }
        
        result_stats_numeric_reactive(result_stats_select)
        
        result_stats_display <- format_numeric_cols(
          result_stats_select,
          digits = 5,
          exclude = "ID"
        )
        result_stats_download(result_stats_display)
        result_stats_select_reactive(result_stats_select)
        result_stats_display_rv(result_stats_display)
        
        
      } else {
        showNotification("No valid statistics to display", type = "warning")
        result_stats_download(NULL)
        result_stats_numeric_reactive(NULL)
      }
    })
    
    
    # ---- Basic stats table output ----
    output$basic_stats_table <- DT::renderDT({
      df <- shiny::req(result_stats_display_rv())
      DT::datatable(
        df,
        extensions = "Buttons",
        options = gs_dt_options(pageLength = 10L),
        rownames = FALSE,
        class = "compact nowrap",
        callback = DT::JS("
      table.columns.adjust();
    ")
      )
    })
    output$download_basic_stats <- downloadHandler(
      filename = function() paste0("basic_statistics_", Sys.Date(), ".csv"),
      content = function(file) {
        df <- shiny::req(result_stats_select_reactive())
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    output$download_basic_stats_txt <- downloadHandler(
      filename = function() paste0("basic_statistics_", Sys.Date(), ".txt"),
      content = function(file) {
        df <- shiny::req(result_stats_select_reactive())
        utils::write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    # =========================================================#
    ### Population-specific stats (DB-native, vectorised) ####
    # =========================================================#
    
    # ---- UI selector choices (no loops, DB query only)
    observe({
      db_ready()
      con <- con_r()
      
      df <- DBI::dbGetQuery(con, sprintf("
    SELECT DISTINCT Population
    FROM %s
    WHERE Population IS NOT NULL
    ORDER BY Population
  ", sql_ident(con, tbl_meta_r())))
      
      shiny::validate(need(nrow(df) > 0, "No populations available in meta table yet."))
      choices <- df$Population
      
      updateSelectInput(session, "selected_pop_overall",
                        choices = choices,
                        selected = choices[1])
    })
    
    # ---- Table: per-locus stats for selected population
    output$basic_stats_by_pop_selected <- DT::renderDT({
      db_ready()
      con  <- con_r()
      base <- base_r()
      
      shiny::req(input$selected_pop_overall)
      pop_name <- input$selected_pop_overall
      
      df <- duck_pop_stats_by_pop_one(
        con       = con,
        pop_name  = pop_name,
        tbl_hf    = tbl_hf_r(),
        tbl_meta  = tbl_meta_r(),
        base      = base,
        missing_code = 0L
      )
      
      if (is.null(df) || nrow(df) == 0) {
        return(
          DT::datatable(
            data.frame(Message = paste("No data available for population:", pop_name)),
            extensions = "Buttons",
            options = gs_dt_options(pageLength = 10L),
            rownames = FALSE,
            class = "compact nowrap",
            callback = DT::JS("table.columns.adjust();")
          )
        )
      }
      
      df_display <- df
      num_cols <- vapply(df_display, is.numeric, logical(1))
      df_display[num_cols] <- lapply(df_display[num_cols], round, 5)
      
      DT::datatable(
        df_display,
        extensions = "Buttons",
        options = gs_dt_options(pageLength = 10L),
        rownames = FALSE,
        class = "compact nowrap",
        callback = DT::JS("table.columns.adjust();")
      )
    })
    
    # ---- Table: overall-by-pop (Ho/Hs), DB-native
    output$overall_by_pop <- renderTable({
      db_ready()
      con  <- con_r()
      base <- base_r()
      
      df <- duck_pop_stats_overall(
        con       = con,
        tbl_hf    = tbl_hf_r(),
        tbl_meta  = tbl_meta_r(),
        base      = base,
        missing_code = 0L
      )
      
      if (is.null(df) || nrow(df) == 0) {
        return(data.frame(Message = "No population data"))
      }
      
      # format for display only
      df$Ho <- round(df$Ho, 4)
      df$Hs <- round(df$Hs, 4)
      df$`Fis (WC)` <- round(df$`Fis (WC)`, 4)
      
      df
    }, rownames = FALSE, digits = 4)
    
    # =========================================================#
    ## Download handlers (population section) ####
    # =========================================================#
    
    # 1) Population-specific (per locus) stats
    output$download_pop_stats <- downloadHandler(
      filename = function() paste0("pop_stats_", input$selected_pop_overall, "_", Sys.Date(), ".csv"),
      content = function(file) {
        db_ready()
        con  <- con_r()
        base <- base_r()
        shiny::req(input$selected_pop_overall)
        
        df <- duck_pop_stats_by_pop_one(
          con       = con,
          pop_name  = input$selected_pop_overall,
          tbl_hf    = tbl_hf_r(),
          tbl_meta  = tbl_meta_r(),
          base      = base,
          missing_code = 0L
        )
        
        if (is.null(df) || nrow(df) == 0) df <- data.frame(Message = "No data available")
        write.csv(df, file, row.names = FALSE)
      }
    )
    
    output$download_pop_stats_txt <- downloadHandler(
      filename = function() paste0("pop_stats_", input$selected_pop_overall, "_", Sys.Date(), ".txt"),
      content = function(file) {
        db_ready()
        con  <- con_r()
        base <- base_r()
        shiny::req(input$selected_pop_overall)
        
        df <- duck_pop_stats_by_pop_one(
          con       = con,
          pop_name  = input$selected_pop_overall,
          tbl_hf    = tbl_hf_r(),
          tbl_meta  = tbl_meta_r(),
          base      = base,
          missing_code = 0L
        )
        
        if (is.null(df) || nrow(df) == 0) df <- data.frame(Message = "No data available")
        write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    
    # 2) Overall by population (Ho/Hs/Fis Nei)
    output$download_overall_by_pop <- downloadHandler(
      filename = function() paste0("overall_by_population_", Sys.Date(), ".csv"),
      content = function(file) {
        db_ready()
        con  <- con_r()
        base <- base_r()
        
        df <- duck_pop_stats_overall(
          con       = con,
          tbl_hf    = tbl_hf_r(),
          tbl_meta  = tbl_meta_r(),
          base      = base,
          missing_code = 0L
        )
        
        if (is.null(df) || nrow(df) == 0) df <- data.frame(Message = "No population data")
        write.csv(df, file, row.names = FALSE)
      }
    )
    
    output$download_overall_by_pop_txt <- downloadHandler(
      filename = function() paste0("overall_by_population_", Sys.Date(), ".txt"),
      content = function(file) {
        db_ready()
        con  <- con_r()
        base <- base_r()
        
        df <- duck_pop_stats_overall(
          con       = con,
          tbl_hf    = tbl_hf_r(),
          tbl_meta  = tbl_meta_r(),
          base      = base,
          missing_code = 0L
        )
        
        if (is.null(df) || nrow(df) == 0) df <- data.frame(Message = "No population data")
        write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    
    output$gene_diversity_table <- DT::renderDT({
      df <- shiny::req(hs_by_pop_wide_r())
      
      df_disp <- df
      num_cols <- names(df_disp)[vapply(df_disp, is.numeric, logical(1))]
      df_disp[num_cols] <- lapply(df_disp[num_cols], round, 3)
      
      rng <- range(unlist(df[num_cols], use.names = FALSE), na.rm = TRUE)
      if (!all(is.finite(rng)) || diff(rng) == 0) rng <- c(0, 1)
      
      DT::datatable(
        df_disp,
        extensions = "Buttons",
        options = gs_dt_options(pageLength = 10L),
        rownames = FALSE,
        class = "compact nowrap",
        callback = DT::JS("table.columns.adjust();")
      ) %>%
        DT::formatStyle(
          columns = num_cols,
          backgroundColor = DT::styleColorBar(rng, "lightblue"),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center",
          color = "black"
        )
    })
    
    output$download_gene_diversity <- downloadHandler(
      filename = function() paste0("gene_diversity_hs_by_pop_", Sys.Date(), ".csv"),
      content = function(file) {
        df <- shiny::req(hs_by_pop_wide_r())
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    output$download_gene_diversity_txt <- downloadHandler(
      filename = function() paste0("gene_diversity_hs_by_pop_", Sys.Date(), ".txt"),
      content = function(file) {
        df <- shiny::req(hs_by_pop_wide_r())
        utils::write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    
    
    
    
    
    
    # ==================================== FIS SECTION ANALYSIS ===============================================
    fis_context <- reactive({
      level <- input$analysis_level
      if (is.null(level) || is.na(level) || level == "") level <- "By Locus"
      
      is_pop <- identical(level, "By Population")
      
      list(
        level  = level,
        is_pop = is_pop,
        # display label only (NOT a column name in df anymore)
        id_label = if (is_pop) "Population" else "Locus",
        title = if (is_pop) {
          "Bootstrap-based FIS inference across populations"
        } else {
          "Bootstrap-based FIS inference across genetic loci"
        }
      )
    })
    
    run_fis_by_pop <- function(
    n_perm,
    n_boot,
    conf_level,
    missing_code = 0L
    ) {
      mat  <- hf_mat_r()
      base <- base_r()
      con  <- con_r()
      
      pop_df <- DBI::dbGetQuery(con, sprintf("
    SELECT DISTINCT Population
    FROM %s
    WHERE Population IS NOT NULL
    ORDER BY Population
  ", sql_ident(con, tbl_meta_r())))
      
      pop_codes <- as.character(sort(unique(mat[, 1])))
      pop_names <- as.character(pop_df$Population[seq_along(pop_codes)])
      
      pop_lookup <- stats::setNames(pop_names, pop_codes)
      
      # 2) Observed population-wise FIS — WC84 ratio-of-sums (C++)
      obs_pop <- wc_fis_by_pop_wc84(
        dat     = mat,
        pop_col = 0L,
        base    = base
      )
      obs_codes <- names(obs_pop)
      if (is.null(obs_codes) || !length(obs_codes)) {
        obs_codes <- pop_codes
        names(obs_pop) <- obs_codes
      } else {
        obs_codes <- as.character(obs_codes)
      }
      
      names(obs_pop) <- unname(pop_lookup[obs_codes])
      
      # 3) Bootstrap individuals within populations (C++)
      boot_mat <- boot_indiv_wc_fis_by_pop(
        mat     = mat,
        pop_col = 0L,
        NAcode  = as.integer(missing_code),
        B       = as.integer(n_boot),
        base    = base
      )
      # rename bootstrap columns to pop_names
      if (!is.null(colnames(boot_mat))) {
        boot_codes <- as.character(colnames(boot_mat))
        colnames(boot_mat) <- unname(pop_lookup[boot_codes])
      }
      
      alpha <- (1 - conf_level) / 2
      boot_mean <- colMeans(boot_mat, na.rm = TRUE)
      ci_l <- apply(boot_mat, 2, stats::quantile, probs = alpha,     na.rm = TRUE, type = 7)
      ci_u <- apply(boot_mat, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE, type = 7)
      
      # 4) Permutation test (C++) + p-values (two-sided abs)
      perm_res <- NULL
      pvals <- rep(NA_real_, length(pop_names))
      names(pvals) <- pop_names
      
      if (!is.null(n_perm) && n_perm > 0) {
        perm_res <- batch_permute_wc_fis_by_pop(
          dat            = mat,
          pop_col_1based = 1L,
          base           = as.integer(base),
          B              = as.integer(n_perm)
        )
        
        # rename permutation columns to pop_names
        if (!is.null(colnames(perm_res))) {
          perm_codes <- as.character(colnames(perm_res))
          colnames(perm_res) <- unname(pop_lookup[perm_codes])
        }
        
        common_names <- intersect(names(obs_pop), colnames(perm_res))
        pvals <- rep(NA_real_, length(names(obs_pop)))
        names(pvals) <- names(obs_pop)
        
        pvals[common_names] <- vapply(common_names, function(pn) {
          permj <- perm_res[, pn]
          permj <- permj[is.finite(permj)]
          obsj  <- obs_pop[pn]
          if (!is.finite(obsj) || length(permj) == 0) return(NA_real_)
          ge <- sum(abs(permj) >= abs(obsj))
          (ge + 1) / (length(permj) + 1)
        }, numeric(1))
        
      }
      
      # 5) Overall row (mean across populations)
      obs_overall <- mean(obs_pop, na.rm = TRUE)
      
      boot_overall <- rowMeans(boot_mat, na.rm = TRUE)
      boot_overall <- boot_overall[is.finite(boot_overall)]
      
      if (length(boot_overall) > 0) {
        overall_boot_mean <- mean(boot_overall, na.rm = TRUE)
        overall_ci_l <- as.numeric(stats::quantile(boot_overall, probs = alpha,     na.rm = TRUE, type = 7))
        overall_ci_u <- as.numeric(stats::quantile(boot_overall, probs = 1 - alpha, na.rm = TRUE, type = 7))
      } else {
        overall_boot_mean <- NA_real_
        overall_ci_l <- NA_real_
        overall_ci_u <- NA_real_
      }
      
      overall_p <- NA_real_
      if (!is.null(perm_res)) {
        perm_overall <- rowMeans(perm_res, na.rm = TRUE)
        perm_overall <- perm_overall[is.finite(perm_overall)]
        if (length(perm_overall) > 0 && is.finite(obs_overall)) {
          ge <- sum(abs(perm_overall) >= abs(obs_overall))
          overall_p <- (ge + 1) / (length(perm_overall) + 1)
        }
      }
      
      # 6) Final table
      final_df <- data.frame(
        ID           = pop_names,
        Observed_FIS = as.numeric(obs_pop[pop_names]),
        Boot_Mean    = as.numeric(boot_mean[pop_names]),
        CI_L         = as.numeric(ci_l[pop_names]),
        CI_U         = as.numeric(ci_u[pop_names]),
        P_value      = as.numeric(pvals[pop_names]),
        stringsAsFactors = FALSE
      )
      
      overall_row <- data.frame(
        ID           = "Overall",
        Observed_FIS = obs_overall,
        Boot_Mean    = overall_boot_mean,
        CI_L         = overall_ci_l,
        CI_U         = overall_ci_u,
        P_value      = overall_p,
        stringsAsFactors = FALSE
      )
      
      final_df <- rbind(final_df, overall_row)
      
      list(
        final_table          = final_df,
        permutation_results  = perm_res,
        bootstrap_results    = boot_mat,
        metadata = list(
          id_col     = "Population",
          n_perm     = n_perm,
          n_boot     = n_boot,
          conf_level = conf_level,
          base       = base,
          pop_names  = pop_names
        )
      )
    }
    
    
    run_fis_by_locus <- function(
    n_perm = 1000,
    n_boot = 1000,
    conf_level = 0.95,
    missing_code = 0L
    ) {
      # 1) Source from DuckDB (Design 1)
      mat  <- hf_mat_r()
      base <- base_r()

      # Safety: integer matrix + pop codes in col 1 (R is 1-based)
      mat <- as.matrix(mat)
      storage.mode(mat) <- "integer"
      stopifnot(is.integer(mat))
      stopifnot(ncol(mat) >= 2L)
      stopifnot(all(mat[, 1] > 0, na.rm = TRUE))
      stopifnot(!is.na(base), base > 0L)
      
      # Locus names: all columns except population code column
      locus_names <- colnames(mat)[-1]
      if (is.null(locus_names) || !length(locus_names)) {
        locus_names <- paste0("L", seq_len(ncol(mat) - 1L))
      }
      
      # 2) Observed FIS (C++)
      observed_fis <- fis_wc_cpp(mat, base = as.integer(base))$FIS
      # 2b) Correct overall FIS - WC84 ratio-of-sums: sum(B) / sum(B+C) across loci
      # This is NOT the mean of per-locus FIS values.
      fis_overall_obs_stats <- observed_wc84_stats_cpp(
        dat            = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base)
      )
      fis_obs_overall <- as.numeric(fis_overall_obs_stats$FIS_overall_ratio_of_sums)
      # 3) Bootstrap individuals within populations (C++)
      boot_mat <- boot_indiv_wc_fis(
        mat     = mat,
        pop_col = 0L,
        NAcode  = as.integer(missing_code),
        B       = as.integer(n_boot),
        base    = as.integer(base)
      )
      # 3b) Bootstrap POPULATIONS (pop-block) - captures uncertainty from the
      # sampling of populations, not just individuals within populations.
      # This is the same resampling scheme used for FST.
      boot_pop_mat <- boot_popblock_wc_fis(
        mat     = mat,
        pop_col = 0L,
        NAcode  = as.integer(missing_code),
        B       = as.integer(n_boot),
        base    = as.integer(base)
      )
      # 4) Bootstrap summaries
      summary_list <- summarize_fis_results(
        boot = boot_mat,
        conf = conf_level
      )
      # 4b) Pop-block bootstrap summaries
      summary_pop_list <- summarize_fis_results(
        boot = boot_pop_mat,
        conf = conf_level
      )
      # 5) Build final dataframe (Observed + CI from individual bootstrap, matching Genetix)
      final_df <- create_results_dataframe(
        obs         = observed_fis,
        sum         = summary_list,
        locus_names = locus_names
      )

      # ---- Canonicalise identifier column to ID
      id_src <- NULL
      if ("Locus" %in% names(final_df)) {
        id_src <- "Locus"
      } else if ("ID" %in% names(final_df)) {
        id_src <- "ID"
      } else if (length(locus_names) == nrow(final_df) - 1L) {
        final_df <- tibble::rownames_to_column(final_df, var = "Locus")
        id_src <- "Locus"
      } else {
        stop(
          "run_fis_by_locus(): final_df has no identifier column. Columns are: ",
          paste(names(final_df), collapse = ", ")
        )
      }
      final_df <- dplyr::rename(final_df, ID = dplyr::all_of(id_src))
      # Add Boot_Mean from individual bootstrap.
      final_df$Boot_Mean <- c(
        as.numeric(summary_list$mean),
        as.numeric(summary_list$overall_mean)
      )
      # Override the Overall row with the correct ratio-of-sums observed and CI.
      overall_row_idx <- which(final_df$ID == "Overall")
      if (length(overall_row_idx) == 1L) {
        final_df$Observed_FIS[overall_row_idx] <- fis_obs_overall
        final_df$CI_L[overall_row_idx]         <- summary_list$overall_ci_lower
        final_df$CI_U[overall_row_idx]         <- summary_list$overall_ci_upper
        final_df$Boot_Mean[overall_row_idx]    <- summary_list$overall_mean
      }
      # 6) Permutation test p-values (two-sided abs)
      perm_res <- NULL
      
      if (!is.null(n_perm) && n_perm > 0) {
        
        perm_res <- batch_permute_wc_fis(
          dat            = mat,
          pop_col_1based = 1L,
          base           = as.integer(base),
          B              = as.integer(n_perm)
        )
        
        B_eff <- nrow(perm_res)
        perm_res_loci <- perm_res[, seq_len(ncol(perm_res) - 1L), drop = FALSE]

        pvals_loci <- vapply(seq_along(observed_fis), function(j) {
          permj <- perm_res_loci[, j]
          permj <- permj[is.finite(permj)]
          obsj  <- observed_fis[j]
          if (!is.finite(obsj) || length(permj) == 0) return(NA_real_)
          ge <- sum(abs(permj) >= abs(obsj))
          (ge + 1) / (length(permj) + 1)
        }, numeric(1))

        # Last column of perm_res is the ratio-of-sums overall FIS per replicate
        perm_overall <- perm_res[, ncol(perm_res)]
        perm_overall <- perm_overall[is.finite(perm_overall)]
        
        p_overall <- if (is.finite(fis_obs_overall) && length(perm_overall) > 0) {
          ge <- sum(abs(perm_overall) >= abs(fis_obs_overall))
          (ge + 1) / (length(perm_overall) + 1)
        } else {
          NA_real_
        }
        
        final_df$P_value <- c(pvals_loci, p_overall)
        
        
      } else {
        final_df$P_value <- rep(NA_real_, nrow(final_df))
      }
      
      list(
        final_table             = final_df,
        permutation_results     = perm_res,
        bootstrap_results       = boot_mat,        # individual bootstrap (within pops)
        bootstrap_pop_results   = boot_pop_mat,    # pop-block bootstrap
        ci_indiv = list(          # CIs from individual bootstrap
          mean  = summary_list$mean,
          ci_lo = summary_list$ci_lower,
          ci_hi = summary_list$ci_upper,
          overall_mean  = summary_list$overall_mean,
          overall_ci_lo = summary_list$overall_ci_lower,
          overall_ci_hi = summary_list$overall_ci_upper
        ),
        ci_pop = list(            # CIs from pop-block bootstrap
          mean  = summary_pop_list$mean,
          ci_lo = summary_pop_list$ci_lower,
          ci_hi = summary_pop_list$ci_upper,
          overall_mean  = summary_pop_list$overall_mean,
          overall_ci_lo = summary_pop_list$overall_ci_lower,
          overall_ci_hi = summary_pop_list$overall_ci_upper
        ),
        metadata = list(
          id_col      = "Locus",
          n_perm      = n_perm,
          n_boot      = n_boot,
          conf_level  = conf_level,
          base        = as.integer(base),
          locus_names = locus_names
        )
      )
    }
    
    
    observeEvent(input$Run_FIS_Analysis, {
      db_ready()
      
      if (input$n_perm < 10 || input$n_boot < 10) {
        showNotification(
          "Number of permutations and bootstrap replicates should be at least 10",
          type = "warning"
        )
      }
      
      waiter <- Waiter$new(
        id    = session$ns("fis_results_table"),
        html  = spin_3(),
        color = transparent(0.7)
      )
      waiter$show()
      on.exit(waiter$hide(), add = TRUE)
      
      tryCatch({
        start_time <- Sys.time()
        shinyWidgets::updateProgressBar(session, "fis_progress", value = 5)
        
        level <- if (is.null(input$analysis_level)) "By Locus" else input$analysis_level
        
        results <- if (identical(level, "By Population")) {
          run_fis_by_pop(
            n_perm       = input$n_perm,
            n_boot       = input$n_boot,
            conf_level   = input$conf_level,
            missing_code = 0L
          )
        } else {
          run_fis_by_locus(
            n_perm       = input$n_perm,
            n_boot       = input$n_boot,
            conf_level   = input$conf_level,
            missing_code = 0L
          )
        }
        
        shinyWidgets::updateProgressBar(session, "fis_progress", value = 100)
        
        fis_boot_timing(round(difftime(Sys.time(), start_time, units = "secs"), 1))
        fis_boot_results(results)
        perm_results(results$permutation_results)

        # Per-allele F-stats (FIS, FST, FIT via WC84 components)
        allele_mat  <- as.matrix(hf_mat_r())
        storage.mode(allele_mat) <- "integer"
        allele_base <- as.integer(base_r())
        fis_allele_results(wc84_per_allele_fstats_cpp(
          dat          = allele_mat,
          pop_col      = 0L,
          base         = allele_base,
          missing_code = 0L
        ))

        showNotification("Bootstrap FIS analysis completed successfully!", type = "message")

      }, error = function(e) {
        fis_boot_results(NULL)
        perm_results(NULL)
        fis_allele_results(NULL)
        showNotification(paste("Error in bootstrap analysis:", e$message), type = "error")
      })
    })
    
    
    ## FIS value boxes ----
    ### Global FIS ----
    output$global_fis_box <- renderValueBox({
      shiny::req(fis_boot_results())
      df <- fis_boot_results()$final_table
      shiny::validate(shiny::need("ID" %in% names(df), "FIS results malformed: missing ID column."))
      
      fs <- df %>%
        dplyr::filter(ID == "Overall") %>%
        dplyr::pull(Observed_FIS)
      
      fs <- fs[1]
      display <- ifelse(is.na(fs), "NA",
                        ifelse(abs(fs) < 0.0001, "\u2248 0.0000", format(round(fs, 4), nsmall = 4)))
      
      color <- ifelse(is.na(fs), "light-blue",
                      ifelse(fs > 0.1, "maroon", ifelse(fs > 0.05, "orange", "aqua")))
      valueBox(
        value = display,
        subtitle = HTML("<small>FIS<br>global</small>"),
        color = color,
        icon = icon("dna")
      )
    })
    
    ### Global p-value ----
    output$global_pvalue_box <- renderValueBox({
      shiny::req(fis_boot_results())
      df <- fis_boot_results()$final_table
      shiny::validate(shiny::need("ID" %in% names(df), "FIS results malformed: missing ID column."))
      
      p <- df %>%
        dplyr::filter(ID == "Overall") %>%
        dplyr::pull(P_value)
      
      p <- p[1]
      display <- ifelse(is.na(p), "NA",
                        ifelse(p < 0.0001, "< 0.0001",
                               ifelse(p < 0.001, "< 0.001", format(round(p, 4), nsmall = 4))))
      
      color <- ifelse(is.na(p), "light-blue",
                      ifelse(p < 0.001, "maroon", ifelse(p < 0.05, "orange", "aqua")))
      
      valueBox(
        value = display,
        subtitle = HTML("<small>Global <i>p</i>-value<br>Bilateral test</small>"),
        color = color,
        icon = icon("balance-scale")
      )
    })
    
    
    ### Significant loci ----
    output$significant_loci_box <- renderValueBox({
      shiny::req(fis_boot_results())
      ctx <- fis_context()
      
      df <- fis_boot_results()$final_table
      shiny::validate(shiny::need(all(c("ID", "P_value") %in% names(df)),
                    "FIS results malformed: missing ID or P_value."))
      
      df2 <- df %>% dplyr::filter(ID != "Overall")
      
      sig   <- df2 %>% dplyr::filter(!is.na(P_value), P_value < 0.05) %>% nrow()
      total <- nrow(df2)
      pct   <- ifelse(total > 0, round(100 * sig / total, 1), 0)
      
      color <- ifelse(sig > 0, "yellow", "aqua")
      
      valueBox(
        value = paste0(sig, " / ", total),
        subtitle = HTML(paste0(
          "<small>Significant ", ctx$id_label, "s<br><small>", pct, "% of total</small></small>"
        )),
        color = color,
        icon = icon("vial")
      )
    })
    
    ### Computation time ----
    output$analysis_time_box <- renderValueBox({
      shiny::req(fis_boot_timing())
      sec <- fis_boot_timing()
      display <- ifelse(sec < 60, paste0(sec, " s"),
                        paste0(round(sec / 60, 1), " min"))
      
      valueBox(
        value = display,
        subtitle = HTML("<small>Computation Time<br>Permutation + Bootstrap</small>"),
        color = "light-blue",
        icon = icon("hourglass-half")
      )
    })
    
    ## FIS result table ----
    output$fis_results_table <- DT::renderDT({
      shiny::req(fis_boot_results())
      ctx <- fis_context()
      
      df <- fis_boot_results()$final_table
      shiny::validate(shiny::need("ID" %in% names(df), "FIS results malformed: missing ID column."))
      
      # consistent column order
      df <- df %>%
        dplyr::select(dplyr::any_of(c("ID", "Observed_FIS", "Boot_Mean", "P_value", "CI_L", "CI_U")))
      
      pretty_names <- c(
        ID           = ctx$id_label,
        Observed_FIS = "Observed FIS",
        Boot_Mean    = "Bootstrap mean",
        P_value      = "P-value",
        CI_L         = "CI lower",
        CI_U         = "CI upper"
      )
      
      DT::datatable(
        df,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = c("copy"),
          pageLength = 10,
          scrollX = TRUE
        ),
        rownames = FALSE,
        colnames = unname(pretty_names[names(df)])
      ) %>%
        DT::formatRound(
          columns = intersect(c("Observed_FIS", "Boot_Mean", "P_value", "CI_L", "CI_U"), names(df)),
          digits = 4
        )
    })
    
    
    ### Download FIS table (TXT)
    output$download_fis_table_txt <- downloadHandler(
      filename = function() paste("fis_results_", Sys.Date(), ".txt", sep = ""),
      content  = function(file) {
        shiny::req(fis_boot_results())
        write.table(fis_boot_results()$final_table, file,
                    sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    ### Download FIS table (CSV)
    output$download_fis_table <- downloadHandler(
      filename = function() paste("fis_results_", Sys.Date(), ".csv", sep = ""),
      content = function(file) {
        shiny::req(fis_boot_results())
        write.csv(
          fis_boot_results()$final_table,
          file,
          row.names = FALSE
        )
      }
    )
    
    ## FIS visualization ----
    make_fis_plot <- function() {
      shiny::req(fis_boot_results())
      ctx <- fis_context()
      
      df <- fis_boot_results()$final_table
      shiny::validate(shiny::need(all(c("ID", "Observed_FIS", "CI_L", "CI_U") %in% names(df)),
                    "FIS results malformed: missing required columns."))
      
      overall_row <- df %>% dplyr::filter(ID == "Overall")
      df          <- df %>% dplyr::filter(ID != "Overall")

      if (nrow(df) == 0) {
        return(ggplot() + labs(title = "No FIS data available") + theme_minimal())
      }

      p <- ggplot(df, aes(x = ID, y = Observed_FIS)) +
        geom_point(size = 3.5, color = "#ff9800") +
        geom_errorbar(aes(ymin = CI_L, ymax = CI_U), width = 0.2,
                      color = "#ff9800", linewidth = 0.9) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        labs(
          title = ctx$title,
          x = ctx$id_label,
          y = "FIS (heterozygote deficit/excess; W&C estimator)"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(face = "bold", hjust = 0.5)
        )

      # Add Overall FIS as a labelled dotted horizontal line
      if (nrow(overall_row) == 1L && is.finite(overall_row$Observed_FIS)) {
        ov <- overall_row$Observed_FIS
        p <- p +
          geom_hline(yintercept = ov, linetype = "dotted",
                     color = "#2c3e50", linewidth = 1) +
          ggplot2::annotate("text", x = Inf, y = ov,
                            label = sprintf("Overall FIS = %.4f", ov),
                            hjust = 1.05, vjust = -0.4,
                            color = "#2c3e50", size = 3.5, fontface = "italic")
      }
      p
    }
    
    output$fis_plot <- renderPlot({ make_fis_plot() 
    })
    
    ### Download FIS plot ----
    
    output$download_fis_plot <- downloadHandler(
      filename = function() paste0("fis_plot_", Sys.Date(), ".png"),
      content = function(file) {
        p <- make_fis_plot()
        ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
      }
    )
    
    
    
    ## ---- Locus × Population cross-table ----
    fis_locus_pop_r <- reactiveVal(NULL)

    observeEvent(input$run_fis_locus_pop, {
      db_ready()
      shiny::req(hf_mat_r(), base_r(), con_r())

      n_perm   <- isolate(input$fis_lp_n_perm)
      base     <- base_r()
      mat      <- hf_mat_r()
      con      <- con_r()

      mat <- as.matrix(mat); storage.mode(mat) <- "integer"
      locus_names <- colnames(mat)[-1L]
      L           <- length(locus_names)

      pop_df <- DBI::dbGetQuery(con, sprintf(
        "SELECT DISTINCT Population FROM %s WHERE Population IS NOT NULL ORDER BY Population",
        sql_ident(con, tbl_meta_r())
      ))
      pop_codes <- as.character(sort(unique(mat[, 1])))
      pop_names <- as.character(pop_df$Population[seq_along(pop_codes)])
      pop_lookup <- stats::setNames(pop_names, pop_codes)

      np      <- length(pop_names)
      fis_m   <- matrix(NA_real_, nrow = L, ncol = np, dimnames = list(locus_names, pop_names))
      pval_m  <- matrix(NA_real_, nrow = L, ncol = np, dimnames = list(locus_names, pop_names))

      withProgress(message = "Computing FIS by locus \u00d7 population\u2026", value = 0, {
        for (pi in seq_along(pop_codes)) {
          incProgress(1 / np, detail = pop_names[pi])
          code <- as.integer(pop_codes[pi])
          pname <- pop_names[pi]

          # Sub-matrix for this population only; recode pop to 1
          sub <- mat[mat[, 1] == code, , drop = FALSE]
          sub[, 1] <- 1L

          if (nrow(sub) < 2L) next

          # Observed WC84 FIS per locus
          res <- fis_wc_cpp(sub, base = as.integer(base))
          fis_m[, pname] <- as.numeric(res$FIS)

          # Permutation p-values per locus
          if (!is.null(n_perm) && n_perm > 0L) {
            perm <- batch_permute_wc_fis(
              dat            = sub,
              pop_col_1based = 1L,
              base           = as.integer(base),
              B              = as.integer(n_perm)
            )
            perm_loci <- perm[, seq_len(ncol(perm) - 1L), drop = FALSE]
            for (li in seq_len(L)) {
              obs_l <- fis_m[li, pname]
              if (!is.finite(obs_l)) next
              permj <- perm_loci[, li]
              permj <- permj[is.finite(permj)]
              if (length(permj) == 0L) next
              ge <- sum(abs(permj) >= abs(obs_l))
              pval_m[li, pname] <- (ge + 1L) / (length(permj) + 1L)
            }
          }
        }
      })

      fis_locus_pop_r(list(fis = fis_m, pval = pval_m,
                           pop_names = pop_names, locus_names = locus_names))
    })

    output$fis_locus_pop_obs <- DT::renderDT({
      shiny::req(fis_locus_pop_r())
      df <- as.data.frame(round(fis_locus_pop_r()$fis, 4))
      df <- tibble::rownames_to_column(df, "Locus")
      DT::datatable(df, rownames = FALSE, filter = "top",
                    extensions = "Buttons",
                    options = list(pageLength = 25, scrollX = TRUE,
                                   dom = "Bfrtip", buttons = c("copy"))) %>%
        DT::formatRound(columns = fis_locus_pop_r()$pop_names, digits = 4)
    })

    output$fis_locus_pop_pval <- DT::renderDT({
      shiny::req(fis_locus_pop_r())
      df <- as.data.frame(round(fis_locus_pop_r()$pval, 4))
      df <- tibble::rownames_to_column(df, "Locus")
      DT::datatable(df, rownames = FALSE, filter = "top",
                    extensions = "Buttons",
                    options = list(pageLength = 25, scrollX = TRUE,
                                   dom = "Bfrtip", buttons = c("copy"))) %>%
        DT::formatRound(columns = fis_locus_pop_r()$pop_names, digits = 4) %>%
        DT::formatStyle(
          columns    = fis_locus_pop_r()$pop_names,
          background = DT::styleInterval(c(0.01, 0.05),
                                         c("#f8d7da", "#fff3cd", "white"))
        )
    })

    output$download_fis_locus_pop <- downloadHandler(
      filename = function() paste0("fis_locus_by_pop_", Sys.Date(), ".csv"),
      content  = function(file) {
        shiny::req(fis_locus_pop_r())
        r    <- fis_locus_pop_r()
        fis  <- as.data.frame(round(r$fis,  4))
        pval <- as.data.frame(round(r$pval, 4))
        fis  <- tibble::rownames_to_column(fis,  "Locus")
        pval <- tibble::rownames_to_column(pval, "Locus")
        names(pval)[-1] <- paste0(names(pval)[-1], "_pval")
        write.csv(merge(fis, pval, by = "Locus", sort = FALSE), file, row.names = FALSE)
      }
    )
    output$download_fis_locus_pop_txt <- downloadHandler(
      filename = function() paste0("fis_locus_by_pop_", Sys.Date(), ".txt"),
      content  = function(file) {
        shiny::req(fis_locus_pop_r())
        r    <- fis_locus_pop_r()
        fis  <- as.data.frame(round(r$fis,  4))
        pval <- as.data.frame(round(r$pval, 4))
        fis  <- tibble::rownames_to_column(fis,  "Locus")
        pval <- tibble::rownames_to_column(pval, "Locus")
        names(pval)[-1] <- paste0(names(pval)[-1], "_pval")
        write.table(merge(fis, pval, by = "Locus", sort = FALSE),
                    file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )

    ## Per-allele F-stats table ----
    output$fis_allele_table <- DT::renderDT({
      shiny::req(fis_allele_results())
      df <- fis_allele_results()
      DT::datatable(
        df,
        rownames  = FALSE,
        filter    = "top",
        extensions = "Buttons",
        options   = list(
          pageLength = 25,
          scrollX    = TRUE,
          dom        = "Bfrtip",
          buttons    = list("colvis")
        )
      ) %>%
        DT::formatRound(columns = c("Freq", "FIS", "FST", "FIT"), digits = 4) %>%
        DT::formatStyle(
          "FIS",
          backgroundColor = DT::styleInterval(
            c(-0.1, 0.1, 0.3),
            c("#d4edda", "#fff3cd", "#f8d7da", "#721c24")
          )
        ) %>%
        DT::formatStyle(
          "FST",
          backgroundColor = DT::styleInterval(
            c(0.05, 0.15, 0.25),
            c("#d4edda", "#fff3cd", "#f8d7da", "#721c24")
          )
        )
    })

    ### Download per-allele table ----
    output$download_fis_allele_table <- downloadHandler(
      filename = function() paste0("fis_per_allele_", Sys.Date(), ".csv"),
      content  = function(file) {
        shiny::req(fis_allele_results())
        utils::write.csv(fis_allele_results(), file, row.names = FALSE)
      }
    )

    output$download_fis_allele_table_txt <- downloadHandler(
      filename = function() paste0("fis_per_allele_", Sys.Date(), ".txt"),
      content  = function(file) {
        shiny::req(fis_allele_results())
        utils::write.table(fis_allele_results(), file,
                           row.names = FALSE, sep = "\t", quote = FALSE)
      }
    )

    ## Compute per-allele F-statistics (independent of bootstrap analyses) ----
    observeEvent(input$compute_allele_fstats, {
      shiny::req(hf_mat_r(), base_r())
      tryCatch({
        allele_mat  <- as.matrix(hf_mat_r())
        storage.mode(allele_mat) <- "integer"
        allele_base <- as.integer(base_r())
        fis_allele_results(wc84_per_allele_fstats_cpp(
          dat          = allele_mat,
          pop_col      = 0L,
          base         = allele_base,
          missing_code = 0L
        ))
        showNotification("Per-allele F-statistics computed successfully!", type = "message")
      }, error = function(e) {
        fis_allele_results(NULL)
        showNotification(paste("Error computing per-allele F-statistics:", e$message), type = "error")
      })
    })

    # ==================================== FIT SECTION ANALYSIS ===============================================
    ## FIT Analysis reactives ----
    fit_boot_results <- reactiveVal(NULL)
    fit_boot_timing  <- reactiveVal(NULL)
    fit_perm_results <- reactiveVal(NULL)
    ## FIT boot and perm analysis ----
    run_fit_analysis <- function(
    n_perm = 1000,
    n_boot = 1000,
    conf_level = 0.95,
    missing_code = 0L,
    cpp_file = "src/fit_permute_bootstrap_wc.cpp",
    cpp_verbose = FALSE
    ) {
      
      
      # DB-first source
      db_ready()
      mat  <- hf_mat_r()
      base <- base_r()
      
      mat <- as.matrix(mat)
      storage.mode(mat) <- "integer"
      shiny::validate(
        shiny::need(is.integer(mat), "hf_mat_r() must return an integer matrix"),
        shiny::need(ncol(mat) >= 2L, "hf_mat_r() must be pop   >= 1 locus"),
        shiny::need(all(mat[, 1] > 0, na.rm = TRUE), "Population codes must be positive integers (1..K)"),
        shiny::need(isTRUE(is.finite(base)) && base > 1L, "base_r() returned invalid base")
      )
      
      locus_names <- colnames(mat)[-1L]
      if (is.null(locus_names) || length(locus_names) == 0L) {
        locus_names <- paste0("L", seq_len(ncol(mat) - 1L))
      }
      
      # observed_wc84_stats_cpp: WC84 per-locus FST/FIT/HI/HS for the FIT bootstrap pipeline.
      # (observed_wc84_stats_cpp from the OpenMP file is the full-stat version used by the FST section.)
      obs_stats <- observed_wc84_stats_cpp(
        dat            = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base)
      )

      observed_fit <- obs_stats$FIT
      names(observed_fit) <- locus_names
      
      boot_res <- boot_wc84_stats_popblock_cpp(
        mat_int        = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base),
        B              = as.integer(n_boot),
        conf_level     = conf_level,
        seed           = 1L,
        n_threads      = 0L
      )
      
      CI_FIT <- boot_res$CI_FIT
      ci_lower <- as.numeric(CI_FIT["lo", ])
      ci_upper <- as.numeric(CI_FIT["hi", ])
      
      fit_boot_mat <- boot_res$FIT_boot
      boot_mean <- as.numeric(colMeans(fit_boot_mat, na.rm = TRUE))
      perm_res <- NULL
      perm_fit_mat <- NULL
      
      pvals_fit <- rep(NA_real_, length(locus_names))
      names(pvals_fit) <- locus_names
      
      if (!is.null(n_perm) && n_perm > 0) {
        
        # FSTAT's "Randomising alleles overall samples" - global allele shuffle
        perm_res <- batch_permute_fit_global(
          dat            = mat,
          pop_col_1based = 1L,
          missing_code   = as.integer(missing_code),
          base           = as.integer(base),
          B              = as.integer(n_perm)
        )
        
        perm_fit_mat <- perm_res$FIT_perm
        colnames(perm_fit_mat) <- locus_names
        
        pvals_fit <- perm_res$p_FIT
        names(pvals_fit) <- locus_names
      }
      
      res_loci <- data.frame(
        ID          = locus_names,
        Observed_FIT = as.numeric(observed_fit),
        Boot_Mean    = boot_mean,
        CI_L         = ci_lower,
        CI_U         = ci_upper,
        P_value      = as.numeric(pvals_fit),
        stringsAsFactors = FALSE
      )
      
      # ---- Overall (FIT) ----
      alpha <- (1 - conf_level) / 2
      
      # WC84 ratio-of-sums overall FIT - NOT the mean of per-locus FIT values
      overall_obs <- as.numeric(obs_stats$FIT_overall_ratio_of_sums)
      
      # bootstrap overall distribution (mean across loci per bootstrap replicate)
      boot_overall <- rowMeans(fit_boot_mat, na.rm = TRUE)
      boot_overall <- boot_overall[is.finite(boot_overall)]
      
      # overall bootstrap mean consistent with locus-level Boot_Mean
      overall_boot_mean <- mean(boot_mean, na.rm = TRUE)
      
      if (length(boot_overall) < 10) {
        overall_ci_l <- NA_real_
        overall_ci_u <- NA_real_
      } else {
        overall_ci_l <- as.numeric(stats::quantile(boot_overall, probs = alpha,     na.rm = TRUE, type = 7))
        overall_ci_u <- as.numeric(stats::quantile(boot_overall, probs = 1 - alpha, na.rm = TRUE, type = 7))
      }
      
      # overall p-value from permutation (two-sided abs), if available
      overall_p <- NA_real_
      if (!is.null(perm_res)) {
        # p-value already computed in C++ using ratio-of-sums overall FIT
        overall_p <- as.numeric(perm_res$p_FIT_overall)
      }
      
      overall_row <- data.frame(
        ID          = "Overall",
        Observed_FIT = overall_obs,
        Boot_Mean    = overall_boot_mean,
        CI_L         = overall_ci_l,
        CI_U         = overall_ci_u,
        P_value      = overall_p,
        stringsAsFactors = FALSE
      )
      
      
      final_df <- rbind(res_loci, overall_row)
      
      list(
        final_table = final_df,
        observed_fit = observed_fit,
        permutation_results = perm_fit_mat,
        bootstrap_results = list(
          FIT_boot = fit_boot_mat,
          CI_FIT   = CI_FIT
        ),
        metadata = list(
          n_loci = length(locus_names),
          n_permutations = n_perm,
          n_bootstrap = n_boot,
          conf_level = conf_level,
          base = base,
          locus_names = locus_names,
          pval_method_FIT = if (!is.null(perm_res)) perm_res$pval_method_FIT else NA_character_
        )
      )
    }
    
    ## FIT observeEvent bootstrap and permutation (button) ====
    
    observeEvent(input$Run_FIT_Analysis, {
      db_ready()
      
      if (input$n_perm_fit < 10 || input$n_boot_fit < 10) {
        showNotification(
          "Number of permutations and bootstrap replicates should be at least 10",
          type = "warning"
        )
      }
      
      waiter <- Waiter$new(
        id    = c("fit_results_table", "fit_permutation_plot", "fit_bootstrap_plot"),
        html  = spin_3(),
        color = transparent(0.7)
      )
      waiter$show()
      on.exit(waiter$hide(), add = TRUE)
      
      tryCatch({
        
        start_time <- Sys.time()
        shinyWidgets::updateProgressBar(session, "fit_progress", value = 5)
        
        shinyWidgets::updateProgressBar(session, "fit_progress", value = 12)
        
        # Run analysis
        results <- run_fit_analysis(
          n_perm     = input$n_perm_fit,
          n_boot     = input$n_boot_fit,
          conf_level = input$conf_level_fit,
          missing_code   = 0L,
          cpp_file       = "src/fit_permute_bootstrap_wc.cpp",
          cpp_verbose    = FALSE
        )
        
        shinyWidgets::updateProgressBar(session, "fit_progress", value = 100)
        
        fit_boot_timing(round(difftime(Sys.time(), start_time, units = "secs"), 1))
        fit_boot_results(results)
        
        # If you use perm_results() elsewhere for plotting, keep this:
        fit_perm_results(results$permutation_results)
        
        
        showNotification("Bootstrap FIT analysis completed successfully!", type = "message")
        
      }, error = function(e) {
        showNotification(paste("Error in FIT analysis:", e$message), type = "error")
        fit_boot_results(NULL)
        fit_perm_results(NULL)
      })
    })
    
    
    ## ===== FIT value boxes =====
    ### Global FIT =====
    output$global_fit_box <- renderValueBox({
      shiny::req(fit_boot_results())
      
      ft <- fit_boot_results()$final_table %>%
        dplyr::filter(ID == "Overall") %>%
        dplyr::pull(Observed_FIT)
      
      display <- ifelse(abs(ft) < 0.0001, "\u2248 0.0000", format(round(ft, 4), nsmall = 4))
      
      # Same threshold scheme as FIS (tune if you want FIT-specific cutoffs)
      color <- ifelse(ft > 0.1, "maroon", ifelse(ft > 0.05, "orange", "aqua"))
      
      valueBox(
        value = display,
        subtitle = HTML("<small>FIT<br>global</small>"),
        color = color,
        icon = icon("dna")
      )
    })
    
    ### Global p-value =====
    output$global_fit_pvalue_box <- renderValueBox({
      shiny::req(fit_boot_results())

      df <- fit_boot_results()$final_table
      p <- df$P_value[df$ID == "Overall"][1]

      display <- if (is.na(p)) {
        "N/A"
      } else if (p < 0.0001) {
        "< 0.0001"
      } else if (p < 0.001) {
        "< 0.001"
      } else {
        format(round(p, 4), nsmall = 4)
      }

      color <- if (is.na(p)) {
        "red"
      } else if (p < 0.001) {
        "red"
      } else if (p < 0.05) {
        "yellow"
      } else {
        "green"
      }

      valueBox(
        value = display,
        subtitle = HTML("<small>Global <i>p</i>-value<br>Two-sided |FIT| permutation</small>"),
        color = color,
        icon = icon("balance-scale"),
        width = NULL
      )
    })
    ### Significant loci =====
    output$significant_loci_fit_box <- renderValueBox({
      shiny::req(fit_boot_results())
      
      fit_data <- fit_boot_results()$final_table |>
        dplyr::filter(ID != "Overall")
      
      total_loci <- nrow(fit_data)
      
      if (total_loci > 0 && "P_value" %in% names(fit_data)) {
        sig_loci <- sum(fit_data$P_value < 0.05, na.rm = TRUE)
        pct <- round(100 * sig_loci / total_loci, 1)
      } else {
        sig_loci <- 0
        pct <- 0
      }
      
      color <- if (sig_loci > 0) "yellow" else "aqua"
      
      valueBox(
        value = paste0(sig_loci, " / ", total_loci),
        subtitle = HTML(paste0("<small>Significant loci (p&lt;0.05)<br>", pct, "% of total</small>")),
        color = color,
        icon = icon("vial"),
        width = NULL
      )
    })
    ### Computation time
    output$analysis_time_fit_box <- renderValueBox({
      shiny::req(fit_boot_timing())
      
      time_sec <- fit_boot_timing()
      time_display <- ifelse(time_sec < 60, 
                             paste0(time_sec, " s"), 
                             paste0(round(time_sec / 60, 1), " min"))
      
      valueBox(
        value = time_display,
        subtitle = HTML("<small>Computation Time<br>FIT Analysis</small>"),
        color = "light-blue",
        icon = icon("clock"),
        width = NULL
      )
    })
    
    ## FIT Table results 
    output$fit_results_table <- DT::renderDT({
      shiny::req(fit_boot_results())
      
      df <- fit_boot_results()$final_table %>%
        dplyr::select(dplyr::any_of(c("ID","Observed_FIT","Boot_Mean","P_value","CI_L","CI_U")))
      
      # force Overall last
      if ("Overall" %in% df$ID) {
        df <- rbind(
          df[df$ID != "Overall", , drop = FALSE],
          df[df$ID == "Overall", , drop = FALSE]
        )
      }
      
      pretty <- c(
        ID           = "Locus",
        Observed_FIT = "Observed FIT",
        Boot_Mean    = "Bootstrap mean",
        P_value      = "P-value",
        CI_L         = "CI lower",
        CI_U         = "CI upper"
      )
      
      DT::datatable(
        df,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = c("copy"),
          pageLength = 15,
          scrollX = TRUE
        ),
        rownames = FALSE,
        colnames = unname(pretty[names(df)])
      ) %>%
        DT::formatRound(
          columns = intersect(c("Observed_FIT","Boot_Mean","P_value","CI_L","CI_U"), names(df)),
          digits = 4
        )
    })
    
    
    
    ## FIT visualization ====
    output$fit_plot <- renderPlot({
      shiny::req(fit_boot_results())
      
      df <- fit_boot_results()$final_table %>%
        dplyr::filter(ID != "Overall")
      
      # robust guards
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        return(ggplot2::ggplot() +
                 ggplot2::labs(title = "No FIT data available") +
                 ggplot2::theme_minimal())
      }
      
      # mark significant loci if P_value exists
      if ("P_value" %in% names(df)) {
        df <- df %>%
          dplyr::mutate(Significant = !is.na(P_value) & P_value < 0.05)
      } else {
        df$Significant <- FALSE
      }
      
      ggplot2::ggplot(df, ggplot2::aes(x = ID, y = Observed_FIT)) +
        ggplot2::geom_point(ggplot2::aes(shape = Significant), size = 3, color = "#ff9800") +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = CI_L, ymax = CI_U), width = 0.2, color = "#ff9800") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::labs(
          title = "FIT estimates with confidence intervals",
          x = "Locus",
          y = "FIT estimate",
          shape = "p < 0.05"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          plot.title  = ggplot2::element_text(face = "bold", hjust = 0.5)
        )
    })
    
    
    
    # ===== FIT DOWNLOAD HANDLERS =====
    output$download_fit_table <- downloadHandler(
      filename = function() {
        paste("fit_results_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        shiny::req(fit_boot_results())
        write.csv(fit_boot_results()$final_table, file, row.names = FALSE)
      }
    )
    
    output$download_fit_table_txt <- downloadHandler(
      filename = function() {
        paste("fit_results_", Sys.Date(), ".txt", sep = "")
      },
      content = function(file) {
        shiny::req(fit_boot_results())
        write.table(fit_boot_results()$final_table, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    
    output$download_fit_plot <- downloadHandler(
      filename = function() {
        paste("fit_plot_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        shiny::req(fit_boot_results())
        
        df <- fit_boot_results()$final_table %>%
          dplyr::filter(ID != "Overall")
        
        if (nrow(df) == 0) {
          p <- ggplot() + labs(title = "No FIT data available") + theme_minimal()
        } else {
          p <- ggplot(df, aes(x = ID, y = Observed_FIT)) +
            geom_point(size = 3, color = "#ff9800") +
            geom_errorbar(aes(ymin = CI_L, ymax = CI_U), width = 0.2, color = "#ff9800") +
            geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
            labs(
              title = "FIT Estimates with Confidence Intervals",
              x = "Locus",
              y = "FIT Estimate"
            ) +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(face = "bold", hjust = 0.5)
            )
        }
        
        ggsave(file, plot = p, width = 12, height = 6, dpi = 300)
      }
    )
    
    
    
    # ==================================== FST SECTION ANALYSIS ===============================================
    
    # hard guard: fail early if C++ funcs are missing
    stopifnot(
      exists("observed_wc84_stats_cpp",   mode = "function"),  # OpenMP full-stat version (FST section)
      exists("batch_permute_wc84_fst_parallel", mode = "function"),
      exists("boot_popblock_wc84_parallel",     mode = "function")
    )
    
    
    ## FST Analysis reactives ----
    fst_boot_results <- reactiveVal(NULL)
    fst_boot_timing  <- reactiveVal(NULL)
    fst_perm_results <- reactiveVal(NULL)
    fst_parallel_meta <- reactiveVal(NULL)
    
    run_bootstrap_fst_analysis <- function(n_perm, n_boot, conf_level, missing_code = 0L) {
      db_ready()
      mat  <- hf_mat_r()
      base <- base_r()
      
      mat <- as.matrix(mat)
      storage.mode(mat) <- "integer"
      shiny::validate(
        shiny::need(is.integer(mat), "hf_mat_r() must return an integer matrix"),
        shiny::need(ncol(mat) >= 2L, "shiny::need pop + at least 1 locus."),
        shiny::need(all(mat[,1] > 0, na.rm = TRUE), "Population codes must be positive integers (1..K)"),
        shiny::need(isTRUE(is.finite(base)) && base > 1L, "Invalid base from params."),
        shiny::need(isTRUE(is.finite(conf_level)) && conf_level > 0 && conf_level < 1, "conf_level must be in (0,1)."),
        shiny::need(isTRUE(is.finite(n_boot)) && n_boot >= 10, "n_boot must be >= 10.")
      )
      
      
      # locus names (everything except pop column)
      loci_names <- colnames(mat)[-1L]
      if (is.null(loci_names) || length(loci_names) == 0L) {
        loci_names <- paste0("L", seq_len(ncol(mat) - 1L))
      }
      
      # ---------------------------#
      # 1) Observed (WC84) from C++
      # ---------------------------#
      obs_res <- .step("observed_wc84_stats_cpp()", observed_wc84_stats_cpp(
        dat            = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base)
      ))
      
      fst_obs_vec <- obs_res$FST
      fst_obs_by_locus <- as.numeric(fst_obs_vec)
      
      loc_obs <- as.character(obs_res$locus_names %||% loci_names)
      names(fst_obs_by_locus) <- loc_obs
      
      # overall FST (WC84 multilocus) = sum(a)/sum(a+b+c)
      fst_obs_overall <- as.numeric(obs_res$FST_overall_ratio_of_sums)
      
      # ---------------------------#
      # 1b) Locus bootstrap (resample loci with replacement)
      # ---------------------------#
      locus_comp <- wc84_locus_components_cpp(
        dat            = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base)
      )
      locus_boot_tbl <- locus_bootstrap_wc84_cpp(
        A          = locus_comp$A,
        Bv         = locus_comp$B,
        C          = locus_comp$C,
        HS         = locus_comp$HS,
        HT         = locus_comp$HT,
        B_reps     = as.integer(n_boot),
        conf_level = conf_level,
        seed       = .seed(),
        n_threads  = .n_threads()
      )

      # ---------------------------#
      # 2) Permutation (single call, fast) from C++
      # ---------------------------#
      perm_res <- NULL
      perm_mat <- NULL
      pvals_by_locus <- rep(NA_real_, length(loci_names))
      pval_overall <- NA_real_
      
      if (!is.null(n_perm) && n_perm > 0) {
        perm_res <-  .step("batch_permute_wc84_fst_auto()", batch_permute_wc84_fst_auto(
          dat            = mat,
          pop_col_1based = 1L,
          missing_code   = as.integer(missing_code),
          base           = as.integer(base),
          B              = as.integer(n_perm),
          n_threads      = .n_threads(),
          seed           = .seed(),
          pval_method    = "greater",              # FST: usually one-sided
          perm_scheme    = "permute_pop_labels",    # breaks pop structure (H0: no structure)
          debug         = FALSE
        ))

        fst_parallel_meta(attr(perm_res, "parallel") %||% NULL)

        perm_mat <- perm_res$FST_perm
        # Ensure colnames align with loci (prefer locus_names from C++ when present)
        if (!is.null(perm_res$locus_names)) colnames(perm_mat) <- as.character(perm_res$locus_names)
        
        pvals_by_locus <- as.numeric(perm_res$p_G)
        if (!is.null(perm_res$locus_names)) names(pvals_by_locus) <- as.character(perm_res$locus_names)
        
        pval_overall <- as.numeric(perm_res$p_FST_overall)
      }
      
      # ---------------------------#
      # 4) Bootstrap (pop blocks) for overall CI (reflected/basic)
      # ---------------------------#
      boot_pop_res <- .step("boot_popblock_wc84_fst_auto()",boot_popblock_wc84_fst_auto(
        mat            = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base),
        B              = as.integer(n_boot),
        n_threads      = .n_threads(),
        seed           = .seed(),
        debug          = FALSE
      ))

      fst_parallel_meta(attr(boot_pop_res, "parallel") %||% fst_parallel_meta())

      perm_parallel_meta <- if (!is.null(perm_res)) attr(perm_res, "parallel") else NULL
      boot_parallel_meta <- attr(boot_pop_res, "parallel")

      # Pop-block bootstrap matrix (per-locus replicates)
      loc <- as.character(boot_pop_res$locus_names)
      
      boot_pop_mat <- boot_pop_res$FST_boot
      hs_boot_mat  <- boot_pop_res$HS_boot
      ht_boot_mat  <- boot_pop_res$HT_boot
      
      colnames(boot_pop_mat) <- loc
      colnames(hs_boot_mat)  <- loc
      colnames(ht_boot_mat)  <- loc
      
      # Overall bootstrap vector MUST be WC84 multilocus ratio-of-sums.
      boot_overall <- boot_pop_res$FST_overall_boot
      if (is.null(boot_overall) || !is.numeric(boot_overall)) {
        stop("C++ boot_popblock_wc84_fst() did not return numeric FST_overall_boot. Recompile the correct fst_permute_bootstrap.cpp.")
      }
      
      # Summarise using C++ (NORMAL bootstrap CI, FSTAT-like)
      loc <- as.character(boot_pop_res$locus_names)
      
      fst_obs_vec_aligned <- as.numeric(fst_obs_vec)
      names(fst_obs_vec_aligned) <- as.character(obs_res$locus_names %||% loc)
      fst_obs_vec_aligned <- fst_obs_vec_aligned[match(loc, names(fst_obs_vec_aligned))]
      
      sum_pop <- .step("summarize_boot_ci(FST)",summarize_boot_ci(
        boot_mat     = boot_pop_mat,
        obs          = fst_obs_vec_aligned,
        obs_overall  = fst_obs_overall,
        boot_overall = boot_overall,
        confidence   = conf_level
      ))
      
      # =========================#
      # HS / HT CI + means
      # =========================#
      
      # Per-locus bootstrap matrices
      hs_boot_mat <- boot_pop_res$HS_boot
      ht_boot_mat <- boot_pop_res$HT_boot
      
      # Make sure colnames are consistent
      colnames(hs_boot_mat) <- as.character(boot_pop_res$locus_names)
      colnames(ht_boot_mat) <- as.character(boot_pop_res$locus_names)
      
      # Observed vectors (named by locus)
      loc <- as.character(boot_pop_res$locus_names)
      
      hs_obs_vec <- as.numeric(boot_pop_res$HS_obs)
      names(hs_obs_vec) <- loc
      hs_obs_vec <- hs_obs_vec[match(loc, names(hs_obs_vec))]
      
      ht_obs_vec <- as.numeric(boot_pop_res$HT_obs)
      names(ht_obs_vec) <- loc
      ht_obs_vec <- ht_obs_vec[match(loc, names(ht_obs_vec))]
      
      
      # Overall observed + overall bootstrap vectors (from C++)
      hs_obs_overall <- as.numeric(boot_pop_res$HS_overall_obs)
      ht_obs_overall <- as.numeric(boot_pop_res$HT_overall_obs)
      
      hs_overall_boot <- boot_pop_res$HS_overall_boot
      ht_overall_boot <- boot_pop_res$HT_overall_boot
      
      # ── Population-block (subsamples) bootstrap CI ──────────────────────────
      sum_hs <- .step("summarize_boot_ci(HS)", summarize_boot_ci(
        boot_mat     = hs_boot_mat,
        obs          = hs_obs_vec,
        obs_overall  = hs_obs_overall,
        boot_overall = hs_overall_boot,
        confidence   = conf_level
      ))

      sum_ht <- .step("summarize_boot_ci(HT)", summarize_boot_ci(
        boot_mat     = ht_boot_mat,
        obs          = ht_obs_vec,
        obs_overall  = ht_obs_overall,
        boot_overall = ht_overall_boot,
        confidence   = conf_level
      ))

      # ── Individual bootstrap CI for HS (resample individuals within pops) ──
      indiv_boot_res <- .step("boot_indiv_hs_cpp()", boot_indiv_hs_cpp(
        dat            = mat,
        pop_col_1based = 1L,
        missing_code   = as.integer(missing_code),
        base           = as.integer(base),
        B              = as.integer(n_boot),
        seed           = .seed(),
        n_threads      = .n_threads()
      ))
      hs_indiv_boot_mat <- indiv_boot_res$HS_boot
      colnames(hs_indiv_boot_mat) <- as.character(indiv_boot_res$locus_names)
      hs_indiv_overall_boot <- indiv_boot_res$HS_overall_boot

      sum_hs_indiv <- .step("summarize_boot_ci(HS indiv)", summarize_boot_ci(
        boot_mat     = hs_indiv_boot_mat,
        obs          = hs_obs_vec,
        obs_overall  = hs_obs_overall,
        boot_overall = hs_indiv_overall_boot,
        confidence   = conf_level
      ))

      loc_fst <- loc
      loc_hs  <- loc
      loc_ht  <- loc

      if (!is.null(names(fst_obs_by_locus))) {
        fst_obs_by_locus <- fst_obs_by_locus[match(loc_fst, names(fst_obs_by_locus))]
        names(fst_obs_by_locus) <- loc_fst
      }
      if (!is.null(names(pvals_by_locus))) {
        pvals_by_locus <- pvals_by_locus[match(loc_fst, names(pvals_by_locus))]
        names(pvals_by_locus) <- loc_fst
      }

      # ── Helper: build per-locus + Overall HS table for one bootstrap mode ──
      .hs_boot_tbl <- function(sum_obj, indiv_boot_mat, indiv_boot_overall,
                               obs_vec, obs_overall, loc_names) {
        boot_se_loci   <- apply(indiv_boot_mat, 2, sd, na.rm = TRUE)
        boot_se_overall <- sd(as.numeric(indiv_boot_overall), na.rm = TRUE)
        per_locus <- data.frame(
          ID          = loc_names,
          Observed_HS = as.numeric(obs_vec[loc_names]),
          Boot_Mean   = as.numeric(sum_obj$mean),
          Boot_SE     = as.numeric(boot_se_loci),
          CI_L        = as.numeric(sum_obj$ci_lo),
          CI_U        = as.numeric(sum_obj$ci_hi),
          stringsAsFactors = FALSE
        )
        overall <- data.frame(
          ID          = "Overall",
          Observed_HS = as.numeric(obs_overall),
          Boot_Mean   = as.numeric(sum_obj$overall_mean),
          Boot_SE     = as.numeric(boot_se_overall),
          CI_L        = as.numeric(sum_obj$overall_ci_lo),
          CI_U        = as.numeric(sum_obj$overall_ci_hi),
          stringsAsFactors = FALSE
        )
        rbind(per_locus, overall)
      }

      # Table 1: HS by individuals
      hs_indiv_tbl <- .hs_boot_tbl(sum_hs_indiv,
                                    hs_indiv_boot_mat, hs_indiv_overall_boot,
                                    hs_obs_vec, hs_obs_overall, loc_hs)

      # Table 2: HS by populations (block bootstrap)
      hs_pop_tbl <- .hs_boot_tbl(sum_hs,
                                  hs_boot_mat, hs_overall_boot,
                                  hs_obs_vec, hs_obs_overall, loc_hs)

      # Table 3: HS by loci — Overall only (from locus_boot_tbl)
      hs_locus_tbl <- tryCatch({
        row <- locus_boot_tbl[locus_boot_tbl$Statistic == "HS", ,drop = FALSE]
        data.frame(
          Statistic   = "HS",
          Observed    = as.numeric(row$Observed),
          Boot_Mean   = as.numeric(row$Boot_Mean),
          Boot_SE     = as.numeric(row$SE),
          CI_L        = as.numeric(row$CI_L),
          CI_U        = as.numeric(row$CI_U),
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        data.frame(Statistic="HS", Observed=hs_obs_overall,
                   Boot_Mean=NA_real_, Boot_SE=NA_real_,
                   CI_L=NA_real_, CI_U=NA_real_,
                   stringsAsFactors=FALSE)
      })

      # HT tables (populations + locus bootstrap, unchanged logic)
      ht_pop_se_loci   <- apply(ht_boot_mat, 2, sd, na.rm = TRUE)
      ht_pop_se_overall <- sd(as.numeric(ht_overall_boot), na.rm = TRUE)
      ht_tbl <- data.frame(
        ID          = loc_ht,
        Observed_HT = as.numeric(ht_obs_vec[loc_ht]),
        Boot_Mean   = as.numeric(sum_ht$mean),
        Boot_SE     = as.numeric(ht_pop_se_loci),
        CI_L        = as.numeric(sum_ht$ci_lo),
        CI_U        = as.numeric(sum_ht$ci_hi),
        stringsAsFactors = FALSE
      )
      ht_overall_row <- data.frame(
        ID          = "Overall",
        Observed_HT = as.numeric(ht_obs_overall),
        Boot_Mean   = as.numeric(sum_ht$overall_mean),
        Boot_SE     = as.numeric(ht_pop_se_overall),
        CI_L        = as.numeric(sum_ht$overall_ci_lo),
        CI_U        = as.numeric(sum_ht$overall_ci_hi),
        stringsAsFactors = FALSE
      )
      ht_final <- rbind(ht_tbl, ht_overall_row)

      # Table 4: HS per population — observed + individual bootstrap CI
      # Bootstrap unit: individuals resampled with replacement within each population.
      # Observed HS per population = mean of per-locus WC84 HS (same formula as
      # hs_by_pop_locus_from_mat), averaged across loci for that population.
      .hs_one_pop_locus <- function(g, base) {
        ok <- is.finite(g) & g > 0L
        if (sum(ok) <= 1L) return(NA_real_)
        g  <- g[ok]; a1 <- g %/% base; a2 <- g %% base
        ok2 <- a1 > 0L & a2 > 0L
        if (sum(ok2) <= 1L) return(NA_real_)
        a1 <- a1[ok2]; a2 <- a2[ok2]; n <- length(a1)
        cnt <- table(c(a1, a2)); p <- as.numeric(cnt) / (2 * n)
        (2 * n / (2 * n - 1)) * (1 - sum(p^2))
      }

      pop_codes_pp  <- as.integer(mat[, 1])
      pop_levels_pp <- attr(mat, "pop_levels")
      pops_pp <- sort(unique(pop_codes_pp[is.finite(pop_codes_pp) & pop_codes_pp > 0L]))
      n_loci_pp <- ncol(mat) - 1L

      set.seed(.seed())
      hs_per_pop_rows <- lapply(pops_pp, function(pp) {
        idx <- which(pop_codes_pp == pp)
        n_pp <- length(idx)
        pop_name <- if (!is.null(pop_levels_pp) && pp >= 1L && pp <= length(pop_levels_pp))
          as.character(pop_levels_pp[[pp]]) else as.character(pp)

        hs_obs_loci <- vapply(seq_len(n_loci_pp), function(j)
          .hs_one_pop_locus(as.integer(mat[idx, j + 1L]), base), numeric(1))
        hs_obs_mean <- mean(hs_obs_loci, na.rm = TRUE)

        boot_means <- vapply(seq_len(n_boot), function(b) {
          idx_b <- idx[sample.int(n_pp, n_pp, replace = TRUE)]
          hs_b  <- vapply(seq_len(n_loci_pp), function(j)
            .hs_one_pop_locus(as.integer(mat[idx_b, j + 1L]), base), numeric(1))
          mean(hs_b, na.rm = TRUE)
        }, numeric(1))

        alpha <- (1 - conf_level) / 2
        data.frame(
          Population  = pop_name,
          Observed_HS = hs_obs_mean,
          Boot_Mean   = mean(boot_means, na.rm = TRUE),
          Boot_SE     = sd(boot_means, na.rm = TRUE),
          CI_L        = as.numeric(quantile(boot_means, alpha,     na.rm = TRUE)),
          CI_U        = as.numeric(quantile(boot_means, 1 - alpha, na.rm = TRUE)),
          N_loci      = as.integer(sum(!is.na(hs_obs_loci))),
          stringsAsFactors = FALSE
        )
      })
      hs_per_pop_tbl <- if (length(hs_per_pop_rows) > 0)
        do.call(rbind, hs_per_pop_rows)
      else
        data.frame(Population=character(), Observed_HS=numeric(),
                   Boot_Mean=numeric(), Boot_SE=numeric(),
                   CI_L=numeric(), CI_U=numeric(),
                   N_loci=integer(), stringsAsFactors=FALSE)

      # Backward-compat alias used by value box and download handler
      hs_final <- hs_pop_tbl
      
      boot_mean <- as.numeric(sum_pop$mean)
      ci_lower  <- as.numeric(sum_pop$ci_lo)
      ci_upper  <- as.numeric(sum_pop$ci_hi)
      
      names(boot_mean) <- loc_fst
      names(ci_lower)  <- loc_fst
      names(ci_upper)  <- loc_fst
      
      # Overall outputs
      overall_boot_mean <- as.numeric(sum_pop$overall_mean)
      overall_ci_l      <- as.numeric(sum_pop$overall_ci_lo)
      overall_ci_u      <- as.numeric(sum_pop$overall_ci_hi)
      
      # Optional: median from the bootstrap matrix (if you still want it)
      boot_median <- apply(boot_pop_mat, 2, median, na.rm = TRUE)
      
      
      # ---------------------------#
      # 6) Final results table
      # ---------------------------#
      loc <- loc_fst
      if (is.null(loc) || length(loc) == 0L) loc <- loci_names
      
      obs_vec <- fst_obs_by_locus
      if (!is.null(names(obs_vec))) obs_vec <- obs_vec[match(loc, names(obs_vec))]
      
      p_vec <- pvals_by_locus
      if (!is.null(names(p_vec))) p_vec <- p_vec[match(loc, names(p_vec))]
      
      res_tbl <- data.frame(
        ID          = loc,
        Observed_FST = as.numeric(obs_vec),
        Boot_Mean    = as.numeric(boot_mean[loc]),
        Boot_Median  = as.numeric(boot_median[loc]),
        P_value      = as.numeric(p_vec),
        CI_L         = as.numeric(ci_lower[loc]),
        CI_U         = as.numeric(ci_upper[loc]),
        stringsAsFactors = FALSE
      )
      
      # overall_boot_mean already computed from v_overall above (more coherent)
      # if it was NA (too few finite), fallback:
      # Option 1: define an overall "bootstrap median" coherently from boot_overall
      overall_boot_median <- if (exists("boot_overall") && length(boot_overall) > 0) {
        median(boot_overall, na.rm = TRUE)
      } else {
        NA_real_
      }
      
      overall_row <- data.frame(
        ID          = "Overall",
        Observed_FST = fst_obs_overall,
        Boot_Mean    = overall_boot_mean,
        Boot_Median  = overall_boot_median,
        P_value      = pval_overall,
        CI_L         = overall_ci_l,
        CI_U         = overall_ci_u,
        stringsAsFactors = FALSE
      )
      
      # force same column order as res_tbl
      overall_row <- overall_row[, names(res_tbl), drop = FALSE]
      
      final_results <- rbind(res_tbl, overall_row)
      
      list(
        final_table      = final_results,
        locus_boot_table = locus_boot_tbl,
        hs_table     = hs_final,
        hs_indiv_tbl = hs_indiv_tbl,
        hs_pop_tbl   = hs_pop_tbl,
        hs_locus_tbl = hs_locus_tbl,
        hs_per_pop_tbl = hs_per_pop_tbl,
        ht_table    = ht_final,
        observed_fst = fst_obs_by_locus,
        permutation_results = perm_mat,
        boot_parallel_meta = boot_parallel_meta,
        perm_parallel_meta = perm_parallel_meta,
        bootstrap_results = list(
          population = boot_pop_res,
          population_matrix = boot_pop_mat,
          overall_boot = boot_overall,
          hs = sum_hs,
          ht = sum_ht
        ),
        metadata = list(
          parallel_perm = perm_parallel_meta,
          parallel_boot = boot_parallel_meta,
          requested_threads = as.integer(.n_threads()),
          n_loci = length(loci_names),
          n_permutations = n_perm,
          n_bootstrap = n_boot,
          conf_level = conf_level,
          base = base,
          loci_names = loci_names,
          pval_method = if (!is.null(perm_res)) perm_res$pval_method else NA_character_
        )
        
      )
    }
    
    
    ## Run FST bootstrap and permutation (button) ----
    observeEvent(input$run_FST_Analysis, {
      
      db_ready()
      
      if (input$n_boot_fst < 10) {
        showNotification("Number of bootstrap replicates should be at least 10", type = "warning")
        return(NULL)
      }
      
      waiter <- Waiter$new(
        id    = c(session$ns("fst_results_table"), session$ns("fst_plot")),
        html  = spin_3(),
        color = transparent(0.7)
      )
      waiter$show()
      on.exit(waiter$hide(), add = TRUE)
      
      tryCatch({
        start_time <- Sys.time()
        
        shinyWidgets::updateProgressBar(session, "fst_progress", value = 15)
        
        results <- run_bootstrap_fst_analysis(
          n_perm         = input$n_perm_fst,
          n_boot         = input$n_boot_fst,
          conf_level     = input$conf_level_fst,
          missing_code   = 0L
        )
        
        shinyWidgets::updateProgressBar(session, "fst_progress", value = 100)
        
        duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
        fst_boot_timing(duration)
        fst_boot_results(results)
        
        showNotification(
          paste("FST analysis completed successfully! Time:", duration, "seconds"),
          type = "message"
        )
        
      }, error = function(e) {
        fst_boot_results(NULL)
        fst_boot_timing(NULL)
        showNotification(paste("Error in FST analysis:", e$message), type = "error")
      })
    })

    ## Run FST bootstrap and permutation (button from Genetic diversities tab) ----
    observeEvent(input$run_FST_Analysis_div, {

      db_ready()

      if (input$n_boot_fst_div < 10) {
        showNotification("Number of bootstrap replicates should be at least 10", type = "warning")
        return(NULL)
      }

      waiter <- Waiter$new(
        id    = c(session$ns("hs_indiv_table"), session$ns("hs_pop_table"),
                  session$ns("hs_locus_table"), session$ns("ht_results_table")),
        html  = spin_3(),
        color = transparent(0.7)
      )
      waiter$show()
      on.exit(waiter$hide(), add = TRUE)

      tryCatch({
        start_time <- Sys.time()

        shinyWidgets::updateProgressBar(session, "fst_progress_div", value = 15)

        results <- run_bootstrap_fst_analysis(
          n_perm         = input$n_perm_fst_div,
          n_boot         = input$n_boot_fst_div,
          conf_level     = if (!is.null(input$conf_level_fst_div)) input$conf_level_fst_div else 0.95,
          missing_code   = 0L
        )

        shinyWidgets::updateProgressBar(session, "fst_progress_div", value = 100)

        duration <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
        fst_boot_timing(duration)
        fst_boot_results(results)

        showNotification(
          paste("FST analysis completed successfully! Time:", duration, "seconds"),
          type = "message"
        )

      }, error = function(e) {
        fst_boot_results(NULL)
        fst_boot_timing(NULL)
        showNotification(paste("Error in FST analysis:", e$message), type = "error")
      })
    })

    ## ===== FST value boxes =====
    
    ### Global FST ----
    output$global_fst_box <- renderValueBox({
      shiny::req(fst_boot_results())
      
      fst <- fst_boot_results()$final_table %>%
        dplyr::filter(ID == "Overall") %>%
        dplyr::pull(Observed_FST)
      
      display <- ifelse(is.na(fst), "N/A",
                        ifelse(abs(fst) < 0.0001, "\u2248 0.0000", format(round(fst, 4), nsmall = 4)))
      
      # Thresholds are a choice; you can keep yours.
      color <- if (is.na(fst)) {
        "red"
      } else if (fst > 0.25) {
        "red"
      } else if (fst > 0.15) {
        "yellow"
      } else if (fst > 0.05) {
        "green"
      } else {
        "aqua"
      }
      
      valueBox(
        value = display,
        color = color,
        subtitle = HTML("<small>FST<br>global</small>"),
        icon = icon("globe-americas"),
        width = NULL
      )
    })
    
    ### Global p-value ----
    output$global_fst_pvalue_box <- renderValueBox({
      shiny::req(fst_boot_results())
      
      p <- fst_boot_results()$final_table %>%
        dplyr::filter(ID == "Overall") %>%
        dplyr::pull(P_value)
      
      p <- p[1]
      
      display <- if (is.na(p)) {
        "N/A"
      } else if (p < 0.0001) {
        "< 0.0001"
      } else if (p < 0.001) {
        "< 0.001"
      } else {
        format(round(p, 4), nsmall = 4)
      }
      
      color <- if (is.na(p)) {
        "red"
      } else if (p < 0.001) {
        "red"
      } else if (p < 0.05) {
        "yellow"
      } else {
        "green"
      }
      
      valueBox(
        value = display,
        subtitle = HTML("<small>Global <i>p</i>-value<br>Bilateral test</small>"),
        color = color,
        icon = icon("balance-scale"),
        width = NULL
      )
    })
    
    ### Significant loci ----
    output$significant_loci_fst_box <- renderValueBox({
      shiny::req(fst_boot_results())
      
      fst_data <- fst_boot_results()$final_table %>%
        dplyr::filter(ID != "Overall")
      
      total_loci <- nrow(fst_data)
      
      if (total_loci > 0 && "P_value" %in% names(fst_data)) {
        sig_loci <- sum(!is.na(fst_data$P_value) & fst_data$P_value < 0.05)
        pct <- round(100 * sig_loci / total_loci, 1)
      } else {
        sig_loci <- 0
        pct <- 0
      }
      
      color <- if (sig_loci > 0) "yellow" else "aqua"
      
      valueBox(
        value = paste0(sig_loci, " / ", total_loci),
        subtitle = HTML(paste0("<small>Significant loci (p&lt;0.05)<br>", pct, "% of total</small>")),
        color = color,
        icon = icon("vial"),
        width = NULL
      )
    })
    
    ### Computation time ----
    output$analysis_time_fst_box <- renderValueBox({
      shiny::req(fst_boot_timing())
      res <- fst_boot_results()
      
      time_sec <- fst_boot_timing()
      time_display <- ifelse(is.na(time_sec), "N/A",
                             ifelse(time_sec < 60,
                                    paste0(time_sec, " s"),
                                    paste0(round(time_sec / 60, 1), " min")))
      
      # Prefer meta returned by wrappers (tells you if parallel was actually used)
      meta <- NULL
      if (is.list(res)) {
        meta <- res$boot_parallel_meta %||% res$perm_parallel_meta %||% (res$metadata$parallel_boot %||% res$metadata$parallel_perm)
      }
      thr <- suppressWarnings(as.integer(meta$requested_threads %||% (res$metadata$requested_threads %||% NA_integer_)))
      if (!is.finite(thr) || length(thr) != 1L || thr < 1L) thr <- NA_integer_
      
      used_par <- isTRUE(meta$used_parallel)
      backend  <- meta$backend %||% NA_character_
      
      thr_label <- if (is.na(thr)) {
        "Threads used: unknown"
      } else if (isTRUE(used_par)) {
        paste0("Threads used: ", thr, " (parallel)")
      } else {
        paste0("Threads used: ", thr, " (serial)")
      }
      
      valueBox(
        value = time_display,
        subtitle = HTML(paste0("<small>Computation Time<br>", thr_label, "</small>")),
        color = "light-blue",
        icon = icon("clock"),
        width = NULL
      )
    })
    
    ## ===== Diversities tab value boxes (global_fst, global_hs, global_ht, time) =====

    output$global_fst_div_box <- renderValueBox({
      res <- fst_boot_results()
      shiny::req(!is.null(res), !is.null(res$final_table))
      fst <- res$final_table %>%
        dplyr::filter(ID == "Overall") %>%
        dplyr::pull(Observed_FST)
      fst <- if (length(fst) == 0) NA_real_ else fst[[1]]
      display <- if (is.na(fst)) "N/A" else format(round(fst, 4), nsmall = 4)
      color <- if (is.na(fst)) "aqua" else if (fst > 0.15) "maroon" else if (fst > 0.05) "orange" else "aqua"
      valueBox(value = display, subtitle = HTML("<small>Global FST<br>Population subdivision</small>"),
               color = color, icon = icon("sitemap"), width = NULL)
    })

    output$global_hs_box <- renderValueBox({
      res <- fst_boot_results()
      shiny::req(!is.null(res), !is.null(res$hs_table))
      hs <- res$hs_table %>%
        dplyr::filter(ID == "Overall") %>%
        dplyr::pull(Observed_HS)
      hs <- if (length(hs) == 0) NA_real_ else hs[[1]]
      display <- if (is.na(hs)) "N/A" else format(round(hs, 4), nsmall = 4)
      color <- if (is.na(hs)) "aqua" else if (hs > 0.5) "maroon" else if (hs > 0.2) "orange" else "aqua"
      valueBox(value = display, subtitle = HTML("<small>Global HS<br>Within-pop gene diversity</small>"),
               color = color, icon = icon("dna"), width = NULL)
    })

    output$global_ht_box <- renderValueBox({
      res <- fst_boot_results()
      shiny::req(!is.null(res), !is.null(res$ht_table))
      ht <- res$ht_table %>%
        dplyr::filter(ID == "Overall") %>%
        dplyr::pull(Observed_HT)
      ht <- if (length(ht) == 0) NA_real_ else ht[[1]]
      display <- if (is.na(ht)) "N/A" else format(round(ht, 4), nsmall = 4)
      color <- if (is.na(ht)) "aqua" else if (ht > 0.5) "maroon" else if (ht > 0.2) "orange" else "aqua"
      valueBox(value = display, subtitle = HTML("<small>Global HT<br>Total gene diversity</small>"),
               color = color, icon = icon("globe"), width = NULL)
    })

    output$analysis_time_div_box <- renderValueBox({
      shiny::req(fst_boot_timing())
      time_sec <- fst_boot_timing()
      time_display <- if (is.na(time_sec)) "N/A" else if (time_sec < 60)
        paste0(time_sec, " s") else paste0(round(time_sec / 60, 1), " min")
      valueBox(value = time_display,
               subtitle = HTML("<small>Computation Time<br>(shared with Subdivision)</small>"),
               color = "light-blue", icon = icon("clock"), width = NULL)
    })

    ## ===== Locus bootstrap summary table =====
    output$locus_boot_table <- DT::renderDT({
      res <- fst_boot_results()
      shiny::req(is.list(res), !is.null(res$locus_boot_table))
      df <- res$locus_boot_table
      pretty <- c(
        Statistic = "Statistic",
        Observed  = "Observed",
        Boot_Mean = "Bootstrap mean",
        SE        = "Bootstrap SE",
        CI_L      = "CI lower",
        CI_U      = "CI upper"
      )
      DT::datatable(
        df,
        extensions = "Buttons",
        options = list(
          dom = "t",
          pageLength = 10,
          scrollX = TRUE
        ),
        rownames = FALSE,
        colnames = unname(pretty[names(df)])
      ) %>%
        DT::formatRound(
          columns = intersect(
            c("Observed", "Boot_Mean", "SE", "CI_L", "CI_U"), names(df)),
          digits = 4
        )
    })

    ## ===== FST, HT, HS result tables =====
    ### HT #####
    output$ht_results_table <- DT::renderDT({
      res <- fst_boot_results()
      shiny::req(is.list(res), !is.null(res$ht_table))

      df <- res$ht_table %>%
        dplyr::select(dplyr::any_of(c(
          "ID","Observed_HT",
          "Subsamp_CI_L","Subsamp_CI_U",
          "Locus_CI_L","Locus_CI_U"
        )))

      if ("Overall" %in% df$ID)
        df <- rbind(df[df$ID != "Overall",, drop=FALSE],
                    df[df$ID == "Overall",, drop=FALSE])

      pretty_names <- c(
        ID           = "Locus",
        Observed_HT  = "Observed HT",
        Subsamp_CI_L = "Populations CI lower",
        Subsamp_CI_U = "Populations CI upper",
        Locus_CI_L   = "Locus bootstrap CI lower",
        Locus_CI_U   = "Locus bootstrap CI upper"
      )

      DT::datatable(df, extensions = "Buttons",
        options = list(dom="Bfrtip", buttons=c("copy"), pageLength=15, scrollX=TRUE),
        rownames = FALSE,
        colnames = unname(pretty_names[names(df)])
      ) %>%
        DT::formatRound(
          columns = intersect(c("Observed_HT","Subsamp_CI_L","Subsamp_CI_U",
                                "Locus_CI_L","Locus_CI_U"), names(df)),
          digits = 4
        )
    })
    
    ### HS #####
    # Helper: render one of the 3 HS tables (all share the same column structure)
    .render_hs_tbl <- function(df) {
      shiny::req(is.data.frame(df), nrow(df) > 0)
      if ("ID" %in% names(df) && "Overall" %in% df$ID)
        df <- rbind(df[df$ID != "Overall",, drop=FALSE],
                    df[df$ID == "Overall",, drop=FALSE])
      num_cols <- intersect(c("Observed_HS","Boot_Mean","Boot_SE","CI_L","CI_U"), names(df))
      pretty <- c(ID="Locus", Observed_HS="Observed HS",
                  Boot_Mean="Bootstrap mean", Boot_SE="Bootstrap SE",
                  CI_L="CI lower", CI_U="CI upper")
      DT::datatable(df, rownames=FALSE, extensions="Buttons",
        options=list(pageLength=15, scrollX=TRUE, dom="Bfrtip", buttons=c("copy")),
        colnames=unname(pretty[names(df)])
      ) %>% DT::formatRound(columns=num_cols, digits=4)
    }

    output$hs_indiv_table <- DT::renderDT({
      res <- fst_boot_results()
      shiny::req(is.list(res), !is.null(res$hs_indiv_tbl))
      .render_hs_tbl(res$hs_indiv_tbl)
    })

    output$hs_pop_table <- DT::renderDT({
      res <- fst_boot_results()
      shiny::req(is.list(res), !is.null(res$hs_pop_tbl))
      .render_hs_tbl(res$hs_pop_tbl)
    })

    output$hs_locus_table <- DT::renderDT({
      res <- fst_boot_results()
      shiny::req(is.list(res), !is.null(res$hs_locus_tbl))
      df <- res$hs_locus_tbl
      num_cols <- intersect(c("Observed","Boot_Mean","Boot_SE","CI_L","CI_U"), names(df))
      pretty <- c(Statistic="Statistic", Observed="Observed HS",
                  Boot_Mean="Bootstrap mean", Boot_SE="Bootstrap SE",
                  CI_L="CI lower", CI_U="CI upper")
      DT::datatable(df, rownames=FALSE, extensions="Buttons",
        options=list(pageLength=5, scrollX=TRUE, dom="Bfrtip", buttons=c("copy")),
        colnames=unname(pretty[names(df)])
      ) %>% DT::formatRound(columns=num_cols, digits=4)
    })
    
    output$hs_per_pop_table <- DT::renderDT({
      res <- fst_boot_results()
      shiny::req(is.list(res), !is.null(res$hs_per_pop_tbl))
      df <- res$hs_per_pop_tbl
      num_cols <- intersect(c("Observed_HS","Boot_Mean","Boot_SE","CI_L","CI_U"), names(df))
      pretty <- c(Population="Population",
                  Observed_HS="Observed HS", Boot_Mean="Bootstrap mean",
                  Boot_SE="Bootstrap SE", CI_L="CI lower", CI_U="CI upper",
                  N_loci="N loci")
      DT::datatable(df, rownames=FALSE, extensions="Buttons",
        options=list(pageLength=25, scrollX=TRUE, dom="Bfrtip", buttons=c("copy")),
        colnames=unname(pretty[names(df)])
      ) %>% DT::formatRound(columns=num_cols, digits=4)
    })

    ### FST #####
    output$fst_results_table <- DT::renderDT({
      res <- fst_boot_results()
      shiny::req(is.list(res), !is.null(res$final_table))
      
      df <- res$final_table %>%
        dplyr::select(dplyr::any_of(c("ID","Observed_FST","Boot_Mean","P_value","CI_L","CI_U")))
      
      pretty_names <- c(
        ID           = "Locus",
        Observed_FST = "Observed FST",
        Boot_Mean    = "Bootstrap mean",
        P_value      = "P-value",
        CI_L         = "CI lower",
        CI_U         = "CI upper"
      )
      DT::datatable(
        df,
        extensions = "Buttons",
        
        options = list(
          dom = "Bfrtip",
          buttons = c("copy"),
          pageLength = 15,
          scrollX = TRUE
        ),
        rownames = FALSE,
        colnames = unname(pretty_names[names(df)])
      ) %>%
        DT::formatRound(columns = intersect(c("Observed_FST","Boot_Mean","P_value","CI_L","CI_U"), names(df)), digits = 4)
    })
    
    
    ## ===== FST, HT, HS  plots =====

    # Helper: build a combined per-locus + Overall plot for HS or HT.
    # Mirrors FST plot style: size=3, width=0.2, no size/linewidth aesthetics.
    .diversity_plot <- function(full_df, obs_col, ci_l_col, ci_u_col, y_label, title) {
      df_loci   <- full_df %>% dplyr::filter(ID != "Overall")
      df_overall <- full_df %>% dplyr::filter(ID == "Overall")

      if (nrow(df_loci) == 0)
        return(ggplot() + labs(title = paste("No", y_label, "data available")) + theme_minimal())

      locus_levels <- c(df_loci$ID, if (nrow(df_overall) > 0) "Overall")
      df_loci   <- df_loci   %>% dplyr::mutate(ID = factor(ID, levels = locus_levels))
      df_overall <- df_overall %>% dplyr::mutate(ID = factor(ID, levels = locus_levels))

      p <- ggplot(mapping = aes(x = ID, y = .data[[obs_col]])) +
        # per-locus: same style as FST plot
        geom_point(data = df_loci, size = 3, color = "#3498db", shape = 16) +
        geom_errorbar(data = df_loci,
                      aes(ymin = .data[[ci_l_col]], ymax = .data[[ci_u_col]]),
                      width = 0.2, color = "#3498db") +
        labs(title = title, x = "Locus", y = y_label) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title  = element_text(face = "bold", hjust = 0.5))

      # Overall: same size and linewidth, distinct by colour + shape only
      if (nrow(df_overall) > 0 &&
          is.finite(df_overall[[obs_col]][1]) &&
          ci_l_col %in% names(df_overall) &&
          ci_u_col %in% names(df_overall) &&
          is.finite(df_overall[[ci_l_col]][1])) {
        p <- p +
          geom_point(data = df_overall, size = 3, color = "#e74c3c", shape = 18) +
          geom_errorbar(data = df_overall,
                        aes(ymin = .data[[ci_l_col]], ymax = .data[[ci_u_col]]),
                        width = 0.2, color = "#e74c3c")
      }
      p
    }

    ### HT ####
    output$ht_plot <- renderPlot({
      res <- fst_boot_results()
      shiny::req(is.list(res), !is.null(res$ht_table))
      df <- res$ht_table
      ci_l <- if ("Subsamp_CI_L" %in% names(df)) "Subsamp_CI_L" else "CI_L"
      ci_u <- if ("Subsamp_CI_U" %in% names(df)) "Subsamp_CI_U" else "CI_U"
      .diversity_plot(df, "Observed_HT", ci_l, ci_u, "HT",
                      "HT per locus \u2014 populations bootstrap CI")
    })

    ### HS ####
    output$hs_plot <- renderPlot({
      res <- fst_boot_results()
      shiny::req(is.list(res), !is.null(res$hs_pop_tbl))
      .diversity_plot(res$hs_pop_tbl, "Observed_HS", "CI_L", "CI_U", "HS",
                      "HS per locus \u2014 populations block bootstrap CI")
    })
    ### FST ####
    output$fst_plot <- renderPlot({
      res <- fst_boot_results()
      shiny::req(is.list(res), !is.null(res$final_table))
      
      df <- res$final_table
      
      # make sure Overall is last
      loci_only <- df$ID[df$ID != "Overall"]
      df$ID <- factor(df$ID, levels = c(sort(unique(loci_only)), "Overall"))
      
      df <- df %>%
        dplyr::mutate(Significant = !is.na(P_value) & P_value < 0.05)
      
      ggplot(df, aes(x = ID, y = Observed_FST)) +
        geom_point(aes(shape = Significant), size = 3, color = "#3498db") +
        geom_errorbar(aes(ymin = CI_L, ymax = CI_U), width = 0.2, color = "#3498db") +
        # geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        labs(
          title = "FST estimates with confidence intervals",
          x = "Locus",
          y = "FST estimate",
          shape = "p < 0.05"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(face = "bold", hjust = 0.5)
        )
    })
    
    
    ## ===== FST, HT, HS  download handlers =====
    ### HT ####
    output$download_ht_table <- downloadHandler(
      filename = function() sprintf("ht_results_%s.csv", Sys.Date()),
      content = function(file) {
        res <- fst_boot_results()
        shiny::req(is.list(res), !is.null(res$ht_table))
        utils::write.csv(res$ht_table, file, row.names = FALSE)
      }
    )
    output$download_ht_table_txt <- downloadHandler(
      filename = function() sprintf("ht_results_%s.txt", Sys.Date()),
      content = function(file) {
        res <- fst_boot_results()
        shiny::req(is.list(res), !is.null(res$ht_table))
        utils::write.table(res$ht_table, file,
                           sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    output$download_ht_plot <- downloadHandler(
      filename = function() sprintf("ht_plot_%s.png", Sys.Date()),
      content = function(file) {
        res <- fst_boot_results()
        shiny::req(is.list(res), !is.null(res$ht_table))
        
        df <- res$ht_table %>% dplyr::filter(ID != "Overall")
        shiny::req(is.data.frame(df), nrow(df) > 0)
        
        p <- ggplot(df, aes(x = ID, y = Observed_HT)) +
          geom_point(size = 3, color = "#3498db") +
          geom_errorbar(aes(ymin = CI_L, ymax = CI_U), width = 0.2, color = "#3498db") +
          labs(
            title = "HT estimates with confidence intervals",
            x = "Locus",
            y = "HT"
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title  = element_text(face = "bold", hjust = 0.5)
          )
        
        ggsave(file, plot = p, width = 12, height = 6, dpi = 300)
      }
    )
    ### HS ####
    output$download_hs_table <- downloadHandler(
      filename = function() sprintf("hs_pop_results_%s.csv", Sys.Date()),
      content = function(file) {
        res <- fst_boot_results()
        shiny::req(is.list(res), !is.null(res$hs_pop_tbl))
        utils::write.csv(res$hs_pop_tbl, file, row.names = FALSE)
      }
    )
    output$download_hs_table_txt <- downloadHandler(
      filename = function() sprintf("hs_pop_results_%s.txt", Sys.Date()),
      content = function(file) {
        res <- fst_boot_results()
        shiny::req(is.list(res), !is.null(res$hs_pop_tbl))
        utils::write.table(res$hs_pop_tbl, file, sep="\t", row.names=FALSE, quote=FALSE)
      }
    )
    output$download_hs_plot <- downloadHandler(
      filename = function() sprintf("hs_plot_%s.png", Sys.Date()),
      content = function(file) {
        res <- fst_boot_results()
        shiny::req(is.list(res), !is.null(res$hs_pop_tbl))
        p <- .diversity_plot(res$hs_pop_tbl, "Observed_HS", "CI_L", "CI_U", "HS",
                             "HS per locus \u2014 populations block bootstrap CI")
        ggsave(file, plot = p, width = 12, height = 6, dpi = 300)
      }
    )
    ### FST ####
    output$download_fst_table <- downloadHandler(
      filename = function() paste0("fst_results_", Sys.Date(), ".csv"),
      content = function(file) {
        shiny::req(fst_boot_results())
        write.csv(fst_boot_results()$final_table, file, row.names = FALSE)
      }
    )
    output$download_fst_table_txt <- downloadHandler(
      filename = function() paste0("fst_results_", Sys.Date(), ".txt"),
      content = function(file) {
        shiny::req(fst_boot_results())
        write.table(fst_boot_results()$final_table, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    output$download_fst_plot <- downloadHandler(
      filename = function() paste0("fst_plot_", Sys.Date(), ".png"),
      content = function(file) {
        shiny::req(fst_boot_results())
        
        df <- fst_boot_results()$final_table
        shiny::req(is.data.frame(df), nrow(df) > 0)
        
        loci_only <- df$ID[df$ID != "Overall"]
        df$ID <- factor(df$ID, levels = c(sort(unique(loci_only)), "Overall"))
        
        df <- df %>% dplyr::mutate(Significant = !is.na(P_value) & P_value < 0.05)
        
        p <- ggplot(df, aes(x = ID, y = Observed_FST)) +
          geom_point(aes(shape = Significant), size = 3, color = "#3498db") +
          geom_errorbar(aes(ymin = CI_L, ymax = CI_U), width = 0.2, color = "#3498db") +
          # geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
          labs(
            title = "FST estimates with confidence intervals",
            x = "Locus",
            y = "FST estimate",
            shape = "p < 0.05"
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title  = element_text(face = "bold", hjust = 0.5)
          )
        
        ggsave(file, plot = p, width = 12, height = 6, dpi = 300)
      }
    )

    ## -- Testing outputs ------------------------------------------------------
    ### Testing local panmixia (FIS permutation) ----
    output$fis_pval_testing <- DT::renderDT({
      shiny::req(fis_boot_results())
      df <- fis_boot_results()$final_table
      shiny::validate(shiny::need(all(c("ID", "Observed_FIS", "P_value") %in% names(df)),
                                  "FIS results malformed."))
      df[, c("ID", "Observed_FIS", "P_value")]
    },
    options = list(pageLength = 25, scrollX = TRUE),
    rownames = FALSE,
    caption = "FIS permutation test \u2014 H0: local panmixia"
    )

    ### Testing global panmixia (FIT permutation) ----
    output$fit_pval_testing <- DT::renderDT({
      shiny::req(fit_boot_results())
      df <- fit_boot_results()$final_table
      shiny::validate(shiny::need(all(c("ID", "Observed_FIT", "P_value") %in% names(df)),
                                  "FIT results malformed."))
      df[, c("ID", "Observed_FIT", "P_value")]
    },
    options = list(pageLength = 25, scrollX = TRUE),
    rownames = FALSE,
    caption = "FIT permutation test \u2014 H0: global panmixia"
    )

    ### Testing subdivision (FST permutation) ----
    output$fst_pval_testing <- DT::renderDT({
      shiny::req(fst_boot_results())
      df <- fst_boot_results()$final_table
      shiny::validate(shiny::need(all(c("ID", "Observed_FST", "P_value") %in% names(df)),
                                  "FST results malformed."))
      df[, c("ID", "Observed_FST", "P_value")]
    },
    options = list(pageLength = 25, scrollX = TRUE),
    rownames = FALSE,
    caption = "FST permutation test \u2014 H0: no subdivision"
    )

  })
}