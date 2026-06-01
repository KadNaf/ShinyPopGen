# helper.R

# Shared CSS + JS injected by every general-stats module UI
gs_head <- function() {
  tagList(
    tags$style(HTML("
    .box-title { font-weight: bold; font-size: 16px; }
    .success-box { border-left: 4px solid #28a745; }
    .warning-box { border-left: 4px solid #ffc107; }
    .info-box { border-left: 4px solid #17a2b8; }
    .value-box { font-size: 24px; font-weight: bold; text-align: center; }
    .shiny-notification { position: fixed; top: 20px; right: 20px; }
    .section-title {
      color: #333a43;
      text-align: center;
      margin: 30px 0 20px 0;
      padding: 10px;
      border-bottom: 2px solid #CEB175;
    }
    .function-box { margin-bottom: 15px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
    .btn-function { width: 100%; margin-top: 10px; }
    h4, h5 { font-weight: bold; }
    .value-box .small-box h3 { font-size: 28px !important; font-weight: bold; }
    .value-box .small-box p  { font-size: 14px !important; margin-bottom: 0; }
    .progress-bar { border-radius: 10px; height: 20px; transition: width 0.5s ease; }
    .progress     { height: 20px; margin-bottom: 10px; }
    .status-indicator { transition: all 0.3s ease; }
    .status-indicator:hover { transform: translateY(-2px); box-shadow: 0 4px 8px rgba(0,0,0,0.1); }
    .small-box h3, .small-box p { color: #000000 !important; }
    div.dataTables_length { margin-bottom: 10px; }
    div.dt-buttons { margin-top: 6px; margin-bottom: 8px; }
    table.dataTable th, table.dataTable td { white-space: nowrap; }
  ")),
    tags$script(HTML("
    $(document).on('shown.bs.tab', 'a[data-toggle=\"tab\"]', function () {
      $($.fn.dataTable.tables(true)).DataTable().columns.adjust();
    });
    $(window).on('resize', function () {
      $($.fn.dataTable.tables(true)).DataTable().columns.adjust();
    });
  "))
  )
}

# Compact gradient banner replacing h2 section-title in each module
module_banner <- function(icon_name, title, subtitle, accent = "#6B64EF") {
  # Unique gradient ID per call: avoids cross-SVG id collision in the DOM
  gid <- paste0("spg-g", format(Sys.time(), "%H%M%S"), sample.int(9999, 1))
  shiny::div(
    class = "spg-module-banner",
    style = paste0("border-bottom-color:", accent, ";"),
    # Left: module icon + text
    shiny::div(class = "spg-banner-icon", shiny::icon(icon_name)),
    shiny::div(
      style = "flex:1; min-width:0;",
      shiny::tags$h2(title, class = "spg-banner-title"),
      shiny::tags$p(subtitle, class = "spg-banner-subtitle")
    ),
    # Right: ShinyPopGen brand SVG text + circular logo
    shiny::div(
      style = "flex-shrink:0; display:flex; flex-direction:row; align-items:center; gap:12px;",
      shiny::HTML(paste0(
        '<svg viewBox="0 0 220 115" height="88" xmlns="http://www.w3.org/2000/svg" aria-label="ShinyPopGen">',
        '<defs>',
        '<linearGradient id="', gid, '" x1="0" y1="0" x2="1" y2="0">',
        '<stop offset="0%" stop-color="#8F86FF"/>',
        '<stop offset="100%" stop-color="#5AA7FF"/>',
        '</linearGradient>',
        '</defs>',
        '<text x="2" y="40" fill="#F4F6FF" font-size="36" font-family="Inter,Segoe UI,Roboto,sans-serif" font-weight="300" letter-spacing="-0.5">Shiny</text>',
        '<text x="2" y="78" fill="url(#', gid, ')" font-size="38" font-family="Inter,Segoe UI,Roboto,sans-serif" font-weight="500" letter-spacing="-0.8">PopGen</text>',
        '<line x1="2" y1="88" x2="198" y2="88" stroke="#7074D8" stroke-width="1"/>',
        '<text x="2" y="106" fill="#A8ACF8" font-size="9" font-family="Inter,Segoe UI,Roboto,sans-serif" font-weight="400" letter-spacing="2.5">POPULATION GENETICS</text>',
        '</svg>'
      )),
      shiny::tags$img(
        src   = "spg_www/shinypopgen_logo.svg",
        height = "72px",
        alt   = "ShinyPopGen",
        style = "opacity:0.88; filter:drop-shadow(0 2px 10px rgba(0,0,0,0.5));"
      )
    )
  )
}

.time_it <- function(label, expr) {
  expr <- substitute(expr)
  eval(expr, envir = parent.frame())
}

# Import and format data (server_import_data tab) ====

# NULL-only coalesce (SAFE for any object)
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Scalar coalesce: NULL / length0 / scalar NA / scalar "" only
`%|||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) return(y)
  
  # only treat scalar NA / "" as missing
  if (length(x) == 1L) {
    if (is.na(x)) return(y)
    if (is.character(x) && !nzchar(x)) return(y)
  }
  
  x
}

.merge_paired_preview <- function(df, allele_cols) {
  # allele_cols: c("B12","B12_1","C07","C07_1",...)
  loci <- unique(sub("(_1|\\.[0-9]+)$", "", allele_cols))
  
  pick_b <- function(locus) {
    cands <- c(paste0(locus, "_1"), paste0(locus, ".", 1:9))
    hit <- cands[cands %in% names(df)]
    if (length(hit)) hit[1] else NA_character_
  }
  
  for (locus in loci) {
    a <- locus
    b <- pick_b(locus)
    if (!a %in% names(df) || is.na(b) || !b %in% names(df)) next
    
    df[[locus]] <- paste0(df[[a]], "/", df[[b]])
    df[[a]] <- paste0(df[[a]], "/", df[[b]])
    df[[b]] <- NULL
  }
  
  df
}

.get_locus_cols_from_marker_cols <- function(marker_cols) {
  marker_cols <- as.character(marker_cols)
  has_pairs <- any(grepl("(_1|\\.[0-9]+)$", marker_cols))
  if (!has_pairs) return(marker_cols)
  
  base <- unique(sub("(_1|\\.[0-9]+)$", "", marker_cols))
  
  # keep only loci that have an a1 column AND at least one a2 candidate
  keep <- base[base %in% marker_cols & vapply(base, function(b) {
    any(c(paste0(b, "_1"), paste0(b, ".", 1:9)) %in% marker_cols)
  }, logical(1))]
  
  keep
}

.build_hf_from_params <- function(con, rv, missing_code_raw, batch_size = 3000L) {
  req(con)
  
  # ---- marker_cols_raw (required for paired build)
  marker_json <- tryCatch(
    DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='marker_cols_raw'")$value[1],
    error = function(e) NA_character_
  )
  
  marker_cols_raw <- if (!is.na(marker_json) && nzchar(marker_json)) {
    tryCatch(jsonlite::fromJSON(marker_json), error = function(e) character(0))
  } else character(0)
  
  if (!length(marker_cols_raw)) stop("No marker_cols_raw found in DuckDB params.")
  
  # ---- base
  base <- tryCatch(
    DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='base'")$value[1],
    error = function(e) NA_character_
  )
  base <- suppressWarnings(as.integer(base))
  if (!is.finite(base) || base <= 0L) stop("Invalid base in DuckDB params.")
  
  # ---- genotype_format
  fmt <- .duckdb_get_param(con, "genotype_format", default = "auto")
  if (identical(fmt, "auto")) {
    # safest default: if raw has suffixes -> paired else string
    fmt <- if (any(grepl("(_1|\\.[0-9]+)$", marker_cols_raw))) "paired" else "string"
  }
  
  .time_it("HF build (chunked)", {
    .duckdb_build_hf_from_raw_chunked(
      con              = con,
      tbl_raw          = rv$tbl_raw,
      tbl_meta         = rv$tbl_meta,
      tbl_hf           = rv$tbl_hf,
      marker_cols      = marker_cols_raw,  # <--- IMPORTANT
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




.duckdb_set_params <- function(con, params, tbl_params = "params",
                               json_keys = c("marker_cols_raw", "locus_cols")){
  stopifnot(DBI::dbIsValid(con))
  stopifnot(is.list(params) || is.environment(params))
  
  # coerce to named list
  if (is.environment(params)) params <- as.list(params)
  if (length(params) == 0) return(invisible(TRUE))
  if (is.null(names(params)) || any(!nzchar(names(params)))) {
    stop("`.duckdb_set_params()` expects a *named* list.")
  }
  
  tbl_sql <- .sql_ident(tbl_params)
  
  # ensure params table exists
  DBI::dbExecute(con, sprintf("
    CREATE TABLE IF NOT EXISTS %s (
      key   VARCHAR PRIMARY KEY,
      value VARCHAR
    );
  ", tbl_sql))
  
  keys <- names(params)
  
  # serialize values (JSON for marker_cols; scalars as plain text; other vectors/lists as JSON)
  vals <- Map(function(k, x) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) return(NA_character_)
    
    # force JSON for selected keys (marker_cols)
    if (k %in% json_keys) {
      return(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"))
    }
    
    # scalar -> plain string
    if (length(x) == 1 && !is.list(x)) {
      return(as.character(x))
    }
    
    # non-scalar -> JSON (safe default)
    jsonlite::toJSON(x, auto_unbox = TRUE, null = "null")
  }, keys, unname(params))
  
  df <- data.frame(
    key   = as.character(keys),
    value = as.character(unlist(vals, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  
  # parameterised upsert via temporary table
  tmp <- paste0("tmp_params_", as.integer(stats::runif(1, 1, 1e9)))
  tmp_sql <- .sql_ident(tmp)
  
  DBI::dbWriteTable(con, tmp, df, temporary = TRUE, overwrite = TRUE)
  
  DBI::dbExecute(con, sprintf("
    INSERT INTO %s AS p
    SELECT key, value FROM %s
    ON CONFLICT(key) DO UPDATE SET value = excluded.value;
  ", tbl_sql, tmp_sql))
  
  invisible(TRUE)
}



infer_base_from_marker_strings <- function(marker_df, sep = "/", sample_n = 5000L, default_base = 1000L) {
  # marker_df: data.frame of character-ish genotype strings
  x <- unlist(marker_df, use.names = FALSE)
  x <- x[!is.na(x)]
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  if (!length(x)) return(as.integer(default_base))
  
  # sample to keep it cheap
  if (length(x) > sample_n) x <- sample(x, sample_n)
  
  # keep only strings with sep and digits
  x <- x[grepl("/", x, fixed = TRUE)]
  if (!length(x)) return(as.integer(default_base))
  
  a <- suppressWarnings(as.integer(sub("/.*$", "", x)))
  b <- suppressWarnings(as.integer(sub("^.*/", "", x)))
  v <- max(c(a, b), na.rm = TRUE)
  if (!is.finite(v) || v <= 0) return(as.integer(default_base))
  
  # base must be > max allele. your previous logic used 10^width, which is fine.
  width <- nchar(as.character(v))
  as.integer(10L ^ width)
}


parse_col_index_ranges <- function(x, n_max) {
  if (is.null(x) || !nzchar(trimws(x))) return(integer(0))
  
  x <- gsub("\\s+", "", x)
  chunks <- unlist(strsplit(x, ",", fixed = TRUE))
  out <- integer(0)
  
  for (ch in chunks) {
    if (!nzchar(ch)) next
    
    if (grepl("^[0-9]+([:-])[0-9]+$", ch)) {
      ab <- unlist(strsplit(ch, "[:-]"))
      a <- suppressWarnings(as.integer(ab[1]))
      b <- suppressWarnings(as.integer(ab[2]))
      if (!is.na(a) && !is.na(b)) out <- c(out, seq.int(min(a,b), max(a,b)))
    } else if (grepl("^[0-9]+$", ch)) {
      a <- suppressWarnings(as.integer(ch))
      if (!is.na(a)) out <- c(out, a)
    }
  }
  
  out <- unique(out)
  out <- out[out >= 1 & out <= n_max]
  sort(out)
}

.compress_idx_ranges <- function(idxs) {
  idxs <- sort(unique(as.integer(idxs)))
  idxs <- idxs[!is.na(idxs)]
  if (!length(idxs)) return("")
  runs <- split(idxs, cumsum(c(1, diff(idxs) != 1)))
  pieces <- vapply(runs, function(r) {
    if (length(r) == 1) as.character(r)
    else paste0(min(r), "-", max(r))
  }, character(1))
  paste(pieces, collapse = ",")
}


.duckdb_build_hf_from_raw_chunked <- function(
    con,
    tbl_raw,
    tbl_meta = NULL,
    tbl_hf = "hf",
    marker_cols,
    genotype_format = c("auto","paired","string"),
    base,
    missing_code_raw = NULL,
    missing_gt = 0L,
    batch_size = 10000L,
    attach_pop_code = TRUE
) {
  stopifnot(DBI::dbIsValid(con))
  stopifnot(DBI::dbExistsTable(con, tbl_raw))
  stopifnot(length(marker_cols) > 0)
  stopifnot(is.finite(base) && base > 0)
  
  genotype_format <- match.arg(genotype_format)
  
  # ---- create output table
  DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS %s;", .sql_ident(tbl_hf)))
  DBI::dbExecute(con, sprintf(
    "CREATE TABLE %s (indiv_id BIGINT, locus_id VARCHAR, gt INTEGER);",
    .sql_ident(tbl_hf)
  ))
  
  tbl_raw_sql <- .sql_ident(tbl_raw)
  tbl_hf_sql  <- .sql_ident(tbl_hf)
  
  # ---- explicit format control (no more relying on list_length())
  has_pairs <- switch(
    genotype_format,
    paired = TRUE,
    string = FALSE,
    auto   = any(grepl("(_1|\\.[0-9]+)$", marker_cols))
  )
  
  # treat these strings as missing (string-format only)
  miss_clause <- ""
  if (!is.null(missing_code_raw) && nzchar(trimws(as.character(missing_code_raw)))) {
    miss_sql <- DBI::dbQuoteString(con, as.character(missing_code_raw))
    miss_clause <- sprintf("WHEN gt_str = %s THEN %d", miss_sql, as.integer(missing_gt))
  }
  
  if (isTRUE(has_pairs)) {
    # =========================================================#
    # PAIRED FORMAT: two allele columns per locus (e.g. B12 + B12_1)
    # =========================================================#
    loci_all <- unique(sub("(_1|\\.[0-9]+)$", "", marker_cols))
    n_loci   <- length(loci_all)
    batches  <- split(seq_len(n_loci), ceiling(seq_len(n_loci) / batch_size))
    
    pick_b <- function(locus, marker_cols) {
      cands <- c(paste0(locus, "_1"), paste0(locus, ".", 1:9))
      hit <- cands[cands %in% marker_cols]
      if (length(hit) == 0) NA_character_ else hit[1]
    }
    
    for (idx in batches) {
      loci   <- loci_all[idx]
      n_batch <- length(loci)
      
      a_cols <- loci
      b_cols <- vapply(loci, pick_b, character(1), marker_cols = marker_cols)
      
      # guard
      missing_a <- setdiff(a_cols, marker_cols)
      missing_b <- setdiff(b_cols, marker_cols)
      if (length(missing_a) || length(missing_b) || anyNA(b_cols)) {
        stop(
          "Paired allele format detected, but missing expected columns:\n",
          if (length(missing_a)) paste0("  missing a1: ", paste(missing_a, collapse = ", "), "\n") else "",
          if (length(missing_b) || anyNA(b_cols)) paste0("  missing a2: ", paste(unique(c(missing_b, b_cols[is.na(b_cols)])), collapse = ", "), "\n") else ""
        )
      }
      
      a_ident <- vapply(a_cols, function(x) as.character(DBI::dbQuoteIdentifier(con, x)), character(1))
      b_ident <- vapply(b_cols, function(x) as.character(DBI::dbQuoteIdentifier(con, x)), character(1))
      
      loci_sql <- paste(vapply(loci, function(x) DBI::dbQuoteString(con, x), character(1)), collapse = ", ")
      
      list_a <- paste(sprintf("TRY_CAST(%s AS INTEGER)", a_ident), collapse = ", ")
      list_b <- paste(sprintf("TRY_CAST(%s AS INTEGER)", b_ident), collapse = ", ")
      
      sql <- sprintf(
        "
      INSERT INTO %s
      WITH t AS (
        SELECT rowid AS indiv_id
        FROM %s
      ),
      lists AS (
        SELECT
          t.indiv_id,
          [%s] AS locus_list,
          [%s] AS a1_list,
          [%s] AS a2_list
        FROM %s r
        JOIN t ON r.rowid = t.indiv_id
      ),
      u AS (
        SELECT
          indiv_id,
          list_extract(locus_list, i) AS locus_id,
          list_extract(a1_list,   i) AS a1,
          list_extract(a2_list,   i) AS a2
        FROM lists,
        range(1, %d + 1) r(i)
      )
      SELECT
        indiv_id,
        locus_id,
        CASE
          WHEN a1 IS NULL OR a2 IS NULL THEN %d
          WHEN a1 <= 0 OR a2 <= 0 THEN %d
          ELSE a1 * %d + a2
        END AS gt
      FROM u;
      ",
        tbl_hf_sql,
        tbl_raw_sql,
        loci_sql,
        list_a,
        list_b,
        tbl_raw_sql,
        as.integer(n_batch),
        as.integer(missing_gt),
        as.integer(missing_gt),
        as.integer(base)
      )
      
      DBI::dbExecute(con, sql)
    }
    
  } else {
    # =========================================================#
    # STRING FORMAT: one column per locus like "195/197"
    # =========================================================#
    
    n <- length(marker_cols)
    batches <- split(seq_len(n), ceiling(seq_len(n) / batch_size))
    
    for (idx in batches) {
      cols <- marker_cols[idx]
      
      cols_ident  <- vapply(cols, function(x) as.character(DBI::dbQuoteIdentifier(con, x)), character(1))
      select_cols <- paste(cols_ident, collapse = ", ")
      unpivot_in  <- paste(cols_ident, collapse = ", ")
      
      sql <- sprintf(
        "
        INSERT INTO %s
        WITH t AS (
          SELECT rowid AS indiv_id, %s
          FROM %s
        ),
        u AS (
          SELECT
            indiv_id,
            locus_id,
            gt_str,
            TRY_CAST(split_part(gt_str,'/',1) AS INTEGER) AS a1,
            TRY_CAST(split_part(gt_str,'/',2) AS INTEGER) AS a2
          FROM t
          UNPIVOT (gt_str FOR locus_id IN (%s))
        )
        SELECT
          indiv_id,
          locus_id,
          CASE
            WHEN gt_str IS NULL OR gt_str = '' THEN %d
            %s
            WHEN a1 IS NULL OR a2 IS NULL THEN %d
            WHEN a1 <= 0 OR a2 <= 0 THEN %d
            ELSE a1 * %d + a2
          END AS gt
        FROM u;
        ",
        tbl_hf_sql,
        select_cols,
        tbl_raw_sql,
        unpivot_in,
        as.integer(missing_gt),
        miss_clause,
        as.integer(missing_gt),
        as.integer(missing_gt),
        as.integer(base)
      )
      
      DBI::dbExecute(con, sql)
    }
  }
  
  # ---- attach pop_code (optional)
  if (isTRUE(attach_pop_code) &&
      !is.null(tbl_meta) &&
      DBI::dbExistsTable(con, tbl_meta) &&
      "Population" %in% DBI::dbListFields(con, tbl_meta)) {
    
    hf_cols <- DBI::dbListFields(con, tbl_hf)
    if (!("pop_code" %in% hf_cols)) {
      DBI::dbExecute(con, sprintf("ALTER TABLE %s ADD COLUMN pop_code INTEGER;", tbl_hf_sql))
    }
    
    DBI::dbExecute(con, sprintf("
      UPDATE %s h
      SET pop_code = d.dense_id
      FROM (
        SELECT individual AS indiv_id,
               dense_rank() OVER (ORDER BY Population) AS dense_id
        FROM %s
      ) d
      WHERE h.indiv_id = d.indiv_id;
    ", tbl_hf_sql, .sql_ident(tbl_meta)))
  }
  
  invisible(TRUE)
}



.duckdb_tune_for_big_import <- function(con) {
  # use all logical cores for import; Shiny sessions don't compete here
  n_threads <- max(1L, parallel::detectCores(logical = TRUE))

  DBI::dbExecute(con, sprintf("PRAGMA threads=%d;", n_threads))
  DBI::dbExecute(con, "PRAGMA enable_progress_bar=false;")
  DBI::dbExecute(con, "SET preserve_insertion_order=false;")
  DBI::dbExecute(con, "SET memory_limit='16GB';")
  # parallel CSV reader (DuckDB >= 0.8; silently ignored on older versions)
  tryCatch(
    DBI::dbExecute(con, "SET enable_parallel_csv_reader=true;"),
    error = function(e) invisible(NULL)
  )
  invisible(TRUE)
}

populate_manual_from_detection <- function(session, det, colnames_all) {
  req(session, det, colnames_all)
  
  cn <- setNames(colnames_all, colnames_all)
  
  pop_sel <- if (!is.null(det$population) && nzchar(det$population) && det$population %in% colnames_all) det$population else ""
  lat_sel <- if (!is.null(det$latitude)  && nzchar(det$latitude)  && det$latitude  %in% colnames_all) det$latitude  else ""
  lon_sel <- if (!is.null(det$longitude) && nzchar(det$longitude) && det$longitude %in% colnames_all) det$longitude else ""
  
  opts <- list(placeholder = "select", maxOptions = 1000)
  
  updateSelectizeInput(session, "pop_data",      choices = cn, selected = pop_sel, server = TRUE, options = opts)
  updateSelectizeInput(session, "latitude_data", choices = cn, selected = lat_sel, server = TRUE, options = opts)
  updateSelectizeInput(session, "longitude_data",choices = cn, selected = lon_sel, server = TRUE, options = opts)
  
  updateTextInput(session, "col_ranges_data", value = det$marker_range %||% "")
  
  if (!is.null(det$metadata_cols) && length(det$metadata_cols)) {
    meta_idx <- match(det$metadata_cols, colnames_all)
    meta_idx <- meta_idx[!is.na(meta_idx)]
    updateTextInput(session, "metadata_ranges", value = .compress_idx_ranges(meta_idx))
  } else {
    updateTextInput(session, "metadata_ranges", value = "")
  }
  
  invisible(TRUE)
}



reset_downstream_state <- function(rv) {
  rv$preview_raw  <- NULL
  rv$preview_meta <- NULL
  rv$det <- NULL
  rv$colnames_all <- NULL
  
  if (!is.null(rv$con)) {
    try(.duckdb_clear_params(rv$con), silent = TRUE)
  }
  
  if (!is.null(rv$populationsLL_grouped) && is.function(rv$populationsLL_grouped)) {
    rv$populationsLL_grouped(NULL)
  }
  
  # IMPORTANT: do NOT re-render the widget; just clear layers
  leaflet::leafletProxy("map") %>%
    leaflet::clearMarkers() %>%
    leaflet::clearShapes()
}

.duckdb_import_raw <- function(con, tbl_raw, file_path, sep, header) {
  stopifnot(!is.null(con), nzchar(tbl_raw))
  if (!nzchar(file_path) || !file.exists(file_path)) {
    stop("duckdb import: file_path missing or does not exist")
  }
  
  DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS %s;", .sql_ident(tbl_raw)))
  
  q_file <- DBI::dbQuoteString(con, normalizePath(file_path, winslash = "/"))
  q_sep  <- DBI::dbQuoteString(con, sep)
  hdr    <- if (isTRUE(header)) "true" else "false"
  
  sql <- sprintf(
    "CREATE TABLE %s AS
     SELECT * FROM read_csv_auto(%s,
        delim=%s,
        header=%s,
        all_varchar=true,
        ignore_errors=true,
        parallel=true
     );",
    .sql_ident(tbl_raw), q_file, q_sep, hdr
  )

  # parallel=true is a DuckDB >= 0.9 hint; fall back silently if unsupported
  tryCatch(
    DBI::dbExecute(con, sql),
    error = function(e) {
      sql_fallback <- sprintf(
        "CREATE TABLE %s AS
         SELECT * FROM read_csv_auto(%s,
            delim=%s,
            header=%s,
            all_varchar=true,
            ignore_errors=true
         );",
        .sql_ident(tbl_raw), q_file, q_sep, hdr
      )
      DBI::dbExecute(con, sql_fallback)
    }
  )
  TRUE
}

.duckdb_build_meta <- function(con, tbl_raw, tbl_meta, colnames_all, pop_data,
                               latitude_data, longitude_data, metadata_ranges) {
  n_all <- length(colnames_all)
  
  meta_idx   <- parse_col_index_ranges(metadata_ranges, n_max = n_all)
  meta_names <- if (length(meta_idx)) colnames_all[meta_idx] else character(0)
  
  meta_names <- unique(meta_names)
  meta_names <- setdiff(meta_names, c(pop_data, latitude_data, longitude_data))
  
  # quote column identifiers
  q <- function(x) paste0('"', gsub('"', '""', x), '"')
  
  sel <- c(sprintf('%s AS Population', q(pop_data)))
  
  has_lat <- !is.null(latitude_data) && nzchar(latitude_data) && latitude_data %in% colnames_all
  has_lon <- !is.null(longitude_data) && nzchar(longitude_data) && longitude_data %in% colnames_all
  
  # robust casts: junk -> NULL instead of error
  if (has_lat) sel <- c(sel, sprintf('TRY_CAST(%s AS DOUBLE) AS Latitude',  q(latitude_data)))
  if (has_lon) sel <- c(sel, sprintf('TRY_CAST(%s AS DOUBLE) AS Longitude', q(longitude_data)))
  
  if (length(meta_names)) sel <- c(sel, sprintf('%s', q(meta_names)))
  
  sql <- sprintf(
    "CREATE OR REPLACE TABLE %s AS
     SELECT %s
     FROM %s;",
    .sql_ident(tbl_meta),
    paste(sel, collapse = ", "),
    .sql_ident(tbl_raw)
  )
  
  DBI::dbExecute(con, sql)
  
  invisible(list(meta_names = meta_names, has_gps = has_lat && has_lon))
}

.duckdb_get_params <- function(con) {
  if (is.null(con) || !DBI::dbExistsTable(con, "params")) return(list())
  
  df <- DBI::dbGetQuery(con, "SELECT key, value FROM params;")
  if (!nrow(df)) return(list())
  
  out <- as.list(stats::setNames(df$value, df$key))
  
  # auto-decode JSON arrays that we stored for vectors
  for (k in names(out)) {
    v <- out[[k]]
    if (!is.character(v) || length(v) != 1) next
    
    vv <- trimws(v)
    
    # ---- FIX: guard NA / empty
    if (is.na(vv) || !nzchar(vv)) next
    
    # ---- FIX: safer array detection
    if (startsWith(vv, "[") && endsWith(vv, "]")) {
      out[[k]] <- tryCatch(jsonlite::fromJSON(vv), error = function(e) v)
    }
  }
  
  out
}


.duckdb_clear_params <- function(con) {
  if (is.null(con)) return(invisible(FALSE))
  ok <- tryCatch(DBI::dbExistsTable(con, "params"), error = function(e) FALSE)
  if (isTRUE(ok)) DBI::dbExecute(con, "DELETE FROM params;")
  invisible(TRUE)
}

.duckdb_get_param <- function(con, key, default = NA_character_) {
  if (!DBI::dbExistsTable(con, "params")) return(default)
  x <- tryCatch(
    DBI::dbGetQuery(con, "SELECT value FROM params WHERE key = ?", params = list(key))$value[1],
    error = function(e) NA_character_
  )
  if (is.na(x) || !nzchar(x)) default else x
}

.duckdb_get_param_json <- function(con, key) {
  x <- .duckdb_get_param(con, key, default = NA_character_)
  if (is.na(x) || !nzchar(x)) return(character(0))
  tryCatch(jsonlite::fromJSON(x), error = function(e) character(0))
}

.find_longest_contiguous_block <- function(idxs) {
  if (!length(idxs)) return(integer(0))
  idxs <- sort(unique(idxs))
  runs <- split(idxs, cumsum(c(1, diff(idxs) != 1)))
  runs[[which.max(lengths(runs))]]
}

.normalize_name <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(x))
}

# very small typo-tolerant match: allow one of these to match after normalization
.is_popcode_name <- function(nm) {
  nm0 <- .normalize_name(nm)
  nm0 %in% c(
    "pop", "popid", "popcode", "populationid", "population_id",
    "siteid", "localityid", "locationid"
  )
}

.is_mostly_numeric <- function(x, min_frac = 0.8) {
  xx <- suppressWarnings(as.numeric(as.character(x)))
  mean(!is.na(xx)) >= min_frac
}

.is_lat_name <- function(nm) {
  nm <- .normalize_name(nm)
  nm %in% c("lat","latitude","latdd","latitudedd","decimallatitude","y","ywgs84","coordy","northing")
}

.is_lon_name <- function(nm) {
  nm <- .normalize_name(nm)
  nm %in% c("lon","lng","long","longitude","londd","longitudedd","decimallongitude","x","xwgs84","coordx","easting")
}

# numeric plausibility check
.is_coord_numeric <- function(x, kind = c("lat","lon"), min_frac = 0.8, require_decimal = FALSE) {
  kind <- match.arg(kind)
  
  s <- as.character(x)
  s <- trimws(s)
  
  xx <- suppressWarnings(as.numeric(s))
  ok <- !is.na(xx)
  if (mean(ok) < min_frac) return(FALSE)
  
  # optional: reject integer-only columns (HostInd etc.)
  if (require_decimal) {
    dec_frac <- mean(grepl("\\.", s[ok]))
    if (dec_frac < min_frac) return(FALSE)
  }
  
  if (kind == "lat") return(mean(xx[ok] >= -90 & xx[ok] <= 90) >= min_frac)
  if (kind == "lon") return(mean(xx[ok] >= -180 & xx[ok] <= 180) >= min_frac)
  FALSE
}

# genotype-like content detection
.is_genotype_col <- function(x, missing_info, min_frac = 0.5) {
  s <- as.character(x)
  s <- s[!is.na(s)]
  if (!length(s)) return(FALSE)
  
  s <- trimws(s)
  s <- s[s != ""]
  if (!length(s)) return(FALSE)
  
  miss_set <- missing_info$miss_set
  if (!is.null(miss_set) && length(miss_set)) {
    s <- s[!tolower(s) %in% tolower(miss_set)]
  }
  if (!length(s)) return(FALSE)
  
  pat_sep <- "^[0-9]+\\s*[/\\-_\\|:]\\s*[0-9]+$"
  if (mean(grepl(pat_sep, s)) >= min_frac) return(TRUE)
  
  pat_int <- "^[0-9]+$"
  frac_int <- mean(grepl(pat_int, s))
  
  if (frac_int >= 0.9) {
    xi <- suppressWarnings(as.integer(s[grepl(pat_int, s)]))
    xi <- xi[!is.na(xi)]
    if (!length(xi)) return(FALSE)
    
    # exclude obvious IDs: too many unique values (e.g., HostInd 1..n)
    uniq_ratio <- length(unique(xi)) / length(xi)
    if (uniq_ratio > 0.8) return(FALSE)
    
    # keep your existing “not GPS/UTM-like” guardrails
    if (max(xi, na.rm = TRUE) > 5000) return(FALSE)
    if (max(nchar(as.character(xi)), na.rm = TRUE) > 4) return(FALSE)
    
    return(TRUE)
  }
  FALSE
}

.is_idlike_name <- function(nm) {
  nm0 <- .normalize_name(nm)
  grepl(
    "(^id$|^ind$|^indiv|^individual|hostind|host|specimen|sample|code|barcode|tag|uid|uuid|run|rep|plate|well|popcode|popid)",
    nm0
  )
}

.is_idlike_content <- function(x, min_unique_ratio = 0.8, min_int_frac = 0.9) {
  s <- as.character(x)
  s <- trimws(s)
  s <- s[!is.na(s) & s != ""]
  if (!length(s)) return(FALSE)
  
  # only consider if mostly integer-like
  if (mean(grepl("^[0-9]+$", s)) < min_int_frac) return(FALSE)
  
  xi <- suppressWarnings(as.integer(s))
  xi <- xi[!is.na(xi)]
  if (!length(xi)) return(FALSE)
  
  uniq_ratio <- length(unique(xi)) / length(xi)
  uniq_ratio >= min_unique_ratio
}

.is_pop_name <- function(x) {
  if (length(x) != 1L || is.na(x)) return(FALSE)
  nm <- tolower(trimws(as.character(x)))
  
  # common variants in your app / popgen datasets
  patterns <- c(
    "^pop$",
    "^population$",
    "^populations$",
    "^pop_id$",
    "^popid$",
    "^deme$",
    "^site$",
    "^locality$",
    "^location$",
    "^sampling_site$",
    "^samplinglocation$",
    "^group$",
    "^cluster$",
    "^subpop$",
    "^subpopulation$",
    "^strata$",
    "^stratum$",
    "^region$",
    "^zone$"
  )
  
  any(vapply(patterns, function(p) grepl(p, nm, perl = TRUE), logical(1)))
}

.sql_ident <- function(x) {
  x <- gsub('"', '""', x, fixed = TRUE)
  paste0('"', x, '"')
}

update_metadata_choices <- function(session, colnames_all) {
  req(session)
  if (is.null(colnames_all) || !length(colnames_all)) return(invisible(FALSE))
  
  cn   <- setNames(colnames_all, colnames_all)
  opts <- list(placeholder = "select", maxOptions = 1000)
  
  updateSelectizeInput(session, "pop_data",
                       choices = cn, selected = "", server = TRUE, options = opts)
  
  updateSelectizeInput(session, "latitude_data",
                       choices = cn, selected = "", server = TRUE, options = opts)
  
  updateSelectizeInput(session, "longitude_data",
                       choices = cn, selected = "", server = TRUE, options = opts)
  
  updateTextInput(session, "metadata_ranges", value = "")
  
  invisible(TRUE)
}


detect_columns_auto <- function(df, missing_info) {
  stopifnot(is.data.frame(df))
  cn <- names(df)
  
  # --- population
  pop_candidates <- which(vapply(cn, .is_pop_name, logical(1)))
  
  pop_col <- NA_character_
  if (length(pop_candidates)) {
    nonnum <- pop_candidates[!vapply(pop_candidates, function(i) .is_mostly_numeric(df[[i]]), logical(1))]
    if (length(nonnum)) {
      pop_col <- cn[nonnum[1]]
    } else {
      pop_col <- cn[pop_candidates[1]]
    }
  }
  # --- GPS
  lat_candidates <- which(vapply(cn, .is_lat_name, logical(1)))
  lon_candidates <- which(vapply(cn, .is_lon_name, logical(1)))
  
  lat_col <- NA_character_
  lon_col <- NA_character_
  
  # 1) Name-based selection (wins)
  if (length(lat_candidates)) {
    for (i in lat_candidates) {
      if (.is_coord_numeric(df[[i]], "lat")) { lat_col <- cn[i]; break }
    }
  }
  if (length(lon_candidates)) {
    for (i in lon_candidates) {
      if (.is_coord_numeric(df[[i]], "lon")) { lon_col <- cn[i]; break }
    }
  }
  
  # 2) Fallback scan ONLY if still missing after name-based
  if (is.na(lat_col)) {
    for (i in seq_along(cn)) {
      if (.is_coord_numeric(df[[i]], "lat", require_decimal = TRUE)) { lat_col <- cn[i]; break }
    }
  }
  if (is.na(lon_col)) {
    for (i in seq_along(cn)) {
      if (.is_coord_numeric(df[[i]], "lon", require_decimal = TRUE)) { lon_col <- cn[i]; break }
    }
  }
  
  # --- markers (content-driven)
  marker_flags <- vapply(df, .is_genotype_col, logical(1), missing_info = missing_info)
  marker_idxs <- which(marker_flags)
  
  # do not allow population/GPS as markers
  marker_idxs <- setdiff(marker_idxs, match(pop_col, cn))
  marker_idxs <- setdiff(marker_idxs, match(lat_col, cn))
  marker_idxs <- setdiff(marker_idxs, match(lon_col, cn))
  
  # --- force ID/code columns out of markers (send to metadata instead)
  idlike_by_name    <- which(vapply(cn, .is_idlike_name, logical(1)))
  idlike_by_content <- which(vapply(df, .is_idlike_content, logical(1)))
  
  idlike_idxs <- sort(unique(c(idlike_by_name, idlike_by_content)))
  
  # never remove the chosen population column, even if it looks "id-like"
  idlike_idxs <- setdiff(idlike_idxs, match(pop_col, cn))
  
  marker_idxs <- setdiff(marker_idxs, idlike_idxs)
  
  
  
  # --- force population-code columns (numeric) out of markers
  popcode_idxs <- which(vapply(cn, .is_popcode_name, logical(1)))
  
  # keep only those that are mostly numeric
  popcode_idxs <- popcode_idxs[vapply(popcode_idxs, function(i) .is_mostly_numeric(df[[i]]), logical(1))]
  
  # never remove the chosen population column
  popcode_idxs <- setdiff(popcode_idxs, match(pop_col, cn))
  
  marker_idxs <- setdiff(marker_idxs, popcode_idxs)
  
  marker_cols <- cn[marker_idxs]
  
  # --- metadata
  # GPS columns (lat_col, lon_col) are kept in meta_cols so they appear
  # in the metadata range display as a contiguous block. .duckdb_build_meta()
  # already excludes them from the generic metadata SELECT to avoid duplication.
  special <- c(pop_col, marker_cols)
  special <- special[!is.na(special)]
  meta_cols <- setdiff(cn, special)
  
  block <- .find_longest_contiguous_block(marker_idxs)
  marker_range <- if (length(block)) sprintf("%d:%d", min(block), max(block)) else NA_character_
  
  list(
    population    = pop_col,
    latitude      = lat_col,
    longitude     = lon_col,
    marker_cols   = marker_cols,
    marker_range  = marker_range,
    metadata_cols = meta_cols
  )
}

normalize_missing_code <- function(x) {
  # Default
  if (is.null(x) || length(x) == 0) x <- ""
  s <- trimws(as.character(x)[1])
  
  # Empty -> default
  if (identical(s, "")) s <- "0"
  
  # If user already provided a diploid-like "A/B" (or "A-B", "A_B"), keep it and canonicalise separator
  if (grepl("[/_-]", s)) {
    s2 <- gsub("[-_]", "/", s)
    parts <- unlist(strsplit(s2, "/", fixed = TRUE))
    parts <- trimws(parts)
    parts <- parts[parts != ""]
    if (length(parts) == 1) parts <- c(parts, parts)
    if (length(parts) >= 2) {
      a <- parts[1]; b <- parts[2]
      miss_gt <- paste0(a, "/", b)
    } else {
      miss_gt <- "0/0"
    }
  } else {
    # Single token -> treat as homozygous missing
    miss_gt <- paste0(s, "/", s)
  }
  
  # Build a robust missing set
  # - include both "/" and "-" and "_" versions
  # - include the raw input (as typed) and trimmed variants
  miss_set <- unique(c(
    miss_gt,
    gsub("/", "-", miss_gt),
    gsub("/", "_", miss_gt),
    s,
    gsub("[-_]", "/", s)
  ))
  
  # Also treat common NA-like strings as missing (user may type NA or leave actual NA)
  miss_set <- unique(c(miss_set, "NA", "NaN", "N/A", "na", "nan", "n/a"))
  
  list(miss_gt = miss_gt, miss_set = miss_set)
}

# ========================= server_general_stats ===============================
.step <- function(label, expr) {
  force(expr)
}

# Convenience wrapper: pick parallel if available
.pkg_has_fun <- function(fname) {
  exists(fname, mode = "function", inherits = TRUE)
}

.normalize_threads <- function(n_threads) {
  nt <- suppressWarnings(as.integer(n_threads))
  if (!is.finite(nt) || is.na(nt)) nt <- 1L
  max(1L, nt)
}

# Convenience wrapper: pick parallel if available + self-report
batch_permute_wc84_fst_auto <- function(dat,
                                        pop_col_1based = 1,
                                        missing_code = 0,
                                        base = 1000,
                                        B = 999,
                                        n_threads = 1,
                                        seed = 1,
                                        pval_method = "two_sided_abs",
                                        perm_scheme = "within_pop_alleles",
                                        debug = FALSE) {
  
  # normalize threads
  nt <- suppressWarnings(as.integer(n_threads))
  if (!is.finite(nt) || is.na(nt) || nt < 1L) nt <- 1L
  
  has_cpp <- exists("batch_permute_wc84_fst_parallel", mode = "function")
  if (!has_cpp) stop("Missing C++ symbol: batch_permute_wc84_fst_parallel()")
  
  use_par <- (nt > 1L)
  if (isTRUE(debug)) {
    msg <- if (use_par) "[perm FST] using PARALLEL backend" else "[perm FST] using 1-thread backend"
    message(sprintf("%s (requested n_threads=%d)", msg, nt))
  }
  
  res <- batch_permute_wc84_fst_parallel(
    dat            = dat,
    pop_col_1based = pop_col_1based,
    missing_code   = missing_code,
    base           = base,
    B              = B,
    n_threads      = nt,
    seed           = seed,
    pval_method    = pval_method,
    perm_scheme    = perm_scheme
  )
  
  attr(res, "parallel") <- list(
    used_parallel     = isTRUE(use_par),
    requested_threads = nt,
    backend           = "cpp_parallel",
    parallel_symbol   = "batch_permute_wc84_fst_parallel"
  )
  res
}

boot_indiv_wc84_fst_auto <- function(mat,
                                     pop_col_1based = 1,
                                     missing_code = 0,
                                     base = 1000,
                                     B = 1000,
                                     n_threads = 1,
                                     seed = 1,
                                     debug = FALSE) {
  nt <- .normalize_threads(n_threads)
  
  use_par <- (nt > 1L) && .pkg_has_fun("boot_indiv_wc84_fst_parallel")
  if (isTRUE(debug)) {
    msg <- if (use_par) "[boot indiv FST] using PARALLEL backend" else "[boot indiv FST] using SERIAL backend"
    message(sprintf("%s (requested n_threads=%d)", msg, nt))
  }
  
  res <- if (use_par) {
    boot_indiv_wc84_fst_parallel(
      mat,
      pop_col_1based = pop_col_1based,
      missing_code   = missing_code,
      base           = base,
      B              = B,
      n_threads      = nt,
      seed           = seed
    )
  } else {
    boot_indiv_wc84_fst(
      mat,
      pop_col_1based = pop_col_1based,
      missing_code   = missing_code,
      base           = base,
      B              = B
    )
  }
  
  attr(res, "parallel") <- list(
    used_parallel     = isTRUE(use_par),
    requested_threads = nt,
    backend           = if (use_par) "cpp_parallel" else "cpp_serial",
    parallel_symbol   = if (use_par) "boot_indiv_wc84_fst_parallel" else NA_character_
  )
  res
}

boot_popblock_wc84_fst_auto <- function(mat,
                                        pop_col_1based = 1,
                                        missing_code = 0,
                                        base = 1000,
                                        B = 1000,
                                        n_threads = 1,
                                        seed = 1,
                                        debug = FALSE) {
  nt <- .normalize_threads(n_threads)
  
  has_cpp <- .pkg_has_fun("boot_popblock_wc84_parallel")
  if (!has_cpp) stop("Missing C++ symbol: boot_popblock_wc84_parallel()")
  
  use_par <- (nt > 1L)
  if (isTRUE(debug)) {
    msg <- if (use_par) "[boot popblock] using PARALLEL backend" else "[boot popblock] using 1-thread backend"
    message(sprintf("%s (requested n_threads=%d)", msg, nt))
  }
  
  res <- boot_popblock_wc84_parallel(
    mat            = mat,
    pop_col_1based = pop_col_1based,
    missing_code   = missing_code,
    base           = base,
    B              = B,
    n_threads      = nt,
    seed           = seed
  )
  
  attr(res, "parallel") <- list(
    used_parallel     = isTRUE(use_par),
    requested_threads = nt,
    backend           = "cpp_parallel",
    parallel_symbol   = "boot_popblock_wc84_parallel"
  )
  res
}

.hash_key <- function(x) {
  # Stable hash; digest is common in Shiny stacks.
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(x, algo = "xxhash64"))
  }
  # fallback: not cryptographic, but stable-ish across session
  paste0(sprintf("%08x", as.integer(abs(stats::runif(1) * 2^31))), collapse = "")
}

.clamp01 <- function(x) pmin(pmax(x, 0), 1)

.compute_basic_stats <- function(mat, base, k) {
  stopifnot(is.matrix(mat), ncol(mat) >= 2L, base > 1L)

  loci <- colnames(mat)[-1L]

  # --- WC84 stats (per locus, Weir & Cockerham 1984) ---
  # observed_wc84_stats_cpp returns:
  #   FST, FIT, FIS  : WC84 F-statistics
  #   HS             : n_i-weighted within-population gene diversity (WC84)
  #   HT             : total gene diversity from pooled allele counts (WC84)
  st <- observed_wc84_stats_cpp(
    dat            = mat,
    pop_col_1based = 1L,
    missing_code   = 0L,
    base           = base
  )

  FST_wc <- as.numeric(st$FST)
  FIT_wc <- as.numeric(st$FIT)
  FIS_wc <- as.numeric(st$FIS)

  # WC84's own HS and HT — n_i-weighted, consistent with the WC84 F-statistic estimators.
  # These are used exclusively for Fst' (W&C) = FST_WC / FST_max_WC (Meirmans 2006).
  Hs_wc <- as.numeric(st$HS)   # n_i-weighted HS
  Ht_wc <- as.numeric(st$HT)   # pooled-count HT

  # --- Nei heterozygosity estimators (Ho / Hs / Ht, Nei 1977/1987) ---
  # These use the unweighted (equal-population) convention and are kept
  # for Ho, the displayed Hs/Ht columns, GST and GST''.
  # NOTE: Nei's HS and WC84's HS diverge for unequal sample sizes.
  nei <- nei_het_stats_cpp(
    dat            = mat,
    pop_col_1based = 1L,
    missing_code   = 0L,
    base           = base
  )

  Ho   <- as.numeric(nei$Ho)
  Hs   <- as.numeric(nei$Hs)   # Nei's HS — used for GST / GST'' / display
  Ht   <- as.numeric(nei$Ht)   # Nei's HT — used for GST / GST'' / display
  k_eff <- if (!is.null(nei$k_eff)) as.integer(nei$k_eff) else rep.int(as.integer(k), length(Ht))

  # Nei Fis/Fst from heterozygosities
  Fis_nei <- rep(NA_real_, length(Hs))

  ok_fis  <- is.finite(Ho) & is.finite(Hs) & Hs > 0
  Fis_nei[ok_fis] <- 1 - Ho[ok_fis] / Hs[ok_fis]

  Fst_nei <- rep(NA_real_, length(Ht))
  ok_fst  <- is.finite(Ht) & is.finite(Hs) & Ht > 0
  Fst_nei[ok_fst] <- (Ht[ok_fst] - Hs[ok_fst]) / Ht[ok_fst]

  # GST (Nei): 1 - Hs/Ht
  GST <- rep(NA_real_, length(Ht))
  ok_gst <- is.finite(Hs) & is.finite(Ht) & Ht > 0
  GST[ok_gst] <- 1 - Hs[ok_gst] / Ht[ok_gst]

  # --- Meirmans recoding for FSTmax ---------------------------------
  pop <- as.integer(mat[, 1])
  loci_idx <- seq.int(2L, ncol(mat))

  maps <- vector("list", length(loci_idx))
  names(maps) <- colnames(mat)[loci_idx]
  max_new_allele <- 0L

  # Pass 1: build allele maps per locus x population
  for (jj in seq_along(loci_idx)) {
    j <- loci_idx[jj]

    g <- as.integer(mat[, j])

    ok_gt <- is.finite(g) & g != 0L & g > 0L
    a1_all <- g %/% base
    a2_all <- g %%  base

    # valid diploid packed genotypes only
    ok <- ok_gt & a1_all > 0L & a2_all > 0L

    map_j <- list()

    if (any(ok)) {
      a1 <- a1_all[ok]
      a2 <- a2_all[ok]
      p  <- pop[ok]

      next_id <- 1L
      pop_ids <- sort(unique(p))

      for (pp in pop_ids) {
        idxp <- (p == pp)
        alleles_pp <- sort(unique(c(a1[idxp], a2[idxp])))

        if (length(alleles_pp) == 0L) next

        new_ids <- as.integer(seq.int(next_id, length.out = length(alleles_pp)))
        names(new_ids) <- as.character(alleles_pp)

        map_j[[as.character(pp)]] <- new_ids
        next_id <- next_id + length(alleles_pp)
      }

      max_new_allele <- max(max_new_allele, next_id - 1L)
    }

    maps[[jj]] <- map_j
  }

  new_base <- as.integer(max(2L, max_new_allele + 1L))

  # Pass 2: apply recoding
  mat_recoded <- mat
  storage.mode(mat_recoded) <- "integer"

  for (jj in seq_along(loci_idx)) {
    j <- loci_idx[jj]

    g <- as.integer(mat[, j])

    ok_gt <- is.finite(g) & g != 0L & g > 0L
    a1 <- g %/% base
    a2 <- g %%  base

    # valid diploid packed genotypes only
    ok <- ok_gt & a1 > 0L & a2 > 0L

    g_new <- rep.int(0L, length(g))

    if (any(ok)) {
      pop_ids <- sort(unique(pop[ok]))

      for (pp in pop_ids) {
        idxp <- which(ok & pop == pp)
        mp <- maps[[jj]][[as.character(pp)]]
        if (is.null(mp) || length(idxp) == 0L) next

        a1_new <- unname(mp[as.character(a1[idxp])])
        a2_new <- unname(mp[as.character(a2[idxp])])

        if (anyNA(a1_new) || anyNA(a2_new)) {
          bad_a1 <- unique(a1[idxp][is.na(a1_new)])
          bad_a2 <- unique(a2[idxp][is.na(a2_new)])
          stop(
            paste0(
              "Internal recoding error after filtering invalid packed genotypes. ",
              "locus=", colnames(mat)[j],
              ", pop=", pp,
              ", missing map a1={", paste(bad_a1, collapse = ","), "}",
              ", missing map a2={", paste(bad_a2, collapse = ","), "}"
            )
          )
        }

        lo <- pmin(a1_new, a2_new)
        hi <- pmax(a1_new, a2_new)

        g_new[idxp] <- as.integer(lo * new_base + hi)
      }
    }

    # invalid packed genotypes stay as missing (0)
    mat_recoded[, j] <- g_new
  }

  # ----------------------------------------------------------------
  # FST_max (WC84) on Meirmans-recoded data
  # Meirmans (2006): recode alleles so every population has a
  # completely private set of alleles (no sharing across pops).
  # WC84 FST on the recoded matrix = maximum attainable FST
  # given the observed within-population diversity (HS).
  # ----------------------------------------------------------------
  st_max <- observed_wc84_stats_cpp(
    dat            = mat_recoded,
    pop_col_1based = 1L,
    missing_code   = 0L,
    base           = new_base          # <-- recoded base, NOT original base
  )

  Fst_max <- as.numeric(st_max$FST)
  if (!is.null(st_max$locus_names) && length(st_max$locus_names) == length(Fst_max)) {
    names(Fst_max) <- as.character(st_max$locus_names)
    Fst_max <- Fst_max[match(loci, names(Fst_max))]
  }

  # ----------------------------------------------------------------
  # F'ST (W&C, Meirmans 2006) — empirical standardisation
  #
  #   F'ST(Meirmans) = FST_WC84(obs) / FST_WC84(Meirmans-recoded)
  #
  # FST_max is obtained by running WC84 on the recoded matrix above.
  # ----------------------------------------------------------------
  Fst_prime_meirmans <- rep(NA_real_, length(FST_wc))
  ok_fp <- is.finite(FST_wc) & is.finite(Fst_max) & Fst_max > 0
  Fst_prime_meirmans[ok_fp] <- FST_wc[ok_fp] / Fst_max[ok_fp]
  # ----------------------------------------------------------------
  # F'ST (Hedrick) — Thierry/FSTAT-compatible implementation
  #
  #   FST'(Hedrick) = FST / (1 - Hs)
  #
  # using WC84 FST in the numerator and Nei's unweighted Hs in the
  # denominator.
  # ----------------------------------------------------------------
  Fst_prime_hedrick <- rep(NA_real_, length(FST_wc))
  ok_fph <- is.finite(FST_wc) & is.finite(Hs) & (1 - Hs) > 0
  Fst_prime_hedrick[ok_fph] <- FST_wc[ok_fph] / (1 - Hs[ok_fph])


  # GST'' (FSTAT): k(Ht-Hs)/((kHt-Hs)(1-Hs))
  GST2 <- rep(NA_real_, length(Ht))
  ok_gst2 <- is.finite(Hs) & is.finite(Ht) & is.finite(k_eff) & k_eff > 1 &
    (1 - Hs) != 0 & (k_eff * Ht - Hs) != 0
  GST2[ok_gst2] <- (k_eff[ok_gst2] * (Ht[ok_gst2] - Hs[ok_gst2])) /
    ((k_eff[ok_gst2] * Ht[ok_gst2] - Hs[ok_gst2]) * (1 - Hs[ok_gst2]))

  out <- data.frame(
    ID                   = loci,
    Ho                   = Ho,
    Hs                   = Hs,
    Ht                   = Ht,
    `Fit (W&C)`          = .clamp01(FIT_wc),
    `Fis (W&C)`          = .clamp01(FIS_wc),
    `Fst (W&C)`          = .clamp01(FST_wc),
    `Fis (Nei)`          = .clamp01(Fis_nei),
    `Fst (Nei)`          = .clamp01(Fst_nei),
    `Fst-max (Meirmans)` = .clamp01(Fst_max),
    `Fst' (Meirmans)`    = .clamp01(Fst_prime_meirmans),
    `Fst' (Hedrick)`     = .clamp01(Fst_prime_hedrick),
    GST                  = .clamp01(GST),
    `GST''`              = .clamp01(GST2),
    stringsAsFactors = FALSE
  )
  names(out) <- c(
    "ID",
    "Ho",
    "Hs",
    "Ht",
    "Fit (W&C)",
    "Fis (W&C)",
    "Fst (W&C)",
    "Fis (Nei)",
    "Fst (Nei)",
    "Fst-max (Meirmans)",
    "Fst' (Meirmans)",
    "Fst' (Hedrick)",
    "GST",
    "GST''"
  )

  # --- WC84 overall ratio-of-sums ---
  FST_wc_overall <- .clamp01(as.numeric(st$FST_overall_ratio_of_sums))
  FIT_wc_overall <- .clamp01(as.numeric(st$FIT_overall_ratio_of_sums))
  FIS_wc_overall <- .clamp01(as.numeric(st$FIS_overall_ratio_of_sums))

  # --- Overall H (Nei) — for display columns + GST / GST'' ---
  Ho_overall    <- mean(Ho,   na.rm = TRUE)
  Hs_overall    <- mean(Hs,   na.rm = TRUE)   # Nei's HS
  Ht_overall    <- mean(Ht,   na.rm = TRUE)   # Nei's HT


  # Use a global k for overall derived metrics
  k_global <- as.integer(k)

  # Nei overall
  FisN_overall <- if (
    is.finite(Ho_overall) && is.finite(Hs_overall) && Hs_overall > 0
  ) {
    1 - Ho_overall / Hs_overall
  } else {
    NA_real_
  }

  FstN_overall <- if (
    is.finite(Ht_overall) && is.finite(Hs_overall) && Ht_overall > 0
  ) {
    (Ht_overall - Hs_overall) / Ht_overall
  } else {
    NA_real_
  }

  # Recompute overall derived metrics from global Ho/Hs/Ht
  GST_overall <- if (
    is.finite(Hs_overall) && is.finite(Ht_overall) && Ht_overall > 0
  ) {
    1 - Hs_overall / Ht_overall
  } else {
    NA_real_
  }

  Fst_max_overall <- as.numeric(st_max$FST_overall_ratio_of_sums)

  Fst_prime_meirmans_overall <- if (
    is.finite(FST_wc_overall) &&
    is.finite(Fst_max_overall) &&
    Fst_max_overall > 0
  ) {
    FST_wc_overall / Fst_max_overall
  } else {
    NA_real_
  }
  
  Fst_prime_hedrick_overall <- if (
    is.finite(FST_wc_overall) &&
    is.finite(Hs_overall) &&
    (1 - Hs_overall) > 0
  ) {
    FST_wc_overall / (1 - Hs_overall)
  } else {
    NA_real_
  }
  
  GST2_overall <- if (
    is.finite(Hs_overall) &&
    is.finite(Ht_overall) &&
    is.finite(k_global) &&
    k_global > 1 &&
    (1 - Hs_overall) != 0 &&
    (k_global * Ht_overall - Hs_overall) != 0
  ) {
    (k_global * (Ht_overall - Hs_overall)) /
      ((k_global * Ht_overall - Hs_overall) * (1 - Hs_overall))
  } else {
    NA_real_
  }

  overall_row <- data.frame(
    ID                   = "Overall",
    Ho                   = Ho_overall,
    Hs                   = Hs_overall,
    Ht                   = Ht_overall,
    `Fit (W&C)`          = FIT_wc_overall,
    `Fis (W&C)`          = FIS_wc_overall,
    `Fst (W&C)`          = FST_wc_overall,
    `Fis (Nei)`          = .clamp01(FisN_overall),
    `Fst (Nei)`          = .clamp01(FstN_overall),
    `Fst-max (Meirmans)` = .clamp01(Fst_max_overall),
    `Fst' (Meirmans)`    = .clamp01(Fst_prime_meirmans_overall),
    `Fst' (Hedrick)`     = .clamp01(Fst_prime_hedrick_overall),
    GST                  = .clamp01(GST_overall),
    `GST''`              = .clamp01(GST2_overall),
    stringsAsFactors = FALSE
  )

  names(overall_row) <- names(out)
  rbind(out, overall_row)
}




.duckdb_write_table <- function(con, tbl, df, overwrite = TRUE) {
  stopifnot(DBI::dbIsValid(con))
  if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
  DBI::dbWriteTable(con, name = tbl, value = df, overwrite = overwrite)
  invisible(TRUE)
}

sql_ident <- function(con, x) {
  as.character(DBI::dbQuoteIdentifier(con, x))
}

duck_tbl_exists <- function(con, tbl) {
  isTRUE(tryCatch(DBI::dbExistsTable(con, tbl), error = function(e) FALSE))
}

duck_list_fields <- function(con, tbl) {
  tryCatch(DBI::dbListFields(con, tbl), error = function(e) character(0))
}

duck_count_rows <- function(con, tbl) {
  q <- sprintf("SELECT COUNT(*) AS n FROM %s", sql_ident(con, tbl))
  DBI::dbGetQuery(con, q)$n[[1]]
}

duck_pop_sizes <- function(con, tbl_meta) {
  if (!duck_tbl_exists(con, tbl_meta)) return(data.frame())
  q <- sprintf("
    SELECT Population, COUNT(*) AS Sample_Size
    FROM %s
    GROUP BY Population
    ORDER BY Sample_Size DESC, Population
  ", sql_ident(con, tbl_meta))
  DBI::dbGetQuery(con, q)
}

duck_n_pop <- function(con, tbl_meta) {
  if (!duck_tbl_exists(con, tbl_meta)) return(0L)
  q <- sprintf("SELECT COUNT(DISTINCT Population) AS n FROM %s", sql_ident(con, tbl_meta))
  as.integer(DBI::dbGetQuery(con, q)$n[[1]])
}

# ---- Gene diversity (Hs) per Population x Locus (DB-first)
duck_hs_by_pop_locus_long <- function(con,
                                      tbl_hf   = "hf",
                                      tbl_meta = "meta",
                                      base,
                                      missing_code = 0L) {
  stopifnot(DBI::dbIsValid(con))
  stopifnot(is.numeric(base), length(base) == 1L, base > 1)
  
  q <- sprintf("
    WITH g0 AS (
      SELECT
        m.Population AS Population,
        h.locus_id   AS Locus,
        h.gt         AS gt,
        CAST(floor(h.gt / %d) AS INTEGER) AS a1,
        CAST(h.gt %% %d AS INTEGER)       AS a2
      FROM %s h
      JOIN %s m
        ON m.individual = h.indiv_id
      WHERE h.gt IS NOT NULL
        AND h.gt <> %d
        AND h.gt > 0
        AND m.Population IS NOT NULL
    ),
    g AS (
      SELECT *
      FROM g0
      WHERE a1 > 0
        AND a2 > 0
    ),
    per_pop_locus AS (
      SELECT
        Population,
        Locus,
        COUNT(*)::DOUBLE AS n
      FROM g
      GROUP BY Population, Locus
    ),
    allele_counts AS (
      SELECT Population, Locus, a1 AS allele, COUNT(*)::DOUBLE AS c
      FROM g
      GROUP BY Population, Locus, a1
      UNION ALL
      SELECT Population, Locus, a2 AS allele, COUNT(*)::DOUBLE AS c
      FROM g
      GROUP BY Population, Locus, a2
    ),
    allele_summed AS (
      SELECT Population, Locus, allele, SUM(c)::DOUBLE AS c
      FROM allele_counts
      GROUP BY Population, Locus, allele
    )
    SELECT
      s.Locus,
      s.Population,
      (2.0*p.n / (2.0*p.n - 1.0)) * (1.0 - SUM(POWER(s.c / (2.0*p.n), 2))) AS Hs   -- unbiased
    FROM allele_summed s
    JOIN per_pop_locus p
      USING (Population, Locus)
    WHERE p.n > 1
    GROUP BY s.Locus, s.Population, p.n
    ORDER BY s.Locus, s.Population
  ",
               as.integer(base), as.integer(base),
               sql_ident(con, tbl_hf),
               sql_ident(con, tbl_meta),
               as.integer(missing_code)
  )
  
  DBI::dbGetQuery(con, q)
}




format_numeric_cols <- function(df, digits = 5, exclude = "Locus") {
  for (nm in names(df)) {
    if (nm %in% exclude) next
    x <- df[[nm]]
    if (is.factor(x)) x <- as.character(x)
    if (!is.numeric(x)) {
      suppressWarnings(x_num <- as.numeric(x))
    } else {
      x_num <- x
    }
    if (is.numeric(x_num) && !all(is.na(x_num))) {
      df[[nm]] <- sprintf(paste0("%.", digits, "f"), x_num)
    } else {
      df[[nm]] <- x
    }
  }
  df
}


#-------------------------------------------------#
## basic stats ####
#-------------------------------------------------#
# ---- DB-native population stats (NO genotype strings) -----------------------

duck_pop_stats_overall <- function(con, tbl_hf = "hf", tbl_meta = "meta",
                                   base, missing_code = 0L) {
  stopifnot(is.numeric(base), length(base) == 1L, base > 1)
  
  q <- sprintf("
    WITH g0 AS (
      SELECT
        m.Population AS Population,
        h.locus_id   AS Locus,
        CAST(floor(h.gt / %d) AS INTEGER) AS a1,
        CAST(h.gt %% %d AS INTEGER)       AS a2
      FROM %s h
      JOIN %s m
        ON m.individual = h.indiv_id
      WHERE h.gt IS NOT NULL
        AND h.gt <> %d
        AND h.gt > 0
    ),
    g AS (
      SELECT * FROM g0 WHERE a1 > 0 AND a2 > 0
    ),
    per_pop_locus AS (
      SELECT
        Population,
        Locus,
        COUNT(*)::DOUBLE AS n,
        SUM(CASE WHEN a1 <> a2 THEN 1 ELSE 0 END)::DOUBLE AS n_het
      FROM g
      GROUP BY Population, Locus
    ),
    allele_counts AS (
      SELECT Population, Locus, a1 AS allele, COUNT(*)::DOUBLE AS c
      FROM g
      GROUP BY Population, Locus, a1
      UNION ALL
      SELECT Population, Locus, a2 AS allele, COUNT(*)::DOUBLE AS c
      FROM g
      GROUP BY Population, Locus, a2
    ),
    allele_summed AS (
      SELECT
        Population, Locus, allele,
        SUM(c)::DOUBLE AS c
      FROM allele_counts
      GROUP BY Population, Locus, allele
    ),
    hs_by_pop_locus AS (
      SELECT
        s.Population,
        s.Locus,
        CASE
          WHEN p.n > 1 THEN ((2.0 * p.n) / (2.0 * p.n - 1.0)) *
               (1.0 - SUM(POWER(c / (2.0 * p.n), 2)))
          ELSE NULL
        END AS Hs
      FROM allele_summed s
      JOIN per_pop_locus p
        USING (Population, Locus)
      WHERE p.n > 0
      GROUP BY s.Population, s.Locus, p.n
    ),
    locus_stats AS (
      SELECT
        p.Population,
        p.Locus,
        p.n       AS n,
        p.n_het   AS n_het,
        (p.n_het / NULLIF(p.n, 0)) AS Ho,
        h.Hs AS Hs
      FROM per_pop_locus p
      LEFT JOIN hs_by_pop_locus h
        USING (Population, Locus)
    ),
    overall AS (
      SELECT
        Population,
        AVG(Ho) AS Ho,
        AVG(Hs) AS Hs,
        CASE
          WHEN SUM(
            CASE WHEN Hs IS NOT NULL AND Hs > 0 THEN n * Hs ELSE 0.0 END
          ) = 0 THEN NULL
          ELSE 1.0 - SUM(
            CASE WHEN Hs IS NOT NULL AND Hs > 0 THEN n_het ELSE 0.0 END
          ) / SUM(
            CASE WHEN Hs IS NOT NULL AND Hs > 0 THEN n * Hs ELSE 0.0 END
          )
        END AS Fis_WC
      FROM locus_stats
      GROUP BY Population
    )
    SELECT
      Population,
      Ho,
      Hs,
      Fis_WC AS \"Fis (WC)\"
    FROM overall
    ORDER BY Population
  ", as.integer(base), as.integer(base),
               DBI::dbQuoteIdentifier(con, tbl_hf),
               DBI::dbQuoteIdentifier(con, tbl_meta),
               as.integer(missing_code))
  
  DBI::dbGetQuery(con, q)
}

duck_pop_stats_by_pop_one <- function(con, pop_name,
                                      tbl_hf = "hf", tbl_meta = "meta",
                                      base, missing_code = 0L) {
  stopifnot(is.character(pop_name), length(pop_name) == 1L)
  
  q <- sprintf("
    WITH g0 AS (
      SELECT
        h.locus_id AS Locus,
        CAST(floor(h.gt / %d) AS INTEGER) AS a1,
        CAST(h.gt %% %d AS INTEGER)       AS a2
      FROM %s h
      JOIN %s m
        ON m.individual = h.indiv_id
      WHERE m.Population = ?
        AND h.gt IS NOT NULL
        AND h.gt <> %d
        AND h.gt > 0
    ),
    g AS (
      SELECT * FROM g0 WHERE a1 > 0 AND a2 > 0
    ),
    per_locus AS (
      SELECT
        Locus,
        COUNT(*)::DOUBLE AS n,
        SUM(CASE WHEN a1 <> a2 THEN 1 ELSE 0 END)::DOUBLE AS n_het
      FROM g
      GROUP BY Locus
    ),
    allele_counts AS (
      SELECT Locus, a1 AS allele, COUNT(*)::DOUBLE AS c
      FROM g
      GROUP BY Locus, a1
      UNION ALL
      SELECT Locus, a2 AS allele, COUNT(*)::DOUBLE AS c
      FROM g
      GROUP BY Locus, a2
    ),
    allele_summed AS (
      SELECT Locus, allele, SUM(c)::DOUBLE AS c
      FROM allele_counts
      GROUP BY Locus, allele
    ),
    hs_by_locus AS (
      SELECT
        s.Locus,
        CASE
          WHEN p.n > 1 THEN ((2.0 * p.n) / (2.0 * p.n - 1.0)) *
               (1.0 - SUM(POWER(c / (2.0 * p.n), 2)))
          ELSE NULL
        END AS Hs
      FROM allele_summed s
      JOIN per_locus p USING (Locus)
      WHERE p.n > 0
      GROUP BY s.Locus, p.n
    )
    SELECT
      p.Locus,
      (p.n_het / NULLIF(p.n, 0)) AS Ho,
      h.Hs AS Hs,
      CASE
        WHEN h.Hs IS NULL OR h.Hs <= 0 THEN NULL
        ELSE 1.0 - (p.n_het / p.n) / h.Hs
      END AS \"Fis (WC)\"
    FROM per_locus p
    LEFT JOIN hs_by_locus h USING (Locus)
    ORDER BY p.Locus
  ", as.integer(base), as.integer(base),
               DBI::dbQuoteIdentifier(con, tbl_hf),
               DBI::dbQuoteIdentifier(con, tbl_meta),
               as.integer(missing_code))
  
  DBI::dbGetQuery(con, q, params = list(pop_name))
}


compute_pop_stats_from_mat <- function(mat, base, pop_code) {
  
  # mat structure:
  # col 1 = pop
  # col 2..L+1 = packed genotypes
  
  pop_vec <- mat[, 1]
  geno    <- mat[, -1, drop = FALSE]
  
  idx <- which(pop_vec == pop_code)
  if (length(idx) == 0L) return(NULL)
  
  G <- geno[idx, , drop = FALSE]
  
  L <- ncol(G)
  
  Ho  <- numeric(L)
  Hs  <- numeric(L)
  Fis <- numeric(L)
  
  for (j in seq_len(L)) {
    
    g <- G[, j]
    g <- g[g != 0L]  # remove missing
    
    if (length(g) == 0L) {
      Ho[j]  <- NA_real_
      Hs[j]  <- NA_real_
      Fis[j] <- NA_real_
      next
    }
    
    # decode
    a1 <- g %/% base
    a2 <- g %%  base
    
    # Ho = proportion heterozygotes
    Ho[j] <- mean(a1 != a2)
    
    # allele vector
    alleles <- c(a1, a2)
    
    # allele frequencies
    tab <- table(alleles)
    p   <- tab / sum(tab)
    
    Hs[j] <- 1 - sum(p^2)
    
    Fis[j] <- if (Hs[j] > 0) 1 - Ho[j] / Hs[j] else NA_real_
  }
  
  data.frame(
    Locus = colnames(mat)[-1L],
    Ho = Ho,
    Hs = Hs,
    `Fis (Nei)` = Fis,
    stringsAsFactors = FALSE
  )
}

run_basic_wc84 <- function(mat, base) {
  
  res <- wc84_components_fst(
    dat          = mat,
    pop_col      = 0L,
    missing_code = 0L,
    base         = base
  )
  
  # assume res returns list with per_locus and overall
  
  per_locus <- as.data.frame(res$per_locus)
  
  colnames(per_locus) <- c(
    "Locus",
    "Fit (W&C)",
    "Fst (W&C)",
    "Fis (W&C)"
  )
  
  overall <- data.frame(
    Locus = "Overall",
    "Fit (W&C)" = res$overall_fit,
    "Fst (W&C)" = res$overall_fst,
    "Fis (W&C)" = res$overall_fis
  )
  
  rbind(per_locus, overall)
}

summarize_boot_ci <- function(boot_mat,
                              obs,
                              obs_overall = NA_real_,
                              boot_overall = NULL,
                              confidence = 0.95) {
  stopifnot(is.matrix(boot_mat), is.numeric(confidence), confidence > 0, confidence < 1)
  
  # Align obs to boot_mat columns if names exist
  if (!is.null(colnames(boot_mat))) {
    loc <- colnames(boot_mat)
    
    if (!is.null(names(obs))) {
      obs <- as.numeric(obs)[match(loc, names(obs))]
      names(obs) <- loc
    } else {
      obs <- as.numeric(obs)
      if (length(obs) == length(loc)) names(obs) <- loc
    }
  } else {
    obs <- as.numeric(obs)
  }
  
  alpha <- (1 - confidence) / 2
  probs <- c(alpha, 1 - alpha)
  
  # per-locus summaries
  boot_mean <- apply(boot_mat, 2, function(x) mean(x, na.rm = TRUE))
  ci <- t(apply(boot_mat, 2, function(x) {
    if (all(is.na(x))) return(c(NA_real_, NA_real_))
    as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
  }))
  
  ci_lo <- ci[, 1]
  ci_hi <- ci[, 2]
  
  # overall summaries (if provided)
  overall_mean  <- NA_real_
  overall_ci_lo <- NA_real_
  overall_ci_hi <- NA_real_
  
  if (!is.null(boot_overall)) {
    bo <- as.numeric(boot_overall)
    if (!all(is.na(bo))) {
      overall_mean <- mean(bo, na.rm = TRUE)
      q <- as.numeric(stats::quantile(bo, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
      overall_ci_lo <- q[1]
      overall_ci_hi <- q[2]
    }
  }
  
  list(
    mean           = boot_mean,
    ci_lo          = ci_lo,
    ci_hi          = ci_hi,
    overall_mean   = overall_mean,
    overall_ci_lo  = overall_ci_lo,
    overall_ci_hi  = overall_ci_hi,
    obs            = obs,
    obs_overall    = as.numeric(obs_overall),
    confidence     = confidence
  )
}

# ================== LINKAGE DISEQUILIBRIUM (LD tab) ===========================

#' Linkage disequilibrium (LD) helpers: contingency tables, G-test, permutations
#'
#' This block implements a simple LD test per population by:
#'  1) building haplotype contingency tables for all locus pairs,
#'  2) computing a likelihood-ratio (G) statistic for each table,
#'  3) estimating p-values by permutation (randomisation) within population.
#'
#' The code assumes:
#'  - `data$Population` defines populations (one value per sample/row),
#'  - loci columns contain haplotype/genotype strings (e.g. "92/100"),
#'  - missing genotypes are coded as "0/0" when `include_missing = FALSE`.
#'
#' Main functions
#' - `create_contingency_tables(data, loci, include_missing)`:
#'     For each population and each pair of loci, returns a contingency table
#'     (rows = locus1 haplotypes, cols = locus2 haplotypes). Optionally drops
#'     rows where either locus is "0/0".
#'
#' - `calculate_g_stat(contingency_table)`:
#'     Computes expected counts under independence and the G statistic:
#'       G = 2 * sum(O * log(O/E))   (only for cells with O > 0)
#'     Returns both `expected` and `g_stat`.
#'
#' - `add_g_stats(contingency_tables)`:
#'     Wraps each (pop, pair) table into a list containing:
#'       observed contingency table, expected table, and G statistic.
#'     Pair names are sanitised to use "." separators.
#'
#' - `randomized_g_stats(data, loci, n_simulations, calculate_g_stat, include_missing)`:
#'     Permutation test per population and locus pair:
#'       independently shuffles locus1 and locus2 columns within population,
#'       recomputes G, and stores the simulated G distribution.
#'     Parallelised over populations (foreach/doParallel).
#'
#' P-values and reporting
#' - `calculate_pvalues(observed_g_stats, simulated_g_stats, epsilon)`:
#'     For each (pop, pair), estimates p = mean(G_sim >= G_obs - epsilon),
#'     dropping NA simulations.
#'
#' - `calculate_global_pvalues(observed_g_stats, simulated_g_stats)`:
#'     Computes a crude “global” p-value per locus pair by comparing the mean
#'     observed G across populations to the pooled simulated G values.
#'
#' - `create_summary_table(pvalues, global_pvalues)`:
#'     Builds a data.frame with one row per locus pair, one column per population
#'     (p-values), plus a `Global_P_Value` column.
#'
#' @name ld_helpers
NULL



# Function to create contingency tables for each population
create_contingency_tables <- function(data, loci, include_missing = TRUE) {
  populations <- unique(data$Population)
  
  contingency_tables <- lapply(populations, function(pop) {
    pop_data <- data[data$Population == pop, ]
    locus_pairs <- combn(loci, 2, simplify = FALSE)
    
    contingency_list <- setNames(
      lapply(locus_pairs, function(pair) {
        locus1 <- pair[1]
        locus2 <- pair[2]
        
        if (!include_missing) {
          pop_data <- pop_data[pop_data[[locus1]] != "0/0" & pop_data[[locus2]] != "0/0", ]
        }
        
        haplotype_data <- data.frame(
          Locus1_haplotype = pop_data[[locus1]],
          Locus2_haplotype = pop_data[[locus2]]
        )
        
        table(haplotype_data$Locus1_haplotype, haplotype_data$Locus2_haplotype)
      }),
      sapply(locus_pairs, function(pair) paste(pair[1], pair[2], sep = "-"))
    )
    return(contingency_list)
  })
  names(contingency_tables) <- populations
  return(contingency_tables)
}

# Function to calculate G-statistic
calculate_g_stat <- function(contingency_table) {
  nt <- sum(contingency_table)
  row_sum <- rowSums(contingency_table)
  col_sum <- colSums(contingency_table)
  expected <- outer(row_sum, col_sum) / nt
  non_zero <- contingency_table > 0
  observed_non_zero <- contingency_table[non_zero]
  expected_non_zero <- expected[non_zero]
  g_stat <- 2 * sum(observed_non_zero * log(observed_non_zero / expected_non_zero), na.rm = TRUE)
  return(list(expected = expected, g_stat = g_stat))
}

# Function to add G-statistics to contingency tables
add_g_stats <- function(contingency_tables) {
  lapply(contingency_tables, function(pop_tables) {
    setNames(
      lapply(names(pop_tables), function(pair_name) {
        contingency_table <- pop_tables[[pair_name]]
        g_stat <- calculate_g_stat(contingency_table)
        list(
          contingency_table = contingency_table,
          expected_contingency_table = g_stat$expected,
          g_stat = g_stat$g_stat
        )
      }),
      gsub("[-`]+", ".", names(pop_tables))
    )
  })
}

# Function to calculate p-values
calculate_pvalues <- function(observed_g_stats, simulated_g_stats, epsilon = 1e-10) {
  results <- lapply(names(observed_g_stats), function(pop) {
    observed <- observed_g_stats[[pop]]
    simulated <- simulated_g_stats[[pop]]
    setNames(
      lapply(names(observed), function(pair) {
        observed_g <- observed[[pair]]$g_stat
        simulated_g <- simulated[[pair]][!is.na(simulated[[pair]])]
        p_value <- if (length(simulated_g) > 0) mean(simulated_g >= (observed_g - epsilon)) else NaN
        list(observed_g_stat = observed_g, p_value = p_value)
      }),
      names(observed)
    )
  })
  names(results) <- names(observed_g_stats)
  return(results)
}

# Function to calculate global p-values
calculate_global_pvalues <- function(observed_g_stats, simulated_g_stats) {
  locus_pairs <- unique(unlist(lapply(observed_g_stats, names)))
  sapply(locus_pairs, function(pair) {
    g_obs <- unlist(lapply(observed_g_stats, function(pop) pop[[pair]]$g_stat))
    g_sim <- unlist(lapply(simulated_g_stats, function(pop) pop[[pair]]))
    mean(g_sim >= mean(g_obs, na.rm = TRUE))
  })
}

# Function to create summary table
create_summary_table <- function(pvalues, global_pvalues) {
  all_pairs <- unique(unlist(lapply(pvalues, names)))
  summary_table <- data.frame(Locus_Pair = all_pairs)
  for (pop in names(pvalues)) {
    summary_table[[pop]] <- sapply(all_pairs, function(pair) pvalues[[pop]][[pair]]$p_value)
  }
  summary_table$Global_P_Value <- sapply(all_pairs, function(pair) global_pvalues[pair])
  return(summary_table)
}

# Function to generate randomized G-statistics
randomized_g_stats <- function(data, loci, n_simulations, calculate_g_stat, include_missing = TRUE) {
  workers <- parallel::detectCores() - 1
  cl <- makeCluster(workers)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("calculate_g_stat"), envir = environment())
  
  populations <- unique(data$Population)
  locus_pairs <- combn(loci, 2, simplify = FALSE)
  
  results <- foreach(pop = populations, .combine = 'c', .packages = 'dplyr') %dopar% {
    pop_data <- data[data$Population == pop, ]
    pop_results <- setNames(vector("list", length(locus_pairs)), sapply(locus_pairs, function(pair) gsub("[-`]+", ".", paste(pair[1], pair[2], sep = "."))))
    
    for (pair in locus_pairs) {
      locus1 <- pair[1]
      locus2 <- pair[2]
      g_stats <- numeric(n_simulations)
      for (i in 1:n_simulations) {
        randomized_data <- pop_data
        randomized_data[[locus1]] <- sample(pop_data[[locus1]])
        randomized_data[[locus2]] <- sample(pop_data[[locus2]])
        
        if (!include_missing) {
          randomized_data <- randomized_data[randomized_data[[locus1]] != "0/0" & randomized_data[[locus2]] != "0/0", ]
        }
        
        if (nrow(randomized_data) > 0) {
          contingency_table <- table(randomized_data[[locus1]], randomized_data[[locus2]])
          g_stats[i] <- calculate_g_stat(contingency_table)$g_stat
        }
      }
      pop_results[[gsub("[-`]+", ".", paste(locus1, locus2, sep = "."))]] <- g_stats
    }
    list(setNames(list(pop_results), pop))
  }
  
  stopCluster(cl)
  return(do.call(c, results))
}