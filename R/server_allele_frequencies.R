server_allele_frequencies <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    `%||%` <- function(a, b) if (!is.null(a)) a else b

    safe_choice <- function(x, default = "all") {
      if (is.null(x) || length(x) == 0L || identical(x, "") || all(is.na(x))) default
      else as.character(x[[1]])
    }

    sql_id  <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))
    sql_str <- function(con, x) as.character(DBI::dbQuoteString(con, x))

    # ── Reactive plumbing ──────────────────────────────────────────────────
    db_tick    <- reactive({ rv$db_tick })
    con_r      <- reactive({ req(rv$con); rv$con })
    tbl_meta_r <- reactive({ rv$tbl_meta %||% "meta" })

    tbl_hf_r <- reactive({
      con <- con_r()
      if (exists("duck_tbl_exists",    mode = "function", inherits = TRUE) &&
          exists(".duckdb_get_params", mode = "function", inherits = TRUE) &&
          duck_tbl_exists(con, "params")) {
        p <- .duckdb_get_params(con); return(as.character(p$tbl_hf %||% "hf"))
      }
      "hf"
    })

    db_ready <- reactive({
      db_tick(); con <- con_r()
      shiny::req(isTRUE(rv$db_ready))
      shiny::validate(
        shiny::need(DBI::dbExistsTable(con, tbl_meta_r()), "DuckDB meta table missing."),
        shiny::need(DBI::dbExistsTable(con, tbl_hf_r()),   "DuckDB hf table missing.")
      )
      TRUE
    })

    base_r <- reactive({
      db_ready()
      b <- rv$base_af %||% rv$base %||% rv$base_r %||% rv$genotype_base
      b <- suppressWarnings(as.integer(b))
      if (length(b) == 1L && is.finite(b) && b > 1L) return(as.integer(b))
      con <- con_r()
      if (DBI::dbExistsTable(con, "params") &&
          exists(".duckdb_get_params", mode = "function", inherits = TRUE)) {
        p <- .duckdb_get_params(con)
        b <- suppressWarnings(as.integer(
          p$base %||% p$base_scalar_full %||% p$base_scalar_preview))
        if (length(b) == 1L && is.finite(b) && b > 1L) return(as.integer(b))
      }
      1000L
    })

    hf_schema_r <- reactive({
      db_ready(); con <- con_r()
      info <- DBI::dbGetQuery(con,
        sprintf("PRAGMA table_info(%s)", DBI::dbQuoteIdentifier(con, tbl_hf_r())))
      cols <- info$name
      if (all(c("individual","locus","g") %in% cols))
        return(list(ind_col="individual", locus_col="locus",    gt_col="g"))
      if (all(c("indiv_id","locus_id","gt") %in% cols))
        return(list(ind_col="indiv_id",   locus_col="locus_id", gt_col="gt"))
      shiny::validate(shiny::need(FALSE,
        "hf must contain (individual,locus,g) or (indiv_id,locus_id,gt)."))
    })

    meta_schema_r <- reactive({
      db_ready(); con <- con_r()
      info <- DBI::dbGetQuery(con,
        sprintf("PRAGMA table_info(%s)", DBI::dbQuoteIdentifier(con, tbl_meta_r())))
      cols    <- info$name
      ind_col <- if ("individual" %in% cols) "individual"
                 else if ("indiv_id" %in% cols) "indiv_id"
                 else shiny::validate(shiny::need(FALSE,
                   "No individual column in meta."))
      pop_col <- c("Population","population","pop","pop_code")[
        c("Population","population","pop","pop_code") %in% cols][1]
      shiny::validate(shiny::need(!is.na(pop_col),
        "No population column in meta."))
      list(ind_col=ind_col, pop_col=pop_col)
    })

    # ── Locus order (first physical appearance) ────────────────────────────
    locus_order_cte <- function(con, hf_tbl_q, hl_q)
      sprintf("locus_order AS (
  SELECT CAST(%s AS VARCHAR) AS _lo_marker, MIN(rowid) AS _lo_rank
  FROM %s GROUP BY CAST(%s AS VARCHAR))", hl_q, hf_tbl_q, hl_q)

    # ── Population / marker lists ──────────────────────────────────────────
    pops_r <- reactive({
      db_ready(); con <- con_r(); ms <- meta_schema_r()
      as.character(DBI::dbGetQuery(con, sprintf(
        "SELECT DISTINCT CAST(%s AS VARCHAR) AS p FROM %s
         WHERE %s IS NOT NULL ORDER BY p",
        sql_id(con,ms$pop_col), sql_id(con,tbl_meta_r()),
        sql_id(con,ms$pop_col)))$p)
    })

    markers_r <- reactive({
      db_ready(); con <- con_r(); hs <- hf_schema_r()
      hf_q <- sql_id(con,tbl_hf_r()); hl_q <- sql_id(con,hs$locus_col)
      as.character(DBI::dbGetQuery(con, sprintf("
        WITH %s
        SELECT DISTINCT CAST(%s AS VARCHAR) AS Marker, lo._lo_rank
        FROM %s h LEFT JOIN locus_order lo
          ON CAST(%s AS VARCHAR)=lo._lo_marker
        ORDER BY lo._lo_rank ASC",
        locus_order_cte(con,hf_q,hl_q), hl_q, hf_q, hl_q))$Marker)
    })

    observe({
      pops <- pops_r(); markers <- markers_r()
      # existing tab selects
      updateSelectInput(session,"selected_population",
        choices=c("All populations"="all", stats::setNames(pops,pops)), selected="all")
      updateSelectizeInput(session,"selected_marker",
        choices=c("All markers"="all", stats::setNames(markers,markers)),
        selected="all", server=TRUE)
      updateSelectInput(session,"selected_population_subsamples",
        choices=c("All populations"="all", stats::setNames(pops,pops)), selected="all")
      # new fstat tab selects
      updateSelectInput(session,"fstat_population",
        choices=c("All populations"="all", stats::setNames(pops,pops)), selected="all")
      updateSelectizeInput(session,"fstat_marker",
        choices=c("All markers"="all", stats::setNames(markers,markers)),
        selected="all", server=TRUE)
    })

    selected_population_r            <- reactive(safe_choice(input$selected_population,"all"))
    selected_marker_r                <- reactive(safe_choice(input$selected_marker,"all"))
    selected_population_subsamples_r <- reactive(
      safe_choice(input$selected_population_subsamples,"all"))
    fstat_population_r               <- reactive(safe_choice(input$fstat_population,"all"))
    fstat_marker_r                   <- reactive(safe_choice(input$fstat_marker,"all"))

    analysis_ready_r <- reactive({
      req(input$update_analysis > 0L); db_ready(); TRUE })

    fstat_ready_r <- reactive({
      req(input$update_fstat > 0L); db_ready(); TRUE })

    # ── Missing data by pop × locus ────────────────────────────────────────
    missing_by_pop_locus_r <- reactive({
      db_ready(); con <- con_r(); hs <- hf_schema_r(); ms <- meta_schema_r()
      hf_q <- sql_id(con,tbl_hf_r()); meta_q <- sql_id(con,tbl_meta_r())
      hi_q <- sql_id(con,hs$ind_col); hl_q <- sql_id(con,hs$locus_col)
      hg_q <- sql_id(con,hs$gt_col);  mi_q <- sql_id(con,ms$ind_col)
      pop_q <- sql_id(con,ms$pop_col)
      DBI::dbGetQuery(con, sprintf("
        WITH %s,
        base AS (
          SELECT CAST(m.%s AS VARCHAR) AS Population,
                 CAST(h.%s AS VARCHAR) AS Marker,
            COUNT(*) AS Sample_Size,
            SUM(CASE WHEN h.%s IS NULL OR h.%s<=0 THEN 1 ELSE 0 END) AS Missing_Data,
            SUM(CASE WHEN h.%s IS NOT NULL AND h.%s>0 THEN 1 ELSE 0 END) AS Genotyped_Data,
            SUM(CASE WHEN h.%s IS NULL OR h.%s<=0 THEN 1 ELSE 0 END)*1.0/COUNT(*)
              AS Missing_Proportion
          FROM %s h INNER JOIN %s m
            ON CAST(h.%s AS VARCHAR)=CAST(m.%s AS VARCHAR)
          WHERE m.%s IS NOT NULL GROUP BY m.%s, h.%s)
        SELECT b.* FROM base b
        LEFT JOIN locus_order lo ON b.Marker=lo._lo_marker
        ORDER BY b.Population, lo._lo_rank ASC",
        locus_order_cte(con,hf_q,hl_q),
        pop_q,hl_q, hg_q,hg_q, hg_q,hg_q, hg_q,hg_q,
        hf_q,meta_q, hi_q,mi_q, pop_q,pop_q,hl_q))
    })

    # ── Existing allele-frequency helpers ──────────────────────────────────
    build_af_by_pop_sql <- function(sel_pop="all", sel_mark="all") {
      con <- con_r(); hs <- hf_schema_r(); ms <- meta_schema_r()
      base <- as.integer(base_r())
      hf_q <- sql_id(con,tbl_hf_r()); meta_q <- sql_id(con,tbl_meta_r())
      hi_q <- sql_id(con,hs$ind_col); hl_q <- sql_id(con,hs$locus_col)
      hg_q <- sql_id(con,hs$gt_col);  mi_q <- sql_id(con,ms$ind_col)
      pop_q <- sql_id(con,ms$pop_col)
      ew <- character(0)
      if (!identical(sel_pop,"all"))
        ew <- c(ew, sprintf("CAST(m.%s AS VARCHAR)=%s",pop_q,sql_str(con,sel_pop)))
      if (!identical(sel_mark,"all"))
        ew <- c(ew, sprintf("CAST(h.%s AS VARCHAR)=%s",hl_q,sql_str(con,sel_mark)))
      ew_sql <- if(length(ew)) paste0("\n        AND ",paste(ew,collapse="\n        AND ")) else ""
      sprintf("
    WITH %s,
    valid_gt AS (
      SELECT CAST(m.%s AS VARCHAR) AS Population,
             CAST(h.%s AS VARCHAR) AS Marker, h.%s AS gt
      FROM %s h INNER JOIN %s m ON CAST(h.%s AS VARCHAR)=CAST(m.%s AS VARCHAR)
      WHERE m.%s IS NOT NULL AND h.%s IS NOT NULL AND h.%s>0%s),
    alleles AS (
      SELECT Population,Marker,CAST(FLOOR(gt/%d) AS BIGINT) AS Allele FROM valid_gt
      UNION ALL
      SELECT Population,Marker,CAST(gt%%%d AS BIGINT) AS Allele FROM valid_gt),
    counts AS (
      SELECT Population,Marker,Allele,COUNT(*) AS Count
      FROM alleles GROUP BY Population,Marker,Allele)
    SELECT c.Population,c.Marker,CAST(c.Allele AS VARCHAR) AS Allele,c.Count,
      SUM(c.Count) OVER(PARTITION BY c.Population,c.Marker) AS N_total,
      c.Count*1.0/SUM(c.Count) OVER(PARTITION BY c.Population,c.Marker) AS Frequency
    FROM counts c
    LEFT JOIN locus_order lo ON c.Marker=lo._lo_marker
    ORDER BY c.Population,lo._lo_rank ASC,c.Allele",
        locus_order_cte(con,hf_q,hl_q),
        pop_q,hl_q,hg_q, hf_q,meta_q,hi_q,mi_q,
        pop_q,hg_q,hg_q, ew_sql, base,base)
    }

    build_af_global_sql <- function(sel_mark="all") {
      con <- con_r(); hs <- hf_schema_r(); base <- as.integer(base_r())
      hf_q <- sql_id(con,tbl_hf_r()); hl_q <- sql_id(con,hs$locus_col)
      hg_q <- sql_id(con,hs$gt_col)
      ew_sql <- if(!identical(sel_mark,"all"))
        sprintf("\n        AND CAST(h.%s AS VARCHAR)=%s",hl_q,sql_str(con,sel_mark)) else ""
      sprintf("
    WITH %s,
    valid_gt AS (
      SELECT CAST(h.%s AS VARCHAR) AS Marker, h.%s AS gt
      FROM %s h WHERE h.%s IS NOT NULL AND h.%s>0%s),
    alleles AS (
      SELECT Marker,CAST(FLOOR(gt/%d) AS BIGINT) AS Allele FROM valid_gt
      UNION ALL
      SELECT Marker,CAST(gt%%%d AS BIGINT) AS Allele FROM valid_gt),
    counts AS (SELECT Marker,Allele,COUNT(*) AS Count FROM alleles GROUP BY Marker,Allele)
    SELECT c.Marker,CAST(c.Allele AS VARCHAR) AS Allele,c.Count,
      SUM(c.Count) OVER(PARTITION BY c.Marker) AS N_total,
      c.Count*1.0/SUM(c.Count) OVER(PARTITION BY c.Marker) AS Frequency
    FROM counts c
    LEFT JOIN locus_order lo ON c.Marker=lo._lo_marker
    ORDER BY lo._lo_rank ASC,c.Allele",
        locus_order_cte(con,hf_q,hl_q),
        hl_q,hg_q,hf_q,hg_q,hg_q,ew_sql,base,base)
    }

    run_af_by_pop <- function(sel_pop="all",sel_mark="all") {
      db_ready(); con <- con_r()
      out <- DBI::dbGetQuery(con,build_af_by_pop_sql(sel_pop,sel_mark))
      if(nrow(out)>0L) out$Frequency <- round(as.numeric(out$Frequency),4)
      out
    }
    run_af_global <- function(sel_mark="all") {
      db_ready(); con <- con_r()
      out <- DBI::dbGetQuery(con,build_af_global_sql(sel_mark))
      if(nrow(out)>0L) out$Frequency <- round(as.numeric(out$Frequency),4)
      out
    }

    af_table_by_pop_r <- reactive({
      analysis_ready_r()
      run_af_by_pop(selected_population_r(),selected_marker_r())
    })
    af_table_global_r <- reactive({
      analysis_ready_r()
      run_af_global(selected_marker_r())
    })

    # ── Genetic diversity ──────────────────────────────────────────────────
    build_genetic_diversity_sql <- function() {
      con <- con_r(); hs <- hf_schema_r(); ms <- meta_schema_r()
      base <- as.integer(base_r())
      hf_q <- sql_id(con,tbl_hf_r()); meta_q <- sql_id(con,tbl_meta_r())
      hi_q <- sql_id(con,hs$ind_col); hl_q <- sql_id(con,hs$locus_col)
      hg_q <- sql_id(con,hs$gt_col);  mi_q <- sql_id(con,ms$ind_col)
      pop_q <- sql_id(con,ms$pop_col)
      sprintf("
    WITH %s,
    valid_gt AS (
      SELECT CAST(m.%s AS VARCHAR) AS Population,
             CAST(h.%s AS VARCHAR) AS Marker, h.%s AS gt
      FROM %s h INNER JOIN %s m ON CAST(h.%s AS VARCHAR)=CAST(m.%s AS VARCHAR)
      WHERE m.%s IS NOT NULL AND h.%s IS NOT NULL AND h.%s>0),
    het AS (SELECT Population,Marker,COUNT(*) AS N_genotypes,
      SUM(CASE WHEN CAST(FLOOR(gt/%d) AS BIGINT)<>CAST(gt%%%d AS BIGINT)
               THEN 1 ELSE 0 END) AS N_heterozygotes
      FROM valid_gt GROUP BY Population,Marker),
    alleles AS (
      SELECT Population,Marker,CAST(FLOOR(gt/%d) AS BIGINT) AS Allele FROM valid_gt
      UNION ALL
      SELECT Population,Marker,CAST(gt%%%d AS BIGINT) AS Allele FROM valid_gt),
    counts AS (SELECT Population,Marker,Allele,COUNT(*) AS Count
      FROM alleles GROUP BY Population,Marker,Allele),
    freqs AS (SELECT Population,Marker,Allele,
      Count*1.0/SUM(Count) OVER(PARTITION BY Population,Marker) AS Frequency FROM counts),
    diversity AS (SELECT Population,Marker,
      COUNT(*) AS Na,
      1.0/SUM(Frequency*Frequency) AS Ne,
      1.0-SUM(Frequency*Frequency) AS He
      FROM freqs GROUP BY Population,Marker)
    SELECT d.Population,d.Marker,d.Na,
      ROUND(d.Ne,3) AS Ne,ROUND(d.He,3) AS He,
      ROUND(CASE WHEN h.N_genotypes>0
        THEN h.N_heterozygotes*1.0/h.N_genotypes ELSE NULL END,3) AS Ho,
      ROUND(CASE WHEN d.He>0 AND h.N_genotypes>0
        THEN (d.He-h.N_heterozygotes*1.0/h.N_genotypes)/d.He ELSE NULL END,3) AS Fis,
      COALESCE(h.N_genotypes,0) AS N_genotypes
    FROM diversity d
    LEFT JOIN het h USING(Population,Marker)
    LEFT JOIN locus_order lo ON d.Marker=lo._lo_marker
    ORDER BY d.Population,lo._lo_rank ASC",
        locus_order_cte(con,hf_q,hl_q),
        pop_q,hl_q,hg_q, hf_q,meta_q,hi_q,mi_q,
        pop_q,hg_q,hg_q, base,base, base,base)
    }

    genetic_diversity_r <- reactive({
      analysis_ready_r(); con <- con_r()
      DBI::dbGetQuery(con,build_genetic_diversity_sql())
    })

    # ── Marker details ─────────────────────────────────────────────────────
    build_marker_details_sql <- function() {
      con <- con_r(); hs <- hf_schema_r(); ms <- meta_schema_r()
      base <- as.integer(base_r())
      hf_q <- sql_id(con,tbl_hf_r()); meta_q <- sql_id(con,tbl_meta_r())
      hi_q <- sql_id(con,hs$ind_col); hl_q <- sql_id(con,hs$locus_col)
      hg_q <- sql_id(con,hs$gt_col);  mi_q <- sql_id(con,ms$ind_col)
      pop_q <- sql_id(con,ms$pop_col)
      sprintf("
    WITH %s,
    valid_gt AS (
      SELECT CAST(h.%s AS VARCHAR) AS Marker, h.%s AS gt
      FROM %s h WHERE h.%s IS NOT NULL AND h.%s>0),
    alleles AS (
      SELECT Marker,CAST(FLOOR(gt/%d) AS BIGINT) AS Allele FROM valid_gt
      UNION ALL
      SELECT Marker,CAST(gt%%%d AS BIGINT) AS Allele FROM valid_gt),
    allele_counts AS (SELECT Marker,Allele,COUNT(*) AS Count FROM alleles
      GROUP BY Marker,Allele),
    allele_summary AS (SELECT Marker,COUNT(*) AS N_alleles,SUM(Count) AS Total_observations,
      string_agg(CAST(Allele AS VARCHAR),', ' ORDER BY Allele) AS Alleles
      FROM allele_counts GROUP BY Marker),
    missing_by_pop AS (
      SELECT CAST(m.%s AS VARCHAR) AS Population,
             CAST(h.%s AS VARCHAR) AS Marker,
        SUM(CASE WHEN h.%s IS NULL OR h.%s<=0 THEN 1 ELSE 0 END)*1.0/COUNT(*)
          AS Missing_Proportion
      FROM %s h INNER JOIN %s m ON CAST(h.%s AS VARCHAR)=CAST(m.%s AS VARCHAR)
      WHERE m.%s IS NOT NULL GROUP BY m.%s,h.%s),
    missing_summary AS (SELECT Marker,AVG(Missing_Proportion) AS Avg_Missing_Proportion
      FROM missing_by_pop GROUP BY Marker)
    SELECT a.Marker,a.N_alleles,a.Total_observations,
      ROUND(COALESCE(m.Avg_Missing_Proportion,0),3) AS Avg_Missing_Proportion,a.Alleles
    FROM allele_summary a
    LEFT JOIN missing_summary m USING(Marker)
    LEFT JOIN locus_order lo ON a.Marker=lo._lo_marker
    ORDER BY lo._lo_rank ASC",
        locus_order_cte(con,hf_q,hl_q),
        hl_q,hg_q,hf_q,hg_q,hg_q, base,base,
        pop_q,hl_q,hg_q,hg_q,
        hf_q,meta_q,hi_q,mi_q,
        pop_q,pop_q,hl_q)
    }

    marker_details_r <- reactive({
      analysis_ready_r(); con <- con_r()
      DBI::dbGetQuery(con,build_marker_details_sql())
    })

    # ══════════════════════════════════════════════════════════════════════
    # NEW — Fstat-style table reactive
    # ══════════════════════════════════════════════════════════════════════
    fstat_long_r <- reactive({
      fstat_ready_r()
      con <- con_r(); hs <- hf_schema_r(); ms <- meta_schema_r()
      base <- as.integer(base_r())
      hf_q <- sql_id(con,tbl_hf_r()); meta_q <- sql_id(con,tbl_meta_r())
      hi_q <- sql_id(con,hs$ind_col); hl_q <- sql_id(con,hs$locus_col)
      hg_q <- sql_id(con,hs$gt_col);  mi_q <- sql_id(con,ms$ind_col)
      pop_q <- sql_id(con,ms$pop_col)

      pop_f  <- if (!identical(fstat_population_r(),"all"))
        sprintf("AND CAST(m.%s AS VARCHAR)=%s",pop_q,sql_str(con,fstat_population_r())) else ""
      mark_f <- if (!identical(fstat_marker_r(),"all"))
        sprintf("AND CAST(h.%s AS VARCHAR)=%s",hl_q,sql_str(con,fstat_marker_r())) else ""

      DBI::dbGetQuery(con, sprintf("
WITH %s,
valid_gt AS (
  SELECT CAST(m.%s AS VARCHAR) AS Population,
         CAST(h.%s AS VARCHAR) AS Marker, h.%s AS gt
  FROM %s h INNER JOIN %s m ON CAST(h.%s AS VARCHAR)=CAST(m.%s AS VARCHAR)
  WHERE m.%s IS NOT NULL AND h.%s IS NOT NULL AND h.%s>0 %s %s),
allele_universe AS (
  SELECT DISTINCT Marker, CAST(FLOOR(gt/%d) AS BIGINT) AS Allele FROM valid_gt
  UNION
  SELECT DISTINCT Marker, CAST(gt%%%d AS BIGINT) AS Allele FROM valid_gt),
populations AS (SELECT DISTINCT Population FROM valid_gt),
grid AS (
  SELECT p.Population, a.Marker, a.Allele
  FROM populations p CROSS JOIN allele_universe a),
alleles_split AS (
  SELECT Population,Marker,CAST(FLOOR(gt/%d) AS BIGINT) AS Allele FROM valid_gt
  UNION ALL
  SELECT Population,Marker,CAST(gt%%%d AS BIGINT) AS Allele FROM valid_gt),
counts AS (
  SELECT Population,Marker,Allele,COUNT(*) AS n
  FROM alleles_split GROUP BY Population,Marker,Allele),
totals AS (
  SELECT Population,Marker,SUM(n) AS N_alleles
  FROM counts GROUP BY Population,Marker),
freqs AS (
  SELECT g.Population,g.Marker,g.Allele,
    COALESCE(c.n*1.0/NULLIF(t.N_alleles,0),0.0) AS Frequency,
    COALESCE(c.n,0) AS Count, COALESCE(t.N_alleles,0) AS N_alleles
  FROM grid g
  LEFT JOIN counts c ON g.Population=c.Population AND g.Marker=c.Marker AND g.Allele=c.Allele
  LEFT JOIN totals t ON g.Population=t.Population AND g.Marker=t.Marker),
div AS (
  SELECT Population,Marker,
    COUNT(CASE WHEN Frequency>0 THEN 1 END) AS Na,
    1.0/SUM(Frequency*Frequency) AS Ne,
    1.0-SUM(Frequency*Frequency) AS He
  FROM freqs WHERE Frequency>0 GROUP BY Population,Marker),
het AS (
  SELECT Population,Marker,COUNT(*) AS N_gt,
    SUM(CASE WHEN CAST(FLOOR(gt/%d) AS BIGINT)<>CAST(gt%%%d AS BIGINT)
             THEN 1 ELSE 0 END) AS N_het
  FROM valid_gt GROUP BY Population,Marker),
sample_info AS (
  SELECT CAST(m.%s AS VARCHAR) AS Population,
         CAST(h.%s AS VARCHAR) AS Marker,
    SUM(CASE WHEN h.%s IS NOT NULL AND h.%s>0 THEN 1 ELSE 0 END) AS N_genotyped,
    SUM(CASE WHEN h.%s IS NULL OR h.%s<=0    THEN 1 ELSE 0 END) AS N_missing
  FROM %s h INNER JOIN %s m ON CAST(h.%s AS VARCHAR)=CAST(m.%s AS VARCHAR)
  WHERE m.%s IS NOT NULL %s %s
  GROUP BY m.%s,h.%s)
SELECT f.Population,f.Marker,CAST(f.Allele AS VARCHAR) AS Allele,
  ROUND(f.Frequency,4) AS Frequency, f.Count, f.N_alleles,
  s.N_genotyped, s.N_missing,
  d.Na, ROUND(d.Ne,3) AS Ne, ROUND(d.He,3) AS He,
  ROUND(CASE WHEN h.N_gt>0 THEN h.N_het*1.0/h.N_gt ELSE NULL END,3) AS Ho,
  ROUND(CASE WHEN d.He>0 AND h.N_gt>0
    THEN (d.He-h.N_het*1.0/h.N_gt)/d.He ELSE NULL END,3) AS Fis
FROM freqs f
LEFT JOIN sample_info s ON f.Population=s.Population AND f.Marker=s.Marker
LEFT JOIN div         d ON f.Population=d.Population AND f.Marker=d.Marker
LEFT JOIN het         h ON f.Population=h.Population AND f.Marker=h.Marker
LEFT JOIN locus_order lo ON f.Marker=lo._lo_marker
ORDER BY lo._lo_rank ASC, f.Population, f.Allele",
        locus_order_cte(con,hf_q,hl_q),
        pop_q,hl_q,hg_q, hf_q,meta_q,hi_q,mi_q,
        pop_q,hg_q,hg_q, pop_f,mark_f,
        base,base, base,base, base,base,
        pop_q,hl_q, hg_q,hg_q, hg_q,hg_q,
        hf_q,meta_q,hi_q,mi_q, pop_q,pop_f,mark_f, pop_q,hl_q))
    })

    # ── Pivot to wide display ──────────────────────────────────────────────
    fstat_wide_r <- reactive({
      long <- req(fstat_long_r())
      if (nrow(long)==0L) return(NULL)

      pops    <- unique(long$Population)
      markers <- unique(long$Marker)

      # Global allele frequencies
      gc  <- aggregate(Count~Marker+Allele, data=long, FUN=sum)
      gt  <- aggregate(Count~Marker,        data=long, FUN=sum)
      names(gt)[2] <- "Total"
      gc  <- merge(gc,gt,by="Marker")
      gc$Global_Freq <- gc$Count/gc$Total

      rows <- list()

      for (loc in markers) {
        lng  <- long[long$Marker==loc,,drop=FALSE]
        alleles <- sort(unique(as.numeric(lng$Allele)))

        mk_row <- function(type, label, vals, global_val=NA_character_) {
          r <- c(Locus="", Row_label=label, Row_type=type,
                 stats::setNames(sapply(pops, function(p) {
                   v <- vals[[p]]
                   if (is.null(v)||is.na(v)) "" else as.character(v)
                 }), pops),
                 Global=if(is.na(global_val)) "" else as.character(global_val))
          r
        }

        # N genotyped
        ng <- stats::setNames(lapply(pops, function(p) {
          x <- lng$N_genotyped[lng$Population==p]; if(length(x)) x[1] else NA}), pops)
        gng <- sum(sapply(ng, function(v) if(is.null(v)||is.na(v)) 0 else v))
        r <- mk_row("stat","N genotyped",ng,gng)
        r["Locus"] <- loc   # locus name on first row only
        rows[[length(rows)+1]] <- r

        # N missing
        nm <- stats::setNames(lapply(pops, function(p) {
          x <- lng$N_missing[lng$Population==p]; if(length(x)) x[1] else NA}), pops)
        gnm <- sum(sapply(nm, function(v) if(is.null(v)||is.na(v)) 0 else v))
        rows[[length(rows)+1]] <- mk_row("stat","N missing",nm,gnm)

        # Allele rows (with zero-fill)
        for (al in alleles) {
          al_c <- as.character(al)
          fp <- stats::setNames(lapply(pops, function(p) {
            f <- lng$Frequency[lng$Population==p & lng$Allele==al_c]
            if(length(f)) sprintf("%.4f",f[1]) else "0.0000"
          }), pops)
          gcr   <- gc[gc$Marker==loc & gc$Allele==al_c,,drop=FALSE]
          gfreq <- if(nrow(gcr)) sprintf("%.4f",gcr$Global_Freq[1]) else "0.0000"
          rows[[length(rows)+1]] <- mk_row("allele",al_c,fp,gfreq)
        }

        # Diversity stats
        for (stat in c("Na","Ne","He","Ho","Fis")) {
          sv <- stats::setNames(lapply(pops, function(p) {
            x <- lng[[stat]][lng$Population==p]
            if(!length(x)||is.na(x[1])) return(NA)
            if(stat=="Na") as.character(x[1]) else sprintf("%.3f",x[1])
          }), pops)
          rows[[length(rows)+1]] <- mk_row("div_stat",stat,sv,NA_character_)
        }
      }

      df <- do.call(rbind, lapply(rows, function(r)
        as.data.frame(t(r), stringsAsFactors=FALSE)))
      attr(df,"pops") <- pops
      df
    })

    # ── Fstat DT render ────────────────────────────────────────────────────
    output$fstat_table <- DT::renderDT({
      wide <- req(fstat_wide_r())
      pops <- attr(wide,"pops") %||%
        setdiff(names(wide),c("Locus","Row_label","Row_type","Global"))

      display   <- wide[,c("Locus","Row_label",pops,"Global"),drop=FALSE]
      col_labels <- c("Locus","Allele / stat",pops,"Global")

      idx_stat <- which(wide$Row_type=="stat")     - 1L
      idx_div  <- which(wide$Row_type=="div_stat") - 1L

      DT::datatable(
        display,
        rownames=FALSE, colnames=col_labels,
        selection="none",
        class="compact hover stripe",
        options=list(
          pageLength=200, scrollX=TRUE, ordering=FALSE, dom="lrtip",
          rowCallback=DT::JS(sprintf("
function(row,data,index){
  var si=[%s], di=[%s];
  if(si.indexOf(index)>-1){
    $('td',row).css({'font-size':'11px','color':'#6b7280','background':'#f9fafb'});
    $('td:eq(1)',row).css('font-style','italic');
    if(data[1]==='N missing'){
      $('td',row).slice(2).each(function(){
        var v=parseInt($(this).text());
        if(!isNaN(v)&&v>0)
          $(this).css({'color':'#854F0B','background':'#FAEEDA'});
      });
    }
  }
  if(di.indexOf(index)>-1){
    $('td',row).css({'font-size':'11px','color':'#374151','background':'#f3f4f6'});
    $('td:eq(1)',row).css({'font-weight':'500','font-style':'normal'});
  }
  if(data[0]!==''){
    $('td:eq(0)',row).css({'font-weight':'600','font-size':'12px'});
    $('td',row).css('border-top','2px solid #d1d5db');
  }
  $('td',row).slice(2).each(function(){
    if($(this).text()==='0.0000') $(this).css('color','#d1d5db');
  });
}",
            paste(idx_stat,collapse=","),
            paste(idx_div, collapse=","))))
      )
    }, server=TRUE)

    # ── Value boxes ────────────────────────────────────────────────────────
    n_individuals_r <- reactive({
      db_ready(); con <- con_r(); ms <- meta_schema_r()
      DBI::dbGetQuery(con, sprintf(
        "SELECT COUNT(DISTINCT CAST(%s AS VARCHAR)) AS n FROM %s WHERE %s IS NOT NULL",
        sql_id(con,ms$ind_col), sql_id(con,tbl_meta_r()),
        sql_id(con,ms$ind_col)))$n[[1]]
    })
    n_populations_r <- reactive({ length(pops_r()) })
    n_markers_r     <- reactive({ length(markers_r()) })

    summary_trigger_r <- eventReactive(input$generate_summary, {
      list(md=missing_by_pop_locus_r(),
           n_ind=n_individuals_r(), n_pop=n_populations_r(), n_mark=n_markers_r())
    }, ignoreInit=TRUE)

    output$vb_individuals <- renderUI({
      s <- summary_trigger_r()
      tags$div(class="af-vbox-val", if(!is.null(s)) s$n_ind else "\u2014")
    })
    output$vb_populations <- renderUI({
      s <- summary_trigger_r()
      tags$div(class="af-vbox-val", if(!is.null(s)) s$n_pop else "\u2014")
    })
    output$vb_markers <- renderUI({
      s <- summary_trigger_r()
      tags$div(class="af-vbox-val", if(!is.null(s)) s$n_mark else "\u2014")
    })
    output$vb_missing <- renderUI({
      s <- summary_trigger_r()
      if(is.null(s)||nrow(s$md)==0)
        return(tags$div(class="af-vbox-val","\u2014"))
      prop  <- round(sum(s$md$Missing_Data)/sum(s$md$Sample_Size),3)
      color <- if(prop>.20)"#A32D2D" else if(prop>.10)"#854F0B" else "#3B6D11"
      tags$div(class="af-vbox-val",style=paste0("color:",color,";"),prop)
    })

    # ── Summary tables ─────────────────────────────────────────────────────
    output$populations_Summary_table <- renderTable({
      # Return a placeholder until the user clicks "Generate summary"
      s <- summary_trigger_r()
      if (is.null(s)) {
        return(data.frame(
          Message = "Click 'Generate data summary' to populate",
          stringsAsFactors = FALSE
        ))
      }
      md <- s$md
      if (!is.data.frame(md) || nrow(md) == 0) {
        return(data.frame(Message = "No data available", stringsAsFactors = FALSE))
      }
      sp <- aggregate(
        cbind(Sample_Size, Missing_Data, Genotyped_Data) ~ Population,
        data = md,
        FUN  = sum
      )
      sp$Missing_Proportion <- round(sp$Missing_Data / sp$Sample_Size, 3)
      sp <- sp[, c("Population","Sample_Size","Genotyped_Data",
                   "Missing_Data","Missing_Proportion")]
      # Force integers for counts, keep proportion as numeric
      sp[["Sample_Size"]]    <- as.integer(sp[["Sample_Size"]])
      sp[["Genotyped_Data"]] <- as.integer(sp[["Genotyped_Data"]])
      sp[["Missing_Data"]]   <- as.integer(sp[["Missing_Data"]])
      names(sp) <- c("Population", "N", "Genotyped", "Missing", "Miss. %")
      sp
    }, striped = TRUE, hover = TRUE, spacing = "xs", width = "100%",
       digits = 3, colnames = TRUE, align = "lrrrr")

    output$data_summary_by_pop_locus <- DT::renderDT({
      s <- req(summary_trigger_r()); md <- s$md
      if(nrow(md)==0) return(DT::datatable(data.frame(Message="No data available")))
      DT::datatable(
        md[,c("Population","Marker","Sample_Size","Genotyped_Data",
              "Missing_Data","Missing_Proportion")],
        options=list(pageLength=15,scrollX=TRUE), rownames=FALSE,
        colnames=c("Population","Marker","Sample size","Genotyped",
                   "Missing","Missing proportion")
      ) %>% DT::formatRound("Missing_Proportion",4)
    })

    output$data_summary_by_locus_mean <- DT::renderDT({
      s <- req(summary_trigger_r()); md <- s$md
      if(nrow(md)==0) return(DT::datatable(data.frame(Message="No data available")))
      md$Marker <- factor(md$Marker,levels=markers_r())
      sl <- aggregate(cbind(Sample_Size,Missing_Data,Genotyped_Data)~Marker,data=md,FUN=mean)
      sl <- sl[order(sl$Marker),]; sl$Marker <- as.character(sl$Marker)
      sl$Missing_Proportion <- sl$Missing_Data/sl$Sample_Size
      DT::datatable(
        sl[,c("Marker","Sample_Size","Genotyped_Data","Missing_Data","Missing_Proportion")],
        options=list(pageLength=15,scrollX=TRUE), rownames=FALSE,
        colnames=c("Marker","Avg sample size","Avg genotyped",
                   "Avg missing","Avg missing proportion")
      ) %>% DT::formatRound(c("Sample_Size","Genotyped_Data",
                               "Missing_Data","Missing_Proportion"),3)
    })

    subsamples_data <- reactive({
      md <- missing_by_pop_locus_r()
      if(!identical(selected_population_subsamples_r(),"all"))
        md <- md[md$Population==selected_population_subsamples_r(),,drop=FALSE]
      md
    })

    output$data_summary_by_Subsamples_locus_sum <- DT::renderDT({
      val <- subsamples_data()
      if(nrow(val)==0)
        return(DT::datatable(data.frame(Message="No data for selected population")))
      val$Marker <- factor(val$Marker,levels=markers_r())
      sl <- aggregate(cbind(Sample_Size,Missing_Data,Genotyped_Data)~Marker,data=val,FUN=sum)
      sl <- sl[order(sl$Marker),]; sl$Marker <- as.character(sl$Marker)
      sl$Missing_Proportion <- sl$Missing_Data/sl$Sample_Size
      has_pop <- !identical(selected_population_subsamples_r(),"all")
      if(has_pop){
        sl$Population <- selected_population_subsamples_r()
        sl <- sl[,c("Population","Marker","Sample_Size","Genotyped_Data",
                    "Missing_Data","Missing_Proportion")]
      }
      DT::datatable(sl,options=list(pageLength=15,scrollX=TRUE),rownames=FALSE,
        colnames=if(has_pop)
          c("Population","Locus","Total sample size","Total genotyped",
            "Total missing","Proportion missing")
        else c("Locus","Total sample size","Total genotyped",
               "Total missing","Proportion missing")
      ) %>% DT::formatRound(c("Sample_Size","Genotyped_Data",
                               "Missing_Data","Missing_Proportion"),3)
    })

    output$data_summary_by_locus_sum <- DT::renderDT({
      s <- req(summary_trigger_r()); md <- s$md
      if(nrow(md)==0) return(DT::datatable(data.frame(Message="No data available")))
      md$Marker <- factor(md$Marker,levels=markers_r())
      sl <- aggregate(cbind(Sample_Size,Missing_Data,Genotyped_Data)~Marker,data=md,FUN=sum)
      sl <- sl[order(sl$Marker),]; sl$Marker <- as.character(sl$Marker)
      sl$Missing_Proportion <- sl$Missing_Data/sl$Sample_Size
      DT::datatable(
        sl[,c("Marker","Sample_Size","Genotyped_Data","Missing_Data","Missing_Proportion")],
        options=list(pageLength=15,scrollX=TRUE),rownames=FALSE,
        colnames=c("Locus","Total sample size","Total genotyped",
                   "Total missing","Proportion missing")
      ) %>% DT::formatRound(c("Sample_Size","Genotyped_Data",
                               "Missing_Data","Missing_Proportion"),3)
    })

    output$data_summary_by_pop_sum <- DT::renderDT({
      s <- req(summary_trigger_r()); md <- s$md
      if(nrow(md)==0) return(DT::datatable(data.frame(Message="No data available")))
      sp <- aggregate(cbind(Sample_Size,Missing_Data,Genotyped_Data)~Population,
                      data=md,FUN=sum)
      sp$Missing_Proportion <- sp$Missing_Data/sp$Sample_Size
      DT::datatable(
        sp[,c("Population","Sample_Size","Genotyped_Data","Missing_Data","Missing_Proportion")],
        options=list(pageLength=15,scrollX=TRUE),rownames=FALSE,
        colnames=c("Population","Total sample size","Total genotyped",
                   "Total missing","Proportion missing")
      ) %>% DT::formatRound(c("Sample_Size","Genotyped_Data",
                               "Missing_Data","Missing_Proportion"),3)
    })

    output$data_summary_global <- DT::renderDT({
      s <- req(summary_trigger_r()); md <- s$md
      if(nrow(md)==0) return(DT::datatable(data.frame(Message="No data available")))
      DT::datatable(data.frame(
        Total_Sample_Size  = sum(md$Sample_Size),
        Total_Genotyped    = sum(md$Genotyped_Data),
        Total_Missing      = sum(md$Missing_Data),
        Proportion_Missing = round(sum(md$Missing_Data)/sum(md$Sample_Size),4)),
        options=list(dom="t"),rownames=FALSE,
        colnames=c("Total sample size","Total genotyped",
                   "Total missing","Proportion missing"))
    })

    # ── Existing frequency DTs ─────────────────────────────────────────────
    output$allele_freq_by_pop <- DT::renderDT({
      d <- af_table_by_pop_r()
      shiny::validate(shiny::need(nrow(d) > 0, "No data available for selected filters."))
      DT::datatable(
        d,
        options  = list(pageLength = 15, scrollX = TRUE, dom = "lftip"),
        rownames = FALSE
      ) %>%
        DT::formatStyle("Frequency",
          backgroundColor = DT::styleInterval(
            c(0.1, 0.5), c("#ffebee","#fff3e0","#e8f5e8")))
    }, server = TRUE)

    output$allele_freq_global <- DT::renderDT({
      d <- af_table_global_r()
      shiny::validate(shiny::need(nrow(d) > 0, "No data available for selected filters."))
      DT::datatable(
        d,
        options  = list(pageLength = 15, scrollX = TRUE, dom = "lftip"),
        rownames = FALSE
      ) %>%
        DT::formatStyle("Frequency",
          backgroundColor = DT::styleInterval(
            c(0.1, 0.5), c("#ffebee","#fff3e0","#e8f5e8")))
    }, server = TRUE)

    output$quick_stats_pop <- renderText({
      analysis_ready_r()
      ds <- run_af_by_pop(selected_population_r(),"all")
      pn <- if(selected_population_r()!="all") selected_population_r() else "All populations"
      if(nrow(ds)>0&&!all(is.na(ds$Frequency)))
        paste0("Population: ",pn,"\nTotal alleles: ",nrow(ds),
               "\nAverage frequency: ",round(mean(ds$Frequency,na.rm=TRUE),3),
               "\nMaximum frequency: ",round(max(ds$Frequency,na.rm=TRUE),3),
               "\nMinimum frequency: ",round(min(ds$Frequency,na.rm=TRUE),3))
      else "No data available"
    })

    output$quick_stats_global <- renderText({
      analysis_ready_r()
      ds <- run_af_global("all")
      if(nrow(ds)>0&&!all(is.na(ds$Frequency)))
        paste0("Global analysis\nTotal alleles: ",nrow(ds),
               "\nPolymorphic markers: ",length(unique(ds$Marker)),
               "\nAverage frequency: ",round(mean(ds$Frequency,na.rm=TRUE),3))
      else "No data available"
    })

    # ── Plots ──────────────────────────────────────────────────────────────
    .zissou1 <- c("#3B9AB2","#78B7C5","#EBCC2A","#E1AF00","#F21A00")
    .zpal <- function(n) colorRampPalette(.zissou1)(max(n, 5L))

    output$allele_freq_plot <- renderPlot({
      analysis_ready_r(); pd <- af_table_by_pop_r()
      shiny::validate(shiny::need(nrow(pd)>0,"No data available."))
      pd$Marker <- factor(pd$Marker,levels=markers_r())
      sm <- selected_marker_r()
      grps <- length(unique(interaction(pd$Population, pd$Marker)))
      ggplot(pd,aes(x=Allele,y=Frequency,fill=interaction(Population,Marker)))+
        geom_col(position="dodge")+
        facet_wrap(~Marker,scales="free_x")+
        labs(title=paste("Allele frequencies",if(sm!="all") paste("- Marker:",sm) else "- All markers"),
             x="Alleles",y="Frequency",fill="Population-Marker")+
        theme_minimal()+theme(axis.text.x=element_text(angle=45,hjust=1))+
        scale_fill_manual(values=.zpal(grps))
    })

    output$missing_data_plot <- renderPlot({
      md <- missing_by_pop_locus_r()
      if(nrow(md)>0){
        md$Marker <- factor(md$Marker,levels=markers_r())
        n_pops <- length(unique(md$Population))
        ggplot(md,aes(x=Marker,y=Missing_Proportion,fill=Population))+
          geom_col(position="dodge")+
          labs(title="Proportion of missing data",x="Markers",y="Proportion missing")+
          theme_minimal()+theme(axis.text.x=element_text(angle=45,hjust=1))+
          scale_fill_manual(values=.zpal(n_pops))+
          geom_hline(yintercept=0.20,linetype="dashed",color="#E1AF00",alpha=0.8)
      } else {
        ggplot()+annotate("text",x=1,y=1,label="No missing data available",size=6)+theme_void()
      }
    })

    output$allele_richness_plot <- renderPlot({
      dd <- genetic_diversity_r()
      shiny::validate(shiny::need(nrow(dd)>0,"No diversity data available."))
      dd$Marker <- factor(dd$Marker,levels=markers_r())
      n_mrk <- length(unique(dd$Marker))
      ggplot(dd,aes(x=Population,y=Na,fill=Marker))+
        geom_col(position="dodge")+
        labs(title="Allelic richness by population and marker",
             x="Population",y="Number of alleles (Na)")+
        theme_minimal()+theme(axis.text.x=element_text(angle=45,hjust=1))+
        scale_fill_manual(values=.zpal(n_mrk))
    })

    # ── Marker details DT ──────────────────────────────────────────────────
    output$marker_details <- DT::renderDT({
      ms <- marker_details_r()
      shiny::validate(shiny::need(nrow(ms)>0,"No data available."))
      DT::datatable(ms,options=list(pageLength=10,scrollX=TRUE),rownames=FALSE,
        colnames=c("Marker","Number of alleles","Total observations",
                   "Average missing proportion","Observed alleles")
      ) %>% DT::formatStyle("Avg_Missing_Proportion",
        backgroundColor=DT::styleInterval(c(0.10,0.20),c("#e8f5e8","#fff3e0","#ffebee")))
    },server=TRUE)

    # ── Complete matrix ────────────────────────────────────────────────────
    complete_allele_matrix <- eventReactive(input$generate_complete_matrix,{
      db_ready()
      showNotification("Generating complete allele matrix...",
        type="message",duration=NULL,id="matrix_calc")
      on.exit(removeNotification(id="matrix_calc"),add=TRUE)
      con <- con_r(); hs <- hf_schema_r(); ms <- meta_schema_r()
      base <- as.integer(base_r())
      hf_q <- sql_id(con,tbl_hf_r()); meta_q <- sql_id(con,tbl_meta_r())
      hi_q <- sql_id(con,hs$ind_col); hl_q <- sql_id(con,hs$locus_col)
      hg_q <- sql_id(con,hs$gt_col);  mi_q <- sql_id(con,ms$ind_col)
      pop_q <- sql_id(con,ms$pop_col)
      out <- DBI::dbGetQuery(con, sprintf("
    WITH %s,
    valid_gt AS (
      SELECT CAST(m.%s AS VARCHAR) AS Population,
             CAST(h.%s AS VARCHAR) AS Marker, h.%s AS gt
      FROM %s h INNER JOIN %s m ON CAST(h.%s AS VARCHAR)=CAST(m.%s AS VARCHAR)
      WHERE m.%s IS NOT NULL AND h.%s IS NOT NULL AND h.%s>0),
    populations AS (SELECT DISTINCT Population FROM valid_gt),
    alleles AS (
      SELECT Population,Marker,CAST(FLOOR(gt/%d) AS BIGINT) AS Allele FROM valid_gt
      UNION ALL
      SELECT Population,Marker,CAST(gt%%%d AS BIGINT) AS Allele FROM valid_gt),
    alleles_by_marker AS (SELECT DISTINCT Marker,Allele FROM alleles),
    pop_marker_grid AS (
      SELECT p.Population,a.Marker,a.Allele
      FROM populations p CROSS JOIN alleles_by_marker a),
    counts AS (SELECT Population,Marker,Allele,COUNT(*) AS Count
      FROM alleles GROUP BY Population,Marker,Allele),
    totals AS (SELECT Population,Marker,COUNT(*) AS N_total
      FROM alleles GROUP BY Population,Marker)
    SELECT g.Population,g.Marker,CAST(g.Allele AS VARCHAR) AS Allele,
      COALESCE(c.Count*1.0/NULLIF(t.N_total,0),0.0) AS Frequency,
      COALESCE(c.Count,0) AS Count,COALESCE(t.N_total,0) AS N_total,
      CASE WHEN COALESCE(c.Count,0)>0 THEN 'Present' ELSE 'Absent' END AS Status
    FROM pop_marker_grid g
    LEFT JOIN counts c ON g.Population=c.Population AND g.Marker=c.Marker AND g.Allele=c.Allele
    LEFT JOIN totals t ON g.Population=t.Population AND g.Marker=t.Marker
    LEFT JOIN locus_order lo ON g.Marker=lo._lo_marker
    ORDER BY lo._lo_rank ASC,g.Population,g.Allele",
        locus_order_cte(con,hf_q,hl_q),
        pop_q,hl_q,hg_q, hf_q,meta_q,hi_q,mi_q,
        pop_q,hg_q,hg_q, base,base))
      showNotification(paste0("Matrix generated: ",nrow(out)," rows"),
        type="message",duration=5)
      out
    },ignoreInit=TRUE)

    output$complete_freq_matrix <- DT::renderDT({
      m <- complete_allele_matrix()
      if(is.null(m)) return(DT::datatable(
        data.frame(Message="Click 'Generate Complete Matrix' to start"),
        options=list(dom="t"),rownames=FALSE,colnames=""))
      m$Frequency <- round(m$Frequency,4)
      DT::datatable(m,
        options=list(pageLength=20,scrollX=TRUE,order=list(list(0,"asc"))),
        rownames=FALSE,filter="top",
        colnames=c("Population","Locus","Allele","Frequency","Count","N total","Status")
      ) %>%
        DT::formatStyle("Frequency",backgroundColor=DT::styleInterval(
          c(0.001,0.1,0.5),c("#f5f5f5","#ffebee","#fff3e0","#e8f5e8"))) %>%
        DT::formatStyle("Status",backgroundColor=DT::styleEqual(
          c("Present","Absent"),c("#d4edda","#f8d7da")))
    })

    # ── Downloads ──────────────────────────────────────────────────────────
    output$download_fstat_csv <- downloadHandler(
      filename=function() paste0("allele_freq_by_locus_",Sys.Date(),".csv"),
      content=function(file){
        wide <- fstat_wide_r(); if(is.null(wide)) return(NULL)
        pops <- attr(wide,"pops")
        write.csv(wide[,c("Locus","Row_label",pops,"Global")],file,row.names=FALSE)
      })
    output$download_fstat_txt <- downloadHandler(
      filename=function() paste0("allele_freq_by_locus_",Sys.Date(),".txt"),
      content=function(file){
        wide <- fstat_wide_r(); if(is.null(wide)) return(NULL)
        pops <- attr(wide,"pops")
        write.table(wide[,c("Locus","Row_label",pops,"Global")],
                    file,sep="\t",row.names=FALSE,quote=FALSE)
      })
    output$download_freq_pop_csv <- downloadHandler(
      filename=function() paste0("frequencies_by_population_",Sys.Date(),".csv"),
      content=function(file) write.csv(run_af_by_pop("all","all"),file,row.names=FALSE))
    output$download_freq_pop_txt <- downloadHandler(
      filename=function() paste0("frequencies_by_population_",Sys.Date(),".txt"),
      content=function(file) write.table(run_af_by_pop("all","all"),file,sep="\t",row.names=FALSE,quote=FALSE))
    output$download_freq_global_csv <- downloadHandler(
      filename=function() paste0("frequencies_global_",Sys.Date(),".csv"),
      content=function(file) write.csv(run_af_global("all"),file,row.names=FALSE))
    output$download_freq_global_txt <- downloadHandler(
      filename=function() paste0("frequencies_global_",Sys.Date(),".txt"),
      content=function(file) write.table(run_af_global("all"),file,sep="\t",row.names=FALSE,quote=FALSE))
    output$download_summary_pop_locus <- downloadHandler(
      filename=function() paste0("summary_pop_marker_",Sys.Date(),".csv"),
      content=function(file) write.csv(missing_by_pop_locus_r(),file,row.names=FALSE))
    output$download_summary_pop_locus_txt <- downloadHandler(
      filename=function() paste0("summary_pop_marker_",Sys.Date(),".txt"),
      content=function(file) write.table(missing_by_pop_locus_r(),file,sep="\t",row.names=FALSE,quote=FALSE))
    output$download_summary_locus_mean <- downloadHandler(
      filename=function() paste0("summary_marker_mean_",Sys.Date(),".csv"),
      content=function(file){
        md <- missing_by_pop_locus_r(); md$Marker <- factor(md$Marker,levels=markers_r())
        sl <- aggregate(cbind(Sample_Size,Missing_Data,Genotyped_Data)~Marker,data=md,FUN=mean)
        sl <- sl[order(sl$Marker),]; sl$Marker <- as.character(sl$Marker)
        sl$Missing_Proportion <- sl$Missing_Data/sl$Sample_Size
        write.csv(sl,file,row.names=FALSE)})
    output$download_summary_locus_mean_txt <- downloadHandler(
      filename=function() paste0("summary_marker_mean_",Sys.Date(),".txt"),
      content=function(file){
        md <- missing_by_pop_locus_r(); md$Marker <- factor(md$Marker,levels=markers_r())
        sl <- aggregate(cbind(Sample_Size,Missing_Data,Genotyped_Data)~Marker,data=md,FUN=mean)
        sl <- sl[order(sl$Marker),]; sl$Marker <- as.character(sl$Marker)
        sl$Missing_Proportion <- sl$Missing_Data/sl$Sample_Size
        write.table(sl,file,sep="\t",row.names=FALSE,quote=FALSE)})
    output$download_summary_Subsamples_locus_sum <- downloadHandler(
      filename=function(){
        suf <- if(!identical(selected_population_subsamples_r(),"all"))
          paste0("_",selected_population_subsamples_r()) else "_all"
        paste0("summary_subsamples",suf,"_",Sys.Date(),".csv")},
      content=function(file){
        val <- subsamples_data(); if(nrow(val)==0){write.csv(data.frame(Message="No data"),file,row.names=FALSE);return()}
        val$Marker <- factor(val$Marker,levels=markers_r())
        sl <- aggregate(cbind(Sample_Size,Missing_Data,Genotyped_Data)~Marker,data=val,FUN=sum)
        sl <- sl[order(sl$Marker),]; sl$Marker <- as.character(sl$Marker)
        sl$Missing_Proportion <- sl$Missing_Data/sl$Sample_Size
        write.csv(sl,file,row.names=FALSE)})
    output$download_summary_Subsamples_locus_sum_txt <- downloadHandler(
      filename=function(){
        suf <- if(!identical(selected_population_subsamples_r(),"all"))
          paste0("_",selected_population_subsamples_r()) else "_all"
        paste0("summary_subsamples",suf,"_",Sys.Date(),".txt")},
      content=function(file){
        val <- subsamples_data(); if(nrow(val)==0){write.table(data.frame(Message="No data"),file,sep="\t",row.names=FALSE,quote=FALSE);return()}
        val$Marker <- factor(val$Marker,levels=markers_r())
        sl <- aggregate(cbind(Sample_Size,Missing_Data,Genotyped_Data)~Marker,data=val,FUN=sum)
        sl <- sl[order(sl$Marker),]; sl$Marker <- as.character(sl$Marker)
        sl$Missing_Proportion <- sl$Missing_Data/sl$Sample_Size
        write.table(sl,file,sep="\t",row.names=FALSE,quote=FALSE)})
    output$download_summary_locus_sum <- downloadHandler(
      filename=function() paste0("summary_marker_sum_",Sys.Date(),".csv"),
      content=function(file){
        md <- missing_by_pop_locus_r(); md$Marker <- factor(md$Marker,levels=markers_r())
        sl <- aggregate(cbind(Sample_Size,Missing_Data,Genotyped_Data)~Marker,data=md,FUN=sum)
        sl <- sl[order(sl$Marker),]; sl$Marker <- as.character(sl$Marker)
        sl$Missing_Proportion <- sl$Missing_Data/sl$Sample_Size
        write.csv(sl,file,row.names=FALSE)})
    output$download_summary_locus_sum_txt <- downloadHandler(
      filename=function() paste0("summary_marker_sum_",Sys.Date(),".txt"),
      content=function(file){
        md <- missing_by_pop_locus_r(); md$Marker <- factor(md$Marker,levels=markers_r())
        sl <- aggregate(cbind(Sample_Size,Missing_Data,Genotyped_Data)~Marker,data=md,FUN=sum)
        sl <- sl[order(sl$Marker),]; sl$Marker <- as.character(sl$Marker)
        sl$Missing_Proportion <- sl$Missing_Data/sl$Sample_Size
        write.table(sl,file,sep="\t",row.names=FALSE,quote=FALSE)})
    output$download_summary_pop_sum <- downloadHandler(
      filename=function() paste0("summary_population_sum_",Sys.Date(),".csv"),
      content=function(file){
        md <- missing_by_pop_locus_r()
        sp <- aggregate(cbind(Sample_Size,Missing_Data,Genotyped_Data)~Population,data=md,FUN=sum)
        sp$Missing_Proportion <- sp$Missing_Data/sp$Sample_Size
        write.csv(sp,file,row.names=FALSE)})
    output$download_summary_pop_sum_txt <- downloadHandler(
      filename=function() paste0("summary_population_sum_",Sys.Date(),".txt"),
      content=function(file){
        md <- missing_by_pop_locus_r()
        sp <- aggregate(cbind(Sample_Size,Missing_Data,Genotyped_Data)~Population,data=md,FUN=sum)
        sp$Missing_Proportion <- sp$Missing_Data/sp$Sample_Size
        write.table(sp,file,sep="\t",row.names=FALSE,quote=FALSE)})
    output$download_summary_global <- downloadHandler(
      filename=function() paste0("summary_global_",Sys.Date(),".csv"),
      content=function(file){
        md <- missing_by_pop_locus_r()
        shiny::validate(shiny::need(nrow(md)>0,"No data"))
        write.csv(data.frame(Total_Sample_Size=sum(md$Sample_Size),
          Total_Genotyped=sum(md$Genotyped_Data),Total_Missing=sum(md$Missing_Data),
          Proportion_Missing=round(sum(md$Missing_Data)/sum(md$Sample_Size),4)),
          file,row.names=FALSE)})
    output$download_summary_global_txt <- downloadHandler(
      filename=function() paste0("summary_global_",Sys.Date(),".txt"),
      content=function(file){
        md <- missing_by_pop_locus_r()
        validate(need(nrow(md)>0,"No data"))
        write.table(data.frame(Total_Sample_Size=sum(md$Sample_Size),
          Total_Genotyped=sum(md$Genotyped_Data),Total_Missing=sum(md$Missing_Data),
          Proportion_Missing=round(sum(md$Missing_Data)/sum(md$Sample_Size),4)),
          file,sep="\t",row.names=FALSE,quote=FALSE)})
    output$download_complete_matrix_csv <- downloadHandler(
      filename=function() paste0("complete_allele_matrix_",Sys.Date(),".csv"),
      content=function(file){
        m <- complete_allele_matrix()
        if(is.null(m)){showNotification("Generate matrix first!",type="error");return(NULL)}
        write.csv(m,file,row.names=FALSE)})
    output$download_complete_matrix_txt <- downloadHandler(
      filename=function() paste0("complete_allele_matrix_",Sys.Date(),".txt"),
      content=function(file){
        m <- complete_allele_matrix()
        if(is.null(m)){showNotification("Generate matrix first!",type="error");return(NULL)}
        write.table(m,file,sep="\t",row.names=FALSE,quote=FALSE)})

    output$download_allele_plots <- downloadHandler(
      filename = function() paste0("allele_plots_", Sys.Date(), ".png"),
      content = function(file) {
        req(analysis_ready_r())
        pd <- af_table_by_pop_r()
        md <- missing_by_pop_locus_r()
        dd <- genetic_diversity_r()
        zp <- function(n) colorRampPalette(c("#3B9AB2","#78B7C5","#EBCC2A","#E1AF00","#F21A00"))(max(n,5L))

        p1 <- if (nrow(pd) > 0) {
          pd$Marker <- factor(pd$Marker, levels = markers_r())
          ggplot2::ggplot(pd, ggplot2::aes(x=Allele, y=Frequency,
                           fill=interaction(Population,Marker))) +
            ggplot2::geom_col(position="dodge") +
            ggplot2::facet_wrap(~Marker, scales="free_x") +
            ggplot2::labs(title="Allele frequencies", x="Allele", y="Frequency", fill="Pop-Marker") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1)) +
            ggplot2::scale_fill_manual(values=zp(length(unique(interaction(pd$Population,pd$Marker)))))
        } else ggplot2::ggplot() + ggplot2::annotate("text",x=1,y=1,label="No data") + ggplot2::theme_void()

        p2 <- if (nrow(md) > 0) {
          md$Marker <- factor(md$Marker, levels = markers_r())
          ggplot2::ggplot(md, ggplot2::aes(x=Marker,y=Missing_Proportion,fill=Population)) +
            ggplot2::geom_col(position="dodge") +
            ggplot2::labs(title="Missing data proportion", x="Markers", y="Proportion missing") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1)) +
            ggplot2::scale_fill_manual(values=zp(length(unique(md$Population)))) +
            ggplot2::geom_hline(yintercept=0.20, linetype="dashed", color="#E1AF00", alpha=0.8)
        } else ggplot2::ggplot() + ggplot2::annotate("text",x=1,y=1,label="No missing data") + ggplot2::theme_void()

        p3 <- if (nrow(dd) > 0) {
          dd$Marker <- factor(dd$Marker, levels = markers_r())
          ggplot2::ggplot(dd, ggplot2::aes(x=Population,y=Na,fill=Marker)) +
            ggplot2::geom_col(position="dodge") +
            ggplot2::labs(title="Allelic richness", x="Population", y="Number of alleles (Na)") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1)) +
            ggplot2::scale_fill_manual(values=zp(length(unique(dd$Marker))))
        } else ggplot2::ggplot() + ggplot2::annotate("text",x=1,y=1,label="No data") + ggplot2::theme_void()

        png(file, width=1400, height=1800, res=120)
        gridExtra::grid.arrange(p1, p2, p3, ncol=1)
        dev.off()
      }
    )

    outputOptions(output,"allele_freq_plot",    suspendWhenHidden=FALSE)
    outputOptions(output,"missing_data_plot",   suspendWhenHidden=FALSE)
    outputOptions(output,"allele_richness_plot",suspendWhenHidden=FALSE)
  })
}