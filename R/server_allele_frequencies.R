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

    # ── Locus order ────────────────────────────────────────────────────────
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
      updateSelectInput(session, "fstat_population",
        choices  = c("All populations"="all", stats::setNames(pops,pops)),
        selected = "all")
      updateSelectizeInput(session, "fstat_marker",
        choices  = c("All markers"="all", stats::setNames(markers,markers)),
        selected = "all", server=TRUE)
    })

    fstat_population_r <- reactive(safe_choice(input$fstat_population, "all"))
    fstat_marker_r     <- reactive(safe_choice(input$fstat_marker,     "all"))

    fstat_ready_r <- reactive({
      req(input$update_fstat > 0L); db_ready(); TRUE })

    # ── Missing data (value box only) ──────────────────────────────────────
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

    # ── Fstat long reactive ────────────────────────────────────────────────
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

      gc  <- aggregate(Count~Marker+Allele, data=long, FUN=sum)
      gt  <- aggregate(Count~Marker,        data=long, FUN=sum)
      names(gt)[2] <- "Total"
      gc  <- merge(gc,gt,by="Marker")
      gc$Global_Freq <- gc$Count/gc$Total

      rows <- list()

      for (loc in markers) {
        lng     <- long[long$Marker==loc,,drop=FALSE]
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

        ng <- stats::setNames(lapply(pops, function(p) {
          x <- lng$N_genotyped[lng$Population==p]; if(length(x)) x[1] else NA}), pops)
        gng <- sum(sapply(ng, function(v) if(is.null(v)||is.na(v)) 0 else v))
        r <- mk_row("stat","N genotyped",ng,gng)
        r["Locus"] <- loc
        rows[[length(rows)+1]] <- r

        nm <- stats::setNames(lapply(pops, function(p) {
          x <- lng$N_missing[lng$Population==p]; if(length(x)) x[1] else NA}), pops)
        gnm <- sum(sapply(nm, function(v) if(is.null(v)||is.na(v)) 0 else v))
        rows[[length(rows)+1]] <- mk_row("stat","N missing",nm,gnm)

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

      display    <- wide[,c("Locus","Row_label",pops,"Global"),drop=FALSE]
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

    # ── Downloads fstat ────────────────────────────────────────────────────
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
  })
}