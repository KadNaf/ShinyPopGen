# module/server_null_alleles.R

server_null_alleles <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    # ── Helpers ────────────────────────────────────────────────────────────
    `%||%` <- function(a, b) if (!is.null(a)) a else b

    safe_choice <- function(x, default = "all") {
      if (is.null(x) || length(x) == 0L || identical(x, "") || all(is.na(x))) default
      else as.character(x[[1]])
    }

    sql_id  <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))
    sql_str <- function(con, x) as.character(DBI::dbQuoteString(con, x))

    # Stable input ID for a locus treatment selector (no spaces / special chars)
    treat_id <- function(loc) paste0("treat_", gsub("[^A-Za-z0-9]", "_", loc))

    # ── DB plumbing ────────────────────────────────────────────────────────
    db_tick    <- reactive({ rv$db_tick })
    con_r      <- reactive({ req(rv$con); rv$con })
    tbl_meta_r <- reactive({ rv$tbl_meta %||% "meta" })

    tbl_hf_r <- reactive({
      con <- con_r()
      if (exists("duck_tbl_exists",    mode = "function", inherits = TRUE) &&
          exists(".duckdb_get_params", mode = "function", inherits = TRUE) &&
          duck_tbl_exists(con, "params")) {
        p <- .duckdb_get_params(con)
        return(as.character(p$tbl_hf %||% "hf"))
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
                 else shiny::validate(shiny::need(FALSE, "No individual column in meta."))
      pop_col <- c("Population","population","pop","pop_code")[
        c("Population","population","pop","pop_code") %in% cols][1]
      shiny::validate(shiny::need(!is.na(pop_col), "No population column in meta."))
      list(ind_col=ind_col, pop_col=pop_col)
    })

    # Every query using this CTE must also add:
    #   LEFT JOIN locus_order lo ON CAST(h.<locus_col> AS VARCHAR) = lo._lo_marker
    locus_order_cte <- function(con, hf_tbl_q, hl_q)
      sprintf("locus_order AS (
  SELECT CAST(%s AS VARCHAR) AS _lo_marker, MIN(rowid) AS _lo_rank
  FROM %s GROUP BY CAST(%s AS VARCHAR))", hl_q, hf_tbl_q, hl_q)

    # ── Marker / population lists ──────────────────────────────────────────
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
        FROM %s h
        LEFT JOIN locus_order lo ON CAST(%s AS VARCHAR) = lo._lo_marker
        ORDER BY lo._lo_rank ASC",
        locus_order_cte(con,hf_q,hl_q), hl_q, hf_q, hl_q))$Marker)
    })

    observe({
      markers <- markers_r(); pops <- pops_r()
      updateSelectInput(session, "t1_locus",
        choices  = c("All loci"="all", stats::setNames(markers,markers)),
        selected = "all")
      updateSelectInput(session, "t1_pop",
        choices  = c("All populations"="all", stats::setNames(pops,pops)),
        selected = "all")
      updateSelectInput(session, "t2_locus",
        choices  = c("All loci"="all", stats::setNames(markers,markers)),
        selected = "all")
    })

    # ── Per-locus treatment selector UI ───────────────────────────────────
    # KEY FIX: use session$ns() inside renderUI to namespace dynamically generated IDs
    output$locus_treatment_ui <- renderUI({
      ns_fn   <- session$ns          # session$ns() is the correct way inside renderUI
      markers <- markers_r()
      if (length(markers) == 0L) return(tags$p("No markers loaded yet."))

      items <- lapply(markers, function(loc) {
        tags$div(class = "na-treat-item",
          tags$div(class = "na-treat-lbl", loc),
          selectInput(
            inputId  = ns_fn(treat_id(loc)),   # correctly namespaced
            label    = NULL,
            choices  = c(
              "999999 \u2014 null homozygote"          = "null_homo",
              "000000 \u2014 absent / PCR failure"     = "absent"
            ),
            selected = "null_homo",
            width    = "100%"
          )
        )
      })
      tags$div(class = "na-treat-grid", items)
    })

    # Named vector: locus_name → "null_homo" | "absent"
    # Inside moduleServer, input$<id> does NOT need ns() — Shiny handles it automatically
    locus_treatments_r <- reactive({
      markers <- markers_r()
      treats  <- sapply(markers, function(loc) {
        val <- input[[treat_id(loc)]]    # no ns() needed inside moduleServer
        if (is.null(val) || !val %in% c("null_homo","absent")) "null_homo" else val
      })
      stats::setNames(treats, markers)
    })

    # ══════════════════════════════════════════════════════════════════════
    # EM NULL-HOMO model — Chapuis & Estoup (2007) / FreeNA exact
    #
    # Missing (gt=0|NA) = H_00 = null homozygotes (coded 999999 in Genepop).
    # N = efpop (total individuals — null homos included in denominator).
    # r update: r_new = Σ_a [r/(p_a + 2r)] * (H_aa/N)  +  H_00/N
    # ══════════════════════════════════════════════════════════════════════
    em_null_homo <- function(gt_vec, base) {
      efpop     <- length(gt_vec)
      null_mask <- is.na(gt_vec) | gt_vec <= 0L
      H_00      <- sum(null_mask)
      N         <- efpop

      if (N == 0L) return(list(rd=0.0, n=efpop))
      valid_gt <- gt_vec[!null_mask]
      if (length(valid_gt) == 0L)
        return(list(rd=ifelse(N>0L, H_00/N, 0.0), n=efpop))

      a1 <- floor(valid_gt / base)
      a2 <- valid_gt %% base
      all_alleles <- sort(unique(c(a1, a2)))
      all_alleles <- all_alleles[all_alleles >= 0L]
      if (length(all_alleles) == 0L) return(list(rd=0.0, n=efpop))

      n_valid  <- N - H_00
      genefreq <- sapply(all_alleles, function(a)
        (sum(a1==a) + sum(a2==a)) / (2L * n_valid))

      r    <- if (H_00 > 0L) sqrt(H_00 / N) else sqrt(1.0 / (N + 1.0))
      H_ii <- sapply(all_alleles, function(a) sum(a1==a & a2==a))
      H_iX <- sapply(all_alleles, function(a) sum((a1==a & a2!=a)|(a2==a & a1!=a)))
      hotot <- sum(H_ii)

      # cpt=0 initialisation — exact FreeNA Pascal formula
      p <- numeric(length(all_alleles))
      for (ai in seq_along(all_alleles)) {
        if (genefreq[ai] <= 0) { p[ai] <- 0.0; next }
        ii <- H_ii[ai]; jj <- H_iX[ai]
        if (H_00 > 0L) {
          X <- H_00 + hotot - ii + (N - H_00 - hotot) - jj; Y <- N
        } else {
          X <- 1.0 + hotot - ii + (N - hotot) - jj;         Y <- N + 1.0
        }
        p[ai] <- 1.0 - sqrt(max(0.0, X / Y))
      }

      # EM iterations — denominator = N (total)
      for (iter in seq_len(5000L)) {
        new_p <- numeric(length(all_alleles))
        ri <- 0.0; re <- 0L
        for (ai in seq_along(all_alleles)) {
          if (genefreq[ai] <= 0) { new_p[ai] <- 0.0; next }
          pa <- p[ai]; denom <- pa + 2.0 * r
          if (denom <= 0) { new_p[ai] <- 0.0; next }
          p_new     <- (pa + r)/denom * (H_ii[ai]/N) + H_iX[ai]/(2.0*N)
          ri        <- ri + r/denom * (H_ii[ai]/N)
          new_p[ai] <- p_new
          if (abs(p_new - pa) > 1e-6) re <- re + 1L
        }
        r_new <- ri + H_00 / N
        if (abs(r_new - r) > 1e-6) re <- re + 1L
        p <- new_p; r <- max(0.0, r_new)
        if (re == 0L) break
      }
      list(rd=r, n=efpop)
    }

    # ══════════════════════════════════════════════════════════════════════
    # EM ABSENT model — Chapuis & Estoup (2007) / FreeNA exact
    #
    # Missing (gt=0|NA) = absent / PCR failure (coded 000000 in Genepop).
    # N = n_valid (valid genotypes only — absents excluded from denominator).
    # r update: r_new = Σ_a [r/(p_a + 2r)] * (H_aa/N)   (no null homo term)
    # ══════════════════════════════════════════════════════════════════════
    em_absent <- function(gt_vec, base) {
      efpop     <- length(gt_vec)
      null_mask <- is.na(gt_vec) | gt_vec <= 0L
      n_absent  <- sum(null_mask)
      N         <- efpop - n_absent          # valid genotypes only

      if (N == 0L) return(list(rd=0.0, n=efpop))
      valid_gt <- gt_vec[!null_mask]
      if (length(valid_gt) == 0L) return(list(rd=0.0, n=efpop))

      a1 <- floor(valid_gt / base)
      a2 <- valid_gt %% base
      all_alleles <- sort(unique(c(a1, a2)))
      all_alleles <- all_alleles[all_alleles >= 0L]
      if (length(all_alleles) == 0L) return(list(rd=0.0, n=efpop))

      genefreq <- sapply(all_alleles, function(a)
        (sum(a1==a) + sum(a2==a)) / (2L * N))

      r    <- sqrt(1.0 / (N + 1.0))
      H_ii <- sapply(all_alleles, function(a) sum(a1==a & a2==a))
      H_iX <- sapply(all_alleles, function(a) sum((a1==a & a2!=a)|(a2==a & a1!=a)))
      hotot <- sum(H_ii)

      # cpt=0 init — absent branch (no null homos)
      p <- numeric(length(all_alleles))
      for (ai in seq_along(all_alleles)) {
        if (genefreq[ai] <= 0) { p[ai] <- 0.0; next }
        ii <- H_ii[ai]; jj <- H_iX[ai]
        X <- 1.0 + hotot - ii + (N - hotot) - jj; Y <- N + 1.0
        p[ai] <- 1.0 - sqrt(max(0.0, X / Y))
      }

      # EM iterations — denominator = N (n_valid)
      for (iter in seq_len(5000L)) {
        new_p <- numeric(length(all_alleles))
        ri <- 0.0; re <- 0L
        for (ai in seq_along(all_alleles)) {
          if (genefreq[ai] <= 0) { new_p[ai] <- 0.0; next }
          pa <- p[ai]; denom <- pa + 2.0 * r
          if (denom <= 0) { new_p[ai] <- 0.0; next }
          p_new     <- (pa + r)/denom * (H_ii[ai]/N) + H_iX[ai]/(2.0*N)
          ri        <- ri + r/denom * (H_ii[ai]/N)
          new_p[ai] <- p_new
          if (abs(p_new - pa) > 1e-6) re <- re + 1L
        }
        r_new <- ri             # no null homo contribution
        if (abs(r_new - r) > 1e-6) re <- re + 1L
        p <- new_p; r <- max(0.0, r_new)
        if (re == 0L) break
      }
      list(rd=r, n=efpop)
    }

    # ══════════════════════════════════════════════════════════════════════
    # Fetch genotypes from DuckDB and dispatch EM per locus treatment
    # ══════════════════════════════════════════════════════════════════════
    fetch_and_run_em <- function(sel_locus="all", sel_pop="all") {
      db_ready()
      con   <- con_r(); hs <- hf_schema_r(); ms <- meta_schema_r()
      base  <- as.integer(base_r())
      hf_q  <- sql_id(con,tbl_hf_r());  meta_q <- sql_id(con,tbl_meta_r())
      hi_q  <- sql_id(con,hs$ind_col);  hl_q   <- sql_id(con,hs$locus_col)
      hg_q  <- sql_id(con,hs$gt_col);   mi_q   <- sql_id(con,ms$ind_col)
      pop_q <- sql_id(con,ms$pop_col)

      filters <- character(0)
      if (!identical(sel_locus,"all"))
        filters <- c(filters,
          sprintf("CAST(h.%s AS VARCHAR)=%s", hl_q, sql_str(con,sel_locus)))
      if (!identical(sel_pop,"all"))
        filters <- c(filters,
          sprintf("CAST(m.%s AS VARCHAR)=%s", pop_q, sql_str(con,sel_pop)))
      w_extra <- if (length(filters))
        paste0(" AND ", paste(filters, collapse=" AND ")) else ""

      # LEFT JOIN locus_order so ORDER BY lo._lo_rank works
      sql <- sprintf("
        WITH %s
        SELECT
          CAST(m.%s AS VARCHAR) AS Population,
          CAST(h.%s AS VARCHAR) AS Marker,
          h.%s                  AS gt
        FROM %s h
        INNER JOIN %s m
          ON CAST(h.%s AS VARCHAR) = CAST(m.%s AS VARCHAR)
        LEFT JOIN locus_order lo
          ON CAST(h.%s AS VARCHAR) = lo._lo_marker
        WHERE m.%s IS NOT NULL%s
        ORDER BY lo._lo_rank ASC, Population",
        locus_order_cte(con,hf_q,hl_q),
        pop_q, hl_q, hg_q,
        hf_q,
        meta_q, hi_q, mi_q,
        hl_q,
        pop_q, w_extra)

      raw <- DBI::dbGetQuery(con, sql)
      if (nrow(raw) == 0L) return(data.frame())

      treatments   <- locus_treatments_r()
      locus_levels <- markers_r()
      combos       <- unique(raw[,c("Population","Marker"),drop=FALSE])

      results <- vector("list", nrow(combos))
      for (i in seq_len(nrow(combos))) {
        pop_i  <- combos$Population[i]
        mark_i <- combos$Marker[i]
        gts    <- raw$gt[raw$Population==pop_i & raw$Marker==mark_i]

        treat <- treatments[mark_i]
        if (is.na(treat) || length(treat)==0L) treat <- "null_homo"

        em <- if (identical(as.character(treat),"absent"))
          em_absent(gts, base)
        else
          em_null_homo(gts, base)

        n_exp <- em$n * (em$rd^2)
        results[[i]] <- data.frame(
          Locus        = mark_i,
          Population   = pop_i,
          p_nulls      = round(em$rd,        5),
          N            = as.integer(em$n),
          N_exp_blanks = round(n_exp,         9),
          p_nulls_x_N  = round(em$rd * em$n, 5),
          stringsAsFactors = FALSE
        )
      }

      out <- do.call(rbind, results)
      if (!is.null(locus_levels) && length(locus_levels)) {
        out$Locus <- factor(out$Locus, levels=locus_levels)
        out <- out[order(out$Locus, out$Population),]
        out$Locus <- as.character(out$Locus)
      }
      out
    }

    # ── Ready guards ───────────────────────────────────────────────────────
    t1_ready_r <- reactive({ req(input$run_t1 > 0L); db_ready(); TRUE })
    t2_ready_r <- reactive({ req(input$run_t2 > 0L); db_ready(); TRUE })

    # ── Tab 1 data ─────────────────────────────────────────────────────────
    t1_data_r <- reactive({
      t1_ready_r()
      withProgress(message="Running EM algorithm (FreeNA)...", value=0.2, {
        d <- fetch_and_run_em(
          sel_locus = safe_choice(input$t1_locus,"all"),
          sel_pop   = safe_choice(input$t1_pop,  "all"))
        setProgress(1); d
      })
    })

    # ── Tab 2 global summary ───────────────────────────────────────────────
    # Verified formulas (match reference output exactly):
    #
    #   Av(N_exp_blanks) = SUM_pops(N_i * p_i^2)          [total expected null homos]
    #   Av(p_nulls)      = SUM(N_i * p_i) / SUM(N_i)      [N-weighted mean]
    #   N_tot            = total individuals (all pops)
    #   N_blanks         = observed missing genotypes (gt=0, all pops)
    #   f(expBlanks)     = Av(N_exp_blanks) / N_tot
    #   p_nulls          = Av(p_nulls)  [same value, N-weighted mean]
    #   p-value          = REMOVED (not applicable)
    t2_data_r <- reactive({
      t2_ready_r()
      withProgress(message="Computing global summary...", value=0.2, {
        sel_loc <- safe_choice(input$t2_locus,"all")

        # Run EM for all pops (needed for aggregation)
        long <- fetch_and_run_em(sel_locus=sel_loc, sel_pop="all")
        if (nrow(long)==0L) return(data.frame())

        # Observed missing counts per locus (all pops combined)
        db_ready()
        con   <- con_r(); hs <- hf_schema_r(); ms <- meta_schema_r()
        hf_q  <- sql_id(con,tbl_hf_r());  meta_q <- sql_id(con,tbl_meta_r())
        hi_q  <- sql_id(con,hs$ind_col);  hl_q   <- sql_id(con,hs$locus_col)
        hg_q  <- sql_id(con,hs$gt_col);   mi_q   <- sql_id(con,ms$ind_col)
        pop_q <- sql_id(con,ms$pop_col)

        lf_extra <- if (!identical(sel_loc,"all"))
          sprintf(" AND CAST(h.%s AS VARCHAR)=%s", hl_q, sql_str(con,sel_loc)) else ""

        # LEFT JOIN locus_order so ORDER BY _lo_rank works
        obs <- DBI::dbGetQuery(con, sprintf("
          WITH %s
          SELECT
            CAST(h.%s AS VARCHAR) AS Marker,
            COUNT(*) AS N_tot,
            SUM(CASE WHEN h.%s IS NULL OR h.%s <= 0 THEN 1 ELSE 0 END) AS N_blanks,
            MIN(lo._lo_rank) AS _lo_rank
          FROM %s h
          INNER JOIN %s m
            ON CAST(h.%s AS VARCHAR) = CAST(m.%s AS VARCHAR)
          LEFT JOIN locus_order lo
            ON CAST(h.%s AS VARCHAR) = lo._lo_marker
          WHERE m.%s IS NOT NULL%s
          GROUP BY CAST(h.%s AS VARCHAR)
          ORDER BY _lo_rank ASC",
          locus_order_cte(con,hf_q,hl_q),
          hl_q,
          hg_q, hg_q,
          hf_q,
          meta_q, hi_q, mi_q,
          hl_q,
          pop_q, lf_extra,
          hl_q))

        locus_levels <- markers_r()
        loci_in_long <- if (!is.null(locus_levels) && length(locus_levels))
          locus_levels[locus_levels %in% unique(long$Locus)]
        else unique(long$Locus)

        rows <- lapply(loci_in_long, function(loc) {
          sub     <- long[long$Locus==loc,,drop=FALSE]
          if (nrow(sub)==0L) return(NULL)

          obs_row  <- obs[obs$Marker==loc,,drop=FALSE]
          n_tot    <- if (nrow(obs_row)) as.integer(obs_row$N_tot[1])    else sum(sub$N)
          n_blanks <- if (nrow(obs_row)) as.integer(obs_row$N_blanks[1]) else NA_integer_

          # Av(N_exp_blanks) = SUM(N_i * p_i^2)
          av_n_exp <- sum(sub$N * (sub$p_nulls^2), na.rm=TRUE)

          # Av(p_nulls) = N-weighted mean = SUM(N_i * p_i) / SUM(N_i)
          vidx  <- !is.na(sub$p_nulls)
          av_p  <- if (any(vidx) && sum(sub$N[vidx])>0)
            sum(sub$p_nulls[vidx] * sub$N[vidx]) / sum(sub$N[vidx])
          else NA_real_

          # f(expBlanks) = Av(N_exp_blanks) / N_tot
          f_exp <- if (!is.na(av_n_exp) && n_tot>0) av_n_exp / n_tot else NA_real_

          data.frame(
            Locus        = loc,
            Av_N_exp     = round(av_n_exp, 9),   # SUM(N_i * p_i^2)
            Av_p_nulls   = round(av_p,     9),   # N-weighted mean
            N_tot        = n_tot,
            N_blanks     = n_blanks,
            f_expBlanks  = round(f_exp,    9),   # Av_N_exp / N_tot
            p_nulls      = round(av_p,     9),   # = Av_p_nulls
            stringsAsFactors = FALSE
          )
        })

        setProgress(1)
        do.call(rbind, Filter(Negate(is.null), rows))
      })
    })

    # ── Value boxes ────────────────────────────────────────────────────────
    output$vb_loci <- renderUI({
      tryCatch(tags$span(length(markers_r())),
               error=function(e) tags$span("\u2014"))
    })

    output$vb_pops <- renderUI({
      tryCatch(tags$span(length(pops_r())),
               error=function(e) tags$span("\u2014"))
    })

    output$vb_n <- renderUI({
      tryCatch({
        db_ready(); con <- con_r(); ms <- meta_schema_r()
        n <- DBI::dbGetQuery(con, sprintf(
          "SELECT COUNT(DISTINCT CAST(%s AS VARCHAR)) AS n FROM %s WHERE %s IS NOT NULL",
          sql_id(con,ms$ind_col), sql_id(con,tbl_meta_r()),
          sql_id(con,ms$ind_col)))$n[[1]]
        tags$span(n)
      }, error=function(e) tags$span("\u2014"))
    })

    output$vb_avg_null <- renderUI({
      tryCatch({
        d <- t1_data_r()
        if (nrow(d)==0||all(is.na(d$p_nulls))) return(tags$span("\u2014"))
        v   <- round(mean(d$p_nulls, na.rm=TRUE), 4)
        col <- if(v>.20)"#9d174d" else if(v>.10)"#854d0e" else "#166534"
        tags$span(style=paste0("color:",col,";"), v)
      }, error=function(e) tags$span("\u2014"))
    })

    output$vb_max_null <- renderUI({
      tryCatch({
        d <- t1_data_r()
        if (nrow(d)==0||all(is.na(d$p_nulls))) return(tags$span("\u2014"))
        v   <- round(max(d$p_nulls, na.rm=TRUE), 4)
        col <- if(v>.30)"#9d174d" else if(v>.15)"#854d0e" else "#166534"
        tags$span(style=paste0("color:",col,";"), v)
      }, error=function(e) tags$span("\u2014"))
    })

    # ── Tab 1 DT ───────────────────────────────────────────────────────────
    output$dt_t1 <- DT::renderDT({
      d <- t1_data_r()
      shiny::validate(shiny::need(nrow(d)>0,
        "No data. Select parameters and click Compute."))

      disp        <- d
      names(disp) <- c("Locus names","Farm","p_nulls","N",
                       "N_exp_blanks","p_nulls\u00d7N")

      DT::datatable(disp,
        rownames=FALSE,
        options=list(
          pageLength=20, scrollX=TRUE, dom="lftip",
          columnDefs=list(list(className="dt-right", targets=2:5))),
        class="compact hover stripe"
      ) |>
        DT::formatRound("p_nulls",        5) |>
        DT::formatRound("N_exp_blanks",   9) |>
        DT::formatRound("p_nulls\u00d7N", 5) |>
        DT::formatStyle("p_nulls",
          backgroundColor=DT::styleInterval(
            c(0.05,0.10,0.20,0.30),
            c("#f0fdf4","#dcfce7","#fefce8","#fff7ed","#fef2f2"))) |>
        DT::formatStyle("Locus names", fontWeight="600", color="#0f172a") |>
        DT::formatStyle("Farm",        color="#475569")
    }, server=TRUE)

    # ── Tab 2 DT ───────────────────────────────────────────────────────────
    output$dt_t2 <- DT::renderDT({
      d <- t2_data_r()
      shiny::validate(shiny::need(nrow(d)>0,
        "No data. Select parameters and click Compute."))

      disp        <- d
      names(disp) <- c("Locus names",
                       "Av(N_exp_blanks)",  # SUM(N_i * p_i^2)
                       "Av(p_nulls)",       # N-weighted mean
                       "N_tot","N_blanks",
                       "f(expBlanks)",      # Av(N_exp_blanks) / N_tot
                       "p_nulls")           # = Av(p_nulls)

      DT::datatable(disp,
        rownames=FALSE,
        options=list(
          pageLength=20, scrollX=TRUE, dom="lftip",
          columnDefs=list(list(className="dt-right", targets=1:6))),
        class="compact hover stripe"
      ) |>
        DT::formatRound("Av(N_exp_blanks)", 9) |>
        DT::formatRound("Av(p_nulls)",      9) |>
        DT::formatRound("f(expBlanks)",     9) |>
        DT::formatRound("p_nulls",          9) |>
        DT::formatStyle("p_nulls",
          backgroundColor=DT::styleInterval(
            c(0.05,0.10,0.20),
            c("#f0fdf4","#dcfce7","#fefce8","#fef2f2"))) |>
        DT::formatStyle("Locus names", fontWeight="600", color="#0f172a")
    }, server=TRUE)

    # ── Downloads Tab 1 ────────────────────────────────────────────────────
    output$dl_t1_csv <- downloadHandler(
      filename=function() paste0("null_allele_per_pop_locus_",Sys.Date(),".csv"),
      content=function(file) {
        d <- t1_data_r(); if(nrow(d)==0) return(invisible(NULL))
        names(d) <- c("Locus_names","Farm","p_nulls","N","N_exp_blanks","p_nulls_x_N")
        write.csv(d, file, row.names=FALSE)
      }
    )
    output$dl_t1_txt <- downloadHandler(
      filename=function() paste0("null_allele_per_pop_locus_",Sys.Date(),".txt"),
      content=function(file) {
        d <- t1_data_r(); if(nrow(d)==0) return(invisible(NULL))
        names(d) <- c("Locus_names","Farm","p_nulls","N","N_exp_blanks","p_nulls_x_N")
        write.table(d, file, sep="\t", row.names=FALSE, quote=FALSE)
      }
    )

    # ── Downloads Tab 2 ────────────────────────────────────────────────────
    output$dl_t2_csv <- downloadHandler(
      filename=function() paste0("null_allele_global_",Sys.Date(),".csv"),
      content=function(file) {
        d <- t2_data_r(); if(nrow(d)==0) return(invisible(NULL))
        names(d) <- c("Locus_names","Av_N_exp_blanks","Av_p_nulls",
                      "N_tot","N_blanks","f_expBlanks","p_nulls")
        write.csv(d, file, row.names=FALSE)
      }
    )
    output$dl_t2_txt <- downloadHandler(
      filename=function() paste0("null_allele_global_",Sys.Date(),".txt"),
      content=function(file) {
        d <- t2_data_r(); if(nrow(d)==0) return(invisible(NULL))
        names(d) <- c("Locus_names","Av_N_exp_blanks","Av_p_nulls",
                      "N_tot","N_blanks","f_expBlanks","p_nulls")
        write.table(d, file, sep="\t", row.names=FALSE, quote=FALSE)
      }
    )

  })
}