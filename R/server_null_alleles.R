# module/server_null_alleles.R
# Null allele frequency estimation (EM), FST-ENA correction, DCSE-INA genetic distance
# Bootstrap 95% CI: over loci (global + per-locus) and over individuals (sub-samples)
#
# CORRECTIONS in this version:
#   1. cs_distance(): renamed pi/pj -> fi/fj to avoid collision with base::pi.
#   2. compute_dc_pairwise() — DCSE-INA: now follows Pascal prod_CS_correction EXACTLY:
#      freq_ina = c(pfreq, null=rd) even when pfreq is empty (pop with only absent/null geno).
#      Check ni_ci>0 && ni_cj>0 only (not length(pfreq)>0), mirroring Pascal's ajustement_r.
#   3. precompute_locus_stats_dc() — same fix for INA branch.
#   4. boot_tbl(): FIX locus name display — "Locus","Pop1","Pop2" all treated as character
#      labels, not coerced to numeric. Uses seq_len(nrow(d)) loop instead of apply()
#      to preserve column types correctly.
#   5. Per-locus tables in boot_loci_fst_pair() and boot_loci_dc_pair(): locus names
#      preserved as character, not factor-coerced.
#
# References:
#   Dempster, Laird & Rubin (1977)  — EM algorithm
#   Chapuis & Estoup (2007)         — FreeNA: ENA and INA corrections
#   Weir (1996)                     — FST following Genepop method
#   Cavalli-Sforza & Edwards (1967) — Chord genetic distance (DCSE)

server_null_alleles <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    # ── Helpers ────────────────────────────────────────────────────────────────
    `%||%` <- function(a, b) if (!is.null(a)) a else b
    safe_choice <- function(x, default = "all") {
      if (is.null(x) || length(x) == 0L || identical(x, "") || all(is.na(x))) default
      else as.character(x[[1]])
    }
    sql_id   <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))
    sql_str  <- function(con, x) as.character(DBI::dbQuoteString(con, x))
    treat_id <- function(loc) paste0("treat_", gsub("[^A-Za-z0-9]", "_", loc))

    # ── DB plumbing ────────────────────────────────────────────────────────────
    db_tick    <- reactive({ rv$db_tick })
    con_r      <- reactive({ req(rv$con); rv$con })
    tbl_meta_r <- reactive({ rv$tbl_meta %||% "meta" })

    tbl_hf_r <- reactive({
      con <- con_r()
      if (exists("duck_tbl_exists", mode = "function", inherits = TRUE) &&
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
                 else shiny::validate(shiny::need(FALSE, "No individual column found in meta."))
      pop_col <- c("Population","population","pop","pop_code")[
        c("Population","population","pop","pop_code") %in% cols][1]
      shiny::validate(shiny::need(!is.na(pop_col), "No population column found in meta."))
      list(ind_col = ind_col, pop_col = pop_col)
    })

    locus_order_cte <- function(con, hf_tbl_q, hl_q)
      sprintf("locus_order AS (
  SELECT CAST(%s AS VARCHAR) AS _lo_marker, MIN(rowid) AS _lo_rank
  FROM %s GROUP BY CAST(%s AS VARCHAR))", hl_q, hf_tbl_q, hl_q)

    # ── Marker / population lists ──────────────────────────────────────────────
    pops_r <- reactive({
      db_ready(); con <- con_r(); ms <- meta_schema_r()
      as.character(DBI::dbGetQuery(con, sprintf(
        "SELECT DISTINCT CAST(%s AS VARCHAR) AS p FROM %s WHERE %s IS NOT NULL ORDER BY p",
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
      updateSelectInput(session,"t1_locus",choices=c("All loci"="all",stats::setNames(markers,markers)),selected="all")
      updateSelectInput(session,"t1_pop",  choices=c("All populations"="all",stats::setNames(pops,pops)),selected="all")
      updateSelectInput(session,"t2_locus",choices=c("All loci"="all",stats::setNames(markers,markers)),selected="all")
      updateSelectInput(session,"fl_locus",choices=c("All loci"="all",stats::setNames(markers,markers)),selected="all")
      updateSelectInput(session,"fl_pop1", choices=c("All pairs"="all",stats::setNames(pops,pops)),selected="all")
      updateSelectInput(session,"fl_pop2", choices=c("All pairs"="all",stats::setNames(pops,pops)),selected="all")
    })

    # ── Per-locus treatment selector ───────────────────────────────────────────
    output$locus_treatment_ui <- renderUI({
      ns_fn <- session$ns; markers <- markers_r()
      if (length(markers) == 0L) return(tags$p("No markers loaded yet."))
      items <- lapply(markers, function(loc)
        tags$div(class="na-treat-item",
          tags$div(class="na-treat-lbl",loc),
          selectInput(inputId=ns_fn(treat_id(loc)),label=NULL,
            choices=c("999999 \u2014 null homozygote"="null_homo",
                      "000000 \u2014 absent / PCR failure"="absent"),
            selected="null_homo",width="100%")))
      tags$div(class="na-treat-grid",items)
    })

    locus_treatments_r <- reactive({
      markers <- markers_r()
      treats  <- sapply(markers, function(loc) {
        val <- input[[treat_id(loc)]]
        if (is.null(val) || !val %in% c("null_homo","absent")) "null_homo" else val
      })
      stats::setNames(treats, markers)
    })

    # ══════════════════════════════════════════════════════════════════════════
    #  EM ALGORITHM — exact translation of rDempster_per_locus (Pascal FreeNA)
    # ══════════════════════════════════════════════════════════════════════════
    em_freena <- function(gt_vec, base, treat = "null_homo") {
      efpop      <- length(gt_vec)
      absent_msk <- is.na(gt_vec) | gt_vec <= 0L
      n_absent   <- sum(absent_msk)
      valid_gt   <- gt_vec[!absent_msk]

      empty <- list(rd=0.0, pfreq=numeric(0), genefreq_obs=numeric(0),
                    H_ii=numeric(0), H_iX=numeric(0), N=0L, efpop=efpop,
                    n_absent=n_absent, n_null_homo=0L, alleles=integer(0), n_valid_geno=0L)

      if (length(valid_gt) == 0L) return(empty)

      a1_all <- floor(valid_gt / base)
      a2_all <- valid_gt %% base
      null_code     <- if (base >= 1000L) 999L else 99L
      null_homo_msk <- (a1_all == null_code) & (a2_all == null_code)
      n_null_homo   <- sum(null_homo_msk)

      valid_a1 <- a1_all[!null_homo_msk]
      valid_a2 <- a2_all[!null_homo_msk]
      alleles  <- sort(unique(c(valid_a1, valid_a2)))
      alleles  <- alleles[alleles >= 0L & alleles != null_code]

      N <- efpop - n_absent
      if (N == 0L || length(alleles) == 0L) {
        empty$N <- N; empty$n_null_homo <- n_null_homo; return(empty)
      }

      n_valid_geno <- N - n_null_homo
      genefreq_obs <- sapply(alleles, function(a)
        (sum(valid_a1==a) + sum(valid_a2==a)) / (2L * n_valid_geno))
      H_ii  <- sapply(alleles, function(a) sum(valid_a1==a & valid_a2==a))
      H_iX  <- sapply(alleles, function(a)
        sum((valid_a1==a & valid_a2!=a) | (valid_a2==a & valid_a1!=a)))
      hotot <- sum(H_ii)

      rd <- if (treat == "null_homo" && n_null_homo > 0L)
              sqrt(n_null_homo / N) else sqrt(1.0 / (N + 1.0))

      p <- numeric(length(alleles))
      for (ai in seq_along(alleles)) {
        if (genefreq_obs[ai] <= 0) { p[ai] <- 0.0; next }
        ii <- H_ii[ai]; jj <- H_iX[ai]
        if (treat == "null_homo" && n_null_homo > 0L) {
          X <- n_null_homo + hotot - ii + (N - n_null_homo - hotot) - jj; Y <- N
        } else {
          X <- 1.0 + hotot - ii + (N - hotot) - jj; Y <- N + 1.0
        }
        p[ai] <- 1.0 - sqrt(max(0.0, X / Y))
      }

      for (iter in seq_len(5000L)) {
        new_p <- numeric(length(alleles)); rdi <- 0.0; re <- 0L
        for (ai in seq_along(alleles)) {
          if (genefreq_obs[ai] <= 0) { new_p[ai] <- 0.0; next }
          pa <- p[ai]; denom <- pa + 2.0 * rd
          if (denom <= 0) { new_p[ai] <- 0.0; next }
          p_new     <- (pa + rd) / denom * (H_ii[ai] / N) + H_iX[ai] / (2.0 * N)
          rdi       <- rdi + rd / denom * (H_ii[ai] / N)
          new_p[ai] <- p_new
          if (abs(p_new - pa) > 1e-6) re <- re + 1L
        }
        rd_new <- if (treat == "null_homo") rdi + (2.0 * n_null_homo) / (2.0 * N) else rdi
        if (abs(rd_new - rd) > 1e-6) re <- re + 1L
        p <- new_p; rd <- max(0.0, rd_new)
        if (re == 0L) break
      }

      a_chr <- as.character(alleles)
      list(rd=rd, pfreq=stats::setNames(p,a_chr),
           genefreq_obs=stats::setNames(genefreq_obs,a_chr),
           H_ii=stats::setNames(H_ii,a_chr), H_iX=stats::setNames(H_iX,a_chr),
           N=N, efpop=efpop, n_absent=n_absent, n_null_homo=n_null_homo,
           alleles=alleles, n_valid_geno=n_valid_geno)
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  WEIR (1996) FST COMPONENTS
    # ══════════════════════════════════════════════════════════════════════════
    weir_components_allele <- function(pop_list, use_corr = FALSE) {
      r     <- length(pop_list)
      N_tot <- sum(sapply(pop_list, `[[`, "ni"))
      N2    <- sum(sapply(pop_list, function(p) p$ni^2))
      if (N_tot == 0L || r < 2L) return(list(s2P=0.0, s2I=0.0, s2G=0.0))
      nc <- (N_tot - N2/N_tot) / (r - 1)
      if (nc <= 0 || N_tot - r <= 0) return(list(s2P=0.0, s2I=0.0, s2G=0.0))
      snA  <- sum(sapply(pop_list, `[[`, "nA"))
      s2A  <- sum(sapply(pop_list, function(p) if (p$ni>0) p$nA^2/(2*p$ni) else 0.0))
      sAA  <- if (use_corr) sum(sapply(pop_list,`[[`,"AA_corr"))
              else           sum(sapply(pop_list,`[[`,"AA"))
      MSG  <- (0.5*snA - sAA) / N_tot
      MSI  <- (0.5*snA + sAA - s2A) / (N_tot - r)
      MSP  <- (s2A - 0.5*snA^2/N_tot) / (r - 1)
      list(s2P=(MSP-MSI)/(2*nc), s2I=0.5*(MSI-MSG), s2G=MSG)
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  CAVALLI-SFORZA & EDWARDS (1967) CHORD DISTANCE
    #  FIX: fi/fj instead of pi/pj to avoid collision with base::pi
    #  Pascal: 2/pi * sqrt(2*(1-CSprod))  where pi = 3.14159...
    # ══════════════════════════════════════════════════════════════════════════
    cs_distance <- function(freq_i, freq_j) {
      alleles <- union(names(freq_i), names(freq_j))
      csprod  <- 0.0
      for (a in alleles) {
        fi <- freq_i[a]; fj <- freq_j[a]
        fi <- if (is.null(fi) || is.na(fi)) 0.0 else as.numeric(fi)
        fj <- if (is.null(fj) || is.na(fj)) 0.0 else as.numeric(fj)
        if (fi > 0 && fj > 0) csprod <- csprod + sqrt(fi * fj)
      }
      if (csprod > 1.0) return(NA_real_)   # dcnotapploc in Pascal
      (2.0 / base::pi) * sqrt(2.0 * (1.0 - csprod))
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  BUILD INA FREQUENCY VECTOR — Pascal: prod_CS_correction / ajustement_r
    #  nallc = nall + 1; corrdgenefreq[iloc, ipop, nallc-1] = rd[iloc, ipop]
    #  This must be done even when pfreq is empty (pop with only absent/null genos)
    #  => always append null allele with freq = rd, regardless of pfreq length
    # ══════════════════════════════════════════════════════════════════════════
    make_ina_freq <- function(em) {
      # Pascal ajustement_r: corrdgenefreq[..., nallc-1] := rd[iloc, ipop]
      # Works even if pfreq is empty (only null homos in this pop)
      c(em$pfreq, `__null__` = em$rd)
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  FETCH ALL GENOTYPES → raw_data_r and em_results_r
    # ══════════════════════════════════════════════════════════════════════════
    raw_data_r <- reactive({
      db_ready()
      con   <- con_r(); hs <- hf_schema_r(); ms <- meta_schema_r()
      hf_q  <- sql_id(con,tbl_hf_r());  meta_q <- sql_id(con,tbl_meta_r())
      hi_q  <- sql_id(con,hs$ind_col);  hl_q   <- sql_id(con,hs$locus_col)
      hg_q  <- sql_id(con,hs$gt_col);   mi_q   <- sql_id(con,ms$ind_col)
      pop_q <- sql_id(con,ms$pop_col)
      sql <- sprintf("
        WITH %s
        SELECT
          CAST(m.%s AS VARCHAR) AS Population,
          CAST(m.%s AS VARCHAR) AS Individual,
          CAST(h.%s AS VARCHAR) AS Marker,
          h.%s                  AS gt
        FROM %s h
        INNER JOIN %s m ON CAST(h.%s AS VARCHAR) = CAST(m.%s AS VARCHAR)
        LEFT JOIN locus_order lo ON CAST(h.%s AS VARCHAR) = lo._lo_marker
        WHERE m.%s IS NOT NULL
        ORDER BY lo._lo_rank ASC, Population, Individual",
        locus_order_cte(con,hf_q,hl_q),
        pop_q, sql_id(con,ms$ind_col), hl_q, hg_q,
        hf_q, meta_q, hi_q, mi_q, hl_q, pop_q)
      DBI::dbGetQuery(con, sql)
    })

    em_results_r <- reactive({
      db_ready()
      raw        <- raw_data_r()
      if (nrow(raw) == 0L) return(list())
      base       <- as.integer(base_r())
      markers    <- markers_r(); pops <- pops_r()
      treatments <- locus_treatments_r()
      em_res <- list()
      for (loc in markers) {
        em_res[[loc]] <- list()
        treat <- as.character(treatments[loc] %||% "null_homo")
        for (pop in pops) {
          gts <- raw$gt[raw$Marker == loc & raw$Population == pop]
          em_res[[loc]][[pop]] <-
            if (length(gts) == 0L)
              list(rd=0.0, pfreq=numeric(0), genefreq_obs=numeric(0),
                   H_ii=numeric(0), H_iX=numeric(0), N=0L, efpop=0L,
                   n_absent=0L, n_null_homo=0L, alleles=integer(0), n_valid_geno=0L)
            else em_freena(gts, base, treat)
        }
      }
      em_res
    })

    # ══════════════════════════════════════════════════════════════════════════
    #  TABS 1 & 2 — p_nulls per locus x population
    # ══════════════════════════════════════════════════════════════════════════
    fetch_and_run_em_simple <- function(sel_locus = "all", sel_pop = "all") {
      db_ready()
      con   <- con_r(); hs <- hf_schema_r(); ms <- meta_schema_r()
      base  <- as.integer(base_r())
      hf_q  <- sql_id(con,tbl_hf_r());  meta_q <- sql_id(con,tbl_meta_r())
      hi_q  <- sql_id(con,hs$ind_col);  hl_q   <- sql_id(con,hs$locus_col)
      hg_q  <- sql_id(con,hs$gt_col);   mi_q   <- sql_id(con,ms$ind_col)
      pop_q <- sql_id(con,ms$pop_col)
      filters <- character(0)
      if (!identical(sel_locus,"all"))
        filters <- c(filters, sprintf("CAST(h.%s AS VARCHAR)=%s",hl_q,sql_str(con,sel_locus)))
      if (!identical(sel_pop,"all"))
        filters <- c(filters, sprintf("CAST(m.%s AS VARCHAR)=%s",pop_q,sql_str(con,sel_pop)))
      w_extra <- if (length(filters)) paste0(" AND ",paste(filters,collapse=" AND ")) else ""
      sql <- sprintf("
        WITH %s
        SELECT CAST(m.%s AS VARCHAR) AS Population, CAST(h.%s AS VARCHAR) AS Marker, h.%s AS gt
        FROM %s h INNER JOIN %s m ON CAST(h.%s AS VARCHAR) = CAST(m.%s AS VARCHAR)
        LEFT JOIN locus_order lo ON CAST(h.%s AS VARCHAR) = lo._lo_marker
        WHERE m.%s IS NOT NULL%s ORDER BY lo._lo_rank ASC, Population",
        locus_order_cte(con,hf_q,hl_q), pop_q, hl_q, hg_q,
        hf_q, meta_q, hi_q, mi_q, hl_q, pop_q, w_extra)
      raw <- DBI::dbGetQuery(con, sql)
      if (nrow(raw) == 0L) return(data.frame())
      treatments   <- locus_treatments_r()
      locus_levels <- markers_r()
      combos       <- unique(raw[,c("Population","Marker"),drop=FALSE])
      results      <- vector("list", nrow(combos))
      for (i in seq_len(nrow(combos))) {
        pop_i <- combos$Population[i]; mark_i <- combos$Marker[i]
        gts   <- raw$gt[raw$Population==pop_i & raw$Marker==mark_i]
        treat <- as.character(treatments[mark_i] %||% "null_homo")
        em    <- em_freena(gts, base, treat)
        n_exp <- em$N * (em$rd^2)
        results[[i]] <- data.frame(
          Locus=mark_i, Population=pop_i,
          p_nulls=round(em$rd,5), N=as.integer(em$N),
          N_exp_blanks=round(n_exp,9), p_nulls_x_N=round(em$rd*em$N,5),
          stringsAsFactors=FALSE)
      }
      out <- do.call(rbind, results)
      if (!is.null(locus_levels) && length(locus_levels)) {
        out$Locus <- factor(out$Locus, levels=locus_levels)
        out <- out[order(out$Locus, out$Population),]
        out$Locus <- as.character(out$Locus)
      }
      out
    }

    t1_ready_r <- reactive({ req(input$run_t1 > 0L); db_ready(); TRUE })
    t2_ready_r <- reactive({ req(input$run_t2 > 0L); db_ready(); TRUE })

    t1_data_r <- reactive({
      t1_ready_r()
      withProgress(message="Running EM algorithm (FreeNA)...", value=0.2, {
        d <- fetch_and_run_em_simple(
          sel_locus=safe_choice(input$t1_locus,"all"),
          sel_pop  =safe_choice(input$t1_pop,  "all"))
        setProgress(1); d
      })
    })

    t2_data_r <- reactive({
      t2_ready_r()
      withProgress(message="Computing global summary...", value=0.2, {
        sel_loc <- safe_choice(input$t2_locus,"all")
        long    <- fetch_and_run_em_simple(sel_locus=sel_loc, sel_pop="all")
        if (nrow(long)==0L) return(data.frame())
        db_ready()
        con   <- con_r(); hs <- hf_schema_r(); ms <- meta_schema_r()
        hf_q  <- sql_id(con,tbl_hf_r()); meta_q <- sql_id(con,tbl_meta_r())
        hi_q  <- sql_id(con,hs$ind_col); hl_q   <- sql_id(con,hs$locus_col)
        hg_q  <- sql_id(con,hs$gt_col);  mi_q   <- sql_id(con,ms$ind_col)
        pop_q <- sql_id(con,ms$pop_col)
        lf_extra <- if (!identical(sel_loc,"all"))
          sprintf(" AND CAST(h.%s AS VARCHAR)=%s",hl_q,sql_str(con,sel_loc)) else ""
        obs <- DBI::dbGetQuery(con, sprintf("
          WITH %s
          SELECT CAST(h.%s AS VARCHAR) AS Marker, COUNT(*) AS N_tot,
            SUM(CASE WHEN h.%s IS NULL OR h.%s <= 0 THEN 1 ELSE 0 END) AS N_blanks,
            MIN(lo._lo_rank) AS _lo_rank
          FROM %s h INNER JOIN %s m ON CAST(h.%s AS VARCHAR)=CAST(m.%s AS VARCHAR)
          LEFT JOIN locus_order lo ON CAST(h.%s AS VARCHAR)=lo._lo_marker
          WHERE m.%s IS NOT NULL%s GROUP BY CAST(h.%s AS VARCHAR) ORDER BY _lo_rank ASC",
          locus_order_cte(con,hf_q,hl_q), hl_q, hg_q, hg_q,
          hf_q, meta_q, hi_q, mi_q, hl_q, pop_q, lf_extra, hl_q))
        locus_levels <- markers_r()
        loci_in_long <- if (!is.null(locus_levels)&&length(locus_levels))
          locus_levels[locus_levels %in% unique(long$Locus)] else unique(long$Locus)
        rows <- lapply(loci_in_long, function(loc) {
          sub      <- long[long$Locus==loc,,drop=FALSE]
          if (nrow(sub)==0L) return(NULL)
          obs_row  <- obs[obs$Marker==loc,,drop=FALSE]
          n_tot    <- if (nrow(obs_row)) as.integer(obs_row$N_tot[1])    else sum(sub$N)
          n_blanks <- if (nrow(obs_row)) as.integer(obs_row$N_blanks[1]) else NA_integer_
          av_n_exp <- sum(sub$N*(sub$p_nulls^2), na.rm=TRUE)
          vidx     <- !is.na(sub$p_nulls)
          av_p     <- if (any(vidx)&&sum(sub$N[vidx])>0)
            sum(sub$p_nulls[vidx]*sub$N[vidx])/sum(sub$N[vidx]) else NA_real_
          f_exp    <- if (!is.na(av_n_exp)&&n_tot>0) av_n_exp/n_tot else NA_real_
          data.frame(Locus=loc, Av_N_exp=round(av_n_exp,9), Av_p_nulls=round(av_p,9),
                     N_tot=n_tot, N_blanks=n_blanks, f_expBlanks=round(f_exp,9),
                     p_nulls=round(av_p,9), stringsAsFactors=FALSE)
        })
        setProgress(1); do.call(rbind, Filter(Negate(is.null), rows))
      })
    })

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 3 — GLOBAL FST (raw + ENA, multilocus)
    # ══════════════════════════════════════════════════════════════════════════
    compute_fst_global <- function(em_res) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]])
      s1 <- s3 <- s1c <- s3c <- 0.0
      rows <- vector("list", length(markers))
      for (li in seq_along(markers)) {
        loc    <- markers[li]; em_loc <- em_res[[loc]]
        alleles_obs <- sort(unique(unlist(lapply(em_loc, function(e) e$alleles))))
        ni_raw  <- sapply(pops, function(p) { e<-em_loc[[p]]; max(0L,e$efpop-e$n_absent-e$n_null_homo) })
        ni_corr <- sapply(pops, function(p) { e<-em_loc[[p]]; max(0L,e$efpop-e$n_absent) })
        r_raw <- sum(ni_raw>0L); r_corr <- sum(ni_corr>0L)
        N_raw <- sum(ni_raw); N2_raw <- sum(ni_raw^2)
        N_cor <- sum(ni_corr); N2_cor <- sum(ni_corr^2)
        nc_raw  <- if (N_raw>0&&r_raw>1)  (N_raw -N2_raw /N_raw ) /(r_raw -1) else 0.0
        nc_corr <- if (N_cor>0&&r_corr>1) (N_cor -N2_cor /N_cor ) /(r_corr-1) else 0.0
        s1l <- s3l <- s1lc <- s3lc <- 0.0
        for (a in alleles_obs) {
          a_chr <- as.character(a)
          pop_raw <- lapply(pops, function(p) {
            e <- em_loc[[p]]; ni <- max(0L,e$efpop-e$n_absent-e$n_null_homo)
            pf <- if (a_chr %in% names(e$genefreq_obs)) e$genefreq_obs[[a_chr]] else 0.0
            AA <- if (a_chr %in% names(e$H_ii)) e$H_ii[[a_chr]] else 0L
            list(ni=ni, nA=pf*2L*ni, AA=AA, AA_corr=AA)
          })
          cmp <- weir_components_allele(pop_raw, use_corr=FALSE)
          s1l <- s1l+cmp$s2P; s3l <- s3l+cmp$s2P+cmp$s2I+cmp$s2G
          pop_ena <- lapply(pops, function(p) {
            e <- em_loc[[p]]; ni <- max(0L,e$efpop-e$n_absent)
            pf <- if (a_chr %in% names(e$pfreq)) e$pfreq[[a_chr]] else 0.0
            AA <- if (a_chr %in% names(e$H_ii)) e$H_ii[[a_chr]] else 0L
            d  <- pf+2.0*e$rd; AAc <- if (AA>0&&d>0) AA*(pf/d) else 0.0
            list(ni=ni, nA=pf*2L*ni, AA=AA, AA_corr=AAc)
          })
          cmpc <- weir_components_allele(pop_ena, use_corr=TRUE)
          s1lc <- s1lc+cmpc$s2P; s3lc <- s3lc+cmpc$s2P+cmpc$s2I+cmpc$s2G
        }
        fst_loc  <- if (s3l  != 0) s1l /s3l  else NA_real_
        fst_locc <- if (s3lc != 0) s1lc/s3lc else NA_real_
        if (!is.na(fst_loc)  && nc_raw >0) { s1  <- s1 +s1l *nc_raw;  s3  <- s3 +s3l *nc_raw  }
        if (!is.na(fst_locc) && nc_corr>0) { s1c <- s1c+s1lc*nc_corr; s3c <- s3c+s3lc*nc_corr }
        rows[[li]] <- data.frame(Locus=loc,
          FST_raw=round(fst_loc,6), FST_ENA=round(fst_locc,6),
          Delta_FST=round(fst_locc-fst_loc,6),
          N_pops_eff_raw=r_raw, N_pops_eff_ENA=r_corr, stringsAsFactors=FALSE)
      }
      list(global_raw=if(s3>0)s1/s3 else NA_real_,
           global_ena=if(s3c>0)s1c/s3c else NA_real_,
           per_locus=do.call(rbind,rows))
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 4 — PAIRWISE FST (raw + ENA)
    # ══════════════════════════════════════════════════════════════════════════
    compute_fst_pairwise <- function(em_res) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]]); n_pops <- length(pops)
      if (n_pops < 2L) return(list(matrix_raw=NULL,matrix_ena=NULL,long=data.frame()))
      s12p  <- matrix(0.0,n_pops,n_pops,dimnames=list(pops,pops))
      s32p  <- matrix(0.0,n_pops,n_pops,dimnames=list(pops,pops))
      s12pc <- matrix(0.0,n_pops,n_pops,dimnames=list(pops,pops))
      s32pc <- matrix(0.0,n_pops,n_pops,dimnames=list(pops,pops))
      for (loc in markers) {
        em_loc <- em_res[[loc]]
        alleles_obs <- sort(unique(unlist(lapply(em_loc, function(e) e$alleles))))
        for (ii in seq_len(n_pops-1L)) {
          for (jj in seq(ii+1L,n_pops)) {
            pi_n <- pops[ii]; pj_n <- pops[jj]
            ei <- em_loc[[pi_n]]; ej <- em_loc[[pj_n]]
            ni_ri <- max(0L,ei$efpop-ei$n_absent-ei$n_null_homo)
            ni_rj <- max(0L,ej$efpop-ej$n_absent-ej$n_null_homo)
            ni_ci <- max(0L,ei$efpop-ei$n_absent)
            ni_cj <- max(0L,ej$efpop-ej$n_absent)
            N_r   <- ni_ri+ni_rj; N2_r <- ni_ri^2+ni_rj^2
            N_c   <- ni_ci+ni_cj; N2_c <- ni_ci^2+ni_cj^2
            nc_r  <- if (N_r>0&&ni_ri>0&&ni_rj>0) (N_r-N2_r/N_r) else 0.0
            nc_c  <- if (N_c>0&&ni_ci>0&&ni_cj>0) (N_c-N2_c/N_c) else 0.0
            for (a in alleles_obs) {
              a_chr <- as.character(a)
              if (nc_r > 0) {
                pd <- list(
                  list(ni=ni_ri, nA=(if(a_chr %in% names(ei$genefreq_obs)) ei$genefreq_obs[[a_chr]] else 0.0)*2L*ni_ri,
                       AA=if(a_chr %in% names(ei$H_ii)) ei$H_ii[[a_chr]] else 0L, AA_corr=0.0),
                  list(ni=ni_rj, nA=(if(a_chr %in% names(ej$genefreq_obs)) ej$genefreq_obs[[a_chr]] else 0.0)*2L*ni_rj,
                       AA=if(a_chr %in% names(ej$H_ii)) ej$H_ii[[a_chr]] else 0L, AA_corr=0.0))
                cmp <- weir_components_allele(pd, use_corr=FALSE)
                s12p[ii,jj] <- s12p[ii,jj]+cmp$s2P*nc_r
                s32p[ii,jj] <- s32p[ii,jj]+(cmp$s2P+cmp$s2I+cmp$s2G)*nc_r
              }
              if (nc_c > 0) {
                pf_i <- if(a_chr %in% names(ei$pfreq)) ei$pfreq[[a_chr]] else 0.0
                pf_j <- if(a_chr %in% names(ej$pfreq)) ej$pfreq[[a_chr]] else 0.0
                AA_i <- if(a_chr %in% names(ei$H_ii)) ei$H_ii[[a_chr]] else 0L
                AA_j <- if(a_chr %in% names(ej$H_ii)) ej$H_ii[[a_chr]] else 0L
                di   <- pf_i+2.0*ei$rd; dj <- pf_j+2.0*ej$rd
                pdc  <- list(
                  list(ni=ni_ci, nA=pf_i*2L*ni_ci, AA=AA_i, AA_corr=if(AA_i>0&&di>0) AA_i*(pf_i/di) else 0.0),
                  list(ni=ni_cj, nA=pf_j*2L*ni_cj, AA=AA_j, AA_corr=if(AA_j>0&&dj>0) AA_j*(pf_j/dj) else 0.0))
                cmpc <- weir_components_allele(pdc, use_corr=TRUE)
                s12pc[ii,jj] <- s12pc[ii,jj]+cmpc$s2P*nc_c
                s32pc[ii,jj] <- s32pc[ii,jj]+(cmpc$s2P+cmpc$s2I+cmpc$s2G)*nc_c
              }
            }
          }
        }
      }
      mat_raw <- matrix(NA_real_,n_pops,n_pops,dimnames=list(pops,pops))
      mat_ena <- matrix(NA_real_,n_pops,n_pops,dimnames=list(pops,pops))
      for (ii in seq_len(n_pops-1L)) for (jj in seq(ii+1L,n_pops)) {
        mat_raw[jj,ii] <- if (s32p[ii,jj] >0) s12p[ii,jj] /s32p[ii,jj]  else NA_real_
        mat_ena[jj,ii] <- if (s32pc[ii,jj]>0) s12pc[ii,jj]/s32pc[ii,jj] else NA_real_
      }
      long_rows <- list()
      for (ii in seq_len(n_pops-1L)) for (jj in seq(ii+1L,n_pops))
        long_rows[[length(long_rows)+1L]] <- data.frame(
          Pop1=pops[ii], Pop2=pops[jj],
          FST_raw=round(mat_raw[jj,ii],6), FST_ENA=round(mat_ena[jj,ii],6),
          Delta_FST=round(mat_ena[jj,ii]-mat_raw[jj,ii],6), stringsAsFactors=FALSE)
      list(matrix_raw=mat_raw, matrix_ena=mat_ena, long=do.call(rbind,long_rows))
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 5 — PAIRWISE DCSE (raw + INA)
    #  FIX: Pascal prod_CS_correction — ajustement_r ALWAYS appends rd as null
    #  allele state, even when pfreq is empty. Check only ni_ci>0 && ni_cj>0.
    # ══════════════════════════════════════════════════════════════════════════
    compute_dc_pairwise <- function(em_res) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]]); n_pops <- length(pops)
      if (n_pops < 2L) return(list(matrix_raw=NULL,matrix_ina=NULL,long=data.frame()))
      dc_sum_raw <- matrix(0.0,n_pops,n_pops,dimnames=list(pops,pops))
      dc_sum_ina <- matrix(0.0,n_pops,n_pops,dimnames=list(pops,pops))
      nloc_eff   <- matrix(length(markers),n_pops,n_pops,dimnames=list(pops,pops))
      nloc_eff_c <- matrix(length(markers),n_pops,n_pops,dimnames=list(pops,pops))
      for (loc in markers) {
        em_loc <- em_res[[loc]]
        for (ii in seq_len(n_pops-1L)) {
          for (jj in seq(ii+1L,n_pops)) {
            ei <- em_loc[[pops[ii]]]; ej <- em_loc[[pops[jj]]]
            # Raw DCSE: genefreq_obs, null homos excluded — Pascal prod_CS
            ni_ri <- ei$efpop - ei$n_absent - ei$n_null_homo
            ni_rj <- ej$efpop - ej$n_absent - ej$n_null_homo
            if (ni_ri > 0L && ni_rj > 0L && length(ei$genefreq_obs) > 0 && length(ej$genefreq_obs) > 0) {
              d_raw <- cs_distance(ei$genefreq_obs, ej$genefreq_obs)
              if (!is.na(d_raw)) dc_sum_raw[jj,ii] <- dc_sum_raw[jj,ii] + d_raw
              else nloc_eff[jj,ii] <- nloc_eff[jj,ii] - 1L
            } else {
              nloc_eff[jj,ii] <- nloc_eff[jj,ii] - 1L
            }
            # INA DCSE: corrdgenefreq + rd as null allele — Pascal prod_CS_correction
            # FIX: only check ni_ci>0 && ni_cj>0, NOT length(pfreq)>0
            # ajustement_r always sets corrdgenefreq[...,nallc-1] = rd[iloc,ipop]
            ni_ci <- ei$efpop - ei$n_absent
            ni_cj <- ej$efpop - ej$n_absent
            if (ni_ci > 0L && ni_cj > 0L) {
              freq_ina_i <- make_ina_freq(ei)   # c(pfreq, __null__=rd)
              freq_ina_j <- make_ina_freq(ej)
              d_ina <- cs_distance(freq_ina_i, freq_ina_j)
              if (!is.na(d_ina)) dc_sum_ina[jj,ii] <- dc_sum_ina[jj,ii] + d_ina
              else nloc_eff_c[jj,ii] <- nloc_eff_c[jj,ii] - 1L
            } else {
              nloc_eff_c[jj,ii] <- nloc_eff_c[jj,ii] - 1L
            }
          }
        }
      }
      mat_raw <- matrix(NA_real_,n_pops,n_pops,dimnames=list(pops,pops))
      mat_ina <- matrix(NA_real_,n_pops,n_pops,dimnames=list(pops,pops))
      for (ii in seq_len(n_pops-1L)) for (jj in seq(ii+1L,n_pops)) {
        mat_raw[jj,ii] <- if (nloc_eff[jj,ii]  >0L) dc_sum_raw[jj,ii]/nloc_eff[jj,ii]   else NA_real_
        mat_ina[jj,ii] <- if (nloc_eff_c[jj,ii]>0L) dc_sum_ina[jj,ii]/nloc_eff_c[jj,ii] else NA_real_
      }
      long_rows <- list()
      for (ii in seq_len(n_pops-1L)) for (jj in seq(ii+1L,n_pops))
        long_rows[[length(long_rows)+1L]] <- data.frame(
          Pop1=pops[ii], Pop2=pops[jj],
          DCSE_raw=round(mat_raw[jj,ii],6), DCSE_INA=round(mat_ina[jj,ii],6),
          Delta_DCSE=round(mat_ina[jj,ii]-mat_raw[jj,ii],6), stringsAsFactors=FALSE)
      list(matrix_raw=mat_raw, matrix_ina=mat_ina, long=do.call(rbind,long_rows))
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  TAB 6 — PER-LOCUS FST FOR EACH PAIR
    # ══════════════════════════════════════════════════════════════════════════
    compute_fst_per_locus_pair <- function(em_res, sel_locus="all",
                                           sel_pop1="all", sel_pop2="all") {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]])
      if (!identical(sel_locus,"all")) markers <- markers[markers==sel_locus]
      all_pairs <- combn(pops,2,simplify=FALSE)
      pairs <- if (!identical(sel_pop1,"all")&&!identical(sel_pop2,"all"))
                 list(c(sel_pop1,sel_pop2))
               else if (!identical(sel_pop1,"all")) Filter(function(x) sel_pop1 %in% x, all_pairs)
               else if (!identical(sel_pop2,"all")) Filter(function(x) sel_pop2 %in% x, all_pairs)
               else all_pairs
      rows <- list()
      for (loc in markers) {
        em_loc <- em_res[[loc]]
        alleles_obs <- sort(unique(unlist(lapply(em_loc, function(e) e$alleles))))
        for (pair in pairs) {
          pi_n <- pair[1]; pj_n <- pair[2]
          if (!pi_n %in% pops || !pj_n %in% pops) next
          ei <- em_loc[[pi_n]]; ej <- em_loc[[pj_n]]
          ni_ri <- max(0L,ei$efpop-ei$n_absent-ei$n_null_homo)
          ni_rj <- max(0L,ej$efpop-ej$n_absent-ej$n_null_homo)
          ni_ci <- max(0L,ei$efpop-ei$n_absent)
          ni_cj <- max(0L,ej$efpop-ej$n_absent)
          N_r   <- ni_ri+ni_rj; N_c <- ni_ci+ni_cj
          nc_r  <- if (N_r>0&&ni_ri>0&&ni_rj>0) (N_r-(ni_ri^2+ni_rj^2)/N_r) else 0.0
          nc_c  <- if (N_c>0&&ni_ci>0&&ni_cj>0) (N_c-(ni_ci^2+ni_cj^2)/N_c) else 0.0
          s1_r <- s3_r <- s1_c <- s3_c <- 0.0
          for (a in alleles_obs) {
            a_chr <- as.character(a)
            if (nc_r > 0) {
              pd <- list(
                list(ni=ni_ri, nA=(if(a_chr %in% names(ei$genefreq_obs)) ei$genefreq_obs[[a_chr]] else 0.0)*2L*ni_ri,
                     AA=if(a_chr %in% names(ei$H_ii)) ei$H_ii[[a_chr]] else 0L, AA_corr=0.0),
                list(ni=ni_rj, nA=(if(a_chr %in% names(ej$genefreq_obs)) ej$genefreq_obs[[a_chr]] else 0.0)*2L*ni_rj,
                     AA=if(a_chr %in% names(ej$H_ii)) ej$H_ii[[a_chr]] else 0L, AA_corr=0.0))
              cmp <- weir_components_allele(pd, use_corr=FALSE)
              s1_r <- s1_r+cmp$s2P*nc_r; s3_r <- s3_r+(cmp$s2P+cmp$s2I+cmp$s2G)*nc_r
            }
            if (nc_c > 0) {
              pf_i <- if(a_chr %in% names(ei$pfreq)) ei$pfreq[[a_chr]] else 0.0
              pf_j <- if(a_chr %in% names(ej$pfreq)) ej$pfreq[[a_chr]] else 0.0
              AA_i <- if(a_chr %in% names(ei$H_ii)) ei$H_ii[[a_chr]] else 0L
              AA_j <- if(a_chr %in% names(ej$H_ii)) ej$H_ii[[a_chr]] else 0L
              di   <- pf_i+2.0*ei$rd; dj <- pf_j+2.0*ej$rd
              pdc  <- list(
                list(ni=ni_ci, nA=pf_i*2L*ni_ci, AA=AA_i, AA_corr=if(AA_i>0&&di>0) AA_i*(pf_i/di) else 0.0),
                list(ni=ni_cj, nA=pf_j*2L*ni_cj, AA=AA_j, AA_corr=if(AA_j>0&&dj>0) AA_j*(pf_j/dj) else 0.0))
              cmpc <- weir_components_allele(pdc, use_corr=TRUE)
              s1_c <- s1_c+cmpc$s2P*nc_c; s3_c <- s3_c+(cmpc$s2P+cmpc$s2I+cmpc$s2G)*nc_c
            }
          }
          rows[[length(rows)+1L]] <- data.frame(
            Locus=loc, Pop1=pi_n, Pop2=pj_n,
            FST_raw=round(if(s3_r!=0)s1_r/s3_r else NA_real_,6),
            FST_ENA=round(if(s3_c!=0)s1_c/s3_c else NA_real_,6),
            Delta=round((if(s3_c!=0)s1_c/s3_c else NA_real_)-(if(s3_r!=0)s1_r/s3_r else NA_real_),6),
            N_i_raw=ni_ri, N_j_raw=ni_rj, N_i_ENA=ni_ci, N_j_ENA=ni_cj,
            stringsAsFactors=FALSE)
        }
      }
      if (length(rows)==0L) return(data.frame())
      do.call(rbind,rows)
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  PRE-COMPUTE PER-LOCUS STATS (for fast vectorised bootstrap over loci)
    # ══════════════════════════════════════════════════════════════════════════
    precompute_locus_stats_fst_global <- function(em_res) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]])
      s1v <- s3v <- s1cv <- s3cv <- numeric(length(markers))
      for (li in seq_along(markers)) {
        loc <- markers[li]; em_loc <- em_res[[loc]]
        alleles_obs <- sort(unique(unlist(lapply(em_loc, function(e) e$alleles))))
        for (a in alleles_obs) {
          a_chr <- as.character(a)
          pop_raw <- lapply(pops, function(p) {
            e <- em_loc[[p]]; ni <- max(0L,e$efpop-e$n_absent-e$n_null_homo)
            pf <- if (a_chr %in% names(e$genefreq_obs)) e$genefreq_obs[[a_chr]] else 0.0
            AA <- if (a_chr %in% names(e$H_ii)) e$H_ii[[a_chr]] else 0L
            list(ni=ni, nA=pf*2L*ni, AA=AA, AA_corr=AA)
          })
          cmp <- weir_components_allele(pop_raw, use_corr=FALSE)
          s1v[li]  <- s1v[li] +cmp$s2P; s3v[li]  <- s3v[li] +cmp$s2P+cmp$s2I+cmp$s2G
          pop_ena <- lapply(pops, function(p) {
            e <- em_loc[[p]]; ni <- max(0L,e$efpop-e$n_absent)
            pf <- if (a_chr %in% names(e$pfreq)) e$pfreq[[a_chr]] else 0.0
            AA <- if (a_chr %in% names(e$H_ii)) e$H_ii[[a_chr]] else 0L
            d  <- pf+2.0*e$rd; AAc <- if (AA>0&&d>0) AA*(pf/d) else 0.0
            list(ni=ni, nA=pf*2L*ni, AA=AA, AA_corr=AAc)
          })
          cmpc <- weir_components_allele(pop_ena, use_corr=TRUE)
          s1cv[li] <- s1cv[li]+cmpc$s2P; s3cv[li] <- s3cv[li]+cmpc$s2P+cmpc$s2I+cmpc$s2G
        }
      }
      list(s1=s1v, s3=s3v, s1c=s1cv, s3c=s3cv, markers=markers)
    }

    precompute_locus_stats_fst_pair <- function(em_res) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]])
      pairs   <- combn(pops,2,simplify=FALSE); n_pairs <- length(pairs)
      s1_raw <- matrix(0.0,n_pairs,length(markers))
      s3_raw <- matrix(0.0,n_pairs,length(markers))
      s1_ena <- matrix(0.0,n_pairs,length(markers))
      s3_ena <- matrix(0.0,n_pairs,length(markers))
      for (li in seq_along(markers)) {
        loc <- markers[li]; em_loc <- em_res[[loc]]
        alleles_obs <- sort(unique(unlist(lapply(em_loc, function(e) e$alleles))))
        for (pi in seq_along(pairs)) {
          pi_n <- pairs[[pi]][1]; pj_n <- pairs[[pi]][2]
          ei <- em_loc[[pi_n]]; ej <- em_loc[[pj_n]]
          ni_ri <- max(0L,ei$efpop-ei$n_absent-ei$n_null_homo)
          ni_rj <- max(0L,ej$efpop-ej$n_absent-ej$n_null_homo)
          ni_ci <- max(0L,ei$efpop-ei$n_absent)
          ni_cj <- max(0L,ej$efpop-ej$n_absent)
          N_r   <- ni_ri+ni_rj; N_c <- ni_ci+ni_cj
          nc_r  <- if (N_r>0&&ni_ri>0&&ni_rj>0) (N_r-(ni_ri^2+ni_rj^2)/N_r) else 0.0
          nc_c  <- if (N_c>0&&ni_ci>0&&ni_cj>0) (N_c-(ni_ci^2+ni_cj^2)/N_c) else 0.0
          for (a in alleles_obs) {
            a_chr <- as.character(a)
            if (nc_r > 0) {
              pd <- list(
                list(ni=ni_ri, nA=(if(a_chr %in% names(ei$genefreq_obs)) ei$genefreq_obs[[a_chr]] else 0.0)*2L*ni_ri,
                     AA=if(a_chr %in% names(ei$H_ii)) ei$H_ii[[a_chr]] else 0L, AA_corr=0.0),
                list(ni=ni_rj, nA=(if(a_chr %in% names(ej$genefreq_obs)) ej$genefreq_obs[[a_chr]] else 0.0)*2L*ni_rj,
                     AA=if(a_chr %in% names(ej$H_ii)) ej$H_ii[[a_chr]] else 0L, AA_corr=0.0))
              cmp <- weir_components_allele(pd, use_corr=FALSE)
              s1_raw[pi,li] <- s1_raw[pi,li]+cmp$s2P*nc_r
              s3_raw[pi,li] <- s3_raw[pi,li]+(cmp$s2P+cmp$s2I+cmp$s2G)*nc_r
            }
            if (nc_c > 0) {
              pf_i <- if(a_chr %in% names(ei$pfreq)) ei$pfreq[[a_chr]] else 0.0
              pf_j <- if(a_chr %in% names(ej$pfreq)) ej$pfreq[[a_chr]] else 0.0
              AA_i <- if(a_chr %in% names(ei$H_ii)) ei$H_ii[[a_chr]] else 0L
              AA_j <- if(a_chr %in% names(ej$H_ii)) ej$H_ii[[a_chr]] else 0L
              di   <- pf_i+2.0*ei$rd; dj <- pf_j+2.0*ej$rd
              pdc  <- list(
                list(ni=ni_ci, nA=pf_i*2L*ni_ci, AA=AA_i, AA_corr=if(AA_i>0&&di>0) AA_i*(pf_i/di) else 0.0),
                list(ni=ni_cj, nA=pf_j*2L*ni_cj, AA=AA_j, AA_corr=if(AA_j>0&&dj>0) AA_j*(pf_j/dj) else 0.0))
              cmpc <- weir_components_allele(pdc, use_corr=TRUE)
              s1_ena[pi,li] <- s1_ena[pi,li]+cmpc$s2P*nc_c
              s3_ena[pi,li] <- s3_ena[pi,li]+(cmpc$s2P+cmpc$s2I+cmpc$s2G)*nc_c
            }
          }
        }
      }
      list(s1_raw=s1_raw, s3_raw=s3_raw, s1_ena=s1_ena, s3_ena=s3_ena,
           markers=markers, pairs=pairs)
    }

    # FIX: precompute_locus_stats_dc — use make_ina_freq() (ni_ci>0 check only)
    precompute_locus_stats_dc <- function(em_res) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]])
      pairs   <- combn(pops,2,simplify=FALSE); n_pairs <- length(pairs)
      dc_raw_m <- matrix(NA_real_, n_pairs, length(markers))
      dc_ina_m <- matrix(NA_real_, n_pairs, length(markers))
      for (li in seq_along(markers)) {
        loc <- markers[li]; em_loc <- em_res[[loc]]
        for (pi in seq_along(pairs)) {
          ei <- em_loc[[pairs[[pi]][1]]]; ej <- em_loc[[pairs[[pi]][2]]]
          ni_ri <- ei$efpop-ei$n_absent-ei$n_null_homo
          ni_rj <- ej$efpop-ej$n_absent-ej$n_null_homo
          ni_ci <- ei$efpop-ei$n_absent; ni_cj <- ej$efpop-ej$n_absent
          # Raw: require observed alleles
          if (ni_ri>0 && ni_rj>0 && length(ei$genefreq_obs)>0 && length(ej$genefreq_obs)>0)
            dc_raw_m[pi,li] <- cs_distance(ei$genefreq_obs, ej$genefreq_obs)
          # INA: FIX — only check ni_ci>0, always include rd as null allele
          if (ni_ci>0 && ni_cj>0)
            dc_ina_m[pi,li] <- cs_distance(make_ina_freq(ei), make_ina_freq(ej))
        }
      }
      list(dc_raw=dc_raw_m, dc_ina=dc_ina_m, markers=markers, pairs=pairs)
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  BOOTSTRAP OVER LOCI — fully vectorised
    # ══════════════════════════════════════════════════════════════════════════
    boot_loci_fst_global <- function(lstat, nboot) {
      L   <- length(lstat$markers)
      idx <- matrix(sample.int(L, L*nboot, replace=TRUE), nrow=nboot)
      S1  <- matrix(lstat$s1[idx],  nrow=nboot); S3  <- matrix(lstat$s3[idx],  nrow=nboot)
      S1c <- matrix(lstat$s1c[idx], nrow=nboot); S3c <- matrix(lstat$s3c[idx], nrow=nboot)
      fst_raw_boot <- ifelse(rowSums(S3) >0, rowSums(S1) /rowSums(S3),  NA_real_)
      fst_ena_boot <- ifelse(rowSums(S3c)>0, rowSums(S1c)/rowSums(S3c), NA_real_)
      # Per-locus jackknife (leave-one-out)
      fst_obs_all <- if (sum(lstat$s3c)>0) sum(lstat$s1c)/sum(lstat$s3c) else NA_real_
      per_locus <- do.call(rbind, lapply(seq_len(L), function(l) {
        s1c_jk <- sum(lstat$s1c[-l]); s3c_jk <- sum(lstat$s3c[-l])
        fst_jk <- if (s3c_jk>0) s1c_jk/s3c_jk else NA_real_
        pseudo  <- L * fst_obs_all - (L-1) * fst_jk
        data.frame(
          Locus        = lstat$markers[l],
          FST_ENA_obs  = round(if(lstat$s3c[l]>0) lstat$s1c[l]/lstat$s3c[l] else NA_real_, 6),
          JK_pseudoval = round(pseudo, 6),
          stringsAsFactors=FALSE)
      }))
      list(raw=quantile(fst_raw_boot,c(0.025,0.5,0.975),na.rm=TRUE),
           ena=quantile(fst_ena_boot,c(0.025,0.5,0.975),na.rm=TRUE),
           per_locus=per_locus, dist_ena=fst_ena_boot)
    }

    # FIX: per_locus data.frame — Locus stored as character, not coerced
    boot_loci_fst_pair <- function(lstat, nboot) {
      L       <- length(lstat$markers)
      n_pairs <- nrow(lstat$s1_raw)
      idx     <- matrix(sample.int(L, L*nboot, replace=TRUE), nrow=nboot)
      results <- vector("list", n_pairs)
      for (pi in seq_len(n_pairs)) {
        s1r <- lstat$s1_raw[pi,]; s3r <- lstat$s3_raw[pi,]
        s1e <- lstat$s1_ena[pi,]; s3e <- lstat$s3_ena[pi,]
        RS1r <- rowSums(matrix(s1r[idx],nrow=nboot)); RS3r <- rowSums(matrix(s3r[idx],nrow=nboot))
        RS1e <- rowSums(matrix(s1e[idx],nrow=nboot)); RS3e <- rowSums(matrix(s3e[idx],nrow=nboot))
        boot_raw <- ifelse(RS3r>0, RS1r/RS3r, NA_real_)
        boot_ena <- ifelse(RS3e>0, RS1e/RS3e, NA_real_)
        fst_ena_obs <- if (sum(s3e)>0) sum(s1e)/sum(s3e) else NA_real_
        fst_raw_obs <- if (sum(s3r)>0) sum(s1r)/sum(s3r) else NA_real_
        # FIX: explicit character() for Locus column — avoid apply() coercion
        per_loc <- data.frame(
          Locus          = as.character(lstat$markers),
          Pop1           = lstat$pairs[[pi]][1],
          Pop2           = lstat$pairs[[pi]][2],
          FST_ENA_locus  = round(ifelse(s3e>0, s1e/s3e, NA_real_), 6),
          FST_raw_locus  = round(ifelse(s3r>0, s1r/s3r, NA_real_), 6),
          stringsAsFactors=FALSE)
        results[[pi]] <- list(
          summary=data.frame(
            Pop1=lstat$pairs[[pi]][1], Pop2=lstat$pairs[[pi]][2],
            FST_ENA_obs=round(fst_ena_obs,6),
            CI_lo_loci  =round(quantile(boot_ena,0.025,na.rm=TRUE),6),
            Median_loci =round(quantile(boot_ena,0.500,na.rm=TRUE),6),
            CI_hi_loci  =round(quantile(boot_ena,0.975,na.rm=TRUE),6),
            FST_raw_obs =round(fst_raw_obs,6),
            CI_lo_raw   =round(quantile(boot_raw,0.025,na.rm=TRUE),6),
            CI_hi_raw   =round(quantile(boot_raw,0.975,na.rm=TRUE),6),
            stringsAsFactors=FALSE),
          per_locus=per_loc)
      }
      list(summary  =do.call(rbind,lapply(results,`[[`,"summary")),
           per_locus=do.call(rbind,lapply(results,`[[`,"per_locus")))
    }

    # FIX: per_locus data.frame — Locus as character, use make_ina_freq for INA
    boot_loci_dc_pair <- function(lstat, nboot) {
      L       <- length(lstat$markers)
      n_pairs <- nrow(lstat$dc_raw)
      idx     <- matrix(sample.int(L, L*nboot, replace=TRUE), nrow=nboot)
      results <- vector("list", n_pairs)
      for (pi in seq_len(n_pairs)) {
        dr <- lstat$dc_raw[pi,]; di <- lstat$dc_ina[pi,]
        boot_raw <- rowMeans(matrix(dr[idx],nrow=nboot), na.rm=TRUE)
        boot_ina <- rowMeans(matrix(di[idx],nrow=nboot), na.rm=TRUE)
        # FIX: explicit character() for Locus column
        per_loc <- data.frame(
          Locus          = as.character(lstat$markers),
          Pop1           = lstat$pairs[[pi]][1],
          Pop2           = lstat$pairs[[pi]][2],
          DCSE_INA_locus = round(di, 6),
          DCSE_raw_locus = round(dr, 6),
          stringsAsFactors=FALSE)
        results[[pi]] <- list(
          summary=data.frame(
            Pop1=lstat$pairs[[pi]][1], Pop2=lstat$pairs[[pi]][2],
            DCSE_INA_obs =round(mean(di,na.rm=TRUE),6),
            CI_lo_loci   =round(quantile(boot_ina,0.025,na.rm=TRUE),6),
            Median_loci  =round(quantile(boot_ina,0.500,na.rm=TRUE),6),
            CI_hi_loci   =round(quantile(boot_ina,0.975,na.rm=TRUE),6),
            DCSE_raw_obs =round(mean(dr,na.rm=TRUE),6),
            CI_lo_raw    =round(quantile(boot_raw,0.025,na.rm=TRUE),6),
            CI_hi_raw    =round(quantile(boot_raw,0.975,na.rm=TRUE),6),
            stringsAsFactors=FALSE),
          per_locus=per_loc)
      }
      list(summary  =do.call(rbind,lapply(results,`[[`,"summary")),
           per_locus=do.call(rbind,lapply(results,`[[`,"per_locus")))
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  BOOTSTRAP OVER INDIVIDUALS
    # ══════════════════════════════════════════════════════════════════════════
    boot_indiv_fst_pair <- function(raw_df, em_res, base, treatments, nboot) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]])
      pairs   <- combn(pops,2,simplify=FALSE)
      idx_by_pop_loc <- list()
      for (pop in pops) for (loc in markers)
        idx_by_pop_loc[[paste0(pop,"___",loc)]] <-
          which(raw_df$Population==pop & raw_df$Marker==loc)
      inds_by_pop <- lapply(pops, function(p) unique(raw_df$Individual[raw_df$Population==p]))
      names(inds_by_pop) <- pops
      boot_fst <- matrix(NA_real_, nboot, length(pairs))
      for (b in seq_len(nboot)) {
        em_b <- list()
        for (loc in markers) {
          em_b[[loc]] <- list()
          treat <- as.character(treatments[loc] %||% "null_homo")
          for (pop in pops) {
            inds <- inds_by_pop[[pop]]
            if (length(inds)==0L) { em_b[[loc]][[pop]] <- em_res[[loc]][[pop]]; next }
            resampled_inds <- sample(inds, length(inds), replace=TRUE)
            base_rows <- idx_by_pop_loc[[paste0(pop,"___",loc)]]
            ind_col   <- raw_df$Individual[base_rows]
            idx_rows  <- unlist(lapply(resampled_inds, function(ind) base_rows[ind_col==ind]))
            gts <- raw_df$gt[idx_rows]
            em_b[[loc]][[pop]] <- if (length(gts)==0L) em_res[[loc]][[pop]]
                                  else em_freena(gts, base, treat)
          }
        }
        pw <- compute_fst_pairwise(em_b)
        for (pi in seq_along(pairs))
          boot_fst[b,pi] <- pw$matrix_ena[pairs[[pi]][2], pairs[[pi]][1]]
      }
      obs_long <- compute_fst_pairwise(em_res)$long
      results  <- vector("list", length(pairs))
      for (pi in seq_along(pairs)) {
        v <- boot_fst[,pi]; p1 <- pairs[[pi]][1]; p2 <- pairs[[pi]][2]
        obs_row <- obs_long[obs_long$Pop1==p1 & obs_long$Pop2==p2,,drop=FALSE]
        fst_obs <- if (nrow(obs_row)>0) obs_row$FST_ENA[1] else NA_real_
        results[[pi]] <- data.frame(Pop1=p1, Pop2=p2,
          FST_ENA_obs =round(fst_obs,6),
          CI_lo_indiv =round(quantile(v,0.025,na.rm=TRUE),6),
          Median_indiv=round(quantile(v,0.500,na.rm=TRUE),6),
          CI_hi_indiv =round(quantile(v,0.975,na.rm=TRUE),6),
          stringsAsFactors=FALSE)
      }
      do.call(rbind, results)
    }

    boot_indiv_dc_pair <- function(raw_df, em_res, base, treatments, nboot) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]])
      pairs   <- combn(pops,2,simplify=FALSE)
      idx_by_pop_loc <- list()
      for (pop in pops) for (loc in markers)
        idx_by_pop_loc[[paste0(pop,"___",loc)]] <-
          which(raw_df$Population==pop & raw_df$Marker==loc)
      inds_by_pop <- lapply(pops, function(p) unique(raw_df$Individual[raw_df$Population==p]))
      names(inds_by_pop) <- pops
      boot_dc <- matrix(NA_real_, nboot, length(pairs))
      for (b in seq_len(nboot)) {
        em_b <- list()
        for (loc in markers) {
          em_b[[loc]] <- list()
          treat <- as.character(treatments[loc] %||% "null_homo")
          for (pop in pops) {
            inds <- inds_by_pop[[pop]]
            if (length(inds)==0L) { em_b[[loc]][[pop]] <- em_res[[loc]][[pop]]; next }
            resampled_inds <- sample(inds, length(inds), replace=TRUE)
            base_rows <- idx_by_pop_loc[[paste0(pop,"___",loc)]]
            ind_col   <- raw_df$Individual[base_rows]
            idx_rows  <- unlist(lapply(resampled_inds, function(ind) base_rows[ind_col==ind]))
            gts <- raw_df$gt[idx_rows]
            em_b[[loc]][[pop]] <- if (length(gts)==0L) em_res[[loc]][[pop]]
                                  else em_freena(gts, base, treat)
          }
        }
        dc_b <- compute_dc_pairwise(em_b)
        for (pi in seq_along(pairs))
          boot_dc[b,pi] <- dc_b$matrix_ina[pairs[[pi]][2], pairs[[pi]][1]]
      }
      obs_long <- compute_dc_pairwise(em_res)$long
      results  <- vector("list", length(pairs))
      for (pi in seq_along(pairs)) {
        v <- boot_dc[,pi]; p1 <- pairs[[pi]][1]; p2 <- pairs[[pi]][2]
        obs_row <- obs_long[obs_long$Pop1==p1 & obs_long$Pop2==p2,,drop=FALSE]
        dc_obs  <- if (nrow(obs_row)>0) obs_row$DCSE_INA[1] else NA_real_
        results[[pi]] <- data.frame(Pop1=p1, Pop2=p2,
          DCSE_INA_obs=round(dc_obs,6),
          CI_lo_indiv =round(quantile(v,0.025,na.rm=TRUE),6),
          Median_indiv=round(quantile(v,0.500,na.rm=TRUE),6),
          CI_hi_indiv =round(quantile(v,0.975,na.rm=TRUE),6),
          stringsAsFactors=FALSE)
      }
      do.call(rbind, results)
    }

    # ── Event-triggered reactives ──────────────────────────────────────────────
    em_r <- reactive({
      db_ready()
      withProgress(message="EM FreeNA — computing null allele frequencies...", value=0.1,
        { res <- em_results_r(); setProgress(1); res })
    })
    fst_global_r <- eventReactive(input$run_fst_global, {
      req(length(em_r())>0)
      withProgress(message="Computing global FST (ENA)...", value=0.2,
        { res <- compute_fst_global(em_r()); setProgress(1); res })
    })
    fst_pair_r <- eventReactive(input$run_fst_pair, {
      req(length(em_r())>0)
      withProgress(message="Computing pairwise FST (ENA)...", value=0.2,
        { res <- compute_fst_pairwise(em_r()); setProgress(1); res })
    })
    dc_r <- eventReactive(input$run_dc, {
      req(length(em_r())>0)
      withProgress(message="Computing pairwise DCSE (INA)...", value=0.2,
        { res <- compute_dc_pairwise(em_r()); setProgress(1); res })
    })
    fst_locus_r <- eventReactive(input$run_fst_locus, {
      req(length(em_r())>0)
      withProgress(message="Computing per-locus FST per pair...", value=0.2, {
        res <- compute_fst_per_locus_pair(em_r(),
          sel_locus=safe_choice(input$fl_locus,"all"),
          sel_pop1 =safe_choice(input$fl_pop1, "all"),
          sel_pop2 =safe_choice(input$fl_pop2, "all"))
        setProgress(1); res })
    })
    boot_fst_global_r <- eventReactive(input$run_boot_fst_global, {
      req(length(em_r())>0)
      nboot <- max(999L,min(9999L,as.integer(input$run_boot_fst_global_nboot %||% 5000L)))
      withProgress(message=sprintf("Bootstrap FST-ENA global — %d replicates...",nboot),value=0.1,{
        lstat <- precompute_locus_stats_fst_global(em_r()); setProgress(0.4)
        res   <- boot_loci_fst_global(lstat,nboot); setProgress(1); res
      })
    })
    boot_fst_pair_r <- eventReactive(input$run_boot_fst_pair, {
      req(length(em_r())>0)
      nboot <- max(999L,min(9999L,as.integer(input$run_boot_fst_pair_nboot %||% 5000L)))
      btype <- input$run_boot_fst_pair_type %||% "both_boot"
      withProgress(message=sprintf("Bootstrap pairwise FST-ENA — %d replicates...",nboot),value=0.1,{
        res <- list()
        if (btype %in% c("loci","both_boot")) {
          lstat <- precompute_locus_stats_fst_pair(em_r()); setProgress(0.35)
          res$loci <- boot_loci_fst_pair(lstat,nboot); setProgress(0.6)
        }
        if (btype %in% c("indiv","both_boot")) {
          raw_df <- raw_data_r()
          res$indiv <- boot_indiv_fst_pair(raw_df,em_r(),as.integer(base_r()),
                                           locus_treatments_r(),nboot)
          setProgress(0.95)
        }
        setProgress(1); res
      })
    })
    boot_dc_r <- eventReactive(input$run_boot_dc, {
      req(length(em_r())>0)
      nboot <- max(999L,min(9999L,as.integer(input$run_boot_dc_nboot %||% 5000L)))
      btype <- input$run_boot_dc_type %||% "both_boot"
      withProgress(message=sprintf("Bootstrap pairwise DCSE-INA — %d replicates...",nboot),value=0.1,{
        res <- list()
        if (btype %in% c("loci","both_boot")) {
          lstat <- precompute_locus_stats_dc(em_r()); setProgress(0.35)
          res$loci <- boot_loci_dc_pair(lstat,nboot); setProgress(0.6)
        }
        if (btype %in% c("indiv","both_boot")) {
          raw_df <- raw_data_r()
          res$indiv <- boot_indiv_dc_pair(raw_df,em_r(),as.integer(base_r()),
                                          locus_treatments_r(),nboot)
          setProgress(0.95)
        }
        setProgress(1); res
      })
    })

    # ── Value boxes ────────────────────────────────────────────────────────────
    output$vb_loci <- renderUI({ tryCatch(tags$span(length(markers_r())),error=function(e)tags$span("\u2014")) })
    output$vb_pops <- renderUI({ tryCatch(tags$span(length(pops_r())),  error=function(e)tags$span("\u2014")) })
    output$vb_n    <- renderUI({
      tryCatch({
        db_ready(); con <- con_r(); ms <- meta_schema_r()
        n <- DBI::dbGetQuery(con, sprintf(
          "SELECT COUNT(DISTINCT CAST(%s AS VARCHAR)) AS n FROM %s WHERE %s IS NOT NULL",
          sql_id(con,ms$ind_col),sql_id(con,tbl_meta_r()),sql_id(con,ms$ind_col)))$n[[1]]
        tags$span(n)
      }, error=function(e) tags$span("\u2014"))
    })
    output$vb_avg_null <- renderUI({
      tryCatch({
        d <- t1_data_r()
        if (nrow(d)==0||all(is.na(d$p_nulls))) return(tags$span("\u2014"))
        v   <- round(mean(d$p_nulls,na.rm=TRUE),4)
        col <- if(v>.20)"#9d174d" else if(v>.10)"#854d0e" else "#166534"
        tags$span(style=paste0("color:",col,";"),v)
      }, error=function(e) tags$span("\u2014"))
    })
    output$vb_max_null <- renderUI({
      tryCatch({
        d <- t1_data_r()
        if (nrow(d)==0||all(is.na(d$p_nulls))) return(tags$span("\u2014"))
        v   <- round(max(d$p_nulls,na.rm=TRUE),4)
        col <- if(v>.30)"#9d174d" else if(v>.15)"#854d0e" else "#166534"
        tags$span(style=paste0("color:",col,";"),v)
      }, error=function(e) tags$span("\u2014"))
    })
    output$vb_fst_ena <- renderUI({
      tryCatch({
        r <- fst_global_r(); v <- round(r$global_ena,4)
        col <- if(!is.na(v)&&v>.15)"#9d174d" else if(!is.na(v)&&v>.05)"#854d0e" else "#166534"
        tags$span(style=paste0("color:",col,";"),if(is.na(v))"\u2014" else v)
      }, error=function(e) tags$span("\u2014"))
    })

    # ── Tab 1 & 2 DT ──────────────────────────────────────────────────────────
    output$dt_t1 <- DT::renderDT({
      d <- t1_data_r()
      shiny::validate(shiny::need(nrow(d)>0,"No data. Select parameters and click Compute."))
      disp <- d; names(disp) <- c("Locus names","Farm","p_nulls","N","N_exp_blanks","p_nulls\u00d7N")
      DT::datatable(disp,rownames=FALSE,
        options=list(pageLength=20,scrollX=TRUE,dom="lftip",columnDefs=list(list(className="dt-right",targets=2:5))),
        class="compact hover stripe") |>
        DT::formatRound("p_nulls",5)|>DT::formatRound("N_exp_blanks",9)|>DT::formatRound("p_nulls\u00d7N",5)|>
        DT::formatStyle("p_nulls",backgroundColor=DT::styleInterval(c(0.05,0.10,0.20,0.30),c("#f0fdf4","#dcfce7","#fefce8","#fff7ed","#fef2f2")))|>
        DT::formatStyle("Locus names",fontWeight="600",color="#0f172a")|>
        DT::formatStyle("Farm",color="#475569")
    }, server=TRUE)

    output$dt_t2 <- DT::renderDT({
      d <- t2_data_r()
      shiny::validate(shiny::need(nrow(d)>0,"No data. Select parameters and click Compute."))
      disp <- d; names(disp) <- c("Locus names","Av(N_exp_blanks)","Av(p_nulls)","N_tot","N_blanks","f(expBlanks)","p_nulls")
      DT::datatable(disp,rownames=FALSE,
        options=list(pageLength=20,scrollX=TRUE,dom="lftip",columnDefs=list(list(className="dt-right",targets=1:6))),
        class="compact hover stripe") |>
        DT::formatRound("Av(N_exp_blanks)",9)|>DT::formatRound("Av(p_nulls)",9)|>
        DT::formatRound("f(expBlanks)",9)|>DT::formatRound("p_nulls",9)|>
        DT::formatStyle("p_nulls",backgroundColor=DT::styleInterval(c(0.05,0.10,0.20),c("#f0fdf4","#dcfce7","#fefce8","#fef2f2")))|>
        DT::formatStyle("Locus names",fontWeight="600",color="#0f172a")
    }, server=TRUE)

    # ── Tab 3 DT + bootstrap ───────────────────────────────────────────────────
    output$dt_fst_global <- DT::renderDT({
      r <- fst_global_r(); d <- r$per_locus
      shiny::validate(shiny::need(nrow(d)>0,"No data. Click Compute."))
      summary_row <- data.frame(
        Locus=paste0("[Multilocus  Raw=",round(r$global_raw,6),"  |  ENA=",round(r$global_ena,6),"]"),
        FST_raw=r$global_raw, FST_ENA=r$global_ena, Delta_FST=r$global_ena-r$global_raw,
        N_pops_eff_raw=NA_integer_, N_pops_eff_ENA=NA_integer_, stringsAsFactors=FALSE)
      disp <- rbind(summary_row,d)
      names(disp) <- c("Locus","Raw FST","FST-ENA","\u0394FST (ENA\u2212raw)","N eff. pops (raw)","N eff. pops (ENA)")
      DT::datatable(disp,rownames=FALSE,
        options=list(pageLength=25,scrollX=TRUE,dom="lftip",columnDefs=list(list(className="dt-right",targets=1:5))),
        class="compact hover stripe") |>
        DT::formatRound("Raw FST",6)|>DT::formatRound("FST-ENA",6)|>
        DT::formatRound("\u0394FST (ENA\u2212raw)",6)|>
        DT::formatStyle("FST-ENA",backgroundColor=DT::styleInterval(c(0.05,0.15,0.25),c("#f0fdf4","#dcfce7","#fefce8","#fef2f2")))|>
        DT::formatStyle("Locus",fontWeight="600",color="#0f172a")
    }, server=TRUE)

    output$ui_boot_fst_global <- renderUI({
      r <- boot_fst_global_r()
      shiny::validate(shiny::need(!is.null(r),"Run bootstrap first."))
      fst_obs_global <- tryCatch(fst_global_r()$global_ena, error=function(e) NA_real_)
      fst_raw_global <- tryCatch(fst_global_r()$global_raw, error=function(e) NA_real_)
      # Per-locus table using seq_len to preserve Locus names
      d    <- r$per_locus
      rows <- paste(sapply(seq_len(nrow(d)), function(i)
        sprintf('<tr><td class="locus-label">%s</td><td>%s</td><td>%s</td></tr>',
                htmltools::htmlEscape(as.character(d$Locus[i])),
                d$FST_ENA_obs[i], d$JK_pseudoval[i])),
        collapse="")
      tags$div(class="na-boot-result",
        tags$strong("Bootstrap over loci \u2014 Global FST-ENA (95% CI, percentile method)"),
        tags$br(),
        sprintf("FST-ENA obs. : %.6f", fst_obs_global), tags$br(),
        sprintf("95%% CI       : [ %.6f  \u2013  %.6f ]   (median %.6f)", r$ena[1], r$ena[3], r$ena[2]),
        tags$br(), tags$br(),
        tags$strong("Raw FST (for comparison):"), tags$br(),
        sprintf("FST raw obs. : %.6f", fst_raw_global), tags$br(),
        sprintf("95%% CI       : [ %.6f  \u2013  %.6f ]   (median %.6f)", r$raw[1], r$raw[3], r$raw[2]),
        tags$br(), tags$br(),
        tags$strong("Per-locus jackknife pseudo-values (leave-one-out):"),
        tags$div(class="na-matrix-wrap",
          HTML(paste0(
            '<table class="na-matrix" style="width:auto">',
            '<thead><tr><th>Locus</th><th>FST-ENA obs.</th><th>JK pseudo-value</th></tr></thead><tbody>',
            rows, '</tbody></table>')))
      )
    })

    # ── Matrix renderer ────────────────────────────────────────────────────────
    render_matrix_html <- function(mat, fmt=6,
                                   thr =c(0.05,0.15,0.25),
                                   clrs=c("#f0fdf4","#dcfce7","#fefce8","#fef2f2")) {
      pops <- rownames(mat); n <- length(pops)
      cell <- function(i,j) {
        if (i==j) return('<td class="diag">\u2014</td>')
        if (i<j)  return('<td class="upper">\u00b7</td>')
        v <- mat[i,j]
        if (is.na(v)) return('<td style="color:#94a3b8;">NA</td>')
        bg <- clrs[findInterval(v,thr)+1L]
        sprintf('<td style="background:%s;">%s</td>',bg,round(v,fmt))
      }
      thead <- paste0('<tr><th></th>',paste(sprintf('<th>%s</th>',pops[-n]),collapse=""),'</tr>')
      tbody <- paste(sapply(seq_len(n),function(i){
        if(i==1L) return("")
        paste0('<tr><td class="pop-label">',pops[i],'</td>',
               paste(sapply(seq_len(n),function(j)cell(i,j)),collapse=""),'</tr>')
      }),collapse="")
      HTML(sprintf('<div class="na-matrix-wrap"><table class="na-matrix"><thead>%s</thead><tbody>%s</tbody></table></div>',
                   thead,tbody))
    }

    # ── FIX: boot_tbl — use seq_len(nrow(d)) to preserve column types ─────────
    # "Locus","Pop1","Pop2" treated as character labels (not coerced to numeric)
    # numeric columns formatted with formatC
    boot_tbl <- function(d, cols, col_labels, char_cols = c("Locus","Pop1","Pop2")) {
      rows_html <- sapply(seq_len(nrow(d)), function(i) {
        cells <- paste(sapply(cols, function(cn) {
          val <- d[[cn]][i]
          if (cn %in% char_cols) {
            sprintf('<td class="locus-label">%s</td>', htmltools::htmlEscape(as.character(val)))
          } else {
            num <- suppressWarnings(as.numeric(val))
            cell_txt <- if (is.na(num)) "NA" else formatC(num, digits=6, format="f")
            sprintf('<td>%s</td>', cell_txt)
          }
        }), collapse="")
        paste0("<tr>", cells, "</tr>")
      })
      HTML(paste0(
        '<table class="na-matrix" style="width:100%"><thead><tr>',
        paste(sprintf("<th>%s</th>", col_labels), collapse=""),
        '</tr></thead><tbody>',
        paste(rows_html, collapse=""),
        '</tbody></table>'))
    }

    # ── Tab 4 pairwise FST ─────────────────────────────────────────────────────
    output$ui_fst_pair_matrix <- renderUI({
      r <- fst_pair_r(); typ <- input$fst_pair_type
      shiny::validate(shiny::need(!is.null(r$matrix_raw),"Click Compute."))
      if (identical(typ,"both")) tags$div(
        tags$p(tags$strong("Raw FST")), render_matrix_html(r$matrix_raw), tags$br(),
        tags$p(tags$strong("FST-ENA")), render_matrix_html(r$matrix_ena))
      else if (identical(typ,"raw")) render_matrix_html(r$matrix_raw)
      else render_matrix_html(r$matrix_ena)
    })

    output$dt_fst_pair <- DT::renderDT({
      r <- fst_pair_r(); d <- r$long
      shiny::validate(shiny::need(nrow(d)>0,"No data. Click Compute."))
      names(d) <- c("Pop 1","Pop 2","Raw FST","FST-ENA","\u0394FST (ENA\u2212raw)")
      DT::datatable(d,rownames=FALSE,
        options=list(pageLength=20,scrollX=TRUE,dom="lftip",columnDefs=list(list(className="dt-right",targets=2:4))),
        class="compact hover stripe") |>
        DT::formatRound("Raw FST",6)|>DT::formatRound("FST-ENA",6)|>
        DT::formatRound("\u0394FST (ENA\u2212raw)",6)|>
        DT::formatStyle("FST-ENA",backgroundColor=DT::styleInterval(c(0.05,0.15,0.25),c("#f0fdf4","#dcfce7","#fefce8","#fef2f2")))
    }, server=TRUE)

    output$ui_boot_fst_pair <- renderUI({
      r <- boot_fst_pair_r()
      shiny::validate(shiny::need(!is.null(r)&&length(r)>0,"Run bootstrap first."))
      parts <- list()
      if (!is.null(r$loci)) {
        parts$summ <- tagList(
          tags$p(tags$strong("\u25b6 Bootstrap over loci \u2014 pairwise FST-ENA (95% CI)")),
          tags$div(class="na-matrix-wrap",
            boot_tbl(r$loci$summary,
              c("Pop1","Pop2","FST_ENA_obs","CI_lo_loci","Median_loci","CI_hi_loci","FST_raw_obs","CI_lo_raw","CI_hi_raw"),
              c("Pop 1","Pop 2","FST-ENA obs.","CI lo (loci)","Median (loci)","CI hi (loci)","Raw FST obs.","CI lo raw","CI hi raw"))),
          tags$br(),
          tags$p(tags$strong("Per-locus observed FST-ENA and Raw FST")),
          tags$div(class="na-matrix-wrap",
            boot_tbl(r$loci$per_locus,
              c("Locus","Pop1","Pop2","FST_ENA_locus","FST_raw_locus"),
              c("Locus","Pop 1","Pop 2","FST-ENA locus","Raw FST locus")))
        )
      }
      if (!is.null(r$indiv)) {
        parts$indiv <- tagList(tags$br(),
          tags$p(tags$strong("\u25b6 Bootstrap over individuals \u2014 pairwise FST-ENA (95% CI)")),
          tags$div(class="na-matrix-wrap",
            boot_tbl(r$indiv,
              c("Pop1","Pop2","FST_ENA_obs","CI_lo_indiv","Median_indiv","CI_hi_indiv"),
              c("Pop 1","Pop 2","FST-ENA obs.","CI lo (indiv)","Median (indiv)","CI hi (indiv)")))
        )
      }
      tags$div(class="na-boot-result", do.call(tagList, parts))
    })

    # ── Tab 5 pairwise DCSE ────────────────────────────────────────────────────
    output$ui_dc_matrix <- renderUI({
      r <- dc_r(); typ <- input$dc_type
      shiny::validate(shiny::need(!is.null(r$matrix_raw),"Click Compute."))
      thr <- c(0.1,0.25,0.4); clrs <- c("#eff6ff","#dbeafe","#fef9c3","#fef2f2")
      if (identical(typ,"both")) tags$div(
        tags$p(tags$strong("Raw DCSE")),  render_matrix_html(r$matrix_raw,thr=thr,clrs=clrs), tags$br(),
        tags$p(tags$strong("DCSE-INA")), render_matrix_html(r$matrix_ina,thr=thr,clrs=clrs))
      else if (identical(typ,"raw")) render_matrix_html(r$matrix_raw,thr=thr,clrs=clrs)
      else render_matrix_html(r$matrix_ina,thr=thr,clrs=clrs)
    })

    output$dt_dc <- DT::renderDT({
      r <- dc_r(); d <- r$long
      shiny::validate(shiny::need(nrow(d)>0,"No data. Click Compute."))
      names(d) <- c("Pop 1","Pop 2","Raw DCSE","DCSE-INA","\u0394DCSE (INA\u2212raw)")
      DT::datatable(d,rownames=FALSE,
        options=list(pageLength=20,scrollX=TRUE,dom="lftip",columnDefs=list(list(className="dt-right",targets=2:4))),
        class="compact hover stripe") |>
        DT::formatRound("Raw DCSE",6)|>DT::formatRound("DCSE-INA",6)|>
        DT::formatRound("\u0394DCSE (INA\u2212raw)",6)
    }, server=TRUE)

    output$ui_boot_dc <- renderUI({
      r <- boot_dc_r()
      shiny::validate(shiny::need(!is.null(r)&&length(r)>0,"Run bootstrap first."))
      parts <- list()
      if (!is.null(r$loci)) {
        parts$summ <- tagList(
          tags$p(tags$strong("\u25b6 Bootstrap over loci \u2014 pairwise DCSE-INA (95% CI)")),
          tags$div(class="na-matrix-wrap",
            boot_tbl(r$loci$summary,
              c("Pop1","Pop2","DCSE_INA_obs","CI_lo_loci","Median_loci","CI_hi_loci","DCSE_raw_obs","CI_lo_raw","CI_hi_raw"),
              c("Pop 1","Pop 2","DCSE-INA obs.","CI lo (loci)","Median (loci)","CI hi (loci)","Raw DCSE obs.","CI lo raw","CI hi raw"))),
          tags$br(),
          tags$p(tags$strong("Per-locus observed DCSE-INA and Raw DCSE")),
          tags$div(class="na-matrix-wrap",
            boot_tbl(r$loci$per_locus,
              c("Locus","Pop1","Pop2","DCSE_INA_locus","DCSE_raw_locus"),
              c("Locus","Pop 1","Pop 2","DCSE-INA locus","Raw DCSE locus")))
        )
      }
      if (!is.null(r$indiv)) {
        parts$indiv <- tagList(tags$br(),
          tags$p(tags$strong("\u25b6 Bootstrap over individuals \u2014 pairwise DCSE-INA (95% CI)")),
          tags$div(class="na-matrix-wrap",
            boot_tbl(r$indiv,
              c("Pop1","Pop2","DCSE_INA_obs","CI_lo_indiv","Median_indiv","CI_hi_indiv"),
              c("Pop 1","Pop 2","DCSE-INA obs.","CI lo (indiv)","Median (indiv)","CI hi (indiv)")))
        )
      }
      tags$div(class="na-boot-result", do.call(tagList, parts))
    })

    # ── Tab 6 per-locus FST ────────────────────────────────────────────────────
    output$dt_fst_locus <- DT::renderDT({
      d <- fst_locus_r()
      shiny::validate(shiny::need(nrow(d)>0,"No data. Click Compute."))
      names(d) <- c("Locus","Pop 1","Pop 2","Raw FST","FST-ENA","\u0394FST","N_i raw","N_j raw","N_i ENA","N_j ENA")
      DT::datatable(d,rownames=FALSE,
        options=list(pageLength=25,scrollX=TRUE,dom="lftip",columnDefs=list(list(className="dt-right",targets=3:9))),
        class="compact hover stripe") |>
        DT::formatRound("Raw FST",6)|>DT::formatRound("FST-ENA",6)|>DT::formatRound("\u0394FST",6)|>
        DT::formatStyle("FST-ENA",backgroundColor=DT::styleInterval(c(0.05,0.15,0.25),c("#f0fdf4","#dcfce7","#fefce8","#fef2f2")))|>
        DT::formatStyle("Locus",fontWeight="600",color="#0f172a")
    }, server=TRUE)

    # ── Download handlers ──────────────────────────────────────────────────────
    make_dl <- function(data_fn, base_name, col_nms=NULL) {
      mk <- function(ext, write_fn)
        downloadHandler(
          filename=function() paste0(base_name,"_",Sys.Date(),".",ext),
          content=function(file) {
            d <- data_fn(); if (is.null(d)||nrow(d)==0L) return(invisible(NULL))
            if (!is.null(col_nms)) names(d) <- col_nms
            write_fn(d,file)
          })
      list(csv=mk("csv",function(d,f) write.csv(d,f,row.names=FALSE)),
           txt=mk("txt",function(d,f) write.table(d,f,sep="\t",row.names=FALSE,quote=FALSE)))
    }

    dl1 <- make_dl(function()t1_data_r(),"null_allele_per_pop_locus",
      c("Locus_names","Farm","p_nulls","N","N_exp_blanks","p_nulls_x_N"))
    output$dl_t1_csv <- dl1$csv; output$dl_t1_txt <- dl1$txt

    dl2 <- make_dl(function()t2_data_r(),"null_allele_global",
      c("Locus_names","Av_N_exp_blanks","Av_p_nulls","N_tot","N_blanks","f_expBlanks","p_nulls"))
    output$dl_t2_csv <- dl2$csv; output$dl_t2_txt <- dl2$txt

    dl3 <- make_dl(function()fst_global_r()$per_locus,"fst_global_ENA",
      c("Locus","FST_raw","FST_ENA","Delta_FST","N_pops_eff_raw","N_pops_eff_ENA"))
    output$dl_fst_global_csv <- dl3$csv; output$dl_fst_global_txt <- dl3$txt

    output$dl_fst_pair_csv <- downloadHandler(
      filename=function() paste0("fst_pairwise_ENA_",Sys.Date(),".csv"),
      content=function(file) { r<-fst_pair_r(); mat<-round(r$matrix_ena,6)
        write.csv(cbind(Population=rownames(mat),as.data.frame(mat)),file,row.names=FALSE) })
    output$dl_fst_pair_txt <- downloadHandler(
      filename=function() paste0("fst_pairwise_ENA_",Sys.Date(),".txt"),
      content=function(file) { r<-fst_pair_r(); mat<-round(r$matrix_ena,6)
        write.table(cbind(Population=rownames(mat),as.data.frame(mat)),file,sep="\t",row.names=FALSE,quote=FALSE) })
    dl4l <- make_dl(function()fst_pair_r()$long,"fst_pairwise_long_ENA",c("Pop1","Pop2","FST_raw","FST_ENA","Delta_FST"))
    output$dl_fst_pair_long_csv <- dl4l$csv; output$dl_fst_pair_long_txt <- dl4l$txt

    output$dl_dc_csv <- downloadHandler(
      filename=function() paste0("dcse_pairwise_INA_",Sys.Date(),".csv"),
      content=function(file) { r<-dc_r(); mat<-round(r$matrix_ina,6)
        write.csv(cbind(Population=rownames(mat),as.data.frame(mat)),file,row.names=FALSE) })
    output$dl_dc_txt <- downloadHandler(
      filename=function() paste0("dcse_pairwise_INA_",Sys.Date(),".txt"),
      content=function(file) { r<-dc_r(); mat<-round(r$matrix_ina,6)
        write.table(cbind(Population=rownames(mat),as.data.frame(mat)),file,sep="\t",row.names=FALSE,quote=FALSE) })
    dl5l <- make_dl(function()dc_r()$long,"dcse_pairwise_long_INA",c("Pop1","Pop2","DCSE_raw","DCSE_INA","Delta_DCSE"))
    output$dl_dc_long_csv <- dl5l$csv; output$dl_dc_long_txt <- dl5l$txt

    dl6 <- make_dl(function()fst_locus_r(),"fst_per_locus_pair_ENA",
      c("Locus","Pop1","Pop2","FST_raw","FST_ENA","Delta_FST","N_i_raw","N_j_raw","N_i_ENA","N_j_ENA"))
    output$dl_fst_locus_csv <- dl6$csv; output$dl_fst_locus_txt <- dl6$txt

    # Bootstrap downloads
    dl_b3 <- make_dl(function() {
      r <- boot_fst_global_r()
      rbind(data.frame(Type="FST_ENA",CI_lo=r$ena[1],Median=r$ena[2],CI_hi=r$ena[3]),
            data.frame(Type="FST_raw",CI_lo=r$raw[1],Median=r$raw[2],CI_hi=r$raw[3]))
    }, "bootstrap_fst_global_ENA")
    output$dl_boot_fst_global_csv <- dl_b3$csv; output$dl_boot_fst_global_txt <- dl_b3$txt

    dl_b4 <- make_dl(function() {
      r <- boot_fst_pair_r()
      parts <- list()
      if (!is.null(r$loci))  parts <- c(parts, list(r$loci$summary), list(r$loci$per_locus))
      if (!is.null(r$indiv)) parts <- c(parts, list(r$indiv))
      tryCatch(Reduce(function(a,b) merge(a,b,all=TRUE), parts), error=function(e) parts[[1]])
    }, "bootstrap_fst_pair_ENA")
    output$dl_boot_fst_pair_csv <- dl_b4$csv; output$dl_boot_fst_pair_txt <- dl_b4$txt

    dl_b5 <- make_dl(function() {
      r <- boot_dc_r()
      parts <- list()
      if (!is.null(r$loci))  parts <- c(parts, list(r$loci$summary), list(r$loci$per_locus))
      if (!is.null(r$indiv)) parts <- c(parts, list(r$indiv))
      tryCatch(Reduce(function(a,b) merge(a,b,all=TRUE), parts), error=function(e) parts[[1]])
    }, "bootstrap_dcse_pair_INA")
    output$dl_boot_dc_csv <- dl_b5$csv; output$dl_boot_dc_txt <- dl_b5$txt

  }) # end moduleServer
}