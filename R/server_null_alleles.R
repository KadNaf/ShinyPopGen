# module/server_null_alleles.R
# Null allele frequency estimation (EM), FST-ENA, DCSE-INA

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
    treat_id <- function(loc) paste0("coding_", gsub("[^A-Za-z0-9]", "_", loc))

    ci_bounds <- function(alpha) {
      lo <- alpha / 2
      hi <- 1 - alpha / 2
      list(lo = lo, hi = hi,
           label = paste0(round((1 - alpha) * 100, 3), "% CI"))
    }

    fmt6 <- function(x) if (is.na(x)) "NA" else formatC(as.numeric(x), digits = 6, format = "f")

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
      updateSelectInput(session,"fl_locus",choices=c("All loci"="all",stats::setNames(markers,markers)),selected="all")
      updateSelectInput(session,"fl_pop1", choices=c("All pairs"="all",stats::setNames(pops,pops)),selected="all")
      updateSelectInput(session,"fl_pop2", choices=c("All pairs"="all",stats::setNames(pops,pops)),selected="all")
    })

    # ── Per-locus coding UI — radio buttons, default = absent (000000) ─────────
    output$locus_coding_ui <- renderUI({
      ns_fn   <- session$ns
      markers <- markers_r()
      if (length(markers) == 0L) return(tags$p("No markers loaded yet."))
      items <- lapply(markers, function(loc) {
        tags$div(class = "na-locus-item",
          tags$div(class = "na-locus-name", loc),
          radioButtons(
            inputId  = ns_fn(treat_id(loc)),
            label    = NULL,
            choices  = c(
              "000000 — absent / PCR failure" = "absent",
              "999999 — null homozygote"      = "null_homo"
            ),
            selected = "absent",
            inline   = TRUE
          )
        )
      })
      tags$div(class = "na-locus-grid", items)
    })

    locus_treatments_r <- reactive({
      markers <- markers_r()
      treats  <- sapply(markers, function(loc) {
        val <- input[[treat_id(loc)]]
        if (is.null(val) || !val %in% c("null_homo","absent")) "absent" else val
      })
      stats::setNames(treats, markers)
    })

    # ══════════════════════════════════════════════════════════════════════════
    #  EM ALGORITHM
    # ══════════════════════════════════════════════════════════════════════════
    em_freena <- function(gt_vec, base, treat = "absent") {
      efpop      <- length(gt_vec)
      absent_msk <- is.na(gt_vec) | gt_vec <= 0L
      n_absent   <- sum(absent_msk)
      valid_gt   <- gt_vec[!absent_msk]

      empty <- list(rd=0.0, pfreq=numeric(0), genefreq_obs=numeric(0),
                    H_ii=numeric(0), H_iX=numeric(0), N=0L, efpop=efpop,
                    n_absent=n_absent, n_null_homo=0L, alleles=integer(0),
                    n_valid_geno=0L)

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
      if (n_valid_geno == 0L) {
        empty$N <- N; empty$n_null_homo <- n_null_homo; return(empty)
      }

      genefreq_obs <- sapply(alleles, function(a)
        (sum(valid_a1==a) + sum(valid_a2==a)) / (2L * n_valid_geno))
      H_ii  <- sapply(alleles, function(a) sum(valid_a1==a & valid_a2==a))
      H_iX  <- sapply(alleles, function(a)
        sum((valid_a1==a & valid_a2!=a) | (valid_a2==a & valid_a1!=a)))
      hotot <- sum(H_ii)

      rd <- if (treat == "null_homo" && n_null_homo > 0L)
              sqrt(n_null_homo / N)
            else
              sqrt(1.0 / (N + 1.0))

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
        rd_new <- if (treat == "null_homo")
                    rdi + (2.0 * n_null_homo) / (2.0 * N)
                  else
                    rdi
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

    # ── Weir (1996) FST components ─────────────────────────────────────────────
    weir_components_allele <- function(pop_list, use_corr = FALSE) {
      r     <- length(pop_list)
      N_tot <- sum(sapply(pop_list, `[[`, "ni"))
      N2    <- sum(sapply(pop_list, function(p) p$ni^2))
      if (N_tot == 0L || r < 2L) return(list(s2P=0.0, s2I=0.0, s2G=0.0))
      nc <- (N_tot - N2/N_tot) / (r - 1)
      if (nc <= 0 || N_tot - r <= 0) return(list(s2P=0.0, s2I=0.0, s2G=0.0))
      snA  <- sum(sapply(pop_list, `[[`, "nA"))
      s2A  <- sum(sapply(pop_list, function(p) if(p$ni>0) p$nA^2/(2*p$ni) else 0.0))
      sAA  <- if (use_corr) sum(sapply(pop_list,`[[`,"AA_corr"))
              else           sum(sapply(pop_list,`[[`,"AA"))
      MSG  <- (0.5*snA - sAA) / N_tot
      MSI  <- (0.5*snA + sAA - s2A) / (N_tot - r)
      MSP  <- (s2A - 0.5*snA^2/N_tot) / (r - 1)
      list(s2P=(MSP-MSI)/(2*nc), s2I=0.5*(MSI-MSG), s2G=MSG)
    }

    # ── CS distance ────────────────────────────────────────────────────────────
    cs_distance <- function(freq_i, freq_j) {
      alleles <- union(names(freq_i), names(freq_j))
      csprod  <- 0.0
      for (a in alleles) {
        fi <- freq_i[a]; fj <- freq_j[a]
        fi <- if (is.null(fi)||is.na(fi)) 0.0 else as.numeric(fi)
        fj <- if (is.null(fj)||is.na(fj)) 0.0 else as.numeric(fj)
        if (fi > 0 && fj > 0) csprod <- csprod + sqrt(fi * fj)
      }
      if (csprod > 1.0) return(NA_real_)
      (2.0 / base::pi) * sqrt(2.0 * (1.0 - csprod))
    }

    make_ina_freq <- function(em) c(em$pfreq, `__null__` = em$rd)

    # ══════════════════════════════════════════════════════════════════════════
    #  FETCH ALL GENOTYPES
    # ══════════════════════════════════════════════════════════════════════════
    raw_data_r <- reactive({
      db_ready()
      con   <- con_r(); hs <- hf_schema_r(); ms <- meta_schema_r()
      hf_q  <- sql_id(con,tbl_hf_r());  meta_q <- sql_id(con,tbl_meta_r())
      hi_q  <- sql_id(con,hs$ind_col);  hl_q   <- sql_id(con,hs$locus_col)
      hg_q  <- sql_id(con,hs$gt_col);   mi_q   <- sql_id(con,ms$ind_col)
      pop_q <- sql_id(con,ms$pop_col)
      DBI::dbGetQuery(con, sprintf("
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
        hf_q, meta_q, hi_q, mi_q, hl_q, pop_q))
    })

    # ══════════════════════════════════════════════════════════════════════════
    #  MAIN COMPUTATIONS
    # ══════════════════════════════════════════════════════════════════════════
    
    compute_fst_global_full <- function(em_res) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]])
      s1_vec <- s3_vec <- s1c_vec <- s3c_vec <- numeric(length(markers))
      rows <- vector("list", length(markers))
      for (li in seq_along(markers)) {
        loc <- markers[li]; em_loc <- em_res[[loc]]
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
        s1_vec[li]  <- s1l *nc_raw;  s3_vec[li]  <- s3l *nc_raw
        s1c_vec[li] <- s1lc*nc_corr; s3c_vec[li] <- s3lc*nc_corr
        rows[[li]] <- data.frame(Locus=loc,
          FST_raw=round(fst_loc,6), FST_ENA=round(fst_locc,6),
          Delta_FST=round(fst_locc-fst_loc,6),
          N_pops_raw=r_raw, N_pops_ENA=r_corr, stringsAsFactors=FALSE)
      }
      s1 <- sum(s1_vec); s3 <- sum(s3_vec)
      s1c <- sum(s1c_vec); s3c <- sum(s3c_vec)
      list(global_raw=if(s3>0)s1/s3 else NA_real_,
           global_ena=if(s3c>0)s1c/s3c else NA_real_,
           per_locus=do.call(rbind,rows),
           s1_vec=s1_vec, s3_vec=s3_vec,
           s1c_vec=s1c_vec, s3c_vec=s3c_vec,
           markers=markers)
    }

    compute_fst_pairwise <- function(em_res) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]]); n_pops <- length(pops)
      if (n_pops < 2L) return(list(matrix_raw=NULL,matrix_ena=NULL,long=data.frame(),
                                   s1_raw=NULL,s3_raw=NULL,s1_ena=NULL,s3_ena=NULL,pairs=NULL))
      pairs   <- combn(pops,2,simplify=FALSE); n_pairs <- length(pairs)
      s12p  <- matrix(0.0,n_pops,n_pops,dimnames=list(pops,pops))
      s32p  <- matrix(0.0,n_pops,n_pops,dimnames=list(pops,pops))
      s12pc <- matrix(0.0,n_pops,n_pops,dimnames=list(pops,pops))
      s32pc <- matrix(0.0,n_pops,n_pops,dimnames=list(pops,pops))
      s1_raw_ml <- matrix(0.0,n_pairs,length(markers))
      s3_raw_ml <- matrix(0.0,n_pairs,length(markers))
      s1_ena_ml <- matrix(0.0,n_pairs,length(markers))
      s3_ena_ml <- matrix(0.0,n_pairs,length(markers))

      for (li in seq_along(markers)) {
        loc <- markers[li]; em_loc <- em_res[[loc]]
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
            pi_idx <- which(sapply(pairs, function(p) p[1]==pi_n && p[2]==pj_n))
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
                if (length(pi_idx)==1) {
                  s1_raw_ml[pi_idx,li] <- s1_raw_ml[pi_idx,li]+cmp$s2P*nc_r
                  s3_raw_ml[pi_idx,li] <- s3_raw_ml[pi_idx,li]+(cmp$s2P+cmp$s2I+cmp$s2G)*nc_r
                }
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
                if (length(pi_idx)==1) {
                  s1_ena_ml[pi_idx,li] <- s1_ena_ml[pi_idx,li]+cmpc$s2P*nc_c
                  s3_ena_ml[pi_idx,li] <- s3_ena_ml[pi_idx,li]+(cmpc$s2P+cmpc$s2I+cmpc$s2G)*nc_c
                }
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
      list(matrix_raw=mat_raw, matrix_ena=mat_ena, long=do.call(rbind,long_rows),
           s1_raw=s1_raw_ml, s3_raw=s3_raw_ml,
           s1_ena=s1_ena_ml, s3_ena=s3_ena_ml,
           pairs=pairs, markers=markers)
    }

    compute_dc_pairwise <- function(em_res) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]]); n_pops <- length(pops)
      if (n_pops < 2L) return(list(matrix_raw=NULL,matrix_ina=NULL,long=data.frame(),
                                   dc_raw=NULL,dc_ina=NULL,pairs=NULL))
      pairs   <- combn(pops,2,simplify=FALSE); n_pairs <- length(pairs)
      dc_sum_raw <- matrix(0.0,n_pops,n_pops,dimnames=list(pops,pops))
      dc_sum_ina <- matrix(0.0,n_pops,n_pops,dimnames=list(pops,pops))
      nloc_eff   <- matrix(length(markers),n_pops,n_pops,dimnames=list(pops,pops))
      nloc_eff_c <- matrix(length(markers),n_pops,n_pops,dimnames=list(pops,pops))
      dc_raw_ml  <- matrix(NA_real_,n_pairs,length(markers))
      dc_ina_ml  <- matrix(NA_real_,n_pairs,length(markers))
      for (li in seq_along(markers)) {
        loc <- markers[li]; em_loc <- em_res[[loc]]
        for (ii in seq_len(n_pops-1L)) {
          for (jj in seq(ii+1L,n_pops)) {
            ei <- em_loc[[pops[ii]]]; ej <- em_loc[[pops[jj]]]
            ni_ri <- ei$efpop-ei$n_absent-ei$n_null_homo
            ni_rj <- ej$efpop-ej$n_absent-ej$n_null_homo
            ni_ci <- ei$efpop-ei$n_absent; ni_cj <- ej$efpop-ej$n_absent
            pi_idx <- which(sapply(pairs, function(p) p[1]==pops[ii] && p[2]==pops[jj]))
            if (ni_ri>0L&&ni_rj>0L&&length(ei$genefreq_obs)>0&&length(ej$genefreq_obs)>0) {
              d_raw <- cs_distance(ei$genefreq_obs, ej$genefreq_obs)
              if (!is.na(d_raw)) {
                dc_sum_raw[jj,ii] <- dc_sum_raw[jj,ii]+d_raw
                if (length(pi_idx)==1) dc_raw_ml[pi_idx,li] <- d_raw
              } else nloc_eff[jj,ii] <- nloc_eff[jj,ii]-1L
            } else nloc_eff[jj,ii] <- nloc_eff[jj,ii]-1L
            if (ni_ci>0L&&ni_cj>0L) {
              d_ina <- cs_distance(make_ina_freq(ei), make_ina_freq(ej))
              if (!is.na(d_ina)) {
                dc_sum_ina[jj,ii] <- dc_sum_ina[jj,ii]+d_ina
                if (length(pi_idx)==1) dc_ina_ml[pi_idx,li] <- d_ina
              } else nloc_eff_c[jj,ii] <- nloc_eff_c[jj,ii]-1L
            } else nloc_eff_c[jj,ii] <- nloc_eff_c[jj,ii]-1L
          }
        }
      }
      mat_raw <- matrix(NA_real_,n_pops,n_pops,dimnames=list(pops,pops))
      mat_ina <- matrix(NA_real_,n_pops,n_pops,dimnames=list(pops,pops))
      for (ii in seq_len(n_pops-1L)) for (jj in seq(ii+1L,n_pops)) {
        mat_raw[jj,ii] <- if(nloc_eff[jj,ii]  >0L) dc_sum_raw[jj,ii]/nloc_eff[jj,ii]   else NA_real_
        mat_ina[jj,ii] <- if(nloc_eff_c[jj,ii]>0L) dc_sum_ina[jj,ii]/nloc_eff_c[jj,ii] else NA_real_
      }
      long_rows <- list()
      for (ii in seq_len(n_pops-1L)) for (jj in seq(ii+1L,n_pops))
        long_rows[[length(long_rows)+1L]] <- data.frame(
          Pop1=pops[ii], Pop2=pops[jj],
          DCSE_raw=round(mat_raw[jj,ii],6), DCSE_INA=round(mat_ina[jj,ii],6),
          Delta_DCSE=round(mat_ina[jj,ii]-mat_raw[jj,ii],6), stringsAsFactors=FALSE)
      list(matrix_raw=mat_raw, matrix_ina=mat_ina, long=do.call(rbind,long_rows),
           dc_raw=dc_raw_ml, dc_ina=dc_ina_ml, pairs=pairs, markers=markers)
    }

    compute_per_locus_pair <- function(em_res) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]])
      rows_fst <- rows_dc <- list()
      for (loc in markers) {
        em_loc <- em_res[[loc]]
        alleles_obs <- sort(unique(unlist(lapply(em_loc, function(e) e$alleles))))
        for (ii in seq_len(length(pops)-1L)) {
          for (jj in seq(ii+1L,length(pops))) {
            pi_n <- pops[ii]; pj_n <- pops[jj]
            ei <- em_loc[[pi_n]]; ej <- em_loc[[pj_n]]
            ni_ri <- max(0L,ei$efpop-ei$n_absent-ei$n_null_homo)
            ni_rj <- max(0L,ej$efpop-ej$n_absent-ej$n_null_homo)
            ni_ci <- max(0L,ei$efpop-ei$n_absent)
            ni_cj <- max(0L,ej$efpop-ej$n_absent)
            N_r <- ni_ri+ni_rj; N_c <- ni_ci+ni_cj
            nc_r <- if(N_r>0&&ni_ri>0&&ni_rj>0) (N_r-(ni_ri^2+ni_rj^2)/N_r) else 0.0
            nc_c <- if(N_c>0&&ni_ci>0&&ni_cj>0) (N_c-(ni_ci^2+ni_cj^2)/N_c) else 0.0
            s1_r<-s3_r<-s1_c<-s3_c<-0.0
            for (a in alleles_obs) {
              a_chr <- as.character(a)
              if (nc_r>0) {
                pd <- list(
                  list(ni=ni_ri,nA=(if(a_chr %in% names(ei$genefreq_obs))ei$genefreq_obs[[a_chr]] else 0.0)*2L*ni_ri,
                       AA=if(a_chr %in% names(ei$H_ii))ei$H_ii[[a_chr]] else 0L,AA_corr=0.0),
                  list(ni=ni_rj,nA=(if(a_chr %in% names(ej$genefreq_obs))ej$genefreq_obs[[a_chr]] else 0.0)*2L*ni_rj,
                       AA=if(a_chr %in% names(ej$H_ii))ej$H_ii[[a_chr]] else 0L,AA_corr=0.0))
                cmp <- weir_components_allele(pd,use_corr=FALSE)
                s1_r<-s1_r+cmp$s2P*nc_r; s3_r<-s3_r+(cmp$s2P+cmp$s2I+cmp$s2G)*nc_r
              }
              if (nc_c>0) {
                pf_i<-if(a_chr %in% names(ei$pfreq))ei$pfreq[[a_chr]] else 0.0
                pf_j<-if(a_chr %in% names(ej$pfreq))ej$pfreq[[a_chr]] else 0.0
                AA_i<-if(a_chr %in% names(ei$H_ii))ei$H_ii[[a_chr]] else 0L
                AA_j<-if(a_chr %in% names(ej$H_ii))ej$H_ii[[a_chr]] else 0L
                di<-pf_i+2.0*ei$rd; dj<-pf_j+2.0*ej$rd
                pdc <- list(
                  list(ni=ni_ci,nA=pf_i*2L*ni_ci,AA=AA_i,AA_corr=if(AA_i>0&&di>0)AA_i*(pf_i/di) else 0.0),
                  list(ni=ni_cj,nA=pf_j*2L*ni_cj,AA=AA_j,AA_corr=if(AA_j>0&&dj>0)AA_j*(pf_j/dj) else 0.0))
                cmpc <- weir_components_allele(pdc,use_corr=TRUE)
                s1_c<-s1_c+cmpc$s2P*nc_c; s3_c<-s3_c+(cmpc$s2P+cmpc$s2I+cmpc$s2G)*nc_c
              }
            }
            rows_fst[[length(rows_fst)+1L]] <- data.frame(
              Locus=loc, Pop1=pi_n, Pop2=pj_n,
              FST_raw=round(if(s3_r!=0)s1_r/s3_r else NA_real_,6),
              FST_ENA=round(if(s3_c!=0)s1_c/s3_c else NA_real_,6),
              stringsAsFactors=FALSE)
            d_raw_l <- if(ni_ri>0&&ni_rj>0&&length(ei$genefreq_obs)>0&&length(ej$genefreq_obs)>0)
              cs_distance(ei$genefreq_obs,ej$genefreq_obs) else NA_real_
            d_ina_l <- if(ni_ci>0&&ni_cj>0)
              cs_distance(make_ina_freq(ei),make_ina_freq(ej)) else NA_real_
            rows_dc[[length(rows_dc)+1L]] <- data.frame(
              Locus=loc, Pop1=pi_n, Pop2=pj_n,
              DCSE_raw=round(d_raw_l,6), DCSE_INA=round(d_ina_l,6),
              stringsAsFactors=FALSE)
          }
        }
      }
      list(fst=do.call(rbind,rows_fst), dc=do.call(rbind,rows_dc))
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  BOOTSTRAP FUNCTIONS
    # ══════════════════════════════════════════════════════════════════════════
    
    boot_loci_global_fst <- function(fst_full, nboot, alpha) {
      L   <- length(fst_full$markers)
      idx <- matrix(sample.int(L, L*nboot, replace=TRUE), nrow=nboot)
      S1  <- matrix(fst_full$s1_vec[idx],  nrow=nboot)
      S3  <- matrix(fst_full$s3_vec[idx],  nrow=nboot)
      S1c <- matrix(fst_full$s1c_vec[idx], nrow=nboot)
      S3c <- matrix(fst_full$s3c_vec[idx], nrow=nboot)
      fst_raw_b <- ifelse(rowSums(S3) >0, rowSums(S1) /rowSums(S3),  NA_real_)
      fst_ena_b <- ifelse(rowSums(S3c)>0, rowSums(S1c)/rowSums(S3c), NA_real_)
      list(
        raw = quantile(fst_raw_b, c(alpha/2, 0.5, 1-alpha/2), na.rm=TRUE),
        ena = quantile(fst_ena_b, c(alpha/2, 0.5, 1-alpha/2), na.rm=TRUE),
        dist = fst_ena_b
      )
    }

    boot_loci_pair_fst <- function(pair_res, nboot, alpha) {
      L       <- length(pair_res$markers)
      n_pairs <- length(pair_res$pairs)
      idx     <- matrix(sample.int(L, L*nboot, replace=TRUE), nrow=nboot)
      results <- vector("list", n_pairs)
      for (pi in seq_len(n_pairs)) {
        s1r <- pair_res$s1_raw[pi,]; s3r <- pair_res$s3_raw[pi,]
        s1e <- pair_res$s1_ena[pi,]; s3e <- pair_res$s3_ena[pi,]
        RS1r <- rowSums(matrix(s1r[idx],nrow=nboot)); RS3r <- rowSums(matrix(s3r[idx],nrow=nboot))
        RS1e <- rowSums(matrix(s1e[idx],nrow=nboot)); RS3e <- rowSums(matrix(s3e[idx],nrow=nboot))
        br <- ifelse(RS3r>0,RS1r/RS3r,NA_real_); be <- ifelse(RS3e>0,RS1e/RS3e,NA_real_)
        results[[pi]] <- data.frame(
          Pop1         = pair_res$pairs[[pi]][1],
          Pop2         = pair_res$pairs[[pi]][2],
          FST_ENA_obs  = round(if(sum(s3e)>0)sum(s1e)/sum(s3e) else NA_real_,6),
          CI_lo_loci   = round(quantile(be,alpha/2,na.rm=TRUE),6),
          Median_loci  = round(quantile(be,0.5,na.rm=TRUE),6),
          CI_hi_loci   = round(quantile(be,1-alpha/2,na.rm=TRUE),6),
          FST_raw_obs  = round(if(sum(s3r)>0)sum(s1r)/sum(s3r) else NA_real_,6),
          CI_lo_raw    = round(quantile(br,alpha/2,na.rm=TRUE),6),
          CI_hi_raw    = round(quantile(br,1-alpha/2,na.rm=TRUE),6),
          stringsAsFactors=FALSE)
      }
      do.call(rbind, results)
    }

    boot_loci_pair_dc <- function(dc_res, nboot, alpha) {
      L       <- length(dc_res$markers)
      n_pairs <- length(dc_res$pairs)
      idx     <- matrix(sample.int(L, L*nboot, replace=TRUE), nrow=nboot)
      results <- vector("list", n_pairs)
      for (pi in seq_len(n_pairs)) {
        dr <- dc_res$dc_raw[pi,]; di <- dc_res$dc_ina[pi,]
        br <- rowMeans(matrix(dr[idx],nrow=nboot),na.rm=TRUE)
        bi <- rowMeans(matrix(di[idx],nrow=nboot),na.rm=TRUE)
        results[[pi]] <- data.frame(
          Pop1         = dc_res$pairs[[pi]][1],
          Pop2         = dc_res$pairs[[pi]][2],
          DCSE_INA_obs = round(mean(di,na.rm=TRUE),6),
          CI_lo_loci   = round(quantile(bi,alpha/2,na.rm=TRUE),6),
          Median_loci  = round(quantile(bi,0.5,na.rm=TRUE),6),
          CI_hi_loci   = round(quantile(bi,1-alpha/2,na.rm=TRUE),6),
          DCSE_raw_obs = round(mean(dr,na.rm=TRUE),6),
          CI_lo_raw    = round(quantile(br,alpha/2,na.rm=TRUE),6),
          CI_hi_raw    = round(quantile(br,1-alpha/2,na.rm=TRUE),6),
          stringsAsFactors=FALSE)
      }
      do.call(rbind, results)
    }

    boot_subsamples_global_fst <- function(raw_df, em_res, base, treatments, nboot, alpha) {
      markers <- names(em_res); pops <- names(em_res[[markers[1]]])
      idx_by_pop_loc <- list()
      for (pop in pops) for (loc in markers)
        idx_by_pop_loc[[paste0(pop,"___",loc)]] <-
          which(raw_df$Population==pop & raw_df$Marker==loc)
      inds_by_pop <- lapply(pops, function(p) unique(raw_df$Individual[raw_df$Population==p]))
      names(inds_by_pop) <- pops
      boot_raw <- boot_ena <- numeric(nboot)
      for (b in seq_len(nboot)) {
        em_b <- list()
        for (loc in markers) {
          em_b[[loc]] <- list()
          treat <- as.character(treatments[loc] %||% "absent")
          for (pop in pops) {
            inds <- inds_by_pop[[pop]]
            if (length(inds)==0L) { em_b[[loc]][[pop]] <- em_res[[loc]][[pop]]; next }
            ri <- sample(inds, length(inds), replace=TRUE)
            br <- idx_by_pop_loc[[paste0(pop,"___",loc)]]
            ic <- raw_df$Individual[br]
            ir <- unlist(lapply(ri, function(x) br[ic==x]))
            gts <- raw_df$gt[ir]
            em_b[[loc]][[pop]] <- if(length(gts)==0L) em_res[[loc]][[pop]]
                                  else em_freena(gts, base, treat)
          }
        }
        fg <- compute_fst_global_full(em_b)
        boot_raw[b] <- fg$global_raw %||% NA_real_
        boot_ena[b] <- fg$global_ena %||% NA_real_
      }
      list(
        raw = quantile(boot_raw, c(alpha/2,0.5,1-alpha/2), na.rm=TRUE),
        ena = quantile(boot_ena, c(alpha/2,0.5,1-alpha/2), na.rm=TRUE)
      )
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  MAIN REACTIVE — single "run_all" button
    # ══════════════════════════════════════════════════════════════════════════
    results_r <- eventReactive(input$run_all, {
      req(db_ready())
      nboot  <- max(100L, min(99999L, as.integer(input$nboot %||% 5000L)))
      alpha  <- as.numeric(input$ci_level %||% "0.05")
      ci     <- ci_bounds(alpha)
      base   <- as.integer(base_r())
      treats <- locus_treatments_r()
      markers <- markers_r(); pops <- pops_r()

      withProgress(message = "Running computations...", value = 0, {

        setProgress(0.03, detail = "Fetching genotypes from database...")
        raw_df <- raw_data_r()
        shiny::validate(shiny::need(nrow(raw_df)>0, "No genotype data found."))

        setProgress(0.08, detail = "EM algorithm (null allele frequencies)...")
        em_res <- list()
        for (loc in markers) {
          em_res[[loc]] <- list()
          treat <- as.character(treats[loc] %||% "absent")
          for (pop in pops) {
            gts <- raw_df$gt[raw_df$Marker==loc & raw_df$Population==pop]
            em_res[[loc]][[pop]] <-
              if (length(gts)==0L)
                list(rd=0.0,pfreq=numeric(0),genefreq_obs=numeric(0),
                     H_ii=numeric(0),H_iX=numeric(0),N=0L,efpop=0L,
                     n_absent=0L,n_null_homo=0L,alleles=integer(0),n_valid_geno=0L)
              else em_freena(gts, base, treat)
          }
        }

        setProgress(0.15, detail = "Global FST and FST-ENA...")
        fst_global <- compute_fst_global_full(em_res)

        setProgress(0.25, detail = "Pairwise FST and FST-ENA...")
        fst_pair <- compute_fst_pairwise(em_res)

        setProgress(0.35, detail = "Pairwise DCSE and DCSE-INA...")
        dc_pair <- compute_dc_pairwise(em_res)

        setProgress(0.45, detail = "Per-locus x pair statistics...")
        per_locus_pair <- compute_per_locus_pair(em_res)

        setProgress(0.50, detail = sprintf("Bootstrap over loci — global FST (%d reps)...", nboot))
        boot_gl_loci <- boot_loci_global_fst(fst_global, nboot, alpha)

        setProgress(0.60, detail = sprintf("Bootstrap over loci — pairwise FST (%d reps)...", nboot))
        boot_pair_fst_loci <- boot_loci_pair_fst(fst_pair, nboot, alpha)

        setProgress(0.70, detail = sprintf("Bootstrap over loci — pairwise DCSE (%d reps)...", nboot))
        boot_pair_dc_loci <- boot_loci_pair_dc(dc_pair, nboot, alpha)

        setProgress(0.78, detail = sprintf("Bootstrap over sub-samples — global FST (%d reps)...", nboot))
        boot_gl_subs <- boot_subsamples_global_fst(raw_df, em_res, base, treats, nboot, alpha)

        setProgress(0.95, detail = "Assembling results...")

        t1_rows <- list()
        for (loc in markers) {
          for (pop in pops) {
            e <- em_res[[loc]][[pop]]
            n_exp <- e$N * (e$rd^2)
            t1_rows[[length(t1_rows)+1L]] <- data.frame(
              Locus=loc, Population=pop,
              Coding=as.character(treats[loc] %||% "absent"),
              p_nulls=round(e$rd,6), N=as.integer(e$N),
              N_exp_blanks=round(n_exp,6),
              stringsAsFactors=FALSE)
          }
        }
        t1 <- do.call(rbind, t1_rows)
        t1$Locus <- factor(t1$Locus, levels=markers)
        t1 <- t1[order(t1$Locus, t1$Population),]
        t1$Locus <- as.character(t1$Locus)

        t2_rows <- lapply(markers, function(loc) {
          sub  <- t1[t1$Locus==loc,,drop=FALSE]
          vidx <- !is.na(sub$p_nulls)
          av_p <- if (any(vidx)&&sum(sub$N[vidx])>0)
            sum(sub$p_nulls[vidx]*sub$N[vidx])/sum(sub$N[vidx]) else NA_real_
          av_n <- sum(sub$N*(sub$p_nulls^2), na.rm=TRUE)
          data.frame(Locus=loc, Coding=as.character(treats[loc] %||% "absent"),
                     Av_p_nulls=round(av_p,6), Av_N_exp=round(av_n,6),
                     N_tot=sum(sub$N), stringsAsFactors=FALSE)
        })
        t2 <- do.call(rbind, t2_rows)

        setProgress(1)

        list(
          t1 = t1, t2 = t2,
          fst_global    = fst_global,
          fst_pair      = fst_pair,
          dc_pair       = dc_pair,
          per_locus_pair = per_locus_pair,
          boot_gl_loci      = boot_gl_loci,
          boot_gl_subs      = boot_gl_subs,
          boot_pair_fst     = boot_pair_fst_loci,
          boot_pair_dc      = boot_pair_dc_loci,
          nboot = nboot, alpha = alpha, ci = ci,
          treats = treats, markers = markers, pops = pops,
          em_res = em_res
        )
      })
    })

    # ══════════════════════════════════════════════════════════════════════════
    #  METADATA HEADER for output files
    # ══════════════════════════════════════════════════════════════════════════
    meta_header <- function(r, file_desc) {
      ci_pct <- paste0(round((1 - r$alpha) * 100, 3), "%")
      treat_summary <- paste(sapply(r$markers, function(loc) {
        cd <- as.character(r$treats[loc] %||% "absent")
        sprintf("%s:%s", loc, if(cd=="absent") "000000" else "999999")
      }), collapse=", ")
      c(
        paste0("# ", file_desc),
        "# Method: Expectation-Maximization (EM) algorithm — Dempster, Laird & Rubin (1977)",
        "# ENA correction (Excluding Null Alleles) — Chapuis & Estoup (2007) / FreeNA",
        "# INA correction (Including Null Alleles) — Chapuis & Estoup (2007) / FreeNA",
        "# FST: Weir (1996) following Genepop method",
        "# DCSE: Cavalli-Sforza & Edwards (1967) chord genetic distance",
        paste0("# Bootstrap replicates: ", r$nboot),
        paste0("# Confidence interval: ", ci_pct, " (alpha = ", r$alpha, ")"),
        paste0("# Locus coding (000000=absent/PCR failure; 999999=null homozygote):"),
        paste0("#   ", treat_summary),
        "#"
      )
    }

    write_with_header <- function(hdr, df, file, sep = ",") {
      writeLines(hdr, con = file)
      write.table(df, file = file, sep = sep, row.names = FALSE,
                  quote = FALSE, append = TRUE,
                  col.names = TRUE)
    }

    half_matrix_txt <- function(df, stat_col, pops, loc) {
      sub  <- df[df$Locus == loc,,drop=FALSE]
      n    <- length(pops)
      lines <- character(0)
      lines <- c(lines, paste0("# Locus: ", loc, "  Statistic: ", stat_col))
      hdr <- paste(c("", pops[-n]), collapse="\t")
      lines <- c(lines, hdr)
      for (i in seq(2, n)) {
        row_vals <- sapply(seq_len(n), function(j) {
          if (j >= i) return("")
          row <- sub[sub$Pop1==pops[j] & sub$Pop2==pops[i],,drop=FALSE]
          if (nrow(row)==0) return("NA")
          v <- row[[stat_col]][1]
          if (is.na(v)) "NA" else as.character(round(v,6))
        })
        lines <- c(lines, paste(c(pops[i], row_vals[-n]), collapse="\t"))
      }
      lines
    }

    # ══════════════════════════════════════════════════════════════════════════
    #  FILE 1 — Null allele frequencies
    # ══════════════════════════════════════════════════════════════════════════
    file1_data <- reactive({
      r <- results_r()
      list(
        header = meta_header(r, "File 1 — Null allele frequencies per locus x population"),
        t1     = r$t1,
        t2     = r$t2
      )
    })

    output$dl_file1_csv <- downloadHandler(
      filename = function() paste0("null_allele_frequencies_", Sys.Date(), ".csv"),
      content  = function(file) {
        d <- file1_data()
        hdr <- c(d$header,
                 "# Section 1: p_nulls per locus x population",
                 "#")
        write_with_header(hdr, d$t1, file, sep = ",")
        write("", file = file, append = TRUE)
        write("# Section 2: N-weighted mean per locus (all populations)",
              file = file, append = TRUE)
        write.table(d$t2, file = file, sep = ",", row.names = FALSE,
                    quote = FALSE, append = TRUE, col.names = TRUE)
      }
    )
    
    output$dl_file1_txt <- downloadHandler(
      filename = function() paste0("null_allele_frequencies_", Sys.Date(), ".txt"),
      content  = function(file) {
        d <- file1_data()
        hdr <- c(d$header,
                 "# Section 1: p_nulls per locus x population", "#")
        write_with_header(hdr, d$t1, file, sep = "\t")
        write("", file = file, append = TRUE)
        write("# Section 2: N-weighted mean per locus", file = file, append = TRUE)
        write.table(d$t2, file = file, sep = "\t", row.names = FALSE,
                    quote = FALSE, append = TRUE, col.names = TRUE)
      }
    )

    # ══════════════════════════════════════════════════════════════════════════
    #  FILE 2 — Global FST & FST-ENA with bootstrap CIs
    # ══════════════════════════════════════════════════════════════════════════
    file2_data <- reactive({
      r <- results_r()
      ci_pct <- paste0(round((1-r$alpha)*100,3),"%")
      
      pl <- r$fst_global$per_locus

      glob <- data.frame(
        Locus           = "GLOBAL_MULTILOCUS",
        FST_raw         = round(r$fst_global$global_raw, 6),
        FST_ENA         = round(r$fst_global$global_ena, 6),
        Delta_FST       = round(r$fst_global$global_ena - r$fst_global$global_raw, 6),
        CI_lo_raw_loci  = round(r$boot_gl_loci$raw[1], 6),
        CI_hi_raw_loci  = round(r$boot_gl_loci$raw[3], 6),
        CI_lo_ENA_loci  = round(r$boot_gl_loci$ena[1], 6),
        CI_hi_ENA_loci  = round(r$boot_gl_loci$ena[3], 6),
        CI_lo_raw_subs  = round(r$boot_gl_subs$raw[1], 6),
        CI_hi_raw_subs  = round(r$boot_gl_subs$raw[3], 6),
        CI_lo_ENA_subs  = round(r$boot_gl_subs$ena[1], 6),
        CI_hi_ENA_subs  = round(r$boot_gl_subs$ena[3], 6),
        N_pops_raw = NA_integer_, N_pops_ENA = NA_integer_,
        stringsAsFactors = FALSE
      )

      pl$CI_lo_raw_loci <- NA_real_; pl$CI_hi_raw_loci <- NA_real_
      pl$CI_lo_ENA_loci <- NA_real_; pl$CI_hi_ENA_loci <- NA_real_
      pl$CI_lo_raw_subs <- NA_real_; pl$CI_hi_raw_subs <- NA_real_
      pl$CI_lo_ENA_subs <- NA_real_; pl$CI_hi_ENA_subs <- NA_real_

      out <- rbind(glob[, names(pl)], pl)
      list(header = meta_header(r, "File 2 — Global FST and FST-ENA with bootstrap CIs"),
           data   = out)
    })

    output$dl_file2_csv <- downloadHandler(
      filename = function() paste0("global_FST_ENA_CI_", Sys.Date(), ".csv"),
      content  = function(file) { d <- file2_data()
        write_with_header(d$header, d$data, file, sep=",") }
    )
    output$dl_file2_txt <- downloadHandler(
      filename = function() paste0("global_FST_ENA_CI_", Sys.Date(), ".txt"),
      content  = function(file) { d <- file2_data()
        write_with_header(d$header, d$data, file, sep="\t") }
    )

    # ══════════════════════════════════════════════════════════════════════════
    #  FILE 3 — Pairwise long format
    # ══════════════════════════════════════════════════════════════════════════
    file3_data <- reactive({
      r <- results_r()
      fst_l <- r$fst_pair$long
      dc_l  <- r$dc_pair$long
      bf    <- r$boot_pair_fst
      bd    <- r$boot_pair_dc

      merged <- merge(fst_l, dc_l, by=c("Pop1","Pop2"), all=TRUE)
      if (!is.null(bf) && nrow(bf)>0)
        merged <- merge(merged, bf[,c("Pop1","Pop2","CI_lo_loci","CI_hi_loci",
                                      "CI_lo_raw","CI_hi_raw")],
                        by=c("Pop1","Pop2"), all.x=TRUE,
                        suffixes=c("","_FST_loci"))
      if (!is.null(bd) && nrow(bd)>0)
        merged <- merge(merged, bd[,c("Pop1","Pop2","CI_lo_loci","CI_hi_loci",
                                      "CI_lo_raw","CI_hi_raw")],
                        by=c("Pop1","Pop2"), all.x=TRUE,
                        suffixes=c("_FST","_DCSE"))
      list(header = meta_header(r, "File 3 — Pairwise statistics (all loci combined), long format"),
           data   = merged)
    })

    output$dl_file3_csv <- downloadHandler(
      filename = function() paste0("pairwise_long_format_", Sys.Date(), ".csv"),
      content  = function(file) { d <- file3_data()
        write_with_header(d$header, d$data, file, sep=",") }
    )
    output$dl_file3_txt <- downloadHandler(
      filename = function() paste0("pairwise_long_format_", Sys.Date(), ".txt"),
      content  = function(file) { d <- file3_data()
        write_with_header(d$header, d$data, file, sep="\t") }
    )

    # ══════════════════════════════════════════════════════════════════════════
    #  FILE 4 — Per-locus half-matrices (TXT only - formatted)
    # ══════════════════════════════════════════════════════════════════════════
    output$dl_file4 <- downloadHandler(
      filename = function() paste0("per_locus_half_matrices_", Sys.Date(), ".txt"),
      content  = function(file) {
        r <- results_r()
        hdr <- meta_header(r, "File 4 — Per-locus half-matrices (FST, FST-ENA, DCSE, DCSE-INA)")
        writeLines(hdr, con=file)
        for (loc in r$markers) {
          for (sc in c("FST_raw","FST_ENA")) {
            ln <- half_matrix_txt(r$per_locus_pair$fst, sc, r$pops, loc)
            write(ln, file=file, append=TRUE)
            write("", file=file, append=TRUE)
          }
          for (sc in c("DCSE_raw","DCSE_INA")) {
            ln <- half_matrix_txt(r$per_locus_pair$dc, sc, r$pops, loc)
            write(ln, file=file, append=TRUE)
            write("", file=file, append=TRUE)
          }
        }
      }
    )

    # ── Download buttons UI ────────────────────────────────────────────────────
    output$ui_dl_file1 <- renderUI({
      req(results_r())
      tags$div(class="na-dl-row",
        downloadButton(ns("dl_file1_csv"), ".csv", class="btn btn-default btn-xs"),
        downloadButton(ns("dl_file1_txt"), ".txt", class="btn btn-default btn-xs"))
    })
    
    output$ui_dl_file2 <- renderUI({
      req(results_r())
      tags$div(class="na-dl-row",
        downloadButton(ns("dl_file2_csv"), ".csv", class="btn btn-default btn-xs"),
        downloadButton(ns("dl_file2_txt"), ".txt", class="btn btn-default btn-xs"))
    })
    
    output$ui_dl_file3 <- renderUI({
      req(results_r())
      tags$div(class="na-dl-row",
        downloadButton(ns("dl_file3_csv"), ".csv", class="btn btn-default btn-xs"),
        downloadButton(ns("dl_file3_txt"), ".txt", class="btn btn-default btn-xs"))
    })
    
    output$ui_dl_file4 <- renderUI({
      req(results_r())
      tags$div(class="na-dl-row",
        downloadButton(ns("dl_file4"), ".txt", class="btn btn-default btn-xs"))
    })

    # ── Run status ─────────────────────────────────────────────────────────────
    output$ui_run_status <- renderUI({
      r <- tryCatch(results_r(), error = function(e) NULL)
      if (is.null(r)) return(NULL)
      ci_pct <- paste0(round((1-r$alpha)*100,3),"%")
      tags$div(class="na-info", style="margin-top:.5rem;",
        icon("check-circle"), " ",
        tags$strong("Computation complete."),
        sprintf(" %d loci \u00b7 %d populations \u00b7 %d replicates \u00b7 %s CI.",
                length(r$markers), length(r$pops), r$nboot, ci_pct),
        " Output files are ready for download above."
      )
    })

    # ── Value boxes ────────────────────────────────────────────────────────────
    output$vb_loci <- renderUI({
      tryCatch(tags$span(length(markers_r())), error=function(e) tags$span("\u2014"))
    })
    output$vb_pops <- renderUI({
      tryCatch(tags$span(length(pops_r())), error=function(e) tags$span("\u2014"))
    })
    output$vb_n <- renderUI({
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
        r <- results_r()
        v <- round(mean(r$t1$p_nulls, na.rm=TRUE), 4)
        col <- if(v>.20)"#9d174d" else if(v>.10)"#854d0e" else "#166534"
        tags$span(style=paste0("color:",col,";"), v)
      }, error=function(e) tags$span("\u2014"))
    })
    output$vb_max_null <- renderUI({
      tryCatch({
        r <- results_r()
        v <- round(max(r$t1$p_nulls, na.rm=TRUE), 4)
        col <- if(v>.30)"#9d174d" else if(v>.15)"#854d0e" else "#166534"
        tags$span(style=paste0("color:",col,";"), v)
      }, error=function(e) tags$span("\u2014"))
    })
    output$vb_fst_ena <- renderUI({
      tryCatch({
        r <- results_r(); v <- round(r$fst_global$global_ena, 4)
        col <- if(!is.na(v)&&v>.15)"#9d174d" else if(!is.na(v)&&v>.05)"#854d0e" else "#166534"
        tags$span(style=paste0("color:",col,";"), if(is.na(v))"\u2014" else v)
      }, error=function(e) tags$span("\u2014"))
    })

    # ── Tab 1: null allele frequencies DTs ────────────────────────────────────
    output$dt_t1 <- DT::renderDT({
      r <- results_r()
      shiny::validate(shiny::need(nrow(r$t1)>0, "No data yet. Click Compute."))
      d <- r$t1; names(d) <- c("Locus","Population","Coding","p_nulls","N","N_exp_blanks")
      DT::datatable(d, rownames=FALSE,
        options=list(pageLength=20,scrollX=TRUE,dom="lftip",
          columnDefs=list(list(className="dt-right",targets=3:5))),
        class="compact hover stripe") |>
        DT::formatRound("p_nulls",6) |> DT::formatRound("N_exp_blanks",6) |>
        DT::formatStyle("p_nulls",backgroundColor=DT::styleInterval(
          c(0.05,0.10,0.20,0.30),c("#f0fdf4","#dcfce7","#fefce8","#fff7ed","#fef2f2"))) |>
        DT::formatStyle("Locus",fontWeight="600",color="#0f172a")
    }, server=TRUE)

    output$dt_t2 <- DT::renderDT({
      r <- results_r()
      shiny::validate(shiny::need(nrow(r$t2)>0, "No data yet. Click Compute."))
      d <- r$t2; names(d) <- c("Locus","Coding","Av(p_nulls)","Av(N_exp_blanks)","N_tot")
      DT::datatable(d, rownames=FALSE,
        options=list(pageLength=20,scrollX=TRUE,dom="lftip",
          columnDefs=list(list(className="dt-right",targets=2:4))),
        class="compact hover stripe") |>
        DT::formatRound("Av(p_nulls)",6) |> DT::formatRound("Av(N_exp_blanks)",6) |>
        DT::formatStyle("Av(p_nulls)",backgroundColor=DT::styleInterval(
          c(0.05,0.10,0.20),c("#f0fdf4","#dcfce7","#fefce8","#fef2f2"))) |>
        DT::formatStyle("Locus",fontWeight="600",color="#0f172a")
    }, server=TRUE)

    # ── Tab 2: FST DTs ─────────────────────────────────────────────────────────
    output$dt_fst_global <- DT::renderDT({
      r <- results_r(); d <- r$fst_global$per_locus
      shiny::validate(shiny::need(nrow(d)>0, "No data yet. Click Compute."))
      glob <- data.frame(
        Locus="[GLOBAL MULTILOCUS]",
        FST_raw=round(r$fst_global$global_raw,6),
        FST_ENA=round(r$fst_global$global_ena,6),
        Delta_FST=round(r$fst_global$global_ena-r$fst_global$global_raw,6),
        N_pops_raw=NA_integer_, N_pops_ENA=NA_integer_, stringsAsFactors=FALSE)
      disp <- rbind(glob, d)
      names(disp) <- c("Locus","Raw FST","FST-ENA","\u0394FST","N pops (raw)","N pops (ENA)")
      DT::datatable(disp, rownames=FALSE,
        options=list(pageLength=25,scrollX=TRUE,dom="lftip",
          columnDefs=list(list(className="dt-right",targets=1:5))),
        class="compact hover stripe") |>
        DT::formatRound("Raw FST",6)|>DT::formatRound("FST-ENA",6)|>
        DT::formatRound("\u0394FST",6)|>
        DT::formatStyle("FST-ENA",backgroundColor=DT::styleInterval(
          c(0.05,0.15,0.25),c("#f0fdf4","#dcfce7","#fefce8","#fef2f2")))|>
        DT::formatStyle("Locus",fontWeight="600",color="#0f172a")
    }, server=TRUE)

    # ── Helper for matrix display ──────────────────────────────────────────────
    render_mat_html <- function(mat, fmt=6,
                                thr =c(0.05,0.15,0.25),
                                clrs=c("#f0fdf4","#dcfce7","#fefce8","#fef2f2")) {
      if (is.null(mat)) return(HTML("<p>No data available</p>"))
      pops <- rownames(mat); n <- length(pops)
      cell <- function(i,j) {
        if (i==j) return('<td class="diag">\u2014</td>')
        if (i<j)  return('<td class="upper">\u00b7</td>')
        v <- mat[i,j]; if (is.na(v)) return('<td style="color:#94a3b8;">NA</td>')
        bg <- clrs[findInterval(v,thr)+1L]
        sprintf('<td style="background:%s;">%s</td>',bg,round(v,fmt))
      }
      thead <- paste0('<tr><th></th>',paste(sprintf('<th>%s</th>',pops[-n]),collapse=""),'</tr>')
      tbody <- paste(sapply(seq_len(n),function(i){
        if(i==1L) return("")
        paste0('<tr><td class="lbl">',pops[i],'</td>',
               paste(sapply(seq_len(n),function(j)cell(i,j)),collapse=""),'</tr>')
      }),collapse="")
      HTML(sprintf('<div class="na-matrix-wrap"><table class="na-matrix"><thead>%s</thead><tbody>%s</tbody></table></div>',
                   thead,tbody))
    }

    boot_tbl <- function(d, cols, col_labels, char_cols = c("Pop1","Pop2","Locus")) {
      if (is.null(d) || nrow(d)==0) return(HTML("<p>No bootstrap data available</p>"))
      rows_html <- sapply(seq_len(nrow(d)), function(i) {
        cells <- paste(sapply(cols, function(cn) {
          val <- d[[cn]][i]
          if (cn %in% char_cols)
            sprintf('<td class="lbl">%s</td>', htmltools::htmlEscape(as.character(val)))
          else {
            num <- suppressWarnings(as.numeric(val))
            sprintf('<td>%s</td>', if(is.na(num)) "NA" else formatC(num,digits=6,format="f"))
          }
        }), collapse="")
        paste0("<tr>",cells,"</tr>")
      })
      HTML(paste0(
        '<table class="na-matrix" style="width:100%"><thead><tr>',
        paste(sprintf("<th>%s</th>",col_labels),collapse=""),
        '</tr></thead><tbody>',paste(rows_html,collapse=""),'</tbody></table>'))
    }

    # ── Global bootstrap CI display ────────────────────────────────────────────
    output$ui_boot_global_fst <- renderUI({
      r <- tryCatch(results_r(), error=function(e) NULL)
      if (is.null(r)) return(tags$p("Run computation first.", style="color:#94a3b8;"))
      ci_pct <- paste0(round((1-r$alpha)*100,3),"%")
      bl <- r$boot_gl_loci; bs <- r$boot_gl_subs
      tags$div(class="na-boot-result",
        tags$strong(sprintf("Global FST-ENA \u2014 observed: %.6f", r$fst_global$global_ena)),
        tags$br(),
        sprintf("%s CI (bootstrap over loci):      [ %.6f  \u2013  %.6f ]", ci_pct, bl$ena[1], bl$ena[3]),
        tags$br(),
        sprintf("%s CI (bootstrap over sub-samples): [ %.6f  \u2013  %.6f ]", ci_pct, bs$ena[1], bs$ena[3]),
        tags$br(), tags$br(),
        tags$strong(sprintf("Global Raw FST \u2014 observed: %.6f", r$fst_global$global_raw)),
        tags$br(),
        sprintf("%s CI (bootstrap over loci):      [ %.6f  \u2013  %.6f ]", ci_pct, bl$raw[1], bl$raw[3]),
        tags$br(),
        sprintf("%s CI (bootstrap over sub-samples): [ %.6f  \u2013  %.6f ]", ci_pct, bs$raw[1], bs$raw[3])
      )
    })

    # ── Pairwise FST matrix ────────────────────────────────────────────────────
    output$ui_fst_pair_matrix <- renderUI({
      r <- tryCatch(results_r(), error=function(e) NULL)
      if (is.null(r)||is.null(r$fst_pair$matrix_raw))
        return(tags$p("Run computation first.", style="color:#94a3b8;"))
      typ <- input$fst_pair_display %||% "both"
      if (identical(typ,"both")) tags$div(
        tags$p(tags$strong("Raw FST")),   render_mat_html(r$fst_pair$matrix_raw), tags$br(),
        tags$p(tags$strong("FST-ENA")),   render_mat_html(r$fst_pair$matrix_ena))
      else if (identical(typ,"raw")) render_mat_html(r$fst_pair$matrix_raw)
      else render_mat_html(r$fst_pair$matrix_ena)
    })

    output$ui_boot_pair_fst <- renderUI({
      r <- tryCatch(results_r(), error=function(e) NULL)
      if (is.null(r)||is.null(r$boot_pair_fst))
        return(tags$p("Run computation first.", style="color:#94a3b8;"))
      ci_pct <- paste0(round((1-r$alpha)*100,3),"%")
      tags$div(class="na-boot-result",
        tags$p(tags$strong(sprintf("Pairwise FST-ENA \u2014 %s CI (bootstrap over loci)", ci_pct))),
        tags$div(class="na-matrix-wrap",
          boot_tbl(r$boot_pair_fst,
            c("Pop1","Pop2","FST_ENA_obs","CI_lo_loci","Median_loci","CI_hi_loci",
              "FST_raw_obs","CI_lo_raw","CI_hi_raw"),
            c("Pop 1","Pop 2","FST-ENA obs.","CI lo","Median","CI hi",
              "Raw FST obs.","CI lo (raw)","CI hi (raw)")))
      )
    })

    # ── Pairwise DCSE matrix + bootstrap ──────────────────────────────────────
    output$ui_dc_matrix <- renderUI({
      r <- tryCatch(results_r(), error=function(e) NULL)
      if (is.null(r)||is.null(r$dc_pair$matrix_raw))
        return(tags$p("Run computation first.", style="color:#94a3b8;"))
      typ <- input$dc_display %||% "both"
      thr <- c(0.1,0.25,0.4); clrs <- c("#eff6ff","#dbeafe","#fef9c3","#fef2f2")
      if (identical(typ,"both")) tags$div(
        tags$p(tags$strong("Raw DCSE")),  render_mat_html(r$dc_pair$matrix_raw,thr=thr,clrs=clrs), tags$br(),
        tags$p(tags$strong("DCSE-INA")), render_mat_html(r$dc_pair$matrix_ina,thr=thr,clrs=clrs))
      else if (identical(typ,"raw")) render_mat_html(r$dc_pair$matrix_raw,thr=thr,clrs=clrs)
      else render_mat_html(r$dc_pair$matrix_ina,thr=thr,clrs=clrs)
    })

    output$ui_boot_pair_dc <- renderUI({
      r <- tryCatch(results_r(), error=function(e) NULL)
      if (is.null(r)||is.null(r$boot_pair_dc))
        return(tags$p("Run computation first.", style="color:#94a3b8;"))
      ci_pct <- paste0(round((1-r$alpha)*100,3),"%")
      tags$div(class="na-boot-result",
        tags$p(tags$strong(sprintf("Pairwise DCSE-INA \u2014 %s CI (bootstrap over loci)", ci_pct))),
        tags$div(class="na-matrix-wrap",
          boot_tbl(r$boot_pair_dc,
            c("Pop1","Pop2","DCSE_INA_obs","CI_lo_loci","Median_loci","CI_hi_loci",
              "DCSE_raw_obs","CI_lo_raw","CI_hi_raw"),
            c("Pop 1","Pop 2","DCSE-INA obs.","CI lo","Median","CI hi",
              "Raw DCSE obs.","CI lo (raw)","CI hi (raw)")))
      )
    })

    # ── Tab 4: per-locus x pair DTs ───────────────────────────────────────────
    output$dt_fst_locus <- DT::renderDT({
      r <- tryCatch(results_r(), error=function(e) NULL)
      shiny::validate(shiny::need(!is.null(r), "Run computation first."))
      d <- r$per_locus_pair$fst
      sl <- safe_choice(input$fl_locus,"all")
      sp1 <- safe_choice(input$fl_pop1,"all"); sp2 <- safe_choice(input$fl_pop2,"all")
      if (!identical(sl,"all"))  d <- d[d$Locus==sl,,drop=FALSE]
      if (!identical(sp1,"all")) d <- d[d$Pop1==sp1|d$Pop2==sp1,,drop=FALSE]
      if (!identical(sp2,"all")) d <- d[d$Pop2==sp2|d$Pop1==sp2,,drop=FALSE]
      shiny::validate(shiny::need(nrow(d)>0,"No data for selected filters."))
      names(d) <- c("Locus","Pop 1","Pop 2","Raw FST","FST-ENA")
      DT::datatable(d, rownames=FALSE,
        options=list(pageLength=25,scrollX=TRUE,dom="lftip",
          columnDefs=list(list(className="dt-right",targets=3:4))),
        class="compact hover stripe") |>
        DT::formatRound("Raw FST",6)|>DT::formatRound("FST-ENA",6)|>
        DT::formatStyle("FST-ENA",backgroundColor=DT::styleInterval(
          c(0.05,0.15,0.25),c("#f0fdf4","#dcfce7","#fefce8","#fef2f2")))|>
        DT::formatStyle("Locus",fontWeight="600",color="#0f172a")
    }, server=TRUE)

    output$dt_dc_locus <- DT::renderDT({
      r <- tryCatch(results_r(), error=function(e) NULL)
      shiny::validate(shiny::need(!is.null(r), "Run computation first."))
      d <- r$per_locus_pair$dc
      sl <- safe_choice(input$fl_locus,"all")
      sp1 <- safe_choice(input$fl_pop1,"all"); sp2 <- safe_choice(input$fl_pop2,"all")
      if (!identical(sl,"all"))  d <- d[d$Locus==sl,,drop=FALSE]
      if (!identical(sp1,"all")) d <- d[d$Pop1==sp1|d$Pop2==sp1,,drop=FALSE]
      if (!identical(sp2,"all")) d <- d[d$Pop2==sp2|d$Pop1==sp2,,drop=FALSE]
      shiny::validate(shiny::need(nrow(d)>0,"No data for selected filters."))
      names(d) <- c("Locus","Pop 1","Pop 2","Raw DCSE","DCSE-INA")
      DT::datatable(d, rownames=FALSE,
        options=list(pageLength=25,scrollX=TRUE,dom="lftip",
          columnDefs=list(list(className="dt-right",targets=3:4))),
        class="compact hover stripe") |>
        DT::formatRound("Raw DCSE",6)|>DT::formatRound("DCSE-INA",6)|>
        DT::formatStyle("Locus",fontWeight="600",color="#0f172a")
    }, server=TRUE)

  })
}