# server_isolation_by_distance.R
# IBD: Rousset (1997) FST/(1-FST) vs geographic distance + Mantel test.
# Based on FreeNA ENA-corrected FST (Chapuis & Estoup 2007)

server_isolation_by_distance <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    # ── Cross-platform parallel helper ────────────────────────────────────────
    .parallel_boot <- function(n_boot, fun) {
      if (.Platform$OS.type == "windows") {
        n_cores <- max(1L, parallel::detectCores() - 1L)
        cl <- parallel::makeCluster(n_cores)
        on.exit(parallel::stopCluster(cl))
        env_funs <- ls(environment(fun), all.names = TRUE)
        parallel::clusterExport(cl, env_funs, envir = environment(fun))
        parallel::clusterExport(cl, "fun", envir = environment())
        results <- parallel::parLapply(cl, seq_len(n_boot), function(i) fun())
      } else {
        n_cores <- max(1L, parallel::detectCores() - 1L)
        results <- parallel::mclapply(seq_len(n_boot), function(i) fun(), 
                                       mc.cores = n_cores)
      }
      results
    }

    # ── Helpers ────────────────────────────────────────────────────────────────
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
    sql_id   <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))
    ci_bounds <- function(alpha) {
      list(lo = alpha/2, hi = 1-alpha/2,
           label = paste0(round((1-alpha)*100, 3), "% CI"))
    }

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
      lat_col <- c("Latitude","latitude")[c("Latitude","latitude") %in% cols][1]
      lon_col <- c("Longitude","longitude")[c("Longitude","longitude") %in% cols][1]
      list(ind_col = ind_col, pop_col = pop_col, lat_col = lat_col, lon_col = lon_col)
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

    # ── GPS coordinates ─────────────────────────────────────────────────────────
    coords_r <- reactive({
      db_ready(); con <- con_r(); ms <- meta_schema_r()
      df <- DBI::dbGetQuery(con, sprintf(
        "SELECT CAST(%s AS VARCHAR) AS Population,
                AVG(CAST(%s AS DOUBLE)) AS Latitude,
                AVG(CAST(%s AS DOUBLE)) AS Longitude
         FROM %s
         WHERE %s IS NOT NULL
           AND %s IS NOT NULL AND %s IS NOT NULL
         GROUP BY %s ORDER BY Population",
        sql_id(con,ms$pop_col),
        sql_id(con,ms$lat_col), sql_id(con,ms$lon_col),
        sql_id(con,tbl_meta_r()),
        sql_id(con,ms$pop_col),
        sql_id(con,ms$lat_col), sql_id(con,ms$lon_col),
        sql_id(con,ms$pop_col)))
      shiny::validate(shiny::need(nrow(df) >= 2L,
        "At least 2 populations with GPS coordinates are required."))
      df
    })

    # ── Haversine distance ─────────────────────────────────────────────────────
    .haversine_km <- function(lat1, lon1, lat2, lon2) {
      R    <- 6371.0
      dlat <- (lat2 - lat1) * pi / 180
      dlon <- (lon2 - lon1) * pi / 180
      a    <- sin(dlat / 2)^2 +
        cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2)^2
      2 * R * asin(sqrt(a))
    }

    .geo_dist_matrix <- function(coords) {
      n   <- nrow(coords)
      mat <- matrix(0.0, n, n,
                    dimnames = list(coords$Population, coords$Population))
      for (i in seq_len(n - 1)) {
        for (j in (i + 1):n) {
          d <- .haversine_km(coords$Latitude[i],  coords$Longitude[i],
                             coords$Latitude[j],  coords$Longitude[j])
          mat[i, j] <- mat[j, i] <- d
        }
      }
      mat
    }

    # ── EM null allele frequency (FreeNA Pascal translation) ──────────────────
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

      # rd initialisation
      rd <- if (treat == "null_homo" && n_null_homo > 0L)
              sqrt(n_null_homo / N)
            else
              sqrt(1.0 / (N + 1.0))

      # corrdgenefreq initialisation (cpt=0 Pascal)
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

      # EM loop
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

    # ── WC84 ENA components (FreeNA Pascal translation) ──────────────────────
    .wc84_ena_2pop <- function(em1, em2) {
      K    <- 2L
      N_k  <- c(em1$N, em2$N)
      if (any(N_k < 2L)) return(NULL)
      n_total <- sum(N_k)
      n_bar   <- n_total / K
      n_c     <- (n_total - sum(N_k^2) / n_total) / (K - 1L)
      if (n_c <= 0 || n_bar <= 1) return(NULL)

      all_alleles <- unique(c(em1$alleles_uniq, em2$alleles_uniq))
      if (length(all_alleles) == 0L) return(NULL)

      a_sum <- 0; b_sum <- 0; c_sum <- 0

      for (allele in all_alleles) {
        p1  <- if (allele %in% names(em1$cq)) em1$cq[[allele]] else 0.0
        p2  <- if (allele %in% names(em2$cq)) em2$cq[[allele]] else 0.0
        pbar <- (N_k[1L] * p1 + N_k[2L] * p2) / n_total
        s2   <- (N_k[1L] * (p1 - pbar)^2 + N_k[2L] * (p2 - pbar)^2) /
                  ((K - 1L) * n_bar)
        h1 <- if (allele %in% names(em1$n_het)) em1$n_het[[allele]] else 0L
        h2 <- if (allele %in% names(em2$n_het)) em2$n_het[[allele]] else 0L
        hbar <- (h1 + h2) / n_total
        term <- pbar * (1 - pbar)
        a_sum <- a_sum + (n_bar / n_c) *
                   (s2 - (1 / (n_bar - 1)) * (term - (K-1)/K * s2 - 0.25 * hbar))
        b_sum <- b_sum + (n_bar / (n_bar - 1)) *
                   (term - (K-1)/K * s2 - (2 * n_bar - 1) / (4 * n_bar) * hbar)
        c_sum <- c_sum + 0.5 * hbar
      }
      list(a = a_sum, b = b_sum, c = c_sum)
    }

    # ── Linearise ──────────────────────────────────────────────────────────────
    .linearise <- function(x) {
      if (is.na(x)) return(NA_real_)
      x <- pmin(pmax(x, 0), 0.9999)
      x / (1 - x)
    }
    .linearise_vec <- function(x) vapply(x, .linearise, numeric(1))

    # ── Pairwise FreeNA FST (ENA corrected) with bootstrap ────────────────────
    .pairwise_freena_fst_ci <- function(hap_df, pop_vector,
                                         null_code = "999999",
                                         n_boot = 500L, conf = 0.95) {
      pops    <- sort(unique(pop_vector))
      loci    <- colnames(hap_df)
      n_pops  <- length(pops)
      alpha   <- (1 - conf) / 2

      .fst_pair <- function(idx1, idx2) {
        a_tot <- 0; b_tot <- 0; c_tot <- 0; ok <- 0L
        for (locus in loci) {
          em1  <- em_freena(hap_df[[locus]][idx1], 1000L, "absent")
          em2  <- em_freena(hap_df[[locus]][idx2], 1000L, "absent")
          comp <- .wc84_ena_2pop(em1, em2)
          if (is.null(comp)) next
          a_tot <- a_tot + comp$a; b_tot <- b_tot + comp$b; c_tot <- c_tot + comp$c
          ok <- ok + 1L
        }
        if (ok == 0L || (a_tot + b_tot + c_tot) == 0) return(NA_real_)
        a_tot / (a_tot + b_tot + c_tot)
      }

      result <- vector("list", n_pops * (n_pops - 1L) / 2L)
      k <- 1L

      for (i in seq_len(n_pops - 1L)) {
        for (j in (i + 1L):n_pops) {
          pi <- pops[i]; pj <- pops[j]
          idx1 <- which(pop_vector == pi)
          idx2 <- which(pop_vector == pj)
          if (length(idx1) < 2L || length(idx2) < 2L) {
            result[[k]] <- data.frame(pop_i = pi, pop_j = pj,
                                       fst = NA_real_, ci_l = NA_real_, ci_u = NA_real_,
                                       stringsAsFactors = FALSE)
            k <- k + 1L; next
          }

          fst_obs <- .fst_pair(idx1, idx2)

          boot_fun <- function() {
            b1 <- sample(idx1, length(idx1), replace = TRUE)
            b2 <- sample(idx2, length(idx2), replace = TRUE)
            .fst_pair(b1, b2)
          }
          
          boot_list <- .parallel_boot(n_boot, boot_fun)
          boot_fst <- unlist(boot_list, use.names = FALSE)
          boot_fst <- boot_fst[is.finite(boot_fst)]
          
          ci_l <- ci_u <- NA_real_
          if (length(boot_fst) > 0L) {
            q    <- quantile(boot_fst, probs = c(alpha, 1 - alpha), na.rm = TRUE)
            ci_l <- q[[1L]]; ci_u <- q[[2L]]
          }
          result[[k]] <- data.frame(pop_i = pi, pop_j = pj,
                                     fst = fst_obs, ci_l = ci_l, ci_u = ci_u,
                                     stringsAsFactors = FALSE)
          k <- k + 1L
        }
      }
      do.call(rbind, result)
    }

    # ── Bootstrap over loci (FreeNA approach) ─────────────────────────────────
    .bootstrap_over_loci <- function(hap_df, pop_vector, null_code = "999999",
                                      n_boot = 1000, conf = 0.95) {
      pops    <- sort(unique(pop_vector))
      loci    <- colnames(hap_df)
      n_loci  <- length(loci)
      n_pops  <- length(pops)
      alpha   <- (1 - conf) / 2

      .fst_loci <- function(idx1, idx2, locus_indices) {
        a_tot <- 0; b_tot <- 0; c_tot <- 0; ok <- 0L
        for (locus in locus_indices) {
          em1  <- em_freena(hap_df[[loci[locus]]][idx1], 1000L, "absent")
          em2  <- em_freena(hap_df[[loci[locus]]][idx2], 1000L, "absent")
          comp <- .wc84_ena_2pop(em1, em2)
          if (is.null(comp)) next
          a_tot <- a_tot + comp$a; b_tot <- b_tot + comp$b; c_tot <- c_tot + comp$c
          ok <- ok + 1L
        }
        if (ok == 0L || (a_tot + b_tot + c_tot) == 0) return(NA_real_)
        a_tot / (a_tot + b_tot + c_tot)
      }

      pop_indices <- lapply(pops, function(p) which(pop_vector == p))
      names(pop_indices) <- pops

      result <- vector("list", n_pops * (n_pops - 1L) / 2L)
      k <- 1L

      for (i in seq_len(n_pops - 1L)) {
        for (j in (i + 1L):n_pops) {
          pi <- pops[i]; pj <- pops[j]
          idx1 <- pop_indices[[pi]]
          idx2 <- pop_indices[[pj]]

          if (length(idx1) < 2L || length(idx2) < 2L) {
            result[[k]] <- data.frame(
              pop_i = pi, pop_j = pj,
              fst = NA_real_, fst_ci_l = NA_real_, fst_ci_u = NA_real_,
              fr = NA_real_, fr_ci_l = NA_real_, fr_ci_u = NA_real_,
              n_loci_used = 0L,
              stringsAsFactors = FALSE
            )
            k <- k + 1L; next
          }

          fst_obs <- .fst_loci(idx1, idx2, seq_len(n_loci))

          boot_fun <- function() {
            locus_idx <- sample(seq_len(n_loci), n_loci, replace = TRUE)
            .fst_loci(idx1, idx2, locus_idx)
          }
          
          boot_list <- .parallel_boot(n_boot, boot_fun)
          boot_fst <- unlist(boot_list, use.names = FALSE)
          boot_fst <- boot_fst[is.finite(boot_fst)]

          fst_ci_l <- fst_ci_u <- NA_real_
          if (length(boot_fst) > 0L) {
            q <- quantile(boot_fst, probs = c(alpha, 1 - alpha), na.rm = TRUE)
            fst_ci_l <- q[[1L]]; fst_ci_u <- q[[2L]]
          }

          fr_obs <- .linearise(fst_obs)
          fr_boot <- .linearise_vec(boot_fst)
          fr_boot <- fr_boot[is.finite(fr_boot)]

          fr_ci_l <- fr_ci_u <- NA_real_
          if (length(fr_boot) > 0L) {
            q <- quantile(fr_boot, probs = c(alpha, 1 - alpha), na.rm = TRUE)
            fr_ci_l <- q[[1L]]; fr_ci_u <- q[[2L]]
          }

          result[[k]] <- data.frame(
            pop_i = pi, pop_j = pj,
            fst = fst_obs, fst_ci_l = fst_ci_l, fst_ci_u = fst_ci_u,
            fr = fr_obs, fr_ci_l = fr_ci_l, fr_ci_u = fr_ci_u,
            n_loci_used = n_loci,
            stringsAsFactors = FALSE
          )
          k <- k + 1L
        }
      }
      do.call(rbind, result)
    }

    # ── Mantel test ────────────────────────────────────────────────────────────
    .mantel_test <- function(fr_vec, dist_vec, n_perm = 9999, use_log = FALSE) {
      x <- if (use_log) log(dist_vec) else dist_vec
      ok <- is.finite(fr_vec) & is.finite(x)
      y  <- fr_vec[ok]; x <- x[ok]
      n  <- length(y)
      if (n < 3L)
        return(list(r = NA_real_, p_value = NA_real_, n_pairs = 0L,
                    slope = NA_real_, intercept = NA_real_, perm_r = NA_real_))

      r_obs <- cor(y, x)
      lm0   <- lm(y ~ x)

      perm_r <- vapply(seq_len(n_perm), function(.i) {
        xp <- sample(x)
        suppressWarnings(cor(y, xp))
      }, numeric(1))
      perm_r  <- perm_r[is.finite(perm_r)]
      p_value <- if (length(perm_r) > 0) mean(perm_r >= r_obs) else NA_real_

      list(r         = r_obs,
           p_value   = p_value,
           n_pairs   = n,
           slope     = unname(coef(lm0)[2L]),
           intercept = unname(coef(lm0)[1L]),
           perm_r    = perm_r)
    }

    # ── Regression parameters ──────────────────────────────────────────────────
    .reg_params <- function(y, x, use_log) {
      xv <- if (use_log) log(x) else x
      ok <- is.finite(y) & is.finite(xv)
      if (sum(ok) < 3L) return(c(intercept = NA_real_, b = NA_real_, 
                                 Nb = NA_real_, Nem = NA_real_))
      cf <- coef(lm(y[ok] ~ xv[ok]))
      b  <- unname(cf[2L])
      list(
        intercept = unname(cf[1L]),
        b         = b,
        Nb        = if (is.finite(b) && b > 0) 1 / b else NA_real_,
        Nem       = if (is.finite(b) && b > 0) 1 / (2 * pi * b) else NA_real_
      )
    }

    # ── Raw genotype data ──────────────────────────────────────────────────────
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
    #  MAIN IBD RESULTS
    # ══════════════════════════════════════════════════════════════════════════
    ibd_results_r <- eventReactive(input$run_ibd, {
      shiny::req(db_ready())
      coords  <- coords_r()
      raw_df  <- raw_data_r()
      markers <- markers_r()
      pops    <- pops_r()
      
      use_log <- (input$model == "2D")
      n_boot_pw <- as.integer(input$n_boot_pw)
      n_boot_loci <- as.integer(input$n_boot_loci)
      alpha  <- as.numeric(input$boot_ci_level %||% "0.05")
      ci     <- ci_bounds(alpha)

      # Build hap_df from raw data
      hap_df <- data.frame(row.names = seq_len(nrow(raw_df)))
      for (loc in markers) {
        gts <- raw_df$gt[raw_df$Marker == loc]
        hap_df[[loc]] <- as.character(gts)
      }
      pop_vector <- raw_df$Population

      # Keep only populations with GPS data
      pops_with_gps <- coords$Population
      keep_pop      <- pop_vector %in% pops_with_gps
      shiny::validate(shiny::need(
        length(unique(pop_vector[keep_pop])) >= 2L,
        "Less than 2 populations have both genotypes and GPS data."))

      hap_sub <- hap_df[keep_pop, , drop = FALSE]
      pop_sub <- pop_vector[keep_pop]
      new_levels <- sort(unique(pop_sub))

      # Pairwise FreeNA ENA-corrected FST + CI (individual bootstrap)
      withProgress(message = "Computing FreeNA ENA-corrected pairwise FST...", value = 0.2, {
        pw <- .pairwise_freena_fst_ci(hap_sub, pop_sub, n_boot = n_boot_pw, conf = 1 - alpha)
      })

      # Bootstrap over loci
      withProgress(message = "Bootstrapping over loci...", value = 0.5, {
        boot_loci <- .bootstrap_over_loci(hap_sub, pop_sub, n_boot = n_boot_loci, conf = 1 - alpha)
      })

      # Geographic distances
      coords_ord <- coords[match(new_levels, coords$Population), ]
      coords_ord <- coords_ord[!is.na(coords_ord$Population), ]
      dist_mat   <- .geo_dist_matrix(coords_ord)

      get_dist <- function(p1, p2) {
        if (p1 %in% rownames(dist_mat) && p2 %in% rownames(dist_mat))
          dist_mat[p1, p2] else NA_real_
      }
      pw$dist_km <- mapply(get_dist, pw$pop_i, pw$pop_j)

      # Linearised FST
      pw$FR   <- .linearise_vec(pw$fst)
      pw$FR_i <- .linearise_vec(pw$ci_l)
      pw$FR_s <- .linearise_vec(pw$ci_u)

      # Merge with bootstrap over loci results
      boot_loci$dist_km <- mapply(get_dist, boot_loci$pop_i, boot_loci$pop_j)
      boot_loci$FR <- .linearise_vec(boot_loci$fst)

      # Regression fits
      reg_avg <- .reg_params(pw$FR,   pw$dist_km, use_log)
      reg_ls  <- .reg_params(pw$FR_s, pw$dist_km, use_log)
      reg_li  <- .reg_params(pw$FR_i, pw$dist_km, use_log)

      # Mantel test (on average F_R)
      withProgress(message = "Running Mantel test...", value = 0.8, {
        mantel <- .mantel_test(pw$FR, pw$dist_km, n_perm = 1000, use_log = use_log)
      })

      list(
        pw        = pw,
        boot_loci = boot_loci,
        dist_mat  = dist_mat,
        reg_avg   = reg_avg,
        reg_ls    = reg_ls,
        reg_li    = reg_li,
        mantel    = mantel,
        use_log   = use_log,
        alpha     = alpha,
        ci        = ci,
        markers   = markers,
        pops      = new_levels,
        n_loci    = length(markers)
      )
    })

    # ── Mantel test results ────────────────────────────────────────────────────
    mantel_results_r <- eventReactive(input$run_mantel, {
      shiny::req(ibd_results_r())
      res <- ibd_results_r()
      n_perm <- as.integer(input$n_perm_mantel)
      use_log <- (input$mantel_model == "2D")
      
      withProgress(message = "Running Mantel test...", value = 0.3, {
        mantel <- .mantel_test(res$pw$FR, res$pw$dist_km, n_perm = n_perm, use_log = use_log)
      })
      
      mantel
    })

    # ── Value boxes ────────────────────────────────────────────────────────────
    output$vb_loci <- renderUI({
      tryCatch(tags$span(length(markers_r())), error=function(e) tags$span("\u2014"))
    })
    output$vb_pops <- renderUI({
      tryCatch(tags$span(length(pops_r())), error=function(e) tags$span("\u2014"))
    })
    output$vb_pairs <- renderUI({
      tryCatch({
        r <- ibd_results_r()
        tags$span(nrow(r$pw))
      }, error=function(e) tags$span("\u2014"))
    })
    output$vb_mantel_r <- renderUI({
      tryCatch({
        r <- ibd_results_r()
        v <- round(r$mantel$r, 4)
        tags$span(if (is.na(v)) "\u2014" else v)
      }, error=function(e) tags$span("\u2014"))
    })
    output$vb_nb <- renderUI({
      tryCatch({
        r <- ibd_results_r()
        v <- r$reg_avg$Nb
        tags$span(if (is.na(v)) "\u2014" else round(v, 1))
      }, error=function(e) tags$span("\u2014"))
    })

    # ── Regression table ───────────────────────────────────────────────────────
    output$reg_table <- DT::renderDT({
      r <- ibd_results_r()
      fmt6 <- function(x) if (is.na(x)) "NA" else formatC(x, format = "f", digits = 4)
      fmt1 <- function(x) if (is.na(x)) "NA" else formatC(x, format = "f", digits = 1)
      df <- data.frame(
        Line      = c("Average", "Upper CI (ls)", "Lower CI (li)"),
        b         = sapply(list(r$reg_avg, r$reg_ls, r$reg_li), function(x) fmt6(x$b)),
        Intercept = sapply(list(r$reg_avg, r$reg_ls, r$reg_li), function(x) fmt6(x$intercept)),
        Nb        = sapply(list(r$reg_avg, r$reg_ls, r$reg_li), function(x) fmt1(x$Nb)),
        Nem       = sapply(list(r$reg_avg, r$reg_ls, r$reg_li), function(x) fmt1(x$Nem)),
        stringsAsFactors = FALSE
      )
      DT::datatable(df, rownames = FALSE,
        options = list(dom = "t", pageLength = 3, ordering = FALSE),
        class   = "compact stripe") %>%
        DT::formatStyle("Line",
          target = "row",
          backgroundColor = DT::styleEqual(
            c("Average", "Upper CI (ls)", "Lower CI (li)"),
            c("#f5f5f5",  "#fff0f0",       "#f0f8ff")))
    })

    output$ui_run_status <- renderUI({
      r <- tryCatch(ibd_results_r(), error = function(e) NULL)
      if (is.null(r)) return(NULL)
      tags$div(class="ibd-info", style="margin-top:.5rem;",
        icon("check-circle"), " ",
        tags$strong("Computation complete."),
        sprintf(" %d loci \u00b7 %d populations \u00b7 %d pairs \u00b7 %s model.",
                r$n_loci, length(r$pops), nrow(r$pw),
                if(r$use_log) "2D (ln)" else "1D (linear)")
      )
    })

    # ── IBD plot ───────────────────────────────────────────────────────────────
    output$ibd_plot <- plotly::renderPlotly({
      r      <- ibd_results_r()
      pw     <- r$pw
      use_log <- r$use_log

      xv      <- if (use_log) log(pw$dist_km) else pw$dist_km
      x_label <- if (use_log) "ln(geographic distance, km)" else "Geographic distance (km)"
      x_seq   <- seq(min(xv, na.rm = TRUE), max(xv, na.rm = TRUE), length.out = 100)

      reg_line <- function(params) {
        if (any(is.na(c(params$intercept, params$b)))) return(NULL)
        data.frame(x = x_seq, y = params$intercept + params$b * x_seq)
      }
      line_avg <- reg_line(r$reg_avg)
      line_ls  <- reg_line(r$reg_ls)
      line_li  <- reg_line(r$reg_li)

      p <- plotly::plot_ly() %>%
        # Error bars (CI segment per point)
        plotly::add_segments(
          data = pw,
          x = ~xv, xend = ~xv, y = ~FR_i, yend = ~FR_s,
          line = list(color = "rgba(100,100,100,0.35)", width = 1),
          showlegend = FALSE, hoverinfo = "none"
        ) %>%
        # Scatter points
        plotly::add_markers(
          data = data.frame(x = xv, y = pw$FR,
                             pop_i = pw$pop_i, pop_j = pw$pop_j,
                             dist  = round(pw$dist_km, 2),
                             fst   = round(pw$fst, 5),
                             fr    = round(pw$FR, 5)),
          x = ~x, y = ~y,
          text = ~paste0(pop_i, " \u2013 ", pop_j,
                         "<br>Dist: ", dist, " km",
                         "<br>F\u209b\u209c-ENA: ", fst,
                         "<br>F\u209b\u209c/(1\u2212F\u209b\u209c): ", fr),
          hoverinfo = "text",
          marker = list(color = "#2CBF9F", size = 7, opacity = 0.85),
          name = "Pairs"
        ) %>%
        { if (!is.null(line_avg))
            plotly::add_lines(., data = line_avg, x = ~x, y = ~y,
              line = list(color = "#333a43", width = 2, dash = "solid"),
              name = paste0("Average  b=", formatC(r$reg_avg$b, format="f", digits=4)))
          else . } %>%
        { if (!is.null(line_ls))
            plotly::add_lines(., data = line_ls, x = ~x, y = ~y,
              line = list(color = "#B40F20", width = 1.5, dash = "dash"),
              name = paste0("Upper CI (ls)  b=", formatC(r$reg_ls$b, format="f", digits=4)))
          else . } %>%
        { if (!is.null(line_li))
            plotly::add_lines(., data = line_li, x = ~x, y = ~y,
              line = list(color = "#3B9AB2", width = 1.5, dash = "dot"),
              name = paste0("Lower CI (li)  b=", formatC(r$reg_li$b, format="f", digits=4)))
          else . } %>%
        plotly::layout(
          xaxis  = list(title = x_label),
          yaxis  = list(title = "F\u209b\u209c / (1 \u2212 F\u209b\u209c)"),
          title  = list(
            text = paste0("IBD (Rousset 1997) \u2014 Mantel r = ",
                          round(r$mantel$r, 4),
                          ", p = ",
                          formatC(r$mantel$p_value, format = "f", digits = 4)),
            font = list(size = 13)),
          legend = list(x = 0.02, y = 0.98, bgcolor = "rgba(255,255,255,0.8)",
                        bordercolor = "#ddd", borderwidth = 1),
          font   = list(family = "Helvetica Neue, Segoe UI, Arial"),
          margin = list(t = 55)
        )
    })

    # ── Pairwise table with bootstrap CIs ─────────────────────────────────────
    output$pairwise_table <- DT::renderDT({
      r <- ibd_results_r()
      pw <- r$pw
      boot <- r$boot_loci
      
      # Merge individual and locus bootstrap results
      df <- data.frame(
        Pop1 = pw$pop_i,
        Pop2 = pw$pop_j,
        Dist_km = round(pw$dist_km, 2),
        FST_ENA = round(pw$fst, 5),
        FST_CI_low_ind = round(pw$ci_l, 5),
        FST_CI_high_ind = round(pw$ci_u, 5),
        FST_CI_low_loci = round(boot$fst_ci_l, 5),
        FST_CI_high_loci = round(boot$fst_ci_u, 5),
        FR = round(pw$FR, 5),
        FR_CI_low_ind = round(pw$FR_i, 5),
        FR_CI_high_ind = round(pw$FR_s, 5),
        FR_CI_low_loci = round(boot$fr_ci_l, 5),
        FR_CI_high_loci = round(boot$fr_ci_u, 5),
        N_loci = boot$n_loci_used,
        stringsAsFactors = FALSE
      )
      
      DT::datatable(df, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15, dom = "lrtip"),
        class = "compact stripe hover",
        colnames = c("Pop 1", "Pop 2", "Dist (km)",
                     "FST-ENA", "FST CI low (ind)", "FST CI high (ind)",
                     "FST CI low (loci)", "FST CI high (loci)",
                     "FR", "FR CI low (ind)", "FR CI high (ind)",
                     "FR CI low (loci)", "FR CI high (loci)",
                     "N loci")
      ) %>%
        DT::formatStyle("FST_ENA",
          backgroundColor = DT::styleInterval(
            c(0.05, 0.15, 0.25),
            c("#d4edda", "#fff3cd", "#f8d7da", "#c3002f22")))
    })

    output$ui_dl_pairwise <- renderUI({
      req(ibd_results_r())
      ns_local <- session$ns
      tags$div(class="ibd-dl-row",
        downloadButton(ns_local("dl_pairwise"), "Download table", class = "btn btn-default btn-sm"))
    })

    output$dl_pairwise <- downloadHandler(
      filename = function() paste0("IBD_pairwise_", Sys.Date(), ".csv"),
      content  = function(file) {
        r <- ibd_results_r()
        pw <- r$pw
        boot <- r$boot_loci
        df <- data.frame(
          Pop1 = pw$pop_i, Pop2 = pw$pop_j,
          Dist_km = round(pw$dist_km, 4),
          FST_ENA = round(pw$fst, 6),
          FST_CI_low_ind = round(pw$ci_l, 6),
          FST_CI_high_ind = round(pw$ci_u, 6),
          FST_CI_low_loci = round(boot$fst_ci_l, 6),
          FST_CI_high_loci = round(boot$fst_ci_u, 6),
          FR = round(pw$FR, 6),
          FR_CI_low_ind = round(pw$FR_i, 6),
          FR_CI_high_ind = round(pw$FR_s, 6),
          FR_CI_low_loci = round(boot$fr_ci_l, 6),
          FR_CI_high_loci = round(boot$fr_ci_u, 6),
          N_loci = boot$n_loci_used,
          stringsAsFactors = FALSE)
        write.csv(df, file, row.names = FALSE)
      }
    )

    # ── Mantel test outputs ────────────────────────────────────────────────────
    output$ui_mantel_status <- renderUI({
      m <- tryCatch(mantel_results_r(), error = function(e) NULL)
      if (is.null(m)) return(NULL)
      tags$div(class="ibd-boot-result",
        tags$strong("Mantel test complete."),
        tags$br(),
        sprintf("r = %.6f, p = %.6f, n = %d pairs", m$r, m$p_value, m$n_pairs)
      )
    })

    # Update Mantel results display
    observe({
      m <- tryCatch(mantel_results_r(), error = function(e) NULL)
      if (!is.null(m)) {
        # Update text values using JavaScript
        session$sendCustomMessage("updateMantelResults", list(
          r = round(m$r, 6),
          p = m$p_value,
          n = m$n_pairs,
          b = round(m$slope, 6),
          intercept = round(m$intercept, 6),
          nb = if (is.finite(m$slope) && m$slope > 0) round(1/m$slope, 1) else NA
        ))
      }
    })

    # ── Mantel plot ────────────────────────────────────────────────────────────
    output$mantel_plot <- plotly::renderPlotly({
      m <- tryCatch(mantel_results_r(), error = function(e) NULL)
      req(m)
      r <- ibd_results_r()
      
      use_log <- (input$mantel_model == "2D")
      xv <- if (use_log) log(r$pw$dist_km) else r$pw$dist_km
      x_label <- if (use_log) "ln(geographic distance, km)" else "Geographic distance (km)"
      
      # Regression line
      x_seq <- seq(min(xv, na.rm = TRUE), max(xv, na.rm = TRUE), length.out = 100)
      y_seq <- m$intercept + m$slope * x_seq
      
      p <- plotly::plot_ly() %>%
        plotly::add_markers(
          x = xv,
          y = r$pw$FR,
          text = ~paste0(r$pw$pop_i, " \u2013 ", r$pw$pop_j,
                         "<br>Dist: ", round(r$pw$dist_km, 2), " km",
                         "<br>FR: ", round(r$pw$FR, 5)),
          hoverinfo = "text",
          marker = list(color = "#2CBF9F", size = 9, opacity = 0.85),
          name = "Pairs"
        ) %>%
        plotly::add_lines(
          x = x_seq,
          y = y_seq,
          line = list(color = "#B40F20", width = 2, dash = "solid"),
          name = paste0("Regression (b = ", round(m$slope, 6), ")")
        ) %>%
        plotly::layout(
          xaxis = list(title = x_label),
          yaxis = list(title = "FR = FST / (1 - FST)"),
          title = list(
            text = paste0("Mantel test \u2014 r = ", round(m$r, 4),
                          ", p = ", formatC(m$p_value, format = "f", digits = 4)),
            font = list(size = 13)),
          legend = list(x = 0.02, y = 0.98),
          margin = list(t = 55)
        )
      p
    })

    # ── Mantel histogram ──────────────────────────────────────────────────────
    output$mantel_hist <- plotly::renderPlotly({
      m <- tryCatch(mantel_results_r(), error = function(e) NULL)
      req(m)
      
      p <- plotly::plot_ly() %>%
        plotly::add_histogram(
          x = m$perm_r,
          nbinsx = 50,
          marker = list(color = "#94a3b8", line = list(color = "#64748b", width = 0.5)),
          name = "Permutations"
        ) %>%
        plotly::add_markers(
          x = m$r,
          y = 0,
          marker = list(color = "#B40F20", size = 12, symbol = "x"),
          name = paste0("Observed r = ", round(m$r, 4))
        ) %>%
        plotly::layout(
          xaxis = list(title = "Correlation coefficient (r)"),
          yaxis = list(title = "Frequency"),
          title = list(
            text = paste0("Permutation distribution \u2014 p = ", 
                          formatC(m$p_value, format = "f", digits = 4)),
            font = list(size = 13)),
          legend = list(x = 0.02, y = 0.98),
          margin = list(t = 55)
        )
      p
    })

    output$ui_dl_mantel <- renderUI({
      req(mantel_results_r())
      ns_local <- session$ns
      tags$div(class="ibd-dl-row",
        downloadButton(ns_local("dl_mantel"), "Download Mantel results", 
                       class = "btn btn-default btn-sm"))
    })

    output$dl_mantel <- downloadHandler(
      filename = function() paste0("Mantel_test_", Sys.Date(), ".csv"),
      content  = function(file) {
        m <- mantel_results_r()
        r <- ibd_results_r()
        df <- data.frame(
          Statistic = c("Correlation (r)", "p-value", "N pairs",
                        "Slope (b)", "Intercept", "Nb = 1/b"),
          Value = c(m$r, m$p_value, m$n_pairs,
                    m$slope, m$intercept,
                    if (is.finite(m$slope) && m$slope > 0) 1/m$slope else NA),
          stringsAsFactors = FALSE
        )
        write.csv(df, file, row.names = FALSE)
      }
    )

    # ── JavaScript to update Mantel results ───────────────────────────────────
    output$mantel_r_text <- renderUI({
      m <- tryCatch(mantel_results_r(), error = function(e) NULL)
      if (is.null(m)) return(tags$span("\u2014"))
      tags$span(style = "color:#0f172a;", round(m$r, 6))
    })
    output$mantel_p_text <- renderUI({
      m <- tryCatch(mantel_results_r(), error = function(e) NULL)
      if (is.null(m)) return(tags$span("\u2014"))
      col <- if (!is.na(m$p_value) && m$p_value < 0.05) "#166534" else "#854d0e"
      tags$span(style = paste0("color:", col, ";"), 
                if (is.na(m$p_value)) "\u2014" else formatC(m$p_value, format = "f", digits = 4))
    })
    output$mantel_n_text <- renderUI({
      m <- tryCatch(mantel_results_r(), error = function(e) NULL)
      if (is.null(m)) return(tags$span("\u2014"))
      tags$span(m$n_pairs)
    })
    output$mantel_b_text <- renderUI({
      m <- tryCatch(mantel_results_r(), error = function(e) NULL)
      if (is.null(m)) return(tags$span("\u2014"))
      tags$span(round(m$slope, 6))
    })
    output$mantel_intercept_text <- renderUI({
      m <- tryCatch(mantel_results_r(), error = function(e) NULL)
      if (is.null(m)) return(tags$span("\u2014"))
      tags$span(round(m$intercept, 6))
    })
    output$mantel_nb_text <- renderUI({
      m <- tryCatch(mantel_results_r(), error = function(e) NULL)
      if (is.null(m)) return(tags$span("\u2014"))
      nb <- if (is.finite(m$slope) && m$slope > 0) 1/m$slope else NA
      tags$span(if (is.na(nb)) "\u2014" else round(nb, 1))
    })

  })
}