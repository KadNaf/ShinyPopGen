# server_isolation_by_distance.R
# IBD: Rousset (1997) FST/(1-FST) vs geographic distance + Mantel test.
# Bootstrap confidence intervals over loci (FreeNA approach).

# Source of truth: DuckDB tables raw (string genotypes) + meta (GPS).
# Pairwise FST via FreeNA ENA correction (Chapuis & Estoup 2007):
#   - EM null allele frequency per locus per population
#   - ENA-corrected allele frequencies (cq from EM) → WC84 FST
#   - 95% CI by individual bootstrap within populations
# Three regression lines: average (F_R), upper CI (F_R_s), lower CI (F_R_i).
# b = slope; Nb = 1/b; Nem = 1/(2*pi*b).
# Mantel test on average F_R vs ln(distance).

server_isolation_by_distance <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    # ── Helpers ────────────────────────────────────────────────────────────────
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
    safe_choice <- function(x, default = "all") {
      if (is.null(x) || length(x) == 0L || identical(x, "") || all(is.na(x))) default
      else as.character(x[[1]])
    }
    sql_id   <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))
    sql_str  <- function(con, x) as.character(DBI::dbQuoteString(con, x))

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

    db_ready <- reactive({
      db_tick(); con <- con_r()
      shiny::req(isTRUE(rv$db_ready))
      shiny::validate(
        shiny::need(DBI::dbExistsTable(con, tbl_meta_r()), "DuckDB meta table missing.")
      )
      TRUE
    })

    # ── Population GPS centroids ──────────────────────────────────────────────
    coords_r <- reactive({
      db_ready()
      con  <- con_r()
      cols <- DBI::dbGetQuery(con, sprintf(
        "SELECT column_name FROM information_schema.columns WHERE table_name = '%s'",
        tbl_meta_r()))$column_name
      shiny::validate(shiny::need(
        all(c("Latitude", "Longitude") %in% cols),
        "No GPS data found. Re-import your dataset and assign Latitude/Longitude columns."))
      df <- DBI::dbGetQuery(con, sprintf(
        "SELECT Population,
                AVG(CAST(Latitude  AS DOUBLE)) AS Latitude,
                AVG(CAST(Longitude AS DOUBLE)) AS Longitude
         FROM %s
         WHERE Population IS NOT NULL
           AND Latitude   IS NOT NULL AND Longitude IS NOT NULL
         GROUP BY Population ORDER BY Population",
        sql_id(con, tbl_meta_r())))
      shiny::validate(shiny::need(nrow(df) >= 2L,
        "At least 2 populations with GPS coordinates are required."))
      df
    })

    # ── Raw string genotypes ──────────────────────────────────────────────────
    raw_genos_r <- reactive({
      db_ready()
      con     <- con_r()
      tbl_raw <- rv$tbl_raw %||% "raw"
      shiny::validate(shiny::need(
        DBI::dbExistsTable(con, tbl_raw),
        "Raw genotype table not found in DuckDB."))

      ok_par <- tryCatch(DBI::dbExistsTable(con, "params"), error = function(e) FALSE)
      shiny::validate(shiny::need(ok_par, "params table not found."))

      marker_json <- tryCatch(
        DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='marker_cols_raw'")$value[1L],
        error = function(e) NA_character_)
      marker_cols_raw <- if (!is.na(marker_json) && nzchar(marker_json))
        tryCatch(jsonlite::fromJSON(marker_json), error = function(e) character(0))
      else character(0)
      shiny::validate(shiny::need(length(marker_cols_raw) > 0L,
        "No marker_cols_raw found in DuckDB params."))

      geno_fmt <- tryCatch(
        DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='genotype_format'")$value[1L],
        error = function(e) NA_character_)
      if (is.na(geno_fmt) || !nzchar(geno_fmt))
        geno_fmt <- if (any(grepl("(_1|\\.[0-9]+)$", marker_cols_raw))) "paired" else "string"

      keep     <- unique(marker_cols_raw)
      keep_sql <- paste(vapply(keep, function(x)
        as.character(DBI::dbQuoteIdentifier(con, x)), character(1L)), collapse = ", ")
      raw_df   <- as.data.frame(
        DBI::dbGetQuery(con, sprintf("SELECT rowid AS individual, %s FROM %s",
                                     keep_sql,
                                     as.character(DBI::dbQuoteIdentifier(con, tbl_raw)))),
        stringsAsFactors = FALSE)
      shiny::validate(shiny::need(nrow(raw_df) > 0L, "No rows in raw table."))

      meta_pop <- DBI::dbGetQuery(con, sprintf(
        "SELECT individual, Population FROM %s WHERE Population IS NOT NULL",
        sql_id(con, tbl_meta_r())))
      raw_df$Population <- meta_pop$Population[match(raw_df$individual, meta_pop$individual)]
      pop_vector <- as.character(raw_df$Population)

      pick_b <- function(locus, nms) {
        cands <- c(paste0(locus, "_1"), paste0(locus, "_2"),
                   paste0(locus, ".", 1:9), paste0(locus, "_", 1:9))
        hit   <- cands[cands %in% nms]; if (length(hit)) hit[1L] else NA_character_
      }
      if (identical(geno_fmt, "paired")) {
        nms  <- names(raw_df)
        loci <- unique(sub("(_1|_2|\\.[0-9]+)$", "", marker_cols_raw))
        hap_df <- data.frame(row.names = seq_len(nrow(raw_df)))
        for (locus in loci) {
          b <- pick_b(locus, nms)
          if (!locus %in% nms || is.na(b) || !b %in% nms) next
          a_val <- as.character(raw_df[[locus]])
          b_val <- as.character(raw_df[[b]])
          a_val[is.na(a_val) | trimws(a_val) == ""] <- "0"
          b_val[is.na(b_val) | trimws(b_val) == ""] <- "0"
          already <- grepl("/", a_val, fixed = TRUE) | grepl("-", a_val, fixed = TRUE)
          hap_df[[locus]] <- ifelse(already, a_val, paste0(a_val, "/", b_val))
        }
      } else {
        hap_df <- as.data.frame(raw_df[, marker_cols_raw, drop = FALSE],
                                 stringsAsFactors = FALSE)
        for (j in seq_along(hap_df)) {
          x <- as.character(hap_df[[j]])
          x[is.na(x) | trimws(x) == ""] <- "0/0"
          hap_df[[j]] <- x
        }
      }
      shiny::validate(shiny::need(ncol(hap_df) > 0L,
        "No locus columns could be reconstructed from raw table."))
      list(hap_df = hap_df, pop_vector = pop_vector, markers = colnames(hap_df))
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

    # ── EM null allele frequency ──────────────────────────────────────────────
    .ibd_em_null <- function(g_vec, null_code = "999999", miss_code = "0",
                              tol = 1e-6, max_iter = 10000L) {
      g <- as.character(g_vec)
      g <- g[!is.na(g) & nzchar(trimws(g))]

      parse2 <- function(gg) {
        if (grepl("/", gg, fixed = TRUE)) strsplit(gg, "/", fixed = TRUE)[[1L]]
        else if (grepl("-", gg, fixed = TRUE)) strsplit(gg, "-", fixed = TRUE)[[1L]]
        else c(gg, gg)
      }

      n_nullhom <- 0L; allele_list <- character(0); genos_obs <- list()

      for (gg in g) {
        al <- parse2(trimws(gg))
        a1 <- trimws(al[1L]); a2 <- trimws(al[2L])
        if ((a1 %in% c("0", "", miss_code)) && (a2 %in% c("0", "", miss_code))) next
        if (a1 == null_code && a2 == null_code) { n_nullhom <- n_nullhom + 1L; next }
        if (a1 == null_code || a2 == null_code) next
        allele_list <- c(allele_list, a1, a2)
        genos_obs[[length(genos_obs) + 1L]] <- c(a1, a2)
      }

      N <- n_nullhom + length(genos_obs)
      empty_r <- list(r = 0.0, cq = setNames(numeric(0), character(0)),
                      N = 0L, n_het = setNames(integer(0), character(0)),
                      alleles_uniq = character(0))
      if (N == 0L)  return(empty_r)
      if (length(genos_obs) == 0L)
        return(list(r = 1.0, cq = setNames(numeric(0), character(0)),
                    N = N, n_het = setNames(integer(0), character(0)),
                    alleles_uniq = character(0)))

      alleles_uniq <- sort(unique(allele_list)); A <- length(alleles_uniq)
      N_obs <- length(genos_obs)
      al_table <- table(factor(allele_list, levels = alleles_uniq))
      genefreq  <- as.numeric(al_table) / (2.0 * N_obs)
      names(genefreq) <- alleles_uniq

      n_hom <- setNames(integer(A), alleles_uniq)
      n_het <- setNames(integer(A), alleles_uniq)
      for (pair in genos_obs) {
        a1 <- pair[1L]; a2 <- pair[2L]
        if (a1 == a2) n_hom[a1] <- n_hom[a1] + 1L
        else { n_het[a1] <- n_het[a1] + 1L; n_het[a2] <- n_het[a2] + 1L }
      }

      rd  <- if (n_nullhom > 0L) sqrt(n_nullhom / N) else sqrt(1.0 / (N + 1L))
      cq  <- setNames(numeric(A), alleles_uniq)
      for (k in seq_len(A)) {
        a <- alleles_uniq[k]
        if (genefreq[k] <= 0) next
        ii <- n_hom[a]; jj <- n_het[a]
        cq[k] <- if (n_nullhom > 0L) 1 - sqrt((n_nullhom + N - ii - jj) / N)
                  else                1 - sqrt((1 + N - ii - jj) / (N + 1L))
        cq[k] <- min(max(cq[k], 1e-10), 1 - 1e-10)
      }

      old_cq <- cq
      for (iter in seq_len(max_iter)) {
        rdi <- 0.0
        for (k in seq_len(A)) {
          a <- alleles_uniq[k]
          if (genefreq[k] <= 0) next
          ii <- n_hom[a]; jj <- n_het[a]; cq_old <- cq[k]
          denom <- cq_old + 2 * rd
          if (denom > 0) {
            cq[k] <- ((cq_old + rd) / denom) * (ii / N) + jj / (2 * N)
            rdi    <- rdi + (rd / denom) * (ii / N)
          }
          cq[k] <- min(max(cq[k], 1e-10), 1 - 1e-10)
        }
        rd_new     <- min(max(rdi + n_nullhom / N, 1e-10), 1 - 1e-10)
        rd_change  <- abs(rd_new - rd)
        cq_change  <- if (A > 0) max(abs(cq - old_cq), na.rm = TRUE) else 0
        old_cq <- cq; rd <- rd_new
        if (rd_change < tol && cq_change < tol) break
      }

      list(r = rd, cq = cq, N = N, n_het = n_het, alleles_uniq = alleles_uniq)
    }

    # ── WC84 ENA components ───────────────────────────────────────────────────
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

    # ── Pairwise FreeNA FST + CI ─────────────────────────────────────────────
    .pairwise_freena_fst_ci <- function(hap_df, pop_vector,
                                         null_code = "999999",
                                         n_boot = 500L, conf = 0.95) {
      pops    <- sort(unique(pop_vector))
      loci    <- colnames(hap_df)
      n_pops  <- length(pops)
      alpha   <- (1 - conf) / 2
      n_cores <- max(1L, parallel::detectCores() - 1L)

      .fst_pair <- function(idx1, idx2) {
        a_tot <- 0; b_tot <- 0; c_tot <- 0; ok <- 0L
        for (locus in loci) {
          em1  <- .ibd_em_null(hap_df[[locus]][idx1], null_code)
          em2  <- .ibd_em_null(hap_df[[locus]][idx2], null_code)
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

          boot_list <- parallel::mclapply(seq_len(n_boot), function(.b) {
            b1 <- sample(idx1, length(idx1), replace = TRUE)
            b2 <- sample(idx2, length(idx2), replace = TRUE)
            .fst_pair(b1, b2)
          }, mc.cores = n_cores)
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

    # ── Bootstrap over loci ────────────────────────────────────────────────────
    .bootstrap_over_loci <- function(hap_df, pop_vector, null_code = "999999",
                                      n_boot = 1000, conf = 0.95) {
      pops    <- sort(unique(pop_vector))
      loci    <- colnames(hap_df)
      n_loci  <- length(loci)
      n_pops  <- length(pops)
      alpha   <- (1 - conf) / 2
      n_cores <- max(1L, parallel::detectCores() - 1L)

      .fst_loci <- function(idx1, idx2, locus_indices) {
        a_tot <- 0; b_tot <- 0; c_tot <- 0; ok <- 0L
        for (locus in locus_indices) {
          em1  <- .ibd_em_null(hap_df[[loci[locus]]][idx1], null_code)
          em2  <- .ibd_em_null(hap_df[[loci[locus]]][idx2], null_code)
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

          boot_list <- parallel::mclapply(seq_len(n_boot), function(.b) {
            locus_idx <- sample(seq_len(n_loci), n_loci, replace = TRUE)
            .fst_loci(idx1, idx2, locus_idx)
          }, mc.cores = n_cores)

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

    .linearise <- function(x) {
      if (is.na(x)) return(NA_real_)
      x <- pmin(pmax(x, 0), 0.9999)
      x / (1 - x)
    }

    .linearise_vec <- function(x) {
      vapply(x, .linearise, numeric(1))
    }

    # ── Mantel test ────────────────────────────────────────────────────────────
    .mantel_test <- function(fr_vec, dist_vec, n_perm = 9999, use_log = FALSE) {
      x <- if (use_log) log(dist_vec) else dist_vec
      ok <- is.finite(fr_vec) & is.finite(x)
      y  <- fr_vec[ok]; x <- x[ok]
      n  <- length(y)
      if (n < 3L)
        return(list(r = NA_real_, p_value = NA_real_, n_pairs = 0L,
                    slope = NA_real_, intercept = NA_real_))

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
           intercept = unname(coef(lm0)[1L]))
    }

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

    # ══════════════════════════════════════════════════════════════════════════
    #  MAIN IBD RESULTS
    # ══════════════════════════════════════════════════════════════════════════
    results_r <- eventReactive(input$run_ibd, {
      shiny::req(db_ready())
      coords  <- coords_r()
      rg      <- raw_genos_r()
      use_log <- (input$model == "2D")
      n_perm  <- as.integer(input$n_perm)
      n_boot  <- as.integer(input$n_boot_pw)

      hap_df     <- rg$hap_df
      pop_vector <- rg$pop_vector
      markers    <- rg$markers

      pops_with_gps <- coords$Population
      keep_pop      <- pop_vector %in% pops_with_gps
      shiny::validate(shiny::need(
        length(unique(pop_vector[keep_pop])) >= 2L,
        "Less than 2 populations have both genotypes and GPS data."))

      hap_sub <- hap_df[keep_pop, , drop = FALSE]
      pop_sub <- pop_vector[keep_pop]
      new_levels <- sort(unique(pop_sub))

      withProgress(message = "Computing FreeNA ENA-corrected pairwise FST...", value = 0.2, {
        pw <- .pairwise_freena_fst_ci(hap_sub, pop_sub, n_boot = n_boot)
      })

      coords_ord <- coords[match(new_levels, coords$Population), ]
      coords_ord <- coords_ord[!is.na(coords_ord$Population), ]
      dist_mat   <- .geo_dist_matrix(coords_ord)

      get_dist <- function(p1, p2) {
        if (p1 %in% rownames(dist_mat) && p2 %in% rownames(dist_mat))
          dist_mat[p1, p2] else NA_real_
      }
      pw$dist_km <- mapply(get_dist, pw$pop_i, pw$pop_j)

      pw$FR   <- .linearise_vec(pw$fst)
      pw$FR_i <- .linearise_vec(pw$ci_l)
      pw$FR_s <- .linearise_vec(pw$ci_u)

      reg_avg <- .reg_params(pw$FR,   pw$dist_km, use_log)
      reg_ls  <- .reg_params(pw$FR_s, pw$dist_km, use_log)
      reg_li  <- .reg_params(pw$FR_i, pw$dist_km, use_log)

      withProgress(message = "Running Mantel test...", value = 0.8, {
        mantel <- .mantel_test(pw$FR, pw$dist_km, n_perm = n_perm, use_log = use_log)
      })

      list(
        pw       = pw,
        dist_mat = dist_mat,
        reg_avg  = reg_avg,
        reg_ls   = reg_ls,
        reg_li   = reg_li,
        mantel   = mantel,
        use_log  = use_log,
        hap_df   = hap_sub,
        pop_vector = pop_sub,
        pops     = new_levels,
        markers  = markers,
        n_loci   = length(markers)
      )
    })

    # ── Bootstrap results ──────────────────────────────────────────────────────
    boot_results_r <- eventReactive(input$run_boot, {
      shiny::req(results_r())
      res <- results_r()
      n_boot <- as.integer(input$n_boot_loci)
      alpha  <- as.numeric(input$boot_ci_level %||% "0.05")

      withProgress(message = "Bootstrapping over loci...", value = 0.3, {
        boot <- .bootstrap_over_loci(res$hap_df, res$pop_vector, 
                                      n_boot = n_boot, conf = 1 - alpha)
      })

      dist_mat <- res$dist_mat
      get_dist <- function(p1, p2) {
        if (p1 %in% rownames(dist_mat) && p2 %in% rownames(dist_mat))
          dist_mat[p1, p2] else NA_real_
      }
      boot$dist_km <- mapply(get_dist, boot$pop_i, boot$pop_j)

      list(
        boot    = boot,
        n_boot  = n_boot,
        alpha   = alpha,
        ci      = ci_bounds(alpha)
      )
    })

    # ── Value boxes ────────────────────────────────────────────────────────────
    output$vb_pops <- renderUI({
      tryCatch({
        r <- results_r()
        n <- length(unique(c(r$pw$pop_i, r$pw$pop_j)))
        tags$span(n)
      }, error = function(e) tags$span("\u2014"))
    })

    output$vb_pairs <- renderUI({
      tryCatch({
        r <- results_r()
        tags$span(nrow(r$pw))
      }, error = function(e) tags$span("\u2014"))
    })

    output$vb_mantel_r <- renderUI({
      tryCatch({
        r <- results_r()
        v <- round(r$mantel$r, 4)
        tags$span(if (is.na(v)) "\u2014" else v)
      }, error = function(e) tags$span("\u2014"))
    })

    output$vb_pval <- renderUI({
      tryCatch({
        r <- results_r()
        v <- r$mantel$p_value
        col <- if (!is.na(v) && v < 0.05) "#166534" else "#854d0e"
        tags$span(style = paste0("color:", col, ";"),
                  if (is.na(v)) "\u2014" else formatC(v, format = "f", digits = 4))
      }, error = function(e) tags$span("\u2014"))
    })

    output$vb_nb <- renderUI({
      tryCatch({
        r <- results_r()
        v <- r$reg_avg$Nb
        tags$span(if (is.na(v)) "\u2014" else round(v, 1))
      }, error = function(e) tags$span("\u2014"))
    })

    # ── Regression table ───────────────────────────────────────────────────────
    output$reg_table <- DT::renderDT({
      r <- results_r()
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
      r <- tryCatch(results_r(), error = function(e) NULL)
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
      r      <- results_r()
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
        plotly::add_segments(
          data = pw,
          x = ~xv, xend = ~xv, y = ~FR_i, yend = ~FR_s,
          line = list(color = "rgba(100,100,100,0.35)", width = 1),
          showlegend = FALSE, hoverinfo = "none"
        ) %>%
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

    # ── FST table ──────────────────────────────────────────────────────────────
    output$fst_table <- DT::renderDT({
      r  <- results_r()
      pw <- r$pw
      df <- data.frame(
        Pop1       = pw$pop_i,
        Pop2       = pw$pop_j,
        Dist_km    = round(pw$dist_km, 2),
        FST_ENA    = round(pw$fst,  5),
        CI_lower   = round(pw$ci_l, 5),
        CI_upper   = round(pw$ci_u, 5),
        FR         = round(pw$FR,   5),
        FR_lower   = round(pw$FR_i, 5),
        FR_upper   = round(pw$FR_s, 5),
        stringsAsFactors = FALSE
      )
      DT::datatable(df, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15,
                       dom = "lrtip"),
        class   = "compact stripe hover",
        colnames = c("Pop 1", "Pop 2", "Dist (km)",
                     "FST-ENA", "CI lower", "CI upper",
                     "FR", "FR lower", "FR upper")
      ) %>%
        DT::formatStyle("FST_ENA",
          backgroundColor = DT::styleInterval(
            c(0.05, 0.15, 0.25),
            c("#d4edda", "#fff3cd", "#f8d7da", "#c3002f22")))
    })

    output$ui_dl_fst <- renderUI({
      req(results_r())
      ns_local <- session$ns
      tags$div(class="ibd-dl-row",
        downloadButton(ns_local("dl_fst"), "Download table", class = "btn btn-default btn-sm"))
    })

    output$dl_fst <- downloadHandler(
      filename = function() paste0("IBD_pairwise_FST_", Sys.Date(), ".csv"),
      content  = function(file) {
        r  <- results_r()
        pw <- r$pw
        df <- data.frame(
          Pop1 = pw$pop_i, Pop2 = pw$pop_j,
          Dist_km = round(pw$dist_km, 4),
          FST_ENA = round(pw$fst, 6), CI_lower = round(pw$ci_l, 6),
          CI_upper = round(pw$ci_u, 6),
          FR = round(pw$FR, 6), FR_lower = round(pw$FR_i, 6),
          FR_upper = round(pw$FR_s, 6),
          stringsAsFactors = FALSE)
        write.csv(df, file, row.names = FALSE)
      }
    )

    # ── Distance table ─────────────────────────────────────────────────────────
    output$dist_table <- DT::renderDT({
      r  <- results_r()
      df <- as.data.frame(round(r$dist_mat, 2))
      df <- cbind(Population = rownames(df), df)
      DT::datatable(df, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15, dom = "t"),
        class   = "compact stripe hover")
    })

    output$ui_dl_dist <- renderUI({
      req(results_r())
      ns_local <- session$ns
      tags$div(class="ibd-dl-row",
        downloadButton(ns_local("dl_dist"), "Download distances", class = "btn btn-default btn-sm"))
    })

    output$dl_dist <- downloadHandler(
      filename = function() paste0("IBD_distances_km_", Sys.Date(), ".csv"),
      content  = function(file) {
        r  <- results_r()
        df <- as.data.frame(round(r$dist_mat, 4))
        df <- cbind(Population = rownames(df), df)
        write.csv(df, file, row.names = FALSE)
      }
    )

    # ── Bootstrap summary boxes ───────────────────────────────────────────────
    output$boot_n_loci <- renderUI({
      r <- tryCatch(results_r(), error = function(e) NULL)
      if (is.null(r)) tags$span("\u2014") else tags$span(r$n_loci)
    })

    output$boot_n_reps <- renderUI({
      b <- tryCatch(boot_results_r(), error = function(e) NULL)
      if (is.null(b)) tags$span("\u2014") else tags$span(b$n_boot)
    })

    output$boot_n_valid <- renderUI({
      b <- tryCatch(boot_results_r(), error = function(e) NULL)
      if (is.null(b)) tags$span("\u2014")
      else tags$span(sum(!is.na(b$boot$fst)))
    })

    output$ui_boot_status <- renderUI({
      b <- tryCatch(boot_results_r(), error = function(e) NULL)
      if (is.null(b)) return(NULL)
      tags$div(class="ibd-boot-result",
        tags$strong(sprintf("Bootstrap complete — %d replicates, %s CI", 
                           b$n_boot, b$ci$label)),
        tags$br(),
        sprintf("%d pairwise estimates with valid CI", sum(!is.na(b$boot$fst_ci_l)))
      )
    })

    # ── Bootstrap table ──────────────────────────────────────────────────────
    output$boot_table <- DT::renderDT({
      b <- boot_results_r()
      boot <- b$boot
      df <- data.frame(
        Pop1 = boot$pop_i,
        Pop2 = boot$pop_j,
        Dist_km = round(boot$dist_km, 2),
        FST_ENA = round(boot$fst, 5),
        FST_CI_low = round(boot$fst_ci_l, 5),
        FST_CI_high = round(boot$fst_ci_u, 5),
        FR = round(boot$fr, 5),
        FR_CI_low = round(boot$fr_ci_l, 5),
        FR_CI_high = round(boot$fr_ci_u, 5),
        N_loci = boot$n_loci_used,
        stringsAsFactors = FALSE
      )
      DT::datatable(df, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15,
                       dom = "lrtip"),
        class = "compact stripe hover",
        colnames = c("Pop 1", "Pop 2", "Dist (km)",
                     "FST-ENA", "FST CI low", "FST CI high",
                     "FR", "FR CI low", "FR CI high",
                     "N loci")
      ) %>%
        DT::formatStyle("FST_ENA",
          backgroundColor = DT::styleInterval(
            c(0.05, 0.15, 0.25),
            c("#d4edda", "#fff3cd", "#f8d7da", "#c3002f22")))
    })

    output$ui_dl_boot <- renderUI({
      req(boot_results_r())
      ns_local <- session$ns
      tags$div(class="ibd-dl-row",
        downloadButton(ns_local("dl_boot"), "Download bootstrap results", 
                       class = "btn btn-default btn-sm"))
    })

    output$dl_boot <- downloadHandler(
      filename = function() paste0("IBD_bootstrap_CI_", Sys.Date(), ".csv"),
      content  = function(file) {
        b <- boot_results_r()
        boot <- b$boot
        df <- data.frame(
          Pop1 = boot$pop_i, Pop2 = boot$pop_j,
          Dist_km = round(boot$dist_km, 4),
          FST_ENA = round(boot$fst, 6), FST_CI_low = round(boot$fst_ci_l, 6),
          FST_CI_high = round(boot$fst_ci_u, 6),
          FR = round(boot$fr, 6), FR_CI_low = round(boot$fr_ci_l, 6),
          FR_CI_high = round(boot$fr_ci_u, 6),
          N_loci = boot$n_loci_used,
          stringsAsFactors = FALSE)
        write.csv(df, file, row.names = FALSE)
      }
    )

    # ── Bootstrap FST plot ───────────────────────────────────────────────────
    output$boot_fst_plot <- plotly::renderPlotly({
      b <- boot_results_r()
      boot <- b$boot
      
      boot <- boot[order(boot$dist_km), ]
      boot$pair <- paste0(boot$pop_i, "–", boot$pop_j)
      
      p <- plotly::plot_ly() %>%
        plotly::add_markers(
          data = boot,
          x = ~pair,
          y = ~fst,
          marker = list(color = "#2CBF9F", size = 8),
          name = "FST-ENA",
          hoverinfo = "text",
          text = ~paste0(pair, "<br>Dist: ", round(dist_km, 1), " km<br>FST: ", round(fst, 4))
        ) %>%
        plotly::add_segments(
          data = boot,
          x = ~pair, xend = ~pair,
          y = ~fst_ci_l, yend = ~fst_ci_u,
          line = list(color = "rgba(44,191,159,0.5)", width = 2),
          showlegend = FALSE,
          hoverinfo = "none"
        ) %>%
        plotly::layout(
          xaxis = list(title = "Population pair", tickangle = -45),
          yaxis = list(title = "FST-ENA (FreeNA corrected)"),
          title = paste0("Bootstrap ", b$ci$label, " for FST-ENA"),
          margin = list(b = 100)
        )
      p
    })

    # ── Bootstrap FR plot ───────────────────────────────────────────────────
    output$boot_fr_plot <- plotly::renderPlotly({
      b <- boot_results_r()
      boot <- b$boot
      
      boot <- boot[order(boot$dist_km), ]
      boot$pair <- paste0(boot$pop_i, "–", boot$pop_j)
      
      p <- plotly::plot_ly() %>%
        plotly::add_markers(
          data = boot,
          x = ~pair,
          y = ~fr,
          marker = list(color = "#3B9AB2", size = 8),
          name = "FR",
          hoverinfo = "text",
          text = ~paste0(pair, "<br>Dist: ", round(dist_km, 1), " km<br>FR: ", round(fr, 4))
        ) %>%
        plotly::add_segments(
          data = boot,
          x = ~pair, xend = ~pair,
          y = ~fr_ci_l, yend = ~fr_ci_u,
          line = list(color = "rgba(59,154,178,0.5)", width = 2),
          showlegend = FALSE,
          hoverinfo = "none"
        ) %>%
        plotly::layout(
          xaxis = list(title = "Population pair", tickangle = -45),
          yaxis = list(title = "FR = FST / (1 - FST)"),
          title = paste0("Bootstrap ", b$ci$label, " for FR"),
          margin = list(b = 100)
        )
      p
    })

  })
}