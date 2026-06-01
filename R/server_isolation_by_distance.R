# server_isolation_by_distance.R
# IBD: Rousset (1997) FST/(1-FST) vs geographic distance + Mantel test.
#
# Source of truth: DuckDB tables raw (string genotypes) + meta (GPS).
# Pairwise FST via FreeNA ENA correction (Chapuis & Estoup 2007):
#   - EM null allele frequency per locus per population
#   - ENA-corrected allele frequencies (cq from EM) → WC84 FST
#   - 95% CI by individual bootstrap within populations
# Three regression lines: average (F_R), upper CI (F_R_s), lower CI (F_R_i).
# b = slope; Nb = 1/b; Nem = 1/(2*pi*b).
# Mantel test on average F_R vs ln(distance).

# ---------------------------------------------------------------------------
# File-local helpers
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# FreeNA ENA-corrected pairwise FST (Chapuis & Estoup 2007)
# ---------------------------------------------------------------------------

# EM null allele frequency — faithful translation of FreeNA Pascal source.
# Returns: list(r, cq, N, n_het, alleles_uniq)
#   r            = converged null allele frequency
#   cq           = named vector: corrected allele freqs (sum to ~1-r)
#   N            = FreeNA N (incl. null hom, excl. missing)
#   n_het        = named int vector: heterozygote counts per allele
#   alleles_uniq = character vector of unique non-null alleles
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
    if (a1 == null_code || a2 == null_code) next  # null het: skip
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

# WC84 a/b/c components for 2 populations using ENA-corrected allele freqs.
# em1, em2: return values from .ibd_em_null()
# Returns: list(a, b, c) or NULL if underdetermined
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
    # H_obs: count of heterozygotes (null hom contribute 0 het)
    h1 <- if (allele %in% names(em1$n_het)) em1$n_het[[allele]] else 0L
    h2 <- if (allele %in% names(em2$n_het)) em2$n_het[[allele]] else 0L
    hbar <- (h1 + h2) / n_total  # denominator = FreeNA N (incl null hom)
    term <- pbar * (1 - pbar)
    a_sum <- a_sum + (n_bar / n_c) *
               (s2 - (1 / (n_bar - 1)) * (term - (K-1)/K * s2 - 0.25 * hbar))
    b_sum <- b_sum + (n_bar / (n_bar - 1)) *
               (term - (K-1)/K * s2 - (2 * n_bar - 1) / (4 * n_bar) * hbar)
    c_sum <- c_sum + 0.5 * hbar
  }
  list(a = a_sum, b = b_sum, c = c_sum)
}

# FreeNA ENA-corrected pairwise FST + 95% CI (individual bootstrap).
# hap_df     : data.frame, rows = individuals, cols = loci (string "a/b" genotypes)
# pop_vector : character vector of population labels (length = nrow(hap_df))
# Bootstrap is parallelised over replicates using parallel::mclapply (ncpu-1 cores).
# Returns data.frame: pop_i, pop_j, fst, ci_l, ci_u
.pairwise_freena_fst_ci <- function(hap_df, pop_vector,
                                     null_code = "999999",
                                     n_boot = 500L, conf = 0.95) {
  pops    <- sort(unique(pop_vector))
  loci    <- colnames(hap_df)
  n_pops  <- length(pops)
  alpha   <- (1 - conf) / 2
  n_cores <- max(1L, parallel::detectCores() - 1L)

  # Compute FST for given individual index sets
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

      # Parallel bootstrap over replicates
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

.linearise <- function(x) {
  # FST / (1 - FST); clamp to [0, 0.9999] before dividing
  x <- pmin(pmax(x, 0), 0.9999)
  x / (1 - x)
}

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

  # Mantel permutation: shuffle x
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

# Fit regression of y ~ x (or log(x)) and return b, Nb, Nem
.reg_params <- function(y, x, use_log) {
  xv <- if (use_log) log(x) else x
  ok <- is.finite(y) & is.finite(xv)
  if (sum(ok) < 3L) return(c(b = NA_real_, Nb = NA_real_, Nem = NA_real_))
  cf <- coef(lm(y[ok] ~ xv[ok]))
  b  <- unname(cf[2L])
  list(
    intercept = unname(cf[1L]),
    b         = b,
    Nb        = if (is.finite(b) && b > 0) 1 / b else NA_real_,
    Nem       = if (is.finite(b) && b > 0) 1 / (2 * pi * b) else NA_real_
  )
}

# ---------------------------------------------------------------------------
# Module server
# ---------------------------------------------------------------------------

server_isolation_by_distance <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    # ── DB truth layer ─────────────────────────────────────────────────────
    `%||%`     <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
    db_tick    <- reactive({ rv$db_tick })
    con_r      <- reactive({ shiny::req(rv$con); rv$con })
    tbl_meta_r <- reactive({ rv$tbl_meta %|||% "meta" })

    db_ready <- reactive({
      db_tick()
      con <- con_r()
      shiny::req(isTRUE(rv$db_ready))
      shiny::validate(
        shiny::need(DBI::dbExistsTable(con, tbl_meta_r()), "DuckDB meta table missing.")
      )
      TRUE
    })

    # ── Population GPS centroids ───────────────────────────────────────────
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
        sql_ident(con, tbl_meta_r())))
      shiny::validate(shiny::need(nrow(df) >= 2L,
        "At least 2 populations with GPS coordinates are required."))
      df
    })

    # ── Raw string genotypes (for FreeNA ENA-corrected FST) ────────────────
    # Reads from the DuckDB raw table, reconstructing "a/b" strings per locus.
    # Returns list(hap_df, pop_vector) identical in structure to null_alleles.
    raw_genos_r <- reactive({
      db_ready()
      con     <- con_r()
      tbl_raw <- rv$tbl_raw %||% "raw"
      shiny::validate(shiny::need(
        DBI::dbExistsTable(con, tbl_raw),
        "Raw genotype table not found in DuckDB. Please re-import the dataset."))

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

      # Population assignment: always from meta (single source of truth).
      # Join raw rows (by rowid = individual) to meta Population column.
      meta_pop <- DBI::dbGetQuery(con, sprintf(
        "SELECT individual, Population FROM %s WHERE Population IS NOT NULL",
        sql_ident(con, tbl_meta_r())))
      raw_df$Population <- meta_pop$Population[match(raw_df$individual, meta_pop$individual)]
      pop_vector <- as.character(raw_df$Population)

      # Reconstruct "a/b" strings per locus
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
      list(hap_df = hap_df, pop_vector = pop_vector)
    })

    # ── Main results (triggered by Run) ────────────────────────────────────
    results_r <- eventReactive(input$run, {
      shiny::req(db_ready())
      coords  <- coords_r()
      rg      <- raw_genos_r()
      use_log <- (input$model == "2D")
      n_perm  <- as.integer(input$n_perm)
      n_boot  <- as.integer(input$n_boot_pw)

      hap_df     <- rg$hap_df
      pop_vector <- rg$pop_vector

      # Keep only populations that have GPS data
      pops_with_gps <- coords$Population
      keep_pop      <- pop_vector %in% pops_with_gps
      shiny::validate(shiny::need(
        length(unique(pop_vector[keep_pop])) >= 2L,
        "Less than 2 populations have both genotypes and GPS data."))

      hap_sub <- hap_df[keep_pop, , drop = FALSE]
      pop_sub <- pop_vector[keep_pop]
      new_levels <- sort(unique(pop_sub))

      # Pairwise FreeNA ENA-corrected FST + CI
      withProgress(message = "Computing FreeNA ENA-corrected pairwise FST...", value = 0.2, {
        pw <- .pairwise_freena_fst_ci(hap_sub, pop_sub, n_boot = n_boot)
      })

      # Pairwise geographic distances
      coords_ord <- coords[match(new_levels, coords$Population), ]
      coords_ord <- coords_ord[!is.na(coords_ord$Population), ]
      dist_mat   <- .geo_dist_matrix(coords_ord)

      # Match pairs to distance matrix
      get_dist <- function(p1, p2) {
        if (p1 %in% rownames(dist_mat) && p2 %in% rownames(dist_mat))
          dist_mat[p1, p2] else NA_real_
      }
      pw$dist_km <- mapply(get_dist, pw$pop_i, pw$pop_j)

      # Linearised FST (F_R, F_R_i, F_R_s)
      pw$FR   <- .linearise(pw$fst)
      pw$FR_i <- .linearise(pw$ci_l)
      pw$FR_s <- .linearise(pw$ci_u)

      # Three regression fits
      reg_avg <- .reg_params(pw$FR,   pw$dist_km, use_log)
      reg_ls  <- .reg_params(pw$FR_s, pw$dist_km, use_log)
      reg_li  <- .reg_params(pw$FR_i, pw$dist_km, use_log)

      # Mantel test (on average F_R)
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
        use_log  = use_log
      )
    })

    # ── Summary boxes ──────────────────────────────────────────────────────
    output$box_npops <- renderValueBox({
      r <- results_r()
      n <- length(unique(c(r$pw$pop_i, r$pw$pop_j)))
      valueBox(n, HTML("Populations<br>with GPS"), icon = icon("map-marker-alt"), color = "teal")
    })
    output$box_npairs <- renderValueBox({
      valueBox(nrow(results_r()$pw), HTML("Population<br>pairs"),
               icon = icon("project-diagram"), color = "blue")
    })
    output$box_mantel_r <- renderValueBox({
      r <- results_r()
      valueBox(round(r$mantel$r, 4), HTML("Mantel r"),
               icon = icon("chart-line"), color = "purple")
    })
    output$box_pval <- renderValueBox({
      pv  <- results_r()$mantel$p_value
      col <- if (!is.na(pv) && pv < 0.05) "green" else "yellow"
      valueBox(
        if (is.na(pv)) "NA" else formatC(pv, format = "f", digits = 4),
        HTML("Mantel p-value<br>(one-sided)"),
        icon = icon("check-circle"), color = col)
    })

    # ── Regression parameters table (b, Nb, Nem) ──────────────────────────
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

    # ── IBD plot: scatter + 3 regression lines ─────────────────────────────
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

      # Build plot
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
                         "<br>F\u209b\u209c(FreeNA): ", fst,
                         "<br>F\u209b\u209c/(1\u2212F\u209b\u209c): ", fr),
          hoverinfo = "text",
          marker = list(color = "#2CBF9F", size = 7, opacity = 0.85),
          name = "Pairs"
        ) %>%

        # Regression line: Average
        { if (!is.null(line_avg))
            plotly::add_lines(., data = line_avg, x = ~x, y = ~y,
              line = list(color = "#333a43", width = 2, dash = "solid"),
              name = paste0("Average  b=", formatC(r$reg_avg$b, format="f", digits=4)))
          else . } %>%

        # Regression line: Upper CI (ls)
        { if (!is.null(line_ls))
            plotly::add_lines(., data = line_ls, x = ~x, y = ~y,
              line = list(color = "#B40F20", width = 1.5, dash = "dash"),
              name = paste0("Upper CI (ls)  b=", formatC(r$reg_ls$b, format="f", digits=4)))
          else . } %>%

        # Regression line: Lower CI (li)
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

    # ── Pairwise FST + linearised values table ─────────────────────────────
    output$fst_table <- DT::renderDT({
      r  <- results_r()
      pw <- r$pw
      df <- data.frame(
        Pop1       = pw$pop_i,
        Pop2       = pw$pop_j,
        Dist_km    = round(pw$dist_km, 2),
        FST_FreeNA = round(pw$fst,  5),
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
                     "FST_FreeNA", "CI lower", "CI upper",
                     "FR", "FR lower", "FR upper")
      ) %>%
        DT::formatStyle("FST_FreeNA",
          backgroundColor = DT::styleInterval(
            c(0.05, 0.15, 0.25),
            c("#d4edda", "#fff3cd", "#f8d7da", "#c3002f22")))
    })

    # ── Pairwise distance matrix ───────────────────────────────────────────
    output$dist_table <- DT::renderDT({
      r  <- results_r()
      df <- as.data.frame(round(r$dist_mat, 2))
      df <- cbind(Population = rownames(df), df)
      DT::datatable(df, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15, dom = "t"),
        class   = "compact stripe hover")
    })

    # ── Downloads ──────────────────────────────────────────────────────────
    output$dl_fst <- downloadHandler(
      filename = function() paste0("IBD_pairwise_FST_", Sys.Date(), ".csv"),
      content  = function(file) {
        r  <- results_r()
        pw <- r$pw
        df <- data.frame(
          Pop1 = pw$pop_i, Pop2 = pw$pop_j,
          Dist_km = round(pw$dist_km, 4),
          FST_FreeNA = round(pw$fst, 6), CI_lower = round(pw$ci_l, 6),
          CI_upper = round(pw$ci_u, 6),
          FR = round(pw$FR, 6), FR_lower = round(pw$FR_i, 6),
          FR_upper = round(pw$FR_s, 6),
          stringsAsFactors = FALSE)
        write.csv(df, file, row.names = FALSE)
      }
    )
    output$dl_dist <- downloadHandler(
      filename = function() paste0("IBD_distances_km_", Sys.Date(), ".csv"),
      content  = function(file) {
        r  <- results_r()
        df <- as.data.frame(round(r$dist_mat, 4))
        df <- cbind(Population = rownames(df), df)
        write.csv(df, file, row.names = FALSE)
      }
    )
  })
}
