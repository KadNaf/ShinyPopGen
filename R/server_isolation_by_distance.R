# server_isolation_by_distance.R
# Faithful R translation of FreeNA_optm2R (Chapuis & Estoup 2007):
#   - EM Dempster null allele frequency
#   - WC84 FST (global + pairwise), with ENA correction
#   - Cavalli-Sforza & Edwards distance, with INA correction
#   - Bootstrap resampling over loci for 95% CI
#   - Mantel test (square or rectangular matrices)

# ---------------------------------------------------------------------------
# File-local helpers — FreeNA algorithms translated from Pascal
# ---------------------------------------------------------------------------

# EM Dempster et al. (1977) null allele frequency per locus per population.
# g_vec: character vector of "a/b" genotypes
# Returns: list(rd = null allele freq, cq = corrected allele freqs (named),
#               nnullhomo = count of null homozygotes, absentgeno = missing,
#               n_eff = effective N excluding null homo + missing,
#               alleles = unique non-null alleles, n_het, n_hom per allele)
.freena_em_null <- function(g_vec, null_code = "999999", miss_code = "0",
                             tol = 1e-6, max_iter = 1000L) {
  g <- as.character(g_vec)
  g <- g[!is.na(g) & nzchar(trimws(g)) & g != "0/0" & g != "0"]

  # Parse alleles
  allele_list <- character(0)
  n_nullhomo <- 0L
  absentgeno <- 0L
  genos_obs <- list()

  for (gg in g) {
    parts <- trimws(strsplit(gg, "/", fixed = TRUE)[[1L]])
    if (length(parts) < 2) { absentgeno <- absentgeno + 1L; next }
    a1 <- parts[1L]; a2 <- parts[2L]
    if ((a1 %in% c("0", "", miss_code)) && (a2 %in% c("0", "", miss_code))) {
      absentgeno <- absentgeno + 1L; next
    }
    if (a1 == null_code && a2 == null_code) { n_nullhomo <- n_nullhomo + 1L; next }
    if (a1 == null_code || a2 == null_code) next  # null het: skip (FreeNA)
    allele_list <- c(allele_list, a1, a2)
    genos_obs[[length(genos_obs) + 1L]] <- c(a1, a2)
  }

  n_eff <- n_nullhomo + length(genos_obs)  # FreeNA N
  if (n_eff == 0L) {
    return(list(rd = 0, cq = numeric(0), nnullhomo = 0L,
                absentgeno = absentgeno, n_eff = 0L,
                alleles = character(0), n_het = integer(0), n_hom = integer(0)))
  }

  alleles <- sort(unique(allele_list))
  A <- length(alleles)
  n_het <- setNames(integer(A), alleles)
  n_hom <- setNames(integer(A), alleles)
  for (pair in genos_obs) {
    if (pair[1L] == pair[2L]) n_hom[pair[1L]] <- n_hom[pair[1L]] + 1L
    else { n_het[pair[1L]] <- n_het[pair[1L]] + 1L; n_het[pair[2L]] <- n_het[pair[2L]] + 1L }
  }

  # Initialization of rd (null allele frequency)
  rd <- if (n_nullhomo > 0L) sqrt(n_nullhomo / n_eff) else sqrt(1 / (n_eff + 1))

  # Initial corrected allele frequencies cq
  cq <- setNames(numeric(A), alleles)
  for (k in seq_len(A)) {
    ii <- n_hom[alleles[k]]; jj <- n_het[alleles[k]]
    hotot <- sum(n_hom)
    cq[k] <- if (n_nullhomo > 0L)
      1 - sqrt((n_nullhomo + hotot - ii + (n_eff - n_nullhomo - hotot) - jj) / n_eff)
    else
      1 - sqrt((1 + hotot - ii + (n_eff - hotot) - jj) / (n_eff + 1))
    cq[k] <- min(max(cq[k], 1e-10), 1 - 1e-10)
  }

  # EM iterations
  for (iter in seq_len(max_iter)) {
    rdi <- 0
    cq_new <- cq
    for (k in seq_len(A)) {
      ii <- n_hom[alleles[k]]; jj <- n_het[alleles[k]]
      denom <- cq[k] + 2 * rd
      if (denom > 0) {
        cq_new[k] <- ((cq[k] + rd) / denom) * (ii / n_eff) + jj / (2 * n_eff)
        rdi <- rdi + (rd / denom) * (ii / n_eff)
      }
      cq_new[k] <- min(max(cq_new[k], 1e-10), 1 - 1e-10)
    }
    rd_new <- rdi + n_nullhomo / n_eff
    if (abs(rd_new - rd) < tol && max(abs(cq_new - cq)) < tol) break
    rd <- rd_new
    cq <- cq_new
  }

  list(rd = rd, cq = cq, nnullhomo = n_nullhomo,
       absentgeno = absentgeno, n_eff = n_eff,
       alleles = alleles, n_het = n_het, n_hom = n_hom)
}

# WC84 FST components per locus (Genepop method).
# em_list: list of .freena_em_null() results per population
# correct: if TRUE, apply ENA correction (uses corrdgenefreq + cAA correction)
# Returns: list(s1l, s3l) — components for this locus
.freena_wc84_locus <- function(em_list, correct = FALSE) {
  npop <- length(em_list)
  if (npop < 2L) return(list(s1l = 0, s3l = 0, npopeff = 0L, ntoteff = 0, ntoteff2 = 0))

  # Collect all alleles
  all_alleles <- unique(unlist(lapply(em_list, function(e) e$alleles)))

  # Effective sample sizes
  ni <- sapply(em_list, function(e) {
    if (correct) e$n_eff + e$nnullhomo - e$absentgeno  # = efpop - absentgeno
    else          e$n_eff                                # = efpop - absentgeno - nnullhomo
  })
  ntoteff  <- sum(ni)
  ntoteff2 <- sum(ni^2)
  npopeff  <- sum(ni > 0)
  nc <- if (ntoteff > 0 && npopeff > 1L)
          (ntoteff - ntoteff2 / ntoteff) / (npopeff - 1L)
        else 0

  s1l <- 0; s3l <- 0

  for (allele in all_alleles) {
    # Allele count nA and observed homozygotes AA
    nA  <- numeric(npop)
    AA  <- numeric(npop)
    for (p in seq_len(npop)) {
      e <- em_list[[p]]
      nn <- 2 * ni[p]
      if (correct) {
        freq_a <- if (allele %in% names(e$cq)) e$cq[[allele]] else 0
        nA[p] <- freq_a * nn
      } else {
        # uncorrected frequency
        n_obs <- e$n_eff - e$nnullhomo  # individuals with non-null genotypes
        freq_a <- if (n_obs > 0 && allele %in% names(e$n_het) && allele %in% names(e$n_hom)) {
          (2 * e$n_hom[[allele]] + e$n_het[[allele]]) / (2 * n_obs)
        } else 0
        nA[p] <- freq_a * nn
      }
      AA[p] <- if (allele %in% names(e$n_hom)) e$n_hom[[allele]] else 0L
      # ENA correction: adjust AA
      if (correct && AA[p] > 0) {
        cq_a <- if (allele %in% names(e$cq)) e$cq[[allele]] else 0
        AA[p] <- AA[p] * (cq_a / (cq_a + 2 * e$rd))
      }
    }

    snA <- sum(nA)
    sAA <- sum(AA)
    s2A <- sum(ifelse(ni > 0, nA^2 / (2 * ni), 0))

    if (ntoteff * nc > 0 && ntoteff > npopeff) {
      MSG <- (0.5 * snA - sAA) / ntoteff
      MSI <- (0.5 * snA + sAA - s2A) / (ntoteff - npopeff)
      MSP <- (s2A - 0.5 * snA^2 / ntoteff) / (npopeff - 1L)
      s2G <- MSG
      s2I <- 0.5 * (MSI - MSG)
      s2P <- (MSP - MSI) / (2 * nc)
      s1l <- s1l + s2P
      s3l <- s3l + s2P + s2I + s2G
    }
  }

  list(s1l = s1l, s3l = s3l, npopeff = npopeff,
       ntoteff = ntoteff, ntoteff2 = ntoteff2)
}

# Cavalli-Sforza & Edwards chord distance per locus.
# em_list: list of .freena_em_null() per population
# correct: if TRUE, use INA correction (include rd as null allele)
.freena_cs_locus <- function(em_list, correct = FALSE) {
  npop <- length(em_list)
  if (npop < 2L) return(list(Dc_mat = matrix(0, npop, npop), nloceff = 0L))

  # Allele freqs per population
  if (correct) {
    freqs <- lapply(em_list, function(e) {
      c(e$cq, `__null__` = e$rd)
    })
  } else {
    freqs <- lapply(em_list, function(e) {
      # uncorrected freqs
      n_obs <- e$n_eff - e$nnullhomo
      if (n_obs <= 0) return(numeric(0))
      sapply(e$alleles, function(a) {
        (2 * e$n_hom[[a]] + e$n_het[[a]]) / (2 * n_obs)
      })
    })
  }

  Dc_mat <- matrix(0, npop, npop)
  nloceff <- 0L

  for (i in seq_len(npop - 1L)) {
    for (j in (i + 1L):npop) {
      fi <- freqs[[i]]; fj <- freqs[[j]]
      # Check non-empty populations
      all_alleles <- union(names(fi), names(fj))
      if (length(all_alleles) == 0L) next

      pi <- sapply(all_alleles, function(a) if (a %in% names(fi)) fi[[a]] else 0)
      pj <- sapply(all_alleles, function(a) if (a %in% names(fj)) fj[[a]] else 0)

      # Normalize
      si <- sum(pi); sj <- sum(pj)
      if (si <= 0 || sj <= 0) next
      pi <- pi / si
      pj <- pj / sj

      csprod <- sum(sqrt(pi * pj))
      if (csprod >= 1) next  # not applicable

      Dc <- (2 / pi) * sqrt(2 * (1 - csprod))
      Dc_mat[i, j] <- Dc_mat[j, i] <- Dc
      nloceff <- nloceff + 1L
    }
  }

  list(Dc_mat = Dc_mat, nloceff = nloceff)
}

# ---------------------------------------------------------------------------
# Full FreeNA pipeline (multi-locus)
# ---------------------------------------------------------------------------

# Main pipeline: computes null allele freqs + pairwise FST + CS + bootstrap CI
# hap_df: data.frame with loci as columns ("a/b" strings)
# pop_vector: population assignment per individual
# Returns a large list with matrices, bootstrap results, log
.freena_pipeline <- function(hap_df, pop_vector,
                              n_boot = 1000L,
                              calc_fst = TRUE, ena_corr = TRUE,
                              calc_cs = TRUE,  ina_corr = TRUE) {
  pops <- sort(unique(pop_vector))
  npop <- length(pops)
  loci <- colnames(hap_df)
  nloc <- length(loci)

  # Storage for per-locus results
  rd_mat <- matrix(NA_real_, nloc, npop,
                    dimnames = list(loci, pops))

  # Pairwise matrices (final across loci)
  fst_mat     <- matrix(0, npop, npop, dimnames = list(pops, pops))
  fst_ena_mat <- matrix(0, npop, npop, dimnames = list(pops, pops))
  cs_mat      <- matrix(0, npop, npop, dimnames = list(pops, pops))
  cs_ina_mat  <- matrix(0, npop, npop, dimnames = list(pops, pops))

  # Per-locus components storage (for bootstrap)
  s1l_fst <- s3l_fst <- array(0, dim = c(nloc, npop, npop),
                                dimnames = list(loci, pops, pops))
  s1l_ena <- s3l_ena <- array(0, dim = c(nloc, npop, npop),
                                dimnames = list(loci, pops, pops))
  Dc_loci <- array(0, dim = c(nloc, npop, npop),
                    dimnames = list(loci, pops, pops))
  Dc_ina_loci <- array(0, dim = c(nloc, npop, npop),
                        dimnames = list(loci, pops, pops))
  nloceff_fst  <- matrix(nloc, npop, npop, dimnames = list(pops, pops))
  nloceff_ena  <- matrix(nloc, npop, npop, dimnames = list(pops, pops))
  nloceff_cs   <- matrix(nloc, npop, npop, dimnames = list(pops, pops))
  nloceff_ina  <- matrix(nloc, npop, npop, dimnames = list(pops, pops))

  # Log buffer
  log_lines <- character(0)
  add_log <- function(...) log_lines <<- c(log_lines, paste0(...))

  add_log("=== FreeNA pipeline ===")
  add_log(sprintf("Individuals: %d, Populations: %d, Loci: %d",
                  nrow(hap_df), npop, nloc))

  # --- Per locus computation ---
  for (loc_idx in seq_len(nloc)) {
    locus <- loci[loc_idx]
    # EM per population
    em_list <- lapply(pops, function(p) {
      .freena_em_null(hap_df[[locus]][pop_vector == p])
    })
    names(em_list) <- pops
    for (p in pops) rd_mat[locus, p] <- em_list[[p]]$rd

    # Pairwise FST components per locus
    if (calc_fst) {
      for (i in seq_len(npop - 1L)) {
        for (j in (i + 1L):npop) {
          # Uncorrected FST (WC84)
          em_ij <- em_list[c(i, j)]
          comp <- .freena_wc84_locus(em_ij, correct = FALSE)
          s1l_fst[loc_idx, i, j] <- s1l_fst[loc_idx, j, i] <- comp$s1l
          s3l_fst[loc_idx, i, j] <- s3l_fst[loc_idx, j, i] <- comp$s3l

          # ENA-corrected FST
          if (ena_corr) {
            comp_ena <- .freena_wc84_locus(em_ij, correct = TRUE)
            s1l_ena[loc_idx, i, j] <- s1l_ena[loc_idx, j, i] <- comp_ena$s1l
            s3l_ena[loc_idx, i, j] <- s3l_ena[loc_idx, j, i] <- comp_ena$s3l
          }
        }
      }
    }

    # Cavalli-Sforza distance per locus
    if (calc_cs) {
      cs_res <- .freena_cs_locus(em_list, correct = FALSE)
      Dc_loci[loc_idx, , ] <- cs_res$Dc_mat
      for (i in seq_len(npop - 1L)) {
        for (j in (i + 1L):npop) {
          if (cs_res$Dc_mat[i, j] == 0 && (i != j)) {
            nloceff_cs[i, j] <- nloceff_cs[j, i] <- nloceff_cs[i, j] - 1L
          }
        }
      }

      if (ina_corr) {
        cs_ina <- .freena_cs_locus(em_list, correct = TRUE)
        Dc_ina_loci[loc_idx, , ] <- cs_ina$Dc_mat
        for (i in seq_len(npop - 1L)) {
          for (j in (i + 1L):npop) {
            if (cs_ina$Dc_mat[i, j] == 0 && (i != j)) {
              nloceff_ina[i, j] <- nloceff_ina[j, i] <- nloceff_ina[i, j] - 1L
            }
          }
        }
      }
    }
  }

  # --- Aggregate across loci ---
  # Global pairwise FST: sum of s1l / sum of s3l over loci
  for (i in seq_len(npop - 1L)) {
    for (j in (i + 1L):npop) {
      s1 <- sum(s1l_fst[, i, j]); s3 <- sum(s3l_fst[, i, j])
      fst_mat[i, j] <- fst_mat[j, i] <- if (s3 > 0) s1 / s3 else NA_real_

      s1 <- sum(s1l_ena[, i, j]); s3 <- sum(s3l_ena[, i, j])
      fst_ena_mat[i, j] <- fst_ena_mat[j, i] <- if (s3 > 0) s1 / s3 else NA_real_

      cs_mat[i, j]     <- cs_mat[j, i]     <- if (nloceff_cs[i, j] > 0)
        sum(Dc_loci[, i, j]) / nloceff_cs[i, j] else NA_real_
      cs_ina_mat[i, j] <- cs_ina_mat[j, i] <- if (nloceff_ina[i, j] > 0)
        sum(Dc_ina_loci[, i, j]) / nloceff_ina[i, j] else NA_real_
    }
  }

  # --- Bootstrap over loci ---
  boot_fst <- boot_fst_ena <- boot_cs <- boot_cs_ina <-
    array(NA_real_, dim = c(npop, npop, n_boot),
          dimnames = list(pops, pops, NULL))

  if (n_boot >= 100L && nloc >= 5L) {
    add_log(sprintf("Bootstrap: %d replicates over %d loci", n_boot, nloc))

    for (b in seq_len(n_boot)) {
      sel <- sample(nloc, nloc, replace = TRUE)

      for (i in seq_len(npop - 1L)) {
        for (j in (i + 1L):npop) {
          if (calc_fst) {
            s1 <- sum(s1l_fst[sel, i, j]); s3 <- sum(s3l_fst[sel, i, j])
            boot_fst[i, j, b] <- boot_fst[j, i, b] <- if (s3 > 0) s1 / s3 else NA

            if (ena_corr) {
              s1 <- sum(s1l_ena[sel, i, j]); s3 <- sum(s3l_ena[sel, i, j])
              boot_fst_ena[i, j, b] <- boot_fst_ena[j, i, b] <- if (s3 > 0) s1 / s3 else NA
            }
          }
          if (calc_cs) {
            eff <- nloceff_cs[i, j]
            boot_cs[i, j, b] <- boot_cs[j, i, b] <- if (eff > 0)
              sum(Dc_loci[sel, i, j]) / eff else NA
            if (ina_corr) {
              eff <- nloceff_ina[i, j]
              boot_cs_ina[i, j, b] <- boot_cs_ina[j, i, b] <- if (eff > 0)
                sum(Dc_ina_loci[sel, i, j]) / eff else NA
            }
          }
        }
      }
    }
  } else {
    add_log("Skipping bootstrap: need nloc >= 5 and n_boot >= 100")
  }

  # --- Compute 95% CI from bootstrap ---
  ci_lo <- function(arr) apply(arr, c(1, 2),
                                function(x) quantile(x, 0.025, na.rm = TRUE))
  ci_hi <- function(arr) apply(arr, c(1, 2),
                                function(x) quantile(x, 0.975, na.rm = TRUE))

  list(
    rd_mat       = rd_mat,
    fst_mat      = fst_mat,
    fst_ena_mat  = fst_ena_mat,
    cs_mat       = cs_mat,
    cs_ina_mat   = cs_ina_mat,
    boot_fst     = boot_fst,
    boot_fst_ena = boot_fst_ena,
    boot_cs      = boot_cs,
    boot_cs_ina  = boot_cs_ina,
    fst_ci_l     = if (!all(is.na(boot_fst)))     ci_lo(boot_fst)     else NULL,
    fst_ci_u     = if (!all(is.na(boot_fst)))     ci_hi(boot_fst)     else NULL,
    fst_ena_ci_l = if (!all(is.na(boot_fst_ena))) ci_lo(boot_fst_ena) else NULL,
    fst_ena_ci_u = if (!all(is.na(boot_fst_ena))) ci_hi(boot_fst_ena) else NULL,
    cs_ci_l      = if (!all(is.na(boot_cs)))      ci_lo(boot_cs)      else NULL,
    cs_ci_u      = if (!all(is.na(boot_cs)))      ci_hi(boot_cs)      else NULL,
    cs_ina_ci_l  = if (!all(is.na(boot_cs_ina)))  ci_lo(boot_cs_ina)  else NULL,
    cs_ina_ci_u  = if (!all(is.na(boot_cs_ina)))  ci_hi(boot_cs_ina)  else NULL,
    pops = pops, loci = loci, nloc = nloc, n_boot = n_boot,
    log = log_lines
  )
}

# ---------------------------------------------------------------------------
# Haversine + geographic distance
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
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      d <- .haversine_km(coords$Latitude[i], coords$Longitude[i],
                         coords$Latitude[j], coords$Longitude[j])
      mat[i, j] <- mat[j, i] <- d
    }
  }
  mat
}

# ---------------------------------------------------------------------------
# Mantel test (supports square or rectangular column-wise format)
# ---------------------------------------------------------------------------

# Parse uploaded file into a symmetric distance matrix.
# If format == "square": standard N×N matrix
# If format == "rectangular": first row = pair labels (PopA-PopB), rows = values;
#   average across rows to get one value per pair, then rebuild square matrix.
.parse_dist_file <- function(file_path, format = "square") {
  raw <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)

  if (format == "square") {
    # First column = row names, rest = matrix
    rn <- as.character(raw[[1L]])
    mat <- as.matrix(raw[, -1L, drop = FALSE])
    rownames(mat) <- rn
    colnames(mat) <- rn
    return(mat)
  }

  # Rectangular: columns are pairs
  # Each column header = "PopA-PopB", rows = multiple distance values
  pair_names <- names(raw)
  pair_names <- pair_names[pair_names != "" & !is.na(pair_names)]

  pops <- unique(unlist(strsplit(pair_names, "-", fixed = FALSE)))
  pops <- pops[pops != "" & !is.na(pops)]

  # Average values per column
  avg_vals <- sapply(pair_names, function(cn) {
    v <- raw[[cn]]
    v <- as.numeric(v)
    mean(v, na.rm = TRUE)
  })
  names(avg_vals) <- pair_names

  # Build square matrix
  n <- length(pops)
  mat <- matrix(NA_real_, n, n, dimnames = list(pops, pops))
  for (pair in pair_names) {
    parts <- strsplit(pair, "-")[[1L]]
    if (length(parts) == 2L) {
      p1 <- parts[1L]; p2 <- parts[2L]
      if (p1 %in% pops && p2 %in% pops) {
        mat[p1, p2] <- mat[p2, p1] <- avg_vals[[pair]]
      }
    }
  }
  diag(mat) <- 0
  mat
}

# Mantel permutation test
.mantel_test <- function(mat1, mat2, n_perm = 9999, method = "pearson") {
  # Extract lower triangles
  ok <- !is.na(mat1) & !is.na(mat2) & lower.tri(mat1)
  v1 <- mat1[ok]; v2 <- mat2[ok]
  n <- length(v1)

  if (n < 3L) {
    return(list(r = NA_real_, p_value = NA_real_, n = n,
                v1 = v1, v2 = v2, perm_r = numeric(0),
                message = "Not enough data points"))
  }

  r_obs <- cor(v1, v2, method = method)

  # Permute rows/columns of mat2 (jointly to preserve symmetry)
  perm_r <- numeric(n_perm)
  for (k in seq_len(n_perm)) {
    idx <- sample(nrow(mat2))
    perm_mat2 <- mat2[idx, idx]
    v2_perm <- perm_mat2[ok]
    perm_r[k] <- cor(v1, v2_perm, method = method)
  }
  perm_r <- perm_r[is.finite(perm_r)]
  p_value <- mean(perm_r >= r_obs)

  list(r = r_obs, p_value = p_value, n = n,
       v1 = v1, v2 = v2, perm_r = perm_r, message = "OK")
}

# ---------------------------------------------------------------------------
# Module server
# ---------------------------------------------------------------------------

server_isolation_by_distance <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── DB truth layer ─────────────────────────────────────────────────────
    `%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
    db_tick    <- reactive({ rv$db_tick })
    con_r      <- reactive({ shiny::req(rv$con); rv$con })
    tbl_meta_r <- reactive({ rv$tbl_meta %|||% "meta" })

    sql_ident <- function(con, x) DBI::dbQuoteIdentifier(con, x)

    db_ready <- reactive({
      db_tick()
      con <- con_r()
      shiny::req(isTRUE(rv$db_ready))
      shiny::validate(
        shiny::need(DBI::dbExistsTable(con, tbl_meta_r()), "DuckDB meta table missing.")
      )
      TRUE
    })

    # ── GPS centroids ──────────────────────────────────────────────────────
    coords_r <- reactive({
      db_ready()
      con <- con_r()
      cols <- DBI::dbGetQuery(con, sprintf(
        "SELECT column_name FROM information_schema.columns WHERE table_name = '%s'",
        tbl_meta_r()))$column_name
      shiny::validate(shiny::need(
        all(c("Latitude", "Longitude") %in% cols),
        "No GPS data found in meta table."))
      df <- DBI::dbGetQuery(con, sprintf(
        "SELECT Population,
                AVG(CAST(Latitude  AS DOUBLE)) AS Latitude,
                AVG(CAST(Longitude AS DOUBLE)) AS Longitude
         FROM %s
         WHERE Population IS NOT NULL
           AND Latitude IS NOT NULL AND Longitude IS NOT NULL
         GROUP BY Population ORDER BY Population",
        as.character(sql_ident(con, tbl_meta_r()))))
      shiny::validate(shiny::need(nrow(df) >= 2L,
        "At least 2 populations with GPS coordinates required."))
      df
    })

    # ── Raw string genotypes ("a/b" per locus) ─────────────────────────────
    raw_genos_r <- reactive({
      db_ready()
      con <- con_r()
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
        "No marker_cols_raw found in params."))

      geno_fmt <- tryCatch(
        DBI::dbGetQuery(con, "SELECT value FROM params WHERE key='genotype_format'")$value[1L],
        error = function(e) NA_character_)
      if (is.na(geno_fmt) || !nzchar(geno_fmt))
        geno_fmt <- if (any(grepl("(_1|\\.[0-9]+)$", marker_cols_raw))) "paired" else "string"

      keep <- unique(marker_cols_raw)
      keep_sql <- paste(vapply(keep, function(x)
        as.character(DBI::dbQuoteIdentifier(con, x)), character(1L)), collapse = ", ")
      raw_df <- as.data.frame(
        DBI::dbGetQuery(con, sprintf("SELECT rowid AS individual, %s FROM %s",
                                     keep_sql,
                                     as.character(DBI::dbQuoteIdentifier(con, tbl_raw)))),
        stringsAsFactors = FALSE)
      shiny::validate(shiny::need(nrow(raw_df) > 0L, "No rows in raw table."))

      meta_pop <- DBI::dbGetQuery(con, sprintf(
        "SELECT individual, Population FROM %s WHERE Population IS NOT NULL",
        as.character(sql_ident(con, tbl_meta_r()))))
      raw_df$Population <- meta_pop$Population[match(raw_df$individual, meta_pop$individual)]
      pop_vector <- as.character(raw_df$Population)

      pick_b <- function(locus, nms) {
        cands <- c(paste0(locus, "_1"), paste0(locus, "_2"),
                   paste0(locus, ".", 1:9), paste0(locus, "_", 1:9))
        hit <- cands[cands %in% nms]; if (length(hit)) hit[1L] else NA_character_
      }
      if (identical(geno_fmt, "paired")) {
        nms <- names(raw_df)
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
        "No locus columns could be reconstructed."))
      list(hap_df = hap_df, pop_vector = pop_vector)
    })

    # ═══════════════════════════════════════════════════════════════════════
    # TAB 1 — Pairwise Distances (FreeNA)
    # ═══════════════════════════════════════════════════════════════════════

    freena_results_r <- eventReactive(input$run_freena, {
      shiny::req(db_ready())
      rg <- raw_genos_r()
      hap_df <- rg$hap_df
      pop_vector <- rg$pop_vector

      # Filter populations with >= 2 individuals
      pop_counts <- table(pop_vector)
      valid_pops <- names(pop_counts[pop_counts >= 2L])
      shiny::validate(shiny::need(
        length(valid_pops) >= 2L,
        "Need at least 2 populations with >= 2 individuals each."))
      keep <- pop_vector %in% valid_pops
      hap_sub <- hap_df[keep, , drop = FALSE]
      pop_sub <- pop_vector[keep]

      shiny::validate(shiny::need(
        ncol(hap_sub) >= 5L,
        "FreeNA bootstrap requires >= 5 loci."))

      n_boot <- as.integer(input$n_boot_loci)

      withProgress(message = "Running FreeNA pipeline...", value = 0, {
        res <- .freena_pipeline(
          hap_sub, pop_sub,
          n_boot   = n_boot,
          calc_fst = input$calc_fst,
          ena_corr = input$ena_corr,
          calc_cs  = input$calc_cs,
          ina_corr = input$ina_corr
        )
        incProgress(1, detail = "Done")
      })

      res$ninds <- sum(keep)
      res
    })

    # Summary boxes
    output$box_ninds <- renderValueBox({
      r <- freena_results_r()
      valueBox(r$ninds, "Individuals", icon = icon("user"), color = "teal")
    })
    output$box_npops <- renderValueBox({
      r <- freena_results_r()
      valueBox(length(r$pops), "Populations", icon = icon("users"), color = "blue")
    })
    output$box_nloci <- renderValueBox({
      r <- freena_results_r()
      valueBox(r$nloc, "Loci", icon = icon("dna"), color = "purple")
    })
    output$box_npairs <- renderValueBox({
      r <- freena_results_r()
      np <- length(r$pops)
      valueBox(np * (np - 1L) / 2L, "Population pairs",
               icon = icon("project-diagram"), color = "orange")
    })

    # Null allele frequency table
    output$rd_table <- DT::renderDT({
      r <- freena_results_r()
      df <- as.data.frame(round(r$rd_mat, 4))
      df <- cbind(Locus = rownames(df), df)
      DT::datatable(df, rownames = FALSE,
        options = list(scrollX = TRUE, scrollY = "240px",
                       pageLength = 50, dom = "rt"),
        class = "compact stripe hover")
    })

    # FST matrix
    output$fst_matrix_table <- DT::renderDT({
      r <- freena_results_r()
      mat <- if (input$fst_matrix_choice == "fst_ena") r$fst_ena_mat else r$fst_mat
      df <- as.data.frame(round(mat, 4))
      df <- cbind(Population = rownames(df), df)
      DT::datatable(df, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 20, dom = "t"),
        class = "compact stripe hover")
    })

    # CS matrix
    output$cs_matrix_table <- DT::renderDT({
      r <- freena_results_r()
      mat <- if (input$cs_matrix_choice == "cs_ina") r$cs_ina_mat else r$cs_mat
      df <- as.data.frame(round(mat, 4))
      df <- cbind(Population = rownames(df), df)
      DT::datatable(df, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 20, dom = "t"),
        class = "compact stripe hover")
    })

    # Detailed pairwise table with CI
    pairwise_df_r <- reactive({
      r <- freena_results_r()
      pops <- r$pops
      np <- length(pops)
      rows <- vector("list", np * (np - 1L) / 2L)
      k <- 1L
      for (i in seq_len(np - 1L)) {
        for (j in (i + 1L):np) {
          row <- data.frame(
            Pop1 = pops[i], Pop2 = pops[j],
            FST      = r$fst_mat[i, j],
            FST_ENA  = r$fst_ena_mat[i, j],
            CS       = r$cs_mat[i, j],
            CS_INA   = r$cs_ina_mat[i, j],
            FST_CI_l = if (!is.null(r$fst_ci_l))     r$fst_ci_l[i, j]     else NA_real_,
            FST_CI_u = if (!is.null(r$fst_ci_u))     r$fst_ci_u[i, j]     else NA_real_,
            ENA_CI_l = if (!is.null(r$fst_ena_ci_l)) r$fst_ena_ci_l[i, j] else NA_real_,
            ENA_CI_u = if (!is.null(r$fst_ena_ci_u)) r$fst_ena_ci_u[i, j] else NA_real_,
            CS_CI_l  = if (!is.null(r$cs_ci_l))      r$cs_ci_l[i, j]      else NA_real_,
            CS_CI_u  = if (!is.null(r$cs_ci_u))      r$cs_ci_u[i, j]      else NA_real_,
            INA_CI_l = if (!is.null(r$cs_ina_ci_l))  r$cs_ina_ci_l[i, j]  else NA_real_,
            INA_CI_u = if (!is.null(r$cs_ina_ci_u))  r$cs_ina_ci_u[i, j]  else NA_real_,
            stringsAsFactors = FALSE
          )
          rows[[k]] <- row
          k <- k + 1L
        }
      }
      do.call(rbind, rows)
    })

    output$pairwise_detail_table <- DT::renderDT({
      df <- pairwise_df_r()
      DT::datatable(
        df[, c("Pop1", "Pop2", "FST", "FST_ENA",
               "ENA_CI_l", "ENA_CI_u",
               "CS", "CS_INA", "INA_CI_l", "INA_CI_u")],
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 20, dom = "lrtip"),
        class = "compact stripe hover",
        colnames = c("Pop 1", "Pop 2",
                     "FST", "FST (ENA)",
                     "ENA CI lo", "ENA CI hi",
                     "CS", "CS (INA)",
                     "INA CI lo", "INA CI hi")) %>%
        DT::formatRound(columns = c("FST", "FST_ENA", "CS", "CS_INA",
                                     "ENA_CI_l", "ENA_CI_u",
                                     "INA_CI_l", "INA_CI_u"),
                         digits = 4)
    })

    # Download handlers
    output$dl_pairwise_csv <- downloadHandler(
      filename = function() paste0("FreeNA_pairwise_", Sys.Date(), ".csv"),
      content  = function(file) write.csv(pairwise_df_r(), file, row.names = FALSE)
    )
    output$dl_freena_log <- downloadHandler(
      filename = function() paste0("FreeNA_log_", Sys.Date(), ".txt"),
      content  = function(file) writeLines(freena_results_r()$log, file)
    )

    # ═══════════════════════════════════════════════════════════════════════
    # TAB 2 — Mantel Test
    # ═══════════════════════════════════════════════════════════════════════

    # Reactive: matrix 1
    mat1_r <- reactive({
      src <- input$mat1_source
      if (src == "tab1_fst") {
        r <- shiny::req(freena_results_r())
        if (input$mat1_fst_choice == "fst_ena") r$fst_ena_mat else r$fst_mat
      } else if (src == "tab1_cs") {
        r <- shiny::req(freena_results_r())
        if (input$mat1_cs_choice == "cs_ina") r$cs_ina_mat else r$cs_mat
      } else {
        shiny::req(input$file_mat1)
        .parse_dist_file(input$file_mat1$datapath, input$mat1_format)
      }
    })

    # Reactive: matrix 2
    mat2_r <- reactive({
      src <- input$mat2_source
      if (src == "gps") {
        coords <- coords_r()
        df <- .geo_dist_matrix(coords)
        if (input$use_log_dist) {
          df[df > 0] <- log(df[df > 0])
          diag(df) <- 0
        }
        df
      } else {
        shiny::req(input$file_mat2)
        .parse_dist_file(input$file_mat2$datapath, input$mat2_format)
      }
    })

    # Run Mantel
    mantel_results_r <- eventReactive(input$run_mantel, {
      m1 <- mat1_r()
      m2 <- mat2_r()

      # Align populations
      common <- intersect(rownames(m1), rownames(m2))
      shiny::validate(shiny::need(
        length(common) >= 3L,
        "Need at least 3 populations present in BOTH matrices."))
      m1 <- m1[common, common]; m2 <- m2[common, common]

      n_perm <- as.integer(input$n_perm_mantel)
      method <- input$mantel_method

      withProgress(message = "Running Mantel test...", value = 0.5, {
        res <- .mantel_test(m1, m2, n_perm = n_perm, method = method)
        incProgress(0.5, detail = "Done")
      })

      res$mat1 <- m1
      res$mat2 <- m2
      res$method <- method
      res$n_perm <- n_perm
      res
    })

    output$box_mantel_r <- renderValueBox({
      r <- mantel_results_r()
      valueBox(if (is.na(r$r)) "NA" else formatC(r$r, format = "f", digits = 4),
               "Mantel r", icon = icon("chart-line"), color = "purple")
    })
    output$box_mantel_p <- renderValueBox({
      r <- mantel_results_r()
      col <- if (!is.na(r$p_value) && r$p_value < 0.05) "green" else "yellow"
      valueBox(if (is.na(r$p_value)) "NA" else formatC(r$p_value, format = "f", digits = 4),
               "P-value", icon = icon("check-circle"), color = col)
    })
    output$box_mantel_n <- renderValueBox({
      r <- mantel_results_r()
      valueBox(r$n, "Pairs", icon = icon("hashtag"), color = "blue")
    })

    output$mantel_summary <- renderPrint({
      r <- mantel_results_r()
      cat("Mantel Test\n")
      cat("===========\n\n")
      cat("Method:      ", r$method, "\n")
      cat("Permutations:", r$n_perm, "\n")
      cat("Pairs:       ", r$n, "\n\n")
      cat(sprintf("Mantel r:    %.4f\n", r$r))
      cat(sprintf("P-value:     %.4f\n\n", r$p_value))
      if (!is.na(r$p_value)) {
        if      (r$p_value < 0.001) cat("Result: Highly significant (p < 0.001) ***\n")
        else if (r$p_value < 0.01)  cat("Result: Very significant (p < 0.01) **\n")
        else if (r$p_value < 0.05)  cat("Result: Significant (p < 0.05) *\n")
        else                        cat("Result: Not significant (p >= 0.05)\n")
      }
    })

    output$mantel_plot <- plotly::renderPlotly({
      r <- mantel_results_r()
      if (is.null(r$v1) || length(r$v1) < 3L) return(NULL)

      df <- data.frame(x = r$v1, y = r$v2)
      fit <- lm(y ~ x, data = df)

      plotly::plot_ly() %>%
        plotly::add_markers(
          data = df, x = ~x, y = ~y,
          marker = list(color = "#2CBF9F", size = 7, opacity = 0.85),
          name = "Pairs",
          hoverinfo = "text",
          text = ~paste0("Mat1: ", round(x, 4), "<br>Mat2: ", round(y, 4))
        ) %>%
        plotly::add_lines(
          data = data.frame(x = df$x, y = fitted(fit)),
          x = ~x, y = ~y,
          line = list(color = "#B40F20", width = 2),
          name = sprintf("Regression (r = %.3f)", r$r)
        ) %>%
        plotly::layout(
          title = list(
            text = sprintf("Mantel r = %.4f, p = %.4f", r$r, r$p_value),
            font = list(size = 14)),
          xaxis = list(title = "Matrix 1 distance"),
          yaxis = list(title = "Matrix 2 distance"),
          legend = list(x = 0.02, y = 0.98, bgcolor = "rgba(255,255,255,0.8)")
        )
    })

  })
}