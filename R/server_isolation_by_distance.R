# server_isolation_by_distance.R
# Two tabs: Pairwise Genetic Distances + Mantel Test

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Cavalli-Sforza & Edwards (1967) chord distance between two populations
# freq1, freq2: named vectors of allele frequencies
.calc_cs_distance <- function(freq1, freq2) {
  all_alleles <- union(names(freq1), names(freq2))
  if (length(all_alleles) == 0) return(NA_real_)
  
  p <- sapply(all_alleles, function(a) ifelse(a %in% names(freq1), freq1[a], 0))
  q <- sapply(all_alleles, function(a) ifelse(a %in% names(freq2), freq2[a], 0))
  
  # Normalized frequencies
  p <- p / sum(p)
  q <- q / sum(q)
  
  # Chord distance
  prod_sum <- sum(sqrt(p * q))
  prod_sum <- min(max(prod_sum, 0), 1)
  
  Dc <- (2 / pi) * sqrt(2 * (1 - prod_sum))
  return(Dc)
}

# Compute allele frequencies for a population at a locus
# genotypes: character vector of "a/b" format
.calc_allele_freq <- function(genotypes) {
  alleles <- c()
  for (g in genotypes) {
    if (is.na(g) || trimws(g) == "" || g == "0/0" || g == "0") next
    parts <- strsplit(g, "/")[[1]]
    if (length(parts) == 2) {
      a1 <- trimws(parts[1])
      a2 <- trimws(parts[2])
      if (a1 != "0" && a1 != "") alleles <- c(alleles, a1)
      if (a2 != "0" && a2 != "") alleles <- c(alleles, a2)
    }
  }
  
  if (length(alleles) == 0) return(setNames(numeric(0), character(0)))
  
  freq_table <- table(alleles)
  freq <- as.numeric(freq_table) / sum(freq_table)
  names(freq) <- names(freq_table)
  return(freq)
}

# Compute pairwise genetic distances with bootstrap over loci
# hap_df: data.frame with loci as columns, "a/b" genotypes
# pop_vector: population assignments
# n_boot: number of bootstrap replicates
# Returns: list with distance matrices and confidence intervals
.pairwise_genetic_distances <- function(hap_df, pop_vector, n_boot = 1000, 
                                        calc_fst = TRUE, calc_cs = TRUE) {
  pops <- sort(unique(pop_vector))
  loci <- colnames(hap_df)
  n_pops <- length(pops)
  n_loci <- length(loci)
  
  if (n_pops < 2) stop("Need at least 2 populations")
  if (n_loci < 2) stop("Need at least 2 loci for bootstrap")
  
  # Initialize result matrices
  cs_dist <- matrix(NA, n_pops, n_pops, dimnames = list(pops, pops))
  fst_dist <- matrix(NA, n_pops, n_pops, dimnames = list(pops, pops))
  
  # Bootstrap storage
  cs_boot <- array(NA, dim = c(n_pops, n_pops, n_boot), 
                   dimnames = list(pops, pops, NULL))
  fst_boot <- array(NA, dim = c(n_pops, n_pops, n_boot),
                    dimnames = list(pops, pops, NULL))
  
  # Function to compute distances for a set of loci
  .compute_distances <- function(locus_subset) {
    cs_mat <- matrix(NA, n_pops, n_pops)
    fst_mat <- matrix(NA, n_pops, n_pops)
    
    for (i in 1:(n_pops - 1)) {
      for (j in (i + 1):n_pops) {
        pop_i <- pops[i]
        pop_j <- pops[j]
        idx_i <- which(pop_vector == pop_i)
        idx_j <- which(pop_vector == pop_j)
        
        # Compute per-locus distances
        cs_locus <- numeric(length(locus_subset))
        fst_locus <- numeric(length(locus_subset))
        
        for (k in seq_along(locus_subset)) {
          locus <- locus_subset[k]
          geno_i <- hap_df[[locus]][idx_i]
          geno_j <- hap_df[[locus]][idx_j]
          
          # Allele frequencies
          freq_i <- .calc_allele_freq(geno_i)
          freq_j <- .calc_allele_freq(geno_j)
          
          if (length(freq_i) > 0 && length(freq_j) > 0) {
            # Cavalli-Sforza distance
            if (calc_cs) {
              cs_locus[k] <- .calc_cs_distance(freq_i, freq_j)
            }
            
            # Simple FST approximation (Nei's Gst)
            if (calc_fst) {
              all_alleles <- union(names(freq_i), names(freq_j))
              p <- sapply(all_alleles, function(a) ifelse(a %in% names(freq_i), freq_i[a], 0))
              q <- sapply(all_alleles, function(a) ifelse(a %in% names(freq_j), freq_j[a], 0))
              
              p_bar <- (p + q) / 2
              H_s <- (sum(p * (1 - p)) + sum(q * (1 - q))) / 2
              H_t <- sum(p_bar * (1 - p_bar))
              
              if (H_t > 0) {
                fst_locus[k] <- (H_t - H_s) / H_t
              } else {
                fst_locus[k] <- 0
              }
            }
          } else {
            cs_locus[k] <- NA
            fst_locus[k] <- NA
          }
        }
        
        # Average over loci
        if (calc_cs) cs_mat[i, j] <- cs_mat[j, i] <- mean(cs_locus, na.rm = TRUE)
        if (calc_fst) fst_mat[i, j] <- fst_mat[j, i] <- mean(fst_locus, na.rm = TRUE)
      }
    }
    
    list(cs = cs_mat, fst = fst_mat)
  }
  
  # Compute observed distances
  obs <- .compute_distances(loci)
  cs_dist <- obs$cs
  fst_dist <- obs$fst
  
  # Bootstrap over loci
  if (n_boot > 0) {
    for (b in 1:n_boot) {
      boot_loci <- sample(loci, n_loci, replace = TRUE)
      boot_res <- .compute_distances(boot_loci)
      cs_boot[, , b] <- boot_res$cs
      fst_boot[, , b] <- boot_res$fst
    }
  }
  
  # Compute confidence intervals
  cs_ci_l <- apply(cs_boot, c(1, 2), function(x) quantile(x, 0.025, na.rm = TRUE))
  cs_ci_u <- apply(cs_boot, c(1, 2), function(x) quantile(x, 0.975, na.rm = TRUE))
  fst_ci_l <- apply(fst_boot, c(1, 2), function(x) quantile(x, 0.025, na.rm = TRUE))
  fst_ci_u <- apply(fst_boot, c(1, 2), function(x) quantile(x, 0.975, na.rm = TRUE))
  
  list(
    cs_dist = cs_dist,
    fst_dist = fst_dist,
    cs_ci_l = cs_ci_l,
    cs_ci_u = cs_ci_u,
    fst_ci_l = fst_ci_l,
    fst_ci_u = fst_ci_u,
    pops = pops,
    n_boot = n_boot
  )
}

# Mantel test function
# matrix1, matrix2: distance matrices (square or rectangular format)
# n_perm: number of permutations
# method: correlation method
.mantel_test <- function(matrix1, matrix2, n_perm = 9999, method = "pearson") {
  # Convert to vectors (lower triangle)
  vec1 <- as.vector(matrix1[lower.tri(matrix1)])
  vec2 <- as.vector(matrix2[lower.tri(matrix2)])
  
  # Remove NAs
  ok <- complete.cases(vec1, vec2)
  vec1 <- vec1[ok]
  vec2 <- vec2[ok]
  
  n <- length(vec1)
  if (n < 3) {
    return(list(r = NA, p_value = NA, n = n, 
                message = "Not enough data points for Mantel test"))
  }
  
  # Observed correlation
  r_obs <- cor(vec1, vec2, method = method)
  
  # Permutation test
  perm_r <- numeric(n_perm)
  for (i in 1:n_perm) {
    perm_vec2 <- sample(vec2)
    perm_r[i] <- cor(vec1, perm_vec2, method = method)
  }
  
  # P-value (one-sided, positive correlation)
  p_value <- mean(perm_r >= r_obs)
  
  list(
    r = r_obs,
    p_value = p_value,
    n = n,
    perm_r = perm_r,
    vec1 = vec1,
    vec2 = vec2
  )
}

# ---------------------------------------------------------------------------
# Module server
# ---------------------------------------------------------------------------

server_isolation_by_distance <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    
    # ── DB truth layer ─────────────────────────────────────────────────────
    `%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
    db_tick <- reactive({ rv$db_tick })
    con_r <- reactive({ shiny::req(rv$con); rv$con })
    tbl_meta_r <- reactive({ rv$tbl_meta %||% "meta" })
    
    db_ready <- reactive({
      db_tick()
      con <- con_r()
      shiny::req(isTRUE(rv$db_ready))
      shiny::validate(
        shiny::need(DBI::dbExistsTable(con, tbl_meta_r()), "DuckDB meta table missing.")
      )
      TRUE
    })
    
    # ── Raw string genotypes ────────────────────────────────────────────────
    raw_genos_r <- reactive({
      db_ready()
      con <- con_r()
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
        sql_ident(con, tbl_meta_r())))
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
        "No locus columns could be reconstructed from raw table."))
      list(hap_df = hap_df, pop_vector = pop_vector)
    })
    
    # ═══════════════════════════════════════════════════════════════════════
    # TAB 1: Pairwise Genetic Distances
    # ═══════════════════════════════════════════════════════════════════════
    
    dist_results_r <- eventReactive(input$run_dist, {
      shiny::req(db_ready())
      rg <- raw_genos_r()
      
      hap_df <- rg$hap_df
      pop_vector <- rg$pop_vector
      
      # Filter populations with enough individuals
      pop_counts <- table(pop_vector)
      valid_pops <- names(pop_counts[pop_counts >= 2])
      keep <- pop_vector %in% valid_pops
      
      shiny::validate(shiny::need(
        length(unique(pop_vector[keep])) >= 2,
        "Need at least 2 populations with ≥2 individuals each."))
      
      hap_sub <- hap_df[keep, , drop = FALSE]
      pop_sub <- pop_vector[keep]
      
      n_boot <- as.integer(input$n_boot_dist)
      calc_fst <- input$calc_fst
      calc_cs <- input$calc_cs
      
      withProgress(message = "Computing pairwise genetic distances...", value = 0, {
        res <- .pairwise_genetic_distances(
          hap_sub, pop_sub, 
          n_boot = n_boot,
          calc_fst = calc_fst,
          calc_cs = calc_cs
        )
        incProgress(1, detail = "Done!")
      })
      
      res
    })
    
    # Summary boxes
    output$box_npops_dist <- renderValueBox({
      r <- dist_results_r()
      valueBox(length(r$pops), "Populations", 
               icon = icon("users"), color = "teal")
    })
    
    output$box_nloci_dist <- renderValueBox({
      rg <- raw_genos_r()
      valueBox(ncol(rg$hap_df), "Loci", 
               icon = icon("dna"), color = "blue")
    })
    
    output$box_npairs_dist <- renderValueBox({
      r <- dist_results_r()
      n_pops <- length(r$pops)
      n_pairs <- n_pops * (n_pops - 1) / 2
      valueBox(n_pairs, "Population pairs", 
               icon = icon("project-diagram"), color = "purple")
    })
    
    output$box_boot_dist <- renderValueBox({
      r <- dist_results_r()
      valueBox(r$n_boot, "Bootstrap replicates", 
               icon = icon("sync"), color = "orange")
    })
    
    # Distance matrix table
    output$dist_matrix_table <- DT::renderDT({
      r <- dist_results_r()
      
      # Choose which distance to display
      if (input$calc_cs) {
        mat <- r$cs_dist
        title <- "Cavalli-Sforza Distance"
      } else if (input$calc_fst) {
        mat <- r$fst_dist
        title <- "FST"
      } else {
        return(NULL)
      }
      
      df <- as.data.frame(round(mat, 4))
      df <- cbind(Population = rownames(df), df)
      
      DT::datatable(df, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15, dom = "t"),
        class = "compact stripe hover",
        caption = title)
    })
    
    # Detailed pairwise table
    output$pairwise_detail_table <- DT::renderDT({
      r <- dist_results_r()
      pops <- r$pops
      n_pops <- length(pops)
      
      rows <- list()
      for (i in 1:(n_pops - 1)) {
        for (j in (i + 1):n_pops) {
          row <- data.frame(
            Population1 = pops[i],
            Population2 = pops[j],
            stringsAsFactors = FALSE
          )
          
          if (input$calc_cs) {
            row$CS_Distance <- r$cs_dist[i, j]
            row$CS_CI_lower <- r$cs_ci_l[i, j]
            row$CS_CI_upper <- r$cs_ci_u[i, j]
          }
          
          if (input$calc_fst) {
            row$FST <- r$fst_dist[i, j]
            row$FST_CI_lower <- r$fst_ci_l[i, j]
            row$FST_CI_upper <- r$fst_ci_u[i, j]
          }
          
          rows[[length(rows) + 1]] <- row
        }
      }
      
      df <- do.call(rbind, rows)
      
      # Round for display
      num_cols <- sapply(df, is.numeric)
      df[num_cols] <- lapply(df[num_cols], function(x) round(x, 4))
      
      DT::datatable(df, rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 20, dom = "lrtip"),
        class = "compact stripe hover") %>%
        DT::formatStyle(columns = names(df)[sapply(df, is.numeric)],
          backgroundColor = DT::styleInterval(
            c(0.05, 0.15, 0.25),
            c("#d4edda", "#fff3cd", "#f8d7da", "#c3002f22")))
    })
    
    # Download handler
    output$dl_pairwise_csv <- downloadHandler(
      filename = function() {
        paste0("pairwise_genetic_distances_", Sys.Date(), ".csv")
      },
      content = function(file) {
        r <- dist_results_r()
        pops <- r$pops
        n_pops <- length(pops)
        
        rows <- list()
        for (i in 1:(n_pops - 1)) {
          for (j in (i + 1):n_pops) {
            row <- data.frame(
              Population1 = pops[i],
              Population2 = pops[j],
              stringsAsFactors = FALSE
            )
            
            if (input$calc_cs) {
              row$CS_Distance <- r$cs_dist[i, j]
              row$CS_CI_lower <- r$cs_ci_l[i, j]
              row$CS_CI_upper <- r$cs_ci_u[i, j]
            }
            
            if (input$calc_fst) {
              row$FST <- r$fst_dist[i, j]
              row$FST_CI_lower <- r$fst_ci_l[i, j]
              row$FST_CI_upper <- r$fst_ci_u[i, j]
            }
            
            rows[[length(rows) + 1]] <- row
          }
        }
        
        df <- do.call(rbind, rows)
        write.csv(df, file, row.names = FALSE)
      }
    )
    
    # ═══════════════════════════════════════════════════════════════════════
    # TAB 2: Mantel Test
    # ═══════════════════════════════════════════════════════════════════════
    
    mantel_results_r <- eventReactive(input$run_mantel, {
      n_perm <- as.integer(input$n_perm_mantel)
      method <- input$mantel_method
      
      if (input$mantel_data_source == "computed") {
        # Use computed pairwise distances + geographic distances
        shiny::req(dist_results_r())
        r <- dist_results_r()
        
        # Genetic distance matrix
        if (input$calc_cs) {
          gen_dist <- r$cs_dist
        } else if (input$calc_fst) {
          gen_dist <- r$fst_dist
        } else {
          shiny::validate("Please compute at least one distance measure")
        }
        
        # Geographic distance matrix
        coords <- coords_r()
        shiny::validate(shiny::need(nrow(coords) >= 2,
          "Need GPS coordinates for geographic distances"))
        
        geo_dist <- .geo_dist_matrix(coords)
        
        # Align populations
        common_pops <- intersect(rownames(gen_dist), rownames(geo_dist))
        shiny::validate(shiny::need(length(common_pops) >= 3,
          "Need at least 3 populations with both genetic and geographic data"))
        
        gen_dist <- gen_dist[common_pops, common_pops]
        geo_dist <- geo_dist[common_pops, common_pops]
        
      } else {
        # Upload files
        shiny::req(input$file_gen_dist)
        shiny::req(input$file_geo_dist)
        
        gen_dist <- read.csv(input$file_gen_dist$datapath, row.names = 1)
        geo_dist <- read.csv(input$file_geo_dist$datapath, row.names = 1)
        
        gen_dist <- as.matrix(gen_dist)
        geo_dist <- as.matrix(geo_dist)
      }
      
      # Run Mantel test
      withProgress(message = "Running Mantel test...", value = 0.5, {
        result <- .mantel_test(gen_dist, geo_dist, n_perm = n_perm, method = method)
        incProgress(0.5, detail = "Done!")
      })
      
      result$gen_dist <- gen_dist
      result$geo_dist <- geo_dist
      result$method <- method
      result$n_perm <- n_perm
      
      result
    })
    
    # Mantel summary boxes
    output$box_mantel_r <- renderValueBox({
      r <- mantel_results_r()
      valueBox(
        if (is.na(r$r)) "NA" else formatC(r$r, format = "f", digits = 4),
        "Mantel r",
        icon = icon("chart-line"), color = "purple"
      )
    })
    
    output$box_mantel_p <- renderValueBox({
      r <- mantel_results_r()
      col <- if (!is.na(r$p_value) && r$p_value < 0.05) "green" else "yellow"
      valueBox(
        if (is.na(r$p_value)) "NA" else formatC(r$p_value, format = "f", digits = 4),
        "P-value",
        icon = icon("check-circle"), color = col
      )
    })
    
    output$box_mantel_n <- renderValueBox({
      r <- mantel_results_r()
      valueBox(r$n, "Data points",
               icon = icon("hashtag"), color = "blue")
    })
    
    # Mantel summary text
    output$mantel_summary <- renderPrint({
      r <- mantel_results_r()
      cat("Mantel Test Results\n")
      cat("===================\n\n")
      cat("Method:", r$method, "\n")
      cat("Permutations:", r$n_perm, "\n")
      cat("Sample size:", r$n, "\n\n")
      cat("Mantel r:", formatC(r$r, format = "f", digits = 4), "\n")
      cat("P-value:", formatC(r$p_value, format = "f", digits = 4), "\n\n")
      
      if (!is.na(r$p_value)) {
        if (r$p_value < 0.001) {
          cat("Result: Highly significant (p < 0.001)\n")
        } else if (r$p_value < 0.01) {
          cat("Result: Very significant (p < 0.01)\n")
        } else if (r$p_value < 0.05) {
          cat("Result: Significant (p < 0.05)\n")
        } else {
          cat("Result: Not significant (p >= 0.05)\n")
        }
      }
    })
    
    # Mantel plot
    output$mantel_plot <- plotly::renderPlotly({
      r <- mantel_results_r()
      
      if (is.null(r$vec1) || is.null(r$vec2)) return(NULL)
      
      df <- data.frame(
        x = r$vec1,
        y = r$vec2
      )
      
      # Fit regression line
      fit <- lm(y ~ x, data = df)
      
      plotly::plot_ly() %>%
        plotly::add_markers(
          data = df, x = ~x, y = ~y,
          marker = list(color = "#2CBF9F", size = 6, opacity = 0.7),
          name = "Data points",
          hoverinfo = "text",
          text = ~paste0("Genetic: ", round(x, 4), 
                        "<br>Geographic: ", round(y, 4))
        ) %>%
        plotly::add_lines(
          data = data.frame(x = df$x, y = fitted(fit)),
          x = ~x, y = ~y,
          line = list(color = "#B40F20", width = 2, dash = "dash"),
          name = paste0("Regression (r = ", formatC(r$r, format = "f", digits = 3), ")")
        ) %>%
        plotly::layout(
          title = list(
            text = paste0("Mantel Test (p = ", formatC(r$p_value, format = "f", digits = 4), ")"),
            font = list(size = 14)
          ),
          xaxis = list(title = "Genetic Distance"),
          yaxis = list(title = "Geographic Distance"),
          legend = list(x = 0.02, y = 0.98, bgcolor = "rgba(255,255,255,0.8)")
        )
    })
    
  })
}