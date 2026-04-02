suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(vegan)
  library(openxlsx)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(tibble)
  library(ggplot2)
  library(readxl)
  library(grid)
})

tsv_path <- "CLR_3254 for dbRDA-single food.tsv"
out_dir <- dirname(tsv_path)

id_col <- "SEQN"
group_col <- "Oral_CDC"
diet_cols_idx <- 3:40
micro_cols_idx <- 41:117

energy_adjust <- TRUE
energy_col_orig <- "Energy_intake_kcal"

perms <- 1001
distance_method <- "euclidean"

label_map <- c(
  "Overall" = "Overall",
  "1" = "Healthy",
  "2" = "Caries",
  "3" = "Periodontitis",
  "4" = "Caries_Periodontitis",
  "777" = "Edentulous"
)

make_clean_names <- function(x) {
  y <- x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    tolower()
  dup <- duplicated(y)
  if (any(dup)) {
    counts <- ave(seq_along(y), y, FUN = seq_along)
    y[dup] <- paste0(y[dup], "_", counts[dup])
  }
  y
}

safe_scale <- function(df_num) {
  sds <- vapply(df_num, sd, numeric(1), na.rm = TRUE)
  keep <- sds > 0 & !is.na(sds)
  df_num <- df_num[, keep, drop = FALSE]
  df_num[] <- scale(df_num, center = TRUE, scale = TRUE)
  df_num
}

drop_linear_dependencies <- function(X) {
  X <- as.data.frame(X)
  if (ncol(X) <= 1) return(X)
  qrX <- qr(as.matrix(X))
  keep_idx <- sort(qrX$pivot[seq_len(qrX$rank)])
  X[, keep_idx, drop = FALSE]
}

`%||%` <- function(x, y) {
  if (is.na(x) || length(x) == 0) y else x
}

prepare_diet_matrix <- function(df, diet_cols, energy_col_orig, energy_adjust = TRUE) {
  energy_col_clean <- make_clean_names(energy_col_orig)
  diet_mat <- df[, diet_cols, drop = FALSE]
  diet_mat[] <- lapply(diet_mat, as.numeric)
  
  if (energy_col_clean %in% colnames(diet_mat)) {
    diet_mat <- diet_mat[, setdiff(colnames(diet_mat), energy_col_clean), drop = FALSE]
  }
  
  if (energy_adjust) {
    if (!energy_col_clean %in% colnames(df)) {
      stop("Energy column not found: ", energy_col_clean)
    }
    
    energy_vec <- as.numeric(df[[energy_col_clean]])
    
    for (v in colnames(diet_mat)) {
      this_y <- as.numeric(diet_mat[[v]])
      if (all(is.na(this_y))) next
      fit <- try(lm(this_y ~ energy_vec), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        diet_mat[[v]] <- resid(fit)
      }
    }
  }
  
  safe_scale(diet_mat)
}

raw <- read_tsv(tsv_path, na = c("", "NA"), show_col_types = FALSE)

orig_names <- names(raw)
clean_names <- make_clean_names(orig_names)
colnames(raw) <- clean_names

write_csv(
  tibble(original = orig_names, cleaned = clean_names),
  file.path(out_dir, "dbrda_name_mapping.csv")
)

id_col_clean <- make_clean_names(id_col)
group_col_clean <- make_clean_names(group_col)
diet_cols <- names(raw)[diet_cols_idx]
micro_cols <- names(raw)[micro_cols_idx]

stopifnot(id_col_clean %in% names(raw), group_col_clean %in% names(raw))
stopifnot(all(diet_cols %in% names(raw)))
stopifnot(all(micro_cols %in% names(raw)))

keep_mask <- complete.cases(raw[, diet_cols])
df <- raw[keep_mask, , drop = FALSE]

diet_mat_all <- prepare_diet_matrix(
  df = df,
  diet_cols = diet_cols,
  energy_col_orig = energy_col_orig,
  energy_adjust = energy_adjust
)

abun <- df[, c(id_col_clean, micro_cols), drop = FALSE]
abun[is.na(abun)] <- 0
comm_mat_all <- as.matrix(abun[, micro_cols, drop = FALSE])
storage.mode(comm_mat_all) <- "double"
rownames(comm_mat_all) <- as.character(df[[id_col_clean]])

groups <- sort(unique(df[[group_col_clean]]))
groups_all <- c("Overall", groups)

run_dbrda_for_group <- function(gval) {
  if (identical(gval, "Overall")) {
    idx <- seq_len(nrow(df))
  } else {
    idx <- which(df[[group_col_clean]] == gval)
  }
  
  n <- length(idx)
  g_label <- dplyr::recode(as.character(gval), !!!label_map, .default = as.character(gval))
  
  if (n < 4) {
    return(list(
      summary = tibble(
        Group = as.character(gval),
        Group_Label = g_label,
        n = n,
        R2_adj = NA_real_,
        F_overall = NA_real_,
        p_overall = NA_real_,
        constrained_axes = NA_integer_,
        note = "Insufficient samples"
      ),
      terms = tibble(
        Group = as.character(gval),
        Term = NA_character_,
        F = NA_real_,
        p = NA_real_,
        R2 = NA_real_
      ),
      sites = tibble(
        Group = as.character(gval),
        SEQN = NA_character_,
        CAP1 = NA_real_,
        CAP2 = NA_real_
      ),
      biplot = tibble(
        Group = as.character(gval),
        Variable = NA_character_,
        CAP1 = NA_real_,
        CAP2 = NA_real_
      )
    ))
  }
  
  comm_g <- comm_mat_all[idx, , drop = FALSE]
  env_g <- diet_mat_all[idx, , drop = FALSE]
  
  sds_g <- vapply(env_g, sd, numeric(1), na.rm = TRUE)
  keep <- sds_g > 0 & !is.na(sds_g)
  env_g <- env_g[, keep, drop = FALSE]
  env_g <- drop_linear_dependencies(env_g)
  
  if (ncol(env_g) < 1) {
    return(list(
      summary = tibble(
        Group = as.character(gval),
        Group_Label = g_label,
        n = n,
        R2_adj = NA_real_,
        F_overall = NA_real_,
        p_overall = NA_real_,
        constrained_axes = NA_integer_,
        note = "No usable dietary constraints"
      ),
      terms = tibble(
        Group = as.character(gval),
        Term = NA_character_,
        F = NA_real_,
        p = NA_real_,
        R2 = NA_real_
      ),
      sites = tibble(
        Group = as.character(gval),
        SEQN = rownames(comm_g),
        CAP1 = NA_real_,
        CAP2 = NA_real_
      ),
      biplot = tibble(
        Group = as.character(gval),
        Variable = NA_character_,
        CAP1 = NA_real_,
        CAP2 = NA_real_
      )
    ))
  }
  
  mod <- capscale(comm_g ~ ., data = as.data.frame(env_g), distance = distance_method)
  
  anov_overall <- tryCatch(anova.cca(mod, permutations = perms), error = function(e) NULL)
  F_overall <- NA_real_
  p_overall <- NA_real_
  
  if (!is.null(anov_overall)) {
    atab <- as.data.frame(anov_overall)
    if (nrow(atab) >= 1 && all(c("F", "Pr(>F)") %in% colnames(atab))) {
      F_overall <- suppressWarnings(as.numeric(atab$F[1]))
      p_overall <- suppressWarnings(as.numeric(atab$`Pr(>F)`[1]))
    }
  }
  
  anov_terms <- tryCatch(anova.cca(mod, permutations = perms, by = "terms"), error = function(e) NULL)
  terms_tbl <- tibble(
    Group = character(0),
    Term = character(0),
    F = numeric(0),
    p = numeric(0),
    R2 = numeric(0)
  )
  
  if (!is.null(anov_terms)) {
    ttab <- as.data.frame(anov_terms)
    ttab$Term <- rownames(ttab)
    
    terms_tbl <- ttab %>%
      transmute(
        Group = as.character(gval),
        Term = Term,
        F = suppressWarnings(as.numeric(F)),
        p = suppressWarnings(as.numeric(`Pr(>F)`))
      )
    
    if (nrow(terms_tbl) > 0) {
      R2_full <- suppressWarnings(RsquareAdj(mod)$r.squared)
      R2_terms <- numeric(nrow(terms_tbl))
      names(R2_terms) <- terms_tbl$Term
      
      for (v in terms_tbl$Term) {
        if (!v %in% colnames(env_g)) {
          R2_terms[v] <- NA_real_
          next
        }
        
        others <- setdiff(colnames(env_g), v)
        
        if (length(others) == 0) {
          R2_terms[v] <- R2_full
        } else {
          env_red <- env_g[, others, drop = FALSE]
          mod_red <- tryCatch(
            capscale(comm_g ~ ., data = as.data.frame(env_red), distance = distance_method),
            error = function(e) NULL
          )
          
          if (is.null(mod_red)) {
            R2_terms[v] <- NA_real_
          } else {
            R2_red <- suppressWarnings(RsquareAdj(mod_red)$r.squared)
            R2_terms[v] <- max(0, R2_full - R2_red)
          }
        }
      }
      
      terms_tbl <- terms_tbl %>%
        mutate(R2 = R2_terms[Term])
    }
  }
  
  R2_adj <- suppressWarnings(RsquareAdj(mod)$adj.r.squared)
  eig_con <- eigenvals(mod, constrained = TRUE)
  n_axes <- length(eig_con)
  
  site_sc <- scores(mod, display = "sites", choices = 1:2, scaling = 1)
  sites_tbl <- tibble(
    Group = as.character(gval),
    SEQN = rownames(comm_g),
    CAP1 = site_sc[, 1],
    CAP2 = site_sc[, 2]
  )
  
  bp_sc <- scores(mod, display = "bp", choices = 1:2, scaling = 1)
  if (is.null(bp_sc) || nrow(bp_sc) == 0) {
    biplot_tbl <- tibble(
      Group = as.character(gval),
      Variable = character(0),
      CAP1 = numeric(0),
      CAP2 = numeric(0)
    )
  } else {
    biplot_tbl <- tibble(
      Group = as.character(gval),
      Variable = rownames(bp_sc),
      CAP1 = bp_sc[, 1],
      CAP2 = bp_sc[, 2]
    )
  }
  
  sum_tbl <- tibble(
    Group = as.character(gval),
    Group_Label = g_label,
    n = n,
    R2_adj = R2_adj,
    F_overall = F_overall,
    p_overall = p_overall,
    constrained_axes = n_axes,
    note = NA_character_
  )
  
  list(summary = sum_tbl, terms = terms_tbl, sites = sites_tbl, biplot = biplot_tbl)
}

n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

res_list <- foreach(g = groups_all, .packages = c("vegan", "dplyr", "stringr", "tibble")) %dopar% {
  run_dbrda_for_group(g)
}

stopCluster(cl)

sum_all <- bind_rows(lapply(res_list, `[[`, "summary")) %>%
  mutate(across(c(R2_adj, F_overall, p_overall), as.numeric)) %>%
  arrange(factor(Group, levels = as.character(groups_all)))

terms_all <- bind_rows(lapply(res_list, `[[`, "terms"))
sites_all <- bind_rows(lapply(res_list, `[[`, "sites"))
biplot_all <- bind_rows(lapply(res_list, `[[`, "biplot"))

diet_cols_used <- colnames(diet_mat_all)

wb <- createWorkbook()

addWorksheet(wb, "dbRDA_Summary")
writeData(wb, "dbRDA_Summary", sum_all)
freezePane(wb, "dbRDA_Summary", firstRow = TRUE)
addFilter(wb, "dbRDA_Summary", row = 1, cols = 1:ncol(sum_all))

addWorksheet(wb, "dbRDA_Terms")
writeData(wb, "dbRDA_Terms", terms_all)
freezePane(wb, "dbRDA_Terms", firstRow = TRUE)
addFilter(wb, "dbRDA_Terms", row = 1, cols = 1:ncol(terms_all))

g_levels <- unique(sum_all$Group)

for (g in g_levels) {
  g_label <- unique(sum_all$Group_Label[sum_all$Group == g])[1]
  st <- dplyr::filter(sites_all, Group == g)
  bp <- dplyr::filter(biplot_all, Group == g)
  
  if (nrow(st) > 0) {
    sheet_sites <- paste0("Sites_", g_label)
    addWorksheet(wb, sheet_sites)
    writeData(wb, sheet_sites, st)
    freezePane(wb, sheet_sites, firstRow = TRUE)
    addFilter(wb, sheet_sites, row = 1, cols = 1:ncol(st))
  }
  
  if (nrow(bp) > 0) {
    sheet_bp <- paste0("Biplot_", g_label)
    addWorksheet(wb, sheet_bp)
    writeData(wb, sheet_bp, bp)
    freezePane(wb, sheet_bp, firstRow = TRUE)
    addFilter(wb, sheet_bp, row = 1, cols = 1:ncol(bp))
  }
}

addWorksheet(wb, "Dietary_Columns_Used")
writeData(wb, "Dietary_Columns_Used", tibble(Dietary_Column = diet_cols_used))
freezePane(wb, "Dietary_Columns_Used", firstRow = TRUE)

out_xlsx <- file.path(
  out_dir,
  if (energy_adjust) {
    "dbRDA_Diet_Microbiome_OralCDC_EnergyAdjusted.xlsx"
  } else {
    "dbRDA_Diet_Microbiome_OralCDC.xlsx"
  }
)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

name_map <- read_csv(file.path(out_dir, "dbrda_name_mapping.csv"), show_col_types = FALSE)
to_original <- setNames(name_map$original, name_map$cleaned)

raw_viz <- read_tsv(tsv_path, na = c("", "NA"), show_col_types = FALSE)
orig_names_viz <- names(raw_viz)
clean_names_viz <- make_clean_names(orig_names_viz)
colnames(raw_viz) <- clean_names_viz

diet_cols_viz <- names(raw_viz)[diet_cols_idx]
micro_cols_viz <- names(raw_viz)[micro_cols_idx]

keep_mask_viz <- complete.cases(raw_viz[, diet_cols_viz])
df_viz <- raw_viz[keep_mask_viz, , drop = FALSE]

diet_mat_viz <- prepare_diet_matrix(
  df = df_viz,
  diet_cols = diet_cols_viz,
  energy_col_orig = energy_col_orig,
  energy_adjust = energy_adjust
)

abun_viz <- df_viz[, c(id_col_clean, micro_cols_viz), drop = FALSE]
abun_viz[is.na(abun_viz)] <- 0
comm_mat_viz <- as.matrix(abun_viz[, micro_cols_viz, drop = FALSE])
storage.mode(comm_mat_viz) <- "double"
rownames(comm_mat_viz) <- as.character(df_viz[[id_col_clean]])

axis_props_for_group <- function(g_code) {
  if (identical(g_code, "Overall")) {
    idx <- seq_len(nrow(df_viz))
  } else {
    idx <- which(df_viz[[group_col_clean]] == g_code)
  }
  
  if (length(idx) < 4) return(c(NA_real_, NA_real_))
  
  comm_g <- comm_mat_viz[idx, , drop = FALSE]
  env_g <- diet_mat_viz[idx, , drop = FALSE]
  
  sds_g <- vapply(env_g, sd, numeric(1), na.rm = TRUE)
  keep <- sds_g > 0 & !is.na(sds_g)
  env_g <- env_g[, keep, drop = FALSE]
  env_g <- drop_linear_dependencies(env_g)
  
  if (ncol(env_g) < 1) return(c(NA_real_, NA_real_))
  
  mod <- capscale(comm_g ~ ., data = as.data.frame(env_g), distance = distance_method)
  eig_con <- eigenvals(mod, constrained = TRUE)
  
  if (length(eig_con) == 0) return(c(NA_real_, NA_real_))
  
  prop <- eig_con / sum(eig_con)
  c(round(100 * prop[1], 1), round(100 * (prop[2] %||% 0), 1))
}

sheets <- excel_sheets(out_xlsx)
groups_code <- c("Overall", "1", "2", "3", "4", "777")
avail_groups <- groups_code[
  paste0("Biplot_", label_map[groups_code]) %in% sheets
]

prep_biplot <- function(g_code) {
  g_label <- label_map[[g_code]]
  sheet_bp <- paste0("Biplot_", g_label)
  
  biplot <- read_excel(out_xlsx, sheet = sheet_bp) %>%
    filter(!is.na(CAP1) & !is.na(CAP2)) %>%
    mutate(orig_label = ifelse(Variable %in% names(to_original), to_original[Variable], Variable))
  
  terms_g <- read_excel(out_xlsx, sheet = "dbRDA_Terms") %>%
    filter(Group == g_code, !is.na(p), p < 0.05)
  
  if (nrow(terms_g) > 0) {
    biplot <- biplot %>% filter(Variable %in% terms_g$Term)
  } else if (nrow(biplot) > 0) {
    biplot <- biplot %>%
      mutate(len = sqrt(CAP1^2 + CAP2^2)) %>%
      slice_max(len, n = min(8, n())) %>%
      select(-len)
  }
  
  list(biplot = biplot, terms = terms_g)
}

out_pdf_labels <- file.path(
  out_dir,
  if (energy_adjust) {
    "dbRDA_Biplot_ArrowsOnly_Labels_EnergyAdjusted.pdf"
  } else {
    "dbRDA_Biplot_ArrowsOnly_Labels.pdf"
  }
)

pdf(out_pdf_labels, width = 6.5, height = 5.5)

label_offset <- 0.03
pad_factor <- 1.18
keep_ratio <- 0.92
plot_margin_mm <- c(t = 6, r = 12, b = 6, l = 12)

for (g_code in avail_groups) {
  dat <- prep_biplot(g_code)
  biplot <- dat$biplot
  
  if (nrow(biplot) > 0) {
    biplot <- biplot %>%
      mutate(
        end_x = CAP1,
        end_y = CAP2,
        len = sqrt(end_x^2 + end_y^2),
        ux = ifelse(len > 0, end_x / len, 0),
        uy = ifelse(len > 0, end_y / len, 0),
        raw_lab_x = end_x + label_offset * ux,
        raw_lab_y = end_y + label_offset * uy
      )
  }
  
  if (nrow(biplot) > 0) {
    max_abs_x <- max(abs(c(biplot$end_x, biplot$raw_lab_x)), na.rm = TRUE)
    max_abs_y <- max(abs(c(biplot$end_y, biplot$raw_lab_y)), na.rm = TRUE)
    x_lim_abs <- max_abs_x * pad_factor
    y_lim_abs <- max_abs_y * pad_factor
    x_lims <- c(-x_lim_abs, x_lim_abs)
    y_lims <- c(-y_lim_abs, y_lim_abs)
    
    biplot <- biplot %>%
      mutate(
        lab_x = sign(raw_lab_x) * pmin(abs(raw_lab_x), keep_ratio * x_lim_abs),
        lab_y = sign(raw_lab_y) * pmin(abs(raw_lab_y), keep_ratio * y_lim_abs)
      )
  } else {
    x_lims <- NULL
    y_lims <- NULL
  }
  
  ap <- axis_props_for_group(g_code)
  x_lab <- if (is.na(ap[1])) "CAP1" else paste0("CAP1 (", sprintf("%.1f", ap[1]), "%)")
  y_lab <- if (is.na(ap[2])) "CAP2" else paste0("CAP2 (", sprintf("%.1f", ap[2]), "%)")
  
  p <- ggplot() +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_segment(
      data = biplot,
      aes(x = 0, y = 0, xend = end_x, yend = end_y),
      arrow = arrow(length = unit(0.24, "cm")),
      linewidth = 0.8,
      alpha = 0.95
    ) +
    geom_text(
      data = biplot,
      aes(x = lab_x, y = lab_y, label = orig_label),
      size = 2.6,
      fontface = "plain"
    ) +
    labs(
      title = paste("dbRDA Arrows Only With Labels:", label_map[[g_code]]),
      x = x_lab,
      y = y_lab,
      subtitle = if (nrow(dat$terms) > 0) {
        "Arrows represent significant dietary variables"
      } else {
        "Arrows represent top-loading dietary variables"
      }
    ) +
    theme_minimal(base_size = 12.5) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(color = "grey30"),
      plot.margin = margin(
        plot_margin_mm["t"],
        plot_margin_mm["r"],
        plot_margin_mm["b"],
        plot_margin_mm["l"],
        unit = "mm"
      )
    )
  
  if (!is.null(x_lims) && !is.null(y_lims) && all(is.finite(c(x_lims, y_lims)))) {
    p <- p + coord_cartesian(xlim = x_lims, ylim = y_lims, clip = "off")
  }
  
  print(p)
}

dev.off()

out_pdf_nolabels <- file.path(
  out_dir,
  if (energy_adjust) {
    "dbRDA_Biplot_ArrowsOnly_NoLabels_EnergyAdjusted.pdf"
  } else {
    "dbRDA_Biplot_ArrowsOnly_NoLabels.pdf"
  }
)

pdf(out_pdf_nolabels, width = 6.5, height = 5.5)

pad_factor_nolab <- 1.18

for (g_code in avail_groups) {
  dat <- prep_biplot(g_code)
  biplot <- dat$biplot
  
  if (nrow(biplot) > 0) {
    biplot <- biplot %>%
      mutate(end_x = CAP1, end_y = CAP2)
  }
  
  if (nrow(biplot) > 0) {
    max_abs_x <- max(abs(biplot$end_x), na.rm = TRUE)
    max_abs_y <- max(abs(biplot$end_y), na.rm = TRUE)
    x_lim_abs <- max_abs_x * pad_factor_nolab
    y_lim_abs <- max_abs_y * pad_factor_nolab
    x_lims <- c(-x_lim_abs, x_lim_abs)
    y_lims <- c(-y_lim_abs, y_lim_abs)
  } else {
    x_lims <- NULL
    y_lims <- NULL
  }
  
  ap <- axis_props_for_group(g_code)
  x_lab <- if (is.na(ap[1])) "CAP1" else paste0("CAP1 (", sprintf("%.1f", ap[1]), "%)")
  y_lab <- if (is.na(ap[2])) "CAP2" else paste0("CAP2 (", sprintf("%.1f", ap[2]), "%)")
  
  p <- ggplot() +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_segment(
      data = biplot,
      aes(x = 0, y = 0, xend = end_x, yend = end_y),
      arrow = arrow(length = unit(0.24, "cm")),
      linewidth = 0.9,
      alpha = 0.95
    ) +
    labs(
      title = paste("dbRDA Arrows Only No Labels:", label_map[[g_code]]),
      x = x_lab,
      y = y_lab,
      subtitle = if (nrow(dat$terms) > 0) {
        "Arrows represent significant dietary variables"
      } else {
        "Arrows represent top-loading dietary variables"
      }
    ) +
    theme_minimal(base_size = 12.5) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(color = "grey30")
    )
  
  if (!is.null(x_lims) && !is.null(y_lims) && all(is.finite(c(x_lims, y_lims)))) {
    p <- p + coord_cartesian(xlim = x_lims, ylim = y_lims, clip = "on")
  }
  
  print(p)
}

dev.off()