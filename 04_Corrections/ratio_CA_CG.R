# ============================================================
# mCA windows by mCG bin + correlation statistics with age
# ============================================================

library(dplyr)
library(ggplot2)
library(readr)
library(broom)
library(tidyr)
library(openxlsx)

# ---------- 1. Paths and parameters ----------
# Uncomment the desired dataset block and comment the rest

# --- EM-seq (human) ---
# data_dir <- "/home/epi-lab/Escritorio/Resultados/CHH/"
# samples <- list(
#   list(file = "35y_cowindows.txt",  age = 35),
#   list(file = "43y_cowindows.txt",  age = 43),
#   list(file = "55y_cowindows.txt",  age = 55),
#   list(file = "67y_cowindows.txt",  age = 67),
#   list(file = "78y_cowindows.txt",  age = 78),
#   list(file = "85y_cowindows.txt",  age = 85),
#   list(file = "103y_cowindows.txt", age = 103)
# )

# --- WGBS (human) ---
 data_dir <- "/home/epi-lab/Escritorio/WGBS/CHH/"
 samples <- list(
   list(file = "cowindows_M9.txt",  age = 5),
   list(file = "cowindows_M8.txt",  age = 12),
   list(file = "cowindows_M7.txt",  age = 16),
   list(file = "cowindows_M6.txt",  age = 25),
   list(file = "cowindows_M5.txt",  age = 49),
   list(file = "cowindows_M4.txt",  age = 64),
   list(file = "cowindows_M3.txt",  age = 77),
   list(file = "cowindows_M2.txt",  age = 77),
   list(file = "cowindows_M1.txt",  age = 81)
 )

# --- Other species ---
# data_dir <- "/home/epi-lab/Escritorio/Animals/CHH/animals/"
# samples <- list(
#   list(file = "cowindows_ornitorrinco.txt", age = 35),
#   list(file = "cowindows_rana.txt",         age = 40),
#   list(file = "cowindows_pollo.txt",        age = 43),
#   list(file = "cowindows_zebrafish.txt",    age = 44),
#   list(file = "cowindows_tiburon.txt",      age = 45)
# )

# --- Chimpanzee ---
# data_dir <- "/home/epi-lab/Escritorio/Animals/CHH/chimpance/"
# samples <- list(
#   list(file = "cowindows_35y.txt", age = 35),
#   list(file = "cowindows_40y.txt", age = 40),
#   list(file = "cowindows_43y.txt", age = 43),
#   list(file = "cowindows_44y.txt", age = 44),
#   list(file = "cowindows_45y.txt", age = 45),
#   list(file = "cowindows_57y.txt", age = 57)
# )

 
 bin_size               <- 10
 mCA_quantile_threshold <- 0.66
 
# ---------- 2. Load and merge all samples ----------
df_all <- lapply(seq_along(samples), function(i) {
  s <- samples[[i]]
  read_delim(file.path(data_dir, s$file), delim = "\t", show_col_types = FALSE) %>%
    mutate(
      chr    = as.character(chr),
      age    = s$age,
      sample = paste0("m", i)
    )
}) %>% bind_rows()

# ---------- 3. Bin mCG values ----------
breaks_vec <- seq(0, 100, by = bin_size)

df_bin <- df_all %>%
  mutate(CG_bin = cut(CG_bvalues, breaks = breaks_vec, include.lowest = TRUE))

# ---------- 4. Background: total windows per bin and sample ----------
df_bg <- df_bin %>%
  group_by(age, sample, CG_bin) %>%
  summarise(n_windows_bin = n(), .groups = "drop")

# ---------- 5. Filter high-mCA windows per sample (above quantile threshold) ----------
df_Qfilter <- df_bin %>%
  group_by(sample) %>%
  filter(CA_bvalues > quantile(CA_bvalues, mCA_quantile_threshold, na.rm = TRUE)) %>%
  ungroup()

# ---------- 6. Count high-mCA windows per bin ----------
df_counts_raw <- df_Qfilter %>%
  group_by(age, sample, CG_bin) %>%
  summarise(n_events = n(), .groups = "drop")

# ---------- 7. Join with background and compute rate ----------
df_counts <- df_counts_raw %>%
  right_join(df_bg, by = c("age", "sample", "CG_bin")) %>%
  mutate(
    n_events    = ifelse(is.na(n_events), 0L, n_events),
    rate_in_bin = n_events / n_windows_bin
  )

# ---------- 8. Compute bin center ----------
df_counts <- df_counts %>%
  mutate(
    bin_index     = as.integer(CG_bin),
    CG_bin_center = (breaks_vec[bin_index] + breaks_vec[bin_index + 1]) / 2
  ) %>%
  dplyr::select(-bin_index)

# ---------- 9. Normalize per sample (sum = 1) ----------
df_counts <- df_counts %>%
  group_by(sample) %>%
  mutate(rate_norm = rate_in_bin / sum(rate_in_bin, na.rm = TRUE)) %>%
  ungroup()

# ---------- 10. Color palette (age gradient) ----------
ages_sorted  <- sort(unique(df_counts$age))
blue_palette <- setNames(
  colorRampPalette(c("#c6dbef", "#035887"))(length(ages_sorted)),
  as.character(ages_sorted)
)

# ============================================================
# MAIN PLOT - normalized rate per sample
# ============================================================

plot_main <- ggplot(df_counts,
                    aes(x     = CG_bin_center,
                        y     = rate_norm,
                        color = as.factor(age),
                        group = sample)) +
  geom_smooth(method = "loess", se = FALSE, span = 1, linewidth = 1) +
  scale_x_continuous(breaks = breaks_vec, limits = c(0, 100)) +
  scale_color_manual(values = blue_palette) +
  labs(
    x     = "mCG (β-value)",
    y     = "mCA/mCG β-value",
    color = "Age\n(years)"
  ) +
  theme_classic(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 18, color = "black"),
    axis.title.y     = element_text(size = 18),
    axis.text.x      = element_text(size = 18, color = "black"),
    axis.title.x     = element_text(size = 18),
    axis.ticks       = element_line(size = 1),
    axis.line        = element_line(size = 1),
    legend.text      = element_text(size = 16),
    legend.title     = element_text(size = 18)
  )
print(plot_main)

# ============================================================
# STATISTICS 1 - Linear regression (LM): rate_norm ~ age per bin
# ============================================================

# --- Slope and intercept ---
df_tidy <- df_counts %>%
  group_by(CG_bin_center) %>%
  filter(n() >= 3) %>%
  group_modify(~ tidy(lm(rate_norm ~ age, data = .x))) %>%
  ungroup() %>%
  filter(term == "age") %>%
  dplyr::select(
    CG_bin_center,
    slope   = estimate,
    SE      = std.error,
    t_stat  = statistic,
    p_value = p.value
  )

# --- R² ---
df_glance <- df_counts %>%
  group_by(CG_bin_center) %>%
  filter(n() >= 3) %>%
  group_modify(~ glance(lm(rate_norm ~ age, data = .x))) %>%
  ungroup() %>%
  dplyr::select(CG_bin_center, R2 = r.squared, R2_adj = adj.r.squared)

# ============================================================
# STATISTICS 2 - Spearman rho: rate_norm ~ age per bin
# ============================================================

df_spearman_by_bin <- df_counts %>%
  group_by(CG_bin_center) %>%
  filter(n() >= 3) %>%
  summarise(
    rho     = cor.test(age, rate_norm, method = "spearman", exact = FALSE)$estimate,
    p_value = cor.test(age, rate_norm, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  )

print(df_spearman_by_bin)

# ============================================================
# STATISTICS 3 - Spearman rho within each sample
# Correlation between CG_bin_center and rate_norm
# ============================================================

df_spearman_per_sample <- df_counts %>%
  group_by(sample, age) %>%
  filter(n() >= 3) %>%
  summarise(
    rho     = cor.test(CG_bin_center, rate_norm, method = "spearman", exact = FALSE)$estimate,
    p_value = cor.test(CG_bin_center, rate_norm, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  arrange(age)

print(df_spearman_per_sample)

# ============================================================
# STATISTICS 4 - Pearson (LM) within each sample
# Correlation between CG_bin_center and rate_norm
# ============================================================

df_pearson_per_sample <- df_counts %>%
  group_by(sample, age) %>%
  filter(n() >= 3) %>%
  group_modify(~ {
    fit <- lm(rate_norm ~ CG_bin_center, data = .x)
    s   <- summary(fit)
    tibble(
      slope   = coef(fit)["CG_bin_center"],
      SE      = s$coefficients["CG_bin_center", "Std. Error"],
      R       = sign(coef(fit)["CG_bin_center"]) * sqrt(s$r.squared),
      R2      = s$r.squared,
      t_stat  = s$coefficients["CG_bin_center", "t value"],
      p_value = s$coefficients["CG_bin_center", "Pr(>|t|)"]
    )
  }) %>%
  arrange(age)

print(df_pearson_per_sample)

# ============================================================
# FOLD CHANGE - bin 5 (0-10) vs bin 95 (90-100) ~ age
# ============================================================

df_fc <- df_counts %>%
  filter(CG_bin_center %in% c(5, 95)) %>%
  dplyr::select(sample, age, CG_bin_center, rate_norm) %>%
  tidyr::pivot_wider(
    names_from   = CG_bin_center,
    values_from  = rate_norm,
    names_prefix = "bin_"
  ) %>%
  mutate(
    FC     = bin_95 / bin_5,
    log2FC = log2(FC)
  ) %>%
  arrange(age)

print(df_fc)

# --- FC ~ age correlation and LM ---
spearman_fc <- cor.test(df_fc$age, df_fc$FC, method = "spearman", exact = FALSE)
fit_fc      <- lm(FC ~ age, data = df_fc)
s_fc        <- summary(fit_fc)

cat("\n--- Spearman FC ~ age ---\n")
cat("rho =", round(spearman_fc$estimate, 3),
    "| p =", round(spearman_fc$p.value, 4), "\n")
cat("\n--- LM slope FC ~ age ---\n")
cat("slope =", round(coef(fit_fc)["age"], 6),
    "| R2 =", round(s_fc$r.squared, 3),
    "| p =", round(s_fc$coefficients["age", "Pr(>|t|)"], 4), "\n")

# ============================================================
# PLOT - FC per sample with regression line
# ============================================================

ggplot(df_fc, aes(x = age, y = FC)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  geom_point(aes(fill = as.factor(age)), shape = 21, size = 4, stroke = 0.5) +
  scale_fill_manual(values = blue_palette) +
  labs(
    x    = "Age (years)",
    y    = "Fold Change (mCA/mCG β-value)",
    fill = "Age\n(years)"
  ) +
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
  theme_classic(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 18, color = "black"),
    axis.title.y     = element_text(size = 18),
    axis.text.x      = element_text(size = 18, color = "black"),
    axis.title.x     = element_text(size = 18),
    axis.ticks       = element_line(size = 1),
    axis.line        = element_line(size = 1),
    legend.text      = element_text(size = 16),
    legend.title     = element_text(size = 18),
    legend.position  = "none"
  )

# ============================================================
# EXPORT ALL STATISTICS TO EXCEL
# ============================================================

# --- Merge LM tables ---
df_lm_by_bin <- df_tidy %>%
  left_join(df_glance, by = "CG_bin_center") %>%
  mutate(R = sign(slope) * sqrt(R2)) %>%
  dplyr::select(CG_bin_center, slope, SE, R, R2, R2_adj, t_stat, p_value)

# --- Styles ---
wb <- createWorkbook()
addWorksheet(wb, "Results")

style_title <- createStyle(
  fontName = "Arial", fontSize = 12,
  fontColour = "#FFFFFF", fgFill = "#1F4E79",
  textDecoration = "bold", halign = "left",
  border = "TopBottomLeftRight", borderColour = "#1F4E79"
)
style_header <- createStyle(
  fontName = "Arial", fontSize = 10,
  fontColour = "#FFFFFF", fgFill = "#2E75B6",
  textDecoration = "bold", halign = "center",
  border = "TopBottomLeftRight", borderColour = "#BFBFBF"
)
style_body <- createStyle(
  fontName = "Arial", fontSize = 10,
  halign = "center", numFmt = "0.000000",
  border = "TopBottomLeftRight", borderColour = "#BFBFBF"
)

# --- Helper function to write each table block ---
write_block <- function(wb, sheet, df, title, start_row) {
  n_rows <- nrow(df)
  n_cols <- ncol(df)
  
  writeData(wb, sheet, x = title, startRow = start_row, startCol = 1)
  mergeCells(wb, sheet, cols = 1:n_cols, rows = start_row)
  addStyle(wb, sheet, style_title, rows = start_row, cols = 1:n_cols, gridExpand = TRUE)
  
  writeData(wb, sheet, x = df,
            startRow = start_row + 1, startCol = 1,
            headerStyle = style_header,
            borders = "all", borderColour = "#BFBFBF")
  
  addStyle(wb, sheet, style_body,
           rows = (start_row + 2):(start_row + 1 + n_rows),
           cols = 1:n_cols, gridExpand = TRUE)
  
  return(start_row + 1 + n_rows + 2)
}

# --- Write all tables sequentially ---
next_row <- 1
next_row <- write_block(wb, "Results", df_lm_by_bin,
                        "1. Pearson (LM): rate_norm ~ age per mCG bin",
                        next_row)
next_row <- write_block(wb, "Results", df_spearman_by_bin,
                        "2. Spearman rho: rate_norm ~ age per mCG bin",
                        next_row)
next_row <- write_block(wb, "Results", df_spearman_per_sample,
                        "3. Spearman rho intra-sample: CG_bin_center ~ rate_norm",
                        next_row)
next_row <- write_block(wb, "Results", df_pearson_per_sample,
                        "4. Pearson (LM) intra-sample: CG_bin_center ~ rate_norm",
                        next_row)

# --- Column widths and save ---
setColWidths(wb, "Results", cols = 1:12, widths = "auto")
out_xlsx <- file.path(data_dir, "statistics_mCA_mCG.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat("Excel saved at:", out_xlsx, "\n")