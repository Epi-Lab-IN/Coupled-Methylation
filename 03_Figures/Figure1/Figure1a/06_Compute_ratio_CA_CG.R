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
