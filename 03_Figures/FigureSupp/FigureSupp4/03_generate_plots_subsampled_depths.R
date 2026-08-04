## 03 — generate plots of shared-depth individual dashboard
##
## Produces the data tables and component/final images directly associated
## with the individual shared-depth dashboard.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

## ---- User configuration: edit paths before running --------------------------
# input data
subbam_metrics_file <- "/path/to/02_results/tables/all_samples_CoM_metrics_by_bam_depth_and_mincov.tsv"
# output directory
outdir <- "/path/to/03_results"
# minimum coverage thresholds
overlay_dashboard_mincov_values <- c(1L, 5L, 10L)
overlay_mincov_colors <- c("1"="#0072B2", "5"="#D55E00", "10"="#009E73")

# sample names and ages
sample_map <- data.table(
  sample=paste0("M", 1:7),
  age=c(35L, 43L, 55, 67L, 78L, 85, 103L)
)

## ---- Directories and run record ---------------------------------------------
table_dir <- file.path(outdir, "tables")
plot_dir <- file.path(outdir, "plots")
log_dir <- file.path(outdir, "logs")
for (d in c(table_dir, plot_dir, log_dir)) dir.create(d, recursive=TRUE, showWarnings=FALSE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_start.txt"))
fwrite(data.table(parameter=c("subbam_metrics_file", "outdir", "overlay_dashboard_mincov_values"),
                  value=c(subbam_metrics_file, outdir, paste(overlay_dashboard_mincov_values, collapse=","))),
       file.path(log_dir, "run_configuration.tsv"), sep="\t")

format_depth_x <- function(x) paste0(formatC(x, format="f", digits=0), "x")
n_expected <- nrow(sample_map)

## ---- Load metrics required for this dashboard --------------------------
if (!file.exists(subbam_metrics_file)) stop("Missing metric table: ", subbam_metrics_file)
metrics <- fread(subbam_metrics_file)
required_columns <- c("sample", "age", "depth_x", "mincov_cg", "mincov_ca", "pair_type", "n_windows", "slope", "mCA_mean")
missing_columns <- setdiff(required_columns, names(metrics))
if (length(missing_columns)) stop("Metric table lacks required columns: ", paste(missing_columns, collapse=", "))

overlay_dt <- metrics[
  pair_type == "matched_min_coverage" &
    mincov_cg == mincov_ca &
    mincov_cg %in% overlay_dashboard_mincov_values &
    is.finite(depth_x) & is.finite(slope) & is.finite(mCA_mean) & is.finite(n_windows)
]
overlay_dt <- copy(overlay_dt)
overlay_dt[, sample := factor(as.character(sample), levels=sample_map$sample)]
overlay_dt[, mincov_factor := factor(as.character(mincov_cg), levels=as.character(overlay_dashboard_mincov_values))]
overlay_dt[, age_numeric := as.numeric(age)]
overlay_dt[, mCA_B_value := mCA_mean / 100]
overlay_dt[, n_CG1CA3_windows := as.numeric(n_windows)]
if (!nrow(overlay_dt)) stop("No finite rows remain for the requested matched minimum-coverage thresholds.")

## Determine depths represented by every sample at every plotted mincov.
depth_coverage <- overlay_dt[, .(n_samples=uniqueN(sample), samples=paste(sort(unique(as.character(sample))), collapse=",")),
                             by=.(mincov_cg, depth_x)][order(mincov_cg, depth_x)]
shared_candidates <- depth_coverage[n_samples == n_expected,
  .(n_mincov_with_full_samples=uniqueN(mincov_cg), mincov_values_with_full_samples=paste(sort(unique(mincov_cg)), collapse=",")),
  by=depth_x]
shared_depths <- shared_candidates[n_mincov_with_full_samples == length(overlay_dashboard_mincov_values), sort(depth_x)]
if (!length(shared_depths)) stop("No depths are shared across all samples and all plotted minimum-coverage thresholds.")
overlay_dt <- overlay_dt[depth_x %in% shared_depths]

## Pearson association across samples at each mincov/depth.
overlay_age_slope <- overlay_dt[, {
  pear <- if (.N >= 3L && stats::sd(slope, na.rm=TRUE) > 0) tryCatch(stats::cor.test(age_numeric, slope, method="pearson"), error=function(e) NULL) else NULL
  .(n_samples=uniqueN(sample), age_slope_pearson_r=if (is.null(pear)) NA_real_ else as.numeric(pear$estimate),
    age_slope_pearson_p=if (is.null(pear)) NA_real_ else pear$p.value)
}, by=.(mincov_cg, mincov_factor, depth_x)][order(mincov_cg, depth_x)]
overlay_age_slope[, sig := fifelse(is.finite(age_slope_pearson_p) & age_slope_pearson_p < 0.05, "p < 0.05", "p >= 0.05")]

## ---- Direct supporting tables ------------------------------------------------
fwrite(overlay_dt, file.path(table_dir, "overlay_dashboard_individual_metrics_sharedDepths.tsv"), sep="\t")
fwrite(depth_coverage, file.path(table_dir, "overlay_dashboard_depth_coverage_by_mincov.tsv"), sep="\t")
fwrite(data.table(depth_x=shared_depths), file.path(table_dir, "overlay_dashboard_shared_depths_used.tsv"), sep="\t")
fwrite(overlay_age_slope, file.path(table_dir, "overlay_dashboard_age_slope_by_mincov_sharedDepths.tsv"), sep="\t")

## ---- Plot setup --------------------------------------------------------------
x_breaks <- sort(shared_depths)
mincov_cols <- overlay_mincov_colors[as.character(overlay_dashboard_mincov_values)]
mincov_cols <- mincov_cols[!is.na(mincov_cols)]
overlay_y_grid_theme <- theme(panel.grid.major.y=element_line(color="grey86", linewidth=0.25),
                              panel.grid.major.x=element_blank(), panel.grid.minor=element_blank())

p_mCA <- ggplot(overlay_dt, aes(x=depth_x, y=mCA_B_value, color=mincov_factor, alpha=age_numeric, group=interaction(sample,mincov_factor))) +
  geom_line(linewidth=0.72) + geom_point(size=1.9) +
  scale_x_continuous(breaks=x_breaks, labels=format_depth_x) +
  scale_y_continuous(labels=percent_format(accuracy=0.1)) +
  scale_color_manual(values=mincov_cols, name="Minimum cytosine coverage") +
  scale_alpha_continuous(name="Age", range=c(0.2,1.0), breaks=sample_map$age) +
  labs(title="Mean mCA B-value", x=NULL, y="Mean mCA (β-value)") + theme_classic(base_size=13) + theme(legend.position="right") + overlay_y_grid_theme

p_slope <- ggplot(overlay_dt, aes(x=depth_x, y=slope, color=mincov_factor, alpha=age_numeric, group=interaction(sample,mincov_factor))) +
  geom_hline(yintercept=0, linetype="dashed", color="grey55") + geom_line(linewidth=0.72) + geom_point(size=1.9) +
  scale_x_continuous(breaks=x_breaks, labels=format_depth_x) +
  scale_color_manual(values=mincov_cols, name="Minimum cytosine coverage") +
  scale_alpha_continuous(name="Age", range=c(0.2,1.0), breaks=sample_map$age) +
  labs(title="CoM slope", x=NULL, y="CoM slope\n(mCA M-value ~ mCG M-value)") + theme_classic(base_size=13) + theme(legend.position="right") + overlay_y_grid_theme

p_windows <- ggplot(overlay_dt, aes(x=depth_x, y=n_CG1CA3_windows, color=mincov_factor, alpha=age_numeric, group=interaction(sample,mincov_factor))) +
  geom_line(linewidth=0.72) + geom_point(size=1.9) +
  scale_x_continuous(breaks=x_breaks, labels=format_depth_x) + scale_y_log10(labels=comma) +
  scale_color_manual(values=mincov_cols, name="Minimum cytosine coverage") +
  scale_alpha_continuous(name="Age", range=c(0.2,1.0), breaks=sample_map$age) +
  labs(title="# of CG1 ∩ CA3 windows", x=NULL, y="# of CG1 ∩ CA3 windows") + theme_classic(base_size=13) + theme(legend.position="right") + overlay_y_grid_theme

p_age_slope <- ggplot(overlay_age_slope, aes(x=depth_x, y=age_slope_pearson_r, color=mincov_factor, group=mincov_factor)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey55") + geom_line(linewidth=0.9) + geom_point(aes(shape=sig), size=2.8, stroke=1.0) +
  scale_x_continuous(breaks=x_breaks, labels=format_depth_x) + coord_cartesian(ylim=c(-1,1)) +
  scale_color_manual(values=mincov_cols, name="Minimum cytosine coverage") +
  scale_shape_manual(values=c("p < 0.05"=8, "p >= 0.05"=1), name="age~slope") +
  labs(title="Age~slope association", x="Nominal sequencing depth", y="age~slope\nPearson r") +
  theme_classic(base_size=13) + theme(legend.position="right") + overlay_y_grid_theme + guides(color="none")

## Component images and target dashboard.
ggsave(file.path(plot_dir, "03_overlay_individual_mCA_sharedDepths.png"), p_mCA, width=9.5, height=5.5, dpi=300)
ggsave(file.path(plot_dir, "03_overlay_individual_CoM_slope_sharedDepths.png"), p_slope, width=9.5, height=5.5, dpi=300)
ggsave(file.path(plot_dir, "03_overlay_individual_CG1capCA3_windows_sharedDepths.png"), p_windows, width=9.5, height=5.5, dpi=300)
ggsave(file.path(plot_dir, "03_overlay_age_slope_sharedDepths.png"), p_age_slope, width=9.5, height=5.5, dpi=300)

dashboard <- (p_mCA / p_slope / p_windows / p_age_slope) + plot_layout(guides="collect") & theme(legend.position="right")
ggsave(file.path(plot_dir, "03_dashboard_overlay_individual_sharedDepths_mCA_slope_windows_ageSlope.png"), dashboard, width=11.8, height=12.8, dpi=300)

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo_end.txt"))
message("03 complete: ", file.path(plot_dir, "03_dashboard_overlay_individual_sharedDepths_mCA_slope_windows_ageSlope.png"))
