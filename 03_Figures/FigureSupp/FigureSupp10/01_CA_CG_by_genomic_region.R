# ========================================
# METHYLATION ANALYSIS BATCH: CG and CA
# Processes all samples in a folder
# Extracts data by genomic region
# ========================================

# Load libraries
library(ggplot2)
library(biomaRt)
library(dplyr)
library(methylKit)
library(genomation)
library(stringr)
library(readr)

# ========================================
# SETTINGS
# ========================================
resultado_dir <- "/home/epi-lab/Escritorio/Resultados/CHH"
database_dir <- "/home/epi-lab/Escritorio/Database"
output_lista <- list()

# ========================================
# LOAD GENE ANNOTATION DATABASE
# ========================================
gene.obj <- readTranscriptFeatures(
  file.path(database_dir, "hg38.bed"),
  remove.unusual = TRUE
)

# ========================================
# LOAD COMMON SITES (FILTER)
# ========================================

# Common CG sites
cg_comun_file <- file.path(resultado_dir, "100bqW_M7_CG_min3")
if (file.exists(cg_comun_file)) {
  cg_comun <- read_delim(cg_comun_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                         quote = "\"",
                         col_types = cols(chr = col_character(), start = col_double(), 
                                          end = col_double(), coverage = col_double(),
                                          numCs = col_double(), numTs = col_double(), 
                                          .default = col_character()))
  cg_comun_GR <- makeGRangesFromDataFrame(cg_comun, keep.extra.columns = TRUE)
} else {
  cg_comun_GR <- NULL
}

# Common CA/CHH sites
ca_comun_file <- file.path(resultado_dir, "CA_comunes_CHH")
if (file.exists(ca_comun_file)) {
  ca_comun <- read_delim(ca_comun_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                         quote = "\"",
                         col_types = cols(chr = col_character(), start = col_double(), 
                                          end = col_double(), coverage = col_double(),
                                          numCs = col_double(), numTs = col_double(), 
                                          .default = col_character()))
  ca_comun_GR <- makeGRangesFromDataFrame(ca_comun, keep.extra.columns = TRUE)
} else {
  ca_comun_GR <- NULL
}

# ========================================
# DETECT SAMPLES AUTOMATICALLY
# ========================================
# Look for CG files (100bqW_M1_CG_min3, 100bqW_M2_CG_min3, etc.)
archivos_cg <- list.files(resultado_dir, pattern = "^100bqW_M[0-9]_CG_min3$")
prefijos <- unique(gsub("^100bqW_", "", gsub("_CG_min3$", "", archivos_cg)))

# Keep only M1 to M7
prefijos <- prefijos[grepl("^M[1-7]$", prefijos)]
prefijos <- sort(prefijos)

# ========================================
# LOAD ENHANCERS (OPTIONAL)
# ========================================
enhancers_file <- file.path(database_dir, "H3K27ac_HS_64.bed")
enhancers_GR <- NULL

if (file.exists(enhancers_file)) {
  enhancers <- read_delim(enhancers_file, delim = "\t", escape_double = FALSE, 
                          col_names = TRUE,
                          trim_ws = TRUE)
  enhancers$chr <- str_replace(enhancers$chr, "chr", "")
  enhancers_GR <- GRanges(seqnames = enhancers$chr, 
                          ranges = IRanges(start = enhancers$start, end = enhancers$end))
}

# Annotate the common sites once for all samples

# Annotate all common CG sites at once
if (!is.null(cg_comun_GR)) {
  ann_comun_cg <- annotateWithGeneParts(cg_comun_GR, gene.obj)
  class_comun_cg <- getMembers(ann_comun_cg)
  
  # Get indices per region
  idx_promo_cg <- which(class_comun_cg[, 1] == 1)
  idx_exon_cg <- which(class_comun_cg[, 1] == 0 & class_comun_cg[, 2] == 1)
  idx_intron_cg <- which(class_comun_cg[, 1] == 0 & class_comun_cg[, 2] == 0 & class_comun_cg[, 3] == 1)
  idx_inter_cg <- which(class_comun_cg[, 1] == 0 & class_comun_cg[, 2] == 0 & class_comun_cg[, 3] == 0)
}

# Annotate all common CA/CHH sites at once
if (!is.null(ca_comun_GR)) {
  ann_comun_ca <- annotateWithGeneParts(ca_comun_GR, gene.obj)
  class_comun_ca <- getMembers(ann_comun_ca)
  
  # Get indices per region
  idx_promo_ca <- which(class_comun_ca[, 1] == 1)
  idx_exon_ca <- which(class_comun_ca[, 1] == 0 & class_comun_ca[, 2] == 1)
  idx_intron_ca <- which(class_comun_ca[, 1] == 0 & class_comun_ca[, 2] == 0 & class_comun_ca[, 3] == 1)
  idx_inter_ca <- which(class_comun_ca[, 1] == 0 & class_comun_ca[, 2] == 0 & class_comun_ca[, 3] == 0)
}

# ========================================
# FUNCTION: process one sample
# ========================================
procesar_muestra <- function(prefijo, resultado_dir, gene.obj, cg_comun_GR = NULL, ca_comun_GR = NULL, enhancers_GR = NULL) {
  
  # Initialize enhancer variables
  enhancers_cg <- data.frame()
  enhancers_ca <- data.frame()
  
  # -------- CG DATA --------
  cpg_file <- file.path(resultado_dir, paste0("100bqW_", prefijo, "_CG_min3"))
  if (!file.exists(cpg_file)) {
    return(NULL)
  }
  
  cpg <- read_delim(cpg_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                    quote = "\"",
                    col_types = cols(chr = col_character(), start = col_double(), 
                                     end = col_double(), coverage = col_double(),
                                     numCs = col_double(), numTs = col_double(), 
                                     .default = col_character()))
  cpg <- makeGRangesFromDataFrame(cpg, keep.extra.columns = TRUE)
  
  # Filter by common CG sites if available
  if (!is.null(cg_comun_GR)) {
    cpg <- cpg[cpg %over% cg_comun_GR]
  }
  
  cpg_df <- as.data.frame(cpg)
  cpg_df$width <- NULL
  colnames(cpg_df)[1] <- "chr"
  cpg_df$Metilacion_CG <- (cpg_df$numCs / cpg_df$coverage) * 100
  cpg_df$mvalues_cg <- log2((cpg_df$numCs + 0.1) / (cpg_df$numTs + 0.1))
  
  # Use precomputed indices (much faster)
  cpg_GR <- as(cpg_df, "GRanges")
  
  if (!is.null(cg_comun_GR) && exists("idx_promo_cg")) {
    hits_promo <- findOverlaps(cpg_GR, cg_comun_GR[idx_promo_cg], select = "first")
    hits_exon <- findOverlaps(cpg_GR, cg_comun_GR[idx_exon_cg], select = "first")
    hits_intron <- findOverlaps(cpg_GR, cg_comun_GR[idx_intron_cg], select = "first")
    hits_inter <- findOverlaps(cpg_GR, cg_comun_GR[idx_inter_cg], select = "first")
    
    # Extract by region
    promotores_cg <- cpg_df[!is.na(hits_promo), ]
    exones_cg <- cpg_df[!is.na(hits_exon), ]
    intrones_cg <- cpg_df[!is.na(hits_intron), ]
    intergenico_cg <- cpg_df[!is.na(hits_inter), ]
  } else {
    # No precomputed annotation, annotate this sample directly
    ann_cpg <- annotateWithGeneParts(cpg_GR, gene.obj)
    class_cpg <- getMembers(ann_cpg)
    
    promotores_cg <- cpg_df[class_cpg[, 1] == 1, ]
    exones_cg <- cpg_df[class_cpg[, 1] == 0 & class_cpg[, 2] == 1, ]
    intrones_cg <- cpg_df[class_cpg[, 1] == 0 & class_cpg[, 2] == 0 & class_cpg[, 3] == 1, ]
    intergenico_cg <- cpg_df[class_cpg[, 1] == 0 & class_cpg[, 2] == 0 & class_cpg[, 3] == 0, ]
  }
  
  # Extract enhancers for CG
  if (!is.null(enhancers_GR)) {
    enhancers_hits <- findOverlaps(cpg_GR, enhancers_GR, select = "first")
    enhancers_cg <- cpg_df[!is.na(enhancers_hits), ]
  }
  
  # -------- CA/CHH DATA --------
  ca_file <- file.path(resultado_dir, paste0(prefijo, "_CA_EMseq_CHH_100bpW"))
  promotores_ca <- exones_ca <- intrones_ca <- intergenico_ca <- NULL
  
  if (file.exists(ca_file)) {
    ca <- read_delim(ca_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                     quote = "\"",
                     col_types = cols(chr = col_character(), start = col_double(), 
                                      end = col_double(), coverage = col_double(),
                                      numCs = col_double(), numTs = col_double(), 
                                      .default = col_character()))
    ca_GR <- makeGRangesFromDataFrame(ca, keep.extra.columns = TRUE)
    
    # Filter by common CA sites if available
    if (!is.null(ca_comun_GR)) {
      ca_GR <- ca_GR[ca_GR %over% ca_comun_GR]
    }
    
    ca <- as.data.frame(ca_GR)
    ca$Metilacion_CA <- (ca$numCs / ca$coverage) * 100
    ca$mvalues_ca <- log2((ca$numCs + 0.1) / (ca$numTs + 0.1))
    
    ca_df <- as.data.frame(ca_GR)
    ca_df$width <- NULL
    colnames(ca_df)[1] <- "chr"
    ca_df$Metilacion_CA <- (ca_df$numCs / ca_df$coverage) * 100
    ca_df$mvalues_ca <- log2((ca_df$numCs + 0.1) / (ca_df$numTs + 0.1))
    
    # Use precomputed indices if they exist
    if (!is.null(ca_comun_GR) && exists("idx_promo_ca")) {
      hits_promo_ca <- findOverlaps(ca_GR, ca_comun_GR[idx_promo_ca], select = "first")
      hits_exon_ca <- findOverlaps(ca_GR, ca_comun_GR[idx_exon_ca], select = "first")
      hits_intron_ca <- findOverlaps(ca_GR, ca_comun_GR[idx_intron_ca], select = "first")
      hits_inter_ca <- findOverlaps(ca_GR, ca_comun_GR[idx_inter_ca], select = "first")
      
      promotores_ca <- ca_df[!is.na(hits_promo_ca), ]
      exones_ca <- ca_df[!is.na(hits_exon_ca), ]
      intrones_ca <- ca_df[!is.na(hits_intron_ca), ]
      intergenico_ca <- ca_df[!is.na(hits_inter_ca), ]
    } else {
      # No precomputed annotation, annotate this sample directly
      ann_ca <- annotateWithGeneParts(ca_GR, gene.obj)
      class_ca <- getMembers(ann_ca)
      
      promotores_ca <- ca_df[class_ca[, 1] == 1, ]
      exones_ca <- ca_df[class_ca[, 1] == 0 & class_ca[, 2] == 1, ]
      intrones_ca <- ca_df[class_ca[, 1] == 0 & class_ca[, 2] == 0 & class_ca[, 3] == 1, ]
      intergenico_ca <- ca_df[class_ca[, 1] == 0 & class_ca[, 2] == 0 & class_ca[, 3] == 0, ]
    }
    
    # Extract enhancers for CA
    enhancers_ca <- data.frame()
    if (!is.null(enhancers_GR)) {
      enhancers_hits_ca <- findOverlaps(ca_GR, enhancers_GR, select = "first")
      enhancers_ca <- ca_df[!is.na(enhancers_hits_ca), ]
    }
  }
  
  # -------- COMPILE RESULTS --------
  # Check that there is data
  if (nrow(promotores_cg) == 0 && nrow(exones_cg) == 0 && 
      nrow(intrones_cg) == 0 && nrow(intergenico_cg) == 0) {
    return(NULL)
  }
  
  resultado <- data.frame(
    Muestra = prefijo,
    Tipo = c("Promotores", "Exones", "Intrones", "Intergénico", "Enhancers"),
    
    Metilacion_CG_mean = c(
      if(nrow(promotores_cg) > 0) mean(promotores_cg$Metilacion_CG) else NA,
      if(nrow(exones_cg) > 0) mean(exones_cg$Metilacion_CG) else NA,
      if(nrow(intrones_cg) > 0) mean(intrones_cg$Metilacion_CG) else NA,
      if(nrow(intergenico_cg) > 0) mean(intergenico_cg$Metilacion_CG) else NA,
      if(nrow(enhancers_cg) > 0) mean(enhancers_cg$Metilacion_CG) else NA
    ),
    
    Conteo_CG = c(
      nrow(promotores_cg), nrow(exones_cg), 
      nrow(intrones_cg), nrow(intergenico_cg),
      nrow(enhancers_cg)
    ),
    
    Metilacion_CA_mean = c(
      if(!is.null(promotores_ca) & nrow(promotores_ca) > 0) mean(promotores_ca$Metilacion_CA) else NA,
      if(!is.null(exones_ca) & nrow(exones_ca) > 0) mean(exones_ca$Metilacion_CA) else NA,
      if(!is.null(intrones_ca) & nrow(intrones_ca) > 0) mean(intrones_ca$Metilacion_CA) else NA,
      if(!is.null(intergenico_ca) & nrow(intergenico_ca) > 0) mean(intergenico_ca$Metilacion_CA) else NA,
      if(!is.null(enhancers_ca) & nrow(enhancers_ca) > 0) mean(enhancers_ca$Metilacion_CA) else NA
    ),
    
    Conteo_CA = c(
      if(!is.null(promotores_ca)) nrow(promotores_ca) else 0,
      if(!is.null(exones_ca)) nrow(exones_ca) else 0,
      if(!is.null(intrones_ca)) nrow(intrones_ca) else 0,
      if(!is.null(intergenico_ca)) nrow(intergenico_ca) else 0,
      if(!is.null(enhancers_ca)) nrow(enhancers_ca) else 0
    )
  )
  
  return(resultado)
}

# ========================================
# RUN BATCH ANALYSIS
# ========================================
resultados_totales <- do.call(
  rbind, 
  Filter(Negate(is.null), lapply(prefijos, procesar_muestra, resultado_dir, gene.obj, cg_comun_GR, ca_comun_GR, enhancers_GR))
)
rownames(resultados_totales) <- NULL

# ========================================
# SUMMARY: RESULTS BY SAMPLE
# ========================================
print(resultados_totales)

# ========================================
# EXPORT RESULTS (OPTIONAL)
# ========================================
write.csv(resultados_totales, 
          "metilacion_emseq_resultados_completos.csv", 
          row.names = FALSE)