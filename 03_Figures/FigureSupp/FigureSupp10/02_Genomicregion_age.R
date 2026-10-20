library(dplyr)
library(tidyr)

# ========================================
# FUNCTION: linear regression of methylation vs age
# per genomic region (Tipo) and sequence context (CG / CA)
# ========================================
age_methylation_regression <- function(meth_data, ages) {
  
  meth_data %>%
    # Add age and reshape CG/CA to long format
    left_join(ages, by = "Muestra") %>%
    select(Muestra, Edad, Tipo,
           CG = Metilacion_CG_mean,
           CA = Metilacion_CA_mean) %>%
    pivot_longer(
      cols = c(CG, CA),
      names_to = "Contexto",
      values_to = "Metilacion"
    ) %>%
    # Keep only complete observations
    filter(!is.na(Edad), !is.na(Metilacion)) %>%
    group_by(Tipo, Contexto) %>%
    # Fit one model per group and extract all statistics at once
    group_modify(~ {
      n <- nrow(.x)
      
      # Not enough points to fit a regression
      if (n < 3) {
        return(tibble(n = n, slope = NA_real_, slope_SE = NA_real_,
                      pvalue = NA_real_, intercept = NA_real_,
                      R_squared = NA_real_, adjusted_R_squared = NA_real_))
      }
      
      fit <- summary(lm(Metilacion ~ Edad, data = .x))
      coefs <- fit$coefficients
      
      tibble(
        n = n,
        slope = coefs["Edad", "Estimate"],
        slope_SE = coefs["Edad", "Std. Error"],
        pvalue = coefs["Edad", "Pr(>|t|)"],
        intercept = coefs["(Intercept)", "Estimate"],
        R_squared = fit$r.squared,
        adjusted_R_squared = fit$adj.r.squared
      )
    }) %>%
    ungroup()
}

# ========================================
# EM-seq
# ========================================
edades_emseq <- tibble(
  Muestra = c("M1", "M2", "M3", "M4", "M5", "M6", "M7"),
  Edad = c(35, 43, 55, 67, 78, 85, 103)
)

resultados_emseq <- age_methylation_regression(
  metilacion_WGBS_region_genomica_resultados_completos,
  edades_emseq
)

write.csv(resultados_emseq,
          "resultados_correlacion_edad_metilacion_CG_CA_EMSEQ.csv",
          row.names = FALSE)

# ========================================
# WGBS
# ========================================
edades_wgbs <- tibble(
  Muestra = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9"),
  Edad = c(81, 77, 77, 64, 49, 25, 16, 12, 5)
)

resultados_wgbs <- age_methylation_regression(
  metilacion_WGBS_region_genomica_resultados_completos,
  edades_wgbs
)

write.csv(resultados_wgbs,
          "resultados_correlacion_edad_metilacion_CG_CA_WGBS.csv",
          row.names = FALSE)