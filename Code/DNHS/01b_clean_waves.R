# DNAm collection age was matched to the corresponding study-wave age to identify the focal wave. 
# Trauma and stress were calculated cumulatively through the focal wave, using baseline totals for each cohort and adding events reported at subsequent waves; when DNAm age matched multiple waves, cumulative values were averaged across matching waves. 
# Smoking was coded as ever smoked (1 if smoking was reported at or before the focal wave; 0 otherwise), with prior smoking status carried forward to Wave 5.

smoke01 <- function(x) case_when(
  x == 1 ~ 1,
  .default = 0
)

ever_smoke <- function(...) {
  as.numeric(any(c(...) == 1, na.rm = TRUE))
}

df_focal <- my_df %>%
  mutate(
    # Cumulative trauma
    trauma_w1 = w1c1_traumanum,
    trauma_w2 = case_when(
      !is.na(w1c1_traumanum) ~ w1c1_traumanum + coalesce(w2c1_traumanum, 0),
      !is.na(w2c2_traumanum) ~ w2c2_traumanum,
      TRUE ~ NA_real_
    ),
    trauma_w3 = trauma_w2 + coalesce(w3_traumanum, 0),
    trauma_w4 = trauma_w3 + coalesce(w4_traumanum, 0),
    trauma_w5 = trauma_w4 + coalesce(w5_traumanum, 0),
    
    # Cumulative stress
    stress_w1 = w1c1_stressnum,
    stress_w2 = case_when(
      !is.na(w1c1_stressnum) ~ w1c1_stressnum + coalesce(w2c1_stressnum, 0),
      !is.na(w2c2_stressnum) ~ w2c2_stressnum,
      TRUE ~ NA_real_
    ),
    stress_w3 = stress_w2 + coalesce(w3_stressnum, 0),
    stress_w4 = stress_w3 + coalesce(w4_stressnum, 0),
    stress_w5 = stress_w4 + coalesce(w5_stressnum, 0)
  ) %>%
  rowwise() %>%
  mutate(
    # Ever smoked by each wave
    smoke_w1 = smoke01(w1c1_CP1A),
    smoke_w2 = case_when(
      !is.na(w1c1_CP1A) ~ ever_smoke(smoke_w1, smoke01(w2c1_CP1A)),
      TRUE ~ smoke01(w2c2_CP1A)
    ),
    smoke_w3 = ever_smoke(smoke_w2, smoke01(w3_CP1A)),
    smoke_w4 = ever_smoke(smoke_w3, smoke01(w4_CP1A)),
    smoke_w5 = smoke_w4,
    
    # Match DNAm age to wave
    n_age_matches = sum(c_across(w1_age_final:w5_age_final) == Age, na.rm = TRUE),
    
    focal_wave = {
      idx <- which(c_across(w1_age_final:w5_age_final) == Age)
      if (length(idx)) paste(idx, collapse = "/") else NA_character_
    },
    
    focal_age = if_else(n_age_matches > 0, Age, NA_real_),
    
    # Mean cumulative trauma across matching waves
    focal_trauma = {
      idx <- which(c_across(w1_age_final:w5_age_final) == Age)
      values <- c_across(trauma_w1:trauma_w5)[idx]
      
      if (!length(values) || all(is.na(values))) NA_real_
      else mean(values, na.rm = TRUE)
    },
    
    # Mean cumulative stress across matching waves
    focal_stress = {
      idx <- which(c_across(w1_age_final:w5_age_final) == Age)
      values <- c_across(stress_w1:stress_w5)[idx]
      
      if (!length(values) || all(is.na(values))) NA_real_
      else mean(values, na.rm = TRUE)
    },
    
    # Ever smoked by focal wave
    focal_smoking = {
      idx <- which(c_across(w1_age_final:w5_age_final) == Age)
      values <- c_across(smoke_w1:smoke_w5)[idx]
      
      if (!length(values)) 0
      else max(values)
    }
  ) %>%
  ungroup()