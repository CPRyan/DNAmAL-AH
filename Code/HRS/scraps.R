
m1 <- lm(
  sle_score1_scaled ~ Female + PAGE + smoker1213 + DNAmAL_r_scaled,
  data = hrs_sle1
)

m2 <- update(
  m1,
  . ~ . - DNAmAL_scaled + methAL_scaled
)

m3 <- update(
  m1,
  . ~ . - DNAmAL_scaled + methALT_scaled
)


# ============================================================
# SLE score 2
# ============================================================

hrs_sle2 <- hrs %>%
  drop_na(
    sle_score2_scaled,
    all_of(covars),
    all_of(dnam_measures)
  )

m4 <- lm(
  sle_score2_scaled ~ Female + PAGE + smoker1213 + DNAmAL_scaled,
  data = hrs_sle2
)

m5 <- update(
  m4,
  . ~ . - DNAmAL_scaled + methAL_scaled
)

m6 <- update(
  m4,
  . ~ . - DNAmAL_scaled + methALT_scaled
)


# ============================================================
# SLE score 3
# ============================================================

hrs_sle3 <- hrs %>%
  drop_na(
    sle_score3_scaled,
    all_of(covars),
    all_of(dnam_measures)
  )

m7 <- lm(
  sle_score3_scaled ~ Female + PAGE + smoker1213 + DNAmAL_scaled,
  data = hrs_sle3
)

m8 <- update(
  m7,
  . ~ . - DNAmAL_scaled + methAL_scaled
)

m9 <- update(
  m7,
  . ~ . - DNAmAL_scaled + methALT_scaled
)



########
# Save
########

# ============================================================
# Save model tables
# ============================================================

# Function
source("Code/HRS/f_write_sjplot_xlsx.R")



pred_labels <- c(
  "(Intercept)",
  "Sex (F)",
  "Age",
  "Smoker (Former)",
  "Smoker (Current)",
  "DNAmAL (scaled)",
  "methAL (scaled)",
  "methALT (scaled)"
)

sle_model_sets <- list(
  "SLE score 1" = list(m1, m2, m3),
  "SLE score 2" = list(m4, m5, m6),
  "SLE score 3" = list(m7, m8, m9)
)

write_sjplot_xlsx(
  model_sets = sle_model_sets,
  output_dir = "Output/HRS/Tables",
  xlsx_name = "SLE_DNAm_model_tables.xlsx",
  dv_label = "SLE Score (scaled)",
  pred_labels = pred_labels
)
