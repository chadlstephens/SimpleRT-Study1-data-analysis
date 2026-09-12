# ==============================================================================
# Study 1: Model 1.1 - Cardiac Phase Effect on RT
# Purpose: Fit the primary mixed-effects model from Section 2.4.3.1.1.
#
# CardiacTim coding CONFIRMED by Dr. Yang (email): 1 = systole, 2 = diastole.
# CardiacPhase is effect-coded systole = -0.5, diastole = +0.5, so a NEGATIVE
# beta1 means RT is faster (lower) in diastole than systole, and a POSITIVE
# beta1 means RT is faster in systole than diastole. Direction is now safe to
# interpret and report.
#
# Run import_clean_study1.R FIRST -- this script assumes `analysis_df` already
# exists in your environment with columns: Subject, RT_clean, CardiacPhase,
# CardiacPhase_raw, and reflects the corrected N=49 sample (subjects 32, 45,
# 48 excluded; 39 included) once you've updated raw_dir accordingly.
# ==============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(MuMIn)

# ---- Sanity check before fitting anything ----
if (!exists("analysis_df")) {
  stop("analysis_df not found. Run import_clean_study1.R first.")
}

cat("Fitting Model 1.1 on N =", n_distinct(analysis_df$Subject), "subjects.\n")
if (n_distinct(analysis_df$Subject) != 49) {
  cat("NOTE: Expected final N is 49 (Yang et al., 2017; subjects 32, 45, 48 excluded, 39 included).",
      "Current N is", n_distinct(analysis_df$Subject),
      "-- check that raw_dir has ace045.txt/ace048.txt removed and ace039.txt added,",
      "then re-run import_clean_study1.R before treating results as final.\n\n")
} else {
  cat("N = 49 confirmed -- matches Yang et al. (2017)'s reported final sample.\n\n")
}


# ==============================================================================
# Step 1: Fit random-intercept-only and random-slope models
# ==============================================================================

m_intercept <- lmer(
  RT_clean ~ CardiacPhase + (1 | Subject),
  data = analysis_df,
  REML = TRUE
)

m_slope <- lmer(
  RT_clean ~ CardiacPhase + (1 + CardiacPhase | Subject),
  data = analysis_df,
  REML = TRUE
)

# ---- Likelihood ratio test to decide which random-effects structure to report ----
# anova() on lmer objects automatically refits with ML for the comparison
lrt_result <- anova(m_intercept, m_slope)
cat("---- Likelihood ratio test: random intercept vs. random slope ----\n")
print(lrt_result)

# Pick the better-fitting model based on the LRT (alpha = .05)
lrt_p <- lrt_result$`Pr(>Chisq)`[2]
use_slope_model <- !is.na(lrt_p) && lrt_p < .05

final_model <- if (use_slope_model) m_slope else m_intercept
cat("\nSelected model:", if (use_slope_model) "random slope (m_slope)" else "random intercept only (m_intercept)", "\n\n")


# ==============================================================================
# Step 2: Report the primary model output
# ==============================================================================

cat("\n================ MODEL 1.1 OUTPUT ================\n")
summary(final_model)

# ---- Fixed effect: CardiacPhase (beta1) ----
fixed_ests <- summary(final_model)$coefficients
beta1      <- fixed_ests["CardiacPhase", "Estimate"]
se1        <- fixed_ests["CardiacPhase", "Std. Error"]
df1        <- fixed_ests["CardiacPhase", "df"]
t1         <- fixed_ests["CardiacPhase", "t value"]
p1         <- fixed_ests["CardiacPhase", "Pr(>|t|)"]

ci <- confint(final_model, parm = "CardiacPhase", method = "Wald")

cat("\n---- beta1 (CardiacPhase) ----\n")
cat("Estimate:  ", round(beta1, 2), "ms\n")
cat("SE:        ", round(se1, 2), "\n")
cat("95% CI:    [", round(ci[1], 2), ",", round(ci[2], 2), "]\n")
cat("t(", round(df1, 1), ") = ", round(t1, 2), ", p = ", format.pval(p1, digits = 3), "\n", sep = "")

# Direction is now interpretable: CardiacPhase is coded systole = -0.5, diastole = +0.5
faster_phase <- if (beta1 < 0) "diastole" else if (beta1 > 0) "systole" else "neither (beta1 = 0)"
cat("\nDirection: beta1", if (beta1 < 0) "< 0" else "> 0",
    "-> RT is numerically faster in", faster_phase, "\n")
cat("(Significance per the p-value above determines whether this difference is reliable.)\n")

# ---- Effect size (approx. Cohen's d from the effect-coded predictor) ----
# Since CardiacPhase spans 1 unit (-0.5 to +0.5), beta1 IS the raw mean difference.
resid_sd <- sigma(final_model)
cohens_d <- beta1 / resid_sd
cat("\nApprox. Cohen's d (beta1 / residual SD):", round(cohens_d, 3), "\n")

# ---- Variance components ----
cat("\n---- Variance components ----\n")
print(VarCorr(final_model))

# ---- R-squared (marginal and conditional) ----
r2_vals <- r.squaredGLMM(final_model)
cat("\n---- R-squared ----\n")
cat("Marginal R2 (fixed effects only):   ", round(r2_vals[1, "R2m"], 4), "\n")
cat("Conditional R2 (fixed + random):    ", round(r2_vals[1, "R2c"], 4), "\n")


# ==============================================================================
# Step 3: Save the model object, tagged with the actual N it was fit on
# ==============================================================================
n_tag <- n_distinct(analysis_df$Subject)
out_file <- paste0("model_1_1_N", n_tag, ".rds")
saveRDS(list(model = final_model, N = n_tag, lrt = lrt_result), out_file)

cat("\n================ END MODEL 1.1 OUTPUT ================\n")
cat("\nSaved to", out_file, "\n")
if (n_tag != 49) {
  cat("Re-run this script once raw_dir reflects the corrected N=49 sample",
      "(subjects 32/45/48 excluded, 39 included) to get the file result you'll report.\n")
} else {
  cat("N = 49 confirmed. This is the sample you'll report.\n")
}