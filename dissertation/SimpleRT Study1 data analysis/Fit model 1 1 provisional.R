# ==============================================================================
# Study 1: Model 1.1 - Cardiac Phase Effect on RT (PROVISIONAL, N=50)
# Purpose: Fit the primary mixed-effects model from Section 2.4.3.1.1 on the
#          current N=50 sample, as a provisional check while awaiting two
#          confirmations from Dr. Yang:
#            1. Which CardiacTim value (1 or 2) is systole vs. diastole
#            2. Which subject ID withdrew (final N should be 49, not 50)
#
# Run import_clean_study1.R FIRST -- this script assumes `analysis_df` already
# exists in your environment with columns: Subject, RT_clean, CardiacPhase,
# CardiacPhase_raw.
#
# *** DO NOT report or interpret the DIRECTION of the beta1 coefficient from
# *** this run. The sign depends entirely on the still-unconfirmed systole/
# *** diastole mapping. Magnitude, significance, and model fit statistics are
# *** unaffected and safe to inspect now.
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
  cat("NOTE: This is a PROVISIONAL fit. Expected final N is 49 (Yang et al., 2017);",
      "current N is", n_distinct(analysis_df$Subject),
      "because the withdrawn-subject exclusion has not yet been applied.\n")
}
cat("NOTE: CardiacPhase direction (systole = -0.5 / diastole = +0.5) is UNCONFIRMED.",
    "Do not interpret the sign of beta1 until Dr. Yang confirms the CardiacTim coding.\n\n")


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

cat("\n================ MODEL 1.1 OUTPUT (PROVISIONAL) ================\n")
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
cat("\n*** Magnitude/significance above are interpretable now. DO NOT state which\n",
    "    phase (systole/diastole) is faster until CardiacTim coding is confirmed. ***\n", sep = "")

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
# Step 3: Save model objects for later comparison against the N=49 fit
# ==============================================================================
saveRDS(list(model = final_model, N = n_distinct(analysis_df$Subject), lrt = lrt_result),
        "model_1_1_provisional_N50.rds")

cat("\n================ END MODEL 1.1 OUTPUT ================\n")
cat("\nSaved to model_1_1_provisional_N50.rds for later comparison against the N=49 fit,\n",
    "once withdrawn_subject_id is confirmed and applied in import_clean_study1.R.\n", sep = "")