# ==============================================================================
# Study 1: Import, Clean, and Diagnose
# Purpose: Stack all 50 raw .txt files, apply the no-face-trial cleaning +
#          Ratcliff (1993) per-subject-per-phase outlier trimming, and inspect
#          diagnostics BEFORE trusting the resulting analysis dataset.
#
# Run check_study1_raw_files.R FIRST and confirm it passes before running this.
# ==============================================================================

library(tidyverse)

# ---- 0. Point this at your raw data folder ------------------------------------
raw_dir <- "/Users/stephens/R/dissertation/data/study1_raw"   # <-- update to your actual path
expected_n_files <- 50          # 52 enrolled, minus 2 excluded for equipment failure

# Sample accounting per Yang et al. (2017): 52 enrolled -> 2 excluded (equipment
# failure) -> 1 withdrew -> 49 analyzed. The 2 equipment-failure exclusions were
# already applied before you received these 50 files. The withdrawn subject's
# file IS still present among the 50 -- exclude it here once the ID is confirmed.
withdrawn_subject_id <- NA_character_   # <-- fill in once Dr. Yang confirms, e.g. "17"


# ==============================================================================
# STEP 3: Import step — stack all 50 files
# ==============================================================================

file_list <- list.files(raw_dir, pattern = "\\.txt$", full.names = TRUE)

cat("Found", length(file_list), "files in", raw_dir, "\n")
if (length(file_list) != expected_n_files) {
  warning("Expected ", expected_n_files, " files, found ", length(file_list),
          ". Re-run check_study1_raw_files.R before proceeding.")
}

read_one <- function(path) {
  # quote = "" disables quote-character interpretation, since embedded XML in
  # Clock.Information (e.g. dt:dt="string") would otherwise confuse read_tsv's
  # quoting logic and misalign column counts for nearly every row.
  # locale(encoding = "UTF-16LE") is required because these E-Prime exports are
  # UTF-16LE encoded, not UTF-8 — confirmed via readr::guess_encoding() and a
  # zero-problems result from problems() after applying this fix.
  read_tsv(
    path,
    skip = 1,
    quote = "",
    locale = locale(encoding = "UTF-16LE"),
    show_col_types = FALSE
  ) %>%
    mutate(source_file = basename(path))
}

raw_all <- map_dfr(file_list, read_one)

cat("\n---- Import summary ----\n")
cat("Total rows imported: ", nrow(raw_all), "\n")
cat("Unique Subjects:     ", n_distinct(raw_all$Subject), "\n")
cat("Unique source files: ", n_distinct(raw_all$source_file), "\n")

if (n_distinct(raw_all$Subject) != expected_n_files) {
  warning("Unique Subject count (", n_distinct(raw_all$Subject),
          ") does not match expected_n_files (", expected_n_files, ").")
}

# ---- Apply the withdrawn-subject exclusion, if known ----
if (is.na(withdrawn_subject_id)) {
  cat("\nNOTE: 'withdrawn_subject_id' is not yet set. All", n_distinct(raw_all$Subject),
      "subjects will be carried through cleaning for now. Do NOT treat any results",
      "as final until this is set and the pipeline is re-run -- the analyzed N",
      "should be 49 per Yang et al. (2017), not", n_distinct(raw_all$Subject), ".\n")
} else {
  if (!withdrawn_subject_id %in% as.character(raw_all$Subject)) {
    warning("withdrawn_subject_id ('", withdrawn_subject_id,
            "') was not found among the imported Subject IDs -- double-check it before proceeding.")
  }
  n_before <- n_distinct(raw_all$Subject)
  raw_all <- raw_all %>% filter(as.character(Subject) != withdrawn_subject_id)
  cat("\nApplied withdrawal exclusion: removed Subject", withdrawn_subject_id,
      "->", n_before, "subjects reduced to", n_distinct(raw_all$Subject), "\n")
}


# ==============================================================================
# STEP 4: Cleaning + diagnostics functions
# ==============================================================================

# ---- 4a. Cleaning function: tags every trial with its disposition ----
clean_study1 <- function(df) {
  df %>%
    filter(`Running[Block]` == "BloackList") %>%                     # drop practice trials
    filter(`Procedure[Trial]` != "TrialProc11") %>%                   # drop catch trials
    filter(`Procedure[Trial]` %in% c("TrialProc1", "TrialProc2")) %>% # no-face trials only
    mutate(
      RT_clean  = coalesce(!!!select(., matches("^Proc\\d{2}End\\.RT$"))),
      ACC_clean = coalesce(!!!select(., matches("^Proc\\d{2}End\\.ACC$"))),
      Miss = ACC_clean == 1 | RT_clean == 0,
      CardiacPhase_raw = CardiacTim,        # VERIFY 1/2 = systole/diastole before proceeding
      CardiacPhase = case_when(
        CardiacPhase_raw == 1 ~ -0.5,       # tentative: systole
        CardiacPhase_raw == 2 ~  0.5,       # tentative: diastole
        TRUE ~ NA_real_
      )
    ) %>%
    group_by(Subject, CardiacPhase_raw) %>%
    mutate(
      cell_n_prefilter = sum(!Miss),                          # non-miss trials available for trimming
      subj_mean = if (cell_n_prefilter[1] >= 10) mean(RT_clean[!Miss], na.rm = TRUE) else NA_real_,
      subj_sd   = if (cell_n_prefilter[1] >= 10) sd(RT_clean[!Miss], na.rm = TRUE)   else NA_real_,
      is_outlier = !Miss & !is.na(subj_sd) & abs(RT_clean - subj_mean) > 3 * subj_sd,
      low_n_flag = cell_n_prefilter < 10,                      # flag degenerate cells instead of silently dropping
      Disposition = case_when(
        low_n_flag           ~ "flagged_low_n",                # kept in output, NOT auto-trimmed
        Miss                 ~ "dropped_miss",
        is_outlier           ~ "dropped_outlier",
        TRUE                 ~ "retained"
      )
    ) %>%
    ungroup()
}

# ---- 4b. Diagnostics: trials in vs. retained per subject x phase ----
trim_diagnostics <- function(tagged_df) {
  tagged_df %>%
    group_by(Subject, CardiacPhase_raw) %>%
    summarise(
      n_total       = n(),
      n_miss        = sum(Disposition == "dropped_miss"),
      n_outlier     = sum(Disposition == "dropped_outlier"),
      n_retained    = sum(Disposition == "retained"),
      low_n_flagged = any(Disposition == "flagged_low_n"),
      pct_retained  = round(100 * n_retained / n_total, 1),
      .groups = "drop"
    ) %>%
    arrange(Subject, CardiacPhase_raw)
}

# ---- 4c. Final analysis-ready dataset (retained trials only) ----
get_clean_trials <- function(tagged_df) {
  tagged_df %>%
    filter(Disposition == "retained") %>%
    select(-cell_n_prefilter, -subj_mean, -subj_sd, -is_outlier, -low_n_flag, -Disposition)
}

# ---- Run the pipeline ----
tagged      <- clean_study1(raw_all)
diag        <- trim_diagnostics(tagged)
analysis_df <- get_clean_trials(tagged)


# ==============================================================================
# STEP 5: Inspect diagnostics BEFORE trusting the analysis dataset
# ==============================================================================

cat("\n================ DIAGNOSTIC REPORT ================\n")

# 5a. Any cells with too few pre-trim trials to compute a stable SD?
low_n_cells <- diag %>% filter(low_n_flagged)
cat("\n-- 5a. Cells flagged for low N (<10 non-miss trials pre-trim) --\n")
if (nrow(low_n_cells) == 0) {
  cat("None. Every subject x phase cell had >=10 trials available for trimming.\n")
} else {
  cat(nrow(low_n_cells), "cell(s) flagged — inspect individually before deciding",
      "whether to exclude the subject/cell or trim anyway with a caveat:\n\n")
  print(low_n_cells)
}

# 5b. Overall retention rate
cat("\n-- 5b. Overall trial retention --\n")
diag %>%
  summarise(
    mean_pct_retained = round(mean(pct_retained), 1),
    min_pct_retained  = min(pct_retained),
    max_pct_retained  = max(pct_retained)
  ) %>%
  print()

# 5c. Subjects/cells with unusually aggressive trimming (<80% retained)
cat("\n-- 5c. Cells with <80% retained (unusually aggressive trimming/missingness) --\n")
aggressive_trim <- diag %>% filter(pct_retained < 80)
if (nrow(aggressive_trim) == 0) {
  cat("None. All cells retained >=80% of trials.\n")
} else {
  cat(nrow(aggressive_trim), "cell(s) below 80% retention — worth checking whether this",
      "reflects genuine noise or a data-quality problem for that subject:\n\n")
  print(aggressive_trim)
}

# 5d. Does the average trials/phase/subject land near the expected ~40?
cat("\n-- 5d. Trials per subject x phase vs. expected ~40 --\n")
diag %>%
  group_by(CardiacPhase_raw) %>%
  summarise(
    mean_retained = round(mean(n_retained), 1),
    min_retained  = min(n_retained),
    max_retained  = max(n_retained),
    .groups = "drop"
  ) %>%
  print()

# 5e. Final N check: how many subjects actually make it into the analysis dataset,
#     and does each contribute both cardiac-phase cells?
cat("\n-- 5e. Subject coverage in the final analysis dataset --\n")
subject_phase_coverage <- analysis_df %>%
  distinct(Subject, CardiacPhase_raw) %>%
  count(Subject, name = "n_phases_present")

incomplete_subjects <- subject_phase_coverage %>% filter(n_phases_present < 2)
cat("Subjects in analysis_df: ", n_distinct(analysis_df$Subject), " (expected ", expected_n_files, ")\n", sep = "")
if (nrow(incomplete_subjects) == 0) {
  cat("PASS: Every subject contributes trials to BOTH cardiac-phase cells.\n")
} else {
  cat("FLAG: The following subject(s) are missing one cardiac-phase cell entirely",
      "after cleaning (e.g. all trials in that cell were trimmed as outliers,",
      "flagged low-N, or missed):\n\n")
  print(incomplete_subjects)
}

cat("\n-- 5f. Sample size accounting vs. Yang et al. (2017) --\n")
n_enrolled   <- 52
n_equipment  <- 2
n_withdrawal <- if (is.na(withdrawn_subject_id)) NA_integer_ else 1
n_current    <- n_distinct(analysis_df$Subject)
n_expected_final <- 49

cat("Enrolled:                     ", n_enrolled, "\n")
cat("Excluded (equipment failure): ", n_equipment, " (subjects 32, 39)\n", sep = "")
cat("Excluded (withdrawal):        ",
    if (is.na(withdrawn_subject_id)) "PENDING -- not yet applied" else paste0("1 (subject ", withdrawn_subject_id, ")"),
    "\n", sep = "")
cat("Subjects in analysis_df now:  ", n_current, "\n")
cat("Expected final N (Yang et al., 2017): ", n_expected_final, "\n")

if (is.na(withdrawn_subject_id)) {
  cat("STATUS: Withdrawal not yet applied -- current N (", n_current,
      ") will NOT match the expected final N (", n_expected_final,
      ") until withdrawn_subject_id is set. Treat all current results as provisional.\n", sep = "")
} else if (n_current == n_expected_final) {
  cat("PASS: Current N matches Yang et al. (2017)'s reported final sample of 49.\n")
} else {
  cat("FLAG: Current N (", n_current, ") does NOT match the expected 49 even after applying",
      " the withdrawal exclusion -- investigate before proceeding.\n", sep = "")
}

cat("\n================ END DIAGNOSTIC REPORT ================\n")
cat("\nReview all sections above. Do not proceed to model fitting until:\n",
    "  - 5a shows no unresolved low-N cells (or you've made a documented decision on them)\n",
    "  - 5c aggressive-trimming cells have been checked for data-quality issues\n",
    "  - 5e shows no subjects missing a full cardiac-phase cell\n",
    "  - You have confirmed the systole/diastole coding for CardiacPhase_raw (see clean_study1())\n",
    "  - withdrawn_subject_id is set and the pipeline re-run, so the final N matches\n",
    "    Yang et al. (2017)'s reported 49\n", sep = "")

# analysis_df is now ready for the Model 1.1 lmer() step, once the above are resolved.