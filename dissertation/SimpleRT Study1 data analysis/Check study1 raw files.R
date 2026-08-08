# ==============================================================================
# Study 1 Raw Data Validation Script
# Purpose: Confirm structural consistency across all 50 raw .txt files BEFORE
#          stacking them with bind_rows()/map_dfr(). Run this first — if either
#          check below fails, do not proceed to the import/cleaning pipeline
#          until the flagged file(s) are resolved.
#          NOTE: Original design N was 52; 2 subjects were excluded due to
#          equipment failure, so 50 files/subjects is the expected count here.
#          NOTE: These E-Prime exports are UTF-16LE encoded (confirmed via
#          guess_encoding() + a zero-problems result from problems()) — all
#          read_tsv() calls below specify locale(encoding = "UTF-16LE") and
#          quote = "" to read them correctly.
# ==============================================================================

library(tidyverse)

# ---- 0. Point this at your raw data folder -----------------------------------
raw_dir <- "/Users/stephens/R/dissertation/data/study1_raw"   # <-- update to your actual path
expected_n_files <- 50          # 52 originally enrolled, minus 2 excluded for equipment failure

# Subject IDs excluded for equipment failure — fill these in with the actual IDs.
# This lets Check 2 confirm the RIGHT two subjects are missing, not just that
# the count happens to be 50 (which could mask a different subject being
# dropped by mistake, e.g. a duplicate file overwrite or a misplaced file).
excluded_subject_ids <- c()   # <-- e.g. c("14", "29")

file_list <- list.files(raw_dir, pattern = "\\.txt$", full.names = TRUE)

cat("Found", length(file_list), "files in", raw_dir, "\n")
if (length(file_list) != expected_n_files) {
  warning("Expected ", expected_n_files, " files, found ", length(file_list),
          ". Check the folder path and file pattern before continuing.")
}

# ==============================================================================
# CHECK 1: Do all 52 files have identical column names, in identical order?
# ==============================================================================

get_header <- function(path) {
  # n_max = 0 reads only the header row (fast; avoids parsing all data)
  # quote = "" disables quote-character interpretation, since embedded XML in
  # Clock.Information (e.g. dt:dt="string") would otherwise confuse the parser.
  # locale(encoding = "UTF-16LE") is required — these E-Prime exports are
  # UTF-16LE encoded, not UTF-8 (confirmed via guess_encoding() + zero-problems
  # result from problems()).
  names(suppressWarnings(read_tsv(
    path, skip = 1, n_max = 0, quote = "",
    locale = locale(encoding = "UTF-16LE"),
    show_col_types = FALSE
  )))
}

headers <- map(file_list, get_header)
names(headers) <- basename(file_list)

# Use the first file's header as the reference to compare all others against
reference_file <- basename(file_list[1])
reference_header <- headers[[1]]

header_check <- tibble(
  file          = names(headers),
  n_columns     = map_int(headers, length),
  same_n_cols   = n_columns == length(reference_header),
  same_names_set = map_lgl(headers, ~ setequal(.x, reference_header)),
  same_order    = map_lgl(headers, ~ identical(.x, reference_header))
)

cat("\n---- CHECK 1: Column name/order consistency ----\n")
cat("Reference file:", reference_file, "with", length(reference_header), "columns\n\n")

if (all(header_check$same_order)) {
  cat("PASS: All", nrow(header_check), "files have identical column names in identical order.\n")
} else {
  cat("FAIL: Column mismatch detected in the following file(s):\n\n")
  print(header_check %>% filter(!same_order))
  
  # For any file with the same columns but different order, show the reordering
  reorder_only <- header_check %>% filter(same_names_set & !same_order)
  if (nrow(reorder_only) > 0) {
    cat("\nFiles with the SAME column set but DIFFERENT order (safe to fix via column selection,\n",
        "but bind_rows() by name would still work correctly — flagged for awareness):\n")
    print(reorder_only$file)
  }
  
  # For any file with genuinely different columns (extra/missing), show the diff
  true_mismatch <- header_check %>% filter(!same_names_set)
  if (nrow(true_mismatch) > 0) {
    cat("\nFiles with GENUINELY DIFFERENT columns (extra and/or missing relative to reference)\n",
        "— these require manual inspection before merging:\n")
    for (f in true_mismatch$file) {
      extra   <- setdiff(headers[[f]], reference_header)
      missing <- setdiff(reference_header, headers[[f]])
      cat("\n  ", f, "\n")
      if (length(extra) > 0)   cat("    Extra columns:   ", paste(extra, collapse = ", "), "\n")
      if (length(missing) > 0) cat("    Missing columns: ", paste(missing, collapse = ", "), "\n")
    }
  }
}

# ==============================================================================
# CHECK 2: 52 unique Subject values, and 530 rows per file
# ==============================================================================

get_subject_and_nrow <- function(path) {
  df <- suppressWarnings(read_tsv(
    path, skip = 1, quote = "",
    locale = locale(encoding = "UTF-16LE"),
    show_col_types = FALSE
  ))
  tibble(
    file        = basename(path),
    n_rows      = nrow(df),
    n_subjects_in_file = n_distinct(df$Subject, na.rm = TRUE),
    subject_id  = paste(unique(df$Subject), collapse = "; ")  # flags files with >1 Subject value
  )
}

structure_check <- map_dfr(file_list, get_subject_and_nrow)

cat("\n---- CHECK 2: Row counts and Subject ID consistency ----\n")

# 2a. Row count check: expect 530 rows per file (10 practice + 520 main trials)
row_count_fail <- structure_check %>% filter(n_rows != 530)
if (nrow(row_count_fail) == 0) {
  cat("PASS: All", nrow(structure_check), "files have exactly 530 rows.\n")
} else {
  cat("FAIL: The following file(s) do NOT have 530 rows:\n\n")
  print(row_count_fail %>% select(file, n_rows))
}

# 2b. Each file should contain exactly ONE unique Subject value
multi_subject_fail <- structure_check %>% filter(n_subjects_in_file != 1)
if (nrow(multi_subject_fail) == 0) {
  cat("\nPASS: Every file contains exactly one unique Subject value.\n")
} else {
  cat("\nFAIL: The following file(s) contain more than one Subject value (or none):\n\n")
  print(multi_subject_fail %>% select(file, n_subjects_in_file, subject_id))
}

# 2c. Across all 52 files, are there exactly 52 UNIQUE Subject values?
#     (i.e., no file is a duplicate/re-export of another subject's data)
all_subject_ids <- structure_check$subject_id
n_unique_subjects <- n_distinct(all_subject_ids)

cat("\nTotal files:", nrow(structure_check),
    "| Unique Subject IDs across files:", n_unique_subjects, "\n")

if (n_unique_subjects != expected_n_files) {
  cat("NOTE: Expected", expected_n_files, "unique subjects (52 enrolled minus 2 excluded",
      "for equipment failure) but found", n_unique_subjects, "- confirm this matches your",
      "current exclusion records.\n")
}

if (n_unique_subjects == nrow(structure_check)) {
  cat("PASS: No duplicate Subject IDs across files.\n")
} else {
  dup_ids <- structure_check %>%
    count(subject_id, name = "n_files") %>%
    filter(n_files > 1)
  cat("FAIL: The following Subject ID(s) appear in more than one file:\n\n")
  print(dup_ids)
  cat("\nFiles sharing a duplicated Subject ID:\n")
  print(structure_check %>% filter(subject_id %in% dup_ids$subject_id) %>%
          select(file, subject_id))
}

# 2d. Confirm the SPECIFIC subjects excluded for equipment failure are the
#     ones actually missing — guards against the count being 50 for the
#     wrong reason (e.g. a different subject's file was accidentally omitted
#     or duplicated, while one of the intended-exclusion IDs slipped through).
if (length(excluded_subject_ids) == 0) {
  cat("\nNOTE: 'excluded_subject_ids' is empty — fill in the two equipment-failure",
      "Subject IDs at the top of this script to enable this check.\n")
} else {
  present_ids <- unique(structure_check$subject_id)
  
  wrongly_present <- intersect(excluded_subject_ids, present_ids)
  if (length(wrongly_present) > 0) {
    cat("\nFAIL: The following excluded Subject ID(s) were found in the raw data folder",
        "(they should have been removed for equipment failure):\n")
    print(wrongly_present)
  } else {
    cat("\nPASS: None of the excluded Subject ID(s) (", paste(excluded_subject_ids, collapse = ", "),
        ") are present among the retained files.\n")
  }
  
  # Also confirm every OTHER expected subject (i.e. not one of the two exclusions)
  # is present, in case a non-excluded subject's file went missing instead.
  # This assumes subject IDs are sequential 1:52 with the two exclusions removed;
  # adjust 'all_expected_ids' if your original numbering scheme differs.
  all_expected_ids <- setdiff(as.character(1:52), excluded_subject_ids)
  unexpectedly_missing <- setdiff(all_expected_ids, present_ids)
  if (length(unexpectedly_missing) > 0) {
    cat("\nFAIL: The following Subject ID(s) were expected (not in the exclusion list)",
        "but are MISSING from the raw data folder:\n")
    print(unexpectedly_missing)
  } else if (length(excluded_subject_ids) > 0) {
    cat("PASS: All non-excluded expected Subject IDs (1-52 minus the 2 exclusions) are present.\n")
  }
}

# ==============================================================================
# Summary table (for a quick full-picture scan)
# ==============================================================================
cat("\n---- Full summary table ----\n")
summary_tbl <- structure_check %>%
  left_join(header_check %>% select(file, same_n_cols, same_order), by = "file")
print(summary_tbl, n = 60)