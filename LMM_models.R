###############################################################################
# Linear mixed-effects model (LMM) for microbiome outcomes
# Random effects: participant nested within city
# Model per feature: y ~ treatment * time + (1|city) + (1|city:participant)
###############################################################################

# ---- Packages ----
# install.packages(c("lme4", "rsq", "readr", "dplyr", "tibble"))
library(lme4)
library(rsq)
library(readr)
library(dplyr)
library(tibble)

# ---- 0) CONFIG 
meta_file <- "Sample_info.txt"        # metadata table
feat_file <- "Microbial_data.txt"     # feature table (columns = features), same row order as metadata OR merge by sample_id

meta_sep <- ""                        # "" = whitespace via read.table; use "," or "\t" if needed
feat_sep <- ""                        # same

# If you need to merge by sample ID, set this to the shared column name in BOTH tables; else NULL assumes row order matches.
sample_id_col <- NULL                 # e.g., "sample_id"

# Required columns in metadata
participant_col <- "ID"
city_col        <- "City"
treat_col       <- "treatment"        # e.g., Control / Intervention
time_col        <- "Time"             # e.g., N1 / N2


# Reference levels (first = reference)
treat_levels <- c("Control", "Intervention")
time_levels  <- c("N1", "N2")

# Output prefix
out_prefix <- "LMM_results_nested_city"

###############################################################################
# 1) READ DATA
###############################################################################
meta <- if (meta_sep == "") read.table(meta_file, header = TRUE, check.names = FALSE) else
  read.delim(meta_file, sep = meta_sep, header = TRUE, check.names = FALSE)

feat <- if (feat_sep == "") read.table(feat_file, header = TRUE, check.names = FALSE) else
  read.delim(feat_file, sep = feat_sep, header = TRUE, check.names = FALSE)


###############################################################################
# 2) ALIGN SAMPLES
###############################################################################
# If sample_id_col is provided, merge by it; otherwise assume row order matches.
if (!is.null(sample_id_col)) {
  stopifnot(sample_id_col %in% colnames(meta), sample_id_col %in% colnames(feat))
  dat <- meta %>% inner_join(feat, by = sample_id_col)
} else {
  if (nrow(meta) != nrow(feat)) stop("meta and feature tables have different number of rows; provide sample_id_col to merge.")
  dat <- bind_cols(meta, feat)
}

###############################################################################
# 3) BASIC CHECKS + FACTOR CODING
###############################################################################
needed <- c(participant_col, city_col, treat_col, time_col)
missing <- setdiff(needed, colnames(dat))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

dat[[participant_col]] <- factor(dat[[participant_col]])
dat[[city_col]]        <- factor(dat[[city_col]])
dat[[treat_col]]       <- factor(as.character(dat[[treat_col]]), levels = treat_levels)
dat[[time_col]]        <- factor(as.character(dat[[time_col]]),  levels = time_levels)

if (anyNA(dat[[treat_col]])) stop("Treatment has values not in treat_levels. Check unique(dat[[treat_col]]).")
if (anyNA(dat[[time_col]]))  stop("Time has values not in time_levels. Check unique(dat[[time_col]]).")

# Interaction grouping factor for nested random effects
dat$city_participant <- interaction(dat[[city_col]], dat[[participant_col]], drop = TRUE)

###############################################################################
# 4) IDENTIFY FEATURE COLUMNS (numeric outcomes only)
###############################################################################
feature_cols <- names(dat)[sapply(dat, is.numeric)]
# If you have numeric metadata columns, exclude them explicitly here:
feature_cols <- setdiff(feature_cols, c())  # e.g., c("age", "BMI")

if (length(feature_cols) == 0) stop("No numeric feature columns found to model.")

###############################################################################
# 5) FIT LMM PER FEATURE + COLLECT RESULTS
###############################################################################
results <- vector("list", length(feature_cols))
names(results) <- feature_cols

for (f in feature_cols) {
  y <- dat[[f]]
  if (anyNA(y)) next 
  
  # Group summaries (means/SD by treatment x time)
  summ <- dat %>%
    group_by(.data[[treat_col]], .data[[time_col]]) %>%
    summarise(mean = mean(.data[[f]], na.rm = TRUE),
              sd   = sd(.data[[f]], na.rm = TRUE),
              .groups = "drop")
  
  # LMM: nested random effects
  form <- as.formula(
    paste0(f, " ~ ", treat_col, " * ", time_col,
           " + (1|", city_col, ") + (1|city_participant)")
  )
  
  m <- lmer(form, data = dat)
  
  # Fixed effects 
  cf <- as.data.frame(coef(summary(m)))
  cf$p_norm <- 2 * (1 - pnorm(abs(cf$`t value`)))
  
  # R2 decomposition
  r2 <- rsq(m, adj = FALSE, type = "v")
  
  # Random intercept variances
  vc <- as.data.frame(VarCorr(m))
  var_city <- vc$vcov[vc$grp == city_col][1]
  var_city_participant <- vc$vcov[vc$grp == "city_participant"][1]
  
  # Helper to extract coefficient rows robustly
  get_row <- function(term) {
    if (!term %in% rownames(cf)) return(c(NA, NA, NA, NA))
    c(cf[term, "Estimate"], cf[term, "Std. Error"], cf[term, "t value"], cf[term, "p_norm"])
  }
  
  # Find terms automatically (works even if factor coding changes)
  treat_term <- grep(paste0("^", treat_col), rownames(cf), value = TRUE)[1]
  time_term  <- grep(paste0("^", time_col),  rownames(cf), value = TRUE)[1]
  int_term   <- grep(":", rownames(cf), value = TRUE)[1]
  
  tr <- get_row(treat_term)
  tm <- get_row(time_term)
  it <- get_row(int_term)
  
  out <- tibble(
    feature = f,
    
    # Means/SDs (kept, since your original output included them)
    mean_control_T0 = summ$mean[summ[[treat_col]] == treat_levels[1] & summ[[time_col]] == time_levels[1]],
    sd_control_T0   = summ$sd  [summ[[treat_col]] == treat_levels[1] & summ[[time_col]] == time_levels[1]],
    mean_control_T1 = summ$mean[summ[[treat_col]] == treat_levels[1] & summ[[time_col]] == time_levels[2]],
    sd_control_T1   = summ$sd  [summ[[treat_col]] == treat_levels[1] & summ[[time_col]] == time_levels[2]],
    mean_case_T0    = summ$mean[summ[[treat_col]] == treat_levels[2] & summ[[time_col]] == time_levels[1]],
    sd_case_T0      = summ$sd  [summ[[treat_col]] == treat_levels[2] & summ[[time_col]] == time_levels[1]],
    mean_case_T1    = summ$mean[summ[[treat_col]] == treat_levels[2] & summ[[time_col]] == time_levels[2]],
    sd_case_T1      = summ$sd  [summ[[treat_col]] == treat_levels[2] & summ[[time_col]] == time_levels[2]],
    
    # Fixed effects
    intercept_est = cf["(Intercept)", "Estimate"],
    intercept_se  = cf["(Intercept)", "Std. Error"],
    intercept_t   = cf["(Intercept)", "t value"],
    intercept_p   = cf["(Intercept)", "p_norm"],
    
    treatment_est = tr[1], treatment_se = tr[2], treatment_t = tr[3], treatment_p = tr[4],
    time_est      = tm[1], time_se      = tm[2], time_t      = tm[3], time_p      = tm[4],
    interaction_est = it[1], interaction_se = it[2], interaction_t = it[3], interaction_p = it[4],
    
    # Model fit
    model_r2  = r2$model,
    fixed_r2  = r2$fixed,
    random_r2 = r2$random,
    var_city = var_city,
    var_city_participant = var_city_participant,
    AIC = AIC(m)
  )
  
  results[[f]] <- out
}

results_df <- bind_rows(results)

###############################################################################
# 7) FDR ADJUSTMENT ACROSS FEATURES (BH)
###############################################################################
# Adjust separately for each focal term across all modeled features
q_treat <- p.adjust(p_treat, method = "BH")
q_time  <- p.adjust(p_time,  method = "BH")
q_int   <- p.adjust(p_int,   method = "BH")

# Map q-values back to results table
results_df <- results_df %>%
  mutate(
    treatment_q    = q_treat[feature],
    time_q         = q_time[feature],
    interaction_q  = q_int[feature]
  )

###############################################################################
# 7) WRITE OUTPUT
###############################################################################
write.csv(results_df, paste0(out_prefix, ".csv"), row.names = FALSE)
message("Done. Wrote: ", paste0(out_prefix, ".csv"))


