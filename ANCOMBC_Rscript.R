############################################################
# ANCOM-BC2
############################################################

library(ANCOMBC)

########################################
# 1) INPUT DATA
#######################################
#Read microbial data
df <- read.delim("microbiomedata.txt", header=TRUE, sep="\t", check.names=FALSE)
otu_mat <- as.matrix(df); storage.mode(otu_mat) <- "numeric"

#Read metadata
meta_df <- read.csv2("metadata.csv", sep=";", check.names=FALSE, stringsAsFactors=FALSE)
names(meta_df) <- trimws(names(meta_df))
names(meta_df) <- sub("^\ufeff", "", names(meta_df))
rownames(meta_df) <- meta_df$sample_id


############################################################
# 1) Align samples + Quality check deadwood
############################################################

deadwood_raw <- "N_deadwood_decayed"

# Ensure sample IDs match between meta and OTU matrix
stopifnot(all(colnames(otu_mat) %in% rownames(meta_df)))

# Reorder meta_df to match otu_mat columns
meta_df <- meta_df[colnames(otu_mat), , drop = FALSE]

# Deadwood numeric
meta_df[[deadwood_raw]] <- as.numeric(as.character(meta_df[[deadwood_raw]]))


############################################################
# 2) Filter rare taxa 
#    Keep taxa present in >=10% samples AND total count >=50
############################################################

prev <- rowMeans(otu_mat > 0)
tot  <- rowSums(otu_mat)

keep_taxa <- (prev >= 0.10) & (tot >= 50)
otu_mat_filt <- otu_mat[keep_taxa, , drop = FALSE]

cat("Taxa kept:", nrow(otu_mat_filt), "of", nrow(otu_mat), "\n")
cat("Samples:", ncol(otu_mat_filt), "\n")

############################################################
# 3) Create deadwood_log1p
############################################################

meta_df$deadwood_log1p <- log1p(meta_df[[deadwood_raw]])
stopifnot(!any(is.na(meta_df$deadwood_log1p)))

############################################################
# 4) Formula + ANCOM-BC2
############################################################

fix_formula <- "deadwood_log1p * N * treatment + City"
rand_formula <- "(1 | ID)"

res_long <- ancombc2(
  data         = otu_mat_filt,
  meta_data    = meta_df,
  fix_formula  = fix_formula,
  rand_formula = rand_formula,
  p_adj_method = "BH",
  prv_cut      = 0,     # already filtered above
  lib_cut      = 0,
  s0_perc      = 0.05,
  struc_zero   = FALSE,
  neg_lb       = TRUE,
  alpha        = 0.05,
  global       = FALSE
)
