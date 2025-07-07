# ---- Load required packages and check installation ----
# Define required packages
required_pkgs <- c("tidyverse", "vroom", "TwoSampleMR", "patchwork", "gt", "meta")
# Check which packages are missing
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
# Stop if any required packages are missing
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}
# Load all required packages
invisible(lapply(required_pkgs, library, character.only = TRUE))

# Authenticate with OpenGWAS (required for some API calls)
ieugwasr::get_opengwas_jwt()

# ---- Read in conversion table for FinnGen outcomes ----
# Read the conversion table for FinnGen outcome names
conversions <- read_tsv("/Volumes/data/backup/data/CRP/finngen/conversions.csv")
conversions  # Print the conversion table

# ---- Define CRP region and load CRP GWAS summary statistics ----
# Define the CRP region (with 300kb padding on each side)
crp_start <- 159682079 - 3e5
crp_end <- 159684379 + 3e5
# Load CRP GWAS summary statistics
crp <- data.table::fread("/Volumes/data/backup/data/CRP/gwas/GCST90029070_buildGRCh37.tsv.gz") 

# ---- Select and clump CRP SNPs in the region ----
# Filter SNPs in the CRP region, select relevant columns, and perform LD clumping
crp_clumped <- crp %>%
  filter(
    chromosome == 1 & 
    (base_pair_location > crp_start) & 
    (base_pair_location < crp_end)
  ) %>%
  select(
    SNP = variant_id, 
    effect_allele.exposure = effect_allele, 
    other_allele.exposure = other_allele, 
    beta.exposure = beta, 
    se.exposure = standard_error, 
    pval.exposure = p_value
  ) %>%
  as_tibble() %>%
  mutate(eaf.exposure = NA) %>%  # EAF not available
  mutate(exposure = "CRP", id.exposure = "CRP") %>%
  mutate(pval.exposure = as.numeric(pval.exposure)) %>%
  filter(pval.exposure < 5e-8) %>%  # Only genome-wide significant SNPs
  mutate(pval = pval.exposure, rsid = SNP) %>%
  ieugwasr::ld_clump_local(
    bfile = "/Users/fh6520/1kg/EUR",
    plink_bin = genetics.binaRies::get_plink_binary(),
    clump_p = 1,
    clump_r2 = 0.01,
    clump_kb = 10000
  )

# View clumped CRP SNPs (for inspection)
crp_clumped %>%
  select(SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, pval.exposure) 

# ---- Meta-analysis function for fixed effects only ----
# Function to perform fixed-effect meta-analysis using the meta package
meta_fixed <- function(df) {
  m1 <- meta::metagen(df$b, df$se)
  # If all new_name are NA, set to NA_character_
  nn <- unique(df$new_name)
  if (length(nn) == 0 || all(is.na(nn))) nn <- NA_character_
  tibble(
    b = m1$TE.fixed,
    se = m1$seTE.fixed,
    pval = m1$pval.fixed,
    lower = m1$lower.fixed,
    upper = m1$upper.fixed,
    het_pval = m1$pval.Q
  )
}

# ---- Get OpenGWAS token and available outcomes ----
# Hard-coded OpenGWAS token (should be kept secure in practice)

# List available outcomes from OpenGWAS (Hamilton F only)
ao <- available_outcomes(opengwas_jwt = token) %>%
  filter(author == "Hamilton F")

# Extract UKB outcome data for clumped CRP SNPs
outcomes_ukb <- extract_outcome_data(crp_clumped$SNP, ao$id, opengwas_jwt = token)

# ---- List of FinnGen files to import ----
# List of FinnGen GWAS files for different outcomes
finngen_files <- c(
  "K11_APPENDACUT_crp_region.tsv.gz",
  "L12_CELLULITIS_crp_region.tsv.gz",
  "K11_CHOLECYST_crp_region.tsv.gz",
  "I9_ENDOCARD_crp_region.tsv.gz",
  "J10_ACUTELOWERNAS_crp_region.tsv.gz",
  "M13_OSTEOMYELITIS_crp_region.tsv.gz",
  "J10_PNEUMONIA_crp_region.tsv.gz",
  "J10_PNEUMOPNEUMO_crp_region.tsv.gz",
  "AB1_OTHER_SEPSIS_crp_region.tsv.gz",
  "AB1_STREPTO_SEPSIS_crp_region.tsv.gz",
  "J10_ACUTEUPPERINFEC_crp_region.tsv.gz",
  "N14_CYSTITIS_crp_region.tsv.gz"
)
# Directory containing FinnGen files
finngen_dir <- "/Volumes/data/backup/data/CRP/tabix/finngen_r12"

# ---- Function to read, TSMR-format, and filter a FinnGen file to clumped SNPs ----
# Reads a FinnGen file, formats it for TwoSampleMR, and filters to clumped SNPs
read_finngen_tsmr <- function(file) {
  path <- file.path(finngen_dir, file)
  df <- data.table::fread(cmd = paste("gzcat", shQuote(path)), data.table = FALSE)
  names(df) <- tolower(names(df))
  out <- tibble::tibble(
    SNP = if ("rsids" %in% names(df)) df$rsids else NA_character_,
    effect_allele.outcome = if ("alt" %in% names(df)) df$alt else NA_character_,
    other_allele.outcome = if ("ref" %in% names(df)) df$ref else NA_character_,
    beta.outcome = if ("beta" %in% names(df)) df$beta else NA_real_,
    se.outcome = if ("sebeta" %in% names(df)) df$sebeta else NA_real_,
    pval.outcome = if ("pval" %in% names(df)) df$pval else NA_real_,
    eaf.outcome = if ("af_alt" %in% names(df)) df$af_alt else NA_real_,
    outcome = gsub("_crp_region.tsv.gz$", "", file),
    id.outcome = file
  ) %>%
    filter(SNP %in% crp_clumped$SNP)
  out
}

# ---- Read and combine all FinnGen outcomes, TSMR-formatted and filtered to clumped SNPs ----
# Read all FinnGen files, format, and combine into one data frame
finngen_outcomes <- purrr::map_dfr(finngen_files, read_finngen_tsmr)

# ---- Read UKB Streptococcus pneumoniae pneumonia outcome ----
# Read UKB GWAS for Streptococcus pneumoniae pneumonia and format for MR
UKB_strep_pneumo <- vroom("/Volumes/data/backup/data/CRP/J13/01_J13.regenie.gz") %>%
  transmute(
    SNP = ID,
    beta.outcome = BETA,
    eaf.outcome = A1FREQ,
    other_allele.outcome = ALLELE0,
    effect_allele.outcome = ALLELE1,
    pval.outcome = 10^-LOG10P,
    se.outcome = SE,
    outcome = "Pneumonia due to Streptococcus pneumoniae",
    id.outcome = outcome
  ) %>%
  filter(SNP %in% crp_clumped$SNP)

# ---- Harmonise exposure and outcome data for MR ----
# Harmonise exposure and outcome data for MR analysis
data <- harmonise_data(
  exposure_dat = crp_clumped,
  outcome_dat = bind_rows(finngen_outcomes, UKB_strep_pneumo, outcomes_ukb), 
  action = 1
)

# ---- Assign harmonised disease names and cohort labels ----
# Assign harmonised disease names and cohort labels for each outcome
data_renamed <- data %>%
  mutate(
    new_name = case_when(
      id.outcome %in% c("ieu-b-4967", "K11_APPENDACUT_crp_region.tsv.gz") ~ "Appendicitis",
      id.outcome %in% c("ieu-b-4970", "L12_CELLULITIS_crp_region.tsv.gz") ~ "Cellulitis",
      id.outcome %in% c("ieu-b-4971", "K11_CHOLECYST_crp_region.tsv.gz") ~ "Cholecystitis",
      id.outcome %in% c("ieu-b-4972", "I9_ENDOCARD_crp_region.tsv.gz") ~ "Endocarditis",
      id.outcome %in% c("ieu-b-4973", "J10_ACUTELOWERNAS_crp_region.tsv.gz") ~ "LRTI",
      id.outcome %in% c("ieu-b-4975", "M13_OSTEOMYELITIS_crp_region.tsv.gz") ~ "Osteomyelitis",
      id.outcome %in% c("ieu-b-4976", "J10_PNEUMONIA_crp_region.tsv.gz") ~ "Pneumonia",
      id.outcome %in% c("J10_PNEUMOPNEUMO_crp_region.tsv.gz", "Pneumonia due to Streptococcus pneumoniae") ~ "Pneumonia (Streptococcus)",
      id.outcome %in% c("ieu-b-4980", "AB1_OTHER_SEPSIS_crp_region.tsv.gz") ~ "Sepsis",
      id.outcome == "AB1_STREPTO_SEPSIS_crp_region.tsv.gz" ~ "Sepsis (Streptococcus)",
      id.outcome %in% c("ieu-b-5063", "J10_ACUTEUPPERINFEC_crp_region.tsv.gz") ~ "URTI",
      id.outcome %in% c("ieu-b-5065", "N14_CYSTITIS_crp_region.tsv.gz") ~ "UTI",
      TRUE ~ id.outcome
    ),
    cohort = case_when(
      str_detect(id.outcome, "_crp_region.tsv.gz$") ~ "FinnGen",
      TRUE ~ "UKB"
    )
  )
data_renamed  # Print the harmonised and labelled data

# ---- Create a lookup table for harmonised names and cohorts ----
# Create a lookup table for new_name, cohort, and id.outcome
link <- data_renamed %>% select(new_name, cohort, id.outcome) %>% distinct()

# ---- Specify MR methods to use ----
# List of MR methods to use in the analysis
mr_methods_to_use <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median")

# ---- Run MR for each study and join harmonised names ----
# Run MR for each outcome and method, join harmonised names, and filter out missing p-values
all_mr_results <- TwoSampleMR::mr(data_renamed, method_list = mr_methods_to_use) %>%
  left_join(link, by = c("id.outcome")) %>%
  filter(!is.na(pval))

# ---- Identify outcomes with both FinnGen and UKB for meta-analysis ----
# Identify outcomes (by new_name and method) with results from both FinnGen and UKB
cohort_res <- all_mr_results %>%
  group_by(new_name, method) %>%
  add_count() %>%
  filter(n == 2) %>%
  ungroup()

# ---- Perform meta-analysis for outcomes with both cohorts ----
# For outcomes with both FinnGen and UKB, perform fixed-effect meta-analysis and combine with cohort results
final_for_plot <- all_mr_results %>%
  group_by(new_name, method) %>%
  add_count() %>%
  filter(n == 2) %>%
  nest() %>%
  mutate(results = map(data, meta_fixed)) %>%
  select(results) %>%
  unnest(results) %>%
  mutate(cohort = "Meta-analysis") %>%
  select(new_name, method, b, se, pval, het_pval, cohort) %>%
  bind_rows(cohort_res) %>%
  select(new_name, method, b, se, pval, het_pval, cohort) 

# ---- Plot forest plot for meta-analysed results ----
all_mr_results %>%
  group_by(new_name, method) %>%
  add_count() %>%
  filter(n == 2) %>%
  nest() %>%
  mutate(results = map(data, meta_fixed)) %>%
  select(results) %>%
  unnest(results) %>%
  mutate(cohort = "Meta-analysis") %>%
  select(new_name, method, b, se, pval, het_pval, cohort) %>%
  filter(method == "Inverse variance weighted") %>%
  arrange(-b) %>%
  ungroup() %>%
  mutate(new_name = if_else(new_name == "Pneumonia (Streptococcus)", "Streptooccal pneumonia", new_name)) %>%
  ggforestplot::forestplot(name = new_name, estimate = b, se = se) +
  xlab("Odds ratio per change in natural log CRP (95% CI)")

# Save the forest plot to PDF
ggsave("R1/Figure1.pdf", height = 8, width = 8)

# ---- Arrange results by p-value ----
# Arrange the final results by p-value (for inspection)
final_for_plot %>% arrange(pval)

# ---- Format results for table output ----
# Format results for table: odds ratio (95% CI), p-value, and output as a gt table
final_for_plot %>%
  mutate(
    column_res = paste0(
      round(exp(b), 2), 
      " (", round(exp(b - 1.96 * se), 2), " ,", round(exp(b + 1.96 * se), 2), 
      "), p = ", pval = signif(pval, 3)
    )
  ) %>%
  select(-b, -se, -pval, -het_pval) %>%
  pivot_wider(values_from = column_res, names_from = method) %>%
  gt::gt()

# ---- Final harmonised data for summary output ----
# This object contains all harmonised SNP-outcome-exposure data, with harmonised names and cohort labels.
final_harmonised_data <- data_renamed

# Write the harmonised data to CSV for record-keeping or further analysis
write_tsv(final_harmonised_data, "R1/harmonised_mr_data.csv")
