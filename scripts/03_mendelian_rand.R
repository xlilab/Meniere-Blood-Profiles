
rm(list = ls())

library(dplyr)
library(TwoSampleMR)
library(ggplot2)

set.seed(0)

# SNPs -> CE (exposure) -> Meniere (outcome) --------------------------------------------------------------
exp_datunclump <- read_exposure_data(
  filename = "GCST90025818_buildGRCh37.tsv.gz",
  clump = FALSE,
  sep= "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col ="effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value")  %>%
  filter(
    pval.exposure < 1e-06 ) %>%
    #pval.exposure < 5e-08 ) %>%
  mutate(
    exposure = "Total free cholesterol levels")
colnames(exp_datunclump)[colnames(exp_datunclump)=="SNP"] = "rsid"
colnames(exp_datunclump)[colnames(exp_datunclump)=="pval.exposure"] = "pval"

exp_dat <- ieugwasr::ld_clump_local(exp_datunclump, clump_kb = 1000,
                                    clump_r2 = 0.001,
                                    clump_p = 1,
                                    bfile = "ref/MAGMA/reference/g1000_eur/g1000_eur",
                                    plink_bin = "plink")

colnames(exp_dat)[colnames(exp_dat)=="pval"] = "pval.exposure"
colnames(exp_dat)[colnames(exp_dat)=="rsid"] = "SNP"

outcome_dat <- read_outcome_data(
  snps = exp_dat$SNP,
  filename = "finngen_R11_H8_MENIERE.gz",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval") %>%
  mutate(
  outcome = "Meniere Disease")

missing <- outcome_dat %>% filter(!SNP %in% exp_dat$SNP) %>% pull(SNP)

exp_out_dat <- harmonise_data(exp_dat, outcome_dat)

mr_pleiotropy_test(exp_out_dat)

# run mr-presso to perform correction of horizontal pleiotropy via outlier removal
exp_out_presso <- run_mr_presso(
  exp_out_dat,
  NbDistribution = 2000,
  SignifThreshold = 0.05)

# main mr-presso results
exp_out_presso[[1]]$`Main MR results`
# global test results
exp_out_presso[[1]]$`MR-PRESSO results`$`Global Test`

# HP outliers
exp_out_dat %>%
  filter(
    row_number() %in% exp_out_presso[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  ) #%>% select(1:9)

# remove HP outliers
exp_out_dat_adj <- exp_out_dat %>%
  filter(
    ! row_number() %in% exp_out_presso[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  )

# re-test for horizontal pleiotropy
mr_pleiotropy_test(exp_out_dat_adj)

# MR
exp_out_dat_adj = exp_out_dat

# perform MR
exp_out_res <- mr(exp_out_dat_adj)

# odds ratio
exp_out_res <- generate_odds_ratios(exp_out_res)

# plot effects
scatterplot = mr_scatter_plot(exp_out_res, exp_out_dat_adj )[[1]] + theme_bw()

# leave one out analysis
exp_out_loo <- mr_leaveoneout(exp_out_dat_adj)

save(exp_out_loo,exp_out_res,exp_out_dat_adj,file = "GCST90025818R11H8MENIERE_1e6_MAGMArefclump.RData")

# SNPs -> Meniere (exposure) -> CE (outcome) ---------------------------------
exp_datunclump <- read_exposure_data(
  filename = "finngen_R11_H8_MENIERE.gz",
  sep = "\t",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval") %>%
  filter(
    pval.exposure < 1e-06
    #pval.exposure < 5e-08
  ) %>%
  mutate( exposure = "Meniere Disease")

colnames(exp_datunclump)[colnames(exp_datunclump)=="SNP"] = "rsid"
colnames(exp_datunclump)[colnames(exp_datunclump)=="pval.exposure"] = "pval"

exp_dat <- ieugwasr::ld_clump_local(exp_datunclump, clump_kb = 1000,
                                    clump_r2 = 0.001,
                                    clump_p = 1,
                                    bfile = "ref/MAGMA/reference/g1000_eur/g1000_eur",
                                    plink_bin = "plink")


colnames(exp_dat)[colnames(exp_dat)=="pval"] = "pval.exposure"
colnames(exp_dat)[colnames(exp_dat)=="rsid"] = "SNP"

# get effects of instruments on outcome
outcome_dat <- read_outcome_data(
  snps = exp_dat$SNP,
  filename = "GCST90025818_buildGRCh37.tsv.gz",
  sep= "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col ="effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
) %>%
  mutate(
    outcome = "Total free cholesterol levels"
  )

missing <- outcome_dat %>% filter(!SNP %in% exp_dat$SNP) %>% pull(SNP)

exp_out_dat <- harmonise_data(exp_dat,outcome_dat)

mr_pleiotropy_test(exp_out_dat)

exp_out_presso <- run_mr_presso(
  exp_out_dat,
  NbDistribution = 2000,
  SignifThreshold = 0.05
)

# main mr-presso results
exp_out_presso[[1]]$`Main MR results`

# global test results
exp_out_presso[[1]]$`MR-PRESSO results`$`Global Test`

# HP outliers
exp_out_dat %>%
  filter(
    row_number() %in% exp_out_presso[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)

exp_out_presso[[1]]$`MR-PRESSO results`$`Distortion Test`

# remove HP outliers
exp_out_dat_adj <- exp_out_dat %>%
  filter(
    ! row_number() %in% exp_out_presso[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  )

# re-test for horizontal pleiotropy
mr_pleiotropy_test(exp_out_dat_adj)

# MR
exp_out_dat_adj = exp_out_dat
exp_out_res <- mr( exp_out_dat_adj)

# odds ratio
exp_out_res <- generate_odds_ratios(exp_out_res)

# plot effects
mr_scatter_plot(exp_out_res, exp_out_dat_adj)[[1]]

# leave one out analysis
exp_out_loo <- mr_leaveoneout(exp_out_dat_adj)

save(exp_out_loo,exp_out_res,exp_out_dat_adj,file = "R11H8MENIEREGCST90025818_1e6_MAGMArefclump.RData")
