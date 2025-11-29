#prepare
library(ieugwasr)
get_opengwas_jwt()
library(TwoSampleMR)
library(MendelianRandomization)
library(dplyr)
library(data.table)
library(MVMR)

#import
ukb<-readr::read_csv("./exposure/BMIadult_ukb.csv")
moba8<-readr::read_csv("./exposure/BMI8year_moba19.csv")
moba1<-readr::read_csv("./exposure/BMI1year_moba19.csv")

#SNP selection
ukb_filtered <- ukb %>% filter(pval < 5e-08)
moba8_filtered <- moba8 %>% filter(pval < 5e-06)
moba1_filtered <- moba1 %>% filter(pval < 5e-06)

ukb_all_snps <- unique(ukb$SNP)
moba8_all_snps <- unique(moba8$SNP)
moba1_all_snps <- unique(moba1$SNP)

ukb_filtered_1 <- ukb_filtered %>% 
  filter(SNP %in% moba8_all_snps & SNP %in% moba1_all_snps)
moba8_filtered_1 <- moba8_filtered %>% 
  filter(SNP %in% ukb_all_snps & SNP %in% moba1_all_snps)
moba1_filtered_1 <- moba1_filtered %>% 
  filter(SNP %in% ukb_all_snps & SNP %in% moba8_all_snps)

#SNP merge
final_common_snps <- unique(c(ukb_filtered_1$SNP, moba8_filtered_1$SNP, moba1_filtered_1$SNP))

#prepare clump
clump_data <- data.frame(
  rsid = final_common_snps,
  pval = sapply(final_common_snps, function(snp) {
    p_vals <- c(
      ifelse(snp %in% ukb_filtered_1$SNP, ukb_filtered_1$pval[ukb_filtered_1$SNP == snp], NA),
      ifelse(snp %in% moba8_filtered_1$SNP, moba8_filtered_1$pval[moba8_filtered_1$SNP == snp], NA),
      ifelse(snp %in% moba1_filtered_1$SNP, moba1_filtered_1$pval[moba1_filtered_1$SNP == snp], NA)
    )
    min(p_vals, na.rm = TRUE)
  })
)
any(is.na(clump_data))

#clump
bfile_path <- "./data_maf0.01_rs_ref/data_maf0.01_rs_ref"
clumped_snps <- ld_clump_local(
  clump_data,
  plink_bin = "./plink/plink.exe",
  bfile = bfile_path,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p = 1
)
write.csv(clumped_snps,"clumped_snps.csv",row.names = FALSE)

#extract SNP parameter
ukb_formal<-read_exposure_data(filename="./exposure/BMIadult_ukb.csv",sep=",",snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele",other_allele_col="other_allele",eaf_col="AF_EUR",pval_col="pval",clump=FALSE)
moba8_formal<-read_exposure_data(filename="./exposure/BMI8year_moba19.csv",sep=",",snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele",eaf_col="EAF",other_allele_col="other_allele",pval_col="pval",clump=FALSE)
moba1_formal<-read_exposure_data(filename="./exposure/BMI1year_moba19.csv",sep=",",snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele",eaf_col="EAF",other_allele_col="other_allele",pval_col="pval",clump=FALSE)

mvmr_exposure_data <- data.frame()

for(snp in clumped_snps$rsid) {
  ukb_snp <- ukb_formal[ukb_formal$SNP == snp, ]
  if(nrow(ukb_snp) > 0) {
    ukb_row <- ukb_snp[1, ]
    ukb_row$id.exposure <- "adult_BMI"
    ukb_row$exposure <- "Adult BMI"
    mvmr_exposure_data <- rbind(mvmr_exposure_data, ukb_row)
  }
  
  moba8_snp <- moba8_formal[moba8_formal$SNP == snp, ]
  if(nrow(moba8_snp) > 0) {
    moba8_row <- moba8_snp[1, ]
    moba8_row$id.exposure <- "child8_BMI"
    moba8_row$exposure <- "Child 8-year BMI"
    mvmr_exposure_data <- rbind(mvmr_exposure_data, moba8_row)
  }
  
  moba1_snp <- moba1_formal[moba1_formal$SNP == snp, ]
  if(nrow(moba1_snp) > 0) {
    moba1_row <- moba1_snp[1, ]
    moba1_row$id.exposure <- "child1_BMI"
    moba1_row$exposure <- "Child 1-year BMI"
    mvmr_exposure_data <- rbind(mvmr_exposure_data, moba1_row)
  }
}

any(is.na(mvmr_exposure_data))

write.csv(mvmr_exposure_data, "mvmr_exposure_data_formatted.csv", row.names = FALSE)

#outcome import
mig_total<-readr::read_csv("./outcome/migraine_total.csv")

#merge
merge<-merge(clumped_snps,mig_total,by.x = "rsid",by.y = "SNP")
write.csv(merge,"merge_exp_migtotal.csv",row.names = FALSE)

#read_outcome_data
outcome_data<-read_outcome_data(snps = clumped_snps$rsid,filename = "merge_exp_migtotal.csv",sep = ",",snp_col = "rsid",beta_col = "beta",se_col = "se",eaf_col="af_alt",effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "pval.y")

#harmonise
data<-mv_harmonise_data(exposure_dat = mvmr_exposure_data, outcome_dat = outcome_data)

#Use MendelianRandomization package
exposure_betas <- data$exposure_beta
exposure_ses <- data$exposure_se
outcome_beta <- data$outcome_beta
outcome_se <- data$outcome_se

#format conversion
mvmr_input <- MendelianRandomization::mr_mvinput(bx = exposure_betas,
                                                 bxse = exposure_ses,
                                                 by = outcome_beta,
                                                 byse = outcome_se,
                                                 exposure = c("Adult_BMI", "Child1_BMI", "Child8_BMI"),
                                                 outcome = "Migraine_total")

#IVW
mvmr_result <- MendelianRandomization::mr_mvivw(mvmr_input, model = "default")

results_df <- data.frame(
  Exposure = c("Adult BMI", "Child 1-year BMI", "Child 8-year BMI"),
  Estimate = mvmr_result@Estimate,
  StdError = mvmr_result@StdError,
  CI_lower = mvmr_result@CILower,
  CI_upper = mvmr_result@CIUpper,
  P_value = mvmr_result@Pvalue
)
View(results_df)

#heterogeneity
q_stat <- mvmr_result@Heter.Stat
print(paste("Q_statistics:", q_stat[1], "P_value:", q_stat[2]))

#horizontal pleiotropy
mvmr_egger <- MendelianRandomization::mr_mvegger(mvmr_input)
print(mvmr_egger)

#Use MVMR package
exposure_betas <- data$exposure_beta
exposure_ses <- data$exposure_se
outcome_beta <- data$outcome_beta
outcome_se <- data$outcome_se

#format conversion
mvmr_data <- format_mvmr(BXGs = exposure_betas,
                         BYG = outcome_beta,
                         seBXGs = exposure_ses,
                         seBYG = outcome_se,
                         RSID = rownames(exposure_betas))

#Conditional F-statistics
mvmr_strength_result <- strength_mvmr(r_input = mvmr_data)

