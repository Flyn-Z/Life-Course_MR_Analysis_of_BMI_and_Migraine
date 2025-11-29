
#prepare
library(ieugwasr)
get_opengwas_jwt()
library(TwoSampleMR)

#clump
ukb<-system.file("GWAS_statistics/BMI_UKB_p5e-8.csv",package="TwoSampleMR")
ukb_exposure<-read_exposure_data(filename=ukb,sep=",",snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele",other_allele_col="other_allele",eaf_col="AF_EUR",pval_col="pval",clump=FALSE)
bfile_path <- "./data_maf0.01_rs_ref/data_maf0.01_rs_ref"
ukb_clump <- ld_clump_local(
  dplyr::tibble(rsid=ukb_exposure$SNP, pval=ukb_exposure$pval.exposure),
  plink_bin = "./plink/plink.exe",
  bfile = bfile_path,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 1
)

moba8<-system.file("GWAS_statistics/BMI_moba8_p5e-6.csv",package="TwoSampleMR")
moba8_exposure<-read_exposure_data(filename=moba8,sep=",",snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele",eaf_col="EAF",other_allele_col="other_allele",pval_col="pval",clump=FALSE)
bfile_path <- "./data_maf0.01_rs_ref/data_maf0.01_rs_ref"
moba8_clump <- ld_clump_local(
  dplyr::tibble(rsid=moba8_exposure$SNP, pval=moba8_exposure$pval.exposure),
  plink_bin = "./plink/plink.exe",
  bfile = bfile_path,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 1
)

moba1<-system.file("GWAS_statistics/BMI_moba1_p5e-6.csv",package="TwoSampleMR")
moba1_exposure<-read_exposure_data(filename=moba1,sep=",",snp_col="SNP",beta_col="beta",se_col="se",effect_allele_col="effect_allele",eaf_col="EAF",other_allele_col="other_allele",pval_col="pval",clump=FALSE)
bfile_path <- "./data_maf0.01_rs_ref/data_maf0.01_rs_ref"
moba1_clump <- ld_clump_local(
  dplyr::tibble(rsid=moba1_exposure$SNP, pval=moba1_exposure$pval.exposure),
  plink_bin = "./plink/plink.exe",
  bfile = bfile_path,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 1
)

#outcome import
mig_total<-readr::read_csv("./outcome/migraine_total.csv")

#merge
merge_ukb_mig<-merge(ukb_clump,mig_total,by.x = "rsid",by.y = "SNP")
write.csv(merge_ukb_mig,"merge_ukb_mig.csv",row.names = FALSE)

merge_moba8_mig<-merge(moba8_clump,mig_total,by.x = "rsid",by.y = "SNP")
write.csv(merge_moba8_mig,"merge_moba8_mig.csv",row.names = FALSE)

merge_moba1_mig<-merge(moba1_clump,mig_total,by.x = "rsid",by.y = "SNP")
write.csv(merge_moba1_mig,"merge_moba1_mig.csv",row.names = FALSE)

#读取结局文件
outcome_data_ukb_mig<-read_outcome_data(snps = ukb_clump$rsid,filename = "merge_ukb_mig.csv",sep = ",",snp_col = "rsid",beta_col = "beta",se_col = "se",eaf_col="af_alt",effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "pval.y")

outcome_data_moba8_mig<-read_outcome_data(snps = moba8_clump$rsid,filename = "merge_moba8_mig.csv",sep = ",",snp_col = "rsid",beta_col = "beta",se_col = "se",eaf_col="af_alt",effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "pval.y")

outcome_data_moba1_mig<-read_outcome_data(snps = moba1_clump$rsid,filename = "merge_moba1_mig.csv",sep = ",",snp_col = "rsid",beta_col = "beta",se_col = "se",eaf_col="af_alt",effect_allele_col = "effect_allele",other_allele_col = "other_allele",pval_col = "pval.y")

#Harmonize
ukb_clump_formal<-merge(ukb_exposure,ukb_clump,by.x = "SNP",by.y = "rsid")
harm_ukb_mig<-harmonise_data(exposure_dat = ukb_clump_formal,outcome_dat = outcome_data_ukb_mig)
harm_ukb_mig<-harm_ukb_mig[harm_ukb_mig$pval.outcome >= 5e-8, ]
write.csv(harm_ukb_mig,"harm_ukb_mig.csv",row.names = FALSE)

moba8_clump_formal<-merge(moba8_exposure,moba8_clump,by.x = "SNP",by.y = "rsid")
harm_moba8_mig<-harmonise_data(exposure_dat = moba8_clump_formal,outcome_dat = outcome_data_moba8_mig)
harm_moba8_mig<-harm_moba8_mig[harm_moba8_mig$pval.outcome >= 5e-8, ]
write.csv(harm_moba8_mig,"harm_moba8_mig.csv",row.names = FALSE)

moba1_clump_formal<-merge(moba1_exposure,moba1_clump,by.x = "SNP",by.y = "rsid")
harm_moba1_mig<-harmonise_data(exposure_dat = moba1_clump_formal,outcome_dat = outcome_data_moba1_mig)
harm_moba1_mig<-harm_moba1_mig[harm_moba1_mig$pval.outcome >= 5e-8, ]
write.csv(harm_moba1_mig,"harm_moba1_mig.csv",row.names = FALSE)

#MR
x1<-generate_odds_ratios(mr_res = mr(harm_ukb_mig))
View(x1)
x2<-generate_odds_ratios(mr_res = mr(harm_moba8_mig))
View(x2)
x3<-generate_odds_ratios(mr_res = mr(harm_moba1_mig))
View(x3)

#heterogeneity
y1<-mr_heterogeneity(harm_ukb_mig)
View(y1)
run_mr_presso(harm_ukb_mig)
y2<-mr_heterogeneity(harm_moba8_mig)
View(y2)
run_mr_presso(harm_moba8_mig)
y3<-mr_heterogeneity(harm_moba1_mig)
View(y3)
run_mr_presso(harm_moba1_mig)

#horizontal pleiotropy
z1<-mr_pleiotropy_test(harm_ukb_mig)
View(z1)
z2<-mr_pleiotropy_test(harm_moba8_mig)
View(z2)
z3<-mr_pleiotropy_test(harm_moba1_mig)
View(z3)

#scattor plot
mr_scatter_plot(mr_results=mr(harm_ukb_mig,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median")),harm_ukb_mig)
mr_scatter_plot(mr_results=mr(harm_moba8_mig,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median")),harm_moba8_mig)
mr_scatter_plot(mr_results=mr(harm_moba1_mig,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median")),harm_moba1_mig)

#funnel plot
mr_funnel_plot(singlesnp_results=mr_singlesnp(harm_ukb_mig))
mr_funnel_plot(singlesnp_results=mr_singlesnp(harm_moba8_mig))
mr_funnel_plot(singlesnp_results=mr_singlesnp(harm_moba1_mig))

#leaveoneout
mr_leaveoneout_plot(leaveoneout_results=mr_leaveoneout(harm_ukb_mig))
mr_leaveoneout_plot(leaveoneout_results=mr_leaveoneout(harm_moba8_mig))
mr_leaveoneout_plot(leaveoneout_results=mr_leaveoneout(harm_moba1_mig))

#Steiger test
harm_ukb_mig$r.exposure <- harm_ukb_mig$beta.exposure / sqrt(harm_ukb_mig$beta.exposure^2 + harm_ukb_mig$se.exposure^2 * (harm_ukb_mig$samplesize.exposure - 2))
harm_ukb_mig$r.outcome <- get_r_from_lor(
  lor = harm_ukb_mig$beta.outcome,
  af = harm_ukb_mig$eaf.outcome,
  ncase = harm_ukb_mig$ncase.outcome,
  ncontrol = harm_ukb_mig$ncontrol.outcome,
  prevalence = harm_ukb_mig$prevalence.outcome
)
steiger_results_ukb_mig <- directionality_test(harm_ukb_mig)
print(steiger_results_ukb_mig)

harm_moba8_mig$r.exposure <- harm_moba8_mig$beta.exposure / sqrt(harm_moba8_mig$beta.exposure^2 + harm_moba8_mig$se.exposure^2 * (harm_moba8_mig$samplesize.exposure - 2))
harm_moba8_mig$r.outcome <- get_r_from_lor(
  lor = harm_moba8_mig$beta.outcome,
  af = harm_moba8_mig$eaf.outcome,
  ncase = harm_moba8_mig$ncase.outcome,
  ncontrol = harm_moba8_mig$ncontrol.outcome,
  prevalence = harm_moba8_mig$prevalence.outcome
)
steiger_results_moba8_mig <- directionality_test(harm_moba8_mig)
print(steiger_results_moba8_mig)

harm_moba1_mig$r.exposure <- harm_moba1_mig$beta.exposure / sqrt(harm_moba1_mig$beta.exposure^2 + harm_moba1_mig$se.exposure^2 * (harm_moba1_mig$samplesize.exposure - 2))
harm_moba1_mig$r.outcome <- get_r_from_lor(
  lor = harm_moba1_mig$beta.outcome,
  af = harm_moba1_mig$eaf.outcome,
  ncase = harm_moba1_mig$ncase.outcome,
  ncontrol = harm_moba1_mig$ncontrol.outcome,
  prevalence = harm_moba1_mig$prevalence.outcome
)
steiger_results_moba1_mig <- directionality_test(harm_moba1_mig)
print(steiger_results_moba1_mig)


