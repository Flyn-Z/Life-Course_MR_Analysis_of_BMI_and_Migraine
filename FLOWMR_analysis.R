library(FLOWMR)
library(dplyr)
library(GRAPPLE)

sel.file <- c("./exposure/BMI_giant17eu.csv", "./exposure/BMIchild_egg15.csv" )
exp.file <- c("./exposure/BMI1year_moba19.csv","./exposure/BMI8year_moba19.csv","./exposure/BMIadult_ukb.csv" )
out.file <- "./outcome/migraine_total.csv"
plink_refdat <- "./data_maf0.01_rs_ref/data_maf0.01_rs_ref"

data.list <- GRAPPLE::getInput(sel.file, exp.file, out.file, plink_refdat, max.p.thres =0.01, plink_exe = "./plink/plink.exe", cal.cor = T)
View(data.list[["data"]])
View(data.list[["cor.mat"]])

dat <-data.list$data
Gamma_hat =rbind(dat$gamma_out1,dat$gamma_exp3,dat$gamma_exp2,dat$gamma_exp1)
Sd_hat = rbind(dat$se_out1,dat$se_exp3,dat$se_exp2,dat$se_exp1)
cor_mat = data.list$cor.mat
result = BayesMediation(Gamma_hat, Sd_hat, cor = cor_mat, inv = TRUE)
View(result[["summary"]])
write.csv(result[["summary"]],"result_BMI_k=4_migraine_total.csv")

path_effect_1 = indirect_effect(result$raw, K = 3, path = c(3,2,1))
print(path_effect)
path_effect_2 = indirect_effect(result$raw, K = 3, path = c(4,3,1))
print(path_effect)
path_effect_3 = indirect_effect(result$raw, K = 3, path = c(4,2,1))
print(path_effect)

ind_effect = indirect_effect(result$raw, K = 4, path = "all")
print(ind_effect)

tot_effect = total_effect(result$raw, K = 4)
print(tot_effect)

FLOWMR::traceplot(result$raw, par = "B", ind = c(1,2))
FLOWMR::traceplot(result$raw, par = "B", ind = c(1,3))
FLOWMR::traceplot(result$raw, par = "B", ind = c(1,4))


