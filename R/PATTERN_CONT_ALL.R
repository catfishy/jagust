library(lme4)
library(coefplot2)
library(ggplot2)
library(lmerTest)
library(pbkrtest)
library(multcomp)
library(contrast)
library(xtable)
library(sjPlot)

# Import data for Ns
df_ef = read.csv('dpgmm_alpha15.37_bilateral_UW_EF__ALL_longdata_continuous.csv')
df_mem = read.csv('dpgmm_alpha15.37_bilateral_UW_MEM__ALL_longdata_continuous.csv')
df_adas = read.csv('dpgmm_alpha15.37_bilateral_ADAScog__ALL_longdata_continuous.csv')
df_avlt = read.csv('dpgmm_alpha15.37_bilateral_AVLT__ALL_longdata_continuous.csv')
df_wmh = read.csv('dpgmm_alpha15.37_bilateral_WMH_percentOfICV__ALL_longdata_continuous.csv')
df_av45 = read.csv('dpgmm_alpha15.37_bilateral_AV45_ALL_longdata_continuous.csv')
df_hcicv = read.csv('dpgmm_alpha15.37_bilateral_FSX_HC_ICV__ALL_longdata_continuous.csv')
df_ef$CORTICAL_SUMMARY_POSITIVE = factor(df_ef$CORTICAL_SUMMARY_POSITIVE)
df_mem$CORTICAL_SUMMARY_POSITIVE = factor(df_mem$CORTICAL_SUMMARY_POSITIVE)
df_adas$CORTICAL_SUMMARY_POSITIVE = factor(df_adas$CORTICAL_SUMMARY_POSITIVE)
df_avlt$CORTICAL_SUMMARY_POSITIVE = factor(df_avlt$CORTICAL_SUMMARY_POSITIVE)
df_wmh$CORTICAL_SUMMARY_POSITIVE = factor(df_wmh$CORTICAL_SUMMARY_POSITIVE)
df_av45$CORTICAL_SUMMARY_POSITIVE = factor(df_av45$CORTICAL_SUMMARY_POSITIVE)
df_hcicv$CORTICAL_SUMMARY_POSITIVE = factor(df_hcicv$CORTICAL_SUMMARY_POSITIVE)
df_ef$RID = factor(df_ef$RID)
df_mem$RID = factor(df_mem$RID)
df_adas$RID = factor(df_adas$RID)
df_avlt$RID = factor(df_avlt$RID)
df_wmh$RID = factor(df_wmh$RID)
df_av45$RID = factor(df_av45$RID)
df_hcicv$RID = factor(df_hcicv$RID)
df_ef$group = factor(df_ef$group)
df_mem$group = factor(df_mem$group)
df_adas$group = factor(df_adas$group)
df_avlt$group = factor(df_avlt$group)
df_av45$group = factor(df_av45$group)
df_wmh$group = factor(df_wmh$group)
df_hcicv$group = factor(df_hcicv$group)
df_ef$diag_prior = factor(df_ef$diag_prior)
df_mem$diag_prior = factor(df_mem$diag_prior)
df_adas$diag_prior = factor(df_adas$diag_prior)
df_avlt$diag_prior = factor(df_avlt$diag_prior)
df_wmh$diag_prior = factor(df_wmh$diag_prior)
df_av45$diag_prior = factor(df_av45$diag_prior)
df_hcicv$diag_prior = factor(df_hcicv$diag_prior)
df_ef$APOE4_BIN = factor(df_ef$APOE4_BIN)
df_mem$APOE4_BIN = factor(df_mem$APOE4_BIN)
df_adas$APOE4_BIN = factor(df_adas$APOE4_BIN)
df_avlt$APOE4_BIN = factor(df_avlt$APOE4_BIN)
df_wmh$APOE4_BIN = factor(df_wmh$APOE4_BIN)
df_av45$APOE4_BIN = factor(df_av45$APOE4_BIN)
df_hcicv$APOE4_BIN = factor(df_hcicv$APOE4_BIN)
df_ef$Gender = factor(df_ef$Gender)
df_mem$Gender = factor(df_mem$Gender)
df_adas$Gender = factor(df_adas$Gender)
df_avlt$Gender = factor(df_avlt$Gender)
df_wmh$Gender = factor(df_wmh$Gender)
df_av45$Gender = factor(df_av45$Gender)
df_hcicv$Gender = factor(df_hcicv$Gender)

# Only keep N/SMC/EMCI/LMCI
valid_diags = c('N','SMC','EMCI','LMCI')
#valid_diags = c('EMCI','LMCI')
#valid_diags = c('N','SMC')
df_mem = df_mem[which(df_mem$diag_prior %in% valid_diags),]
df_ef = df_ef[which(df_ef$diag_prior %in% valid_diags),]
df_avlt = df_avlt[which(df_avlt$diag_prior %in% valid_diags),]
df_adas = df_adas[which(df_adas$diag_prior %in% valid_diags),]
df_wmh = df_wmh[which(df_wmh$diag_prior %in% valid_diags),]
df_av45 = df_av45[which(df_av45$diag_prior %in% valid_diags),]
df_hcicv = df_hcicv[which(df_hcicv$diag_prior %in% valid_diags),]


# LM models        
fm_av45 = lm()


# LME MODELS: Normals 

fm_uwef_csonly = lmer(UW_EF_ ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_ef)
fm_uwmem_csonly = lmer(UW_MEM_ ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_mem)
fm_avlt_csonly = lmer(AVLT. ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_avlt)
fm_adas_csonly = lmer(ADAScog. ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_adas)
fm_wmh_csonly = lmer(WMH_percentOfICV. ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_wmh)
fm_av45_csonly = lmer(AV45 ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_av45)
fm_hcicv_csonly = lmer(100*FSX_HC.ICV_ ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_hcicv)

fm_uwef_cs = lmer(UW_EF_ ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_prior + CORTICAL_SUMMARY_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_ef)
fm_uwmem_cs = lmer(UW_MEM_ ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_prior + CORTICAL_SUMMARY_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_mem)
fm_avlt_cs = lmer(AVLT. ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_prior + CORTICAL_SUMMARY_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_avlt)
fm_adas_cs = lmer(ADAScog. ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_prior + CORTICAL_SUMMARY_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_adas)
fm_wmh_cs = lmer(WMH_percentOfICV. ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_prior + CORTICAL_SUMMARY_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_wmh)
fm_av45_cs = lmer(AV45 ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_prior + CORTICAL_SUMMARY_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_av45)
fm_hcicv_cs = lmer(100*FSX_HC.ICV_ ~ diag_prior + diag_prior:YrsPostBL + CORTICAL_SUMMARY_prior + CORTICAL_SUMMARY_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_hcicv)

fm_uwef = lmer(UW_EF_ ~ diag_prior + diag_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_ef)
fm_uwmem = lmer(UW_MEM_ ~ diag_prior + diag_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_mem)
fm_avlt = lmer(AVLT. ~ diag_prior + diag_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_avlt)
fm_adas = lmer(ADAScog. ~ diag_prior + diag_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_adas)
fm_wmh = lmer(WMH_percentOfICV. ~ diag_prior + diag_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_wmh)
fm_av45 = lmer(AV45 ~ diag_prior + diag_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_av45)
fm_hcicv = lmer(100*FSX_HC.ICV_ ~ diag_prior + diag_prior:YrsPostBL + X1 + X1:YrsPostBL + X4 + X4:YrsPostBL + X22 + X22:YrsPostBL + X25 + X25:YrsPostBL + X5 + X5:YrsPostBL + X3 + X3:YrsPostBL + X6 + X6:YrsPostBL + X9 + X9:YrsPostBL +X0 + X0:YrsPostBL + X34 + X34:YrsPostBL + X2 + X2:YrsPostBL + X19 + X19:YrsPostBL + X29 + X29:YrsPostBL + X56 + X56:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_hcicv)

# MODEL SUMMARIES
uwef_summary = summary(fm_uwef)
uwmem_summary = summary(fm_uwmem)
adas_summary = summary(fm_adas)
avlt_summary = summary(fm_avlt)
wmh_summary = summary(fm_wmh)
hcicv_summary = summary(fm_hcicv)
av45_summary = summary(fm_av45)

uwef_cs_summary = summary(fm_uwef_cs)
uwmem_cs_summary = summary(fm_uwmem_cs)
adas_cs_summary = summary(fm_adas_cs)
avlt_cs_summary = summary(fm_avlt_cs)
wmh_cs_summary = summary(fm_wmh_cs)
hcicv_cs_summary = summary(fm_hcicv_cs)
av45_cs_summary = summary(fm_av45_cs)

# MODEL ANOVAS
fm_avlt_cs = lmer(AVLT. ~ CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE:YrsPostBL + X22 + X22:YrsPostBL + APOE4_BIN + APOE4_BIN:YrsPostBL + Age.AV45 + Gender + Edu..Yrs. + YrsPostBL + (1 + YrsPostBL |RID), df_avlt)
anova(fm_avlt_csonly,fm_avlt_cs)

uwef_anova = anova(fm_uwef_csonly,fm_uwef_cs)
uwmem_anova = anova(fm_uwmem_csonly,fm_uwmem_cs)
adas_anova = anova(fm_adas_csonly,fm_adas_cs)
avlt_anova = anova(fm_avlt_csonly,fm_avlt_cs)
wmh_anova = anova(fm_wmh_csonly,fm_wmh_cs)
hcicv_anova = anova(fm_hcicv_csonly,fm_hcicv_cs)
av45_anova = anova(fm_av45_csonly,fm_av45_cs)

# latex COEFFICIENTS
uwmem_summary_coeff_latex = xtable(uwmem_summary$coefficients)
avlt_summary_coeff_latex = xtable(avlt_summary$coefficients)

# latex ANOVAS
uwmem_anova_latex = xtable(uwmem_anova)
avlt_anova_latex = xtable(avlt_anova)

# PRINT MODEL OUTPUTS
sink('uwef_summary.txt'); print(uwef_summary, correlation=TRUE); sink(file=NULL)
sink('uwmem_summary.txt'); print(uwmem_summary, correlation=TRUE); sink(file=NULL)
sink('adas_summary.txt'); print(adas_summary, correlation=TRUE); sink(file=NULL)
sink('avlt_summary.txt'); print(avlt_summary, correlation=TRUE); sink(file=NULL)
sink('av45_summary.txt'); print(av45_summary, correlation=TRUE); sink(file=NULL)
sink('hcicv_summary.txt'); print(hcicv_summary, correlation=TRUE); sink(file=NULL)
sink('wmh_summary.txt'); print(wmh_summary, correlation=TRUE); sink(file=NULL)

sink('uwef_cs_summary.txt'); print(uwef_cs_summary, correlation=TRUE); sink(file=NULL)
sink('uwmem_cs_summary.txt'); print(uwmem_cs_summary, correlation=TRUE); sink(file=NULL)
sink('adas_cs_summary.txt'); print(adas_cs_summary, correlation=TRUE); sink(file=NULL)
sink('avlt_cs_summary.txt'); print(avlt_cs_summary, correlation=TRUE); sink(file=NULL)
sink('av45_cs_summary.txt'); print(av45_cs_summary, correlation=TRUE); sink(file=NULL)
sink('hcicv_cs_summary.txt'); print(hcicv_cs_summary, correlation=TRUE); sink(file=NULL)
sink('wmh_cs_summary.txt'); print(wmh_cs_summary, correlation=TRUE); sink(file=NULL)

# PRINT MODEL ANOVA
sink('uwef_anova.txt'); print(uwef_anova, correlation=TRUE); sink(file=NULL)
sink('uwmem_anova.txt'); print(uwmem_anova, correlation=TRUE); sink(file=NULL)
sink('adas_anova.txt'); print(adas_anova, correlation=TRUE); sink(file=NULL)
sink('avlt_anova.txt'); print(avlt_anova, correlation=TRUE); sink(file=NULL)
sink('av45_anova.txt'); print(av45_anova, correlation=TRUE); sink(file=NULL)
sink('hcicv_anova.txt'); print(hcicv_anova, correlation=TRUE); sink(file=NULL)
sink('wmh_anova.txt'); print(wmh_anova, correlation=TRUE); sink(file=NULL)

# VISUALIZE MODELS
sjp.int(fm_avlt_cs,type='eff',swapPredictors=TRUE, showCI=TRUE)
sjp.int(fm_av45, type='eff', swapPredictors=TRUE, showCI=TRUE)

# Contrasts
contrasts_all = read.csv('contrasts_all.csv')
test = glht(fm_av45_cs,linfct=data.matrix(contrasts_all))
test_summary = summary(test)
contrasts = test_summary$linfct
coeff = matrix(test_summary$test$coefficients)
colnames(coeff) = 'estimate'
pvalue = matrix(test_summary$test$pvalues)
colnames(pvalue) = 'pvalues'
stderr = matrix(test_summary$test$sigma)
colnames(stderr) = 'stderr'
zvalue = matrix(test_summary$test$tstat)
colnames(zvalue) = 'zvalue'
contrasts = cbind(contrasts,coeff,stderr,zvalue,pvalue)
write.csv(contrasts, file = "av45_contrasts.csv", na = "")