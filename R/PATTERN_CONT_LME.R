library(lme4)
library(coefplot2)
library(ggplot2)
library(lmerTest)
library(pbkrtest)
library(multcomp)
library(contrast)

# Import data for Ns
df_long = read.csv('dpgmm_alpha12.66_bilateral_ADAScog__ALL_longdata.csv')

# Convert to factors
to_factor = c('RID','APOE4_BIN','APOE2_BIN','diag_prior','Gender')
df_long[to_factor] = lapply(df_long[to_factor], as.factor)

# Separate Diag groups
df_long_n = df_long[which(df_long$diag_prior %in% c('N','SMC')),] 
df_long_mci = df_long[which(df_long$diag_prior %in% c('EMCI','LMCI')),] 
df_long_ad = df_long[which(df_long$diag_prior %in% c('AD')),] 

# LME MODELS
fm_cs_n = lmer(ADAScog. ~ CORTICAL_SUMMARY_prior + Age.AV45 + Gender + Edu..Yrs. + diag_prior + APOE4_BIN*YrsPostBL + (1 + YrsPostBL| RID), df_long_n)
fm_pattern_n = lmer(ADAScog. ~ CORTICAL_SUMMARY_prior + X0 + X1 + X3 + X4 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X20 + X21 + X23 + X24 + Age.AV45 + Gender + Edu..Yrs. + diag_prior + APOE4_BIN*YrsPostBL + (1 + YrsPostBL| RID), df_long_n)
fm_cs_mci = lmer(ADAScog. ~ CORTICAL_SUMMARY_prior + Age.AV45 + Gender + Edu..Yrs. + diag_prior + APOE4_BIN*YrsPostBL + (1 + YrsPostBL | RID), df_long_mci)
fm_pattern_mci = lmer(ADAScog. ~ CORTICAL_SUMMARY_prior + X0 + X1 + X3 + X4 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X20 + X21 + X23 + X24 + Age.AV45 + Gender + Edu..Yrs. + diag_prior + APOE4_BIN*YrsPostBL + (1 + YrsPostBL | RID), df_long_mci)
fm_cs_ad = lmer(ADAScog. ~ CORTICAL_SUMMARY_prior + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN*YrsPostBL + (1 + YrsPostBL | RID), df_long_ad)
fm_pattern_ad = lmer(ADAScog. ~ CORTICAL_SUMMARY_prior + X0 + X1 + X3 + X4 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X20 + X21 + X23 + X24 + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN*YrsPostBL + (1 + YrsPostBL | RID), df_long_ad)

# MODEL SUMMARIES
fm_cs_n_summary = summary(fm_cs_n)
fm_pattern_n_summary = summary(fm_pattern_n)
fm_cs_mci_summary = summary(fm_cs_mci)
fm_pattern_mci_summary = summary(fm_pattern_mci)
fm_cs_ad_summary = summary(fm_cs_ad)
fm_pattern_ad_summary = summary(fm_pattern_ad)

# MODEL ANOVAS
fm_cs_n_anova = as.matrix(anova(fm_cs_n))
fm_pattern_n_anova = as.matrix(anova(fm_pattern_n))
fm_cs_mci_anova = as.matrix(anova(fm_cs_mci))
fm_pattern_mci_anova = as.matrix(anova(fm_pattern_mci))
fm_cs_ad_anova = as.matrix(anova(fm_cs_ad))
fm_pattern_ad_anova = as.matrix(anova(fm_pattern_ad))


# PRINT MODEL OUTPUTS
sink('fm_cs_n_summary.txt'); print(fm_cs_n_summary, correlation=TRUE); sink(file=NULL)
sink('fm_pattern_n_summary.txt'); print(fm_pattern_n_summary, correlation=TRUE); sink(file=NULL)
sink('fm_cs_mci_summary.txt'); print(fm_cs_mci_summary, correlation=TRUE); sink(file=NULL)
sink('fm_pattern_mci_summary.txt'); print(fm_pattern_mci_summary, correlation=TRUE); sink(file=NULL)
sink('fm_cs_ad_summary.txt'); print(fm_cs_ad_summary, correlation=TRUE); sink(file=NULL)
sink('fm_pattern_ad_summary.txt'); print(fm_pattern_ad_summary, correlation=TRUE); sink(file=NULL)

# PRINT MODEL COEFFICIENTS
write.csv(fm_cs_n_summary$coefficients, file = "fm_cs_n_coefficients.csv", na = "")
write.csv(fm_pattern_n_summary$coefficients, file = "fm_pattern_n_coefficients.csv", na = "")
write.csv(fm_cs_mci_summary$coefficients, file = "fm_cs_mci_coefficients.csv", na = "")
write.csv(fm_pattern_mci_summary$coefficients, file = "fm_pattern_mci_coefficients.csv", na = "")
write.csv(fm_cs_ad_summary$coefficients, file = "fm_cs_ad_coefficients.csv", na = "")
write.csv(fm_pattern_ad_summary$coefficients, file = "fm_pattern_ad_coefficients.csv", na = "")

# PRINT MODEL ANOVA
write.csv(fm_cs_n_anova, file = "fm_cs_n_anova.csv", na = "")
write.csv(fm_pattern_n_anova, file = "fm_pattern_n_anova.csv", na = "")
write.csv(fm_cs_mci_anova, file = "fm_cs_mci_anova.csv", na = "")
write.csv(fm_pattern_mci_anova, file = "fm_pattern_mci_anova.csv", na = "")
write.csv(fm_cs_ad_anova, file = "fm_cs_ad_anova.csv", na = "")
write.csv(fm_pattern_ad_anova, file = "fm_pattern_ad_anova.csv", na = "")


