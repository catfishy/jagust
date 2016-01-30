library(lme4)
library(coefplot2)
library(ggplot2)
library(lmerTest)
library(pbkrtest)
library(multcomp)
library(contrast)

# Import data for Ns
df_ef_n = read.csv('dpgmm_alpha15.37_bilateral_UW_EF__N_longdata_continuous.csv')
df_mem_n = read.csv('dpgmm_alpha15.37_bilateral_UW_MEM__N_longdata_continuous.csv')
df_adas_n = read.csv('dpgmm_alpha15.37_bilateral_ADAScog__N_longdata_continuous.csv')
df_avlt_n = read.csv('dpgmm_alpha15.37_bilateral_AVLT__N_longdata_continuous.csv')
df_ef_n$RID = factor(df_ef_n$RID)
df_mem_n$RID = factor(df_mem_n$RID)
df_adas_n$RID = factor(df_adas_n$RID)
df_avlt_n$RID = factor(df_avlt_n$RID)
df_ef_n$group = factor(df_ef_n$group)
df_mem_n$group = factor(df_mem_n$group)
df_adas_n$group = factor(df_adas_n$group)
df_avlt_n$group = factor(df_avlt_n$group)
df_ef_n$diag_prior = factor(df_ef_n$diag_prior)
df_mem_n$diag_prior = factor(df_mem_n$diag_prior)
df_adas_n$diag_prior = factor(df_adas_n$diag_prior)
df_avlt_n$diag_prior = factor(df_avlt_n$diag_prior)
df_ef_n$APOE4_BIN = factor(df_ef_n$APOE4_BIN)
df_mem_n$APOE4_BIN = factor(df_mem_n$APOE4_BIN)
df_adas_n$APOE4_BIN = factor(df_adas_n$APOE4_BIN)
df_avlt_n$APOE4_BIN = factor(df_avlt_n$APOE4_BIN)
df_ef_n$Gender = factor(df_ef_n$Gender)
df_mem_n$Gender = factor(df_mem_n$Gender)
df_adas_n$Gender = factor(df_adas_n$Gender)
df_avlt_n$Gender = factor(df_avlt_n$Gender)
df_ef_n$CORTICAL_SUMMARY_POSITIVE = factor(df_ef_n$CORTICAL_SUMMARY_POSITIVE)
df_mem_n$CORTICAL_SUMMARY_POSITIVE = factor(df_mem_n$CORTICAL_SUMMARY_POSITIVE)
df_adas_n$CORTICAL_SUMMARY_POSITIVE = factor(df_adas_n$CORTICAL_SUMMARY_POSITIVE)
df_avlt_n$CORTICAL_SUMMARY_POSITIVE = factor(df_avlt_n$CORTICAL_SUMMARY_POSITIVE)

df_adas_n_positive = df_adas_n[df_adas_n$CORTICAL_SUMMARY_POSITIVE==1,]
df_ef_n_positive = df_ef_n[df_ef_n$CORTICAL_SUMMARY_POSITIVE==1,]
df_mem_n_positive = df_mem_n[df_mem_n$CORTICAL_SUMMARY_POSITIVE==1,]
df_avlt_n_positive = df_avlt_n[df_avlt_n$CORTICAL_SUMMARY_POSITIVE==1,]
df_adas_n_negative = df_adas_n[df_adas_n$CORTICAL_SUMMARY_POSITIVE==0,]
df_ef_n_negative = df_ef_n[df_ef_n$CORTICAL_SUMMARY_POSITIVE==0,]
df_mem_n_negative = df_mem_n[df_mem_n$CORTICAL_SUMMARY_POSITIVE==0,]
df_avlt_n_negative = df_avlt_n[df_avlt_n$CORTICAL_SUMMARY_POSITIVE==0,]

# Import data for MCIs
df_ef_mci = read.csv('dpgmm_alpha15.37_bilateral_UW_EF__MCI_longdata_continuous.csv')
df_mem_mci = read.csv('dpgmm_alpha15.37_bilateral_UW_MEM__MCI_longdata_continuous.csv')
df_adas_mci = read.csv('dpgmm_alpha15.37_bilateral_ADAScog__MCI_longdata_continuous.csv')
df_avlt_mci = read.csv('dpgmm_alpha15.37_bilateral_AVLT__MCI_longdata_continuous.csv')
df_ef_mci$RID = factor(df_ef_mci$RID)
df_mem_mci$RID = factor(df_mem_mci$RID)
df_adas_mci$RID = factor(df_adas_mci$RID)
df_avlt_mci$RID = factor(df_avlt_mci$RID)
df_ef_mci$group = factor(df_ef_mci$group)
df_mem_mci$group = factor(df_mem_mci$group)
df_adas_mci$group = factor(df_adas_mci$group)
df_avlt_mci$group = factor(df_avlt_mci$group)
df_ef_mci$diag_prior = factor(df_ef_mci$diag_prior)
df_mem_mci$diag_prior = factor(df_mem_mci$diag_prior)
df_adas_mci$diag_prior = factor(df_adas_mci$diag_prior)
df_avlt_mci$diag_prior = factor(df_avlt_mci$diag_prior)
df_ef_mci$APOE4_BIN = factor(df_ef_mci$APOE4_BIN)
df_mem_mci$APOE4_BIN = factor(df_mem_mci$APOE4_BIN)
df_adas_mci$APOE4_BIN = factor(df_adas_mci$APOE4_BIN)
df_avlt_mci$APOE4_BIN = factor(df_avlt_mci$APOE4_BIN)
df_ef_mci$Gender = factor(df_ef_mci$Gender)
df_mem_mci$Gender = factor(df_mem_mci$Gender)
df_adas_n$Gender = factor(df_adas_n$Gender)
df_avlt_n$Gender = factor(df_avlt_n$Gender)
df_ef_mci$CORTICAL_SUMMARY_POSITIVE = factor(df_ef_mci$CORTICAL_SUMMARY_POSITIVE)
df_mem_mci$CORTICAL_SUMMARY_POSITIVE = factor(df_mem_mci$CORTICAL_SUMMARY_POSITIVE)
df_adas_mci$CORTICAL_SUMMARY_POSITIVE = factor(df_adas_mci$CORTICAL_SUMMARY_POSITIVE)
df_avlt_mci$CORTICAL_SUMMARY_POSITIVE = factor(df_avlt_mci$CORTICAL_SUMMARY_POSITIVE)

df_adas_mci_positive = df_adas_mci[df_adas_mci$CORTICAL_SUMMARY_POSITIVE==1,]
df_ef_mci_positive = df_ef_mci[df_ef_mci$CORTICAL_SUMMARY_POSITIVE==1,]
df_mem_mci_positive = df_mem_mci[df_mem_mci$CORTICAL_SUMMARY_POSITIVE==1,]
df_avlt_mci_positive = df_avlt_mci[df_avlt_mci$CORTICAL_SUMMARY_POSITIVE==1,]
df_adas_mci_negative = df_adas_mci[df_adas_mci$CORTICAL_SUMMARY_POSITIVE==0,]
df_ef_mci_negative = df_ef_mci[df_ef_mci$CORTICAL_SUMMARY_POSITIVE==0,]
df_mem_mci_negative = df_mem_mci[df_mem_mci$CORTICAL_SUMMARY_POSITIVE==0,]
df_avlt_mci_negative = df_avlt_mci[df_avlt_mci$CORTICAL_SUMMARY_POSITIVE==0,]

# LME MODELS: Normals 
fm_uwef_n_cs = lmer(UW_EF_ ~ CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE*UW_EF_postAV45_ + APOE4_BIN + APOE4_BIN*UW_EF_postAV45_ + Age.AV45 + Gender + Edu..Yrs. + UW_EF_postAV45_ + (1 + UW_EF_postAV45_ | RID), df_ef_n)
fm_uwmem_n_cs = lmer(UW_MEM_ ~ CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE*UW_MEM_postAV45_ + APOE4_BIN + APOE4_BIN*UW_MEM_postAV45_ + Age.AV45 + Gender + Edu..Yrs. + UW_MEM_postAV45_ + (1 + UW_MEM_postAV45_ | RID), df_mem_n)
fm_avlt_n_cs = lmer(AVLT. ~ CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE*TIMEpostAV45_AVLT. + APOE4_BIN + APOE4_BIN*TIMEpostAV45_AVLT. + Age.AV45 + Gender + Edu..Yrs. + TIMEpostAV45_AVLT. + (1 + TIMEpostAV45_AVLT. | RID), df_avlt_n)
fm_adas_n_cs = lmer(ADAScog. ~ CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE*TIMEpostAV45_ADAS. + APOE4_BIN + APOE4_BIN*TIMEpostAV45_ADAS. + Age.AV45 + Gender + Edu..Yrs. + TIMEpostAV45_ADAS. + (1 + TIMEpostAV45_ADAS. | RID), df_adas_n)

fm_uwef_n = lmer(UW_EF_ ~ X1 + X1*UW_EF_postAV45_ + X4 + X4*UW_EF_postAV45_ + X2 + X2*UW_EF_postAV45_ + X19 + X19*UW_EF_postAV45_ + X22 + X22*UW_EF_postAV45_ + X25 + X25*UW_EF_postAV45_ + X29 + X29*UW_EF_postAV45_ + X34 + X34*UW_EF_postAV45_ + X56 + X56*UW_EF_postAV45_ + APOE4_BIN + APOE4_BIN*UW_EF_postAV45_ + Age.AV45 + Gender + Edu..Yrs. + UW_EF_postAV45_ + (1 + UW_EF_postAV45_ | RID), df_ef_n)
fm_uwmem_n = lmer(UW_MEM_ ~ X1 + X1*UW_MEM_postAV45_ + X4 + X4*UW_MEM_postAV45_ + X2 + X2*UW_MEM_postAV45_ + X19 + X19*UW_MEM_postAV45_ + X22 + X22*UW_MEM_postAV45_ + X25 + X25*UW_MEM_postAV45_ + X29 + X29*UW_MEM_postAV45_ + X34 + X34*UW_MEM_postAV45_ + X56 + X56*UW_MEM_postAV45_ + APOE4_BIN + APOE4_BIN*UW_MEM_postAV45_ + Age.AV45 + Gender + Edu..Yrs. + UW_MEM_postAV45_ + (1 + UW_MEM_postAV45_ | RID), df_mem_n)
fm_avlt_n = lmer(AVLT. ~ CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE*TIMEpostAV45_AVLT. + X1 + X1*TIMEpostAV45_AVLT. + X4 + X4*TIMEpostAV45_AVLT. + X2 + X2*TIMEpostAV45_AVLT. + X19 + X19*TIMEpostAV45_AVLT. + X22 + X22*TIMEpostAV45_AVLT. + X25 + X25*TIMEpostAV45_AVLT. + X29 + X29*TIMEpostAV45_AVLT. + X34 + X34*TIMEpostAV45_AVLT. + X56 + X56*TIMEpostAV45_AVLT. + APOE4_BIN + APOE4_BIN*TIMEpostAV45_AVLT. + Age.AV45 + Gender + Edu..Yrs. + TIMEpostAV45_AVLT. + (1 + TIMEpostAV45_AVLT. | RID), df_avlt_n)
fm_adas_n = lmer(ADAScog. ~ X1 + X1*TIMEpostAV45_ADAS. + X4 + X4*TIMEpostAV45_ADAS. + X2 + X2*TIMEpostAV45_ADAS. + X19 + X19*TIMEpostAV45_ADAS. + X22 + X22*TIMEpostAV45_ADAS. + X25 + X25*TIMEpostAV45_ADAS. + X29 + X29*TIMEpostAV45_ADAS. + X34 + X34*TIMEpostAV45_ADAS. + X56 + X56*TIMEpostAV45_ADAS. + APOE4_BIN + APOE4_BIN*TIMEpostAV45_ADAS. + Age.AV45 + Gender + Edu..Yrs. + TIMEpostAV45_ADAS. + (1 + TIMEpostAV45_ADAS. | RID), df_adas_n)

# LME MODELS: MCIS
fm_uwef_mci_cs = lmer(UW_EF_ ~ CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE*UW_EF_postAV45_ + APOE4_BIN + APOE4_BIN*UW_EF_postAV45_ + Age.AV45 + Gender + Edu..Yrs. + UW_EF_postAV45_ + (1 + UW_EF_postAV45_ | RID), df_ef_mci)
fm_uwmem_mci_cs = lmer(UW_MEM_ ~ CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE*UW_MEM_postAV45_ + APOE4_BIN + APOE4_BIN*UW_MEM_postAV45_ + Age.AV45 + Gender + Edu..Yrs. + UW_MEM_postAV45_ + (1 + UW_MEM_postAV45_ | RID), df_mem_mci)
fm_avlt_mci_cs = lmer(AVLT. ~ CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE*TIMEpostAV45_AVLT. + APOE4_BIN + APOE4_BIN*TIMEpostAV45_AVLT. + Age.AV45 + Gender + Edu..Yrs. + TIMEpostAV45_AVLT. + (1 + TIMEpostAV45_AVLT. | RID), df_avlt_mci)
fm_adas_mci_cs = lmer(ADAScog. ~ CORTICAL_SUMMARY_POSITIVE + CORTICAL_SUMMARY_POSITIVE*TIMEpostAV45_ADAS. + APOE4_BIN + APOE4_BIN*TIMEpostAV45_ADAS. + Age.AV45 + Gender + Edu..Yrs. + TIMEpostAV45_ADAS. + (1 + TIMEpostAV45_ADAS. | RID), df_adas_mci)

fm_uwef_mci = lmer(UW_EF_ ~ X1 + X1*UW_EF_postAV45_ + X4 + X4*UW_EF_postAV45_ + X2 + X2*UW_EF_postAV45_ + X19 + X19*UW_EF_postAV45_ + X22 + X22*UW_EF_postAV45_ + X25 + X25*UW_EF_postAV45_ + X29 + X29*UW_EF_postAV45_ + X34 + X34*UW_EF_postAV45_ + X56 + X56*UW_EF_postAV45_ + APOE4_BIN + APOE4_BIN*UW_EF_postAV45_ + Age.AV45 + Gender + Edu..Yrs. + UW_EF_postAV45_ + (1 + UW_EF_postAV45_ | RID), df_ef_mci)
fm_uwmem_mci = lmer(UW_MEM_ ~ X1 + X1*UW_MEM_postAV45_ + X4 + X4*UW_MEM_postAV45_ + X2 + X2*UW_MEM_postAV45_ + X19 + X19*UW_MEM_postAV45_ + X22 + X22*UW_MEM_postAV45_ + X25 + X25*UW_MEM_postAV45_ + X29 + X29*UW_MEM_postAV45_ + X34 + X34*UW_MEM_postAV45_ + X56 + X56*UW_MEM_postAV45_ + APOE4_BIN + APOE4_BIN*UW_MEM_postAV45_ + Age.AV45 + Gender + Edu..Yrs. + UW_MEM_postAV45_ + (1 + UW_MEM_postAV45_ | RID), df_mem_mci)
fm_avlt_mci = lmer(AVLT. ~ X1 + X1*TIMEpostAV45_AVLT. + X4 + X4*TIMEpostAV45_AVLT. + X2 + X2*TIMEpostAV45_AVLT. + X19 + X19*TIMEpostAV45_AVLT. + X22 + X22*TIMEpostAV45_AVLT. + X25 + X25*TIMEpostAV45_AVLT. + X29 + X29*TIMEpostAV45_AVLT. + X34 + X34*TIMEpostAV45_AVLT. + X56 + X56*TIMEpostAV45_AVLT. + APOE4_BIN + APOE4_BIN*TIMEpostAV45_AVLT. + Age.AV45 + Gender + Edu..Yrs. + TIMEpostAV45_AVLT. + (1 + TIMEpostAV45_AVLT. | RID), df_avlt_mci)
fm_adas_mci = lmer(ADAScog. ~ X1 + X1*TIMEpostAV45_ADAS. + X4 + X4*TIMEpostAV45_ADAS. + X2 + X2*TIMEpostAV45_ADAS. + X19 + X19*TIMEpostAV45_ADAS. + X22 + X22*TIMEpostAV45_ADAS. + X25 + X25*TIMEpostAV45_ADAS. + X29 + X29*TIMEpostAV45_ADAS. + X34 + X34*TIMEpostAV45_ADAS. + X56 + X56*TIMEpostAV45_ADAS. + APOE4_BIN + APOE4_BIN*TIMEpostAV45_ADAS. + Age.AV45 + Gender + Edu..Yrs. + TIMEpostAV45_ADAS. + (1 + TIMEpostAV45_ADAS. | RID), df_adas_mci)

# MODEL SUMMARIES
uwef_n_summary = summary(fm_uwef_n)
uwmem_n_summary = summary(fm_uwmem_n)
adas_n_summary = summary(fm_adas_n)
avlt_n_summary = summary(fm_avlt_n)

uwef_n_cs_summary = summary(fm_uwef_n_cs)
uwmem_n_cs_summary = summary(fm_uwmem_n_cs)
adas_n_cs_summary = summary(fm_adas_n_cs)
avlt_n_cs_summary = summary(fm_avlt_n_cs)

uwef_mci_summary = summary(fm_uwef_mci)
uwmem_mci_summary = summary(fm_uwmem_mci)
adas_mci_summary = summary(fm_adas_mci)
avlt_mci_summary = summary(fm_avlt_mci)

uwef_mci_cs_summary = summary(fm_uwef_mci_cs)
uwmem_mci_cs_summary = summary(fm_uwmem_mci_cs)
adas_mci_cs_summary = summary(fm_adas_mci_cs)
avlt_mci_cs_summary = summary(fm_avlt_mci_cs)

# MODEL ANOVAS
uwef_n_anova = as.matrix(anova(fm_uwef_n))
uwmem_n_anova = as.matrix(anova(fm_uwmem_n))
adas_n_anova = as.matrix(anova(fm_adas_n))
avlt_n_anova = as.matrix(anova(fm_avlt_n))
uwef_n_cs_anova = as.matrix(anova(fm_uwef_n_cs))
uwmem_n_cs_anova = as.matrix(anova(fm_uwmem_n_cs))
adas_n_cs_anova = as.matrix(anova(fm_adas_n_cs))
avlt_n_cs_anova = as.matrix(anova(fm_avlt_n_cs))

uwef_mci_anova = as.matrix(anova(fm_uwef_mci))
uwmem_mci_anova = as.matrix(anova(fm_uwmem_mci))
adas_mci_anova = as.matrix(anova(fm_adas_mci))
avlt_mci_anova = as.matrix(anova(fm_avlt_mci))
uwef_mci_cs_anova = as.matrix(anova(fm_uwef_mci_cs))
uwmem_mci_cs_anova = as.matrix(anova(fm_uwmem_mci_cs))
adas_mci_cs_anova = as.matrix(anova(fm_adas_mci_cs))
avlt_mci_cs_anova = as.matrix(anova(fm_avlt_mci_cs))

# PRINT MODEL OUTPUTS
sink('uwef_n_summary.txt'); print(uwef_n_summary, correlation=TRUE); sink(file=NULL)
sink('uwmem_n_summary.txt'); print(uwmem_n_summary, correlation=TRUE); sink(file=NULL)
sink('adas_n_summary.txt'); print(adas_n_summary, correlation=TRUE); sink(file=NULL)
sink('avlt_n_summary.txt'); print(avlt_n_summary, correlation=TRUE); sink(file=NULL)
sink('uwef_n_cs_summary.txt'); print(uwef_n_cs_summary, correlation=TRUE); sink(file=NULL)
sink('uwmem_n_cs_summary.txt'); print(uwmem_n_cs_summary, correlation=TRUE); sink(file=NULL)
sink('adas_n_cs_summary.txt'); print(adas_n_cs_summary, correlation=TRUE); sink(file=NULL)
sink('avlt_n_cs_summary.txt'); print(avlt_n_cs_summary, correlation=TRUE); sink(file=NULL)

sink('uwef_mci_summary.txt'); print(uwef_mci_summary, correlation=TRUE); sink(file=NULL)
sink('uwmem_mci_summary.txt'); print(uwmem_mci_summary, correlation=TRUE); sink(file=NULL)
sink('adas_mci_summary.txt'); print(adas_mci_summary, correlation=TRUE); sink(file=NULL)
sink('avlt_mci_summary.txt'); print(avlt_mci_summary, correlation=TRUE); sink(file=NULL)
sink('uwef_mci_cs_summary.txt'); print(uwef_mci_cs_summary, correlation=TRUE); sink(file=NULL)
sink('uwmem_mci_cs_summary.txt'); print(uwmem_mci_cs_summary, correlation=TRUE); sink(file=NULL)
sink('adas_mci_cs_summary.txt'); print(adas_mci_cs_summary, correlation=TRUE); sink(file=NULL)
sink('avlt_mci_cs_summary.txt'); print(avlt_mci_cs_summary, correlation=TRUE); sink(file=NULL)

# PRINT MODEL COEFFICIENTS
write.csv(uwef_n_summary$coefficients, file = "uwef_n_coefficients.csv", na = "")
write.csv(uwmem_n_summary$coefficients, file = "uwmem_n_coefficients.csv", na = "")
write.csv(adas_n_summary$coefficients, file = "adas_n_coefficients.csv", na = "")
write.csv(avlt_n_summary$coefficients, file = "avlt_n_coefficients.csv", na = "")
write.csv(uwef_n_cs_summary$coefficients, file = "uwef_n_cs_coefficients.csv", na = "")
write.csv(uwmem_n_cs_summary$coefficients, file = "uwmem_n_cs_coefficients.csv", na = "")
write.csv(adas_n_cs_summary$coefficients, file = "adas_n_cs_coefficients.csv", na = "")
write.csv(avlt_n_cs_summary$coefficients, file = "avlt_n_cs_coefficients.csv", na = "")

write.csv(uwef_mci_summary$coefficients, file = "uwef_mci_coefficients.csv", na = "")
write.csv(uwmem_mci_summary$coefficients, file = "uwmem_mci_coefficients.csv", na = "")
write.csv(adas_mci_summary$coefficients, file = "adas_mci_coefficients.csv", na = "")
write.csv(avlt_mci_summary$coefficients, file = "avlt_mci_coefficients.csv", na = "")
write.csv(uwef_mci_cs_summary$coefficients, file = "uwef_mci_cs_coefficients.csv", na = "")
write.csv(uwmem_mci_cs_summary$coefficients, file = "uwmem_mci_cs_coefficients.csv", na = "")
write.csv(adas_mci_cs_summary$coefficients, file = "adas_mci_cs_coefficients.csv", na = "")
write.csv(avlt_mci_cs_summary$coefficients, file = "avlt_mci_cs_coefficients.csv", na = "")

# PRINT MODEL ANOVA
write.csv(uwef_n_anova, file = "uwef_n_anova.csv", na = "")
write.csv(uwmem_n_anova, file = "uwmem_n_anova.csv", na = "")
write.csv(adas_n_anova, file = "adas_n_anova.csv", na = "")
write.csv(avlt_n_anova, file = "avlt_n_anova.csv", na = "")
write.csv(uwef_n_cs_anova, file = "uwef_n_cs_anova.csv", na = "")
write.csv(uwmem_n_cs_anova, file = "uwmem_n_cs_anova.csv", na = "")
write.csv(adas_n_cs_anova, file = "adas_n_cs_anova.csv", na = "")
write.csv(avlt_n_cs_anova, file = "avlt_n_cs_anova.csv", na = "")

write.csv(uwef_mci_anova, file = "uwef_mci_anova.csv", na = "")
write.csv(uwmem_mci_anova, file = "uwmem_mci_anova.csv", na = "")
write.csv(adas_mci_anova, file = "adas_mci_anova.csv", na = "")
write.csv(avlt_mci_anova, file = "avlt_mci_anova.csv", na = "")
write.csv(uwef_mci_cs_anova, file = "uwef_mci_cs_anova.csv", na = "")
write.csv(uwmem_mci_cs_anova, file = "uwmem_mci_cs_anova.csv", na = "")
write.csv(adas_mci_cs_anova, file = "adas_mci_cs_anova.csv", na = "")
write.csv(avlt_mci_cs_anova, file = "avlt_mci_cs_anova.csv", na = "")


# VISUALIZE MODELS

