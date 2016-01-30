fm_uwef_n_cs = lmer(UW_EF_ ~ APOE4_BIN + APOE4_BIN*UW_EF_postAV45_ + Age.AV45 + Gender + Edu..Yrs. + UW_EF_postAV45_ + (1 + UW_EF_postAV45_ | RID), df_ef_n_negative)
fm_uwmem_n_cs = lmer(UW_MEM_ ~ APOE4_BIN + APOE4_BIN*UW_MEM_postAV45_ + Age.AV45 + Gender + Edu..Yrs. + UW_MEM_postAV45_ + (1 + UW_MEM_postAV45_ | RID), df_mem_n_negative)
fm_avlt_n_cs = lmer(AVLT. ~ APOE4_BIN + APOE4_BIN*TIMEpostAV45_AVLT. + Age.AV45 + Gender + Edu..Yrs. + TIMEpostAV45_AVLT. + (1 + TIMEpostAV45_AVLT. | RID), df_avlt_n_negative)
fm_adas_n_cs = lmer(ADAScog. ~ APOE4_BIN + APOE4_BIN*TIMEpostAV45_ADAS. + Age.AV45 + Gender + Edu..Yrs. + TIMEpostAV45_ADAS. + (1 + TIMEpostAV45_ADAS. | RID), df_adas_n_negative)



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

