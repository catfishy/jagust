library(lme4)
library(coefplot2)
library(ggplot2)
library(lmerTest)
library(pbkrtest)
library(multcomp)
library(contrast)
library(xtable)
library(sjPlot)
library(splines)
library(car)
library(stats)

# Import data
df_av45 = read.csv('dpgmm_alpha14.36_bilateral_AV45_ALL_longdata_continuous_slope.csv')

# Convert to factors
df_av45$CORTICAL_SUMMARY_POSITIVE = factor(df_av45$CORTICAL_SUMMARY_POSITIVE)
df_av45$RID = factor(df_av45$RID)
df_av45$group = factor(df_av45$group)
df_av45$diag_prior = factor(df_av45$diag_prior)
df_av45$APOE4_BIN = factor(df_av45$APOE4_BIN)
df_av45$Gender = factor(df_av45$Gender)

# add indicator variables
df_av45$AD = factor(as.numeric(df_av45$diag_prior == 'AD'))

# patterns used
valid_patterns = c(16,4,3,6,19,0,8,7)

# pattern weight models
fm_av45_cs = glm(AD ~ CORTICAL_SUMMARY_prior, family=binomial, df_av45)
fm_av45_pattern = glm(AD ~ X0 + X3 + X4 + X6 + X7 + X8 + X16 + X19, family=binomial, df_av45)


# summaries
fm_av45_cs_summary = summary(fm_av45_cs)
fm_av45_pattern_summary = summary(fm_av45_pattern)

fm_av45_cs_summary$aic
fm_av45_pattern_summary$aic

# anova
fm_av45_nopattern_anova = Anova(fm_av45_nopattern,type='III')
fm_av45_anova = Anova(fm_av45,type='III')
fm_av45_onlypatterns_anova = Anova(fm_av45_onlypatterns,type='III')
