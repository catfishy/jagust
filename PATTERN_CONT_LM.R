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

# Import data for Ns
df_av45 = read.csv('dpgmm_alpha15.37_bilateral_AV45_ALL_longdata_continuous_slope.csv')
df_av45$CORTICAL_SUMMARY_POSITIVE = factor(df_av45$CORTICAL_SUMMARY_POSITIVE)
df_av45$RID = factor(df_av45$RID)
df_av45$group = factor(df_av45$group)
df_av45$diag_prior = factor(df_av45$diag_prior)
df_av45$APOE4_BIN = factor(df_av45$APOE4_BIN)
df_av45$Gender = factor(df_av45$Gender)

# Only keep N/SMC/EMCI/LMCI
valid_diags = c('N','SMC','EMCI','LMCI')
#valid_diags = c('EMCI','LMCI')
#valid_diags = c('N','SMC')
df_av45 = df_av45[which(df_av45$diag_prior %in% valid_diags),]


# PCA on patterns

# LM models
fm_av45_nopattern = lm(AV45_slope ~ diag_prior + poly(CORTICAL_SUMMARY_prior,2) + APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., df_av45)
fm_av45 = lm(AV45_slope ~ diag_prior + poly(CORTICAL_SUMMARY_prior,2) + X1 + X4 + X22 + X25 + X5 + X3 + X6 + X9 + X0 + X34 + X2 + X19 + X29 + X56 + APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., df_av45)

# splined LM models
fm_av45_nopattern_spline = lm(AV45_slope ~ diag_prior + ns(CORTICAL_SUMMARY_prior,df=5) + APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., df_av45)
fm_av45_spline = lm(AV45_slope ~ diag_prior + ns(CORTICAL_SUMMARY_prior,df=5) + X1 + X4 + X22 + X25 + X5 + X3 + X6 + X9 + X0 + X34 + X2 + X19 + X29 + X56 + APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., df_av45)

# summaries
fm_av45_nopattern_summary = summary(fm_av45_nopattern)
fm_av45_summary = summary(fm_av45)
fm_av45_nopattern_spline_summary = summary(fm_av45_nopattern_spline)
fm_av45_spline_summary = summary(fm_av45_spline)

fm_av45_nopattern_summary$r.squared
fm_av45_summary$r.squared
fm_av45_nopattern_spline_summary$r.squared
fm_av45_spline_summary$r.squared

# anova
fm_av45_nopattern_anova = anova(fm_av45_nopattern)
fm_av45_anova = anova(fm_av45)
fm_av45_nopattern_spline_anova = anova(fm_av45_nopattern_spline)
fm_av45_spline_anova = anova(fm_av45_spline)

# anova model comparisons
fm_modelcomparison_anova = anova(fm_av45_nopattern, fm_av45)
fm_modelcomparison_spline_anova = anova(fm_av45_nopattern_spline, fm_av45_spline)

# PRINT MODEL OUTPUTS
sink('av45_lm_nopattern.txt'); print(fm_av45_nopattern_summary, correlation=TRUE); sink(file=NULL)
sink('av45_lm_withpattern.txt'); print(fm_av45_summary, correlation=TRUE); sink(file=NULL)
sink('av45_spline_lm_nopattern.txt'); print(fm_av45_nopattern_spline_summary, correlation=TRUE); sink(file=NULL)
sink('av45_spline_lm_withpattern.txt'); print(fm_av45_spline_summary, correlation=TRUE); sink(file=NULL)

# PRINT MODEL ANOVA
sink('av45_lm_nopattern_anova.txt'); print(fm_av45_nopattern_anova, correlation=TRUE); sink(file=NULL)
sink('av45_lm_withpattern_anova.txt'); print(fm_av45_anova, correlation=TRUE); sink(file=NULL)
sink('av45_spline_lm_nopattern_anova.txt'); print(fm_av45_nopattern_spline_anova, correlation=TRUE); sink(file=NULL)
sink('av45_spline_lm_withpattern_anova.txt'); print(fm_av45_spline_anova, correlation=TRUE); sink(file=NULL)
sink('av45_mc_anova.txt'); print(fm_modelcomparison_anova, correlation=TRUE); sink(file=NULL)
sink('av45_spline_mc_anova.txt'); print(fm_modelcomparison_spline_anova, correlation=TRUE); sink(file=NULL)


# plot
plot(df_av45[,'CORTICAL_SUMMARY_prior'],df_av45[,'AV45_slope'], main='Original Data', xlab='Cortical Summary BL', ylab='Annualized AV45 Slope')
plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45), main='LM Predicted (with patterns)', xlab='Cortical Summary BL', ylab='Annualized AV45 Slope')
plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_nopattern), main='LM Predicted (without patterns)', xlab='Cortical Summary BL', ylab='Annualized AV45 Slope')
plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_spline), main='Cubic Spline LM Predicted (with patterns)', xlab='Cortical Summary BL', ylab='Annualized AV45 Slope')
plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_nopattern_spline), main='Cubic Spline LM Predicted (without patterns)', xlab='Cortical Summary BL', ylab='Annualized AV45 Slope')

# contrasts
contrasts_all = read.csv('contrasts_lm.csv')
test = glht(fm_av45_spline,linfct=data.matrix(contrasts_all))
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
write.csv(contrasts, file = "av45_splinelm_contrasts.csv", na = "")