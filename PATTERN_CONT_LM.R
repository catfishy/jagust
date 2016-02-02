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

# Import data
df_av45 = read.csv('dpgmm_alpha15.37_bilateral_AV45_ALL_longdata_continuous_slope.csv')
df_av45$CORTICAL_SUMMARY_POSITIVE = factor(df_av45$CORTICAL_SUMMARY_POSITIVE)
df_av45$RID = factor(df_av45$RID)
df_av45$group = factor(df_av45$group)
df_av45$diag_prior = factor(df_av45$diag_prior)
df_av45$APOE4_BIN = factor(df_av45$APOE4_BIN)
df_av45$Gender = factor(df_av45$Gender)

# Only keep N/SMC/EMCI/LMCI
#valid_diags = c('N','SMC','EMCI','LMCI')
#valid_diags = c('EMCI','LMCI')
valid_diags = c('N','SMC')
df_av45 = df_av45[which(df_av45$diag_prior %in% valid_diags),]

# only keep negatives
#df_av45 = df_av45[which(df_av45$CORTICAL_SUMMARY_POSITIVE == 0),]

# pattern weight models
fm_av45_nopattern = lm(AV45_slope ~ diag_prior + CORTICAL_SUMMARY_prior*APOE4_BIN +  I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., family=gaussian, df_av45)
fm_av45 = lm(AV45_slope ~ diag_prior + CORTICAL_SUMMARY_prior*APOE4_BIN +  I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + X1*APOE4_BIN + X4*APOE4_BIN + X22*APOE4_BIN + X25*APOE4_BIN + X9*APOE4_BIN + X0*APOE4_BIN + X34*APOE4_BIN + X2*APOE4_BIN + X19*APOE4_BIN + APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., family=gaussian, df_av45)
fm_av45_onlypatterns = lm(AV45_slope ~ diag_prior + X1*APOE4_BIN + X4*APOE4_BIN + X22*APOE4_BIN + X25*APOE4_BIN + X9*APOE4_BIN + X0*APOE4_BIN + X34*APOE4_BIN + X2*APOE4_BIN + X19*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., family=gaussian, df_av45)

# # pattern weight models with binary cs
# fm_av45_nopattern = lm(AV45_slope ~ diag_prior + CORTICAL_SUMMARY_POSITIVE*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., family=gaussian, df_av45)
# fm_av45 = lm(AV45_slope ~ diag_prior + CORTICAL_SUMMARY_POSITIVE*APOE4_BIN + X1*APOE4_BIN + X4*APOE4_BIN + X22*APOE4_BIN + X25*APOE4_BIN + X9*APOE4_BIN + X0*APOE4_BIN + X34*APOE4_BIN + X2*APOE4_BIN + X19*APOE4_BIN + APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., family=gaussian, df_av45)
# fm_av45_onlypatterns = lm(AV45_slope ~ diag_prior + X1*APOE4_BIN + X4*APOE4_BIN + X22*APOE4_BIN + X25*APOE4_BIN + X9*APOE4_BIN + X0*APOE4_BIN + X34*APOE4_BIN + X2*APOE4_BIN + X19*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., family=gaussian, df_av45)
# 

# summaries
fm_av45_nopattern_summary = summary(fm_av45_nopattern)
fm_av45_summary = summary(fm_av45)
fm_av45_onlypatterns_summary = summary(fm_av45_onlypatterns)

fm_av45_nopattern_summary$r.squared
fm_av45_summary$r.squared
fm_av45_onlypatterns_summary$r.squared

fm_av45_nopattern_summary$aic
fm_av45_summary$aic
fm_av45_onlypatterns_summary$aic

# anova
fm_av45_nopattern_anova = anova(fm_av45_nopattern,test='F')
fm_av45_anova = anova(fm_av45,test='F')
fm_av45_onlypatterns_anova = anova(fm_av45_onlypatterns,test='F')

# anova model comparisons
fm_modelcomparison_anova = anova(fm_av45_nopattern, fm_av45,test='LRT')

# PRINT MODEL OUTPUTS
sink('av45_lm_nopattern.txt'); print(fm_av45_nopattern_summary, correlation=TRUE); sink(file=NULL)
sink('av45_lm_withpattern.txt'); print(fm_av45_summary, correlation=TRUE); sink(file=NULL)
sink('av45_lm_onlypattern.txt'); print(fm_av45_onlypatterns_summary, correlation=TRUE); sink(file=NULL)

# PRINT MODEL ANOVA
sink('av45_lm_nopattern_anova.txt'); print(fm_av45_nopattern_anova, correlation=TRUE); sink(file=NULL)
sink('av45_lm_withpattern_anova.txt'); print(fm_av45_anova, correlation=TRUE); sink(file=NULL)
sink('av45_lm_onlypattern_anova.txt'); print(fm_av45_anova, correlation=TRUE); sink(file=NULL)
sink('av45_mc_anova.txt'); print(fm_modelcomparison_anova, correlation=TRUE); sink(file=NULL)

# plot fits
toplot = c(0,1,19,25,22)
for(g in toplot){
  jpeg(paste('fit_original_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],df_av45[,'AV45_slope'], col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('Original Data, Group:',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
  jpeg(paste('fit_withpattern_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45), col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('LM Predicted (with patterns)',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
  #jpeg(''); plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_nopattern), col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('LM Predicted (without patterns)',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
  jpeg(paste('fit_onlypattern_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_onlypatterns), col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('LM Predicted (only patterns)',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
}

# plot residuals
plot(fm_av45$fitted.values, fm_av45$residuals)
plot(fm_av45_nopattern$fitted.values, fm_av45_nopattern$residuals)
plot(fm_av45_onlypatterns$fitted.values, fm_av45_onlypatterns$residuals)


# get coeff names
rownames = row.names(fm_av45_onlypatterns_summary$coefficients)
for(i in rownames){
  print(i)
}

# contrasts
contrasts_all = read.csv('contrasts_lm.csv')
test = glht(fm_av45_onlypatterns,linfct=data.matrix(contrasts_all))
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
write.csv(contrasts, file = "av45_lm_contrasts.csv", na = "")