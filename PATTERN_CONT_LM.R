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

# only keep negatives
#df_av45 = df_av45[which(df_av45$CORTICAL_SUMMARY_POSITIVE == 0),]

# only keep patterns
valid_patterns = c(1,4,22,25,9,0,34,2,19)
#df_av45 = df_av45[which(df_av45$group %in% valid_patterns),]

# pattern weight models
fm_av45_onlycs = lm(AV45_slope ~ CORTICAL_SUMMARY_prior + I(CORTICAL_SUMMARY_prior^2), family=gaussian, df_av45)
fm_av45_nopattern = lm(AV45_slope ~ diag_prior + CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., family=gaussian, df_av45)
fm_av45 = lm(AV45_slope ~ diag_prior + CORTICAL_SUMMARY_prior*APOE4_BIN +  I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + X1*APOE4_BIN + X4*APOE4_BIN + X22*APOE4_BIN + X25*APOE4_BIN + X9*APOE4_BIN + X0*APOE4_BIN + X34*APOE4_BIN + X2*APOE4_BIN + X19*APOE4_BIN + APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., family=gaussian, df_av45)
fm_av45_onlypatterns = lm(AV45_slope ~ diag_prior + X1*APOE4_BIN + X4*APOE4_BIN + X22*APOE4_BIN + X25*APOE4_BIN + X9*APOE4_BIN + X0*APOE4_BIN + X34*APOE4_BIN + X2*APOE4_BIN + X19*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., family=gaussian, df_av45)

# # pattern weight models with all patterns
# fm_av45_nopattern = lm(AV45_slope ~ diag_prior + CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., family=gaussian, df_av45)
# fm_av45 = lm(AV45_slope ~ diag_prior + CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + X5*APOE4_BIN + X3*APOE4_BIN + X6*APOE4_BIN + X29*APOE4_BIN + X56*APOE4_BIN + X1*APOE4_BIN + X4*APOE4_BIN + X22*APOE4_BIN + X25*APOE4_BIN + X9*APOE4_BIN + X0*APOE4_BIN + X34*APOE4_BIN + X2*APOE4_BIN + X19*APOE4_BIN + APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., family=gaussian, df_av45)
# fm_av45_onlypatterns = lm(AV45_slope ~ diag_prior + X5*APOE4_BIN + X3*APOE4_BIN + X6*APOE4_BIN + X29*APOE4_BIN + X56*APOE4_BIN + X1*APOE4_BIN + X4*APOE4_BIN + X22*APOE4_BIN + X25*APOE4_BIN + X9*APOE4_BIN + X0*APOE4_BIN + X34*APOE4_BIN + X2*APOE4_BIN + X19*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., family=gaussian, df_av45)

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
fm_av45_nopattern_anova = Anova(fm_av45_nopattern,type='III')
fm_av45_anova = Anova(fm_av45,type='III')
fm_av45_onlypatterns_anova = Anova(fm_av45_onlypatterns,type='III')
# fm_av45_nopattern_anova = Anova(fm_av45_nopattern,type='III')
# fm_av45_anova = Anova(fm_av45,type='III')
# fm_av45_onlypatterns_anova = Anova(fm_av45_onlypatterns,,type='III')

# anova model comparisons
fm_modelcomparison_anova = anova(fm_av45_nopattern, fm_av45, test='LRT')

# PRINT MODEL OUTPUTS
sink('av45_lm_nopattern.txt'); print(fm_av45_nopattern_summary, correlation=TRUE); sink(file=NULL)
sink('av45_lm_withpattern.txt'); print(fm_av45_summary, correlation=TRUE); sink(file=NULL)
sink('av45_lm_onlypattern.txt'); print(fm_av45_onlypatterns_summary, correlation=TRUE); sink(file=NULL)

# PRINT MODEL ANOVA
sink('av45_lm_nopattern_anova.txt'); print(fm_av45_nopattern_anova, correlation=TRUE); sink(file=NULL)
sink('av45_lm_withpattern_anova.txt'); print(fm_av45_anova, correlation=TRUE); sink(file=NULL)
sink('av45_lm_onlypattern_anova.txt'); print(fm_av45_onlypatterns_anova, correlation=TRUE); sink(file=NULL)
sink('av45_mc_anova.txt'); print(fm_modelcomparison_anova, correlation=TRUE); sink(file=NULL)

# plot fits
toplot = c(0,1,25,22)
labels =
for(g in toplot){
  jpeg(paste('fit_original_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],df_av45[,'AV45_slope'], col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('Original Data, Group:',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
  jpeg(paste('fit_withpattern_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45), col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('LM Predicted (with patterns)',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
  jpeg(paste('fit_nopattern_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_nopattern), col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('LM Predicted (without patterns)',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
  jpeg(paste('fit_onlypattern_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_onlypatterns), col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('LM Predicted (only patterns)',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
}


polfit = function(x) fm_av45_onlycs$coefficients[3]*x^2 + fm_av45_onlycs$coefficients[2]*x + fm_av45_onlycs$coefficients[1]
colors = rainbow(length(toplot))
#colors = topo.colors(length(toplot))

#pdf('fit_original_siggroups.pdf',height=11,width=10)
plot(df_av45[,'CORTICAL_SUMMARY_prior'],df_av45[,'AV45_slope'], pch=4, cex=1, lwd=0.6, main='Significant Pattern Groups', xlab='Baseline Florbetapir Cortical Summary SUVR', ylab='Cortical Summary SUVR Annualized Change')
for(i in 1:length(toplot)){
  g = toplot[i]
  c = colors[i]
  x = df_av45[which(df_av45$group==g),'CORTICAL_SUMMARY_prior']
  y = df_av45[which(df_av45$group==g),'AV45_slope']
  points(x,y, bg=c, pch=21, cex=1.1, lwd=1)
}
curve(polfit,add=T)
legend('topright', legend=sapply(toplot, function(x) paste('Group #',x,sep='')), fill=colors)
#dev.off()

# plot fits
plot(df_av45[,'CORTICAL_SUMMARY_prior'],df_av45[,'AV45_slope'],bg='yellow',pch=21)
points(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45),bg='red',pch=21)
points(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_nopattern),bg='blue',pch=21)
points(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_onlypatterns),bg='green',pch=21)


# plot residuals
plot(1,type-'n')
points(fm_av45_nopattern$fitted.values, fm_av45_nopattern$residuals,bg='blue',pch=21)
points(fm_av45$fitted.values, fm_av45$residuals, bg='red', pch=21)



# get coeff names
rownames = row.names(fm_av45_summary$coefficients)
for(i in rownames){
  print(i)
}

# contrasts
contrasts_all = read.csv('contrasts_lm.csv')
test = glht(fm_av45,linfct=data.matrix(contrasts_all))
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