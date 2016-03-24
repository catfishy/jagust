library(lme4)
library(coefplot2)
library(ggplot2)
library(pbkrtest)
library(multcomp)
library(contrast)
library(xtable)
library(car)
library(stats)

library(gdata)
library(ggbiplot)

isPatternColumn = function(i){
  if (startsWith(i,'X')) return(TRUE) else return(FALSE)
}
isPatternColumn = Vectorize(isPatternColumn)

# Import data
df_av45 = read.csv('dpgmm_alpha13.13_bilateral_UW_MEM__ALL_slopedata.csv')
pattern_columns = Filter(isPatternColumn,names(df_av45))
df_patterns = df_av45[pattern_columns]

# demean pattern variables + binary variables
bin_variables = c('Gender','CORTICAL_SUMMARY_POSITIVE','APOE2_BIN','APOE4_BIN')
pattern_columns = c()
for (i in names(df_av45)){
  if (startsWith(i,'X') | i %in% bin_variables){
    df_av45[,eval(i)] = (df_av45[,eval(i)] * 2) - 1.0
  }
}

# Run PCA
df_patterns.pca = prcomp(df_patterns, center=TRUE, scale.s=TRUE)
df_patterns.pca_rotation = df_patterns.pca$rotation
df_patterns.transformed = predict(df_patterns.pca)
df_av45 = as.data.frame(cbind(as.matrix(df_av45),df_patterns.transformed))


# Convert to non factors to floats
to_factor = c('CORTICAL_SUMMARY_POSITIVE','RID','group','diag_prior','APOE4_BIN','APOE2_BIN','Gender')
for (i in names(df_av45)){
  if (!(i %in% to_factor)){
    df_av45[,eval(i)] = as.numeric(as.character(df_av45[,eval(i)]))
  }
}

# Only keep negatives
# df_av45 = df_av45[which(df_av45$CORTICAL_SUMMARY_POSITIVE == -1),]

# Only keep N/SMC/EMCI/LMCI
valid_diags = c('N','SMC','EMCI','LMCI')
#valid_diags = c('EMCI','LMCI')
#valid_diags = c('N','SMC')
df_av45 = df_av45[which(df_av45$diag_prior %in% valid_diags),]

# pattern weight models
fm_av45_onlycs = lm(UW_MEM__slope ~ CORTICAL_SUMMARY_prior + I(CORTICAL_SUMMARY_prior^2), df_av45)
fm_av45_nopattern = lm(UW_MEM__slope ~ diag_prior + CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs., df_av45)
fm_av45 = lm(UW_MEM__slope ~ diag_prior + 
               CORTICAL_SUMMARY_prior*APOE4_BIN +  
               I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + 
               X0*APOE4_BIN +
               X3*APOE4_BIN +
               X6*APOE4_BIN +
               X7*APOE4_BIN +
               X8*APOE4_BIN +
               X9*APOE4_BIN +
               X12*APOE4_BIN +
               X13*APOE4_BIN +
               X14*APOE4_BIN +
               X17*APOE4_BIN +
               X26*APOE4_BIN +
               X29*APOE4_BIN +
               Age.AV45 + Gender + Edu..Yrs., df_av45)
fm_av45_onlypatterns = lm(UW_MEM__slope ~ diag_prior + 
                            X0*APOE4_BIN +
                            X3*APOE4_BIN +
                            X6*APOE4_BIN +
                            X7*APOE4_BIN +
                            X8*APOE4_BIN +
                            X9*APOE4_BIN +
                            X12*APOE4_BIN +
                            X13*APOE4_BIN +
                            X14*APOE4_BIN +
                            X17*APOE4_BIN +
                            X26*APOE4_BIN +
                            X29*APOE4_BIN +
                            Age.AV45 + Gender + Edu..Yrs., df_av45)
fm_av45_onlypca = lm(UW_MEM__slope ~ diag_prior + 
                       PC1*APOE4_BIN +
                       PC2*APOE4_BIN +
                       PC3*APOE4_BIN +
                       PC4*APOE4_BIN +
                       PC5*APOE4_BIN +
                       PC6*APOE4_BIN +
                       PC7*APOE4_BIN +
                       PC8*APOE4_BIN +
                       PC9*APOE4_BIN +
                       PC10*APOE4_BIN +
                       PC11*APOE4_BIN +
                       Age.AV45 + Gender + Edu..Yrs., df_av45)
fm_av45_pca = lm(UW_MEM__slope ~ diag_prior + 
                   CORTICAL_SUMMARY_prior*APOE4_BIN +  
                   I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + 
                   PC1*APOE4_BIN +
                   PC2*APOE4_BIN +
                   PC3*APOE4_BIN +
                   PC4*APOE4_BIN +
                   PC5*APOE4_BIN +
                   PC6*APOE4_BIN +
                   PC7*APOE4_BIN +
                   PC8*APOE4_BIN +
                   PC9*APOE4_BIN +
                   PC10*APOE4_BIN +
                   PC11*APOE4_BIN +
                   Age.AV45 + Gender + Edu..Yrs., df_av45)
# summaries
fm_av45_nopattern_summary = summary(fm_av45_nopattern)
fm_av45_summary = summary(fm_av45)
fm_av45_onlypatterns_summary = summary(fm_av45_onlypatterns)
fm_av45_onlypca_summary = summary(fm_av45_onlypca)
fm_av45_pca_summary = summary(fm_av45_pca)

fm_av45_nopattern_summary$adj.r.squared
fm_av45_summary$adj.r.squared
fm_av45_onlypatterns_summary$adj.r.squared
fm_av45_onlypca_summary$adj.r.squared
fm_av45_pca_summary$adj.r.squared

fm_av45_nopattern_summary$aic
fm_av45_summary$aic
fm_av45_onlypatterns_summary$aic

# anova
fm_av45_nopattern_anova = Anova(fm_av45_nopattern,type='III')
fm_av45_anova = Anova(fm_av45,type='III')
fm_av45_onlypatterns_anova = Anova(fm_av45_onlypatterns,type='III')
fm_av45_onlypca_anova = Anova(fm_av45_onlypca,type='III')
fm_av45_pca_anova = Anova(fm_av45_pca,type='III')

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
toplot = c(42,19,8,7,1,0,12)
for(g in toplot){
  jpeg(paste('fit_original_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],df_av45[,'AV45_slope'], col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('Original Data, Group:',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
  jpeg(paste('fit_withpattern_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45), col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('LM Predicted (with patterns)',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
  jpeg(paste('fit_nopattern_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_nopattern), col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('LM Predicted (without patterns)',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
  jpeg(paste('fit_onlypattern_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_onlypatterns), col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('LM Predicted (only patterns)',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
}

# plot together
polfit = function(x) fm_av45_onlycs$coefficients[3]*x^2 + fm_av45_onlycs$coefficients[2]*x + fm_av45_onlycs$coefficients[1]
colors = rainbow(length(toplot))
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