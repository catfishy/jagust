lm.testvar = function(var.name) {
  # CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + 
  base_str = paste(target,"~","diag_prior*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
  add_str = lm.addvar(var.name)
  form_base = as.formula(base_str)
  form = as.formula(paste(base_str,add_str))
  fm = lm(form,data=df_av45)
  fm_base = lm(form_base,data=df_av45)
  like = anova(fm_base,fm)
  like.p = like$`Pr(>F)`[2]
}

# PLOTTING
df_av45_apoepos = df_av45[df_av45$APOE4_BIN == 1,]
df_av45_apoeneg = df_av45[df_av45$APOE4_BIN == 0,]

qplot(df_av45_apoepos$CORTICAL_SUMMARY_prior, 
      df_av45_apoepos$CORTICAL_SUMMARY_change, 
      data=df_av45_apoepos, 
      colour=df_av45_apoepos$NSFA_16) + 
  scale_colour_gradient(low="black", high="pink", limits=c(0,2), na.value='black')

qplot(df_av45_apoeneg$CORTICAL_SUMMARY_prior, 
      df_av45_apoeneg$CORTICAL_SUMMARY_change, 
      data=df_av45_apoeneg, 
      colour=df_av45_apoeneg$NSFA_16) + 
  scale_colour_gradient(low="black", high="pink", limits=c(0,2), na.value='black')

qplot(df_av45$CORTICAL_SUMMARY_prior, 
      df_av45$UW_EF_BL_3months, 
      data=df_av45, 
      colour=df_av45$NSFA_24) + 
  scale_colour_gradient(low="black", high="pink", limits=c(0,2),na.value='black')

qplot(fitted(fm_onlypattern), resid(fm_onlypattern)) + geom_hline(yintercept=0)
qplot(df_av45[,target],fitted(fm_onlypattern))
avPlots(fm_onlypattern, id.n=2, id.cex=0.7)



# PRINTING OUTPUTS AND GRAPHING

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
toplot = c(0,1,2,3,14,15)
# plot together
polfit = function(x) fm_av45_onlycs$coefficients[3]*x^2 + fm_av45_onlycs$coefficients[2]*x + fm_av45_onlycs$coefficients[1]
colors = rainbow(length(toplot))
plot(df_av45[,'CORTICAL_SUMMARY_prior'],df_av45[,'CORTICAL_SUMMARY_slope'], pch=4, cex=1, lwd=0.6, main='Significant Pattern Groups', xlab='Baseline Florbetapir Cortical Summary SUVR', ylab='Cortical Summary SUVR Annualized Change')
for(i in 1:length(toplot)){
  g = toplot[i]
  c = colors[i]
  x = df_av45[which(df_av45$group==g),'CORTICAL_SUMMARY_prior']
  y = df_av45[which(df_av45$group==g),'CORTICAL_SUMMARY_slope']
  points(x,y, bg=c, pch=21, cex=1.1, lwd=1)
}
curve(polfit,add=T)
legend('topright', legend=sapply(toplot, function(x) paste('Group #',x,sep='')), fill=colors)


for(g in toplot){
  jpeg(paste('fit_original_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],df_av45[,'CORTICAL_SUMMARY_slope'], col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('Original Data, Group:',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
  jpeg(paste('fit_withpattern_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45), col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('LM Predicted (with patterns)',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
  jpeg(paste('fit_nopattern_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_nopattern), col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('LM Predicted (without patterns)',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
  jpeg(paste('fit_onlypattern_group',g,'.jpeg',sep='')); plot(df_av45[,'CORTICAL_SUMMARY_prior'],predict(fm_av45_onlypatterns), col=ifelse(df_av45[,'group']==g, "red", "black"), main=paste('LM Predicted (only patterns)',g), xlab='Cortical Summary BL', ylab='Annualized AV45 Slope'); dev.off();
}


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