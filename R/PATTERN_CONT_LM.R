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
library(gdata)
library(psych)
library(reshape2)

# # Calculate time to positivity threshold
# df_av45['threshold'] = threshold
# df_av45['threshold_diff'] = df_av45$threshold - df_av45$CORTICAL_SUMMARY_prior
# df_av45['to_threshold'] = df_av45$threshold_diff / df_av45$CORTICAL_SUMMARY_slope
# df_av45[which(df_av45$to_threshold < 0),'to_threshold'] = 0

# CONSTANTS
threshold = 1.27
pattern_prefix = 'NSFA_'
to_factor = c('RID','diag_prior','diag_post','APOE4_BIN','APOE2_BIN','Gender','Diag.AV45_long')
demog_columns = c('RID','APOE4_BIN','diag_prior','Age.AV45','Gender','Edu..Yrs.')
av45_columns = c('CORTICAL_SUMMARY_change','CORTICAL_SUMMARY_prior','CORTICAL_SUMMARY_post')
target = "CORTICAL_SUMMARY_change"

#valid_diags = c('N','SMC','EMCI','LMCI')
valid_diags = c('N','SMC')

#time_col_prefix = 'TIMEpostAV45_ADAS'
#value_col_prefix = 'ADAScog'
time_col_prefix = 'TIMEpostAV45_AVLT.'
value_col_prefix = 'AVLT.'

# FUNCTIONS
isPatternColumn = function(i){
  if (startsWith(i,pattern_prefix)) return(TRUE) else return(FALSE)
}
isPatternColumn = Vectorize(isPatternColumn)

lm.addvar = function(var.name) {
  paste('+',paste(var.name,'*','APOE4_BIN',sep=''))
}

lme.addvar = function(var.name) {
  paste('+',paste(var.name,'*','time',sep=''))
}

to.long = function(df, time_col_prefix, value_col_prefix) {
  # Keep relevant columns
  time_columns = Filter(function(i){startsWith(i,time_col_prefix)}, names(df))
  value_columns = Filter(function(i){startsWith(i,value_col_prefix)}, names(df))
  df = df[c(demog_columns,av45_columns,pattern_columns,time_columns,value_columns)]
  # Convert to long format
  df_time_wide = df[c(demog_columns,av45_columns,pattern_columns,time_columns)]
  colnames(df_time_wide) = gsub(time_col_prefix,'TP',names(df_time_wide))
  df_value_wide = df[c(demog_columns,av45_columns,pattern_columns,value_columns)]
  colnames(df_value_wide) = gsub(value_col_prefix,'TP',names(df_value_wide))
  df_time_long = melt(df_time_wide, 
                      id.vars=c(demog_columns,av45_columns,pattern_columns),
                      measure.vars=Filter(function(x){startsWith(x,'TP')},names(df_time_wide)),
                      variable.name='timepoint',
                      value.name='time')
  df_value_long = melt(df_value_wide,
                       id.vars=c(demog_columns,av45_columns,pattern_columns),
                       measure.vars=Filter(function(x){startsWith(x,'TP')},names(df_value_wide)),
                       variable.name='timepoint',
                       value.name='value')
  merge_on = c(demog_columns,av45_columns,pattern_columns,'timepoint')
  df_long = merge(df_time_long,df_value_long,merge_on)
  df_long[complete.cases(df_long[,names(df_long)]),]
}


# IMPORT
df_av45 = read.csv('pattern_dataset.csv')
df_av45 = df_av45[which(df_av45$diag_prior %in% valid_diags),]
for (i in names(df_av45)){
  if (i %in% to_factor){
    df_av45[,eval(i)] = as.factor(df_av45[,eval(i)])
  }
}
pattern_columns = Filter(isPatternColumn,names(df_av45))
df_long = to.long(df_av45, time_col_prefix, value_col_prefix)

# One by one pattern variable likelihood testing
lme.testvar = function(var.name) {
  base_str = paste('value',"~","CORTICAL_SUMMARY_prior*time + APOE4_BIN*time + diag_prior + Age.AV45 + Gender + Edu..Yrs.")
  add_str = lme.addvar(var.name)
  random_str = '+ (1 + time | RID)'
  form_base = as.formula(paste(base_str,random_str))
  form = as.formula(paste(base_str,add_str,random_str))
  fm = lmer(form,data=df_long)
  fm_base = lmer(form_base,df_long)
  like = anova(fm_base,fm)
  like.p = like$`Pr(>Chisq)`[2]
}
lm.testvar = function(var.name) {
  base_str = paste(target,"~","CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + diag_prior + Age.AV45 + Gender + Edu..Yrs.")
  add_str = lm.addvar(var.name)
  form_base = as.formula(base_str)
  form = as.formula(paste(base_str,add_str))
  fm = lm(form,data=df_av45)
  fm_base = lm(form_base,data=df_av45)
  like = anova(fm_base,fm)
  like.p = like$`Pr(>F)`[2]
}

# LM (Fixed Effects)
like.pvalues = lapply(pattern_columns,lm.testvar)
valid_patterns = pattern_columns[like.pvalues <= 0.05]
form.addons = lapply(valid_patterns,lm.addvar)
base_form = paste(target,"~ diag_prior + Age.AV45 + Gender + Edu..Yrs.")
nopattern_form = paste(target,"~ diag_prior + CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
onlypattern_form = paste(base_form,paste(form.addons,collapse=' '))
full_form = paste(nopattern_form,paste(form.addons,collapse=' '))

fm_base = lm(as.formula(base_form),df_av45)
fm_nopattern = lm(as.formula(nopattern_form),df_av45)
fm_onlypattern = lm(as.formula(onlypattern_form),df_av45)
fm_full = lm(as.formula(full_form),df_av45)

fm_base.summary = summary(fm_base)
fm_nopattern.summary = summary(fm_nopattern)
fm_onlypattern.summary = summary(fm_onlypattern)
fm_full.summary = summary(fm_full)

fm_base.summary$adj.r.squared
fm_nopattern.summary$adj.r.squared
fm_onlypattern.summary$adj.r.squared
fm_full.summary$adj.r.squared

fm_base.anova = Anova(fm_base,type='III')
fm_nopattern.anova = Anova(fm_nopattern,type='III')
fm_onlypattern.anova = Anova(fm_onlypattern,type='III')
fm_full.anova = Anova(fm_full,type='III')

# LME (Mixed Effects)
like.pvalues = lapply(pattern_columns,lme.testvar)
valid_patterns = pattern_columns[like.pvalues <= 0.05]
form.addons = lapply(valid_patterns,lme.addvar)
random_str = '+ (1 + time | RID)'
base_str = "value ~ time + diag_prior + Age.AV45 + Gender + Edu..Yrs."
nopattern_str = "value ~ diag_prior + CORTICAL_SUMMARY_prior*time + APOE4_BIN*time + Age.AV45 + Gender + Edu..Yrs."
base_form = paste(base_str,random_str)
nopattern_form = paste(nopattern_str,random_str)
onlypattern_form = paste(base_str,paste(form.addons,collapse=' '),random_str)
full_form = paste(nopattern_str,paste(form.addons,collapse=' '),random_str)

fm_base = lmer(as.formula(base_form),df_long)
fm_nopattern = lmer(as.formula(nopattern_form),df_long)
fm_onlypattern = lmer(as.formula(onlypattern_form),df_long)
fm_full = lmer(as.formula(full_form),df_long)

fm_base.summary = summary(fm_base)
fm_nopattern.summary = summary(fm_nopattern)
fm_onlypattern.summary = summary(fm_onlypattern)
fm_full.summary = summary(fm_full)

fm_base.summary$AICtab
fm_nopattern.summary$AICtab
fm_onlypattern.summary$AICtab
fm_full.summary$AICtab

fm_base.anova = Anova(fm_base,type='III')
fm_nopattern.anova = Anova(fm_nopattern,type='III')
fm_onlypattern.anova = Anova(fm_onlypattern,type='III')
fm_full.anova = Anova(fm_full,type='III')




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