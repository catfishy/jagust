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
library(piecewiseSEM)
library(LambertW)
library(nnet)
library(DAAG)
library(relaimpo)
library(caret)
library(cvTools)
library(VGAM)
library(lmtest)
library(languageR)
library(stringr)
library(covTest)

source('R/LM_FUNCS.R')

# CONSTANTS
pattern_prefix = 'NSFA_'
to_factor = c('RID','ad_prior','ad_post','positive_prior','positive_post',
              'diag_prior','diag_post','APOE4_BIN','APOE2_BIN','Gender',
              'Diag.AV45','Diag.AV1451','positive_prior','positive_post',
              'AV45_NONTP_wcereb_BIN1.11')
to_standardize = c('CORTICAL_SUMMARY_prior','Age.AV45','Edu..Yrs.')
demog_columns = c('RID','APOE4_BIN','Diag.AV45','Age.AV45','Gender','Edu..Yrs.')
av45_columns = c('CORTICAL_SUMMARY_prior')

#target = "CORTICAL_SUMMARY_change"
#target = "UW_EF_AV45_1"
#target = "UW_EF_slope"
#target = "ADAS_AV45_1"
#target = "ADASslope_postAV45"
#target = "AVLT_AV45_1"
#target = "AVLT_slope_postAV45"
#target = "UW_MEM_AV45_1"
#target = "UW_MEM_slope"
#target = "CSF_ABETA_closest_AV45_1"
#target = "CSF_TAU_closest_AV45_1"
#target = "CSF_PTAU_closest_AV45_1"
target = 'MMSE_AV45_1'

#output_folder = 'R/output/'
output_folder = 'R/output_all_diag/'
output_folder = 'R/output_neg_emci/'

positive_value=1
valid_diags = c('N','SMC','EMCI','LMCI','AD')
#valid_diags = c('N')
#valid_diags = c('EMCI')
#valid_diags = c('LMCI')


#time_col_prefix = 'TIMEpostAV45_ADAS'
#value_col_prefix = 'ADAScog'
#time_col_prefix = 'TIMEpostAV45_AVLT.'
#value_col_prefix = 'AVLT.'
#time_col_prefix = 'WMH_postAV45.'
#value_col_prefix = 'WMH_percentOfICV.'
time_col_prefix = 'UW_EF_postAV45_'
value_col_prefix = 'UW_EF_'

# IMPORT
df_av45 = read.csv('nsfa/av45_pattern_dataset.csv')
non.na = complete.cases(df_av45[,c(demog_columns,av45_columns,target)])
df_av45 = df_av45[non.na,]

# filter by diag or postivity
df_av45 = df_av45[which(df_av45$Diag.AV45 %in% valid_diags),]
df_av45 = df_av45[which(df_av45[,'AV45_NONTP_wcereb_BIN1.11'] == positive_value),]

# make factors
for (i in names(df_av45)){
  if (i %in% to_factor){
    df_av45[,eval(i)] = as.factor(as.character(df_av45[,eval(i)]))
  }
}
pattern_columns = Filter(isPatternColumn,names(df_av45))
naive_columns = Filter(isNaiveColumn,names(df_av45))




# # look at histograms
# for (pcol in pattern_columns) {
#   p = ggplot(df_av45, aes_string(pcol, colour = 'positive_prior')) + geom_histogram(binwidth=0.1)
#   print(p)
# }


# standardize predictors
cross_to_standardize = c(to_standardize,pattern_columns,naive_columns,target)
cross_normalization = preProcess(df_av45[,cross_to_standardize])
df_av45[,cross_to_standardize] = predict(cross_normalization, df_av45[,cross_to_standardize])
# long_to_standardize = c(to_standardize,pattern_columns,'time','value')
# long_normalization = preProcess(df_long[,long_to_standardize])
# df_long[,long_to_standardize] = predict(long_normalization, df_long[,long_to_standardize])

# # Create long dataset
# df_long = to.long(df_av45, time_col_prefix, value_col_prefix)
# df_long$value = Gaussianize(df_long$value, type='hh', method='MLE', return.u=TRUE)

# make response normal
#df_av45[,eval(target)] = Gaussianize(df_av45[,eval(target)], type='hh', method='MLE', return.u=TRUE)

#all.addons = lapply(pattern_columns,lm.addvar)
#naive.addons = lapply(naive_columns,lm.addvar)
all.addons = paste('+',paste(pattern_columns,collapse=' + '))
naive.addons = paste('+',paste(naive_columns,collapse=' + '))
addons_form = str_replace(paste(target,"~",paste(all.addons,collapse=' ')),"\\+ ","")

diag_counts = table(df_av45$Diag.AV45)
if (length((diag_counts[diag_counts>0])) > 1) {
  diag_str = 'Diag.AV45 + '
} else {
  diag_str = ''
}

base_form = paste(target,"~",diag_str,"Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN")
nopattern_form = paste(target,"~",diag_str,"CORTICAL_SUMMARY_prior + I(CORTICAL_SUMMARY_prior^2) + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN")
onlypattern_form = paste(base_form,paste(all.addons,collapse=' '))
pattern_form = paste(target,'~',paste(pattern_columns,collapse=' + '))
full_form = paste(nopattern_form,paste(all.addons,collapse=' '))
naive_form = paste(base_form,paste(naive.addons,collapse=' '))

# LARS lasso

nopattern_x = getxy(nopattern_form,df_av45)
y = as.numeric(df_av45[,target])
nopattern.lars.model = lars(nopattern_x,y,type='lasso')
nopattern.lars.test = covTest(nopattern.lars.model,nopattern_x,y)$results
nopattern.lars.sigcoef.idx = nopattern.lars.test[nopattern.lars.test[,'P-value'] < 0.2 & !is.na(nopattern.lars.test[,'P-value']),'Predictor_Number']
nopattern.lars.coef = coef(nopattern.lars.model, s=which.min(summary(nopattern.lars.model)$Cp), mode='step')
nopattern.lars.r2 = nopattern.lars.model$R2[which.min(summary(nopattern.lars.model)$Cp)]
nopattern.lars.sigcoef = nopattern.lars.coef[nopattern.lars.sigcoef.idx]

onlypattern_x = getxy(onlypattern_form,df_av45)
y = as.numeric(df_av45[,target])
onlypattern.lars.model = lars(onlypattern_x,y,type='lasso')
onlypattern.lars.test = covTest(onlypattern.lars.model,onlypattern_x,y)$results
onlypattern.lars.sigcoef.idx = onlypattern.lars.test[onlypattern.lars.test[,'P-value'] < 0.2 & !is.na(onlypattern.lars.test[,'P-value']),'Predictor_Number']
onlypattern.lars.coef = coef(onlypattern.lars.model, s=which.min(summary(onlypattern.lars.model)$Cp), mode='step')
onlypattern.lars.r2 = onlypattern.lars.model$R2[which.min(summary(onlypattern.lars.model)$Cp)]
onlypattern.lars.sigcoef = onlypattern.lars.coef[onlypattern.lars.sigcoef.idx]

# LASSO Penalized regression
nopattern.lasso.model = run.lasso(nopattern_form,df_av45,'Rsquared')
nopattern.lasso.metric = subset(nopattern.lasso.model$results, fraction == nopattern.lasso.model$bestTune$fraction)
nopattern.lasso.coef = predict.enet(nopattern.lasso.model$finalModel, type='coefficients',s=nopattern.lasso.model$bestTune$fraction, mode='fraction')$coefficients
nopattern.lasso.coef = nopattern.lasso.coef[nopattern.lasso.coef != 0]

onlypattern.lasso.model = run.lasso(onlypattern_form,df_av45,'Rsquared')
onlypattern.lasso.metric = subset(onlypattern.lasso.model$results, fraction == onlypattern.lasso.model$bestTune$fraction)
onlypattern.lasso.coef = predict.enet(onlypattern.lasso.model$finalModel, type='coefficients',s=onlypattern.lasso.model$bestTune$fraction, mode='fraction')$coefficients
onlypattern.lasso.coef = onlypattern.lasso.coef[onlypattern.lasso.coef != 0]

naive.lasso.model = run.lasso(naive_form,df_av45,'RMSE')
naive.lasso.metric = subset(naive.lasso.model$results, fraction == naive.lasso.model$bestTune$fraction)
naive.lasso.coef = predict.enet(naive.lasso.model$finalModel, type='coefficients',s=naive.lasso.model$bestTune$fraction, mode='fraction')$coefficients
naive.lasso.coef = naive.lasso.coef[naive.lasso.coef != 0]

# GLMNET Penalized Regression
nopattern.glmnet.model = run.glmnet(nopattern_form,df_av45,'RMSE')
nopattern.glmnet.metric = subset(nopattern.glmnet.model$results, alpha == nopattern.glmnet.model$bestTune$alpha & lambda == nopattern.glmnet.model$bestTune$lambda)
nopattern.glmnet.coef = predict.glmnet(nopattern.glmnet.model$finalModel,type='coefficients',s=nopattern.glmnet.model$bestTune$lambda)

onlypattern.glmnet.model = run.glmnet(onlypattern_form,df_av45,'RMSE')
onlypattern.glmnet.metric = subset(onlypattern.glmnet.model$results, alpha == onlypattern.glmnet.model$bestTune$alpha & lambda == onlypattern.glmnet.model$bestTune$lambda)
onlypattern.glmnet.coef = predict.glmnet(onlypattern.glmnet.model$finalModel,type='coefficients',s=onlypattern.glmnet.model$bestTune$lambda)

naive.glmnet.model = run.glmnet(naive_form,df_av45,'RMSE')
naive.glmnet.metric = subset(naive.glmnet.model$results, alpha == naive.glmnet.model$bestTune$alpha & lambda == naive.glmnet.model$bestTune$lambda)
naive.glmnet.coef = predict.glmnet(naive.glmnet.model$finalModel,type='coefficients',s=naive.glmnet.model$bestTune$lambda)




# LM RFE

rfe.base = run.rfe(base_form, target, df_av45)
rfe.nopattern = run.rfe(nopattern_form, target, df_av45)
#rfe.onlypattern = run.rfe(onlypattern_form, target, df_av45)
rfe.full = run.rfe(full_form, target, df_av45)

fm_base = rfe.base$fit
fm_nopattern = rfe.nopattern$fit
#fm_onlypattern = rfe.onlypattern$fit
fm_full = rfe.full$fit

fm_base.summary = summary(fm_base)
fm_nopattern.summary = summary(fm_nopattern)
#fm_onlypattern.summary = summary(fm_onlypattern)
fm_full.summary = summary(fm_full)

fm_base.fit = sem.model.fits(fm_base)
fm_nopattern.fit = sem.model.fits(fm_nopattern)
#fm_onlypattern.fit = sem.model.fits(fm_onlypattern)
fm_full.fit = sem.model.fits(fm_full)

fm_base.anova = Anova(fm_base,type='III')
fm_nopattern.anova = Anova(fm_nopattern,type='III')
#fm_onlypattern.anova = Anova(fm_onlypattern,type='III')
fm_full.anova = Anova(fm_full,type='III')


save.printout(paste(output_folder,target,'_fm_base_summary','.txt',sep=''),fm_base.summary)
save.printout(paste(output_folder,target,'_fm_nopattern_summary','.txt',sep=''),fm_nopattern.summary)
save.printout(paste(output_folder,target,'_fm_full_summary','.txt',sep=''),fm_full.summary)
save.printout(paste(output_folder,target,'_fm_base_fit','.txt',sep=''),fm_base.fit)
save.printout(paste(output_folder,target,'_fm_nopattern_fit','.txt',sep=''),fm_nopattern.fit)
save.printout(paste(output_folder,target,'_fm_full_fit','.txt',sep=''),fm_full.fit)
save.printout(paste(output_folder,target,'_fm_base_anova','.txt',sep=''),fm_base.anova)
save.printout(paste(output_folder,target,'_fm_nopattern_anova','.txt',sep=''),fm_nopattern.anova)
save.printout(paste(output_folder,target,'_fm_full_anova','.txt',sep=''),fm_full.anova)
save(rfe.base,file=paste(output_folder,target,'_rfe_base_obj',sep=''))
save(rfe.nopattern,file=paste(output_folder,target,'_rfe_nopattern_obj',sep=''))
save(rfe.full,file=paste(output_folder,target,'_rfe_full_obj',sep=''))

fm_base.plotfn = function() {par(mfrow=c(2,2));plot(fm_base);title("Base Model", outer=T, line=-2);}
fm_nopattern.plotfn = function() {par(mfrow=c(2,2));plot(fm_nopattern);title("No Pattern Model", outer=T, line=-2);}
fm_full.plotfn = function() {par(mfrow=c(2,2));plot(fm_full);title("Full Model", outer=T, line=-2);}
save.plot(paste(output_folder,target,'_fm_base_lmplot.pdf',sep=''), fm_base.plotfn)
save.plot(paste(output_folder,target,'_fm_nopattern_lmplot.pdf',sep=''), fm_nopattern.plotfn)
save.plot(paste(output_folder,target,'_fm_full_lmplot.pdf',sep=''), fm_full.plotfn)

save.plot(paste(output_folder,target,'_fm_base_avplot.pdf',sep=''), function() {avPlots(fm_base, ask=FALSE)})
save.plot(paste(output_folder,target,'_fm_nopattern_avplot.pdf',sep=''), function() {avPlots(fm_nopattern, ask=FALSE)})
save.plot(paste(output_folder,target,'_fm_full_avplot.pdf',sep=''), function() {avPlots(fm_full, ask=FALSE)})

fm_base.summary
fm_nopattern.summary
fm_full.summary

fm_base.fit
fm_nopattern.fit
fm_full.fit


