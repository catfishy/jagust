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
library(scales)

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

# target = "CORTICAL_SUMMARY_change"
# target = "UW_EF_AV45_1"
#target = "UW_EF_slope"
# target = "ADAS_AV45_1"
target = "ADASslope_postAV45"
# target = "AVLT_AV45_1"
#target = "AVLT_slope_postAV45"
# target = "UW_MEM_AV45_1"
# target = "UW_MEM_slope"
# target = "CSF_ABETA_closest_AV45_1"
# target = "UCB_FS_HC.ICV_AV45_1"
# target = "CSF_TAU_closest_AV45_1"
#target = "CSF_PTAU_closest_AV45_1"
# target = 'MMSE_AV45_1'

#output_folder = 'R/output/'
output_folder = 'R/output_all_diag/'
output_folder = 'R/output_neg_emci/'

positive_value=1
valid_diags = c('N','SMC','EMCI','LMCI','AD')
# valid_diags = c('N','SMC')
# valid_diags = c('EMCI')
# valid_diags = c('LMCI')
# valid_diags = c('AD')
# valid_diags = c('EMCI','LMCI','AD')
# valid_diags = c('LMCI','AD')

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

# remove target outliers
# target.mean = mean(df_av45[,target])
# target.sd = sd(df_av45[,target])
# df_av45 = df_av45[df_av45[,target] <= target.mean+target.sd*5,]
# df_av45 = df_av45[df_av45[,target] >= target.mean-target.sd*5,]

# filter by diag or positivity
df_av45 = df_av45[which(df_av45$Diag.AV45 %in% valid_diags),]
df_av45 = df_av45[which(df_av45[,'positive_prior'] == positive_value),]

# # filter by percentage around cutoff
# cutoff = 0.87
# keep = 300
# df_av45$cutoff_diff = abs(df_av45$CORTICAL_SUMMARY_prior-cutoff)
# df_av45 = df_av45[order(df_av45$cutoff_diff, decreasing=FALSE)[1:keep],]


# make factors
for (i in names(df_av45)){
  if (i %in% to_factor){
    df_av45[,eval(i)] = as.factor(as.character(df_av45[,eval(i)]))
  }
}
df_av45$Diag.AV45 = factor(df_av45$Diag.AV45, levels=valid_diags)
pattern_columns = Filter(isPatternColumn,names(df_av45))
naive_columns = Filter(isNaiveColumn,names(df_av45))




# # standardize predictors
# cross_to_standardize = c(to_standardize,pattern_columns,naive_columns,target)
# cross_normalization = preProcess(df_av45[,cross_to_standardize])
# df_av45[,cross_to_standardize] = predict(cross_normalization, df_av45[,cross_to_standardize])



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
# if (length((diag_counts[diag_counts>0])) > 1) {
#   diag_str = 'Diag.AV45 + '
# } else {
#   diag_str = ''
# }
diag_str = ''

base_form = paste(target,"~",diag_str,"Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN")
nopattern_form = paste(target,"~",diag_str,"CORTICAL_SUMMARY_prior + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN")
# I(CORTICAL_SUMMARY_prior^2) + 
onlypattern_form = paste(base_form,paste(all.addons,collapse=' '))
pattern_form = paste(target,'~',paste(pattern_columns,collapse=' + '))
full_form = paste(nopattern_form,paste(all.addons,collapse=' '))
naive_form = paste(base_form,paste(naive.addons,collapse=' '))



# LARS lasso
nopattern_x = getxy(nopattern_form,df_av45)
y = as.numeric(df_av45[,target])
nopattern.lars.model = lars(nopattern_x,y,type='lasso')
nopattern.lars.test = covTest(nopattern.lars.model,nopattern_x,y)$results
nopattern.lars.test = data.frame(nopattern.lars.test[complete.cases(nopattern.lars.test),])
nopattern.lars.coef = coef(nopattern.lars.model, s=which.min(summary(nopattern.lars.model)$Cp), mode='step')
nopattern.lars.test$name = names(nopattern.lars.coef[nopattern.lars.test[,'Predictor_Number']])
nopattern.lars.test$coef = nopattern.lars.coef[nopattern.lars.test[,'Predictor_Number']]
nopattern.lars.sig = nopattern.lars.test[nopattern.lars.test$P.value <= 0.05,]
nopattern.lars.cp = min(summary(nopattern.lars.model)$Cp)
nopattern.lars.r2 = nopattern.lars.model$R2[which.min(summary(nopattern.lars.model)$Cp)]
nopattern.lars.n = attr(summary(nopattern.lars.model)$Cp,'n')
nopattern.lars.p = NROW(nopattern.lars.coef[nopattern.lars.coef != 0])
nopattern.lars.r2adj = r2adj(nopattern.lars.r2,nopattern.lars.n,nopattern.lars.p)
nopattern.lars.nonzero = nopattern.lars.test[nopattern.lars.test$coef != 0,'name']
paste(nopattern.lars.nonzero,collapse=' + ')
paste(target,'~',paste(nopattern.lars.sig[,'name'], collapse=' + '))

onlypattern_x = getxy(onlypattern_form,df_av45)
y = as.numeric(df_av45[,target])
onlypattern.lars.model = lars(onlypattern_x,y,type='lasso')
onlypattern.lars.test = covTest(onlypattern.lars.model,onlypattern_x,y)$results
onlypattern.lars.test = data.frame(onlypattern.lars.test[complete.cases(onlypattern.lars.test),])
onlypattern.lars.coef = coef(onlypattern.lars.model, s=which.min(summary(onlypattern.lars.model)$Cp), mode='step')
onlypattern.lars.test$name = names(onlypattern.lars.coef[onlypattern.lars.test[,'Predictor_Number']])
onlypattern.lars.test$coef = onlypattern.lars.coef[onlypattern.lars.test[,'Predictor_Number']]
onlypattern.lars.sig = onlypattern.lars.test[onlypattern.lars.test$P.value <= 0.05,]
onlypattern.lars.cp = min(summary(onlypattern.lars.model)$Cp)
onlypattern.lars.r2 = onlypattern.lars.model$R2[which.min(summary(onlypattern.lars.model)$Cp)]
onlypattern.lars.n = attr(summary(onlypattern.lars.model)$Cp,'n')
onlypattern.lars.p = NROW(onlypattern.lars.coef[onlypattern.lars.coef != 0])
onlypattern.lars.r2adj = r2adj(onlypattern.lars.r2,onlypattern.lars.n,onlypattern.lars.p)
onlypattern.lars.nonzero = onlypattern.lars.test[onlypattern.lars.test$coef != 0,'name']
paste(onlypattern.lars.nonzero,collapse=' + ')
paste(target,'~',paste(onlypattern.lars.sig[,'name'], collapse=' + '))

full_x = getxy(full_form,df_av45)
y = as.numeric(df_av45[,target])
full.lars.model = lars(full_x,y,type='lasso')
full.lars.test = covTest(full.lars.model,full_x,y)$results
full.lars.test = data.frame(full.lars.test[complete.cases(full.lars.test),])
full.lars.coef = coef(full.lars.model, s=which.min(summary(full.lars.model)$Cp), mode='step')
full.lars.test$name = names(full.lars.coef[full.lars.test[,'Predictor_Number']])
full.lars.test$coef = full.lars.coef[full.lars.test[,'Predictor_Number']]
full.lars.sig = full.lars.test[full.lars.test$P.value <= 0.05,]
full.lars.cp = min(summary(full.lars.model)$Cp)
full.lars.r2 = full.lars.model$R2[which.min(summary(full.lars.model)$Cp)]
full.lars.n = attr(summary(full.lars.model)$Cp,'n')
full.lars.p = NROW(full.lars.coef[full.lars.coef != 0])
full.lars.r2adj = r2adj(full.lars.r2,full.lars.n,full.lars.p)
full.lars.nonzero = full.lars.test[full.lars.test$coef != 0,'name']
paste(full.lars.nonzero,collapse=' + ')
paste(target,'~',paste(full.lars.sig[,'name'], collapse=' + '))


# r2 shrinkage
# test1.form = 'UW_MEM_AV45_1 ~ CORTICAL_SUMMARY_prior + Gender + Edu..Yrs.'
# test2.form = 'UW_MEM_AV45_1 ~ NSFA_6 + NSFA_14 + NSFA_1'
# test3.form = 'UW_MEM_AV45_1 ~ CORTICAL_SUMMARY_prior + NSFA_14 + NSFA_0'
# test1.form = 'CORTICAL_SUMMARY_change ~ CORTICAL_SUMMARY_prior + APOE4_BIN'
# test2.form = 'CORTICAL_SUMMARY_change ~ NSFA_6 + NSFA_5'
# test3.form = 'CORTICAL_SUMMARY_change ~ NSFA_6 + NSFA_5'
# test1.form = 'ADAS_AV45_1 ~ CORTICAL_SUMMARY_prior'
# test2.form = 'ADAS_AV45_1 ~ NSFA_6 + NSFA_14'
# test3.form = 'ADAS_AV45_1 ~ CORTICAL_SUMMARY_prior + NSFA_14 + NSFA_0'
# test1.form = 'AVLT_AV45_1 ~ CORTICAL_SUMMARY_prior + Gender + Edu..Yrs.'
# test2.form = 'AVLT_AV45_1 ~ NSFA_6 + NSFA_14 + Gender + NSFA_1'
# test3.form = 'AVLT_AV45_1 ~ CORTICAL_SUMMARY_prior + NSFA_14 + Gender + NSFA_0'
test1.form = 'ADASslope_postAV45 ~ CORTICAL_SUMMARY_prior'
test2.form = 'ADASslope_postAV45 ~ NSFA_6 + NSFA_24 + NSFA_5'
test3.form = 'ADASslope_postAV45 ~ CORTICAL_SUMMARY_prior'

r2.shrinkage(test1.form, target, df_av45)
r2.shrinkage(test2.form, target, df_av45)
r2.shrinkage(test3.form, target ,df_av45)

toplot = 'NSFA_14'
ggplot(df_av45,aes_string(x='CORTICAL_SUMMARY_prior',y=target,color=toplot)) +
  geom_point(size=3) +
  geom_smooth(method='lm') +
  scale_color_gradient2(limit=c(-2,2), low='#22FF00',mid='white',high='#FF0000',oob=squish) +
  theme(panel.grid=element_blank(), 
        panel.background=element_rect(fill='black'))

avPlots(lm(test1.form,df_av45))
avPlots(lm(test2.form,df_av45))
avPlots(lm(test3.form,df_av45))





# ADAS likelihood test
adas_base_form = paste("ADAS_AV45_1",'~ 1')
adas_braak_form = paste("ADAS_AV45_1",'~',"CORTICAL_SUMMARY_prior + Gender + Edu..Yrs. + APOE4_BIN + Age.AV45")
adas_pattern_form = paste("ADAS_AV45_1",'~',"NSFA_6 + NSFA_14 + NSFA_0 + NSFA_1 + NSFA_21 + NSFA_5 + NSFA_20 + NSFA_24 + Gender + NSFA_4 + NSFA_25 + Edu..Yrs. + NSFA_19 + NSFA_22 + NSFA_9 + NSFA_26 + NSFA_7 + NSFA_17 + NSFA_13 + NSFA_3 + NSFA_28 + NSFA_10 + APOE4_BIN + NSFA_11 + NSFA_15")
adas_base_lm = lm(as.formula(adas_base_form),df_av45)
adas_braak_lm = lm(as.formula(adas_braak_form),df_av45)
adas_pattern_lm = lm(as.formula(adas_pattern_form),df_av45)
adas_braak_anova = anova(adas_base_lm, adas_braak_lm)
adas_pattern_anova = anova(adas_base_lm, adas_pattern_lm)

# AVLT likelihood test
avlt_base_form = paste("AVLT_AV45_1",'~ 1')
avlt_braak_form = paste("AVLT_AV45_1",'~',"CORTICAL_SUMMARY_prior + Gender + Age.AV45 + Edu..Yrs. + APOE4_BIN")
avlt_pattern_form = paste("AVLT_AV45_1",'~',"NSFA_6 + Gender + Age.AV45 + NSFA_0 + NSFA_14 + Edu..Yrs. + NSFA_1 + NSFA_20 + NSFA_21 + APOE4_BIN + NSFA_4 + NSFA_23 + NSFA_25 + NSFA_11 + NSFA_22 + NSFA_5 + NSFA_26 + NSFA_28 + NSFA_3 + NSFA_9 + NSFA_24 + NSFA_15 + NSFA_10 + NSFA_8 + NSFA_12")
avlt_base_lm = lm(as.formula(avlt_base_form),df_av45)
avlt_braak_lm = lm(as.formula(avlt_braak_form),df_av45)
avlt_pattern_lm = lm(as.formula(avlt_pattern_form),df_av45)
avlt_braak_anova = anova(avlt_base_lm, avlt_braak_lm)
avlt_pattern_anova = anova(avlt_base_lm, avlt_pattern_lm)

# AVLT plots
p1 = ggplot(df_av45, aes_string(x='NSFA_6', y='AVLT_AV45_1')) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') +
  theme(plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14),
        axis.text.x=element_text(face='bold', size=14),
        legend.title=element_blank()) +
  ggtitle('NSFA_6 Factor Score vs. AVLT') +
  xlab('NSFA_6 Factor Score') +
  ylab('AVLT')
print(p1)

p2 = ggplot(df_av45, aes_string(x='NSFA_14', y='AVLT_AV45_1')) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') +
  theme(plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14),
        axis.text.x=element_text(face='bold', size=14),
        legend.title=element_blank()) +
  ggtitle('NSFA_14 Factor Score vs. AVLT') +
  xlab('NSFA_14 Factor Score') +
  ylab('AVLT')
print(p2)

p3 = ggplot(df_av45, aes_string(x='CORTICAL_SUMMARY_prior', y='AVLT_AV45_1')) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') +
  theme(plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14),
        axis.text.x=element_text(face='bold', size=14),
        legend.title=element_blank()) +
  ggtitle('Summary SUVR vs. AVLT') +
  xlab('Summary SUVR (whole cereb. ref)') +
  ylab('AVLT')
print(p3)

# Adas plots
p1 = ggplot(df_av45, aes_string(x='NSFA_6', y='ADAS_AV45_1')) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') +
  theme(plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14),
        axis.text.x=element_text(face='bold', size=14),
        legend.title=element_blank()) +
  ggtitle('NSFA_6 Factor Score vs. ADAS-Cog') +
  xlab('NSFA_6 Factor Score') +
  ylab('ADAS-Cog')
print(p1)

p2 = ggplot(df_av45, aes_string(x='NSFA_14', y='ADAS_AV45_1')) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') +
  theme(plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14),
        axis.text.x=element_text(face='bold', size=14),
        legend.title=element_blank()) +
  ggtitle('NSFA_14 Factor Score vs. ADAS-Cog') +
  xlab('NSFA_14 Factor Score') +
  ylab('ADAS-Cog')
print(p2)

p3 = ggplot(df_av45, aes_string(x='CORTICAL_SUMMARY_prior', y='ADAS_AV45_1')) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') +
  theme(plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14),
        axis.text.x=element_text(face='bold', size=14),
        legend.title=element_blank()) +
  ggtitle('Summary SUVR vs. ADAS-Cog') +
  xlab('Summary SUVR (whole cereb. ref)') +
  ylab('ADAS-Cog')
print(p3)



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


