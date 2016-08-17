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
library(Jmisc)
library(lars)
library(covTest)

source('R/LM_FUNCS.R')

# CONSTANTS
pattern_prefix = 'NSFA_'
to_factor = c('RID','ad_prior','ad_post','positive_prior','positive_post',
              'diag_prior','diag_post','APOE4_BIN','APOE2_BIN','Gender',
              'Diag.AV45','Diag.AV1451','positive_prior','positive_post',
              'AV45_NONTP_wcereb_BIN1.11','AV45_NONTP_2_wcereb_BIN1.11',
              'AV45_NONTP_3_wcereb_BIN1.11')
to_standardize = c('Age.AV45','Edu..Yrs.')
demog_columns = c('RID','APOE4_BIN','Diag.AV1451','Age.AV1451','Gender','Edu..Yrs.')
diag_columns = c('Diag.AV45','Diag.AV1451')
# braak_columns = c('AV1451_Braak1_CerebGray_BL',
#                   'AV1451_Braak2_CerebGray_BL',
#                   'AV1451_Braak3_CerebGray_BL',
#                   'AV1451_Braak4_CerebGray_BL',
#                   'AV1451_Braak5_CerebGray_BL',
#                   'AV1451_Braak6_CerebGray_BL')
braak_columns = c('AV1451_Braak12_CerebGray_BL',
                  'AV1451_Braak34_CerebGray_BL',
                  'AV1451_Braak56_CerebGray_BL')

# target = "UW_EF_AV1451_1"
# target = "UW_MEM_AV1451_1"
# target = "ADAS_AV1451_1"
# target = "AVLT_AV1451_1"
target = "AV1451_BL_closest_AV45_wcereb_BIN1.11"

# target = "ADAS_retroslope_AV1451_BL"
# target = "AVLT_retroslope_AV1451_BL"
# target = "UW_MEM_retroslope_AV1451_BL"
# target = "UW_EF_retroslope_AV1451_BL"
# target = "AV1451_BL_closest_AV45_wcereb_retroSlope"

output_folder = 'R/output_av1451/'

valid_diags = c('N','SMC','EMCI','LMCI','AD')
#valid_diags = c('N','SMC','EMCI','LMCI')
#valid_diags = c('N','SMC')
#valid_diags = c('EMCI')
#valid_diags = c('LMCI')

# IMPORT
df_av1451 = read.csv('nsfa/av1451skull_pattern_dataset.csv')
pattern_columns = Filter(isPatternColumn,names(df_av1451))
naive_columns = Filter(isNaiveColumn,names(df_av1451))
non.na = complete.cases(df_av1451[,c(demog_columns,braak_columns,target)])
df_av1451 = df_av1451[non.na,]
for (i in names(df_av1451)){
  if (i %in% to_factor){
    df_av1451[,eval(i)] = as.factor(as.character(df_av1451[,eval(i)]))
  }
}

# remove target outliers
# target.mean = mean(df_av1451[,target])
# target.sd = sd(df_av1451[,target])
# df_av1451 = df_av1451[df_av1451[,target] <= target.mean+target.sd*5,]
# df_av1451 = df_av1451[df_av1451[,target] >= target.mean-target.sd*5,]

# Filter by diag
df_av1451 = df_av1451[which(df_av1451$Diag.AV1451 %in% valid_diags),]

# standardize predictors
# cross_to_standardize = c(to_standardize,pattern_columns,naive_columns,braak_columns,target)
# cross_normalization = preProcess(df_av1451[,cross_to_standardize])
# df_av1451[,cross_to_standardize] = predict(cross_normalization, df_av1451[,cross_to_standardize])
# df_av1451[,target] = scale(df_av1451[,target], center=TRUE, scale=FALSE)

# # look at histograms
# for (pcol in pattern_columns) {
#   p = ggplot(df_av1451, aes_string(pcol)) + geom_histogram(binwidth=0.1)
#   print(p)
# }

# make crossx response normal
#df_av1451[,eval(target)] = Gaussianize(df_av1451[,eval(target)], type='hh', method='MLE', return.u=TRUE)

# Formula setup
# all.addons = lapply(pattern_columns,lm.addvar)
# naive.addons = lapply(naive_columns,lm.addvar)
# braak.addons = lapply(braak_columns,lm.addvar)
all.addons = paste('+',paste(pattern_columns,collapse=' + '))
naive.addons = paste('+',paste(naive_columns,collapse=' + '))
braak.addons = paste('+',paste(braak_columns,collapse=' + '))
patterns_str = paste(all.addons,collapse=' ')
naive_str = paste(naive.addons,collapse=' ')
braak_str = paste(braak.addons,collapse=' ')

# diag_str = 'Diag.AV1451*APOE4_BIN +'
# diag_str = 'Diag.AV1451 +'
diag_str = ''



base_form = paste(target,"~",diag_str,"APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.")
braak_form = paste(target,"~",diag_str,"APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",braak_str)
pattern_form = paste(target,"~",diag_str,"APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",patterns_str)
naive_form = paste(target,"~",diag_str,"APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",naive_str)
full_form = paste(target,"~",diag_str,"APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",patterns_str,braak_str)
# onlypattern_form = str_replace(paste(target,"~",paste(all.addons,collapse=' ')),"\\+ ","")

# LARS lasso
braak_x = getxy(braak_form,df_av1451)
y = as.numeric(df_av1451[,target])
braak.lars.model = lars(braak_x,y,type='lasso')
braak.lars.test = covTest(braak.lars.model,braak_x,y)$results
braak.lars.test = data.frame(braak.lars.test[complete.cases(braak.lars.test),])
braak.lars.coef = coef(braak.lars.model, s=which.min(summary(braak.lars.model)$Cp), mode='step')
braak.lars.test$name = names(braak.lars.coef[braak.lars.test[,'Predictor_Number']])
braak.lars.test$coef = braak.lars.coef[braak.lars.test[,'Predictor_Number']]
braak.lars.sig = braak.lars.test[braak.lars.test$P.value <= 0.05,]
braak.lars.cp = min(summary(braak.lars.model)$Cp)
braak.lars.r2 = braak.lars.model$R2[which.min(summary(braak.lars.model)$Cp)]
braak.lars.n = attr(summary(braak.lars.model)$Cp,'n')
braak.lars.p = NROW(braak.lars.coef[braak.lars.coef != 0])
braak.lars.r2adj = r2adj(braak.lars.r2,braak.lars.n,braak.lars.p)
braak.lars.nonzero = braak.lars.test[braak.lars.test$coef != 0,'name']
paste(braak.lars.nonzero,collapse=' + ')
paste(target,'~',paste(braak.lars.sig[,'name'], collapse=' + '))

pattern_x = getxy(pattern_form,df_av1451)
y = as.numeric(df_av1451[,target])
pattern.lars.model = lars(pattern_x,y,type='lasso')
pattern.lars.test = covTest(pattern.lars.model,pattern_x,y)$results
pattern.lars.test = data.frame(pattern.lars.test[complete.cases(pattern.lars.test),])
pattern.lars.coef = coef(pattern.lars.model, s=which.min(summary(pattern.lars.model)$Cp), mode='step')
pattern.lars.test$name = names(pattern.lars.coef[pattern.lars.test[,'Predictor_Number']])
pattern.lars.test$coef = pattern.lars.coef[pattern.lars.test[,'Predictor_Number']]
pattern.lars.sig = pattern.lars.test[pattern.lars.test$P.value <= 0.05,]
pattern.lars.cp = min(summary(pattern.lars.model)$Cp)
pattern.lars.r2 = pattern.lars.model$R2[which.min(summary(pattern.lars.model)$Cp)]
pattern.lars.n = attr(summary(pattern.lars.model)$Cp,'n')
pattern.lars.p = NROW(pattern.lars.coef[pattern.lars.coef != 0])
pattern.lars.r2adj = r2adj(pattern.lars.r2,pattern.lars.n,pattern.lars.p)
pattern.lars.nonzero = pattern.lars.test[pattern.lars.test$coef != 0,'name']
paste(pattern.lars.nonzero,collapse=' + ')
paste(target,'~',paste(pattern.lars.sig[,'name'], collapse=' + '))

full_x = getxy(full_form,df_av1451)
y = as.numeric(df_av1451[,target])
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
# test1.form = "ADAS_AV1451_1 ~ AV1451_Braak56_CerebGray_BL"
# test2.form = "ADAS_AV1451_1 ~ NSFA_0 + NSFA_3 + NSFA_8"
# test3.form = "ADAS_AV1451_1 ~ NSFA_0 + NSFA_3 + NSFA_8"
# test1.form = "AVLT_AV1451_1 ~ Age.AV1451"
# test2.form = "AVLT_AV1451_1 ~ NSFA_3"
# test3.form = "AVLT_AV1451_1 ~ NSFA_3"
# test1.form = "ADAS_retroslope_AV1451_BL ~ AV1451_Braak56_CerebGray_BL + APOE4_BIN"
# test2.form = "ADAS_retroslope_AV1451_BL ~ NSFA_0"
# test3.form = "ADAS_retroslope_AV1451_BL ~ NSFA_0"
test2.form = "AV1451_BL_closest_AV45_wcereb_retroSlope ~ NSFA_7"

test1_results = c()
test2_results = c()
test3_results = c()
for (i in 1:50) {
  test1_r2 = r2.shrinkage(test1.form, target, df_av1451)[2]
  test2_r2 = r2.shrinkage(test2.form, target, df_av1451)[2]
  test3_r2 = r2.shrinkage(test3.form, target ,df_av1451)[2]
  test1_results = c(test1_results,test1_r2)
  test2_results = c(test2_results,test2_r2)
  test3_results = c(test3_results,test3_r2)
}
mean(test1_results)
mean(test2_results)
mean(test3_results)



# Penalized LM
braak.lasso.model = run.lasso(braak_form,df_av1451,'RMSE')
braak.lasso.metric = subset(braak.lasso.model$results, fraction == braak.lasso.model$bestTune$fraction)
braak.lasso.coef = predict.enet(braak.lasso.model$finalModel, type='coefficients',s=braak.lasso.model$bestTune$fraction, mode='fraction')$coefficients
braak.lasso.coef = braak.lasso.coef[braak.lasso.coef != 0]

pattern.lasso.model = run.lasso(pattern_form,df_av1451,'RMSE')
pattern.lasso.metric = subset(pattern.lasso.model$results, fraction == pattern.lasso.model$bestTune$fraction)
pattern.lasso.coef = predict.enet(pattern.lasso.model$finalModel, type='coefficients',s=pattern.lasso.model$bestTune$fraction, mode='fraction')$coefficients
pattern.lasso.coef = pattern.lasso.coef[pattern.lasso.coef != 0]

naive.lasso.model = run.lasso(naive_form,df_av1451,'RMSE')
naive.lasso.metric = subset(naive.lasso.model$results, fraction == naive.lasso.model$bestTune$fraction)
naive.lasso.coef = predict.enet(naive.lasso.model$finalModel, type='coefficients',s=naive.lasso.model$bestTune$fraction, mode='fraction')$coefficients
naive.lasso.coef = naive.lasso.coef[naive.lasso.coef != 0]

full.lasso.model = run.lasso(full_form,df_av1451,'Rsquared')
full.lasso.metric = subset(full.lasso.model$results, fraction == full.lasso.model$bestTune$fraction)
full.lasso.coef = predict.enet(full.lasso.model$finalModel, type='coefficients',s=full.lasso.model$bestTune$fraction, mode='fraction')$coefficients
full.lasso.coef = full.lasso.coef[full.lasso.coef != 0]




# GLM net

braak.glmnet.model = run.glmnet(braak_form,df_av1451,'Rsquared')
braak.glmnet.metric = subset(braak.glmnet.model$results, alpha == braak.glmnet.model$bestTune$alpha & lambda == braak.glmnet.model$bestTune$lambda)
braak.glmnet.coef = predict.glmnet(braak.glmnet.model$finalModel,type='coefficients',s=braak.glmnet.model$bestTune$lambda)

naive.glmnet.model = run.glmnet(naive_form,df_av1451,'Rsquared')
naive.glmnet.metric = subset(naive.glmnet.model$results, alpha == naive.glmnet.model$bestTune$alpha & lambda == naive.glmnet.model$bestTune$lambda)
naive.glmnet.coef = predict.glmnet(naive.glmnet.model$finalModel,type='coefficients',s=naive.glmnet.model$bestTune$lambda)

pattern.glmnet.model = run.glmnet(pattern_form,df_av1451,'Rsquared')
pattern.glmnet.metric = subset(pattern.glmnet.model$results, alpha == pattern.glmnet.model$bestTune$alpha & lambda == pattern.glmnet.model$bestTune$lambda)
pattern.glmnet.coef = predict.glmnet(pattern.glmnet.model$finalModel,type='coefficients',s=pattern.glmnet.model$bestTune$lambda)

full.glmnet.model = run.glmnet(full_form,df_av1451,'Rsquared')
full.glmnet.metric = subset(full.glmnet.model$results, alpha == full.glmnet.model$bestTune$alpha & lambda == full.glmnet.model$bestTune$lambda)
full.glmnet.coef = predict.glmnet(full.glmnet.model$finalModel,type='coefficients',s=full.glmnet.model$bestTune$lambda)



# RFE


rfe.onlypattern = run.rfe(onlypattern_form, target, df_av1451, 2)
optvars = rfe.onlypattern$optVariables

pattern_form = paste(target,"~",diag_str,"Age.AV1451 + Gender + Edu..Yrs. + APOE4_BIN +",paste(optvars,collapse=' + '))
fm_base = lm(as.formula(base_form),data=df_av1451)
fm_braak = lm(as.formula(braak_form),data=df_av1451)
fm_pattern = lm(as.formula(pattern_form),data=df_av1451)

rfe.base = run.rfe(base_form, target, df_av1451, 2)
rfe.braak = run.rfe(braak_form, target, df_av1451, 2)
rfe.pattern = run.rfe(pattern_form, target, df_av1451, 2)
fm_base = rfe.base$fit
fm_braak = rfe.braak$fit
fm_pattern = rfe.pattern$fit

fm_base.summary = summary(fm_base)
fm_braak.summary = summary(fm_braak)
fm_pattern.summary = summary(fm_pattern)

fm_base.fit = sem.model.fits(fm_base)
fm_braak.fit = sem.model.fits(fm_braak)
fm_pattern.fit = sem.model.fits(fm_pattern)

fm_base.anova = Anova(fm_base,type='III')
fm_braak.anova = Anova(fm_braak,type='III')
fm_pattern.anova = Anova(fm_pattern,type='III')

save.printout(paste(output_folder,target,'_fm_base_summary','.txt',sep=''),fm_base.summary)
save.printout(paste(output_folder,target,'_fm_pattern_summary','.txt',sep=''),fm_pattern.summary)
save.printout(paste(output_folder,target,'_fm_braak_summary','.txt',sep=''),fm_braak.summary)
save.printout(paste(output_folder,target,'_fm_base_fit','.txt',sep=''),fm_base.fit)
save.printout(paste(output_folder,target,'_fm_pattern_fit','.txt',sep=''),fm_pattern.fit)
save.printout(paste(output_folder,target,'_fm_braak_fit','.txt',sep=''),fm_braak.fit)
save.printout(paste(output_folder,target,'_fm_base_anova','.txt',sep=''),fm_base.anova)
save.printout(paste(output_folder,target,'_fm_pattern_anova','.txt',sep=''),fm_pattern.anova)
save.printout(paste(output_folder,target,'_fm_braak_anova','.txt',sep=''),fm_braak.anova)
save(rfe.base,file=paste(output_folder,target,'_rfe_base_obj',sep=''))
save(rfe.pattern,file=paste(output_folder,target,'_rfe_pattern_obj',sep=''))
save(rfe.braak,file=paste(output_folder,target,'_rfe_braak_obj',sep=''))

fm_base.plotfn = function() {par(mfrow=c(2,2));plot(fm_base);title("Base Model", outer=T, line=-2);}
fm_pattern.plotfn = function() {par(mfrow=c(2,2));plot(fm_pattern);title("Pattern Model", outer=T, line=-2);}
fm_braak.plotfn = function() {par(mfrow=c(2,2));plot(fm_braak);title("Braak Model", outer=T, line=-2);}
save.plot(paste(output_folder,target,'_fm_base_lmplot.pdf',sep=''), fm_base.plotfn)
save.plot(paste(output_folder,target,'_fm_pattern_lmplot.pdf',sep=''), fm_pattern.plotfn)
save.plot(paste(output_folder,target,'_fm_braak_lmplot.pdf',sep=''), fm_braak.plotfn)

save.plot(paste(output_folder,target,'_fm_base_avplot.pdf',sep=''), function() {avPlots(fm_base, ask=FALSE)})
save.plot(paste(output_folder,target,'_fm_pattern_avplot.pdf',sep=''), function() {avPlots(fm_pattern, ask=FALSE)})
save.plot(paste(output_folder,target,'_fm_braak_avplot.pdf',sep=''), function() {avPlots(fm_braak, ask=FALSE)})

fm_base.summary
fm_pattern.summary
fm_braak.summary

fm_base.fit
fm_pattern.fit
fm_braak.fit

