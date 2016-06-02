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

source('R/LM_FUNCS.R')

# CONSTANTS
pattern_prefix = 'NSFA_'
to_factor = c('RID','ad_prior','ad_post','positive_prior','positive_post',
              'diag_prior','diag_post','APOE4_BIN','APOE2_BIN','Gender',
              'Diag.AV45_long','positive_prior','positive_post',
              'AV45_NONTP_wcereb_BIN1.11')
to_standardize = c('CORTICAL_SUMMARY_change','CORTICAL_SUMMARY_prior','CORTICAL_SUMMARY_post','Age.AV45','Edu..Yrs.')
demog_columns = c('RID','APOE4_BIN','diag_prior','Age.AV45','Gender','Edu..Yrs.')
av45_columns = c('CORTICAL_SUMMARY_change','CORTICAL_SUMMARY_prior','CORTICAL_SUMMARY_post','positive_prior','positive_post')

#target = "CORTICAL_SUMMARY_change"
#target = "UW_EF_BL_3months"
#target = "UW_EF_slope"
#target = "ADAS_3MTH_AV45"
#target = "ADASslope_postAV45"
#target = "AVLT_AV45_1_3MTHS"
#target = "AVLTslope_postAV45"
#target = "UW_MEM_BL_3months"
#target = "UW_MEM_slope"
target = "CSF_ABETA_closest_AV45_1"
#target = "CSF_ABETA_slope"
#target = "CSF_TAU_closest_AV45_1"
#target = "CSF_TAU_slope"
#target = "CSF_PTAU_closest_AV45_1"
#target = "CSF_PTAU_slope"

#output_folder = 'R/output/'
output_folder = 'R/output_all_diag/'
output_folder = 'R/output_neg_emci/'

positive_value=0
all_diags = c('N','SMC','EMCI','LMCI','AD')
#valid_diags = c('N','SMC','EMCI','LMCI','AD')
#valid_diags = c('N','SMC','EMCI','LMCI')
#valid_diags = c('N','SMC')
valid_diags = c('EMCI')
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
df_av45 = read.csv('nsfa/pattern_dataset.csv')

pattern_columns = Filter(isPatternColumn,names(df_av45))

df_av45 = df_av45[which(df_av45[,'AV45_NONTP_wcereb_BIN1.11'] == positive_value),]
df_av45 = df_av45[which(df_av45$diag_prior %in% valid_diags),]
for (i in names(df_av45)){
  if (i %in% to_factor){
    df_av45[,eval(i)] = as.factor(as.character(df_av45[,eval(i)]))
  }
}
df_long = to.long(df_av45, time_col_prefix, value_col_prefix)

# look at histograms
for (pcol in pattern_columns) {
  p = ggplot(df_av45, aes_string(pcol, colour = 'positive_prior')) + geom_histogram(binwidth=0.1)
  print(p)
}

# look at change
for (pcol in pattern_columns) {
  scan2_col = paste('SCAN2_',pcol,sep='')
  change = (df_av45[,eval(scan2_col)]-df_av45[,eval(pcol)])/(df_av45$AV45_1_2_Diff)
  plot(df_av45[,eval(pcol)],change,xlab=pcol,ylab=pcol)
  #plot(df_av45$CORTICAL_SUMMARY_change,change,xlab='CS Prior',ylab=pcol)
  abline(0,0,col='red')
}


# standardize predictors
# cross_to_standardize = c(to_standardize,pattern_columns,target)
# long_to_standardize = c(to_standardize,pattern_columns,'time','value')
# cross_normalization = preProcess(df_av45[,cross_to_standardize])
# long_normalization = preProcess(df_long[,long_to_standardize])
# df_av45[,cross_to_standardize] = predict(cross_normalization, df_av45[,cross_to_standardize])
# df_long[,long_to_standardize] = predict(long_normalization, df_long[,long_to_standardize])

# make long response normal
df_long$value = Gaussianize(df_long$value, type='hh', method='MLE', return.u=TRUE)

# make crossx response normal
non.na = complete.cases(df_av45[,c(demog_columns,av45_columns,target)])
df_av45 = df_av45[non.na,]
df_av45[,eval(target)] = Gaussianize(df_av45[,eval(target)], type='hh', method='MLE', return.u=TRUE)

# One by one pattern variable likelihood testing
mlm.testvar = function(var.name) {
  #CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + 
  #CORTICAL_SUMMARY_prior*APOE4_BIN + positive_prior*APOE4_BIN
  base_str = paste(target,"~","APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
  add_str = lm.addvar(var.name)
  form_base = as.formula(base_str)
  form = as.formula(paste(base_str,add_str))
  fm = vglm(as.formula(form), family=multinomial(refLevel=1), data=df_av45)
  fm_base = vglm(as.formula(form_base), family=multinomial(refLevel=1), data=df_av45)
  like = VGAM::lrtest(fm, fm_base)
  like.p = like@Body$`Pr(>Chisq)`[2]
}

mlm.cv = function(dataset, form, target) {
  k = 10
  folds = cvFolds(nrow(dataset), K=k)
  holdoutpred = rep(0,nrow(dataset))
  for (i in 1:k) {
    train = dataset[folds$subsets[folds$which != i],]
    validation = dataset[folds$subsets[folds$which == i],]
    newlm = vglm(form, family=multinomial(refLevel=1), data=train)
    newprobs = VGAM::predict(newlm, validation, type='response')
    if (is.null(colnames(newprobs))) {
      newpred = apply(newprobs,1,which.max)
    } else {
      newpred = colnames(newprobs)[apply(newprobs,1,which.max)]
    }
    holdoutpred[folds$subsets[folds$which ==i]] = newpred
  }
  responses = as.numeric(dataset[,eval(target)])
  # find accuracy
  ctab = xtabs(as.formula(paste('~',target,'+ holdoutpred')), data=dataset)
  sum(diag(ctab))/sum(ctab)
}

lme.testvar = function(var.name) {
  #diag_prior*time + CORTICAL_SUMMARY_prior*time + 
  base_str = paste('value',"~","APOE4_BIN*time + Age.AV45 + Gender + Edu..Yrs.")
  add_str = lme.addvar(var.name)
  random_str = '+ (1 + time | RID)'
  form_base = as.formula(paste(base_str,random_str))
  form = as.formula(paste(base_str,add_str,random_str))
  fm = lmer(form,data=df_long)
  fm_base = lmer(form_base,df_long)
  like = anova(fm_base,fm)
  like.p = like$`Pr(>Chisq)`[2]
}

lme.cv = function(dataset, form) {
  k = 20
  subjects = levels(as.factor(dataset$RID))
  folds = cvFolds(length(subjects), K=k)
  holdoutpred = rep(0,nrow(dataset))
  for (i in 1:k) {
    train_subjects = subjects[folds$subsets[folds$which != i]]
    validation_subjects = subjects[folds$subsets[folds$which == i]]
    train = dataset[dataset$RID %in% train_subjects,]
    validation_indices = dataset$RID %in% validation_subjects
    validation = dataset[validation_indices,]
    newlm = lmer(as.formula(form),df_long)
    newpred = predict(newlm, newdata=validation)
    holdoutpred[validation_indices] = newpred
  }
  responses = dataset[,'value']
  rmse(holdoutpred,round(responses))
}

run.rfe = function(form, var.response, dataset) {
  x = as.data.frame(model.matrix(as.formula(form),dataset))[,-1]
  nzv_cols = nearZeroVar(x)
  if (length(nzv_cols) > 0) {
    x = x[, -nzv_cols]
  }
  corr_cols = findCorrelation(cor(x),.8)
  if (length(corr_cols) > 0) {
    x = x[, -corr_cols]
  }

  colnames(x) = lapply(colnames(x), make.names)
  rownames(x) = NULL
  y = as.numeric(df_av45[,var.response])
  
  ctrl = rfeControl(functions = lmFuncs, 
                    method = "repeatedcv", 
                    number = 10,
                    repeats = 5,
                    rerank = TRUE,
                    verbose = FALSE)
  set.seed(1)
  rfe.output = rfe(x, 
                   y, 
                   sizes = c(2:ncol(x)),
                   rfeControl = ctrl,
                   metric = 'Rsquared')
  rfe.output
}


# LM RFE
all.addons = lapply(pattern_columns,lm.addvar)
addons_form = str_replace(paste(target,"~",paste(all.addons,collapse=' ')),"\\+ ","")
# if (sum(df_av45$diag_prior=='SMC') == 0) {
#   diag_str = 'APOE4_BIN +'
# } else {
#   diag_str = 'diag_prior*APOE4_BIN +'
# }
diag_str = 'diag_prior*APOE4_BIN +'

base_form = paste(target,"~",diag_str,"Age.AV45 + Gender + Edu..Yrs.")
nopattern_form = paste(target,"~",diag_str,"CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
onlypattern_form = paste(base_form,paste(all.addons,collapse=' '))
full_form = paste(nopattern_form,paste(all.addons,collapse=' '))

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







# LM STEPAIC
all.addons = lapply(pattern_columns,lm.addvar)
base_form = paste(target,"~ diag_prior*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
nopattern_form = paste(target,"~ diag_prior*APOE4_BIN + CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
onlypattern_form = paste(base_form,paste(all.addons,collapse=' '))

fm_base = lm(as.formula(base_form),data=df_av45)
fm_nopattern = lm(as.formula(nopattern_form),df_av45)
fm_onlypattern = lm(as.formula(onlypattern_form),df_av45)

trControl = trainControl(method = 'LGOCV', number = 10)
fm_base.lmfit = train(as.formula(base_form), data=df_av45, method='lmStepAIC', trControl=trControl)
fm_nopattern.lmfit = train(as.formula(nopattern_form), data=df_av45, method='lmStepAIC', trControl=trControl)
fm_onlypattern.lmfit = train(as.formula(onlypattern_form), data=df_av45, method='lmStepAIC', trControl=trControl)
fm_base.step = fm_base.lmfit$finalModel
fm_nopattern.step = fm_nopattern.lmfit$finalModel
fm_onlypattern.step = fm_onlypattern.lmfit$finalModel

fm_base.step = stepAIC(fm_base, direction='both')
fm_nopattern.step = stepAIC(fm_nopattern, direction='both')
fm_onlypattern.step = stepAIC(fm_onlypattern, direction='both')

fm_base.summary = summary(fm_base.step)
fm_nopattern.summary = summary(fm_nopattern.step)
fm_onlypattern.summary = summary(fm_onlypattern.step)


fm_base.fit = sem.model.fits(fm_base.step)
fm_nopattern.fit = sem.model.fits(fm_nopattern.step)
fm_onlypattern.fit = sem.model.fits(fm_onlypattern.step)

fm_base.anova = Anova(fm_base.step,type='III')
fm_nopattern.anova = Anova(fm_nopattern.step,type='III')
fm_onlypattern.anova = Anova(fm_onlypattern.step,type='III')

fm_base.fit
fm_nopattern.fit
fm_onlypattern.fit

par(mfrow=c(2,2))
plot(fm_base.step)
title("Base Model", outer=T, line=-2)
plot(fm_nopattern.step)
title("No Pattern Model", outer=T, line=-2)
plot(fm_onlypattern.step)
title("Only Pattern Model", outer=T, line=-2)


#params = expand.grid(alpha=seq(0,1,by=0.1),lambda=2^seq(1,-10, by=-0.3))
#lmFit = train(as.formula(onlypattern_form), data=df_av45, method="glmnet", trControl=trControl, tuneGrid=params)

# LM (Fixed Effects)
like.pvalues = lapply(pattern_columns,lm.testvar)
valid_patterns = pattern_columns[like.pvalues <= 0.05]
form.addons = lapply(valid_patterns,lm.addvar)
all.addons = lapply(pattern_columns,lm.addvar)
base_form = paste(target,"~ diag_prior*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
nopattern_form = paste(target,"~ diag_prior*APOE4_BIN + CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
onlypattern_form = paste(base_form,paste(form.addons,collapse=' '))
full_form = paste(nopattern_form,paste(form.addons,collapse=' '))
allpattern_form = paste(base_form,paste(all.addons,collapse=' '))

fm_base = lm(as.formula(base_form),data=df_av45)
fm_nopattern = lm(as.formula(nopattern_form),df_av45)
fm_onlypattern = lm(as.formula(onlypattern_form),df_av45)
fm_full = lm(as.formula(full_form),df_av45)
fm_all = lm(as.formula(allpattern_form),df_av45)

fm_base.summary = summary(fm_base)
fm_nopattern.summary = summary(fm_nopattern)
fm_onlypattern.summary = summary(fm_onlypattern)
fm_full.summary = summary(fm_full)

fm_base.fit = sem.model.fits(fm_base)
fm_nopattern.fit = sem.model.fits(fm_nopattern)
fm_onlypattern.fit = sem.model.fits(fm_onlypattern)
fm_full.fit = sem.model.fits(fm_full)

fm_base.anova = Anova(fm_base,type='III')
fm_nopattern.anova = Anova(fm_nopattern,type='III')
fm_onlypattern.anova = Anova(fm_onlypattern,type='III')
fm_full.anova = Anova(fm_full,type='III')

fm_full.cv = cv.lm(data=df_av45,form.lm=fm_full, m=20, printit=FALSE)
fm_base.cv = cv.lm(data=df_av45,form.lm=fm_base, m=20, printit=FALSE)
fm_nopattern.cv = cv.lm(data=df_av45,form.lm=fm_nopattern, m=20, printit=FALSE)
fm_onlypattern.cv = cv.lm(data=df_av45,form.lm=fm_onlypattern, m=20, printit=FALSE)

attributes(fm_full.cv)$ms
attributes(fm_onlypattern.cv)$ms
attributes(fm_nopattern.cv)$ms
attributes(fm_base.cv)$ms

fm_base.fit
fm_nopattern.fit
fm_onlypattern.fit
fm_full.fit

fm_base.summary$adj.r.squared
fm_nopattern.summary$adj.r.squared
fm_onlypattern.summary$adj.r.squared
fm_full.summary$adj.r.squared

base_form
nopattern_form
onlypattern_form
full_form

anova(fm_full,fm_base)
anova(fm_onlypattern,fm_base)
anova(fm_nopattern,fm_base)


# MLM (multinomial linear regression)
like.pvalues = lapply(pattern_columns,mlm.testvar)
valid_patterns = pattern_columns[like.pvalues <= 0.05]
form.addons = lapply(valid_patterns,lm.addvar)
base_form = paste(target,"~ APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
#CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN
nopattern_form = paste(target,"~ CORTICAL_SUMMARY_prior*APOE4_BIN + positive_prior*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
onlypattern_form = paste(base_form,paste(form.addons,collapse=' '))
full_form = paste(nopattern_form,paste(form.addons,collapse=' '))
null_form = paste(target,'~ 1')

fm_null = vglm(as.formula(null_form), family=multinomial(refLevel=1), data=df_av45)
fm_base = vglm(as.formula(base_form), family=multinomial(refLevel=1), data=df_av45)
fm_nopattern = vglm(as.formula(nopattern_form), family=multinomial(refLevel=1), data=df_av45)
fm_onlypattern = vglm(as.formula(onlypattern_form), family=multinomial(refLevel=1), data=df_av45)
fm_full = vglm(as.formula(full_form), family=multinomial(refLevel=1), data=df_av45)

fm_base.aic = VGAM::AIC(fm_base)
fm_nopattern.aic = VGAM::AIC(fm_nopattern)
fm_onlypattern.aic = VGAM::AIC(fm_onlypattern)
fm_full.aic = VGAM::AIC(fm_full)

fm_base.deviance = VGAM::deviance(fm_base)
fm_nopattern.deviance = VGAM::deviance(fm_nopattern)
fm_onlypattern.deviance = VGAM::deviance(fm_onlypattern)
fm_full.deviance = VGAM::deviance(fm_full)

fm_null.ll = VGAM::logLik(fm_null)
fm_base.ll = VGAM::logLik(fm_base)
fm_nopattern.ll = VGAM::logLik(fm_nopattern)
fm_onlypattern.ll = VGAM::logLik(fm_onlypattern)
fm_full.ll = VGAM::logLik(fm_full)

fm_base.summary = VGAM::summary(fm_base)
fm_nopattern.summary = VGAM::summary(fm_nopattern)
fm_onlypattern.summary = VGAM::summary(fm_onlypattern)
fm_full.summary = VGAM::summary(fm_full)

fm_full.ccr = mlm.cv(df_av45,as.formula(full_form),target)
fm_onlypattern.ccr = mlm.cv(df_av45,as.formula(onlypattern_form),target)
fm_nopattern.ccr = mlm.cv(df_av45,as.formula(nopattern_form),target)
fm_base.ccr = mlm.cv(df_av45,as.formula(base_form),target)


fm_base.aic
fm_nopattern.aic
fm_onlypattern.aic
fm_full.aic

as.vector(1 - (fm_base.ll / fm_null.ll))
as.vector(1 - (fm_nopattern.ll / fm_null.ll))
as.vector(1 - (fm_onlypattern.ll / fm_null.ll))
as.vector(1 - (fm_full.ll / fm_null.ll))

fm_full.ccr
fm_onlypattern.ccr
fm_nopattern.ccr
fm_base.ccr


# LME (Mixed Effects)
like.pvalues = lapply(pattern_columns,lme.testvar)
valid_patterns = pattern_columns[like.pvalues <= 0.05]
form.addons = lapply(valid_patterns,lme.addvar)
random_str = '+ (1 + time | RID)'
#diag_prior*time + 
base_str = "value ~ APOE4_BIN*time + Age.AV45 + Gender + Edu..Yrs."
nopattern_str = "value ~ CORTICAL_SUMMARY_prior*time + APOE4_BIN*time + Age.AV45 + Gender + Edu..Yrs."
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

fm_base.fit = sem.model.fits(fm_base)
fm_nopattern.fit = sem.model.fits(fm_nopattern)
fm_onlypattern.fit = sem.model.fits(fm_onlypattern)
fm_full.fit = sem.model.fits(fm_full)

fm_base.anova = Anova(fm_base,type='III')
fm_nopattern.anova = Anova(fm_nopattern,type='III')
fm_onlypattern.anova = Anova(fm_onlypattern,type='III')
fm_full.anova = Anova(fm_full,type='III')
full_base.anova = anova(fm_full,fm_base)
onlypattern_base.anova = anova(fm_onlypattern,fm_base)
nopattern_base.anova = anova(fm_nopattern,fm_base)

fm_base.fit
fm_nopattern.fit
fm_onlypattern.fit
fm_full.fit

full_base.anova$`Pr(>Chisq)`[2]
onlypattern_base.anova$`Pr(>Chisq)`[2]
nopattern_base.anova$`Pr(>Chisq)`[2]

fm_base.cv = lme.cv(df_long, as.formula(base_form))
fm_nopattern.cv = lme.cv(df_long, as.formula(nopattern_form))
fm_onlypattern.cv = lme.cv(df_long, as.formula(onlypattern_form))
fm_full.cv = lme.cv(df_long, as.formula(full_form))

fm_base.summary$logLik[1]
fm_nopattern.summary$logLik[1]
fm_onlypattern.summary$logLik[1]
fm_full.summary$logLik[1]

fm_base.cv
fm_nopattern.cv
fm_onlypattern.cv
fm_full.cv

base_form
nopattern_form
onlypattern_form
full_form


