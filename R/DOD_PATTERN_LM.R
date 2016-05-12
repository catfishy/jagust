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

# # Calculate time to positivity threshold
# df_av45['threshold'] = threshold
# df_av45['threshold_diff'] = df_av45$threshold - df_av45$CORTICAL_SUMMARY_prior
# df_av45['to_threshold'] = df_av45$threshold_diff / df_av45$CORTICAL_SUMMARY_slope
# df_av45[which(df_av45$to_threshold < 0),'to_threshold'] = 0

# CONSTANTS
pattern_prefix = 'NSFA_'
to_factor = c('SCRNO','ad_prior','ad_post','positive_prior','positive_post',
              'diag_prior','diag_post','APOE4BIN','Sex',
              'Diag_closest_AV45_BL','positive_prior','positive_post',
              'ANTIDEP_USE','SSRI','PTGroup','GroupNum','GroupNum_TBI','GroupNum_PTSD')
to_standardize = c('CORTICAL_SUMMARY_change','CORTICAL_SUMMARY_prior','CORTICAL_SUMMARY_post','Age','Edu')
demog_columns = c('SCRNO','APOE4BIN','diag_prior','Age','Sex','Edu')
av45_columns = c('CORTICAL_SUMMARY_prior','positive_prior')

target = "ADAS_closest_AV45_BL"
target = "ADAS_postAV45_SLOPE"
target = "AVLT_closest_AV45_BL"
target = "AVLT_postAV45_SLOPE"
target = "GD_closest_AV45_BL"
target = "GD_postAV45_SLOPE"
target = "CDR_GLOBAL_closest_AV45_BL"
target = "CDR_GLOBAL_postAV45_SLOPE"


valid_diags = c('N')
#valid_diags = c('MCI')


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

save.printout = function(output_file, obj) {
  sink(output_file); print(obj, correlation=TRUE); sink(file=NULL)
}

save.plot = function(output_file, plot_fn) {
  pdf(file=output_file);plot_fn();dev.off();
}


plot.model = function(model) {
  plot(model)
}

fm.relimp = function(model) {
  calc.relimp(model, type=c('lmg'), rela=TRUE)
}

rmse = function(m, o) {
  sqrt(mean((m-o)^2))
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
df_av45 = read.csv('nsfa/dod_pattern_dataset.csv')
df_av45 = df_av45[which(df_av45$diag_prior %in% valid_diags),]
for (i in names(df_av45)){
  if (i %in% to_factor){
    df_av45[,eval(i)] = as.factor(as.character(df_av45[,eval(i)]))
  }
}
pattern_columns = Filter(isPatternColumn,names(df_av45))

# standardize predictors
# cross_to_standardize = c(to_standardize,pattern_columns,target)
# long_to_standardize = c(to_standardize,pattern_columns,'time','value')
# cross_normalization = preProcess(df_av45[,cross_to_standardize])
# long_normalization = preProcess(df_long[,long_to_standardize])
# df_av45[,cross_to_standardize] = predict(cross_normalization, df_av45[,cross_to_standardize])

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
base_form = paste(target,"~ diag_prior*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
nopattern_form = paste(target,"~ diag_prior*APOE4_BIN + CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
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


save.printout(paste('R/dod_output/',target,'_fm_base_summary','.txt',sep=''),fm_base.summary)
save.printout(paste('R/dod_output/',target,'_fm_nopattern_summary','.txt',sep=''),fm_nopattern.summary)
save.printout(paste('R/dod_output/',target,'_fm_full_summary','.txt',sep=''),fm_full.summary)
save.printout(paste('R/dod_output/',target,'_fm_base_fit','.txt',sep=''),fm_base.fit)
save.printout(paste('R/dod_output/',target,'_fm_nopattern_fit','.txt',sep=''),fm_nopattern.fit)
save.printout(paste('R/dod_output/',target,'_fm_full_fit','.txt',sep=''),fm_full.fit)
save.printout(paste('R/dod_output/',target,'_fm_base_anova','.txt',sep=''),fm_base.anova)
save.printout(paste('R/dod_output/',target,'_fm_nopattern_anova','.txt',sep=''),fm_nopattern.anova)
save.printout(paste('R/dod_output/',target,'_fm_full_anova','.txt',sep=''),fm_full.anova)
save(rfe.base,file=paste('R/dod_output/',target,'_rfe_base_obj',sep=''))
save(rfe.nopattern,file=paste('R/dod_output/',target,'_rfe_nopattern_obj',sep=''))
save(rfe.full,file=paste('R/dod_output/',target,'_rfe_full_obj',sep=''))

fm_base.plotfn = function() {par(mfrow=c(2,2));plot(fm_base);title("Base Model", outer=T, line=-2);}
fm_nopattern.plotfn = function() {par(mfrow=c(2,2));plot(fm_nopattern);title("No Pattern Model", outer=T, line=-2);}
fm_full.plotfn = function() {par(mfrow=c(2,2));plot(fm_full);title("Full Model", outer=T, line=-2);}
save.plot(paste('R/dod_output/',target,'_fm_base_lmplot.pdf',sep=''), fm_base.plotfn)
save.plot(paste('R/dod_output/',target,'_fm_nopattern_lmplot.pdf',sep=''), fm_nopattern.plotfn)
save.plot(paste('R/dod_output/',target,'_fm_full_lmplot.pdf',sep=''), fm_full.plotfn)

save.plot(paste('R/dod_output/',target,'_fm_base_avplot.pdf',sep=''), function() {avPlots(fm_base, ask=FALSE)})
save.plot(paste('R/dod_output/',target,'_fm_nopattern_avplot.pdf',sep=''), function() {avPlots(fm_nopattern, ask=FALSE)})
save.plot(paste('R/dod_output/',target,'_fm_full_avplot.pdf',sep=''), function() {avPlots(fm_full, ask=FALSE)})

fm_base.summary
fm_nopattern.summary
fm_full.summary

fm_base.fit
fm_nopattern.fit
fm_full.fit


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

