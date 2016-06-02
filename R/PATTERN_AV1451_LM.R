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
              'Diag.AV45','Diag.AV1451','positive_prior','positive_post',
              'AV45_NONTP_wcereb_BIN1.11','AV45_NONTP_2_wcereb_BIN1.11',
              'AV45_NONTP_3_wcereb_BIN1.11')
demog_columns = c('RID','APOE4_BIN','Diag.AV1451','Age.AV1451','Gender','Edu..Yrs.')


#target = "UW_EF_AV1451_1"
#target = "UW_MEM_AV1451_1"
#target = "ADAS_AV1451_1"
target = "AVLT_AV1451_1"

output_folder = 'R/output_av1451/'

valid_diags = c('N','SMC','EMCI','LMCI','AD')
#valid_diags = c('N','SMC','EMCI','LMCI')
#valid_diags = c('N','SMC')
#valid_diags = c('EMCI')
#valid_diags = c('LMCI')

# IMPORT
df_av1451 = read.csv('nsfa/av1451_pattern_dataset.csv')
pattern_columns = Filter(isPatternColumn,names(df_av1451))

df_av1451 = df_av1451[which(df_av1451$Diag.AV1451 %in% valid_diags),]
for (i in names(df_av1451)){
  if (i %in% to_factor){
    df_av1451[,eval(i)] = as.factor(as.character(df_av1451[,eval(i)]))
  }
}

# look at histograms
for (pcol in pattern_columns) {
  p = ggplot(df_av1451, aes_string(pcol)) + geom_histogram(binwidth=0.1)
  print(p)
}

# make crossx response normal
non.na = complete.cases(df_av1451[,c(demog_columns,target)])
df_av1451 = df_av1451[non.na,]
df_av1451[,eval(target)] = Gaussianize(df_av1451[,eval(target)], type='hh', method='MLE', return.u=TRUE)


run.rfe = function(form, var.response, dataset) {
  x = as.data.frame(model.matrix(as.formula(form),dataset))[,-1]
  nzv_cols = nearZeroVar(x)
  if (length(nzv_cols) > 0) {
    x = x[, -nzv_cols]
  }
  corr_cols = findCorrelation(cor(x),.9)
  if (length(corr_cols) > 0) {
    x = x[, -corr_cols]
  }

  colnames(x) = lapply(colnames(x), make.names)
  rownames(x) = NULL
  y = as.numeric(df_av1451[,var.response])
  
  ctrl = rfeControl(functions = lmFuncs, 
                    method = "repeatedcv", 
                    number = 10,
                    repeats = 5,
                    rerank = TRUE,
                    verbose = FALSE)
  set.seed(1337)
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

diag_str = 'Diag.AV1451*APOE4_BIN +'
#diag_str = ''

base_form = paste(target,"~",diag_str,"Age.AV1451 + Gender + Edu..Yrs. + APOE4_BIN")
braak_form = paste(target,"~",diag_str,"Age.AV1451 + Gender + Edu..Yrs. + APOE4_BIN + AV1451_PVC_Braak12_CerebGray_BL + AV1451_PVC_Braak34_CerebGray_BL + AV1451_PVC_Braak56_CerebGray_BL")
pattern_form = paste(target,"~",diag_str,"Age.AV1451 + Gender + Edu..Yrs. + APOE4_BIN +",paste(pattern_columns,collapse=' + '))

rfe.base = run.rfe(base_form, target, df_av1451)
rfe.braak = run.rfe(braak_form, target, df_av1451)
rfe.pattern = run.rfe(pattern_form, target, df_av1451)

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
