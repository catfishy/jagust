library(LambertW)
library(lme4)


source('R/LM_FUNCS.R')

pattern_prefix = 'NSFA_'
to_factor = c('RID','ad_prior','ad_post','positive_prior','positive_post','diag_prior','diag_post','APOE4_BIN','APOE2_BIN','Gender','Diag.AV45_long','positive_prior','positive_post')
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
#target = "CSF_ABETA_closest_AV45_1"
#target = "CSF_TAU_closest_AV45_1"
#target = "CSF_PTAU_closest_AV45_1"

#valid_diags = c('N','SMC','EMCI','LMCI','AD')
#valid_diags = c('N','SMC','EMCI','LMCI')
#valid_diags = c('N','SMC')
valid_diags = c('EMCI')
#valid_diags = c('LMCI')

df_av45 = read.csv('nsfa/pattern_dataset.csv')
pattern_columns = Filter(isPatternColumn,names(df_av45))
df_av45 = df_av45[which(df_av45$diag_prior %in% valid_diags),]
for (i in names(df_av45)){
  if (i %in% to_factor){
    df_av45[,eval(i)] = as.factor(as.character(df_av45[,eval(i)]))
  }
}

# Make response normal (TODO: Z-scores)
non.na = complete.cases(df_av45[,c(pattern_columns,demog_columns,av45_columns,target)])
df_av45 = df_av45[non.na,]
df_av45[,eval(target)] = Gaussianize(df_av45[,eval(target)], type='hh', method='MLE', return.u=TRUE)

# Run model
if (length(valid_diags) > 1) {
  diag_str = 'diag_prior*APOE4_BIN +'
} else {
  diag_str = ''
}

#nopattern_form = paste(target,"~",diag_str,"CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
nopattern_form = paste(target,"~",diag_str,"CORTICAL_SUMMARY_prior*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
onlypattern_form = paste(target,"~",diag_str,"NSFA_6*APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")

fm_nopattern = lm(as.formula(nopattern_form),data=df_av45)
fm_onlypattern = lm(as.formula(onlypattern_form),data=df_av45)

fm_nopattern.summary = summary(fm_nopattern)
fm_onlypattern.summary = summary(fm_onlypattern)

fm_nopattern.summary$adj.r.squared
fm_onlypattern.summary$adj.r.squared


