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
              'Diag.AV45_long','positive_prior','positive_post',
              'AV45_NONTP_wcereb_BIN1.11')
to_standardize = c('CORTICAL_SUMMARY_prior','Age.AV45','Edu..Yrs.')
demog_columns = c('RID','APOE4_BIN','Diag.AV45','Age.AV45','Gender','Edu..Yrs.')
av45_columns = c('CORTICAL_SUMMARY_prior')

output_folder = 'R/output_av45_ancova/'

valid_diags = c('N','EMCI','LMCI','AD')

# IMPORT
df_av45 = read.csv('nsfa/av45_pattern_dataset.csv')
df_av45 = df_av45[which(df_av45$Diag.AV45 %in% valid_diags),]
df_av45$Diag.AV45 = factor(df_av45$Diag.AV45, levels=valid_diags)
#non.na = complete.cases(df_av45[,c(demog_columns,av45_columns)])
#df_av45 = df_av45[non.na,]
for (i in names(df_av45)){
  if (i %in% to_factor){
    df_av45[,eval(i)] = as.factor(as.character(df_av45[,eval(i)]))
  }
}


# plot correlation with cortical summary
pvc_r2 = summary(lm(NSFA_6 ~ CORTICAL_SUMMARY_prior, df_av45))$adj.r.squared
nonpvc_r2 = summary(lm(NSFA_6 ~ AV45_NONTP_wcereb, df_av45))$adj.r.squared

p1 = ggplot(df_av45,aes_string(x='AV45_NONTP_wcereb',y='NSFA_6')) +
  geom_point(aes_string(color='Diag.AV45')) + 
  geom_smooth(method='lm') + 
  ylim(-2.5,2.5) +
  theme(plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14),
        axis.text.x=element_text(face='bold', size=14)) +
    ggtitle('Cortical Summary SUVR (whole cereb. ref., Non-PVC) vs. NSFA_6 Factor Score') +
    xlab('Cortical Summary SUVR (whole cereb. ref., Non-PVC)') +
    ylab('NSFA_6 Factor Score') +
    annotate("text",x=1.5,y=-1.5,size=18,label=paste("R2:",round(nonpvc_r2,3)))
print(p1)

p2 = ggplot(df_av45,aes_string(x='CORTICAL_SUMMARY_prior',y='NSFA_6')) +
  geom_point(aes_string(color='Diag.AV45')) + 
  geom_smooth(method='lm') + 
  ylim(-2.5,2.5) +
  theme(plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14),
        axis.text.x=element_text(face='bold', size=14)) +
  ggtitle('Cortical Summary SUVR (whole cereb. ref., PVC) vs. NSFA_6 Factor Score') +
  xlab('Cortical Summary SUVR (whole cereb. ref., PVC)') +
  ylab('NSFA_6 Factor Score') + 
  annotate("text",x=2.5,y=-1.5,size=18,label=paste("R2:",round(pvc_r2,3)))
print(p2)




pattern_columns = Filter(function(i) {startsWith(i,'NSFA_')}, names(df_av45))
scan2_columns = Filter(function(i) {startsWith(i,'SCAN2_NSFA_')}, names(df_av45))
scan3_columns = Filter(function(i) {startsWith(i,'SCAN3_NSFA_')}, names(df_av45))

# pattern_columns = Filter(isPatternColumn,names(df_av45))
# scan2_columns = Filter(isScan2Column,names(df_av45))
# scan3_columns = Filter(isScan3Column,names(df_av45))

# # standardize predictors
# demog_normalization = preProcess(df_av45[,to_standardize])
# pattern_normalization = preProcess(df_av45[,pattern_columns])
# df_av45[,to_standardize] = predict(demog_normalization, df_av45[,to_standardize])
# df_av45[,pattern_columns] = predict(pattern_normalization, df_av45[,pattern_columns])
# # standardize scan 2
# scan2_subset = df_av45[,scan2_columns]
# names(scan2_subset) = gsub("SCAN2_", "", names(scan2_subset))
# scan2_subset = predict(pattern_normalization, scan2_subset)
# names(scan2_subset) = gsub("NSFA_", "SCAN2_NSFA_", names(scan2_subset))
# df_av45[,scan2_columns] = scan2_subset
# # standardize scan 3
# scan3_subset = df_av45[,scan3_columns]
# names(scan3_subset) = gsub("SCAN3_", "", names(scan3_subset))
# scan3_subset = predict(pattern_normalization, scan3_subset)
# names(scan3_subset) = gsub("NSFA_", "SCAN3_NSFA", names(scan3_subset))
# df_av45[,scan3_columns] = scan3_subset

for (pcol in pattern_columns) {
  scan2_col = gsub('NSFA_','SCAN2_NSFA_',pcol)
  scan3_col = gsub('NSFA_','SCAN3_NSFA_',pcol)
  scan2_change = (df_av45[,eval(scan2_col)]-df_av45[,eval(pcol)])/(df_av45$AV45_1_2_Diff)
  scan3_change = (df_av45[,eval(scan3_col)]-df_av45[,eval(pcol)])/(df_av45$AV45_1_3_Diff)
  all_change = scan3_change
  all_change[is.na(all_change)] = scan2_change[is.na(all_change)]
  df_av45[,paste('SCAN2_CHANGE_',pcol,sep='')] = scan2_change
  df_av45[,paste('SCAN3_CHANGE_',pcol,sep='')] = scan3_change
  df_av45[,paste('ALL_CHANGE_',pcol,sep='')] = all_change
}

# plot factor 6 changes
pcol = 'NSFA_6'
p1 = ggplot(df_av45, aes_string(x=pcol, y=paste('ALL_CHANGE_',pcol,sep=''))) +
  geom_point(aes_string(colour ='Diag.AV45')) +
  geom_smooth() +
  theme(plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14),
        axis.text.x=element_text(face='bold', size=14)) + 
  geom_hline(aes(yintercept=0.0)) + 
  ggtitle('Baseline NSFA_6 Score vs. Annualized Change') +
  xlab('Cortical Summary SUVR (whole cereb. ref., PVC)') +
  ylab('Annualized Change')
print(p1)

p2 = ggplot(df_av45, aes_string(x='CORTICAL_SUMMARY_prior', y='CORTICAL_SUMMARY_change')) +
  geom_point(aes_string(colour ='Diag.AV45')) +
  geom_smooth() + 
  theme(plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14),
        axis.text.x=element_text(face='bold', size=14)) +
  geom_hline(aes(yintercept=0.0)) + 
  ggtitle('Baseline Cortical Summary SUVR (whole cereb. ref., PVC) vs. Annualized Change') +
  xlab('Baseline Cortical Summary SUVR (whole cereb. ref., PVC)') +
  ylab('Annualized Change')
print(p2)

p3 = ggplot(df_av45, aes_string(x='AV45_NONTP_wcereb', y='AV45_NONTP_wcereb_Slope')) +
  geom_point(aes_string(colour ='Diag.AV45')) +
  geom_smooth() + 
  theme(plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14),
        axis.text.x=element_text(face='bold', size=14)) +
  geom_hline(aes(yintercept=0.0)) + 
  ggtitle('Baseline Cortical Summary SUVR (whole cereb. ref., Non-PVC) vs. Annualized Change') +
  xlab('Baseline Cortical Summary SUVR (whole cereb. ref., Non-PVC)') +
  ylab('Annualized Change')
print(p3)



# for (pcol in pattern_columns) {
#   scan2_change_col = paste('SCAN2_CHANGE_',pcol,sep='')
#   scan3_change_col = paste('SCAN3_CHANGE_',pcol,sep='')
#   all_change_col = paste('ALL_CHANGE_',pcol,sep='')
#   p = ggplot(df_av45, aes_string(x=pcol, y=all_change_col)) +
#         geom_point(aes_string(x=pcol, colour ='Diag.AV45'), shape=1) +
#         geom_smooth()
#   print(p)
# }


ggplot(df_av45, aes_string(x='CORTICAL_SUMMARY_prior',y='CORTICAL_SUMMARY_change')) + 
  geom_point(aes_string(x='CORTICAL_SUMMARY_prior', colour ='Diag.AV45'), shape=1) + 
  geom_smooth()

