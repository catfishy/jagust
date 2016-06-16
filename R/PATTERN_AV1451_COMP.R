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
              'Diag.AV1451','Diag.AV45','positive_prior','positive_post',
              'AV45_NONTP_wcereb_BIN1.11')
to_standardize = c('Age.AV45','Edu..Yrs.')
demog_columns = c('RID','APOE4_BIN','Diag.AV1451','Age.AV1451','Gender','Edu..Yrs.')
braak_columns = c('AV1451_Braak1_CerebGray_BL',
                  'AV1451_Braak2_CerebGray_BL',
                  'AV1451_Braak3_CerebGray_BL',
                  'AV1451_Braak4_CerebGray_BL',
                  'AV1451_Braak5_CerebGray_BL',
                  'AV1451_Braak6_CerebGray_BL')
braak_cum_columns = c('AV1451_Braak1_CerebGray_BL',
                      'AV1451_Braak2_CerebGray_Cum_BL',
                      'AV1451_Braak3_CerebGray_Cum_BL',
                      'AV1451_Braak4_CerebGray_Cum_BL',
                      'AV1451_Braak5_CerebGray_Cum_BL',
                      'AV1451_Braak6_CerebGray_Cum_BL')


valid_diags = c('N','EMCI','LMCI','AD')

# IMPORT
df_av1451 = read.csv('nsfa/av1451uni_pattern_dataset.csv')
df_av1451 = df_av1451[which(df_av1451$Diag.AV1451 %in% valid_diags),]
df_av1451$Diag.AV1451 = factor(df_av1451$Diag.AV1451, levels=valid_diags)
for (i in names(df_av1451)){
  if (i %in% to_factor){
    df_av1451[,eval(i)] = as.factor(as.character(df_av1451[,eval(i)]))
  }
}
pattern_columns = Filter(function(i) {startsWith(i,'NSFA_')}, names(df_av1451))

# find correlation with braak stages
for (pcol in pattern_columns) {
  for (bcol in braak_columns) {
    r2 = summary(lm(paste(pcol,'~',bcol),df_av1451))$adj.r.squared
    if (r2 > 0.1) {
      print(paste(pcol,bcol,r2))
    }
    
  }
}

for (pcol in pattern_columns) {
  for (bcol in braak_cum_columns) {
    r2 = summary(lm(paste(pcol,'~',bcol),df_av1451))$adj.r.squared
    if (r2 > 0.1) {
      print(paste(pcol,bcol,r2))
    }
  }
}

pcol = 'NSFA_1'
bcol = 'AV1451_Braak3_CerebGray_BL'
r2 = summary(lm(paste(pcol,'~',bcol),df_av1451))$adj.r.squared
p = ggplot(df_av1451,aes_string(x=bcol,y=pcol)) +
  geom_point(aes_string(color='Diag.AV1451')) + 
  geom_smooth(method='lm') + 
  annotate("text",x=1.0,y=0,size=18,label=paste("R2:",round(r2,3)))
print(p)


pvc_r2 = summary(lm(NSFA_6 ~ CORTICAL_SUMMARY_prior, df_av1451))$adj.r.squared
nonpvc_r2 = summary(lm(NSFA_6 ~ AV1451_NONTP_wcereb, df_av1451))$adj.r.squared

p1 = ggplot(df_av1451,aes_string(x='AV1451_NONTP_wcereb',y='NSFA_6')) +
  geom_point(aes_string(color='Diag.AV1451')) + 
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

p2 = ggplot(df_av1451,aes_string(x='CORTICAL_SUMMARY_prior',y='NSFA_6')) +
  geom_point(aes_string(color='Diag.AV1451')) + 
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






