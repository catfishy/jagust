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
              'Diag.AV45','positive_prior','positive_post',
              'AV45_NONTP_wcereb_BIN1.11')
to_standardize = c('CORTICAL_SUMMARY_prior','Age.AV45','Edu..Yrs.')
demog_columns = c('RID','APOE4_BIN','Diag.AV45','Age.AV45','Gender','Edu..Yrs.')
av45_columns = c('CORTICAL_SUMMARY_prior','AV45_NONTP_wcereb')

output_folder = 'R/output_av45_ancova/'

# valid_diags = c('N','EMCI','LMCI','AD')
valid_diags = c('EMCI')

# ancova_factors = 'Age.AV45 + Diag.AV1451 + APOE4_BIN'
# ancova_factors = 'Diag.AV45'
# ancova_factors = 'Diag.AV45:AV45_NONTP_wcereb_BIN1.11'
# ancova_factors = 'AV45_NONTP_wcereb_BIN1.11'
# ancova_factors = 'Age.AV45'
# ancova_factors = 'APOE4_BIN'
# ancova_factors = 'Gender'
# ancova_factors = 'Edu..Yrs.'
ancova_factors = 'AVLT_slope_postAV45'
# ancova_factors = 'ADASslope_postAV45'
# ancova_factors = 'MMSEslope_postAV45'

# IMPORT
df_av45 = read.csv('nsfa/av45_pattern_dataset.csv')
df_av45 = df_av45[which(df_av45$Diag.AV45 %in% valid_diags),]
non.na = complete.cases(df_av45[,c(demog_columns,av45_columns)])
df_av45 = df_av45[non.na,]
for (i in names(df_av45)){
  if (i %in% to_factor){
    df_av45[,eval(i)] = as.factor(as.character(df_av45[,eval(i)]))
  }
}
pattern_columns = Filter(isPatternColumn,names(df_av45))
naive_columns = Filter(isNaiveColumn,names(df_av45))

# standardize predictors
cross_to_standardize = c(to_standardize,pattern_columns,naive_columns)
cross_normalization = preProcess(df_av45[,cross_to_standardize])
df_av45[,cross_to_standardize] = predict(cross_normalization, df_av45[,cross_to_standardize])

# make response normal
#df_av45[,eval(target)] = Gaussianize(df_av45[,eval(target)], type='hh', method='MLE', return.u=TRUE)


#all.addons = lapply(pattern_columns,lm.addvar)
#naive.addons = lapply(naive_columns,lm.addvar)
all.addons = paste('+',paste(pattern_columns,collapse=' + '))
naive.addons = paste('+',paste(naive_columns,collapse=' + '))
pattern_formula = lapply(pattern_columns, function (x) lm.createformula(x,ancova_factors))
av45_formula = lapply(av45_columns, function (x) lm.createformula(x,ancova_factors))
naive_formula = lapply(naive_columns, function (x) lm.createformula(x,ancova_factors))
pvalues = list()
lmformula = list()

all_formula = c(pattern_formula,av45_formula)
all_columns = c(pattern_columns,av45_columns)
for (i in 1:NROW(all_formula)) {
  pattern_str = all_columns[i]
  formula = as.formula(paste(all_formula[i]))
  print(formula)
  cur.lm = lm(formula,df_av45)
  cur.lm.summary = summary(cur.lm)
  cur.lm.coef = cur.lm.summary$coefficients
  fstat = cur.lm.summary$fstatistic
  cur.lm.pvalue = pf(fstat[1],fstat[2],fstat[3],lower.tail=F)
  pvalues[pattern_str] = cur.lm.pvalue
  lmformula[pattern_str] = paste(pattern_formula[i])
  cur.lm.anova = Anova(cur.lm,type='III')
#   save.printout(paste(output_folder,pattern_str,'_',ancova_factors,'_ancova_anova','.txt',sep=''),cur.lm.anova)
#   save.printout(paste(output_folder,pattern_str,'_',ancova_factors,'_ancova_summary','.txt',sep=''),cur.lm.summary)
#   lm.plotfn = function() {par(mfrow=c(2,2));plot(cur.lm);title(pattern_str, outer=T, line=-2);}
#   save.plot(paste(output_folder,pattern_str,'_',ancova_factors,'_lmplot.pdf',sep=''), lm.plotfn)
#   save.plot(paste(output_folder,pattern_str,'_',ancova_factors,'_avplot.pdf',sep=''), function() {avPlots(cur.lm, ask=FALSE)})
}

pvalues.corrected = p.adjust(pvalues,method='bonferroni')
pvalues.sig = pvalues.corrected[pvalues.corrected < 0.05]
pvalues.corrected



to_graph = c('NSFA_6','NSFA_8','NSFA_0')
# Diagnosis graphs
for (pcol in to_graph) {
  p = ggplot(df_av45, aes_string('Diag.AV45',pcol)) +
    geom_violin(trim=TRUE,aes(fill=Diag.AV45),show.legend=FALSE) +
    coord_flip() + 
    theme(axis.title.x=element_text(size=20),
          axis.title.y=element_text(size=20),
          axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14),
          legend.title=element_text(size=14)) + 
    ylab(paste(pcol,'Score')) +
    xlab('Diagnosis')
  print(p)
}

p2 = ggplot(df_av45, aes_string('Diag.AV45','CORTICAL_SUMMARY_prior')) +
  geom_violin(trim=TRUE,aes(fill=Diag.AV45),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab('CS SUVR (whole cereb. ref., PVC)') +
  xlab('Diagnosis')
print(p2)

# Age graphs
for (pcol in to_graph) {
  p = ggplot(df_av45, aes_string(x='Age.AV45',y=pcol)) +
    geom_point(aes_string(color='Diag.AV45'),show.legend=FALSE) + 
    geom_smooth(method='lm') + 
    coord_flip() + 
    theme(axis.title.x=element_text(size=20),
          axis.title.y=element_text(size=20),
          axis.text.x=element_text(size=14), 
          axis.text.y=element_text(size=14)) + 
    ylab(paste(pcol,'Score')) +
    xlab('Age')
  print(p)
}

p2 = ggplot(df_av45, aes_string(x='Age.AV45',y='CORTICAL_SUMMARY_prior')) +
  geom_point(aes_string(color='Diag.AV45'),show.legend=FALSE) + 
  geom_smooth(method='lm') +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14)) + 
  ylab('CS SUVR (whole cereb. ref., PVC)') +
  xlab('Age')
print(p2)

# APOE4 graphs
for (pcol in to_graph) {
  p = ggplot(df_av45, aes_string('APOE4_BIN',pcol)) +
    geom_violin(trim=TRUE,aes(fill=APOE4_BIN),show.legend=FALSE) +
    coord_flip() + 
    theme(axis.title.x=element_text(size=20),
          axis.title.y=element_text(size=20),
          axis.text.x=element_text(size=14), 
          axis.text.y=element_text(size=14),
          legend.title=element_text(size=14)) + 
    ylab(paste(pcol,'Score')) +
    xlab(expression(paste('APOE ',epsilon,'4',collapse='')))
  print(p)
}

p2 = ggplot(df_av45, aes_string('APOE4_BIN','CORTICAL_SUMMARY_prior')) +
  geom_violin(trim=TRUE,aes(fill=APOE4_BIN),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab('CS SUVR (whole cereb. ref., PVC)') +
  xlab(expression(paste('APOE ',epsilon,'4',collapse='')))
print(p2)

# AV45+ graphs
for (pcol in to_graph) {
  p = ggplot(df_av45, aes_string('AV45_NONTP_wcereb_BIN1.11',pcol)) +
    geom_violin(trim=TRUE,aes(fill=AV45_NONTP_wcereb_BIN1.11),show.legend=FALSE) +
    coord_flip() + 
    theme(axis.title.x=element_text(size=20),
          axis.title.y=element_text(size=20),
          axis.text.x=element_text(size=14), 
          axis.text.y=element_text(size=14),
          legend.title=element_text(size=14)) + 
    ylab(paste(pcol,'Score')) +
    xlab('AV45 SUVR Pos.')
  print(p)
}

p2 = ggplot(df_av45, aes_string('AV45_NONTP_wcereb_BIN1.11','CORTICAL_SUMMARY_prior')) +
  geom_violin(trim=TRUE,aes(fill=AV45_NONTP_wcereb_BIN1.11),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab('CS SUVR (whole cereb. ref., PVC)') +
  xlab('AV45 SUVR Pos.')
print(p2)

