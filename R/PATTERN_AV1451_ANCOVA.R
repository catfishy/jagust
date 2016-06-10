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
to_standardize = c('Age.AV45','Edu..Yrs.')
demog_columns = c('RID','APOE4_BIN','Diag.AV45','Age.AV45','Gender','Edu..Yrs.')
braak_columns = c('AV1451_Braak1_CerebGray_BL',
                  'AV1451_Braak2_CerebGray_BL',
                  'AV1451_Braak3_CerebGray_BL',
                  'AV1451_Braak4_CerebGray_BL',
                  'AV1451_Braak5_CerebGray_BL',
                  'AV1451_Braak6_CerebGray_BL')

output_folder = 'R/output_av1451_ancova/'

valid_diags = c('N','EMCI','LMCI')


# IMPORT
df_av1451 = read.csv('nsfa/av1451_pattern_dataset.csv')
df_av1451 = df_av1451[which(df_av1451$Diag.AV45 %in% valid_diags),]
non.na = complete.cases(df_av1451[,c(demog_columns,av45_columns)])
df_av1451 = df_av1451[non.na,]
for (i in names(df_av1451)){
  if (i %in% to_factor){
    df_av1451[,eval(i)] = as.factor(as.character(df_av1451[,eval(i)]))
  }
}
pattern_columns = Filter(isPatternColumn,names(df_av1451))
naive_columns = Filter(isNaiveColumn,names(df_av1451))

# standardize predictors
cross_to_standardize = c(to_standardize,pattern_columns,naive_columns,braak_columns)
cross_normalization = preProcess(df_av1451[,cross_to_standardize])
df_av1451[,cross_to_standardize] = predict(cross_normalization, df_av1451[,cross_to_standardize])

# make response normal
#df_av1451[,eval(target)] = Gaussianize(df_av1451[,eval(target)], type='hh', method='MLE', return.u=TRUE)



ancova_factors = 'APOE4_BIN + Gender + Age.AV1451 + Edu..Yrs. + Diag.AV1451'
pattern_formula = lapply(pattern_columns, function (x) lm.createformula(x,ancova_factors))
naive_formula = lapply(naive_columns, function (x) lm.createformula(x,ancova_factors))
braak_formula = lapply(braak_columns, function (x) lm.createformula(x,ancova_factors))
pvalues = list()

for (i in 1:NROW(pattern_formula)) {
  pattern_str = pattern_columns[i]
  formula = as.formula(paste(pattern_formula[i]))
  cur.lm = lm(formula,df_av1451)
  cur.lm.summary = summary(cur.lm)
  cur.lm.coef = cur.lm.summary$coefficients
  fstat = cur.lm.summary$fstatistic
  cur.lm.pvalue = pf(fstat[1],fstat[2],fstat[3],lower.tail=F)
  pvalues[pattern_str] = cur.lm.pvalue
  cur.lm.anova = Anova(cur.lm,type='III')
  save.printout(paste(output_folder,pattern_str,'_ancova_anova','.txt',sep=''),cur.lm.anova)
  save.printout(paste(output_folder,pattern_str,'_ancova_summary','.txt',sep=''),cur.lm.summary)
  lm.plotfn = function() {par(mfrow=c(2,2));plot(cur.lm);title(pattern_str, outer=T, line=-2);}
  save.plot(paste(output_folder,pattern_str,'_lmplot.pdf',sep=''), lm.plotfn)
  save.plot(paste(output_folder,pattern_str,'_avplot.pdf',sep=''), function() {avPlots(cur.lm, ask=FALSE)})
}

for (i in 1:NROW(braak_formula)) {
  pattern_str = braak_columns[i]
  formula = as.formula(paste(braak_formula[i]))
  cur.lm = lm(formula,df_av1451)
  cur.lm.summary = summary(cur.lm)
  cur.lm.coef = cur.lm.summary$coefficients
  fstat = cur.lm.summary$fstatistic
  cur.lm.pvalue = pf(fstat[1],fstat[2],fstat[3],lower.tail=F)
  pvalues[pattern_str] = cur.lm.pvalue
  cur.lm.anova = Anova(cur.lm,type='III')
  save.printout(paste(output_folder,pattern_str,'_ancova_anova','.txt',sep=''),cur.lm.anova)
  save.printout(paste(output_folder,pattern_str,'_ancova_summary','.txt',sep=''),cur.lm.summary)
  lm.plotfn = function() {par(mfrow=c(2,2));plot(cur.lm);title(pattern_str, outer=T, line=-2);}
  save.plot(paste(output_folder,pattern_str,'_lmplot.pdf',sep=''), lm.plotfn)
  save.plot(paste(output_folder,pattern_str,'_avplot.pdf',sep=''), function() {avPlots(cur.lm, ask=FALSE)})
}

pvalues.corrected = p.adjust(pvalues,method='bonferroni')
pvalues.sig = pvalues.corrected[pvalues.corrected < 0.05]
