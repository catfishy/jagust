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

valid_diags = c('N','SMC','EMCI','LMCI','AD')


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
ancova_factors = 'APOE4_BIN + Gender + Age.AV45 + Edu..Yrs. + Diag.AV45'
pattern_formula = lapply(pattern_columns, function (x) lm.createformula(x,ancova_factors))
naive_formula = lapply(naive_columns, function (x) lm.createformula(x,ancova_factors))
pvalues = list()

for (i in 1:NROW(pattern_formula)) {
  pattern_str = pattern_columns[i]
  formula = as.formula(paste(pattern_formula[i]))
  cur.lm = lm(formula,df_av45)
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
