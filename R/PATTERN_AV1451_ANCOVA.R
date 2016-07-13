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
              'Diag.AV45','Diag.AV1451','positive_prior','positive_post',
              'AV45_NONTP_wcereb_BIN1.11')
to_standardize = c('Age.AV45','Age.AV1451','Edu..Yrs.')
demog_columns = c('RID','APOE4_BIN','Diag.AV45','Diag.AV1451',
                  'Age.AV45','Age.AV1451','Gender','Edu..Yrs.')
# braak_columns = c('AV1451_Braak1_CerebGray_BL',
#                   'AV1451_Braak2_CerebGray_BL',
#                   'AV1451_Braak3_CerebGray_BL',
#                   'AV1451_Braak4_CerebGray_BL',
#                   'AV1451_Braak5_CerebGray_BL',
#                   'AV1451_Braak6_CerebGray_BL',
#                   'AV1451_Braak12_CerebGray_BL',
#                   'AV1451_Braak34_CerebGray_BL',
#                   'AV1451_Braak56_CerebGray_BL')
braak_columns = c('AV1451_Braak12_CerebGray_BL',
                  'AV1451_Braak34_CerebGray_BL',
                  'AV1451_Braak56_CerebGray_BL')


valid_diags = c('N','SMC','EMCI','LMCI','AD')


# IMPORT
output_folder = 'R/output_av1451_ancova/'
df_av1451 = read.csv('nsfa/av1451_pattern_dataset.csv')

# ancova_factors = 'Age.AV1451 + Diag.AV1451 + APOE4_BIN'
ancova_factors = 'Diag.AV1451'
# ancova_factors = 'Age.AV1451'
# ancova_factors = 'APOE4_BIN'
# ancova_factors = 'Gender'
# ancova_factors = 'Edu..Yrs.'

# FILTER
non.na = complete.cases(df_av1451[,c(demog_columns,braak_columns)])
df_av1451 = df_av1451[non.na,]
df_av1451 = df_av1451[which(df_av1451$Diag.AV1451 %in% valid_diags),]
for (i in names(df_av1451)){
  if (i %in% to_factor){
    df_av1451[,eval(i)] = as.factor(as.character(df_av1451[,eval(i)]))
  }
}
df_av1451$Diag.AV1451 = factor(df_av1451$Diag.AV1451, levels=valid_diags)
pattern_columns = Filter(isPatternColumn,names(df_av1451))
naive_columns = Filter(isNaiveColumn,names(df_av1451))

# # standardize predictors
# cross_to_standardize = c(to_standardize,pattern_columns,naive_columns,braak_columns)
# cross_normalization = preProcess(df_av1451[,cross_to_standardize])
# df_av1451[,cross_to_standardize] = predict(cross_normalization, df_av1451[,cross_to_standardize])

# make response normal
#df_av1451[,eval(target)] = Gaussianize(df_av1451[,eval(target)], type='hh', method='MLE', return.u=TRUE)


pattern_formula = lapply(pattern_columns, function (x) lm.createformula(x,ancova_factors))
naive_formula = lapply(naive_columns, function (x) lm.createformula(x,ancova_factors))
braak_formula = lapply(braak_columns, function (x) lm.createformula(x,ancova_factors))
pvalues = list()
lmformula = list()

for (i in 1:NROW(pattern_formula)) {
  pattern_str = pattern_columns[i]
  formula = as.formula(paste(pattern_formula[i]))
  cur.lm = lm(formula,df_av1451)
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

for (i in 1:NROW(braak_formula)) {
  pattern_str = braak_columns[i]
  formula = as.formula(paste(braak_formula[i]))
  cur.lm = lm(formula,df_av1451)
  cur.lm.summary = summary(cur.lm)
  cur.lm.coef = cur.lm.summary$coefficients
  fstat = cur.lm.summary$fstatistic
  cur.lm.pvalue = pf(fstat[1],fstat[2],fstat[3],lower.tail=F)
  pvalues[pattern_str] = cur.lm.pvalue
  lmformula[pattern_str] = paste(braak_formula[i])
  cur.lm.anova = Anova(cur.lm,type='III')
#   save.printout(paste(output_folder,pattern_str,'_',ancova_factors,'_ancova_anova','.txt',sep=''),cur.lm.anova)
#   save.printout(paste(output_folder,pattern_str,'_',ancova_factors,'_ancova_summary','.txt',sep=''),cur.lm.summary)
#   lm.plotfn = function() {par(mfrow=c(2,2));plot(cur.lm);title(pattern_str, outer=T, line=-2);}
#   save.plot(paste(output_folder,pattern_str,'_',ancova_factors,'_lmplot.pdf',sep=''), lm.plotfn)
#   save.plot(paste(output_folder,pattern_str,'_',ancova_factors,'_avplot.pdf',sep=''), function() {avPlots(cur.lm, ask=FALSE)})
}

pvalues.corrected = p.adjust(pvalues,method='bonferroni')
pvalues.sig = pvalues.corrected[pvalues.corrected < 0.05]

pcol = 'NSFA_0'

# Diagnosis graphs
p1 = ggplot(df_av1451, aes_string('Diag.AV1451',pcol)) +
  geom_violin(trim=TRUE,aes(fill=Diag.AV1451),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab(paste(pcol,'Score')) +
  xlab('Diagnosis')

p2 = ggplot(df_av1451, aes_string('Diag.AV1451','AV1451_Braak34_CerebGray_BL')) +
  geom_violin(trim=TRUE,aes(fill=Diag.AV1451),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab('Braak III/IV SUVR (cereb. gray ref.)') +
  xlab('Diagnosis')

p3 = ggplot(df_av1451, aes_string('Diag.AV1451','AV1451_Braak56_CerebGray_BL')) +
  geom_violin(trim=TRUE,aes(fill=Diag.AV1451),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab('Braak V/VI SUVR (cereb. gray ref.)') +
  xlab('Diagnosis')

print(p1)
print(p2)
print(p3)

# Age graphs
p1 = ggplot(df_av1451, aes_string(x='Age.AV1451',y=pcol)) +
  geom_point(aes_string(color='Diag.AV1451'),show.legend=FALSE) + 
  geom_smooth(method='lm') + 
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14)) + 
  ylab(paste(pcol,'Score')) +
  xlab('Age')
p2 = ggplot(df_av1451, aes_string(x='Age.AV1451',y='NSFA_14')) +
  geom_point(aes_string(color='Diag.AV1451'),show.legend=FALSE) + 
  geom_smooth(method='lm') +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14)) + 
  ylab('NSFA_14 Score') +
  xlab('Age')
p3 = ggplot(df_av1451, aes_string(x='Age.AV1451',y='AV1451_Braak34_CerebGray_BL')) +
  geom_point(aes_string(color='Diag.AV1451'),show.legend=FALSE) + 
  geom_smooth(method='lm') + 
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14)) + 
    ylab('Braak III/IV SUVR (cereb. gray ref.)') +
    xlab('Age')
p4 = ggplot(df_av1451, aes_string(x='Age.AV1451',y='AV1451_Braak56_CerebGray_BL')) +
  geom_point(aes_string(color='Diag.AV1451'),show.legend=FALSE) + 
  geom_smooth(method='lm') + 
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14)) + 
  ylab('Braak V/VI SUVR (cereb. gray ref.)') +
  xlab('Age')

print(p1)
print(p2)
print(p3)
print(p4)
  




# APOE4 graphs
pcol = 'NSFA_2'
p1 = ggplot(df_av1451, aes_string('APOE4_BIN',pcol)) +
  geom_violin(trim=TRUE,aes(fill=APOE4_BIN),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab(paste(pcol,'Score')) +
  xlab(expression(paste('APOE ',epsilon,'4',collapse='')))

p2 = ggplot(df_av1451, aes_string('APOE4_BIN','AV1451_Braak34_CerebGray_BL')) +
  geom_violin(trim=TRUE,aes(fill=APOE4_BIN),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab('Braak III/IV SUVR (cereb. gray ref.)') +
  xlab(expression(paste('APOE ',epsilon,'4',collapse='')))

p3 = ggplot(df_av1451, aes_string('APOE4_BIN','AV1451_Braak56_CerebGray_BL')) +
  geom_violin(trim=TRUE,aes(fill=APOE4_BIN),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14), 
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab('Braak V/VI SUVR (cereb. gray ref.)') +
  xlab(expression(paste('APOE ',epsilon,'4',collapse='')))

print(p1)
print(p2)
print(p3)




