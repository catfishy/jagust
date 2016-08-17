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
library(Jmisc)
library(lars)
library(covTest)

source('R/LM_FUNCS.R')

# CONSTANTS
pattern_prefix = 'NSFA_'
to_factor = c('RID','ad_prior','ad_post','positive_prior','positive_post',
              'diag_prior','diag_post','APOE4_BIN','APOE2_BIN','Gender',
              'Diag.AV45','Diag.AV1451','positive_prior','positive_post',
              'AV45_NONTP_wcereb_BIN1.11','AV45_NONTP_2_wcereb_BIN1.11',
              'AV45_NONTP_3_wcereb_BIN1.11','AV1451_BL_closest_AV45_wcereb_BIN1.11')
to_standardize = c('Age.AV45','Edu..Yrs.')
demog_columns = c('RID','APOE4_BIN','Diag.AV1451','Age.AV1451','Gender','Edu..Yrs.')
diag_columns = c('Diag.AV45','Diag.AV1451')
# braak_columns = c('AV1451_Braak1_CerebGray_BL',
#                   'AV1451_Braak2_CerebGray_BL',
#                   'AV1451_Braak3_CerebGray_BL',
#                   'AV1451_Braak4_CerebGray_BL',
#                   'AV1451_Braak5_CerebGray_BL',
#                   'AV1451_Braak6_CerebGray_BL')
braak_columns = c('AV1451_Braak12_CerebGray_BL',
                  'AV1451_Braak34_CerebGray_BL',
                  'AV1451_Braak56_CerebGray_BL')

# target = "UW_EF_AV1451_1"
# target = "UW_MEM_AV1451_1"
# target = "ADAS_AV1451_1"
# target = "AVLT_AV1451_1"
# target = "AV1451_BL_closest_AV45_wcereb_BIN1.11"
target = 'AV1451_BL_closest_AV45_wcereb'

# target = "ADAS_retroslope_AV1451_BL"
# target = "AVLT_retroslope_AV1451_BL"
# target = "UW_MEM_retroslope_AV1451_BL"
# target = "UW_EF_retroslope_AV1451_BL"
# target = "AV1451_BL_closest_AV45_wcereb_retroSlope"

output_folder = 'R/output_av1451/'

valid_diags = c('N','SMC','EMCI','LMCI','AD')
#valid_diags = c('N','SMC','EMCI','LMCI')
#valid_diags = c('N','SMC')
#valid_diags = c('EMCI')
#valid_diags = c('LMCI')

# IMPORT
df_av1451 = read.csv('nsfa/av1451skull_pattern_dataset.csv')
pattern_columns = Filter(isNaiveColumn,names(df_av1451))
# pattern_columns = Filter(isPatternColumn,names(df_av1451))
non.na = complete.cases(df_av1451[,c(demog_columns,braak_columns,target)])
df_av1451 = df_av1451[non.na,]
for (i in names(df_av1451)){
  if (i %in% to_factor){
    df_av1451[,eval(i)] = as.factor(as.character(df_av1451[,eval(i)]))
  }
}

# remove target outliers
# target.mean = mean(df_av1451[,target])
# target.sd = sd(df_av1451[,target])
# df_av1451 = df_av1451[df_av1451[,target] <= target.mean+target.sd*5,]
# df_av1451 = df_av1451[df_av1451[,target] >= target.mean-target.sd*5,]

# Filter by diag
df_av1451 = df_av1451[which(df_av1451$Diag.AV1451 %in% valid_diags),]

# standardize predictors
# cross_to_standardize = c(to_standardize,pattern_columns,naive_columns,braak_columns,target)
# cross_normalization = preProcess(df_av1451[,cross_to_standardize])
# df_av1451[,cross_to_standardize] = predict(cross_normalization, df_av1451[,cross_to_standardize])
# df_av1451[,target] = scale(df_av1451[,target], center=TRUE, scale=FALSE)

# # look at histograms
# for (pcol in pattern_columns) {
#   p = ggplot(df_av1451, aes_string(pcol)) + geom_histogram(binwidth=0.1)
#   print(p)
# }

# make crossx response normal
#df_av1451[,eval(target)] = Gaussianize(df_av1451[,eval(target)], type='hh', method='MLE', return.u=TRUE)

# Formula setup
# all.addons = lapply(pattern_columns,lm.addvar)
# naive.addons = lapply(naive_columns,lm.addvar)
# braak.addons = lapply(braak_columns,lm.addvar)
all.addons = paste('+',paste(pattern_columns,collapse=' + '))
naive.addons = paste('+',paste(naive_columns,collapse=' + '))
braak.addons = paste('+',paste(braak_columns,collapse=' + '))
patterns_str = paste(all.addons,collapse=' ')
naive_str = paste(naive.addons,collapse=' ')
braak_str = paste(braak.addons,collapse=' ')

# diag_str = 'Diag.AV1451*APOE4_BIN +'
# diag_str = 'Diag.AV1451 +'
diag_str = ''



base_form = paste(target,"~",diag_str,"APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.")
braak_form = paste(target,"~",diag_str,"APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",braak_str)
pattern_form = paste(target,"~",diag_str,"APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",patterns_str)
naive_form = paste(target,"~",diag_str,"APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",naive_str)
full_form = paste(target,"~",diag_str,"APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",patterns_str,braak_str)
# onlypattern_form = str_replace(paste(target,"~",paste(all.addons,collapse=' ')),"\\+ ","")

# LARS lasso
braak_x = getxy(braak_form,df_av1451)
y = as.numeric(df_av1451[,target])
braak.lars.model = lars(braak_x,y,type='lasso')
braak.lars.test = covTest(braak.lars.model,braak_x,y)$results
braak.lars.test = data.frame(braak.lars.test[complete.cases(braak.lars.test),])
braak.lars.coef = coef(braak.lars.model, s=which.min(summary(braak.lars.model)$Cp), mode='step')
braak.lars.test$name = names(braak.lars.coef[braak.lars.test[,'Predictor_Number']])
braak.lars.test$coef = braak.lars.coef[braak.lars.test[,'Predictor_Number']]
braak.lars.sig = braak.lars.test[braak.lars.test$P.value <= 0.1,]
braak.lars.cp = min(summary(braak.lars.model)$Cp)
braak.lars.r2 = braak.lars.model$R2[which.min(summary(braak.lars.model)$Cp)]
braak.lars.n = attr(summary(braak.lars.model)$Cp,'n')
braak.lars.p = NROW(braak.lars.coef[braak.lars.coef != 0])
braak.lars.r2adj = r2adj(braak.lars.r2,braak.lars.n,braak.lars.p)
braak.lars.nonzero = braak.lars.test[braak.lars.test$coef != 0,'name']
paste(braak.lars.nonzero,collapse=' + ')
paste(target,'~',paste(braak.lars.sig[,'name'], collapse=' + '))

pattern_x = getxy(pattern_form,df_av1451)
y = as.numeric(df_av1451[,target])
pattern.lars.model = lars(pattern_x,y,type='lasso')
pattern.lars.test = covTest(pattern.lars.model,pattern_x,y)$results
pattern.lars.test = data.frame(pattern.lars.test[complete.cases(pattern.lars.test),])
pattern.lars.coef = coef(pattern.lars.model, s=which.min(summary(pattern.lars.model)$Cp), mode='step')
pattern.lars.test$name = names(pattern.lars.coef[pattern.lars.test[,'Predictor_Number']])
pattern.lars.test$coef = pattern.lars.coef[pattern.lars.test[,'Predictor_Number']]
pattern.lars.sig = pattern.lars.test[pattern.lars.test$P.value <= 0.1,]
pattern.lars.cp = min(summary(pattern.lars.model)$Cp)
pattern.lars.r2 = pattern.lars.model$R2[which.min(summary(pattern.lars.model)$Cp)]
pattern.lars.n = attr(summary(pattern.lars.model)$Cp,'n')
pattern.lars.p = NROW(pattern.lars.coef[pattern.lars.coef != 0])
pattern.lars.r2adj = r2adj(pattern.lars.r2,pattern.lars.n,pattern.lars.p)
pattern.lars.nonzero = pattern.lars.test[pattern.lars.test$coef != 0,'name']
paste(pattern.lars.nonzero,collapse=' + ')
paste(target,'~',paste(pattern.lars.sig[,'name'], collapse=' + '))

full_x = getxy(full_form,df_av1451)
y = as.numeric(df_av1451[,target])
full.lars.model = lars(full_x,y,type='lasso')
full.lars.test = covTest(full.lars.model,full_x,y)$results
full.lars.test = data.frame(full.lars.test[complete.cases(full.lars.test),])
full.lars.coef = coef(full.lars.model, s=which.min(summary(full.lars.model)$Cp), mode='step')
full.lars.test$name = names(full.lars.coef[full.lars.test[,'Predictor_Number']])
full.lars.test$coef = full.lars.coef[full.lars.test[,'Predictor_Number']]
full.lars.sig = full.lars.test[full.lars.test$P.value <= 0.1,]
full.lars.cp = min(summary(full.lars.model)$Cp)
full.lars.r2 = full.lars.model$R2[which.min(summary(full.lars.model)$Cp)]
full.lars.n = attr(summary(full.lars.model)$Cp,'n')
full.lars.p = NROW(full.lars.coef[full.lars.coef != 0])
full.lars.r2adj = r2adj(full.lars.r2,full.lars.n,full.lars.p)
full.lars.nonzero = full.lars.test[full.lars.test$coef != 0,'name']
paste(full.lars.nonzero,collapse=' + ')
paste(target,'~',paste(full.lars.sig[,'name'], collapse=' + '))

# r2 shrinkage
# test1.form = "ADAS_AV1451_1 ~ AV1451_Braak56_CerebGray_BL"
# test2.form = "ADAS_AV1451_1 ~ NSFA_0 + NSFA_3 + NSFA_8"
# test3.form = "ADAS_AV1451_1 ~ NSFA_0 + NSFA_3 + NSFA_8"

# test1.form = "AVLT_AV1451_1 ~ Age.AV1451"
# test2.form = "AVLT_AV1451_1 ~ NSFA_3"
# test3.form = "AVLT_AV1451_1 ~ NSFA_3"

# test1.form = "ADAS_retroslope_AV1451_BL ~ AV1451_Braak56_CerebGray_BL + APOE4_BIN"
# test2.form = "ADAS_retroslope_AV1451_BL ~ NSFA_0"
# test3.form = "ADAS_retroslope_AV1451_BL ~ NSFA_0"

# test2.form = "AV1451_BL_closest_AV45_wcereb_retroSlope ~ NSFA_7"

# test1.form = "AV1451_BL_closest_AV45_wcereb ~ AV1451_Braak12_CerebGray_BL + APOE4_BIN + Age.AV1451 + AV1451_Braak56_CerebGray_BL + Edu..Yrs."
# test2.form = "AV1451_BL_closest_AV45_wcereb ~ NSFA_0 + NSFA_3 + APOE4_BIN + NSFA_1 + Age.AV1451 + NSFA_7 + NSFA_6 + NSFA_10 + NSFA_4 + NSFA_2"
# test3.form = "AV1451_BL_closest_AV45_wcereb ~ AV1451_Braak12_CerebGray_BL + NSFA_0 + NSFA_3 + APOE4_BIN + Age.AV1451 + NSFA_10"

test1.form = "AV1451_BL_closest_AV45_wcereb ~ AV1451_Braak12_CerebGray_BL + Age.AV1451 + AV1451_Braak56_CerebGray_BL"
test2.form = "AV1451_BL_closest_AV45_wcereb ~ SCORE_NSFA_10 + APOE4_BIN + SCORE_NSFA_9"
test3.form = "AV1451_BL_closest_AV45_wcereb ~ SCORE_NSFA_10 + APOE4_BIN"

test1_results = c()
test2_results = c()
test3_results = c()
for (i in 1:50) {
  test1_r2 = r2.shrinkage(test1.form, target, df_av1451)[2]
  test2_r2 = r2.shrinkage(test2.form, target, df_av1451)[2]
  test3_r2 = r2.shrinkage(test3.form, target ,df_av1451)[2]
  test1_results = c(test1_results,test1_r2)
  test2_results = c(test2_results,test2_r2)
  test3_results = c(test3_results,test3_r2)
}
mean(test1_results)
mean(test2_results)
mean(test3_results)



# GLM net
link='binomial'

y = as.numeric(df_av1451[,target])
braak_x = getxy(braak_form,df_av1451)
pattern_x = getxy(pattern_form,df_av1451)
full_x = getxy(full_form,df_av1451)

braak.cvfit = cv.glmnet(braak_x,y,family=link)
pattern.cvfit = cv.glmnet(pattern_x,y,family=link)
full.cvfit = cv.glmnet(full_x,y,family=link)

braak.coef = coef(braak.cvfit,s='lambda.1se')
vars = rownames(braak.coef)[summary(braak.coef)$i[summary(braak.coef)$i != 1]]
paste(target,'~',paste(vars,collapse=' + '))

pattern.coef = coef(pattern.cvfit,s='lambda.1se')
vars = rownames(pattern.coef)[summary(pattern.coef)$i[summary(pattern.coef)$i != 1]]
paste(target,'~',paste(vars,collapse=' + '))

full.coef = coef(full.cvfit,s='lambda.1se')
vars = rownames(full.coef)[summary(full.coef)$i[summary(full.coef)$i != 1]]
paste(target,'~',paste(vars,collapse=' + '))

# form1 = "AV1451_BL_closest_AV45_wcereb_BIN1.11 ~ APOE4_BIN + Age.AV1451 + Gender + AV1451_Braak12_CerebGray_BL + AV1451_Braak56_CerebGray_BL"
# form2 = "AV1451_BL_closest_AV45_wcereb_BIN1.11 ~ APOE4_BIN + Age.AV1451 + Gender + NSFA_0 + NSFA_1 + NSFA_2 + NSFA_3 + NSFA_7 + NSFA_11"
# form3 = "AV1451_BL_closest_AV45_wcereb_BIN1.11 ~ APOE4_BIN + Age.AV1451 + Gender + NSFA_0 + NSFA_7 + AV1451_Braak12_CerebGray_BL + AV1451_Braak56_CerebGray_BL"

form1 = "AV1451_BL_closest_AV45_wcereb_BIN1.11 ~ APOE4_BIN + Age.AV1451 + Gender + AV1451_Braak12_CerebGray_BL + AV1451_Braak56_CerebGray_BL"
form2 = "AV1451_BL_closest_AV45_wcereb_BIN1.11 ~ APOE4_BIN + Age.AV1451 + Gender + SCORE_NSFA_0 + SCORE_NSFA_7 + SCORE_NSFA_10 + SCORE_NSFA_11"
form3 = "AV1451_BL_closest_AV45_wcereb_BIN1.11 ~ APOE4_BIN + Age.AV1451 + Gender + SCORE_NSFA_7 + SCORE_NSFA_10 + AV1451_Braak12_CerebGray_BL + AV1451_Braak56_CerebGray_BL"


summary(glm(form1,data=df_av1451,family=binomial()))
summary(glm(form2,data=df_av1451,family=binomial()))
summary(glm(form3,data=df_av1451,family=binomial()))

