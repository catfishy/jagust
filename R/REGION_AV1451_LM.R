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
to_factor = c('RID','APOE4_BIN','Gender',
              'Diag.AV1451','AV1451_BL_closest_AV45_wcereb_BIN1.11')
demog_columns = c('RID','APOE4_BIN','Diag.AV1451','Age.AV1451','Gender')
diag_columns = c('Diag.AV45','Diag.AV1451')

# target = "ADAS_AV1451_1"
# target = "AVLT_AV1451_1"
# target = "AV1451_BL_closest_AV45_wcereb_BIN1.11"
# target = "AV1451_BL_closest_AV45_wcereb_retroSlope"
# target = 'AV1451_BL_closest_AV45_wcereb'
# target = "ADAS_retroslope_AV1451_BL"
target = "AVLT_retroslope_AV1451_BL"


output_folder = 'R/output_av1451/'

valid_diags = c('N','SMC','EMCI','LMCI','AD')
#valid_diags = c('N','SMC','EMCI','LMCI')
#valid_diags = c('N','SMC')
#valid_diags = c('EMCI')
#valid_diags = c('LMCI')


# IMPORT 
master_csv = 'FDG_AV45_COGdata/FDG_AV45_COGdata_08_23_16.csv'
region_csv = 'datasets/pvc_adni_av1451/tauskullregions_uptake.csv'
df_master = read.csv(master_csv,skip=1)
columns = c('RID',
            'AV1451_CerebGray_BL',
            'AV1451_BL_closest_AV45_wcereb_BIN1.11',
            'AV1451_BL_closest_AV45_wcereb',
            'AV1451_BL_closest_AV45_wcereb_retroSlope',
            'Diag.AV1451',
            'Age.AV1451',
            'APOE4_BIN',
            'Gender',
            'ADAS_AV1451_1','ADAS_retroslope_AV1451_BL',
            'AVLT_AV1451_1','AVLT_retroslope_AV1451_BL')
df_master = df_master[,columns]
df_region = read.csv(region_csv)
df_region$RID = df_region$subject
fs_cols = c("CTX_BANKSSTS","CTX_PRECENTRAL","PALLIDUM",
            "CAUDATE","CEREBGM","CTX_LINGUAL","HEMIWM",
            "THALAMUS_PROPER","ORBITOFR","CTX_ISTHMUSCINGULATE",
            "AMYGDALA","CTX_PARAHIPPOCAMPAL","CTX_INFERIORTEMPORAL", 
            "CTX_LATERALOCCIPITAL","CTX_POSTCENTRAL",
            "CTX_CUNEUS","CTX_INSULA","CTX_ENTORHINAL","PARSFR",
            "CTX_SUPERIORTEMPORAL","CTX_INFERIORPARIETAL", 
            "PUTAMEN","MIDDLEFR","CTX_PERICALCARINE","CTX_SUPRAMARGINAL", 
            "CTX_ROSTRALANTERIORCINGULATE","CTX_PARACENTRAL",
            "CTX_TRANSVERSETEMPORAL","CTX_MIDDLETEMPORAL",
            "CTX_TEMPORALPOLE","CTX_POSTERIORCINGULATE",
            "CTX_SUPERIORPARIETAL","HIPPOCAMPUS",
            "ACCUMBENS_AREA","CTX_FUSIFORM", 
            "CTX_PRECUNEUS","CTX_SUPERIORFRONTAL",
            "CTX_CAUDALANTERIORCINGULATE")
df_region = df_region[,c('RID',fs_cols)]
df = merge(df_master,df_region,by='RID')

# FILTER
non.na = complete.cases(df[,c(demog_columns,target)])
df = df[non.na,]
for (i in names(df)){
  if (i %in% to_factor){
    df[,eval(i)] = as.factor(as.character(df[,eval(i)]))
  }
}

# Filter by diag
df = df[which(df$Diag.AV1451 %in% valid_diags),]

# Only amyloid positives
df = df[which(df$AV1451_BL_closest_AV45_wcereb_BIN1.11 == 1),]

# Formula setup
region.addons = paste('+',paste(fs_cols,collapse=' + '))
region_str = paste(region.addons,collapse=' ')

diag_str = ''

base_form = paste(target,"~",diag_str,"APOE4_BIN + Age.AV1451 + Gender")
region_form = paste(target,"~",diag_str,"APOE4_BIN + Age.AV1451 + Gender",region_str)

# LARS lasso
region_x = getxy(region_form,df)
y = as.numeric(df[,target])
region.lars.model = lars(region_x,y,type='lasso')
region.lars.test = covTest(region.lars.model,region_x,y)$results
region.lars.test = data.frame(region.lars.test[complete.cases(region.lars.test),])
region.lars.coef = coef(region.lars.model, s=which.min(summary(region.lars.model)$Cp), mode='step')
region.lars.test$name = names(region.lars.coef[region.lars.test[,'Predictor_Number']])
region.lars.test$coef = region.lars.coef[region.lars.test[,'Predictor_Number']]
region.lars.sig = region.lars.test[region.lars.test$P.value <= 0.1,]
region.lars.cp = min(summary(region.lars.model)$Cp)
region.lars.r2 = region.lars.model$R2[which.min(summary(region.lars.model)$Cp)]
region.lars.n = attr(summary(region.lars.model)$Cp,'n')
region.lars.p = NROW(region.lars.coef[region.lars.coef != 0])
region.lars.r2adj = r2adj(region.lars.r2,region.lars.n,region.lars.p)
region.lars.nonzero = region.lars.test[region.lars.test$coef != 0,'name']
paste(region.lars.nonzero,collapse=' + ')
paste(target,'~',paste(region.lars.sig[,'name'], collapse=' + '))


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
  test1_r2 = r2.shrinkage(test1.form, target, df)[2]
  test2_r2 = r2.shrinkage(test2.form, target, df)[2]
  test3_r2 = r2.shrinkage(test3.form, target ,df)[2]
  test1_results = c(test1_results,test1_r2)
  test2_results = c(test2_results,test2_r2)
  test3_results = c(test3_results,test3_r2)
}
mean(test1_results)
mean(test2_results)
mean(test3_results)



# GLM net
link='binomial'

y = as.numeric(df[,target])
region_x = getxy(region_form,df)

region.cvfit = cv.glmnet(region_x,y,family=link)

region.coef = coef(region.cvfit,s='lambda.1se')
vars = rownames(region.coef)[summary(region.coef)$i[summary(region.coef)$i != 1]]
paste(target,'~',paste(vars,collapse=' + '))


form1 = "AV1451_BL_closest_AV45_wcereb_BIN1.11 ~ APOE4_BIN + Age.AV1451 + Gender + AV1451_Braak12_CerebGray_BL + AV1451_Braak56_CerebGray_BL"
form2 = "AV1451_BL_closest_AV45_wcereb_BIN1.11 ~ APOE4_BIN + Age.AV1451 + Gender + NSFA_0 + NSFA_1 + NSFA_3 + NSFA_7 + NSFA_11"
form3 = "AV1451_BL_closest_AV45_wcereb_BIN1.11 ~ APOE4_BIN + Age.AV1451 + Gender + NSFA_0 + NSFA_7 + NSFA_11 + AV1451_Braak12_CerebGray_BL + AV1451_Braak56_CerebGray_BL"

# form1 = "AV1451_BL_closest_AV45_wcereb_BIN1.11 ~ APOE4_BIN + Age.AV1451 + Gender + AV1451_Braak12_CerebGray_BL + AV1451_Braak56_CerebGray_BL"
# form2 = "AV1451_BL_closest_AV45_wcereb_BIN1.11 ~ APOE4_BIN + Age.AV1451 + Gender + SCORE_NSFA_10"
# form3 = "AV1451_BL_closest_AV45_wcereb_BIN1.11 ~ APOE4_BIN + Age.AV1451 + Gender + SCORE_NSFA_3 + SCORE_NSFA_7 + SCORE_NSFA_10 + AV1451_Braak12_CerebGray_BL + AV1451_Braak56_CerebGray_BL"


summary(glm(form1,data=df,family=binomial()))
summary(glm(form2,data=df,family=binomial()))
summary(glm(form3,data=df,family=binomial()))

