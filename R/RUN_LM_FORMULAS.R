library(lme4)
library(ggplot2)
library(stats)
library(gdata)
library(reshape2)
library(caret)
library(cvTools)
library(lmtest)
library(stringr)
library(Jmisc)

source('R/LM_FUNCS.R')

# # sparsest local minima
# av45_pos_lm_list = c("regions.ADAS_AV1451_1"="ADAS_AV1451_1 ~ CTX_ISTHMUSCINGULATE + AMYGDALA + CTX_MIDDLETEMPORAL + BRAIN_STEM", "regions.ADAS_retroslope_AV1451_BL"="ADAS_retroslope_AV1451_BL ~ APOE4_BIN + Gender + Edu..Yrs. + CTX_ISTHMUSCINGULATE + AMYGDALA + CTX_LATERALOCCIPITAL + CTX_ENTORHINAL + PUTAMEN + CTX_PARACENTRAL + CTX_MIDDLETEMPORAL + CEREBWM + ACCUMBENS_AREA", "regions.AVLT_AV1451_1"="AVLT_AV1451_1 ~ AMYGDALA", "regions.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ APOE4_BIN + Edu..Yrs. + CTX_ENTORHINAL", "regions.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ Gender + CTX_LINGUAL + AMYGDALA + CTX_ENTORHINAL", "regions.AV1451_BL_closest_AV45_wcereb_retroSlope"="None", "patterns.ADAS_AV1451_1"="ADAS_AV1451_1 ~ NSFA_0 + NSFA_3 + NSFA_4 + NSFA_8", "patterns.ADAS_retroslope_AV1451_BL"="ADAS_retroslope_AV1451_BL ~ NSFA_0", "patterns.AVLT_AV1451_1"="AVLT_AV1451_1 ~ NSFA_0 + NSFA_1 + NSFA_3 + NSFA_4 + NSFA_7 + NSFA_8", "patterns.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ APOE4_BIN", "patterns.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ APOE4_BIN + Gender + NSFA_0 + NSFA_1 + NSFA_3 + NSFA_11", "patterns.AV1451_BL_closest_AV45_wcereb_retroSlope"="None", "full.ADAS_AV1451_1"="ADAS_AV1451_1 ~ AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak4_CerebGray_BL", "full.ADAS_retroslope_AV1451_BL"="ADAS_retroslope_AV1451_BL ~ AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak4_CerebGray_BL", "full.AVLT_AV1451_1"="AVLT_AV1451_1 ~ NSFA_1 + NSFA_3 + NSFA_4 + NSFA_7 + NSFA_8 + AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak2_CerebGray_BL", "full.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ APOE4_BIN + AV1451_PVC_Braak1_CerebGray_BL", "full.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ APOE4_BIN + Gender + NSFA_0 + NSFA_3 + AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak2_CerebGray_BL", "full.AV1451_BL_closest_AV45_wcereb_retroSlope"="None", "braak.ADAS_AV1451_1"="ADAS_AV1451_1 ~ AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak4_CerebGray_BL", "braak.ADAS_retroslope_AV1451_BL"="ADAS_retroslope_AV1451_BL ~ APOE4_BIN + Edu..Yrs. + AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak4_CerebGray_BL", "braak.AVLT_AV1451_1"="AVLT_AV1451_1 ~ AV1451_PVC_Braak1_CerebGray_BL", "braak.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ APOE4_BIN + Edu..Yrs. + AV1451_PVC_Braak1_CerebGray_BL", "braak.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ APOE4_BIN + Gender + AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak2_CerebGray_BL + AV1451_PVC_Braak3_CerebGray_BL", "braak.AV1451_BL_closest_AV45_wcereb_retroSlope"="AV1451_BL_closest_AV45_wcereb_retroSlope ~ Edu..Yrs.", "onlyregions.ADAS_AV1451_1"="ADAS_AV1451_1 ~ CTX_ISTHMUSCINGULATE + AMYGDALA + CTX_LATERALOCCIPITAL + PUTAMEN + CTX_PERICALCARINE + CTX_MIDDLETEMPORAL + CEREBWM + BRAIN_STEM", "onlyregions.ADAS_retroslope_AV1451_BL"="ADAS_retroslope_AV1451_BL ~ AMYGDALA + CTX_ENTORHINAL + CTX_MIDDLETEMPORAL + CEREBWM", "onlyregions.AVLT_AV1451_1"="AVLT_AV1451_1 ~ AMYGDALA", "onlyregions.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ AMYGDALA", "onlyregions.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ CTX_LINGUAL + AMYGDALA + CTX_ENTORHINAL + HIPPOCAMPUS", "onlyregions.AV1451_BL_closest_AV45_wcereb_retroSlope"="None")
# all_subj_lm_list = c("regions.ADAS_AV1451_1"="ADAS_AV1451_1 ~ CTX_ISTHMUSCINGULATE + AMYGDALA + CTX_ENTORHINAL + CTX_MIDDLETEMPORAL + CTX_PRECUNEUS", "regions.ADAS_retroslope_AV1451_BL"="ADAS_retroslope_AV1451_BL ~ APOE4_BIN + Edu..Yrs. + AMYGDALA + CTX_ENTORHINAL + CTX_MIDDLETEMPORAL + CEREBWM", "regions.AVLT_AV1451_1"="AVLT_AV1451_1 ~ Age.AV1451 + Gender + Edu..Yrs. + AMYGDALA + CTX_ENTORHINAL + CEREBWM + BRAIN_STEM", "regions.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ AMYGDALA", "regions.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ APOE4_BIN + Age.AV1451 + AMYGDALA + CTX_ENTORHINAL + CEREBWM", "regions.AV1451_BL_closest_AV45_wcereb_retroSlope"="AV1451_BL_closest_AV45_wcereb_retroSlope ~ APOE4_BIN + Age.AV1451 + Edu..Yrs. + CTX_PARAHIPPOCAMPAL + CTX_CUNEUS + CTX_TRANSVERSETEMPORAL + HIPPOCAMPUS + CTX_CAUDALANTERIORCINGULATE", "patterns.ADAS_AV1451_1"="ADAS_AV1451_1 ~ NSFA_0 + NSFA_3 + NSFA_4 + NSFA_8", "patterns.ADAS_retroslope_AV1451_BL"="ADAS_retroslope_AV1451_BL ~ NSFA_0", "patterns.AVLT_AV1451_1"="AVLT_AV1451_1 ~ NSFA_3 + NSFA_5", "patterns.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ APOE4_BIN + Age.AV1451 + NSFA_0 + NSFA_1", "patterns.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ APOE4_BIN + Age.AV1451 + NSFA_0 + NSFA_1 + NSFA_3 + NSFA_4 + NSFA_6 + NSFA_7 + NSFA_10", "patterns.AV1451_BL_closest_AV45_wcereb_retroSlope"="AV1451_BL_closest_AV45_wcereb_retroSlope ~ NSFA_7", "full.ADAS_AV1451_1"="ADAS_AV1451_1 ~ NSFA_0 + NSFA_3 + NSFA_4 + NSFA_8 + AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak4_CerebGray_BL", "full.ADAS_retroslope_AV1451_BL"="ADAS_retroslope_AV1451_BL ~ APOE4_BIN + NSFA_0 + NSFA_4 + AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak4_CerebGray_BL", "full.AVLT_AV1451_1"="AVLT_AV1451_1 ~ Age.AV1451 + Gender + Edu..Yrs. + NSFA_3 + NSFA_5 + AV1451_PVC_Braak1_CerebGray_BL", "full.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ AV1451_PVC_Braak1_CerebGray_BL", "full.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ APOE4_BIN + Age.AV1451 + NSFA_0 + NSFA_3 + NSFA_4 + AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak2_CerebGray_BL + AV1451_PVC_Braak3_CerebGray_BL", "full.AV1451_BL_closest_AV45_wcereb_retroSlope"="AV1451_BL_closest_AV45_wcereb_retroSlope ~ NSFA_7", "braak.ADAS_AV1451_1"="ADAS_AV1451_1 ~ AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak4_CerebGray_BL", "braak.ADAS_retroslope_AV1451_BL"="ADAS_retroslope_AV1451_BL ~ APOE4_BIN + Edu..Yrs. + AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak4_CerebGray_BL", "braak.AVLT_AV1451_1"="AVLT_AV1451_1 ~ APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs. + AV1451_PVC_Braak1_CerebGray_BL", "braak.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ AV1451_PVC_Braak1_CerebGray_BL", "braak.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ APOE4_BIN + Age.AV1451 + AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak2_CerebGray_BL + AV1451_PVC_Braak3_CerebGray_BL + AV1451_PVC_Braak6_CerebGray_BL", "braak.AV1451_BL_closest_AV45_wcereb_retroSlope"="AV1451_BL_closest_AV45_wcereb_retroSlope ~ APOE4_BIN + Age.AV1451 + Edu..Yrs. + AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak2_CerebGray_BL + AV1451_PVC_Braak6_CerebGray_BL", "onlyregions.ADAS_AV1451_1"="ADAS_AV1451_1 ~ CTX_ISTHMUSCINGULATE + AMYGDALA + CTX_ENTORHINAL + CTX_MIDDLETEMPORAL + CTX_PRECUNEUS", "onlyregions.ADAS_retroslope_AV1451_BL"="ADAS_retroslope_AV1451_BL ~ AMYGDALA + CTX_ENTORHINAL + CTX_MIDDLETEMPORAL", "onlyregions.AVLT_AV1451_1"="AVLT_AV1451_1 ~ AMYGDALA + CTX_ENTORHINAL + CEREBWM + BRAIN_STEM", "onlyregions.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ AMYGDALA", "onlyregions.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ PALLIDUM + AMYGDALA + CTX_ENTORHINAL + CEREBWM", "onlyregions.AV1451_BL_closest_AV45_wcereb_retroSlope"="AV1451_BL_closest_AV45_wcereb_retroSlope ~ CTX_PARAHIPPOCAMPAL + CTX_CUNEUS + CTX_ENTORHINAL + CTX_TRANSVERSETEMPORAL + HIPPOCAMPUS + CTX_CAUDALANTERIORCINGULATE")
# av45_neg_lm_list = c("regions.ADAS_AV1451_1"="ADAS_AV1451_1 ~ Age.AV1451 + Edu..Yrs. + PALLIDUM + CTX_LATERALOCCIPITAL + CTX_ENTORHINAL + PARSFR + CTX_SUPERIORTEMPORAL + PUTAMEN + CTX_ROSTRALANTERIORCINGULATE + CTX_PARACENTRAL + CTX_POSTERIORCINGULATE + CTX_SUPERIORPARIETAL + CEREBWM", "regions.ADAS_retroslope_AV1451_BL"="ADAS_retroslope_AV1451_BL ~ APOE4_BIN + Age.AV1451 + Edu..Yrs. + PALLIDUM + CAUDATE + AMYGDALA + CTX_LATERALOCCIPITAL + CTX_CUNEUS + CTX_SUPERIORTEMPORAL + MIDDLEFR + CTX_TRANSVERSETEMPORAL + CTX_TEMPORALPOLE + CTX_POSTERIORCINGULATE + CTX_SUPERIORPARIETAL + HIPPOCAMPUS + ACCUMBENS_AREA + CTX_PRECUNEUS", "regions.AVLT_AV1451_1"="AVLT_AV1451_1 ~ CTX_ROSTRALANTERIORCINGULATE + CEREBWM", "regions.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ Age.AV1451 + CTX_PARAHIPPOCAMPAL + PUTAMEN + ACCUMBENS_AREA + CTX_CAUDALANTERIORCINGULATE", "regions.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ APOE4_BIN + PUTAMEN", "regions.AV1451_BL_closest_AV45_wcereb_retroSlope"="AV1451_BL_closest_AV45_wcereb_retroSlope ~ APOE4_BIN + ORBITOFR + CTX_SUPRAMARGINAL", "patterns.ADAS_AV1451_1"="ADAS_AV1451_1 ~ NSFA_3 + NSFA_5", "patterns.ADAS_retroslope_AV1451_BL"="", "patterns.AVLT_AV1451_1"="AVLT_AV1451_1 ~ Age.AV1451 + Gender + NSFA_5", "patterns.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ Age.AV1451 + Gender + NSFA_4 + NSFA_7", "patterns.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ APOE4_BIN + Gender + NSFA_0 + NSFA_10", "patterns.AV1451_BL_closest_AV45_wcereb_retroSlope"="AV1451_BL_closest_AV45_wcereb_retroSlope ~ APOE4_BIN", "full.ADAS_AV1451_1"="ADAS_AV1451_1 ~ NSFA_3 + NSFA_5", "full.ADAS_retroslope_AV1451_BL"="", "full.AVLT_AV1451_1"="AVLT_AV1451_1 ~ Age.AV1451 + Gender + NSFA_3 + NSFA_5 + AV1451_PVC_Braak1_CerebGray_BL", "full.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ Age.AV1451 + NSFA_4 + NSFA_7 + AV1451_PVC_Braak2_CerebGray_BL", "full.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ APOE4_BIN + Gender + NSFA_0 + NSFA_10", "full.AV1451_BL_closest_AV45_wcereb_retroSlope"="AV1451_BL_closest_AV45_wcereb_retroSlope ~ APOE4_BIN + NSFA_5 + AV1451_PVC_Braak4_CerebGray_BL + AV1451_PVC_Braak6_CerebGray_BL", "braak.ADAS_AV1451_1"="ADAS_AV1451_1 ~ Age.AV1451", "braak.ADAS_retroslope_AV1451_BL"="", "braak.AVLT_AV1451_1"="AVLT_AV1451_1 ~ APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs. + AV1451_PVC_Braak1_CerebGray_BL + AV1451_PVC_Braak4_CerebGray_BL + AV1451_PVC_Braak6_CerebGray_BL + AV1451_PVC_BraakAll_CerebGray_BL", "braak.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ Age.AV1451 + AV1451_PVC_Braak2_CerebGray_BL", "braak.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ APOE4_BIN + Gender", "braak.AV1451_BL_closest_AV45_wcereb_retroSlope"="AV1451_BL_closest_AV45_wcereb_retroSlope ~ APOE4_BIN + AV1451_PVC_Braak4_CerebGray_BL + AV1451_PVC_Braak6_CerebGray_BL", "onlyregions.ADAS_AV1451_1"="ADAS_AV1451_1 ~ PALLIDUM + CTX_INFERIORTEMPORAL + CTX_LATERALOCCIPITAL + CTX_ENTORHINAL + PARSFR + CTX_SUPERIORTEMPORAL + PUTAMEN + CTX_PARACENTRAL + CTX_POSTERIORCINGULATE + CTX_SUPERIORPARIETAL + CEREBWM", "onlyregions.ADAS_retroslope_AV1451_BL"="", "onlyregions.AVLT_AV1451_1"="AVLT_AV1451_1 ~ CTX_ROSTRALANTERIORCINGULATE + CEREBWM", "onlyregions.AVLT_retroslope_AV1451_BL"="AVLT_retroslope_AV1451_BL ~ CTX_PARAHIPPOCAMPAL + PUTAMEN + ACCUMBENS_AREA + CTX_CAUDALANTERIORCINGULATE", "onlyregions.AV1451_BL_closest_AV45_wcereb"="AV1451_BL_closest_AV45_wcereb ~ PUTAMEN", "onlyregions.AV1451_BL_closest_AV45_wcereb_retroSlope"="AV1451_BL_closest_AV45_wcereb_retroSlope ~ ORBITOFR + CTX_SUPRAMARGINAL")
#   
valid_diags = c('N','SMC','EMCI','LMCI','AD')
demog_columns = c('RID','APOE4_BIN','Diag.AV1451','Age.AV1451','Gender','Edu..Yrs.')
to_factor = c('RID','APOE4_BIN','APOE2_BIN','Gender',
              'Diag.AV45','Diag.AV1451',
              'AV45_NONTP_wcereb_BIN1.11','AV1451_BL_closest_AV45_wcereb_BIN1.11')
to_standardize = c('Age.AV45','Edu..Yrs.','Age.AV1451')
braak_columns = c('AV1451_PVC_Braak1_CerebGray_BL',
                  'AV1451_PVC_Braak2_CerebGray_BL',
                  'AV1451_PVC_Braak3_CerebGray_BL',
                  'AV1451_PVC_Braak4_CerebGray_BL',
                  'AV1451_PVC_Braak5_CerebGray_BL',
                  'AV1451_PVC_Braak6_CerebGray_BL')
all_targets = c("ADAS_AV1451_1",
                "ADAS_retroslope_AV1451_BL",
                "AVLT_AV1451_1",
                "AVLT_retroslope_AV1451_BL",
                'AV1451_BL_closest_AV45_wcereb',
                "AV1451_BL_closest_AV45_wcereb_retroSlope")

df_av1451 = read.csv('nsfa/av1451skull_pattern_dataset.csv')
df_regions = read.csv('datasets/pvc_adni_av1451/tauskullregions_uptake_bilateral.csv')
df_regions = rename(df_regions,c("subject"="RID"))
df_regions$RID = as.factor(df_regions$RID)
all_regions = colnames(df_regions)
all_regions = all_regions[all_regions != 'RID']
df_av1451 = merge(df_av1451,df_regions,by='RID')

df_av1451$Gender = df_av1451$Gender - 1
pattern_columns = Filter(isPatternColumn,names(df_av1451))
non.na = complete.cases(df_av1451[,c(demog_columns,braak_columns)])
df_av1451 = df_av1451[non.na,]
for (i in names(df_av1451)){
  if (i %in% to_factor){
    df_av1451[,eval(i)] = as.factor(as.character(df_av1451[,eval(i)]))
  }
}

df_av1451 = df_av1451[which(df_av1451$Diag.AV1451 %in% valid_diags),]

# Filter by AV45 status
df_av1451_pos = df_av1451[which(df_av1451$AV1451_BL_closest_AV45_wcereb_BIN1.11 == 1),]
df_av1451_neg = df_av1451[which(df_av1451$AV1451_BL_closest_AV45_wcereb_BIN1.11 == 0),]

# Standardize
cross_to_standardize = c(to_standardize,pattern_columns,braak_columns,all_regions,all_targets)
cross_normalization = preProcess(df_av1451[,cross_to_standardize])
df_av1451[,cross_to_standardize] = predict(cross_normalization, df_av1451[,cross_to_standardize])
cross_normalization = preProcess(df_av1451_pos[,cross_to_standardize])
df_av1451_pos[,cross_to_standardize] = predict(cross_normalization, df_av1451_pos[,cross_to_standardize])

# Create regression formulas
# pen_type = 'lars'
pen_type = 'glmnet'

# run_type = 'all_subj'
# df_av1451 = df_av1451

run_type = 'av45_pos'
df_av1451 = df_av1451_pos

# run_type = 'av45_neg'
# df_av1451 = df_av1451_neg

full.coef.csv = paste(pen_type,'_',run_type,'_coef_full.csv',sep='')
braak.coef.csv = paste(pen_type,'_',run_type,'_coef_braak.csv',sep='')
pattern.coef.csv = paste(pen_type,'_',run_type,'_coef_patterns.csv',sep='')
region.coef.csv = paste(pen_type,'_',run_type,'_coef_regions.csv',sep='')
onlyregion.coef.csv = paste(pen_type,'_',run_type,'_coef_onlyregions.csv',sep='')

full.lmforms = create.lm.from.coefcsv(full.coef.csv)
braak.lmforms = create.lm.from.coefcsv(braak.coef.csv)
pattern.lmforms = create.lm.from.coefcsv(pattern.coef.csv)
region.lmforms = create.lm.from.coefcsv(region.coef.csv)
onlyregion.lmforms = create.lm.from.coefcsv(onlyregion.coef.csv)
all.lmforms = c(full=full.lmforms,
                braak=braak.lmforms,
                pattern=pattern.lmforms,
                region=region.lmforms,
                onlyregion=onlyregion.lmforms)

run.lm = function(target, lm_form) {
  result = data.frame()
  non.na = complete.cases(df_av1451[,c(target)])
  df_cur = df_av1451[non.na,]
  fm = lm(lm_form,df_cur)
  summary(fm)
}

# Run linear models
r2_output_file = paste('/usr/local/jagust/',pen_type,'_',run_type,'_r2.csv',sep='')
r2 = data.frame()
for (cur_key in names(all.lmforms)) {
  parts = strsplit(cur_key, "[.]")[[1]]
  target = paste(parts[2:length(parts)],collapse='.')
  if (target == "AV1451_BL_closest_AV45_wcereb_BIN1.11") {
    print("Skipping AV45+/-")
    next
  }
  key = parts[1]
  lm_form = all.lmforms[[cur_key]]
  if (lm_form == "") {
    next
  }
  fm.summary = run.lm(target, lm_form)
  r2.shrink = r2.shrinkage.mean(lm_form, 
                                target, 
                                df_av1451[complete.cases(df_av1451[,c(target)]),])
  r2[cur_key,'shrink.r2'] = r2.shrink
  r2[cur_key,'adj.r2'] = fm.summary$adj.r.squared
  r2[cur_key,'target'] = target
  r2[cur_key,'lm_type'] = key
  outputfile = paste('/usr/local/jagust/',pen_type,'_',run_type,'_',cur_key,'_lm_summary.txt',sep='')
  capture.output(fm.summary, file=outputfile)
}
write.csv(r2,r2_output_file)


# make heatmap
base_size = 16
p = ggplot(all_subj_r2, aes(lmtype, target)) + 
  geom_tile(aes(fill = adj.r2), 
            colour='grey50') + 
  scale_fill_gradient2(low="red", high="blue") + 
  labs(x="",y="") + 
  theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio=1) + 
  theme(panel.border=element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=base_size*0.8),
        axis.text.x = element_text(size=base_size*0.8, 
                                   angle = 300, 
                                   hjust = 0))
print(p)

base_size = 16
p = ggplot(av45_pos_r2, aes(lmtype, target)) + 
  geom_tile(aes(fill = adj.r2), 
            colour='grey50') + 
  scale_fill_gradient2(low="red", high="blue") + 
  labs(x="",y="") + 
  theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio=1) + 
  theme(panel.border=element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=base_size*0.8),
        axis.text.x = element_text(size=base_size*0.8, 
                                   angle = 300, 
                                   hjust = 0))
print(p)

base_size = 16
p = ggplot(av45_neg_r2, aes(lmtype, target)) + 
  geom_tile(aes(fill = adj.r2), 
            colour='grey50') + 
  scale_fill_gradient2(low="red", high="blue") + 
  labs(x="",y="") + 
  theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio=1) + 
  theme(panel.border=element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=base_size*0.8),
        axis.text.x = element_text(size=base_size*0.8, 
                                   angle = 300, 
                                   hjust = 0))
print(p)