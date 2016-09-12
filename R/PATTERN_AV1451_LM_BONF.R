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
library(plyr)

source('R/LM_FUNCS.R')

# CONSTANTS
pattern_prefix = 'NSFA_'
to_factor = c('RID','APOE4_BIN','APOE2_BIN','Gender',
              'Diag.AV45','Diag.AV1451',
              'AV1451_BL_closest_AV45_wcereb_BIN1.11')
to_standardize = c('Age.AV45','Edu..Yrs.','Age.AV1451')
demog_columns = c('RID','APOE4_BIN','Diag.AV1451','Age.AV1451','Gender','Edu..Yrs.')
diag_columns = c('Diag.AV45','Diag.AV1451')
braak_columns = c('AV1451_PVC_Braak1_CerebGray_BL',
                  'AV1451_PVC_Braak2_CerebGray_BL',
                  'AV1451_PVC_Braak3_CerebGray_BL',
                  'AV1451_PVC_Braak4_CerebGray_BL',
                  'AV1451_PVC_Braak5_CerebGray_BL',
                  'AV1451_PVC_Braak6_CerebGray_BL')
braakmax_columns = c('BRAAK1_MAX','BRAAK2_MAX','BRAAK3_MAX',
                     'BRAAK4_MAX','BRAAK5_MAX','BRAAK6_MAX')

output_folder = 'R/output_av1451/'

valid_diags = c('N','SMC','EMCI','LMCI','AD')

# IMPORT
df_av1451 = read.csv('nsfa/av1451skull_pattern_dataset.csv')
df_regions = read.csv('datasets/pvc_adni_av1451/tauskullregions_uptake_bilateral.csv')
df_regions = rename(df_regions,c("subject"="RID"))
df_regions$RID = as.factor(df_regions$RID)
all_regions = colnames(df_regions)
all_regions = all_regions[all_regions != 'RID']
df_av1451 = merge(df_av1451,df_regions,by='RID')

df_max = read.csv('output/08-23-2016/UCBERKELEYAV1451_MAX_08-23-2016_regular_tp.csv')
max_cols = c('RID','BRAAK1','BRAAK2','BRAAK3','BRAAK4','BRAAK5','BRAAK6')
df_max = df_max[,max_cols]
colnames(df_max) = c('RID','BRAAK1_MAX','BRAAK2_MAX','BRAAK3_MAX','BRAAK4_MAX','BRAAK5_MAX','BRAAK6_MAX')
df_av1451 = merge(df_av1451,df_max,by='RID')

# TODO: IMPORT MASK MEAN/MAX
df_mask = read.csv('maskdata_inorm_cerebgm_snorm.csv')

df_av1451$Gender = df_av1451$Gender - 1
pattern_columns = Filter(isPatternColumn,names(df_av1451))
non.na = complete.cases(df_av1451[,c(demog_columns,braak_columns,braakmax_columns)])
df_av1451 = df_av1451[non.na,]
for (i in names(df_av1451)){
  if (i %in% to_factor){
    df_av1451[,eval(i)] = as.factor(as.character(df_av1451[,eval(i)]))
  }
}

# Filter by diag
df_av1451 = df_av1451[which(df_av1451$Diag.AV1451 %in% valid_diags),]

# Filter by AV45 status
df_av1451_pos = df_av1451[which(df_av1451$AV1451_BL_closest_AV45_wcereb_BIN1.11 == 1),]
df_av1451_neg = df_av1451[which(df_av1451$AV1451_BL_closest_AV45_wcereb_BIN1.11 == 0),]

# Choose run

use_min = FALSE

df_av1451 = df_av1451_pos
runname = 'av45_pos'

# df_av1451 = df_av1451_neg
# runname = 'av45_neg'

# df_av1451 = df_av1451
# runname = 'all_subj'

# Standardize variables
cross_to_standardize = c(to_standardize,pattern_columns,
                         braak_columns,braakmax_columns,all_regions)
cross_normalization = preProcess(df_av1451[,cross_to_standardize])
df_av1451[,cross_to_standardize] = predict(cross_normalization, df_av1451[,cross_to_standardize])

# Formula setup
all.addons = paste('+',paste(pattern_columns,collapse=' + '))
braak.addons = paste('+',paste(braak_columns,collapse=' + '))
braakmax.addons = paste('+',paste(braakmax_columns,collapse=' + '))
regions.addons = paste('+',paste(all_regions,collapse=' + '))
patterns_str = paste(all.addons,collapse=' ')
braak_str = paste(braak.addons,collapse=' ')
braakmax_str = paste(braakmax.addons,collapse=' ')
regions_str = paste(regions.addons,collapse=' ')

# all_targets = c("UW_MEM_AV1451_1",
#                 "UW_EF_AV1451_1",
#                 "ADAS_AV1451_1",
#                 "AVLT_AV1451_1",
#                 'AV1451_BL_closest_AV45_wcereb',
#                 "ADAS_retroslope_AV1451_BL",
#                 "AVLT_retroslope_AV1451_BL",
#                 "UW_MEM_retroslope_AV1451_BL",
#                 "UW_EF_retroslope_AV1451_BL",
#                 "AV1451_BL_closest_AV45_wcereb_retroSlope")
all_targets = c("AVLT_AV1451_1",
                "AVLT_retroslope_AV1451_BL",
                "ADAS_AV1451_1",
                "ADAS_retroslope_AV1451_BL",
                'AV1451_BL_closest_AV45_wcereb',
                "AV1451_BL_closest_AV45_wcereb_retroSlope")
# Standardize targets
target_norm = preProcess(df_av1451[,all_targets])
df_av1451[,all_targets] = predict(target_norm, df_av1451[,all_targets])

# add logistic regression targets
all_targets = c(all_targets,'AV1451_BL_closest_AV45_wcereb_BIN1.11')

# penalized regression against cognitive targets
lars_coef_full = data.frame()
lars_coef_fullmax = data.frame()
lars_coef_braak = data.frame()
lars_coef_braakmax = data.frame()
lars_coef_pattern = data.frame()
lars_coef_region = data.frame()
lars_coef_onlyregion = data.frame()
glmnet_coef_full = data.frame()
glmnet_coef_fullmax = data.frame()
glmnet_coef_braak = data.frame()
glmnet_coef_braakmax = data.frame()
glmnet_coef_pattern = data.frame()
glmnet_coef_region = data.frame()
glmnet_coef_onlyregion = data.frame()
for (target in all_targets){
  # remove na target rows
  non.na = complete.cases(df_av1451[,c(target)])
  df_cur = df_av1451[non.na,]
  if (target == 'AV1451_BL_closest_AV45_wcereb_BIN1.11'){
    if (runname != 'all_subj') {
      print("Skipping AV45 +/- target")
      next
    }
    family = 'binomial'
    y = as.factor(df_cur[,target])
    run_lars = FALSE
  } else {
    family = 'gaussian'
    y = as.numeric(df_cur[,target])
    run_lars = TRUE
  }
  
  full_form = paste(target,"~","APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",patterns_str,braak_str)
  fullmax_form = paste(target,"~","APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",patterns_str,braakmax_str)
  braak_form = paste(target,"~","APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",braak_str)
  braakmax_form = paste(target,"~","APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",braakmax_str)
  pattern_form = paste(target,"~","APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",patterns_str)
  region_form = paste(target,'~',"APOE4_BIN + Age.AV1451 + Gender + Edu..Yrs.",regions_str)
  onlyregion_form = paste(target,'~',paste(all_regions,collapse=' + '))
  full_x = getxy(full_form,df_cur,FALSE)
  fullmax_x = getxy(fullmax_form,df_cur,FALSE)
  braak_x = getxy(braak_form,df_cur,FALSE)
  braakmax_x = getxy(braakmax_form,df_cur,FALSE)
  pattern_x = getxy(pattern_form,df_cur,FALSE)
  region_x = getxy(region_form,df_cur,FALSE)
  onlyregion_x = getxy(onlyregion_form,df_cur,FALSE)
  
  print(target)
  glmnet_coef_full = rbind(glmnet_coef_full,get.glmnet.coeff(full_x,y,target,use_min=use_min,family=family))
  glmnet_coef_fullmax = rbind(glmnet_coef_fullmax,get.glmnet.coeff(fullmax_x,y,target,use_min=use_min,family=family))
  glmnet_coef_braak = rbind(glmnet_coef_braak,get.glmnet.coeff(braak_x,y,target,use_min=use_min,family=family))
  glmnet_coef_braakmax = rbind(glmnet_coef_braakmax,get.glmnet.coeff(braakmax_x,y,target,use_min=use_min,family=family))
  glmnet_coef_pattern = rbind(glmnet_coef_pattern,get.glmnet.coeff(pattern_x,y,target,use_min=use_min,family=family))
  glmnet_coef_region = rbind(glmnet_coef_region,get.glmnet.coeff(region_x,y,target,use_min=use_min,family=family))
  glmnet_coef_onlyregion = rbind(glmnet_coef_onlyregion,get.glmnet.coeff(onlyregion_x,y,target,use_min=use_min,family=family))
  if (run_lars) {
    lars_coef_full = rbind(lars_coef_full,get.lars.coeff(full_x,y,target))
    lars_coef_fullmax = rbind(lars_coef_fullmax,get.lars.coeff(fullmax_x,y,target))
    lars_coef_braak = rbind(lars_coef_braak,get.lars.coeff(braak_x,y,target))
    lars_coef_braakmax = rbind(lars_coef_braakmax,get.lars.coeff(braakmax_x,y,target))
    lars_coef_pattern = rbind(lars_coef_pattern,get.lars.coeff(pattern_x,y,target))
    lars_coef_region = rbind(lars_coef_region,get.lars.coeff(region_x,y,target))
    lars_coef_onlyregion = rbind(lars_coef_onlyregion,get.lars.coeff(onlyregion_x,y,target))
  }
}

write.csv(file=paste('lars_',runname,'_coef_full.csv',sep=''),x=lars_coef_full,row.names=TRUE)
write.csv(file=paste('glmnet_',runname,'_coef_full.csv',sep=''),x=glmnet_coef_full,row.names=TRUE)
write.csv(file=paste('lars_',runname,'_coef_fullmax.csv',sep=''),x=lars_coef_fullmax,row.names=TRUE)
write.csv(file=paste('glmnet_',runname,'_coef_fullmax.csv',sep=''),x=glmnet_coef_fullmax,row.names=TRUE)
write.csv(file=paste('lars_',runname,'_coef_braak.csv',sep=''),x=lars_coef_braak,row.names=TRUE)
write.csv(file=paste('glmnet_',runname,'_coef_braak.csv',sep=''),x=glmnet_coef_braak,row.names=TRUE)
write.csv(file=paste('lars_',runname,'_coef_braakmax.csv',sep=''),x=lars_coef_braakmax,row.names=TRUE)
write.csv(file=paste('glmnet_',runname,'_coef_braakmax.csv',sep=''),x=glmnet_coef_braakmax,row.names=TRUE)
write.csv(file=paste('lars_',runname,'_coef_patterns.csv',sep=''),x=lars_coef_pattern,row.names=TRUE)
write.csv(file=paste('glmnet_',runname,'_coef_patterns.csv',sep=''),x=glmnet_coef_pattern,row.names=TRUE)
write.csv(file=paste('lars_',runname,'_coef_regions.csv',sep=''),x=lars_coef_region,row.names=TRUE)
write.csv(file=paste('glmnet_',runname,'_coef_regions.csv',sep=''),x=glmnet_coef_region,row.names=TRUE)
write.csv(file=paste('lars_',runname,'_coef_onlyregions.csv',sep=''),x=lars_coef_onlyregion,row.names=TRUE)
write.csv(file=paste('glmnet_',runname,'_coef_onlyregions.csv',sep=''),x=glmnet_coef_onlyregion,row.names=TRUE)

# make heatmaps
# p = coef.heatplot(lars_coef_full)
# ggsave(paste('lars_',runname,'_coef_full.jpeg',sep=''),plot=p,height=400,width=1400,units='mm',limitsize=FALSE)
# p = coef.heatplot(lars_coef_fullmax)
# ggsave(paste('lars_',runname,'_coef_fullmax.jpeg',sep=''),plot=p,height=400,width=1400,units='mm',limitsize=FALSE)
# p = coef.heatplot(lars_coef_braak)
# ggsave(paste('lars_',runname,'_coef_braak.jpeg',sep=''),plot=p,height=400,width=1000,units='mm',limitsize=FALSE)
# p = coef.heatplot(lars_coef_braakmax)
# ggsave(paste('lars_',runname,'_coef_braakmax.jpeg',sep=''),plot=p,height=400,width=1000,units='mm',limitsize=FALSE)
# p = coef.heatplot(lars_coef_pattern)
# ggsave(paste('lars_',runname,'_coef_patterns.jpeg',sep=''),plot=p,height=400,width=1000,units='mm',limitsize=FALSE)
# p = coef.heatplot(lars_coef_region)
# ggsave(paste('lars_',runname,'_coef_regions.jpeg',sep=''),plot=p,height=400,width=2000,units='mm',limitsize=FALSE)
# p = coef.heatplot(lars_coef_onlyregion)
# ggsave(paste('lars_',runname,'_coef_onlyregions.jpeg',sep=''),plot=p,height=400,width=1800,units='mm',limitsize=FALSE)

p = coef.heatplot(glmnet_coef_full)
ggsave(paste('glmnet_',runname,'_coef_full.jpeg',sep=''),plot=p,height=400,width=1400,units='mm',limitsize=FALSE)
p = coef.heatplot(glmnet_coef_fullmax)
ggsave(paste('glmnet_',runname,'_coef_fullmax.jpeg',sep=''),plot=p,height=400,width=1400,units='mm',limitsize=FALSE)
p = coef.heatplot(glmnet_coef_braak)
ggsave(paste('glmnet_',runname,'_coef_braak.jpeg',sep=''),plot=p,height=400,width=1000,units='mm',limitsize=FALSE)
p = coef.heatplot(glmnet_coef_braakmax)
ggsave(paste('glmnet_',runname,'_coef_braakmax.jpeg',sep=''),plot=p,height=400,width=1000,units='mm',limitsize=FALSE)
p = coef.heatplot(glmnet_coef_pattern)
ggsave(paste('glmnet_',runname,'_coef_patterns.jpeg',sep=''),plot=p,height=400,width=1000,units='mm',limitsize=FALSE)
p = coef.heatplot(glmnet_coef_region)
ggsave(paste('glmnet_',runname,'_coef_regions.jpeg',sep=''),plot=p,height=400,width=2000,units='mm',limitsize=FALSE)
p = coef.heatplot(glmnet_coef_onlyregion)
ggsave(paste('glmnet_',runname,'_coef_onlyregions.jpeg',sep=''),plot=p,height=400,width=1800,units='mm',limitsize=FALSE)
