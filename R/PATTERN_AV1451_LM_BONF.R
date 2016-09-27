library(lme4)
library(coefplot2)
library(ggplot2)
library(stats)
library(reshape2)
library(caret)
library(lmtest)
library(lars)
library(plyr)
library(scales)

source('R/LM_FUNCS.R')

# CONSTANTS
pattern_prefix = 'NSFA_'
to_factor = c('RID','APOE4_BIN','APOE2_BIN','Gender',
              'Diag.AV45','Diag.AV1451',
              'AV1451_BL_closest_AV45_wcereb_BIN1.11',
              'AV45_NONTP_1_wcereb_BIN1.11',
              'AV45_NONTP_2_wcereb_BIN1.11',
              'AV45_NONTP_3_wcereb_BIN1.11')
to_standardize = c('Age.AV45','Edu..Yrs.','Age.AV1451')
demog_columns = c('RID','APOE4_BIN','Diag.AV1451','Age.AV1451','Gender','Edu..Yrs.')
diag_columns = c('Diag.AV45','Diag.AV1451')
braak_columns = c('AV1451_PVC_Braak1_CerebGray_BL',
                  'AV1451_PVC_Braak2_CerebGray_BL',
                  'AV1451_PVC_Braak3_CerebGray_BL',
                  'AV1451_PVC_Braak4_CerebGray_BL',
                  'AV1451_PVC_Braak5_CerebGray_BL',
                  'AV1451_PVC_Braak6_CerebGray_BL',
                  'AV1451_PVC_BraakAll_CerebGray_BL')
braakmax_columns = c('BRAAK1_MAX','BRAAK2_MAX','BRAAK3_MAX',
                     'BRAAK4_MAX','BRAAK5_MAX','BRAAK6_MAX')
mask_columns = c('BRAAKALL_ROI_MASK_MEAN',
                 'BRAAKALL_ROI_MASK_MAX',
                 'AV45_POS_VS_NEG_VOX_MASK_MEAN',
                 'AV45_POS_VS_NEG_VOX_MASK_MAX')

output_folder = 'R/output_av1451/'

valid_diags = c('N','SMC','EMCI','LMCI','AD')

# IMPORT
df_av1451 = read.csv('nsfa/av1451skull_pattern_dataset.csv')
df_regions = read.csv('datasets/pvc_adni_av1451/tauskullregions_uptake_bilateral.csv')
df_regions = rename(df_regions,c("subject"="RID"))
df_regions$RID = as.factor(df_regions$RID)
all_regions = colnames(df_regions)
all_regions = all_regions[!(all_regions %in% c('RID','CEREBGM'))]
df_av1451 = merge(df_av1451,df_regions,by='RID')

df_max = read.csv('output/08-31-2016/UCBERKELEYAV1451_MAX_08-31-2016_regular_tp.csv')
max_cols = c('RID','BRAAK1','BRAAK2','BRAAK3','BRAAK4','BRAAK5','BRAAK6')
df_max = df_max[,max_cols]
colnames(df_max) = c('RID','BRAAK1_MAX','BRAAK2_MAX','BRAAK3_MAX','BRAAK4_MAX','BRAAK5_MAX','BRAAK6_MAX')
df_av1451 = merge(df_av1451,df_max,by='RID')

df_mask = read.csv('maskdata_inorm_cerebgm_snorm.csv')
df_mask$RID = as.factor(df_mask$RID)
df_av1451 = merge(df_av1451,df_mask,by='RID')

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

# df_av1451 = df_av1451_pos
# runname = 'av45_pos'

# df_av1451 = df_av1451_neg
# runname = 'av45_neg'

df_av1451 = df_av1451
runname = 'all_subj'

# Standardize variables
cross_to_standardize = c(to_standardize,pattern_columns,
                         braak_columns,braakmax_columns,all_regions,
                         mask_columns)
cross_normalization = preProcess(df_av1451[,cross_to_standardize])
df_av1451[,cross_to_standardize] = predict(cross_normalization, df_av1451[,cross_to_standardize])


indep_variables = c(braak_columns,
                    braakmax_columns,
                    mask_columns,
                    pattern_columns,
                    all_regions)
dep_variables = c("AVLT_AV1451_1",
                  "AVLT_retroslope_AV1451_BL",
                  "ADAS_AV1451_1",
                  "ADAS_retroslope_AV1451_BL",
                  'AV1451_BL_closest_AV45_wcereb',
                  "AV1451_BL_closest_AV45_wcereb_retroSlope",
                  "CSF_TAU_closest_AV1451_1",
                  "CSF_PTAU_closest_AV1451_1",
                  "CSF_ABETA_closest_AV1451_1",
                  "Age.AV1451")
av45_dep_variables = c("AVLT_AV1451_1",
                       "AVLT_retroslope_AV1451_BL",
                       "ADAS_AV1451_1",
                       "ADAS_retroslope_AV1451_BL",
                       "CSF_TAU_closest_AV1451_1",
                       "CSF_PTAU_closest_AV1451_1",
                       "CSF_ABETA_closest_AV1451_1",
                       "Age.AV1451")
# Standardize targets
target_norm = preProcess(df_av1451[,dep_variables])
df_av1451[,dep_variables] = predict(target_norm, df_av1451[,dep_variables])

results = data.frame()
counter = 0
for (dep_var in dep_variables) {
  for (indep_var in indep_variables) {
    lm_form = paste(dep_var,'~',indep_var,collapse=' ')
    fm = lm(lm_form,df_av1451)
    fm.summary = summary(fm)
    r2 = fm.summary$adj.r.squared
    fstat = fm.summary$fstatistic
    fstat_pvalue = pf(fstat[['value']],fstat[['numdf']],fstat[['dendf']],lower.tail=F)
    coef = fm.summary$coefficients[indep_var,'Estimate']
    tstat = fm.summary$coefficients[indep_var,'t value']
    tstat_pvalue = fm.summary$coefficients[indep_var,"Pr(>|t|)"]
    counter_str = toString(counter)
    results[counter_str,'dep'] = dep_var
    results[counter_str,'indep'] = indep_var
    results[counter_str,'r2'] = r2
    results[counter_str,'coef_estimate'] = coef
    results[counter_str,'fstat_pvalue'] = fstat_pvalue
    results[counter_str,'tstat'] = tstat
    results[counter_str,'tstat_pvalue'] = tstat_pvalue
    counter = counter + 1
  }
}

results_av45 = data.frame()
counter = 0
for (dep_var in av45_dep_variables) {
  for (indep_var in indep_variables) {
    lm_form = paste(dep_var,' ~ AV1451_BL_closest_AV45_wcereb_BIN1.11*',indep_var,sep='')
    fm = lm(lm_form,df_av1451)
    fm.summary = summary(fm)
    r2 = fm.summary$adj.r.squared
    fstat = fm.summary$fstatistic
    fstat_pvalue = pf(fstat[['value']],fstat[['numdf']],fstat[['dendf']],lower.tail=F)
    var_name = paste('AV1451_BL_closest_AV45_wcereb_BIN1.111:',indep_var,sep='')
    coef = fm.summary$coefficients[var_name,'Estimate']
    tstat = fm.summary$coefficients[var_name,'t value']
    tstat_pvalue = fm.summary$coefficients[var_name,"Pr(>|t|)"]
    counter_str = toString(counter)
    results_av45[counter_str,'dep'] = dep_var
    results_av45[counter_str,'indep'] = indep_var
    results_av45[counter_str,'r2'] = r2
    results_av45[counter_str,'coef_estimate'] = coef
    results_av45[counter_str,'fstat_pvalue'] = fstat_pvalue
    results_av45[counter_str,'tstat'] = tstat
    results_av45[counter_str,'tstat_pvalue'] = tstat_pvalue
    counter = counter + 1
  }
}

results_firstav45 = data.frame()
counter = 0
for (dep_var in av45_dep_variables) {
  for (indep_var in indep_variables) {
    lm_form = paste(dep_var,' ~ AV45_NONTP_1_wcereb_BIN1.11*',indep_var,sep='')
    fm = lm(lm_form,df_av1451)
    fm.summary = summary(fm)
    r2 = fm.summary$adj.r.squared
    fstat = fm.summary$fstatistic
    fstat_pvalue = pf(fstat[['value']],fstat[['numdf']],fstat[['dendf']],lower.tail=F)
    var_name = paste('AV45_NONTP_1_wcereb_BIN1.111:',indep_var,sep='')
    coef = fm.summary$coefficients[var_name,'Estimate']
    tstat = fm.summary$coefficients[var_name,'t value']
    tstat_pvalue = fm.summary$coefficients[var_name,"Pr(>|t|)"]
    counter_str = toString(counter)
    results_firstav45[counter_str,'dep'] = dep_var
    results_firstav45[counter_str,'indep'] = indep_var
    results_firstav45[counter_str,'r2'] = r2
    results_firstav45[counter_str,'coef_estimate'] = coef
    results_firstav45[counter_str,'fstat_pvalue'] = fstat_pvalue
    results_firstav45[counter_str,'tstat'] = tstat
    results_firstav45[counter_str,'tstat_pvalue'] = tstat_pvalue
    counter = counter + 1
  }
}

# row order
results$dep = factor(results$dep, levels=dep_variables)
results$indep = factor(results$indep, levels=indep_variables)
results_av45$dep = factor(results_av45$dep, levels=av45_dep_variables)
results_av45$indep = factor(results_av45$indep, levels=indep_variables)
results_firstav45$dep = factor(results_firstav45$dep, levels=av45_dep_variables)
results_firstav45$indep = factor(results_firstav45$indep, levels=indep_variables)

# bonferroni correct
results_uncorrected = data.frame(results)
results_av45_uncorrected = data.frame(results_av45)
results_firstav45_uncorrected = data.frame(results_firstav45)
correct_method = 'bonferroni'
# correct_method = 'holm'

for (dep in dep_variables) {
  results[results$dep == dep,'fstat_pvalue'] = p.adjust(results[results$dep == dep,'fstat_pvalue'],method=correct_method)
  results[results$dep == dep,'tstat_pvalue'] = p.adjust(results[results$dep == dep,'tstat_pvalue'],method=correct_method)
}
for (dep in av45_dep_variables) {
  results_av45[results_av45$dep == dep,'fstat_pvalue'] = p.adjust(results_av45[results_av45$dep == dep,'fstat_pvalue'],method=correct_method)
  results_av45[results_av45$dep == dep,'tstat_pvalue'] = p.adjust(results_av45[results_av45$dep == dep,'tstat_pvalue'],method=correct_method)
}
for (dep in av45_dep_variables) {
  results_firstav45[results_firstav45$dep == dep,'fstat_pvalue'] = p.adjust(results_firstav45[results_firstav45$dep == dep,'fstat_pvalue'],method=correct_method)
  results_firstav45[results_firstav45$dep == dep,'tstat_pvalue'] = p.adjust(results_firstav45[results_firstav45$dep == dep,'tstat_pvalue'],method=correct_method)
}

# create tile plots
base_size = 12

ggplot(results_firstav45_uncorrected, aes(indep, dep)) + 
  geom_tile(aes(fill = tstat_pvalue), 
            colour='grey50') + 
  scale_fill_distiller(palette = "Spectral",
                       limits = c(0,0.1),
                       oob = squish) + 
  labs(x="",y="") + 
  ggtitle("AV1451|Early_AV45+, Pr(>|t|), uncorrected") + 
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

ggplot(results_firstav45_uncorrected, aes(indep, dep)) + 
  geom_tile(aes(fill = coef_estimate), 
            colour='grey50') + 
  scale_fill_distiller(palette = "Spectral",
                       limits = c(-1,1),
                       oob = squish) + 
  labs(x="",y="") + 
  ggtitle("AV1451|Early_AV45+, Coef Est.") + 
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

ggplot(results_av45_uncorrected, aes(indep, dep)) + 
  geom_tile(aes(fill = tstat_pvalue), 
            colour='grey50') + 
  scale_fill_distiller(palette = "Spectral",
                       limits = c(0,0.1),
                       oob = squish) + 
  labs(x="",y="") + 
  ggtitle("AV1451|AV45+, Pr(>|t|), uncorrected") + 
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

ggplot(results_av45_uncorrected, aes(indep, dep)) + 
  geom_tile(aes(fill = coef_estimate), 
            colour='grey50') + 
  scale_fill_distiller(palette = "Spectral",
                       limits = c(-1,1),
                       oob = squish) + 
  labs(x="",y="") + 
  ggtitle("AV1451|AV45+, Coef Est.") + 
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


ggplot(results_uncorrected, aes(indep, dep)) + 
  geom_tile(aes(fill = tstat_pvalue), 
            colour='grey50') + 
  scale_fill_distiller(palette = "Spectral",
                       trans = "log",
                       breaks=c(1e-1,1e-4,1e-7,1e-10),
                       labels=c("1e-1","1e-4","1e-7","1e-10")) + 
  labs(x="",y="") + 
  ggtitle("Pr(>|t|), uncorrected") + 
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

ggplot(results_uncorrected, aes(indep, dep)) + 
  geom_tile(aes(fill = fstat_pvalue), 
            colour='grey50') + 
  scale_fill_distiller(palette = "Spectral",
                       trans = "log",
                       breaks=c(1e-1,1e-4,1e-7,1e-10),
                       labels=c("1e-1","1e-4","1e-7","1e-10")) + 
  labs(x="",y="") + 
  ggtitle("F-statistic p-value, uncorrected") + 
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


ggplot(results, aes(indep, dep)) + 
  geom_tile(aes(fill = tstat_pvalue), 
            colour='grey50') + 
  scale_fill_distiller(palette = "Spectral",
                       trans = "log",
                       breaks=c(1e-1,1e-4,1e-7,1e-10),
                       labels=c("1e-1","1e-4","1e-7","1e-10")) + 
  labs(x="",y="") + 
  ggtitle("Pr(>|t|), bonferroni corrected") + 
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

ggplot(results, aes(indep, dep)) + 
  geom_tile(aes(fill = fstat_pvalue), 
            colour='grey50') + 
  scale_fill_distiller(palette = "Spectral",
                       trans = "log",
                       breaks=c(1e-1,1e-4,1e-7,1e-10),
                       labels=c("1e-1","1e-4","1e-7","1e-10")) + 
  labs(x="",y="") + 
  ggtitle("F-statistic p-value, bonferroni corrected") + 
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

ggplot(results, aes(indep, dep)) + 
  geom_tile(aes(fill = r2), 
            colour='grey50') + 
  scale_fill_distiller(palette = "Spectral") + 
  labs(x="",y="") + 
  ggtitle("Adjusted R-Squared") + 
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

