df_av1451 = read.csv('nsfa/av1451skull_pattern_dataset.csv')

region_df = read.csv('datasets/pvc_adni_av1451/tauskullregions_uptake.csv')
region_df$timepoint = NULL
region_df$subject = as.factor(as.integer(region_df$subject))
non.na = complete.cases(df_av1451[,c('AV1451_BL_closest_AV45_wcereb_BIN1.11','Diag.AV1451')])
df_av1451 = df_av1451[non.na,]

# Regions to test
rois = c('CTX_LH_SUPERIORPARIETAL', 'LEFT_CAUDATE', 'CTX_RH_FUSIFORM', 'CTX_LH_ISTHMUSCINGULATE', 'CTX_RH_CAUDALANTERIORCINGULATE', 'RIGHT_PARSFR', 'HEMIWM', 'CTX_RH_PARACENTRAL', 'CTX_LH_PARAHIPPOCAMPAL', 'CTX_LH_PARACENTRAL', 'CTX_RH_INFERIORTEMPORAL', 'CTX_LH_LATERALOCCIPITAL', 'RIGHT_MIDDLEFR', 'CTX_LH_INFERIORPARIETAL', 'RIGHT_CAUDATE', 'LEFT_ORBITOFR', 'CTX_RH_LINGUAL', 'CTX_RH_ISTHMUSCINGULATE', 'RIGHT_ORBITOFR', 'CTX_RH_MIDDLETEMPORAL', 'CTX_LH_LINGUAL', 'LEFT_PALLIDUM', 'CTX_LH_CAUDALANTERIORCINGULATE', 'CHOROID', 'CEREBGM', 'CTX_LH_TEMPORALPOLE', 'LEFT_AMYGDALA', 'CTX_RH_PRECENTRAL', 'CTX_RH_ROSTRALANTERIORCINGULATE', 'CTX_LH_PRECENTRAL', 'CTX_LH_TRANSVERSETEMPORAL', 'CTX_RH_POSTERIORCINGULATE', 'CTX_LH_MIDDLETEMPORAL', 'CTX_RH_POSTCENTRAL', 'CTX_RH_PERICALCARINE', 'RIGHT_THALAMUS_PROPER', 'CTX_LH_PRECUNEUS', 'CTX_RH_TRANSVERSETEMPORAL', 'CTX_LH_BANKSSTS', 'CTX_LH_SUPERIORTEMPORAL', 'CTX_LH_FUSIFORM', 'CTX_LH_SUPERIORFRONTAL', 'CTX_LH_SUPRAMARGINAL', 'CTX_RH_SUPERIORFRONTAL', 'CTX_LH_POSTERIORCINGULATE', 'LEFT_ACCUMBENS_AREA', 'LEFT_PUTAMEN', 'CTX_LH_ROSTRALANTERIORCINGULATE', 'LEFT_MIDDLEFR', 'LEFT_HIPPOCAMPUS', 'RIGHT_PALLIDUM', 'LEFT_THALAMUS_PROPER', 'CTX_RH_PRECUNEUS', 'CTX_RH_PARAHIPPOCAMPAL', 'CTX_RH_SUPRAMARGINAL', 'CEREBWM', 'RIGHT_AMYGDALA', 'CTX_LH_CUNEUS', 'CTX_LH_POSTCENTRAL', 'CTX_RH_BANKSSTS', 'CTX_RH_ENTORHINAL', 'CTX_LH_PERICALCARINE', 'BRAIN_STEM', 'CTX_RH_TEMPORALPOLE', 'CTX_RH_CUNEUS', 'CTX_RH_INFERIORPARIETAL', 'RIGHT_PUTAMEN', 'CTX_LH_INFERIORTEMPORAL', 'LEFT_PARSFR', 'CTX_LH_INSULA', 'OTHER', 'CTX_LH_ENTORHINAL', 'CTX_RH_LATERALOCCIPITAL', 'CTX_RH_SUPERIORTEMPORAL', 'RIGHT_ACCUMBENS_AREA', 'RIGHT_HIPPOCAMPUS', 'CTX_RH_INSULA', 'CTX_RH_SUPERIORPARIETAL')
region_df[,rois] = region_df[,rois] / region_df$CEREBGM
# amyloid negative normal ADNI subjects
neg_normal = df_av1451[(df_av1451$Diag.AV1451 %in% c('N','SMC')) & (df_av1451$AV1451_BL_closest_AV45_wcereb_BIN1.11 == 0),]$RID
# amyloid positive LMCI/AD ADNI subjects
pos_ad = df_av1451[(df_av1451$Diag.AV1451 %in% c('AD','LMCI')) & (df_av1451$AV1451_BL_closest_AV45_wcereb_BIN1.11 == 1),]$RID

# 2 sample t test for each region
ttest_pvalues = data.frame()
for (roi in rois) {
  null_values = region_df[region_df$subject %in% neg_normal,roi]
  test_values = region_df[region_df$subject %in% pos_ad,roi]
  ftest = var.test(null_values,test_values) # F-test for homoskedasticity
  if (is.nan(ftest$p.value)) {
    ttest_pvalues[roi,'ttest_pvalue'] = NaN
    next
  }
  equal_var = ftest$p.value >= 0.05
  ttest = t.test(null_values,test_values,var.equal=FALSE,paired=FALSE)
  ttest_pvalues[roi,'ttest_pvalue'] = ttest$p.value
}
ttest_pvalues = na.omit(ttest_pvalues)

# Bonferroni correction
ttest_pvalues[,'ttest_pvalues_bonferroni'] = p.adjust(ttest_pvalues$ttest_pvalue,method='bonferroni')
ttest_pvalues[,'ttest_pvalues_holm'] = p.adjust(ttest_pvalues$ttest_pvalue,method='holm')
ttest_pvalues[,'ttest_pvalues_BH'] = p.adjust(ttest_pvalues$ttest_pvalue,method='BH')

write.csv(ttest_pvalues,file="negN_posAD_pvc_region_ttest.csv")
