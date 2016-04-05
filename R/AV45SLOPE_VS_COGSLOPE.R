library(lm4)
library(gdata)
library(psych)
library(reshape2)
library(plyr)

all_content = readLines('FDG_AV45_COGdata/FDG_AV45_COGdata_04_04_16.csv')
skip_first = all_content[-1]
df = read.csv(textConnection(skip_first), header=TRUE, stringsAsFactors=FALSE)

# Convert to factors
to_factor = c('RID',
              'APOE4_BIN',
              'Diag.AV45_long',
              'Gender',
              'AV45_TP_wcereb_BIN1.11',
              'AV45_PVC_CorticalSummary_WholeCereb_1.15_BL')
df[to_factor] = lapply(df[to_factor], as.factor)

value_col = 'ADASslope_postAV45'
#value_col = 'UW_MEM_slope'
#value_col = 'AVLTslope_postAV45'

# Keep relevant columns
demo_columns = c('RID','APOE4_BIN','Diag.AV45_long','Gender','Age.AV45','Edu..Yrs.')
av45_columns = c('AV45_PVC_CorticalSummary_WholeCereb_BL',
                 'AV45_TP_wcereb',
                 'AV45_NONTP_wcereb',
                 'AV45_TP_wcereb_BIN1.11',
                 'AV45_TP_wcereb_Slope_2pts',
                 'AV45_PVC_CorticalSummary_WholeCereb_1.15_BL',
                 'AV45_PVC_CorticalSummary_WholeCereb_slope_2points')
df = df[c(demo_columns,av45_columns,value_col)]
df = df[complete.cases(df[,names(df)]),]
names(df)[names(df)==value_col] <- "value"

# Filter by baseline positivity
df = df[which(df$'AV45_TP_wcereb_BIN1.11' == 1),] 
df = df[which(df$'AV45_PVC_CorticalSummary_WholeCereb_1.15_BL' == 1),]

# Separate Diag groups
df_n = df[which(df$Diag.AV45_long %in% c('N','SMC')),] 
df_emci = df[which(df$Diag.AV45_long %in% c('EMCI')),] 
df_lmci = df[which(df$Diag.AV45_long %in% c('LMCI')),] 
df_ad = df[which(df$Diag.AV45_long %in% c('AD')),] 

# Baseline AV45 vs Cog change lm models
fm_null_tp_n = lm(value ~ AV45_TP_wcereb + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_n)
fm_null_tp_emci = lm(value ~ AV45_TP_wcereb + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_emci)
fm_null_tp_lmci = lm(value ~ AV45_TP_wcereb + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_lmci)
fm_null_tp_ad = lm(value ~ AV45_TP_wcereb + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_ad)

fm_null_pvc_n = lm(value ~ AV45_PVC_CorticalSummary_WholeCereb_BL + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_n)
fm_null_pvc_emci = lm(value ~ AV45_PVC_CorticalSummary_WholeCereb_BL + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_emci)
fm_null_pvc_lmci = lm(value ~ AV45_PVC_CorticalSummary_WholeCereb_BL + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_lmci)
fm_null_pvc_ad = lm(value ~ AV45_PVC_CorticalSummary_WholeCereb_BL + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_ad)

fm_tp_n = lm(value ~ AV45_TP_wcereb + AV45_TP_wcereb_Slope_2pts + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_n)
fm_pvc_n = lm(value ~ AV45_PVC_CorticalSummary_WholeCereb_BL + AV45_PVC_CorticalSummary_WholeCereb_slope_2points + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_n)
fm_tp_emci = lm(value ~ AV45_TP_wcereb + AV45_TP_wcereb_Slope_2pts + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_emci)
fm_pvc_emci = lm(value ~ AV45_PVC_CorticalSummary_WholeCereb_BL + AV45_PVC_CorticalSummary_WholeCereb_slope_2points + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_emci)
fm_tp_lmci = lm(value ~ AV45_TP_wcereb + AV45_TP_wcereb_Slope_2pts + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_lmci)
fm_pvc_lmci = lm(value ~ AV45_PVC_CorticalSummary_WholeCereb_BL + AV45_PVC_CorticalSummary_WholeCereb_slope_2points + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_lmci)
fm_tp_ad = lm(value ~ AV45_TP_wcereb + AV45_TP_wcereb_Slope_2pts + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_ad)
fm_pvc_ad = lm(value ~ AV45_PVC_CorticalSummary_WholeCereb_BL + AV45_PVC_CorticalSummary_WholeCereb_slope_2points + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN, df_ad)

fm_tp_n_anova = anova(fm_null_tp_n,fm_tp_n)
fm_pvc_n_anova = anova(fm_null_pvc_n,fm_pvc_n)
fm_tp_emci_anova = anova(fm_null_tp_emci,fm_tp_emci)
fm_pvc_emci_anova = anova(fm_null_pvc_emci,fm_pvc_emci)
fm_tp_lmci_anova = anova(fm_null_tp_lmci,fm_tp_lmci)
fm_pvc_lmci_anova = anova(fm_null_pvc_lmci,fm_pvc_lmci)
fm_tp_ad_anova = anova(fm_null_tp_ad,fm_tp_ad)
fm_pvc_ad_anova = anova(fm_null_pvc_ad,fm_pvc_ad)

fm_tp_n_anova$`Pr(>F)`[2]
fm_pvc_n_anova$`Pr(>F)`[2]
fm_tp_emci_anova$`Pr(>F)`[2]
fm_pvc_emci_anova$`Pr(>F)`[2]
fm_tp_lmci_anova$`Pr(>F)`[2]
fm_pvc_lmci_anova$`Pr(>F)`[2]
fm_tp_ad_anova$`Pr(>F)`[2]
fm_pvc_ad_anova$`Pr(>F)`[2]

