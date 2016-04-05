library(lme4)
library(gdata)
library(psych)
library(reshape2)

all_content = readLines('FDG_AV45_COGdata/FDG_AV45_COGdata_03_23_16.csv')
skip_first = all_content[-1]
df = read.csv(textConnection(skip_first), header=TRUE, stringsAsFactors=FALSE)

# Convert to factors
to_factor = c('RID','APOE4_BIN','Diag.AV45_long','Gender')
df[to_factor] = lapply(df[to_factor], as.factor)

#time_col_prefix = 'TIMEpostAV45_ADAS'
#value_col_prefix = 'ADAScog'
time_col_prefix = 'TIMEpostAV45_AVLT.'
value_col_prefix = 'AVLT.'

# Keep relevant columns
demo_columns = c('RID','APOE4_BIN','Diag.AV45_long','Gender','Age.AV45','Edu..Yrs.')
av45_columns = c('AV45_PVC_CorticalSummary_WholeCereb_BL','AV45_TP_wcereb','AV45_NONTP_wcereb')
time_columns = Filter(function(i){startsWith(i,time_col_prefix)}, names(df))
value_columns = Filter(function(i){startsWith(i,value_col_prefix)}, names(df))
df = df[c(demo_columns,av45_columns,time_columns,value_columns)]

# Convert to long format
df_time_wide = df[c(demo_columns,av45_columns,time_columns)]
colnames(df_time_wide) = gsub(time_col_prefix,'TP',names(df_time_wide))
df_value_wide = df[c(demo_columns,av45_columns,value_columns)]
colnames(df_value_wide) = gsub(value_col_prefix,'TP',names(df_value_wide))
df_time_long = melt(df_time_wide, 
                    id.vars=c(demo_columns,av45_columns),
                    measure.vars=Filter(function(x){startsWith(x,'TP')},names(df_time_wide)),
                    variable.name='timepoint',
                    value.name='time')
df_value_long = melt(df_value_wide,
                     id.vars=c(demo_columns,av45_columns),
                     measure.vars=Filter(function(x){startsWith(x,'TP')},names(df_value_wide)),
                     variable.name='timepoint',
                     value.name='value')
merge_on = c(demo_columns,av45_columns,'timepoint')
df_long = merge(df_time_long,df_value_long,merge_on)
df_long = df_long[complete.cases(df_long[,names(df_long)]),]

# Separate Diag groups
df_long_n = df_long[which(df_long$Diag.AV45_long %in% c('N','SMC')),] 
df_long_emci = df_long[which(df_long$Diag.AV45_long %in% c('EMCI')),] 
df_long_lmci = df_long[which(df_long$Diag.AV45_long %in% c('LMCI')),] 
df_long_ad = df_long[which(df_long$Diag.AV45_long %in% c('AD')),] 

# Baseline AV45 vs Cog change LME models
fm_null_n = lmer(value ~ Age.AV45 + Gender + Edu..Yrs. + Diag.AV45_long + APOE4_BIN*time + (1 + time | RID), df_long_n)
fm_null_emci = lmer(value ~ Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN*time + (1 + time | RID), df_long_emci)
fm_null_lmci = lmer(value ~ Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN*time + (1 + time | RID), df_long_lmci)
fm_null_ad = lmer(value ~ Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN*time + (1 + time | RID), df_long_ad)

fm_tp_n = lmer(value ~ AV45_TP_wcereb + Age.AV45 + Gender + Edu..Yrs. + Diag.AV45_long + APOE4_BIN*time + (1 + time | RID), df_long_n)
fm_pvc_n = lmer(value ~ AV45_PVC_CorticalSummary_WholeCereb_BL + Age.AV45 + Gender + Edu..Yrs. + Diag.AV45_long + APOE4_BIN*time + (1 + time | RID), df_long_n)
fm_tp_emci = lmer(value ~ AV45_TP_wcereb + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN*time + (1 + time | RID), df_long_emci)
fm_pvc_emci = lmer(value ~ AV45_PVC_CorticalSummary_WholeCereb_BL + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN*time + (1 + time | RID), df_long_emci)
fm_tp_lmci = lmer(value ~ AV45_TP_wcereb + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN*time + (1 + time | RID), df_long_lmci)
fm_pvc_lmci = lmer(value ~ AV45_PVC_CorticalSummary_WholeCereb_BL + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN*time + (1 + time | RID), df_long_lmci)
fm_tp_ad = lmer(value ~ AV45_TP_wcereb + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN*time + (1 + time | RID), df_long_ad)
fm_pvc_ad = lmer(value ~ AV45_PVC_CorticalSummary_WholeCereb_BL + Age.AV45 + Gender + Edu..Yrs. + APOE4_BIN*time + (1 + time | RID), df_long_ad)

fm_tp_n_anova = anova(fm_null_n,fm_tp_n)
fm_pvc_n_anova = anova(fm_null_n,fm_pvc_n)
fm_tp_emci_anova = anova(fm_null_emci,fm_tp_emci)
fm_pvc_emci_anova = anova(fm_null_emci,fm_pvc_emci)
fm_tp_lmci_anova = anova(fm_null_lmci,fm_tp_lmci)
fm_pvc_lmci_anova = anova(fm_null_lmci,fm_pvc_lmci)
fm_tp_ad_anova = anova(fm_null_ad,fm_tp_ad)
fm_pvc_ad_anova = anova(fm_null_ad,fm_pvc_ad)

fm_tp_n_anova$`Pr(>Chisq)`[2]
fm_pvc_n_anova$`Pr(>Chisq)`[2]
fm_tp_emci_anova$`Pr(>Chisq)`[2]
fm_pvc_emci_anova$`Pr(>Chisq)`[2]
fm_tp_lmci_anova$`Pr(>Chisq)`[2]
fm_pvc_lmci_anova$`Pr(>Chisq)`[2]
fm_tp_ad_anova$`Pr(>Chisq)`[2]
fm_pvc_ad_anova$`Pr(>Chisq)`[2]

Age.AV45
Gender
Edu..Yrs.

AV45_PVC_CorticalSummary_WholeCereb_BL	AV45_PVC_CorticalSummary_WholeCereb_Scan2	AV45_PVC_CorticalSummary_WholeCereb_Scan3
AV45_TP_wcereb	AV45_TP_2_wcereb	AV45_TP_3_wcereb
AV45_NONTP_wcereb	AV45_NONTP_2_wcereb	AV45_NONTP_3_wcereb


ADAScog.1
ADAScog.2	
ADAScog.3	
ADAScog.4	
ADAScog.5	
ADAScog.6	
ADAScog.7
ADAScog.8	ADAScog.9	ADAScog.10	ADAScog.11	
TIME_ADAS.1	TIME_ADAS.2	TIME_ADAS.3	TIME_ADAS.4	TIME_ADAS.5	TIME_ADAS.6	TIME_ADAS.7	TIME_ADAS.8	TIME_ADAS.9	TIME_ADAS.10	TIME_ADAS.11	
TIMEreltoAV45_ADAS.1	TIMEreltoAV45_ADAS.2	TIMEreltoAV45_ADAS.3	TIMEreltoAV45_ADAS.4	TIMEreltoAV45_ADAS.5	TIMEreltoAV45_ADAS.6	TIMEreltoAV45_ADAS.7	TIMEreltoAV45_ADAS.8	TIMEreltoAV45_ADAS.9	TIMEreltoAV45_ADAS.10	TIMEreltoAV45_ADAS.11	
TIMEpostAV45_ADAS.1	TIMEpostAV45_ADAS.2	TIMEpostAV45_ADAS.3	TIMEpostAV45_ADAS.4	TIMEpostAV45_ADAS.5	TIMEpostAV45_ADAS.6	TIMEpostAV45_ADAS.7	TIMEpostAV45_ADAS.8	TIMEpostAV45_ADAS.9	TIMEpostAV45_ADAS.10	TIMEpostAV45_ADAS.11