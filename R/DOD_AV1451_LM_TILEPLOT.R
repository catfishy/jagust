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
to_factor = c('SCRNO','APOE4BIN','Sex','PTGroup',
              'Diag_closest_AV45_BL','Diag_closest_AV1451_BL')
to_standardize = c('Age','Edu')
diag_columns = c('Diag.AV45','Diag.AV1451')
braak_columns = c('AV1451_PVC_Braak1_CerebGray_BL',
                  'AV1451_PVC_Braak2_CerebGray_BL',
                  'AV1451_PVC_Braak3_CerebGray_BL',
                  'AV1451_PVC_Braak4_CerebGray_BL',
                  'AV1451_PVC_Braak5_CerebGray_BL',
                  'AV1451_PVC_Braak6_CerebGray_BL',
                  'AV1451_PVC_Braak12_CerebGray_BL',
                  'AV1451_PVC_Braak34_CerebGray_BL',
                  'AV1451_PVC_Braak56_CerebGray_BL',
                  'AV1451_PVC_BraakAll_CerebGray_BL')
valid_diags = c('N','SMC','EMCI','LMCI','AD')

# Choose demographic covariates
# demog_columns = c('SCRNO','APOE4BIN','Sex','Age','Edu')
demog_columns = c('SCRNO','Sex','Age','Edu')

# IMPORT
df_av1451 = read.csv('DOD_DATA/DOD_DATA_11_03_16.csv')
df_regions = read.csv('pvc/pvc_dod_av1451/tauskullregions_suvr_bilateral.csv')
df_regions = rename(df_regions,c("subject"="SCRNO"))
all_regions = colnames(df_regions)
all_regions = all_regions[!(all_regions %in% c('SCRNO','CEREBGM'))]
df_av1451 = merge(df_av1451,df_regions,by='SCRNO')

# Convert variables
non.na = complete.cases(df_av1451[,c(demog_columns,braak_columns)])
df_av1451 = df_av1451[non.na,]
for (i in names(df_av1451)){
  if (i %in% to_factor){
    df_av1451[,eval(i)] = as.factor(as.character(df_av1451[,eval(i)]))
  }
}

# Filter by patient group
valid_ptgroup = c('PTSD','TBI+PTSD')
df_av1451 = df_av1451[which(df_av1451$PTGroup %in% valid_ptgroup),]

# Standardize variables
cross_to_standardize = c(to_standardize,
                         braak_columns,
                         all_regions)
cross_normalization = preProcess(df_av1451[,cross_to_standardize])
df_av1451[,cross_to_standardize] = predict(cross_normalization, df_av1451[,cross_to_standardize])


indep_variables = c(braak_columns,
                    all_regions)
dep_variables = c("AVLT_closest_AV1451_BL",
                  "CAPS_LIFETIME_SCORE",
                  "AV1451_1_closest_AV45_wcereb")

# Standardize targets
target_norm = preProcess(df_av1451[,dep_variables])
df_av1451[,dep_variables] = predict(target_norm, df_av1451[,dep_variables])

results = data.frame()
counter = 0
for (dep_var in dep_variables) {
  for (indep_var in indep_variables) {
    lm_form = paste(dep_var,'~ Age + Edu +',indep_var,collapse=' ')
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
    results[counter_str,'points'] = sum(complete.cases(df_av1451[dep_var]))
    results[counter_str,'r2'] = r2
    results[counter_str,'coef_estimate'] = coef
    results[counter_str,'fstat_pvalue'] = fstat_pvalue
    results[counter_str,'tstat'] = tstat
    results[counter_str,'tstat_pvalue'] = tstat_pvalue
    counter = counter + 1
  }
}

# row order
results$dep = factor(results$dep, levels=dep_variables)
results$indep = factor(results$indep, levels=indep_variables)

# bonferroni correct
results_uncorrected = data.frame(results)
correct_method = 'bonferroni'
# correct_method = 'holm'

for (dep in dep_variables) {
  results[results$dep == dep,'fstat_pvalue'] = p.adjust(results[results$dep == dep,'fstat_pvalue'],method=correct_method)
  results[results$dep == dep,'tstat_pvalue'] = p.adjust(results[results$dep == dep,'tstat_pvalue'],method=correct_method)
}

# create tile plots
base_size = 12

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

