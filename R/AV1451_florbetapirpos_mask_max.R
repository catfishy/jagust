library(ggplot2)

master_csv = 'FDG_AV45_COGdata/FDG_AV45_COGdata_08_23_16.csv'
mask_val_csv = 'florbetapirpos_vs_neg_mask_val.csv'

df_master = read.csv(master_csv,skip=1)
columns = c('RID',
            'AV1451_CerebGray_BL',
            'AV1451_BL_closest_AV45_wcereb_BIN1.11',
            'AV1451_BL_closest_AV45_wcereb',
            'AV1451_BL_closest_AV45_wcereb_retroSlope',
            'AV1451_Braak12_CerebGray_BL',
            'AV1451_Braak34_CerebGray_BL',
            'AV1451_Braak56_CerebGray_BL',
            'Diag.AV1451')
df_master = df_master[,columns]
df_mask_max = read.csv(mask_val_csv)
df = merge(df_master,df_mask_max,by='RID')

to_factor = c('RID','Diag.AV1451','AV1451_BL_closest_AV45_wcereb_BIN1.11')
for (i in names(df)){
  if (i %in% to_factor){
    df[,eval(i)] = as.factor(as.character(df[,eval(i)]))
  }
}
nonna = complete.cases(df[,c('AV1451_BL_closest_AV45_wcereb_BIN1.11')])
df = df[nonna,]

# plot against AV45 SUVR
plot(df$MASK_MAX,df$AV1451_BL_closest_AV45_wcereb)
plot(df$MASK_MEAN,df$AV1451_BL_closest_AV45_wcereb)

# plot against AV45 retro slope
plot(df$MASK_MAX,df$AV1451_BL_closest_AV45_wcereb_retroSlope)
plot(df$MASK_MEAN,df$AV1451_BL_closest_AV45_wcereb_retroSlope)

# violin plot by AV45 positivity
p1 = ggplot(df, aes_string('AV1451_BL_closest_AV45_wcereb_BIN1.11','MASK_MAX')) +
  geom_violin(trim=TRUE,
              na.rm=TRUE,
              aes(fill=AV1451_BL_closest_AV45_wcereb_BIN1.11),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab('MASK Max SUVR') +
  xlab('AV45 Positivity')
print(p1)

p2 = ggplot(df, aes_string('AV1451_BL_closest_AV45_wcereb_BIN1.11','MASK_MEAN')) +
  geom_violin(trim=TRUE,
              na.rm=TRUE,
              aes(fill=AV1451_BL_closest_AV45_wcereb_BIN1.11),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab('MASK Mean SUVR') +
  xlab('AV45 Positivity')
print(p2)
