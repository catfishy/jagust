library(ggplot2)
library(reshape)

## IMPORT DATA ##

# df = read.csv('/usr/local/jagust/AV1451_template_space_hist.csv')
# colnames(df) = gsub("X","",colnames(df))
# colnames(df) = gsub("\\.\\.\\."," - ",colnames(df))

# df = read.csv('/usr/local/jagust/AV1451_template_space_cumulative.csv')
# colnames(df) = gsub("X..",">=",colnames(df))

df_regions = read.csv('datasets/pvc_adni_av1451/tauskullregions_uptake_bilateral.csv')
df_regions = rename(df_regions,c("subject"="RID"))
rownames(df_regions) = df_regions$RID
df_regions$RID = NULL
all_regions = colnames(df_regions)
all_regions = all_regions[!(all_regions %in% c('RID','CEREBGM'))]
df_regions = df_regions[,all_regions]
thresholds = seq(0,4,by=0.1)
df = NULL
for (thresh in thresholds) {
  greater_df = data.frame(rowSums(df_regions >= thresh)/length(all_regions))
  colnames(greater_df) = paste('>=',format(thresh,nsmall=1),sep='')
  if (is.null(df)) {
    df = greater_df
  } else {
    df = cbind(df,greater_df)
  }
}
df$SID = factor(rownames(df))
df$Dataset = 'ADNI'

################

df = rename(df,c("SID"="RID"))
master_df = read.csv('/usr/local/jagust/FDG_AV45_COGdata/FDG_AV45_COGdata_09_19_16.csv', skip=1)
master_df$RID = as.factor(master_df$RID)
master_df = master_df[c('RID','Diag.AV1451','AV1451_BL_closest_AV45_wcereb_BIN1.11')]
valid_diags = c('N','SMC','EMCI','LMCI','AD')
master_df = master_df[master_df$Diag.AV1451 %in% valid_diags,]
master_df$Diag.AV1451 = factor(master_df$Diag.AV1451,levels=valid_diags)

df_adni = df[df$Dataset == 'ADNI',]
df_bacs = df[df$Dataset == 'BACS',]
measures = names(df)
measures = measures[!(measures %in% c('RID','Dataset'))]
df_bacs.melt = melt(df_bacs,id.vars=c('RID'),measure.vars=measures)
df_adni.melt = melt(df_adni,id.vars=c('RID'),measure.vars=measures)
df_adni.melt = merge(df_adni.melt,master_df,by='RID')

ggplot(df_adni.melt, aes(x=variable,y=value, group=Diag.AV1451, color=Diag.AV1451)) + 
  stat_summary(geom='pointrange') +
  stat_summary(geom='line') + 
  ggtitle("ADNI % Voxels above threshold") + 
  theme(panel.border=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0))

measures.imp = c(">=1.0",">=1.1",">=1.2",">=1.3",">=1.4",">=1.5",">=1.6",">=1.7")
for (cur_measure in measures.imp) {
  cur_df = df_adni.melt[df_adni.melt$variable == cur_measure,]
  ggplot(cur_df, aes(x=Diag.AV1451,y=value,color=Diag.AV1451)) +
    geom_point() + 
    stat_summary(geom='crossbar') +
    ggtitle(paste("% Voxels",cur_measure,"by DX")) + 
    ylab(paste("% Voxels",cur_measure))
  ggsave(paste("Voxels",cur_measure,".pdf",sep=''))
}
