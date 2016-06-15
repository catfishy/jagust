

# df_rois = read.csv('nsfa/av1451uni_roi_comparisons.csv')
df_rois = read.csv('nsfa/av45_roi_comparisons.csv')

df_rois.ordered = df_rois[order(-df_rois$varperc),]
best_factors = df_rois.ordered[df_rois.ordered$varperc > 0.05,"FACTOR"]
region_order = c('FRONTAL','PARIETAL','SENSORY','TEMPORAL','OCCIPITAL',
                 'MEDIALOCCIPITAL','CINGULATE','BASALGANGLIA','LIMBIC',
                 'THALAMUS','CEREBRAL_WHITE','CEREBELLUM_GRAY','CEREBELLUM_WHITE')
braak_order = c('BRAAK1','BRAAK2','BRAAK3','BRAAK4','BRAAK5','BRAAK6')

# plot sum of squared loadings
p = ggplot(df_rois,aes(x=FACTOR,y=varperc)) +
  geom_point()
print(p)




# plot regional loading lines
df_rois.long = melt(df_rois[df_rois$FACTOR %in% best_factors,],id.vars='FACTOR',measure.vars=region_order)
p = ggplot(df_rois.long,aes_string(x='variable',y='value', color='FACTOR')) +
  geom_line(aes_string(group='FACTOR')) +
  geom_point() + 
  theme(axis.ticks=element_line(), 
        plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14, angle=30), 
        axis.text.x=element_text(angle=300, face='bold', size=14, hjust=0, vjust=1, color='grey50')) + 
  geom_hline(aes(yintercept=0.0)) + 
  ggtitle('Regional Loadings, by Factor') +
  xlab('Brain regions') +
  ylab('Sum of loading coefficients')
print(p)

# print braak loading lines
df_rois.braak = melt(df_rois[df_rois$FACTOR %in% best_factors,],id.vars='FACTOR',measure.vars=braak_order)
p = ggplot(df_rois.braak,aes_string(x='variable',y='value', color='FACTOR')) +
  geom_line(aes_string(group='FACTOR')) +
  geom_point() + 
  theme(axis.ticks=element_line(), 
        plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14, angle=30), 
        axis.text.x=element_text(angle=300, face='bold', size=14, hjust=0, vjust=1, color='grey50')) + 
  geom_hline(aes(yintercept=0.0)) + 
  ggtitle('Braak Region Loadings, by Factor') +
  xlab('Non-cumulative Braak regions') +
  ylab('Sum of loading coefficients')
print(p)
