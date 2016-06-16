

df_rois = read.csv('nsfa/av1451uni_roi_comparisons.csv')
# df_rois = read.csv('nsfa/av45_roi_comparisons.csv')

df_rois$ASYM_ABS = rowSums(abs(df_rois[, asym_order]))
best_asym = df_rois[order(-df_rois$ASYM_ABS),][1:5,"FACTOR"]
best_factors = df_rois[order(-df_rois$varperc),][1:5,"FACTOR"]
region_order = c('FRONTAL','PARIETAL','SENSORY','TEMPORAL','MEDIAL_TEMPORAL','OCCIPITAL',
                 'MEDIAL_OCCIPITAL','CINGULATE','BASALGANGLIA',
                 'THALAMUS','CEREBRAL_WHITE','CEREBELLUM_GRAY','CEREBELLUM_WHITE')
braak_order = c('BRAAK1','BRAAK2','BRAAK3','BRAAK4','BRAAK5','BRAAK6')
asym_order = c('FRONTAL_ASYM','PARIETAL_ASYM','SENSORY_ASYM','TEMPORAL_ASYM',
               'MEDIAL_TEMPORAL_ASYM','OCCIPITAL_ASYM','MEDIAL_OCCIPITAL_ASYM',
               'CINGULATE_ASYM','BASALGANGLIA_ASYM','THALAMUS_ASYM',
               'CEREBRAL_WHITE_ASYM','CEREBELLUM_GRAY_ASYM','CEREBELLUM_WHITE_ASYM')
factor_levels = df_rois[order(-df_rois$varperc),'FACTOR']
df_rois$FACTOR = factor(df_rois$FACTOR, levels = factor_levels)

# plot sum of squared loadings
p = ggplot(df_rois,aes(x=FACTOR,y=varperc)) +
  geom_point() + 
  theme(axis.ticks=element_line(), 
        plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14, angle=30), 
        axis.text.x=element_text(angle=300, face='bold', size=14, hjust=0, vjust=1, color='grey50')) +
  ggtitle('Sum of Squared Region Loadings, by Factor') +
  xlab('Factors') +
  ylab('Sum of squared region loadings / # regions')
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

# print asym

df_rois.asym = melt(df_rois[df_rois$FACTOR %in% best_asym,],id.vars='FACTOR',measure.vars=asym_order)
p = ggplot(df_rois.asym,aes_string(x='variable',y='value', color='FACTOR')) +
  geom_line(aes_string(group='FACTOR')) +
  geom_point() + 
  theme(axis.ticks=element_line(), 
        plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14, angle=30), 
        axis.text.x=element_text(angle=300, face='bold', size=14, hjust=0, vjust=1, color='grey50')) + 
  geom_hline(aes(yintercept=0.0)) + 
  ggtitle('Regional Asymmetry, by Factor') +
  xlab('Brain Regions') +
  ylab('Asymmetry Index')
print(p)
