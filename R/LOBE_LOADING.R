library(ggplot2)
library(ggrepel)
library(reshape2)

# df_rois = read.csv('nsfa/av1451uni_roi_comparisons.csv')
# df_scores = read.csv('nsfa/av1451uni_pattern_dataset.csv')
# df_loadings = read.csv('nsfa/av1451uni_factor_loadings.csv')

df_rois = read.csv('nsfa/av1451_roi_comparisons.csv')
df_scores = read.csv('nsfa/av1451_pattern_dataset.csv')
df_loadings = read.csv('nsfa/av1451_factor_loadings.csv')

# df_rois = read.csv('nsfa/av45_roi_comparisons.csv')
# df_scores = read.csv('nsfa/av45_pattern_dataset.csv')

pattern_columns = Filter(function(i) {startsWith(i,'NSFA_')}, names(df_scores))

region_order = c('FRONTAL','PARIETAL','SENSORY','TEMPORAL','MEDIAL_TEMPORAL','OCCIPITAL',
                 'MEDIAL_OCCIPITAL','CINGULATE','BASALGANGLIA',
                 'THALAMUS','CEREBRAL_WHITE','CEREBELLUM_GRAY','CEREBELLUM_WHITE')
braak_order = c('BRAAK1','BRAAK2','BRAAK3','BRAAK4','BRAAK5','BRAAK6')
braak_cum_order = c('BRAAK1_CUM','BRAAK2_CUM','BRAAK3_CUM','BRAAK4_CUM','BRAAK5_CUM','BRAAK6_CUM')
neg_asym_order = c('FRONTAL_NEG_ASYM','PARIETAL_NEG_ASYM','SENSORY_NEG_ASYM','TEMPORAL_NEG_ASYM',
               'MEDIAL_TEMPORAL_NEG_ASYM','OCCIPITAL_NEG_ASYM','MEDIAL_OCCIPITAL_NEG_ASYM',
               'CINGULATE_NEG_ASYM','BASALGANGLIA_NEG_ASYM','THALAMUS_NEG_ASYM',
               'CEREBRAL_WHITE_NEG_ASYM','CEREBELLUM_GRAY_NEG_ASYM','CEREBELLUM_WHITE_NEG_ASYM')
pos_asym_order = c('FRONTAL_POS_ASYM','PARIETAL_POS_ASYM','SENSORY_POS_ASYM','TEMPORAL_POS_ASYM',
               'MEDIAL_TEMPORAL_POS_ASYM','OCCIPITAL_POS_ASYM','MEDIAL_OCCIPITAL_POS_ASYM',
               'CINGULATE_POS_ASYM','BASALGANGLIA_POS_ASYM','THALAMUS_POS_ASYM',
               'CEREBRAL_WHITE_POS_ASYM','CEREBELLUM_GRAY_POS_ASYM','CEREBELLUM_WHITE_POS_ASYM')

df_rois$ASYM_ABS = rowSums(abs(df_rois[, c(neg_asym_order,pos_asym_order)]))
patterns_by_score = sort(colSums(abs(df_scores[,pattern_columns])), decreasing=TRUE)

best_asym = df_rois[order(-df_rois$ASYM_ABS),][1:4,"FACTOR"]
best_factors = df_rois[order(-df_rois$varperc),][1:4,"FACTOR"]
best_factor_scores = names(patterns_by_score)[1:4]

factor_levels = df_rois[order(-df_rois$varperc),'FACTOR']
df_rois$FACTOR = factor(df_rois$FACTOR, levels = factor_levels)

# plot sum of squared loadings
p = ggplot(df_rois,aes(x=FACTOR,y=varperc)) +
  geom_point() + 
  theme(axis.ticks=element_line(), 
        plot.title=element_text(size=20),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12, angle=30), 
        axis.text.x=element_text(angle=270, size=12, hjust=0, vjust=1, color='grey50')) +
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
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12, angle=30), 
        axis.text.x=element_text(angle=300, size=12, hjust=0, vjust=1, color='grey50')) + 
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
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12, angle=30), 
        axis.text.x=element_text(angle=300, size=12, hjust=0, vjust=1, color='grey50')) + 
  geom_hline(aes(yintercept=0.0)) + 
  ggtitle('Braak Region Loadings, by Factor') +
  xlab('Non-cumulative Braak regions') +
  ylab('Sum of loading coefficients')
print(p)

# PRINT BRAAK CUM LOADING LINES
df_rois.braakcum = melt(df_rois[df_rois$FACTOR %in% best_factors,],id.vars='FACTOR',measure.vars=braak_cum_order)
p = ggplot(df_rois.braakcum,aes_string(x='variable',y='value', color='FACTOR')) +
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
  xlab('Cumulative Braak regions') +
  ylab('Sum of loading coefficients')
print(p)

# # print asym
# for (i in best_factors) {
#   print(i)
#   cur_neg = melt(df_rois[df_rois$FACTOR == i, neg_asym_order])
#   cur_neg$variable = lapply(cur_neg$variable,function(i) {gsub('_NEG_ASYM','',i)})
#   names(cur_neg) = c('region','value_neg')
#   cur_pos = melt(df_rois[df_rois$FACTOR == i, pos_asym_order])
#   cur_pos$variable = lapply(cur_pos$variable,function(i) {gsub('_POS_ASYM','',i)})
#   names(cur_pos) = c('region','value_pos')
#   cur_points = cur_neg
#   cur_points$value_pos = cur_pos$value_pos
#   cur_points$region = factor(unlist(cur_points$region),levels=region_order)
#   
#   p = ggplot(cur_points,aes_string(x='value_neg',y='value_pos',label='region')) +
#         geom_hline(aes(yintercept=0.0)) +
#         geom_vline(aes(xintercept=0.0)) + 
#         geom_point(size=5,color='black') +
#         geom_label_repel(aes_string(fill='region'),
#                          box.padding = unit(0.25,'lines'),
#                          point.padding = unit(0.5, 'lines'),
#                          fontface='bold',
#                          color='white',
#                          size=8,
#                          show.legend=FALSE,
#                          max.iter = 7e4) + 
#         xlim(-1,1) + 
#         ylim(-1,1) +
#         theme(plot.title=element_text(size=20),
#               axis.title.x=element_text(size=18),
#               axis.title.y=element_text(size=18),
#               axis.text.x=element_text(face='bold', size=14),
#               axis.text.y=element_text(face='bold', size=14)) +
#         ggtitle(paste(i,'Loading Asymmetry')) +
#         xlab('Negative loading asymmetry index') +
#         ylab('Positive loading asymmetry index')
#   print(p)
# }

