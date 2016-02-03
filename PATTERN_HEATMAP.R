library(ggplot2)
library(reshape2)

lobe_keys = c('FRONTAL','PARIETAL','CINGULATE','SENSORY','TEMPORAL','OCCIPITAL',
              'MEDIALOCCIPITAL','BASALGANGLIA','LIMBIC','THALAMUS','CEREBELLUM_GRAY',
              'BRAIN_STEM','CEREBELLUM_WHITE','CEREBRAL_WHITE')


patterns = read.csv('dpgmm_alpha15.37_bilateral_pattern_means.csv')
patterns$pattern = factor(patterns$pattern)
patterns$pattern = with(patterns, paste0(pattern,' (',round(mean_suvr,2),')'))

patterns.m = melt(patterns, id.vars = c("pattern","mean_suvr"))
#patterns.m$lobe_order = match(patterns.m$variable,lobe_keys)
#patterns.m = patterns.m[order(patterns.m$lobe_order),]
patterns.m$variable = as.character(patterns.m$variable)
patterns.m$variable = factor(patterns.m$variable, levels=lobe_keys)
patterns.m$pattern = factor(patterns.m$pattern, levels=(patterns.m$pattern)[order(patterns.m$mean_suvr)])

p = ggplot(patterns.m, aes(variable,pattern)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low='white', high='dodgerblue', name='Total\nUptake\nContrib.') + 
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.ticks=element_line(), 
        plot.title=element_text(size=22),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=12, angle=30), 
        axis.text.x=element_text(angle=330, face='bold', size=10, hjust=0, vjust=1, color='grey50')) + 
  ggtitle('Regional Contributions to Global Uptake by Pattern') +
  xlab('Brain Regions') +
  ylab('Patterns (Mean Cortical Summary SUVR)')
  
print(p)
