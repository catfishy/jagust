library(ggplot2)
library(reshape2)

lobe_keys = c('FRONTAL','PARIETAL','CINGULATE','SENSORY','TEMPORAL','OCCIPITAL',
              'MEDIALOCCIPITAL','BASALGANGLIA','LIMBIC','CEREBGM',
              'BRAIN_STEM','CEREBWM','CEREBRAL_WHITE')

patterns = read.csv('dpgmm_alpha12.66_bilateral_pattern_means.csv')
patterns$pattern = factor(patterns$pattern)
patterns$pattern = with(patterns, paste0(pattern,' (',round(mean_suvr,2),')'))

patterns.m = melt(patterns, id.vars = c("pattern","mean_suvr"))
patterns.m$variable = as.character(patterns.m$variable)
patterns.m$variable = factor(patterns.m$variable, levels=lobe_keys)
patterns.m$pattern = factor(patterns.m$pattern, levels=(patterns.m$pattern)[order(patterns.m$mean_suvr)])

p = ggplot(patterns.m, aes(variable,pattern)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low='white', high='dodgerblue', name='Total\nUptake\nContrib.') + 
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.ticks=element_line(), 
        plot.title=element_text(size=20),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(face='bold', size=14, angle=30), 
        axis.text.x=element_text(angle=300, face='bold', size=14, hjust=0, vjust=1, color='grey50')) + 
  ggtitle('Regional Contributions to Global Florbetapir Uptake, by Pattern') +
  xlab('Brain Regions') +
  ylab('Patterns (Avg. Cortical Summary SUVR)')
  
print(p)
