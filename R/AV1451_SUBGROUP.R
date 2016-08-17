library(mixtools)
library(mclust)
library(ggplot2)
library(ROCR)

df_scores = read.csv('nsfa/av1451skull_pattern_dataset.csv')

valid_diags = c('N','SMC','EMCI','LMCI','AD')
df_scores = df_scores[which(df_scores$Diag.AV1451 %in% valid_diags),]
df_scores$Diag.AV1451 = as.factor(df_scores$Diag.AV1451, levels=valid_diags)
df_scores$AV1451_BL_closest_AV45_wcereb_BIN1.11 = as.factor(df_scores$AV1451_BL_closest_AV45_wcereb_BIN1.11)


target = 'ADAS_AV1451_1'
target.slope = 'ADAS_retroslope_AV1451_BL'

# remove target na
# df_scores = df_scores[complete.cases(df_scores[,c(target)]),]
# df_scores.slope = df_scores[complete.cases(df_scores[,c(target.slope)]),]

# remove target outliers
# target.mean = mean(df_scores[,target])
# target.sd = sd(df_scores[,target])
# target.slope.mean = mean(df_scores.slope[,target.slope])
# target.slope.sd = sd(df_scores.slope[,target.slope])
# df_scores = df_scores[df_scores[,target] <= target.mean+target.sd*5,]
# df_scores = df_scores[df_scores[,target] >= target.mean-target.sd*5,]
# df_scores.slope = df_scores.slope[df_scores.slope[,target.slope] <= target.slope.mean+target.slope.sd*5,]
# df_scores.slope = df_scores.slope[df_scores.slope[,target.slope] >= target.slope.mean-target.slope.sd*5,]

threshold_by = 'AV1451_BL_closest_AV45_wcereb_BIN1.11'
suvr.neg = df_scores[df_scores[,threshold_by] == 0,]
suvr.pos = df_scores[df_scores[,threshold_by] == 1,]
suvr.slope.neg = df_scores.slope[df_scores.slope[,threshold_by] == 0,]
suvr.slope.pos = df_scores.slope[df_scores.slope[,threshold_by] == 1,]

# diag violin plots
p1 = ggplot(df_scores, aes_string('Diag.AV1451','NSFA_0')) +
  geom_violin(trim=TRUE,aes(fill=Diag.AV1451),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab('NSFA_0 Score') +
  xlab('Diagnosis')

p2 = ggplot(df_scores, aes_string('AV1451_BL_closest_AV45_wcereb_BIN1.11','NSFA_0')) +
  geom_violin(trim=TRUE,
              na.rm=TRUE,
              aes(fill=AV1451_BL_closest_AV45_wcereb_BIN1.11),show.legend=FALSE) +
  coord_flip() + 
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        legend.title=element_text(size=14)) + 
  ylab('NSFA_0 Score') +
  xlab('AV45 Positivity')

print(p1)
print(p2)



# Negatives
fm.neg.suvr = summary(lm(target_suvr_form, suvr.neg))
fm.neg.suvr.r2 = fm.neg.suvr$adj.r.squared
fm.neg.suvr.pvalue = pf(fm.neg.suvr$fstatistic['value'],
                        fm.neg.suvr$fstatistic['numdf'],
                        fm.neg.suvr$fstatistic['dendf'],
                        lower.tail=F)
fm.neg.suvr.r2.note = paste("Adjusted~italic(R)^2 ==",round(fm.neg.suvr.r2,3))
fm.neg.suvr.pvalue.note = paste("P ==",round(fm.neg.suvr.pvalue,3))

fm.neg.nsfa = summary(lm(target_nsfa6_form, suvr.neg))
fm.neg.nsfa.r2 = fm.neg.nsfa$adj.r.squared
fm.neg.nsfa.pvalue = pf(fm.neg.nsfa$fstatistic['value'],
                        fm.neg.nsfa$fstatistic['numdf'],
                        fm.neg.nsfa$fstatistic['dendf'],
                        lower.tail=F)
fm.neg.nsfa.r2.note = paste("Adjusted~italic(R)^2 ==",round(fm.neg.nsfa.r2,3))
fm.neg.nsfa.pvalue.note = paste("P ==",round(fm.neg.nsfa.pvalue,3))

fm.neg.slope.suvr = summary(lm(targetslope_suvr_form, suvr.slope.neg))
fm.neg.slope.suvr.r2 = fm.neg.slope.suvr$adj.r.squared
fm.neg.slope.suvr.pvalue = pf(fm.neg.slope.suvr$fstatistic['value'],
                              fm.neg.slope.suvr$fstatistic['numdf'],
                              fm.neg.slope.suvr$fstatistic['dendf'],
                              lower.tail=F)
fm.neg.slope.suvr.r2.note = paste("Adjusted~italic(R)^2 ==",round(fm.neg.slope.suvr.r2,3))
fm.neg.slope.suvr.pvalue.note = paste("P ==",round(fm.neg.slope.suvr.pvalue,3))

fm.neg.slope.nsfa = summary(lm(targetslope_nsfa6_form, suvr.slope.neg))
fm.neg.slope.nsfa.r2 = fm.neg.slope.nsfa$adj.r.squared
fm.neg.slope.nsfa.pvalue = pf(fm.neg.slope.nsfa$fstatistic['value'],
                              fm.neg.slope.nsfa$fstatistic['numdf'],
                              fm.neg.slope.nsfa$fstatistic['dendf'],
                              lower.tail=F)
fm.neg.slope.nsfa.r2.note = paste("Adjusted~italic(R)^2 ==",round(fm.neg.slope.nsfa.r2,3))
fm.neg.slope.nsfa.pvalue.note = paste("P ==",round(fm.neg.slope.nsfa.pvalue,3))


p1 = ggplot(suvr.neg, aes_string(x='NSFA_6', y=target)) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') + 
  ggtitle(paste('NSFA_6 versus',target,'(AV45 SUVR-)')) +
  xlab('NSFA_6 Score') + 
  annotate("text",
           x=max(suvr.neg$NSFA_6),y=max(suvr.neg[,target]),
           hjust=1,size=6,
           label=fm.neg.nsfa.r2.note,parse=TRUE) + 
  annotate("text",
           x=max(suvr.neg$NSFA_6),y=max(suvr.neg[,target]),
           hjust=1,vjust=2.2,size=6,
           label=fm.neg.nsfa.pvalue.note,parse=TRUE)
p2 = ggplot(suvr.slope.neg, aes_string(x='NSFA_6', y=target.slope)) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') + 
  ggtitle(paste('NSFA_6 versus',target.slope,'(AV45 SUVR-)')) +
  xlab('NSFA_6 Score') + 
  annotate("text",
           x=max(suvr.slope.neg$NSFA_6),y=max(suvr.slope.neg[,target.slope]),
           hjust=1,size=6,
           label=fm.neg.slope.nsfa.r2.note,parse=TRUE) + 
  annotate("text",
           x=max(suvr.slope.neg$NSFA_6),y=max(suvr.slope.neg[,target.slope]),
           hjust=1,vjust=2.2,size=6,
           label=fm.neg.slope.nsfa.pvalue.note,parse=TRUE)
p3 = ggplot(suvr.neg, aes_string(x='CORTICAL_SUMMARY_prior', y=target)) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') + 
  ggtitle(paste('AV45 SUVR versus',target,'(AV45 SUVR-)')) +
  xlab('AV45 SUVR') + 
  annotate("text",
           x=max(suvr.neg$CORTICAL_SUMMARY_prior),y=max(suvr.neg[,target]),
           hjust=1,size=6,
           label=fm.neg.suvr.r2.note,parse=TRUE) + 
  annotate("text",
           x=max(suvr.neg$CORTICAL_SUMMARY_prior),y=max(suvr.neg[,target]),
           hjust=1,vjust=2.2,size=6,
           label=fm.neg.suvr.pvalue.note,parse=TRUE)
p4 = ggplot(suvr.slope.neg, aes_string(x='CORTICAL_SUMMARY_prior', y=target.slope)) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') + 
  ggtitle(paste('AV45 SUVR versus',target.slope,'(AV45 SUVR-)')) +
  xlab('AV45 SUVR') + 
  annotate("text",
           x=max(suvr.slope.neg$CORTICAL_SUMMARY_prior),y=max(suvr.slope.neg[,target.slope]),
           hjust=1,size=6,
           label=fm.neg.slope.suvr.r2.note,parse=TRUE) + 
  annotate("text",
           x=max(suvr.slope.neg$CORTICAL_SUMMARY_prior),y=max(suvr.slope.neg[,target.slope]),
           hjust=1,vjust=2.2,size=6,
           label=fm.neg.slope.suvr.pvalue.note,parse=TRUE)
print(p1)
print(p2)
print(p3)
print(p4)

# Positives
fm.pos.suvr = summary(lm(target_suvr_form, suvr.pos))
fm.pos.suvr.r2 = fm.pos.suvr$adj.r.squared
fm.pos.suvr.pvalue = pf(fm.pos.suvr$fstatistic['value'],
                        fm.pos.suvr$fstatistic['numdf'],
                        fm.pos.suvr$fstatistic['dendf'],
                        lower.tail=F)
fm.pos.suvr.r2.note = paste("Adjusted~italic(R)^2 ==",round(fm.pos.suvr.r2,3))
fm.pos.suvr.pvalue.note = paste("P ==",round(fm.pos.suvr.pvalue,3))

fm.pos.nsfa = summary(lm(target_nsfa6_form, suvr.pos))
fm.pos.nsfa.r2 = fm.pos.nsfa$adj.r.squared
fm.pos.nsfa.pvalue = pf(fm.pos.nsfa$fstatistic['value'],
                        fm.pos.nsfa$fstatistic['numdf'],
                        fm.pos.nsfa$fstatistic['dendf'],
                        lower.tail=F)
fm.pos.nsfa.r2.note = paste("Adjusted~italic(R)^2 ==",round(fm.pos.nsfa.r2,3))
fm.pos.nsfa.pvalue.note = paste("P ==",round(fm.pos.nsfa.pvalue,3))

fm.pos.slope.suvr = summary(lm(targetslope_suvr_form, suvr.slope.pos))
fm.pos.slope.suvr.r2 = fm.pos.slope.suvr$adj.r.squared
fm.pos.slope.suvr.pvalue = pf(fm.pos.slope.suvr$fstatistic['value'],
                              fm.pos.slope.suvr$fstatistic['numdf'],
                              fm.pos.slope.suvr$fstatistic['dendf'],
                              lower.tail=F)
fm.pos.slope.suvr.r2.note = paste("Adjusted~italic(R)^2 ==",round(fm.pos.slope.suvr.r2,3))
fm.pos.slope.suvr.pvalue.note = paste("P ==",round(fm.pos.slope.suvr.pvalue,3))

fm.pos.slope.nsfa = summary(lm(targetslope_nsfa6_form, suvr.slope.pos))
fm.pos.slope.nsfa.r2 = fm.pos.slope.nsfa$adj.r.squared
fm.pos.slope.nsfa.pvalue = pf(fm.pos.slope.nsfa$fstatistic['value'],
                              fm.pos.slope.nsfa$fstatistic['numdf'],
                              fm.pos.slope.nsfa$fstatistic['dendf'],
                              lower.tail=F)
fm.pos.slope.nsfa.r2.note = paste("Adjusted~italic(R)^2 ==",round(fm.pos.slope.nsfa.r2,3))
fm.pos.slope.nsfa.pvalue.note = paste("P ==",round(fm.pos.slope.nsfa.pvalue,3))


p1 = ggplot(suvr.pos, aes_string(x='NSFA_6', y=target)) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') + 
  ggtitle(paste('NSFA_6 versus',target,'(AV45 SUVR+)')) +
  xlab('NSFA_6 Score') + 
  annotate("text",
           x=max(suvr.pos$NSFA_6),y=max(suvr.pos[,target]),
           hjust=1,size=6,
           label=fm.pos.nsfa.r2.note,parse=TRUE) + 
  annotate("text",
           x=max(suvr.pos$NSFA_6),y=max(suvr.pos[,target]),
           hjust=1,vjust=2.2,size=6,
           label=fm.pos.nsfa.pvalue.note,parse=TRUE)
p2 = ggplot(suvr.slope.pos, aes_string(x='NSFA_6', y=target.slope)) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') + 
  ggtitle(paste('NSFA_6 versus',target.slope,'(AV45 SUVR+)')) +
  xlab('NSFA_6 Score') + 
  annotate("text",
           x=max(suvr.slope.pos$NSFA_6),y=max(suvr.slope.pos[,target.slope]),
           hjust=1,size=6,
           label=fm.pos.slope.nsfa.r2.note,parse=TRUE) + 
  annotate("text",
           x=max(suvr.slope.pos$NSFA_6),y=max(suvr.slope.pos[,target.slope]),
           hjust=1,vjust=2.2,size=6,
           label=fm.pos.slope.nsfa.pvalue.note,parse=TRUE)
p3 = ggplot(suvr.pos, aes_string(x='CORTICAL_SUMMARY_prior', y=target)) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') + 
  ggtitle(paste('AV45 SUVR versus',target,'(AV45 SUVR+)')) +
  xlab('AV45 SUVR') + 
  annotate("text",
           x=max(suvr.pos$CORTICAL_SUMMARY_prior),y=max(suvr.pos[,target]),
           hjust=1,size=6,
           label=fm.pos.suvr.r2.note,parse=TRUE) + 
  annotate("text",
           x=max(suvr.pos$CORTICAL_SUMMARY_prior),y=max(suvr.pos[,target]),
           hjust=1,vjust=2.2,size=6,
           label=fm.pos.suvr.pvalue.note,parse=TRUE)
p4 = ggplot(suvr.slope.pos, aes_string(x='CORTICAL_SUMMARY_prior', y=target.slope)) +
  geom_point(aes_string(color='Diag.AV45')) +
  geom_smooth(method='lm') + 
  ggtitle(paste('AV45 SUVR versus',target.slope,'(AV45 SUVR+)')) +
  xlab('AV45 SUVR') + 
  annotate("text",
           x=max(suvr.slope.pos$CORTICAL_SUMMARY_prior),y=max(suvr.slope.pos[,target.slope]),
           hjust=1,size=6,
           label=fm.pos.slope.suvr.r2.note,parse=TRUE) + 
  annotate("text",
           x=max(suvr.slope.pos$CORTICAL_SUMMARY_prior),y=max(suvr.slope.pos[,target.slope]),
           hjust=1,vjust=2.2,size=6,
           label=fm.pos.slope.suvr.pvalue.note,parse=TRUE)
print(p1)
print(p2)
print(p3)
print(p4)





# NSFA_6_CLASS performance
pred = prediction(as.integer(df_scores$NSFA_6_CLASS), as.integer(df_scores$AD))
# Recall-Precision curve             
RP.perf = performance(pred, "prec", "rec");
plot(RP.perf);
# ROC curve
ROC.perf = performance(pred, "tpr", "fpr");
plot(ROC.perf);
# ROC area under the curve
auc.tmp = performance(pred,"auc");
auc = as.numeric(auc.tmp@y.values)

# SUVR performance
pred = prediction(as.integer(df_scores$CORTICAL_SUMMARY_prior_BIN1.11), as.integer(df_scores$AD))
# Recall-Precision curve             
RP.perf = performance(pred, "prec", "rec");
plot(RP.perf);
# ROC curve
ROC.perf = performance(pred, "tpr", "fpr");
plot(ROC.perf);
# ROC area under the curve
auc.tmp = performance(pred,"auc");
auc = as.numeric(auc.tmp@y.values)

# Diag bar plots
nsfa.neg = df_scores[df_scores$NSFA_6_CLASS == 0,]
nsfa.pos = df_scores[df_scores$NSFA_6_CLASS == 1,]
suvr.neg = df_scores[df_scores$CORTICAL_SUMMARY_prior_BIN1.11 == 0,]
suvr.pos = df_scores[df_scores$CORTICAL_SUMMARY_prior_BIN1.11 == 1,]
hist(suvr.neg$CORTICAL_SUMMARY_prior)
hist(nsfa.neg$NSFA_6)
p1 = ggplot(nsfa.neg, aes_string('Diag.AV45')) +
  geom_bar() + 
  ggtitle('NSFA NEG')
p2 = ggplot(nsfa.pos, aes_string('Diag.AV45')) +
  geom_bar() + 
  ggtitle('NSFA POS')
p3 = ggplot(suvr.neg, aes_string('Diag.AV45')) +
  geom_bar() + 
  ggtitle('SUVR NEG')
p4 = ggplot(suvr.pos, aes_string('Diag.AV45')) +
  geom_bar() + 
  ggtitle('SUVR POS')
print(p1)
print(p2)
print(p3)
print(p4)

