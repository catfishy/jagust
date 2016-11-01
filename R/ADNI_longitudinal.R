library(ggplot2)

Incr = read.csv("~/Google Drive/ADNI_shared/PROJECT_increasing_amyloid/10.25.16_long.csv")
nnincr = Incr[ which(Incr$Diag.AV45=='N' & Incr$AV45_NONTP_1_BigRef_BIN.79==0), ]
nnincr = nnincr[complete.cases(nnincr$AV45_NONTP_EarlyAccumBigRef_Slope),]

#Create spaghetti plots of AVLT change for increasers vs decreasers
#Include points showing locations of each AVLT visit

ggplot(nnincr,aes(x=AVLT_time, y=AVLT)) +  
  geom_point() + 
  geom_line(aes(color=factor(RID)),size=0.5) + 
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  facet_grid(AV45_NONTP_EarlyAccumBigRef_Slope_INCR ~ .) + 
  theme_bw() + 
  labs(x="Time relative to baseline florbeatpir", y="AVLT") + 
  theme(axis.title = element_text(face="bold", size=18), 
        axis.text = element_text(face="bold", size=16), 
        strip.text = element_text(face="bold", size=16), 
        panel.border = element_rect(colour = "black", size=1.5),
        legend.position="none")

