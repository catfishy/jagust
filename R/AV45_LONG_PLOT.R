library(ggplot2)

data = read.csv("~/Google Drive/ADNI_shared/Andy Documentation/FDG_AV45_COG_data/03.01.16_long3timepts.csv")
data$Diag.AV45_long = factor(data$Diag.AV45_long, levels=c("N","SMC","EMCI","LMCI","AD"))
data$RID = factor(data$RID)
data = data[complete.cases(data[, 'Diag.AV45_long']),]

threetimepts = subset(data, AV45_NONTP_3_BigRef_BIN.79>-1)
twotimepts = subset(data, AV45_NONTP_2_BigRef_BIN.79>-1)

ggplot(twotimepts,aes(x=Age, y=AV45_NONTP_WMratio, color=RID)) + 
  facet_grid(Diag.AV45_long ~ AV45_NONTP_WMratio_Slope_pos) +
  theme(legend.position="none") +
  geom_line() + 
  geom_point()