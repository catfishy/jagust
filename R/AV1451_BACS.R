library(plyr)


df_bacs = read.csv('MegaSpreadSheet_Tau_BACS.csv',stringsAsFactors=FALSE)
bacs_fields = c('NSFA_0','NSFA_3','NSFA_4','Age_at_Tau',
           'Braak12_SUVR_PVCc345','Braak34_SUVR_PVCc346','Braak56_SUVR_PVCc347',
           'global_SUVR_PVC','global_SUVR_MNI','MMSE_tau','group')
# bacs_groups = c(1,2,3,3.5,4,52)
bacs_groups = c(1,2,3,3.5,4,51,52,53,54,55)
for (field in bacs_fields) {
  df_bacs[,field] = as.numeric(df_bacs[,field])
}
df_bacs = df_bacs[df_bacs$LBNL.ID != '',]
df_bacs = df_bacs[df_bacs$group %in% bacs_groups,]
row.names(df_bacs) = df_bacs$LBNL.ID
df_bacs$LBNL.ID = NULL
df_bacs = df_bacs[,c(bacs_fields)]
df_bacs = rename(df_bacs, c("Braak12_SUVR_PVCc345"="AV1451_PVC_Braak12_CerebGray_BL",
                            "Braak34_SUVR_PVCc346"="AV1451_PVC_Braak34_CerebGray_BL",
                            "Braak56_SUVR_PVCc347"="AV1451_PVC_Braak56_CerebGray_BL",
                            "global_SUVR_PVC"="AV1451_PVC_BraakAll_CerebGray_BL",
                            "global_SUVR_MNI"="BRAAKALL_ROI_MASK_MEAN",
                            "Age_at_Tau"="Age.AV1451"))


df_av1451 = read.csv('nsfa/av1451skull_pattern_dataset.csv')
adni_fields = c('NSFA_0','NSFA_3','NSFA_4','Age.AV1451',
                'AV1451_PVC_Braak12_CerebGray_BL',
                'AV1451_PVC_Braak34_CerebGray_BL',
                'AV1451_PVC_Braak56_CerebGray_BL',
                'AV1451_PVC_BraakAll_CerebGray_BL',
                'Diag.AV1451')
row.names(df_av1451) = df_av1451$RID
df_av1451$RID = NULL
df_av1451$Diag.AV1451 = as.character(df_av1451$Diag.AV1451)
df_av1451 = df_av1451[,adni_fields]

df_mask = read.csv('maskdata_inorm_cerebgm_snorm.csv')
row.names(df_mask) = df_mask$RID
df_mask$RID = NULL

df_av1451 = merge(df_av1451,df_mask,by="row.names")
row.names(df_av1451) = df_av1451$Row.names
df_av1451$Row.names = NULL
df_av1451 = df_av1451[,c(adni_fields,'BRAAKALL_ROI_MASK_MEAN')]
df_av1451 = rename(df_av1451, c("Diag.AV1451"="group"))


all_fields = c('NSFA_0',
               'Age.AV1451',
               'AV1451_PVC_BraakAll_CerebGray_BL',
               'BRAAKALL_ROI_MASK_MEAN',
               'group')
df_all = rbind(df_av1451[,all_fields],df_bacs[,all_fields])

# different orderings
print("NSFA0 ordering")
df_nsfa0_order = df_all[complete.cases(df_all$NSFA_0),]
df_nsfa0_order = df_nsfa0_order[order(df_nsfa0_order$NSFA_0),]
paste("'",paste(row.names(df_nsfa0_order),collapse="','"),"'",sep='')

print("BRAAKALL NONPVC ordering")
df_nonpvc_order = df_all[complete.cases(df_all$BRAAKALL_ROI_MASK_MEAN),]
df_nonpvc_order = df_nonpvc_order[order(df_nonpvc_order$BRAAKALL_ROI_MASK_MEAN),]
paste("'",paste(row.names(df_nonpvc_order),collapse="','"),"'",sep='')

print("BRAAKALL PVC ordering")
df_pvc_order = df_all[complete.cases(df_all$AV1451_PVC_BraakAll_CerebGray_BL),]
df_pvc_order = df_pvc_order[order(df_pvc_order$AV1451_PVC_BraakAll_CerebGray_BL),]
paste("'",paste(row.names(df_pvc_order),collapse="','"),"'",sep='')

# determine reference group
young_normal = row.names(df_all[df_all$group == 1,])
needed = 20 - length(young_normal)

lowest_pvcorder = row.names(head(df_pvc_order[df_pvc_order$group != 1,],needed))
lowest_nonpvcorder = row.names(head(df_nonpvc_order[df_nonpvc_order$group != 1,],needed))
lowest_nsfa0order = row.names(head(df_nsfa0_order[df_nsfa0_order$group != 1,],needed))

print("REF NSFA0")
paste("'",paste(c(young_normal,lowest_nsfa0order),collapse="','"),"'",sep='')
print("REF NONPVC")
paste("'",paste(c(young_normal,lowest_nonpvcorder),collapse="','"),"'",sep='')
print("REF PVC")
paste("'",paste(c(young_normal,lowest_pvcorder),collapse="','"),"'",sep='')


df_regions = read.csv('datasets/pvc_adni_av1451/tauskullregions_uptake_bilateral.csv')
row.names(df_regions) = df_regions$subject
df_regions$subject = NULL
regions = colnames(df_regions)
df_av1451_regions = merge(df_regions,df_av1451,by='row.names')
region_with_age = data.frame()
for (region in regions) {
  form = paste('Age.AV1451 ~',region)
  region_with_age[region,'r2'] = summary(lm(form,df_av1451_regions))$adj.r.squared
}



