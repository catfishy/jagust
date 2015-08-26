from utils import *

lut_file = "../FreeSurferColorLUT.txt"

frontal=[1003,1012,1014,1018,1019,1020,1027,1028,1032,2003,2012,2014,2018,2019,2020,2027,2028,2032]
parietal=[1008,1025,1029,1031,2008,2025,2029,2031]
temporal=[1015,1030,2015,2030]
cingulate=[1002,1010,1023,1026,2002,2010,2023,2026]
basal_ganglia=[11,50,12,51,13,52]
occipital=[1011,2011]
left_cerebGM=[8]
right_cerebGM=[47]
left_cerebWM=[7]
right_cerebWM=[46]
brainstem=[16]
hemiWM=[2,41,77,251,252,253,254,255,1004,2004]
non_hemiWM=[85]
other=[0,10,17,18,26,28,30,49,53,54,58,60,62,80,1000,1001,1005,1006,1007,1009,1013,1016,1017,1021,1022,1024,1033,1034,1035,2000,2001,2005,2006,2007,2009,2013,2016,2017,2021,2022,2024,2033,2034,2035]
all_groups=[frontal,parietal,temporal,cingulate,basal_ganglia,occipital,left_cerebGM,right_cerebGM,left_cerebWM,right_cerebWM,brainstem,hemiWM,non_hemiWM,other]
names = ['frontal','parietal','temporal','cingulate','basal_ganglia','occipital','left_cerebGM','right_cerebGM','left_cerebWM','right_cerebWM','brainstem','hemiWM','non_hemiWM','other']
print "\n\nGROUPING 1"
filepath = "grouping_1.mat"
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)

left_frontal=[1003,1012,1014,1018,1019,1020,1027,1028,1032]
right_frontal=[2003,2012,2014,2018,2019,2020,2027,2028,2032]
left_parietal=[1008,1025,1029,1031]
right_parietal=[2008,2025,2029,2031]
left_temporal=[1015,1030]
right_temporal=[2015,2030]
left_cingulate=[1002,1010,1023,1026]
right_cingulate=[2002,2010,2023,2026]
left_basal_ganglia=[11,12,13]
right_basal_ganglia=[50,51,52]
left_occipital=[1011]
right_occipital=[2011]
left_cerebGM=[8]
right_cerebGM=[47]
left_cerebWM=[7]
right_cerebWM=[46]
brainstem=[16]
hemiWM=[2,41,77,251,252,253,254,255,1004,2004]
non_hemiWM=[85]
other=[0,10,17,18,26,28,30,49,53,54,58,60,62,80,1000,1001,1005,1006,1007,1009,1013,1016,1017,1021,1022,1024,1033,1034,1035,2000,2001,2005,2006,2007,2009,2013,2016,2017,2021,2022,2024,2033,2034,2035]
all_groups=[left_frontal,right_frontal,left_parietal,right_parietal,left_temporal,right_temporal,left_cingulate,right_cingulate,
            left_basal_ganglia,right_basal_ganglia,left_occipital,right_occipital,left_cerebGM,right_cerebGM,left_cerebWM,right_cerebWM,brainstem,hemiWM,non_hemiWM,other]
names = ['left_frontal','right_frontal','left_parietal','right_parietal','left_temporal','right_temporal','left_cingulate','right_cingulate',
         'left_basal_ganglia','right_basal_ganglia','left_occipital','right_occipital','left_cerebGM','right_cerebGM','left_cerebWM','right_cerebWM','brainstem','hemiWM','non_hemiWM','other']
print "\n\nGROUPING 2"
filepath = "grouping_2.mat"
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)


frontal=[1003,1012,1014,1018,1019,1020,1027,1028,1032,2003,2012,2014,2018,2019,2020,2027,2028,2032]
parietal=[1008,1025,1029,1031,2008,2025,2029,2031]
temporal=[1015,1030,2015,2030]
cingulate=[1002,1010,1023,1026,2002,2010,2023,2026]
basal_ganglia=[11,50,12,51,13,52]
occipital=[1011,2011,1005,2005,1013,2013,1021,2021]
central_gyri=[1017,1022,1024,2017,2024,2022]
thalamus=[10,49]
left_cerebGM=[8]
right_cerebGM=[47]
left_cerebWM=[7]
right_cerebWM=[46]
brainstem=[16]
hemiWM=[2,41,77,251,252,253,254,255,1004,2004]
non_hemiWM=[85]
hippocampus=[17,53,1006,1009,1016,1033,1034,2006,2009,2016,2033,2034]
other=[0,18,26,28,30,54,58,60,62,80,1000,1001,1007,1035,2000,2001,2007,2035]
all_groups=[frontal,parietal,temporal,cingulate,basal_ganglia,occipital,central_gyri,thalamus,left_cerebGM,right_cerebGM,left_cerebWM,right_cerebWM,brainstem,hemiWM,non_hemiWM,hippocampus,other]
names = ['frontal','parietal','temporal','cingulate','basal_ganglia','occipital','central_gyri','thalamus','left_cerebGM','right_cerebGM','left_cerebWM','right_cerebWM','brainstem','hemiWM','non_hemiWM','hippocampus','other']
print "\n\nGROUPING 3"
filepath = "grouping_3.mat"
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)

left_frontal=[1003,1012,1014,1018,1019,1020,1027,1028,1032]
right_frontal=[2003,2012,2014,2018,2019,2020,2027,2028,2032]
left_parietal=[1008,1025,1029,1031]
right_parietal=[2008,2025,2029,2031]
left_temporal=[1015,1030]
right_temporal=[2015,2030]
left_cingulate=[1002,1010,1023,1026]
right_cingulate=[2002,2010,2023,2026]
left_basal_ganglia=[11,12,13]
right_basal_ganglia=[50,51,52]
left_occipital=[1011,1005,1013,1021]
right_occipital=[2011,2005,2013,2021]
left_central_gyri=[1017,1024,1022]
right_central_gyri=[2017,2024,2022]
left_thalamus=[10]
right_thalamus=[49]
left_cerebGM=[8]
right_cerebGM=[47]
left_cerebWM=[7]
right_cerebWM=[46]
brainstem=[16]
hemiWM=[2,41,77,251,252,253,254,255,1004,2004]
non_hemiWM=[85]
left_hippocampus=[17,1006,1009,1016,1033,1034]
right_hippocampus=[53,2006,2009,2016,2033,2034]
other=[0,18,26,28,30,54,58,60,62,80,1000,1001,1007,1035,2000,2001,2007,2035]
all_groups=[left_frontal,right_frontal,left_parietal,right_parietal,left_temporal,right_temporal,left_cingulate,right_cingulate,left_basal_ganglia,right_basal_ganglia,
            left_occipital,right_occipital,left_central_gyri,right_central_gyri,left_thalamus,right_thalamus,
            left_cerebGM,right_cerebGM,left_cerebWM,right_cerebWM,brainstem,hemiWM,non_hemiWM,left_hippocampus,right_hippocampus,other]
names = ['left_frontal','right_frontal','left_parietal','right_parietal','left_temporal','right_temporal','left_cingulate','right_cingulate','left_basal_ganglia','right_basal_ganglia',
         'left_occipital','right_occipital','left_central_gyri','right_central_gyri','left_thalamus','right_thalamus',
         'left_cerebGM','right_cerebGM','left_cerebWM','right_cerebWM','brainstem','hemiWM','non_hemiWM','left_hippocampus','right_hippocampus','other']
print "\n\nGROUPING 4"
filepath = "grouping_4.mat"
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)

Cluster0=[1027, 2014, 2027, 1014]
Cluster1=[24, 85, 14]
Cluster2=[7, 46, 52, 16, 28, 60, 13]
Cluster3=[18, 2016, 2033, 2000, 1033, 1016, 54, 1000]
Cluster4=[2034, 1034, 1035, 2035]
Cluster5=[44, 5]
Cluster6=[1004, 2004]
Cluster7=[72]
Cluster8=[1026, 2026, 2002, 1002]
Cluster9=[1011, 2005, 1013, 1005, 2013, 2011]
Cluster10=[2, 41]
Cluster11=[2024, 1024, 1022, 2029, 1029, 2022]
Cluster12=[1032, 2032]
Cluster13=[15, 47, 1006, 2006, 8]
Cluster14=[1023, 2010, 2023, 1010, 1025, 2025]
Cluster15=[1021, 2021]
Cluster16=[50, 11, 62, 10, 30, 49]
Cluster17=[1009, 1007, 2009, 2019, 1030, 2015, 1015, 2030, 2007, 1019]
Cluster18=[53, 17, 80]
Cluster19=[253, 254, 252]
Cluster20=[4, 31, 63, 43]
Cluster21=[1001, 2001]
Cluster22=[26, 58]
Cluster23=[251]
Cluster24=[77, 255]
Cluster25=[2020, 1020, 2018, 1012, 2012, 1018]
Cluster26=[51, 12]
Cluster27=[1003, 2028, 1008, 2003, 1031, 2017, 2031, 2008, 1017, 1028]
all_groups=[Cluster0,Cluster1,Cluster2,Cluster3,Cluster4,Cluster5,Cluster6,Cluster7,Cluster8,Cluster9,Cluster10,Cluster11,Cluster12,Cluster13,Cluster14,Cluster15,Cluster16,Cluster17,Cluster18,Cluster19,Cluster20,Cluster21,Cluster22,Cluster23,Cluster24,Cluster25,Cluster26,Cluster27]
names=['Cluster0', 'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4', 'Cluster5', 'Cluster6', 'Cluster7', 'Cluster8', 'Cluster9', 'Cluster10', 'Cluster11', 'Cluster12', 'Cluster13', 'Cluster14', 'Cluster15', 'Cluster16', 'Cluster17', 'Cluster18', 'Cluster19', 'Cluster20', 'Cluster21', 'Cluster22', 'Cluster23', 'Cluster24', 'Cluster25', 'Cluster26', 'Cluster27']
print '\n\nGROUPING KMEANS'
filepath = "grouping_kmeans_26.mat"
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)

Cluster0=[53, 50, 11, 17, 80]
Cluster1=[24, 85, 14]
Cluster2=[1011, 2005, 1013, 1005, 2013, 2011]
Cluster3=[1027, 2020, 2014, 2027, 1020, 2018, 1012, 2012, 1018, 1014]
Cluster4=[15, 47, 1006, 2006, 8]
Cluster5=[2028, 1028, 26, 58]
Cluster6=[62, 10, 30, 49]
Cluster7=[4, 31, 63, 43]
Cluster8=[2024, 1024, 1030, 2030, 1022, 2029, 1029, 2022]
Cluster9=[1009, 2009, 2019, 2015, 1015, 1019]
Cluster10=[1003, 1008, 2003, 1031, 2017, 2031, 2008, 1017]
Cluster11=[44, 5]
Cluster12=[1007, 2034, 1034, 1035, 2007, 2035]
Cluster13=[1021, 2021]
Cluster14=[18, 2016, 2000, 1016, 54, 1000]
Cluster15=[2, 41]
Cluster16=[1001, 2001]
Cluster17=[1023, 1026, 2026, 2010, 2023, 1010, 2002, 1002, 1025, 2025]
Cluster18=[77, 255]
Cluster19=[253, 254, 252]
Cluster20=[51, 12]
Cluster21=[7, 46, 16, 28, 60]
Cluster22=[251]
Cluster23=[72]
Cluster24=[1004, 2004]
Cluster25=[1032, 2032]
Cluster26=[2033, 1033]
Cluster27=[52, 13]
all_groups=[Cluster0,Cluster1,Cluster2,Cluster3,Cluster4,Cluster5,Cluster6,Cluster7,Cluster8,Cluster9,Cluster10,Cluster11,Cluster12,Cluster13,Cluster14,Cluster15,Cluster16,Cluster17,Cluster18,Cluster19,Cluster20,Cluster21,Cluster22,Cluster23,Cluster24,Cluster25,Cluster26,Cluster27]
names=['Cluster0', 'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4', 'Cluster5', 'Cluster6', 'Cluster7', 'Cluster8', 'Cluster9', 'Cluster10', 'Cluster11', 'Cluster12', 'Cluster13', 'Cluster14', 'Cluster15', 'Cluster16', 'Cluster17', 'Cluster18', 'Cluster19', 'Cluster20', 'Cluster21', 'Cluster22', 'Cluster23', 'Cluster24', 'Cluster25', 'Cluster26', 'Cluster27']
print '\n\nGROUPING AGG'
filepath = "grouping_kmeans_55.mat"
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)
