from utils import *

lut_file = "../FreeSurferColorLUT.txt"
groups_included = {}

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
other=[10,17,18,26,28,30,49,53,54,58,60,62,80,1000,1001,1005,1006,1007,1009,1013,1016,1017,1021,1022,1024,1033,1034,1035,2000,2001,2005,2006,2007,2009,2013,2016,2017,2021,2022,2024,2033,2034,2035]
all_groups=[frontal,parietal,temporal,cingulate,basal_ganglia,occipital,left_cerebGM,right_cerebGM,left_cerebWM,right_cerebWM,brainstem,hemiWM,non_hemiWM,other]
names = ['frontal','parietal','temporal','cingulate','basal_ganglia','occipital','left_cerebGM','right_cerebGM','left_cerebWM','right_cerebWM','brainstem','hemiWM','non_hemiWM','other']
print "\n\nGROUPING 1"
filepath = "grouping_1.mat"
groups_included[filepath] = [i for g in all_groups for i in g]
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)
# create input for fake aparc so we can visualize grouping
out_file = '../output/fake_aparc_inputs/GROUP1.mat'
sio.savemat(out_file, {'regions': [{'inds': v, 'val': i} for i,v in enumerate(all_groups)]})


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
other=[10,17,18,26,28,30,49,53,54,58,60,62,80,1000,1001,1005,1006,1007,1009,1013,1016,1017,1021,1022,1024,1033,1034,1035,2000,2001,2005,2006,2007,2009,2013,2016,2017,2021,2022,2024,2033,2034,2035]
all_groups=[left_frontal,right_frontal,left_parietal,right_parietal,left_temporal,right_temporal,left_cingulate,right_cingulate,
            left_basal_ganglia,right_basal_ganglia,left_occipital,right_occipital,left_cerebGM,right_cerebGM,left_cerebWM,right_cerebWM,brainstem,hemiWM,non_hemiWM,other]
names = ['left_frontal','right_frontal','left_parietal','right_parietal','left_temporal','right_temporal','left_cingulate','right_cingulate',
         'left_basal_ganglia','right_basal_ganglia','left_occipital','right_occipital','left_cerebGM','right_cerebGM','left_cerebWM','right_cerebWM','brainstem','hemiWM','non_hemiWM','other']
print "\n\nGROUPING 2"
filepath = "grouping_2.mat"
groups_included[filepath] = [i for g in all_groups for i in g]
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)
# create input for fake aparc so we can visualize grouping
out_file = '../output/fake_aparc_inputs/GROUP2.mat'
sio.savemat(out_file, {'regions': [{'inds': v, 'val': i} for i,v in enumerate(all_groups)]})

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
other=[18,26,28,30,54,58,60,62,80,1000,1001,1007,1035,2000,2001,2007,2035]
all_groups=[frontal,parietal,temporal,cingulate,basal_ganglia,occipital,central_gyri,thalamus,left_cerebGM,right_cerebGM,left_cerebWM,right_cerebWM,brainstem,hemiWM,non_hemiWM,hippocampus,other]
names = ['frontal','parietal','temporal','cingulate','basal_ganglia','occipital','central_gyri','thalamus','left_cerebGM','right_cerebGM','left_cerebWM','right_cerebWM','brainstem','hemiWM','non_hemiWM','hippocampus','other']
print "\n\nGROUPING 3"
filepath = "grouping_3.mat"
groups_included[filepath] = [i for g in all_groups for i in g]
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)
# create input for fake aparc so we can visualize grouping
out_file = '../output/fake_aparc_inputs/GROUP3.mat'
sio.savemat(out_file, {'regions': [{'inds': v, 'val': i} for i,v in enumerate(all_groups)]})

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
other=[18,26,28,30,54,58,60,62,80,1000,1001,1007,1035,2000,2001,2007,2035]
all_groups=[left_frontal,right_frontal,left_parietal,right_parietal,left_temporal,right_temporal,left_cingulate,right_cingulate,left_basal_ganglia,right_basal_ganglia,
            left_occipital,right_occipital,left_central_gyri,right_central_gyri,left_thalamus,right_thalamus,
            left_cerebGM,right_cerebGM,left_cerebWM,right_cerebWM,brainstem,hemiWM,non_hemiWM,left_hippocampus,right_hippocampus,other]
names = ['left_frontal','right_frontal','left_parietal','right_parietal','left_temporal','right_temporal','left_cingulate','right_cingulate','left_basal_ganglia','right_basal_ganglia',
         'left_occipital','right_occipital','left_central_gyri','right_central_gyri','left_thalamus','right_thalamus',
         'left_cerebGM','right_cerebGM','left_cerebWM','right_cerebWM','brainstem','hemiWM','non_hemiWM','left_hippocampus','right_hippocampus','other']
print "\n\nGROUPING 4"
filepath = "grouping_4.mat"
groups_included[filepath] = [i for g in all_groups for i in g]
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)
# create input for fake aparc so we can visualize grouping
out_file = '../output/fake_aparc_inputs/GROUP4.mat'
sio.savemat(out_file, {'regions': [{'inds': v, 'val': i} for i,v in enumerate(all_groups)]})


Cluster0=[1021, 2021]
Cluster1=[2024, 1024, 1022, 2022]
Cluster2=[2034, 1034, 1035, 2035]
Cluster3=[253, 254]
Cluster4=[1003, 2028, 2003, 1028]
Cluster5=[2010, 1010, 1025, 2025]
Cluster6=[1032, 2032]
Cluster8=[1001, 2001]
Cluster9=[2020, 1020, 2018, 1018]
Cluster10=[18, 2000, 54, 1000]
Cluster11=[1009, 2009, 2015, 1015]
Cluster12=[2033, 1033]
Cluster13=[1004, 2004]
Cluster14=[26, 58]
Cluster15=[1011, 2011]
Cluster16=[255]
Cluster17=[2019, 1019]
Cluster18=[1008, 1031, 2031, 2008]
Cluster19=[1030, 2030]
Cluster20=[2002, 1002]
Cluster21=[1007, 2007]
Cluster22=[50, 11]
Cluster23=[2005, 1005]
Cluster24=[1023, 2023]
Cluster25=[2014, 1014]
Cluster26=[1006, 2006]
Cluster28=[51, 12]
Cluster29=[2017, 1017]
Cluster31=[80]
Cluster32=[28, 60]
Cluster33=[16]
Cluster34=[85]
Cluster35=[53, 17]
Cluster36=[52, 13]
Cluster37=[2, 41]
Cluster38=[10, 49]
Cluster39=[47, 8]
Cluster40=[7, 46]
Cluster41=[30]
Cluster42=[77]
Cluster44=[2029, 1029]
Cluster45=[1012, 2012]
Cluster46=[1027, 2027]
Cluster47=[252]
Cluster48=[1013, 2013]
Cluster49=[251]
Cluster50=[62]
Cluster53=[2016, 1016]
Cluster54=[1026, 2026]
all_groups=[Cluster0,Cluster1,Cluster2,Cluster3,Cluster4,Cluster5,Cluster6,Cluster8,
            Cluster9,Cluster10,Cluster11,Cluster12,Cluster13,Cluster14,Cluster15,Cluster16,
            Cluster17,Cluster18,Cluster19,Cluster20,Cluster21,Cluster22,Cluster23,Cluster24,
            Cluster25,Cluster26,Cluster28,Cluster29,Cluster31,Cluster32,Cluster33,
            Cluster34,Cluster35,Cluster36,Cluster37,Cluster38,Cluster39,Cluster40,
            Cluster41,Cluster42,Cluster44,Cluster45,Cluster46,Cluster47,Cluster48,
            Cluster49,Cluster50,Cluster53,Cluster54]
names=['Cluster0', 'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4', 'Cluster5', 'Cluster6', 
       'Cluster8', 'Cluster9', 'Cluster10', 'Cluster11', 'Cluster12', 'Cluster13', 
       'Cluster14', 'Cluster15', 'Cluster16', 'Cluster17', 'Cluster18', 'Cluster19', 'Cluster20', 
       'Cluster21', 'Cluster22', 'Cluster23', 'Cluster24', 'Cluster25', 'Cluster26',
       'Cluster28', 'Cluster29', 'Cluster31', 'Cluster32', 'Cluster33', 'Cluster34', 
       'Cluster35', 'Cluster36', 'Cluster37', 'Cluster38', 'Cluster39', 'Cluster40', 'Cluster41', 
       'Cluster42', 'Cluster44', 'Cluster45', 'Cluster46', 'Cluster47', 'Cluster48', 
       'Cluster49', 'Cluster50', 'Cluster53', 'Cluster54']
print '\n\nGROUPING AGG high'
filepath = "grouping_agghigh.mat"
groups_included[filepath] = [i for g in all_groups for i in g]
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)
# create input for fake aparc so we can visualize grouping
out_file = '../output/fake_aparc_inputs/GROUP_AGGHIGH.mat'
sio.savemat(out_file, {'regions': [{'inds': v, 'val': i} for i,v in enumerate(all_groups)]})

Cluster0=[1023, 1026, 2014, 2026, 2010, 2023, 1012, 1010, 2002, 1002, 2012, 1025, 2025, 1014]
Cluster1=[2, 41, 251]
Cluster2=[1004, 77, 2004, 255]
Cluster3=[18, 2016, 2033, 2000, 1033, 1016, 54, 1000]
Cluster4=[1030, 2030]
Cluster5=[51, 62, 12, 30]
Cluster6=[1027, 1003, 2020, 2027, 1020, 2028, 1008, 2018, 2003, 1031, 2031, 2008, 1028, 1018]
Cluster7=[53, 50, 11, 17, 80, 10, 49]
Cluster8=[1011, 2005, 1013, 1005, 2013, 2011]
Cluster9=[2019, 2015, 1015, 1019]
Cluster10=[2024, 1024, 1022, 2022]
Cluster12=[85]
Cluster14=[1032, 2032]
Cluster16=[47, 8]
Cluster17=[1001, 2001]
Cluster18=[52, 28, 60, 13]
Cluster19=[7, 46]
Cluster20=[1021, 2021]
Cluster22=[253, 254, 252]
Cluster23=[1006, 2006]
Cluster24=[1029, 2029]
Cluster25=[1007, 2007]
Cluster26=[1017, 2017]
Cluster27=[26, 58, 1035, 2035]
Cluster28=[1009, 2009]
Cluster29=[1034, 2034]
Cluster30=[16]
all_groups=[Cluster0,Cluster1,Cluster2,Cluster3,Cluster4,Cluster5,Cluster6,Cluster7,
            Cluster8,Cluster9,Cluster10,Cluster12,Cluster14,
            Cluster16,Cluster17,Cluster18,Cluster19,Cluster20,Cluster22,Cluster23,
            Cluster24,Cluster25,Cluster26,Cluster27,Cluster28,Cluster29,Cluster30]
names=['Cluster0', 'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4', 'Cluster5', 'Cluster6', 'Cluster7', 
       'Cluster8', 'Cluster9', 'Cluster10', 'Cluster12', 'Cluster14', 
       'Cluster16', 'Cluster17', 'Cluster18', 'Cluster19', 'Cluster20', 'Cluster22', 'Cluster23', 
       'Cluster24', 'Cluster25', 'Cluster26', 'Cluster27', 'Cluster28', 'Cluster29', 'Cluster30']
print "\n\nGROUPING AGG low two"
filepath = "grouping_agglowtwo.mat"
groups_included[filepath] = [i for g in all_groups for i in g]
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)
# create input for fake aparc so we can visualize grouping
out_file = '../output/fake_aparc_inputs/GROUP_AGGLOWTWO.mat'
sio.savemat(out_file, {'regions': [{'inds': v, 'val': i} for i,v in enumerate(all_groups)]})


clusters = [('ctx-lh-caudalmiddlefrontal', [1003]),
            ('ctx-lh-lateralorbitofrontal', [1012]),
            ('ctx-lh-medialorbitofrontal', [1014]),
            ('ctx-lh-parsopercularis', [1018]),
            ('ctx-lh-parsorbitalis', [1019]),
            ('ctx-lh-parstriangularis', [1020]),
            ('ctx-lh-rostralmiddlefrontal', [1027]),
            ('ctx-lh-superiorfrontal', [1028]),
            ('ctx-lh-frontalpole', [1032]), 
            ('ctx-rh-caudalmiddlefrontal', [2003]), 
            ('ctx-rh-lateralorbitofrontal', [2012]), 
            ('ctx-rh-medialorbitofrontal', [2014]), 
            ('ctx-rh-parsopercularis', [2018]), 
            ('ctx-rh-parsorbitalis', [2019]), 
            ('ctx-rh-parstriangularis', [2020]), 
            ('ctx-rh-rostralmiddlefrontal', [2027]), 
            ('ctx-rh-superiorfrontal', [2028]), 
            ('ctx-rh-frontalpole', [2032]), 
            ('ctx-lh-caudalanteriorcingulate', [1002]), 
            ('ctx-lh-isthmuscingulate', [1010]), 
            ('ctx-lh-posteriorcingulate', [1023]), 
            ('ctx-lh-rostralanteriorcingulate', [1026]), 
            ('ctx-rh-caudalanteriorcingulate', [2002]), 
            ('ctx-rh-isthmuscingulate', [2010]), 
            ('ctx-rh-posteriorcingulate', [2023]), 
            ('ctx-rh-rostralanteriorcingulate', [2026]), 
            ('ctx-lh-inferiorparietal', [1008]), 
            ('ctx-lh-precuneus', [1025]), 
            ('ctx-lh-superiorparietal', [1029]), 
            ('ctx-lh-supramarginal', [1031]), 
            ('ctx-rh-inferiorparietal', [2008]), 
            ('ctx-rh-precuneus', [2025]), 
            ('ctx-rh-superiorparietal', [2029]), 
            ('ctx-rh-supramarginal', [2031]), 
            ('ctx-lh-middletemporal', [1015]), 
            ('ctx-lh-superiortemporal', [1030]),
            ('ctx-rh-middletemporal', [2015]),
            ('ctx-rh-superiortemporal', [2030]),
            ('left_basal_ganglia', [11,12,13]),
            ('right_basal_ganglia', [50,51,52]),
            ('left_occipital', [1011,1005,1013,1021]),
            ('right_occipital', [2011,2005,2013,2021]),
            ('left_central_gyri', [1017,1024,1022]),
            ('right_central_gyri', [2017,2024,2022]),
            ('left_thalamus', [10]),
            ('right_thalamus', [49]),
            ('left_cerebGM', [8]),
            ('right_cerebGM', [47]),
            ('left_cerebWM', [7]),
            ('right_cerebWM', [46]),
            ('brainstem', [16]),
            ('hemiWM', [2,41,77,251,252,253,254,255,1004,2004]),
            ('non_hemiWM', [85]),
            ('left_hippocampus', [17,1006,1009,1016,1033,1034]),
            ('right_hippocampus', [53,2006,2009,2016,2033,2034]),
            ('other', [2000,1000,60,1035,2007,80,2001,18,2035,62,54,1001,58,1007,28,26,30])]
print "\n\nGROUPING SUMMARY"
filepath = "grouping_summary.mat"
all_groups = [b for a,b in clusters]
names = [a for a,b in clusters]
groups_included[filepath] = [i for g in all_groups for i in g]
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)
# create input for fake aparc so we can visualize grouping
out_file = '../output/fake_aparc_inputs/GROUP_SUMMARY.mat'
sio.savemat(out_file, {'regions': [{'inds': v, 'val': i} for i,v in enumerate(all_groups)]})
print groups_included


clusters = [('ctx-lh-caudalmiddlefrontal', [1003]),
            ('ctx-lh-lateralorbitofrontal', [1012]),
            ('ctx-lh-medialorbitofrontal', [1014]),
            ('ctx-lh-parsopercularis', [1018]),
            ('ctx-lh-parsorbitalis', [1019]),
            ('ctx-lh-parstriangularis', [1020]),
            ('ctx-lh-rostralmiddlefrontal', [1027]),
            ('ctx-lh-superiorfrontal', [1028]),
            ('ctx-lh-frontalpole', [1032]),
            ('ctx-rh-caudalmiddlefrontal', [2003]),
            ('ctx-rh-lateralorbitofrontal', [2012]),
            ('ctx-rh-medialorbitofrontal', [2014]),
            ('ctx-rh-parsopercularis', [2018]),
            ('ctx-rh-parsorbitalis', [2019]),
            ('ctx-rh-parstriangularis', [2020]),
            ('ctx-rh-rostralmiddlefrontal', [2027]),
            ('ctx-rh-superiorfrontal', [2028]),
            ('ctx-rh-frontalpole', [2032]),
            ('ctx-lh-caudalanteriorcingulate', [1002]),
            ('ctx-lh-isthmuscingulate', [1010]),
            ('ctx-lh-posteriorcingulate', [1023]),
            ('ctx-lh-rostralanteriorcingulate', [1026]),
            ('ctx-rh-caudalanteriorcingulate', [2002]),
            ('ctx-rh-isthmuscingulate', [2010]),
            ('ctx-rh-posteriorcingulate', [2023]),
            ('ctx-rh-rostralanteriorcingulate', [2026]),
            ('ctx-lh-inferiorparietal', [1008]),
            ('ctx-lh-precuneus', [1025]),
            ('ctx-lh-superiorparietal', [1029]),
            ('ctx-lh-supramarginal', [1031]),
            ('ctx-rh-inferiorparietal', [2008]),
            ('ctx-rh-precuneus', [2025]),
            ('ctx-rh-superiorparietal', [2029]),
            ('ctx-rh-supramarginal', [2031]),
            ('ctx-lh-middletemporal', [1015]),
            ('ctx-lh-superiortemporal', [1030]),
            ('ctx-rh-middletemporal', [2015]),
            ('ctx-rh-superiortemporal', [2030]),
            ('Left-Caudate', [11]),
            ('Left-Putamen', [12]),
            ('Left-Pallidum', [13]),
            ('Right-Caudate', [50]),
            ('Right-Putamen', [51]),
            ('Right-Pallidum', [52]),
            ('ctx-lh-lateraloccipital', [1011]),
            ('ctx-lh-cuneus', [1005]),
            ('ctx-lh-lingual', [1013]),
            ('ctx-lh-pericalcarine', [1021]),
            ('ctx-rh-lateraloccipital', [2011]),
            ('ctx-rh-cuneus', [2005]),
            ('ctx-rh-lingual', [2013]),
            ('ctx-rh-pericalcarine', [2021]),
            ('ctx-lh-paracentral', [1017]),
            ('ctx-lh-precentral', [1024]),
            ('ctx-lh-postcentral', [1022]),
            ('ctx-rh-paracentral', [2017]),
            ('ctx-rh-precentral', [2024]),
            ('ctx-rh-postcentral', [2022]),
            ('Left-Thalamus-Proper', [10]),
            ('Right-Thalamus-Proper', [49]),
            ('Left-Cerebellum-Cortex', [8]),
            ('Right-Cerebellum-Cortex', [47]),
            ('Left-Cerebellum-White-Matter', [7]),
            ('Right-Cerebellum-White-Matter', [46]),
            ('Brain-Stem', [16]),
            ('Left-Cerebral-White-Matter', [2]),
            ('Right-Cerebral-White-Matter', [41]),
            ('WM-hypointensities', [77]),
            ('CC_Posterior', [251]),
            ('CC_Mid_Posterior', [252]),
            ('CC_Central', [253]),
            ('CC_Mid_Anterior', [254]),
            ('CC_Anterior', [255]),
            ('ctx-lh-corpuscallosum', [1004]),
            ('ctx-rh-corpuscallosum', [2004]),
            ('Optic-Chiasm', [85]),
            ('Left-Hippocampus', [17]),
            ('ctx-lh-entorhinal', [1006]),
            ('ctx-lh-inferiortemporal', [1009]),
            ('ctx-lh-parahippocampal', [1016]),
            ('ctx-lh-temporalpole', [1033]),
            ('ctx-lh-transversetemporal', [1034]),
            ('Right-Hippocampus', [53]),
            ('ctx-rh-entorhinal', [2006]),
            ('ctx-rh-inferiortemporal', [2009]),
            ('ctx-rh-parahippocampal', [2016]),
            ('ctx-rh-temporalpole', [2033]),
            ('ctx-rh-transversetemporal', [2034]),
            ('ctx-rh-unknown', [2000]),
            ('ctx-lh-unknown', [1000]),
            ('Right-VentralDC', [60]),
            ('ctx-lh-insula', [1035]),
            ('ctx-rh-fusiform', [2007]),
            ('non-WM-hypointensities', [80]),
            ('ctx-rh-bankssts', [2001]),
            ('Left-Amygdala', [18]),
            ('ctx-rh-insula', [2035]),
            ('Right-vessel', [62]),
            ('Right-Amygdala', [54]),
            ('ctx-lh-bankssts', [1001]),
            ('Right-Accumbens-area', [58]),
            ('ctx-lh-fusiform', [1007]),
            ('Left-VentralDC', [28]),
            ('Left-Accumbens-area', [26]), 
            ('Left-vessel', [30])]
print "\n\nGROUPING ALLREGIONS"
filepath = "grouping_allregions.mat"
all_groups = [b for a,b in clusters]
names = [a for a,b in clusters]
groups_included[filepath] = [i for g in all_groups for i in g]
createROIGrouping(names, all_groups, filepath)
printROIGrouping(names, all_groups, lut_file)
# create input for fake aparc so we can visualize grouping
out_file = '../output/fake_aparc_inputs/GROUP_ALLREGIONS.mat'
sio.savemat(out_file, {'regions': [{'inds': v, 'val': i} for i,v in enumerate(all_groups)]})
print groups_included
