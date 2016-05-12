import os, sys

import cPickle
import pandas as pd
from sklearn.preprocessing import StandardScaler
import scipy

from utils import importRoussetCSV, bilateralTranslations, importFreesurferLookup, saveFakeAparcInput
from patterns import parseRawDataset, scaleRawInput

def savePatternAsMat(pattern_bl_df, output_file):
    scipy.io.savemat(output_file,{'ID': list(pattern_bl_df.index), 
                                            'Features': list(pattern_bl_df.columns), 
                                            'Y':pattern_bl_df.T.as_matrix()})

def savePatternAsAparc(df, lut_file, bilateral, out_template):
    if bilateral:
        index_lookup = bilateralTranslations(lut_file)
    else:
        lut_table = importFreesurferLookup(lut_file)
        index_lookup = {v.replace('-','_').upper():[k] for k,v in lut_table.iteritems()}
    for colname in df.columns:
        output_name = out_template % colname
        pattern = dict(df[colname])
        saveFakeAparcInput(output_name,pattern,index_lookup)

# SETUP FILES
#master_csv = '../FDG_AV45_COGdata/FDG_AV45_COGdata_04_07_16.csv'
#data_csv = '../datasets/pvc_adni_av45/mostregions_output.csv'
#pattern_mat = '../av45_pattern_bl.mat'
#dod = False

master_csv = '../DOD_DATA/DOD_DATA_05_06_16.csv'
data_csv = '../datasets/pvc_dod_av45/mostregions_output.csv'
pattern_mat = '../dod_av45_pattern_bl.mat'
dod = True

# PATTERN INPUT FILES
# nsfa_activation_csv = '../nsfa/av45_factor_activations.csv'
# nsfa_loading_csv = '../nsfa/av45_factor_loadings.csv'
# model_file = '../dpgmm_alpha12.89_bilateral_spherical_AV45_model_L1.pkl'

nsfa_activation_csv = '../nsfa/dod_av45_factor_activations.csv'
nsfa_loading_csv = '../nsfa/dod_av45_factor_loadings.csv'
model_file = None

# OUTPUT FILES
# output_file = '../nsfa/pattern_dataset.csv'
# nsfa_output_template = "../output/fake_aparc_inputs/nsfa/factor_loading_%s"
# igmm_output_template = "../output/fake_aparc_inputs/igmm/pattern_loading_%s"

output_file = '../nsfa/dod_pattern_dataset.csv'
nsfa_output_template = "../output/fake_aparc_inputs/nsfa/dod_factor_loading_%s"
igmm_output_template = "../output/fake_aparc_inputs/igmm/dod_pattern_loading_%s"


dod = True
bilateral = True
scale_type = 'original'
norm_type = 'L1'
tracer = 'AV45'
lut_file = "../FreeSurferColorLUT.txt"
by_subj, threshold = importRoussetCSV(data_csv, translate_threshold=1.11)

lobe_keys = ['FRONTAL','PARIETAL','CINGULATE','TEMPORAL','OCCIPITAL','MEDIALOCCIPITAL',
             'SENSORY','BASALGANGLIA','LIMBIC','CEREBGM','CEREBWM','CEREBRAL_WHITE','BRAIN_STEM']
pattern_keys = ['BRAIN_STEM', 'CTX_LH_BANKSSTS', 'CTX_LH_CAUDALANTERIORCINGULATE', 'CTX_LH_CAUDALMIDDLEFRONTAL', 'CTX_LH_CUNEUS', 'CTX_LH_ENTORHINAL', 'CTX_LH_FRONTALPOLE', 
                'CTX_LH_FUSIFORM', 'CTX_LH_INFERIORPARIETAL', 'CTX_LH_INFERIORTEMPORAL', 'CTX_LH_INSULA', 'CTX_LH_ISTHMUSCINGULATE', 'CTX_LH_LATERALOCCIPITAL', 
                'CTX_LH_LATERALORBITOFRONTAL', 'CTX_LH_LINGUAL', 'CTX_LH_MEDIALORBITOFRONTAL', 'CTX_LH_MIDDLETEMPORAL', 'CTX_LH_PARACENTRAL', 'CTX_LH_PARAHIPPOCAMPAL', 
                'CTX_LH_PARSOPERCULARIS', 'CTX_LH_PARSORBITALIS', 'CTX_LH_PARSTRIANGULARIS', 'CTX_LH_PERICALCARINE', 'CTX_LH_POSTCENTRAL', 'CTX_LH_POSTERIORCINGULATE', 
                'CTX_LH_PRECENTRAL', 'CTX_LH_PRECUNEUS', 'CTX_LH_ROSTRALANTERIORCINGULATE', 'CTX_LH_ROSTRALMIDDLEFRONTAL', 'CTX_LH_SUPERIORFRONTAL', 'CTX_LH_SUPERIORPARIETAL', 
                'CTX_LH_SUPERIORTEMPORAL', 'CTX_LH_SUPRAMARGINAL', 'CTX_LH_TEMPORALPOLE', 'CTX_LH_TRANSVERSETEMPORAL', 'CTX_RH_BANKSSTS', 'CTX_RH_CAUDALANTERIORCINGULATE', 
                'CTX_RH_CAUDALMIDDLEFRONTAL', 'CTX_RH_CUNEUS', 'CTX_RH_ENTORHINAL', 'CTX_RH_FRONTALPOLE', 'CTX_RH_FUSIFORM', 'CTX_RH_INFERIORPARIETAL', 
                'CTX_RH_INFERIORTEMPORAL', 'CTX_RH_INSULA', 'CTX_RH_ISTHMUSCINGULATE', 'CTX_RH_LATERALOCCIPITAL', 'CTX_RH_LATERALORBITOFRONTAL', 'CTX_RH_LINGUAL', 
                'CTX_RH_MEDIALORBITOFRONTAL', 'CTX_RH_MIDDLETEMPORAL', 'CTX_RH_PARACENTRAL', 'CTX_RH_PARAHIPPOCAMPAL', 'CTX_RH_PARSOPERCULARIS', 'CTX_RH_PARSORBITALIS', 
                'CTX_RH_PARSTRIANGULARIS', 'CTX_RH_PERICALCARINE', 'CTX_RH_POSTCENTRAL', 'CTX_RH_POSTERIORCINGULATE', 'CTX_RH_PRECENTRAL', 'CTX_RH_PRECUNEUS', 
                'CTX_RH_ROSTRALANTERIORCINGULATE', 'CTX_RH_ROSTRALMIDDLEFRONTAL', 'CTX_RH_SUPERIORFRONTAL', 'CTX_RH_SUPERIORPARIETAL', 'CTX_RH_SUPERIORTEMPORAL', 
                'CTX_RH_SUPRAMARGINAL', 'CTX_RH_TEMPORALPOLE', 'CTX_RH_TRANSVERSETEMPORAL', 'LEFT_ACCUMBENS_AREA', 'LEFT_AMYGDALA', 'LEFT_CAUDATE', 'LEFT_CEREBELLUM_CORTEX', 
                'LEFT_CEREBELLUM_WHITE_MATTER', 'LEFT_CEREBRAL_WHITE_MATTER', 'LEFT_HIPPOCAMPUS', 'LEFT_PALLIDUM', 'LEFT_PUTAMEN', 'LEFT_THALAMUS_PROPER', 
                'RIGHT_ACCUMBENS_AREA', 'RIGHT_AMYGDALA', 'RIGHT_CAUDATE', 'RIGHT_CEREBELLUM_CORTEX', 'RIGHT_CEREBELLUM_WHITE_MATTER', 
                'RIGHT_CEREBRAL_WHITE_MATTER', 'RIGHT_HIPPOCAMPUS', 'RIGHT_PALLIDUM', 'RIGHT_PUTAMEN', 'RIGHT_THALAMUS_PROPER']
if bilateral: # truncate keys
    pattern_keys = list(set([_.replace('LH_','').replace('RH_','').replace('RIGHT_','').replace('LEFT_','') for _ in pattern_keys]))

data = parseRawDataset(data_csv, master_csv, pattern_keys, lobe_keys, tracer='AV45', ref_key='WHOLECEREB', norm_type=norm_type, bilateral=bilateral, dod=dod)
pattern_bl_df = data['pattern_bl_df']
pattern_scan2_df = data['pattern_scan2_df']
pattern_scan3_df = data['pattern_scan3_df']
pattern_prior_df = data['pattern_prior_df']
pattern_post_df = data['pattern_post_df']
pattern_change_df = data['pattern_change_df']
uptake_prior_df = data['uptake_prior_df']
uptake_post_df = data['uptake_post_df']
lobes_prior_df = data['lobes_prior_df']
lobes_post_df = data['lobes_post_df']
lobes_change_df = data['lobes_change_df']
result_df = data['result_df']
rchange_df = data['change_df']
pattern_col_order = list(pattern_prior_df.columns)

# SAVE PATTERN MAT FILE (FOR INPUT INTO NSFA)
#savePatternAsMat(pattern_bl_df, pattern_mat)

# Create NSFA patterns df
nsfa_act_df = pd.read_csv(nsfa_activation_csv).T
nsfa_load_df = pd.read_csv(nsfa_loading_csv).T
nsfa_load_df.columns = nsfa_act_df.columns = ['NSFA_%s' % _ for _ in nsfa_act_df.columns]
nsfa_act_df.index = nsfa_act_df.index.astype('int64')

# Create IGMM patterns df
if model_file is not None:
    model = cPickle.load(open(model_file, 'rb'))
    means = model.means_
    bl_patterns_only, scaler = scaleRawInput(pattern_prior_df, scale_type='original')
    proba_df_bl = pd.DataFrame(model.predict_proba(bl_patterns_only)).set_index(bl_patterns_only.index)
    col_probs = proba_df_bl.sum(axis=0)
    valid_groups = list(col_probs[col_probs>0.001].index)
    igmm_load_df = pd.DataFrame({i:means[i] for i in valid_groups})
    igmm_load_df.index = pattern_col_order
    igmm_prob_df = proba_df_bl[valid_groups]
    igmm_prob_df.columns = ['IGMM_%s' % _ for _ in igmm_prob_df.columns]
    igmm_prob_df.index = igmm_prob_df.index.droplevel(1)
else:
    igmm_prob_df = pd.DataFrame()

# Get master data
if dod:
    master_df = pd.read_csv(master_csv, low_memory=False)
    master_df.set_index('SCRNO', inplace=True)
    columns = ['Age','Sex','Edu','APOE4BIN','Diag_closest_AV45_BL',
               'ANTIDEP_USE','SSRI','PTGroup','GroupNum','GroupNum_TBI','GroupNum_PTSD']
    columns += [_ for _ in master_df.columns if _.startswith('ADAS')]
    columns += [_ for _ in master_df.columns if _.startswith('AVLT')]
    columns += [_ for _ in master_df.columns if _.startswith('GD')]
    columns += [_ for _ in master_df.columns if _.startswith('CDR_GLOBAL')]
    columns += [_ for _ in master_df.columns if _.startswith('CDR_SOB')]
    other_df = master_df[columns]
else:
    master_df = pd.read_csv(master_csv, low_memory=False, header=[0,1])
    master_df.columns = master_df.columns.get_level_values(1)
    master_df.set_index('RID', inplace=True)
    columns = ['Age@AV45','Gender','APOE2_BIN','APOE4_BIN','Edu.(Yrs)','Diag@AV45_long','UCB_FS_HC/ICV_slope']
    columns += [_ for _ in master_df.columns if _.startswith('ADAScog.') or _.startswith('TIMEpostAV45_ADAS.')]
    columns += ['ADAS_3MTH_AV45','ADASslope_postAV45']
    columns += [_ for _ in master_df.columns if _.startswith('AVLT.') or _.startswith('TIMEpostAV45_AVLT.')]
    columns += ['AVLT_AV45_1_3MTHS','AVLTslope_postAV45']
    columns += [_ for _ in master_df.columns if _.startswith('WMH_percentOfICV.') or _.startswith('WMH_postAV45.')]
    columns += ['WMH_percentOfICV_AV45_6MTHS','WMH_percentOfICV_slope']
    columns += [_ for _ in master_df.columns if _.startswith('UW_MEM_') or _.startswith('UW_MEM_postAV45_')]
    columns += [_ for _ in master_df.columns if _.startswith('UW_EF_') or _.startswith('UW_EF_postAV45_')]
    other_df = master_df[columns]

print columns

# Apply threshold
result_df['positive_prior'] = (result_df['CORTICAL_SUMMARY_prior'] >= threshold).astype(int)
result_df['positive_post'] = (result_df['CORTICAL_SUMMARY_post'] >= threshold).astype(int)
result_df['ad_prior'] = (result_df['diag_prior'] == 'AD').astype(int)
result_df['ad_post'] = (result_df['diag_post'] == 'AD').astype(int)

# Combine result_df, igmm_prob_df, nsfa_act_df, and other_df
if len(igmm_prob_df.index) > 0:
    combined_df = igmm_prob_df.merge(nsfa_act_df,left_index=True,right_index=True)
else:
    combined_df = nsfa_act_df
combined_df = combined_df.merge(result_df,left_index=True,right_index=True)
combined_df = combined_df.merge(other_df,left_index=True,right_index=True)
if dod:
    combined_df.index.name = 'SCRNO'
else:
    combined_df.index.name = 'RID'

combined_df.to_csv(output_file,index=True)

# # Save loading patterns
savePatternAsAparc(nsfa_load_df, lut_file, bilateral, nsfa_output_template)
# savePatternAsAparc(igmm_load_df, lut_file, bilateral, igmm_output_template)
