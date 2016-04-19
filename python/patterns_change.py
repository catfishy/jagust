from patterns import *

if __name__ == '__main__':
    # SETUP
    bilateral=True
    membership_conf = 0.50
    threshold = 1.15
    scale_type = 'original'
    norm_type = 'L1'
    tracer = 'AV45'

    data_csv = '../datasets/pvc_adni_av45/mostregions_output.csv'
    master_csv = '../FDG_AV45_COGdata/FDG_AV45_COGdata_04_07_16.csv'
    lut_file = "../FreeSurferColorLUT.txt"

    result_keys = ['CORTICAL_SUMMARY_post', 'CORTICAL_SUMMARY_prior', 'CORTICAL_SUMMARY_change', 'diag_prior', 'diag_post']
    lobe_keys = ['FRONTAL','PARIETAL','CINGULATE','TEMPORAL','OCCIPITAL','MEDIALOCCIPITAL',
                 'SENSORY','BASALGANGLIA','LIMBIC','CEREBGM','CEREBWM','CEREBRAL_WHITE','BRAIN_STEM']
    master_keys = ['Age@AV45','Gender','APOE2_BIN','APOE4_BIN','Edu.(Yrs)','SMOKING','DIABETES',
                   'UW_MEM_BL_3months','UW_MEM_slope',
                   'UW_EF_BL_3months','UW_EF_slope',
                   'WMH_percentOfICV_AV45_6MTHS','WMH_percentOfICV_slope',
                   'CSF_TAU_closest_AV45','CSF_TAU_slope',
                   'CSF_ABETA_closest_AV45','CSF_ABETA_slope',
                   'FSX_HC/ICV_BL_3months','FSX_HC/ICV_slope',
                   'FDG_PONS_AV45_6MTHS','FDG_postAV45_slope']
    cat_keys = ['APOE4_BIN','APOE2_BIN','Gender','Handedness','SMOKING','DIABETES']
    summary_keys = ['CORTICAL_SUMMARY_post', 'CORTICAL_SUMMARY_prior', 'CORTICAL_SUMMARY_change']
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

    data = parseRawDataset(data_csv, master_csv, pattern_keys, lobe_keys, tracer='AV45', ref_key='WHOLECEREB', norm_type=norm_type, bilateral=bilateral)

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

    # PRETRAINED
    model_file = '../dpgmm_alpha15.68_bilateral_spherical_AV45_change_model.pkl'
    best_model = cPickle.load(open(model_file, 'rb'))
    alpha = best_model.alpha

    # dump LMM Dataset
    LME_dep_variables = [('ADAScog.','TIMEpostAV45_ADAS.'),   
                         ('UW_MEM_','UW_MEM_postAV45_'),
                         ('UW_EF_','UW_EF_postAV45_'),
                         ('AV45','')]
    for dep_val_var, dep_time_var in LME_dep_variables:
        if dep_val_var == 'AV45':
            saveAV45LMEDataset(best_model, master_csv, pattern_change_df, result_df, bilateral, confident=True)
        else:
            saveLMEDataset(best_model, master_csv, pattern_change_df, result_df, dep_val_var, dep_time_var, bilateral, categorize=False, confident=True)
