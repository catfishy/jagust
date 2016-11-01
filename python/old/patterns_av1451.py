from sklearn.mixture import DPGMM

from patterns import *



def trainFactorModel(pattern_prior_df, pattern_post_df, tracer='AV45'):
    # Scale inputs
    patterns_only, post_patterns_only = scaleInput(pattern_prior_df, pattern_post_df)

    # Choose alpha + do model selection
    components = len(patterns_only.index)
    
    best_model = chooseBestFactorModel(patterns_only, components)

    print best_model
    '''
    alpha = best_model.alpha
    with open("%s_%s_%s_model.pkl" % (generateFileRoot(alpha, model=True),covar_type,tracer), 'wb') as fid:
        cPickle.dump(best_model, fid)   
    '''

def chooseBestFactorModel(patterns, components):
    # Cluster
    models_scores = {}

    for i in range(40):
        g = DPGMM_SCALE(n_components=components,
                        gamma_shape=gamma_shape,
                        gamma_inversescale=gamma_inversescale,
                        covariance_type=covar_type,
                        tol=1e-4,
                        n_iter=1000,
                        params='wmc',
                        init_params='wmc',
                        verbose=False)
        print "Fitting: %s" % i
        g.fit(patterns)
        print g.converged_
        print g.alpha
        y_ = g.predict(patterns)
        print Counter(y_)
        membership = np.zeros((len(patterns.index),components))
        for i,member in enumerate(y_):
            membership[i][member] = 1
        bound = g.lower_bound(patterns, membership)
        score = np.mean(g.score(patterns))
        models_scores[bound] = g

    best_score = max(models_scores.keys())
    best_model = models_scores[best_score]
    return best_model


if __name__ == '__main__':
    # SETUP
    bilateral=False
    lut_file = "../FreeSurferColorLUT.txt"
    membership_conf = 0.50
    threshold = 1.15
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
    data_csv = '../datasets/pvc_adni_av1451/mostregions_output.csv'
    master_csv = '../FDG_AV45_COGdata/FDG_AV45_COGdata_04_07_16.csv'

    data = parseRawDataset(data_csv, master_csv, pattern_keys, lobe_keys, tracer='AV1451', ref_key='WHOLECEREB', bilateral=bilateral)
    pattern_bl_df = data['pattern_bl_df']
    pattern_scan2_df = data['pattern_scan2_df']
    pattern_scan3_df = data['pattern_scan3_df']
    pattern_prior_df = data['pattern_prior_df']
    pattern_post_df = data['pattern_post_df']
    uptake_prior_df = data['uptake_prior_df']
    uptake_post_df = data['uptake_post_df']
    lobes_prior_df = data['lobes_prior_df']
    lobes_post_df = data['lobes_post_df']
    lobes_change_df = data['lobes_change_df']
    result_df = data['result_df']
    rchange_df = data['change_df']
    pattern_col_order = list(pattern_prior_df.columns)

    print pattern_bl_df

    # Calculate
    trainDPGMM(pattern_prior_df, pattern_post_df, result_df, 'spherical', bilateral, tracer='AV1451', shape=1.0, inversescale=1.0)




