import os, sys

import cPickle
import pandas as pd
from sklearn.preprocessing import StandardScaler
import scipy

from utils import *
from patterns import parseRawDataset, scaleRawInput

def savePatternAsMat(pattern_df, output_file):
    scipy.io.savemat(output_file,{'ID': list(pattern_df.index),
                                  'Features': list(pattern_df.columns),
                                  'Y':pattern_df.T.as_matrix()})

def calculateAsymmetry(loading_df):
    lut_file = '../FreeSurferColorLUT.txt'
    lut = importFreesurferLookup(lut_file)
    pattern_regions = list(loading_df.index)
    # create lobe encodings
    lobe_encodings = {}
    for lobe_name, lobe_idx in LOBES.iteritems():
        regions = list(set([lut[_].upper().replace('-','_') for _ in lobe_idx]))
        left_regions = [_ for _ in regions if 'LH' in _ or 'LEFT' in _]
        right_regions = [_ for _ in regions if 'RH' in _ or 'RIGHT' in _]
        left_encoding = np.array([1 if _ in left_regions else 0 for _ in pattern_regions])
        right_encoding = np.array([1 if _ in right_regions else 0 for _ in pattern_regions])
        lobe_encodings[lobe_name] = {'left': left_encoding, 'right': right_encoding}
    # create factor loading encodings
    factor_encodings = {}
    for factor in loading_df.columns:
        encoding = np.array([loading_df.loc[_,factor] for _ in pattern_regions])
        factor_encodings[factor] = encoding
    # compare
    results = []
    for factor, f_encoding in factor_encodings.iteritems():
        asyms = {'FACTOR': factor}
        # normalize
        for lobe_name, encodings in lobe_encodings.iteritems():
            left_encoding = encodings['left']
            right_encoding = encodings['right']
            pos_factor = np.array([max(_,0.0) for _ in f_encoding])
            neg_factor = np.array([abs(min(_,0.0)) for _ in f_encoding])
            left_pos_load = sum(pos_factor*left_encoding)
            right_pos_load = sum(pos_factor*right_encoding)
            left_neg_load = sum(neg_factor*left_encoding)
            right_neg_load = sum(neg_factor*right_encoding)
            if (left_pos_load + right_pos_load) == 0:
                pos_asym = 0.0
            else:
                pos_asym = (left_pos_load - right_pos_load) / (left_pos_load + right_pos_load)
            if (left_neg_load + right_neg_load) == 0:
                neg_asym = 0.0
            else:
                neg_asym = (left_neg_load - right_neg_load) / (left_neg_load + right_neg_load)
            asyms['%s_POS_ASYM' % lobe_name] = pos_asym
            asyms['%s_NEG_ASYM' % lobe_name] = neg_asym
        results.append(asyms)
    df = pd.DataFrame(results).set_index('FACTOR')
    return df


def calculateFactorScores(loading_df, pattern_df, scaler=None):
    # standardize pattern
    if scaler is None:
        scaler = StandardScaler().fit(pattern_df)
    standard_df = pd.DataFrame(scaler.transform(pattern_df))
    standard_df.set_index(pattern_df.index, inplace=True)
    standard_df.columns = pattern_df.columns
    # calculate scores
    scores_df = pd.DataFrame()
    for factor in loading_df.columns:
        coeff = loading_df[factor]
        scores = (standard_df*coeff).sum(axis=1)
        scores_df[factor] = scores
    scores_df.columns = ['SCORE_%s' % _ for _ in scores_df.columns]
    return scores_df, scaler

def compareToLobes(loading_df, bilateral):
    lut_file = '../FreeSurferColorLUT.txt'
    lut = importFreesurferLookup(lut_file)
    pattern_regions = list(loading_df.index)
    # create lobe encodings
    lobe_encodings = {}
    print LOBES
    for lobe_name, lobe_idx in LOBES.iteritems():
        if bilateral:
            regions = set([lut[_].upper().replace('-','_').replace('LH_','').replace('RH_','').replace('RIGHT_','').replace('LEFT_','') for _ in lobe_idx])
        else:
            regions = set([lut[_].upper().replace('-','_') for _ in lobe_idx])
        encoding = np.array([1 if _ in regions else 0 for _ in pattern_regions])
        lobe_encodings[lobe_name] = encoding
    # create factor loading encodings
    factor_encodings = {}
    for factor in loading_df.columns:
        encoding = np.array([loading_df.loc[_,factor] for _ in pattern_regions])
        factor_encodings[factor] = encoding
    # compare
    results = []
    for factor, f_encoding in factor_encodings.iteritems():
        # split
        pos_factor = np.array([max(_,0.0) for _ in f_encoding])
        neg_factor = np.array([abs(min(_,0.0)) for _ in f_encoding])
        pos_sims = {'FACTOR': '%s_POS' % factor}
        neg_sims = {'FACTOR': '%s_NEG' % factor}
        sims = {'FACTOR': factor}
        # normalize
        for lobe_name, encoding in lobe_encodings.iteritems():
            # cross entropy
            pos_cross = sum(pos_factor*encoding)
            neg_cross = sum(neg_factor*encoding)
            cross = pos_cross - neg_cross
            pos_sims[lobe_name] = pos_cross
            neg_sims[lobe_name] = neg_cross
            sims[lobe_name] = cross
        #results += [pos_sims, neg_sims]
        results.append(sims)
    df = pd.DataFrame(results).set_index('FACTOR')
    return df


def compareToBraakStages(loading_df, bilateral):
    lut_file = '../FreeSurferColorLUT.txt'
    lut = importFreesurferLookup(lut_file)
    pattern_regions = list(loading_df.index)
    # create braak encodings
    braak_stages = [BRAAK1,BRAAK2,BRAAK3,BRAAK4,BRAAK5,BRAAK6]
    braak_encodings = {}
    for i, stage in enumerate(braak_stages):
        if bilateral:
            regions = set([lut[_].upper().replace('-','_').replace('LH_','').replace('RH_','').replace('RIGHT_','').replace('LEFT_','') for _ in stage])
        else:
            regions = set([lut[_].upper().replace('-','_') for _ in stage])
        stage_name = 'BRAAK%s' % (i+1)
        encoding = np.array([1 if _ in regions else 0 for _ in pattern_regions])
        braak_encodings[stage_name] = encoding
    # create factor loading encodings
    factor_encodings = {}
    for factor in loading_df.columns:
        encoding = np.array([loading_df.loc[_,factor] for _ in pattern_regions])
        factor_encodings[factor] = encoding
    # compare
    results = []
    for factor, f_encoding in factor_encodings.iteritems():
        # split
        pos_factor = np.array([max(_,0.0) for _ in f_encoding])
        neg_factor = np.array([abs(min(_,0.0)) for _ in f_encoding])
        pos_sims = {'FACTOR': '%s_POS' % factor}
        neg_sims = {'FACTOR': '%s_NEG' % factor}
        sims = {'FACTOR': factor}
        # normalize
        for braak, b_encoding in braak_encodings.iteritems():
            # cross entropy
            pos_cross = sum(pos_factor*b_encoding)
            neg_cross = sum(neg_factor*b_encoding)
            cross = pos_cross - neg_cross
            pos_sims[braak] = pos_cross
            neg_sims[braak] = neg_cross
            sims[braak] = cross
        #results += [pos_sims, neg_sims]
        results.append(sims)
    df = pd.DataFrame(results).set_index('FACTOR')
    return df

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

# FOR ADNI AV45
# master_csv = '../FDG_AV45_COGdata/FDG_AV45_COGdata_06_14_16.csv'
# data_csv = '../datasets/pvc_adni_av45/mostregions_output.csv'
# pattern_mat = '../av45_pattern_bl.mat'
# pattern_mat_2 = '../av45_pattern_scan2.mat'
# pattern_mat_3 = '../av45_pattern_scan3.mat'
# nsfa_activation_csv = '../nsfa/av45_factor_activations.csv'
# nsfa_activation_csv_2 = '../nsfa/av45_factor_activations_scan2.csv'
# nsfa_activation_csv_3 = '../nsfa/av45_factor_activations_scan3.csv'
# nsfa_loading_csv = '../nsfa/av45_factor_loadings.csv'
# nsfa_lambdag_csv = '../nsfa/av45_lambdag.csv'
# # model_file = '../dpgmm_alpha12.89_bilateral_spherical_AV45_model_L1.pkl'
# model_file = None
# output_file = '../nsfa/av45_pattern_dataset.csv'
# topregions_output_file = '../nsfa/av45_top_regions.csv'
# comp_output_file = '../nsfa/av45_roi_comparisons.csv'
# comm_output_file = '../nsfa/av45_communality.csv'
# nsfa_output_template = "../output/fake_aparc_inputs/nsfa/av45_factor_loading_%s"
# igmm_output_template = "../output/fake_aparc_inputs/igmm/av45_pattern_loading_%s"
# dod = False
# bilateral = True

# FOR ADNI AV1451
# master_csv = '../FDG_AV45_COGdata/FDG_AV45_COGdata_06_14_16.csv'
# data_csv = '../datasets/pvc_adni_av1451/mostregions_output.csv'
# pattern_mat = '../av1451_pattern_bl.mat'
# pattern_mat_2 = None
# pattern_mat_3 = None
# nsfa_activation_csv = '../nsfa/av1451_factor_activations.csv'
# nsfa_activation_csv_2 = None
# nsfa_activation_csv_3 = None
# nsfa_loading_csv = '../nsfa/av1451_factor_loadings.csv'
# nsfa_lambdag_csv = '../nsfa/av1451_lambdag.csv'
# model_file = None
# output_file = '../nsfa/av1451_pattern_dataset.csv'
# topregions_output_file = '../nsfa/av1451_top_regions.csv'
# comp_output_file = '../nsfa/av1451_roi_comparisons.csv'
# comm_output_file = '../nsfa/av1451_communality.csv'
# nsfa_output_template = "../output/fake_aparc_inputs/nsfa/av1451_factor_loading_%s"
# igmm_output_template = "../output/fake_aparc_inputs/igmm/av1451_pattern_loading_%s"
# dod = False
# bilateral = True

# FOR ADNI AV1451 UNILATERAL
master_csv = '../FDG_AV45_COGdata/FDG_AV45_COGdata_06_14_16.csv'
data_csv = '../datasets/pvc_adni_av1451/mostregions_output.csv'
pattern_mat = '../av1451uni_pattern_bl.mat'
pattern_mat_2 = None
pattern_mat_3 = None
nsfa_activation_csv = '../nsfa/av1451uni_factor_activations.csv'
nsfa_activation_csv_2 = None
nsfa_activation_csv_3 = None
nsfa_loading_csv = '../nsfa/av1451uni_factor_loadings.csv'
nsfa_lambdag_csv = '../nsfa/av1451uni_lambdag.csv'
model_file = None
output_file = '../nsfa/av1451uni_pattern_dataset.csv'
topregions_output_file = '../nsfa/av1451uni_top_regions.csv'
comp_output_file = '../nsfa/av1451uni_roi_comparisons.csv'
comm_output_file = '../nsfa/av1451uni_communality.csv'
nsfa_output_template = "../output/fake_aparc_inputs/nsfa/av1451uni_factor_loading_%s"
igmm_output_template = "../output/fake_aparc_inputs/igmm/av1451uni_pattern_loading_%s"
dod = False
bilateral = False


# FOR DOD AV45
# master_csv = '../DOD_DATA/DOD_DATA_05_06_16.csv'
# data_csv = '../datasets/pvc_dod_av45/mostregions_output.csv'
# pattern_mat = '../dod_av45_pattern_bl.mat'
# pattern_mat_2 = '../dod_av45_pattern_bl.mat'
# nsfa_activation_csv = '../nsfa/dod_av45_factor_activations.csv'
# nsfa_loading_csv = '../nsfa/dod_av45_factor_loadings.csv'
# model_file = None
# output_file = '../nsfa/dod_pattern_dataset.csv'
# topregions_output_file = '../nsfa/dod_top_regions.csv'
# comp_output_file = '../nsfa/dod_braak_comparisons.csv'
# comm_output_file = '../nsfa/dod_communality.csv'
# nsfa_output_template = "../output/fake_aparc_inputs/nsfa/dod_factor_loading_%s"
# igmm_output_template = "../output/fake_aparc_inputs/igmm/dod_pattern_loading_%s"
# dod = True
# bilateral = True

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

data = parseRawDataset(data_csv, master_csv, pattern_keys, lobe_keys, tracer=tracer, ref_key='WHOLECEREB', norm_type=norm_type, bilateral=bilateral, dod=dod)
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
column_order = pattern_bl_df.columns
pattern_bl_df = pattern_bl_df[column_order]
if pattern_bl_df.index.nlevels > 1:
    pattern_bl_df.index = pattern_bl_df.index.droplevel(1)
pattern_scan2_df = pattern_scan2_df[column_order]
if pattern_scan2_df.index.nlevels > 1:
    pattern_scan2_df.index = pattern_scan2_df.index.droplevel(1)
pattern_scan3_df = pattern_scan3_df[column_order]
if pattern_scan3_df.index.nlevels > 1:
    pattern_scan3_df.index = pattern_scan3_df.index.droplevel(1)

# savePatternAsMat(pattern_bl_df, pattern_mat)
# savePatternAsMat(pattern_scan2_df, pattern_mat_2)
# savePatternAsMat(pattern_scan3_df, pattern_mat_3)
# sys.exit(1)



# Create NSFA patterns df
nsfa_act_df = pd.read_csv(nsfa_activation_csv).T
nsfa_act_df.index = nsfa_act_df.index.astype('int64')
nsfa_load_df = pd.read_csv(nsfa_loading_csv).T
nsfa_lambdag = map(float,list(pd.read_csv(nsfa_lambdag_csv).columns))
# only keep factors with nonzero loadings
columns = ['NSFA_%s' % _ for _ in nsfa_act_df.columns]
nsfa_lambdag_df = pd.DataFrame([nsfa_lambdag],columns=columns).iloc[0]
nsfa_load_df.columns = nsfa_act_df.columns = columns
nsfa_load_df = nsfa_load_df.loc[:,nsfa_load_df.sum() != 0]
nsfa_act_df = nsfa_act_df[nsfa_load_df.columns]
bl_columns = nsfa_act_df.columns
nsfa_lambdag_df = nsfa_lambdag_df[bl_columns]
# add on subsequent timepoints if available
if nsfa_activation_csv_2 is not None:
    nsfa_act_2_df = pd.read_csv(nsfa_activation_csv_2).T
    nsfa_act_2_df.index = nsfa_act_2_df.index.astype('int64')
    nsfa_act_2_df.columns = ['NSFA_%s' % _ for _ in nsfa_act_2_df.columns]
    nsfa_act_2_df = nsfa_act_2_df[bl_columns]
    nsfa_act_2_df.columns = ['SCAN2_%s' % _ for _ in nsfa_act_2_df.columns]
    nsfa_act_df = nsfa_act_df.merge(nsfa_act_2_df, left_index=True, right_index=True, how='outer')
if nsfa_activation_csv_3 is not None:
    nsfa_act_3_df = pd.read_csv(nsfa_activation_csv_3).T
    nsfa_act_3_df.index = nsfa_act_3_df.index.astype('int64')
    columns = ['NSFA_%s' % _ for _ in nsfa_act_3_df.columns]
    nsfa_act_3_df.columns = columns
    nsfa_act_3_df = nsfa_act_3_df[bl_columns]
    nsfa_act_3_df.columns = ['SCAN3_%s' % _ for _ in nsfa_act_3_df.columns]
    nsfa_act_df = nsfa_act_df.merge(nsfa_act_3_df, left_index=True, right_index=True, how='outer')


# get naive ratios
uptake_prior_df.index = uptake_prior_df.index.droplevel(1)
scores_df, scaler = calculateFactorScores(nsfa_load_df, pattern_bl_df, scaler=None)
if not pattern_scan2_df.empty:
    scores_scan2_df, scaler = calculateFactorScores(nsfa_load_df, pattern_scan2_df, scaler=scaler)
    scores_scan2_df.columns = [_.replace('SCORE','SCORE_SCAN2') for _ in scores_scan2_df.columns]
    scores_df = scores_df.merge(scores_scan2_df, left_index=True, right_index=True, how='outer')
if not pattern_scan3_df.empty:
    scores_scan3_df, scaler = calculateFactorScores(nsfa_load_df, pattern_scan3_df, scaler=scaler)
    scores_scan3_df.columns = [_.replace('SCORE','SCORE_SCAN3') for _ in scores_scan3_df.columns]
    scores_df = scores_df.merge(scores_scan3_df, left_index=True, right_index=True, how='outer')

# Look at variance explained + communality
num_vars = nsfa_load_df.shape[1]
#nsfa_standard_load_df = nsfa_load_df/np.sqrt(nsfa_lambdag_df)
nsfa_standard_load_df = nsfa_load_df
ss = np.square(nsfa_standard_load_df).sum(axis=0)/num_vars
comm = np.square(nsfa_standard_load_df).sum(axis=1)
comm.sort_values(ascending=False,inplace=True)
comm.to_csv(comm_output_file)

# Compare to braak STAGES
braak_comp_df = compareToBraakStages(nsfa_load_df, bilateral)
lobe_comp_df = compareToLobes(nsfa_load_df, bilateral)
asym_df = calculateAsymmetry(nsfa_load_df)
comp_df = braak_comp_df.merge(lobe_comp_df,left_index=True,right_index=True)
comp_df = comp_df.merge(asym_df,left_index=True,right_index=True)
comp_df['varperc'] = ss
comp_df.sort_values('varperc', ascending=False, inplace=True)
comp_df.to_csv(comp_output_file)

# Create top ten region lists for each factor
toppos_df = pd.DataFrame()
topneg_df = pd.DataFrame()
for col in nsfa_load_df.columns:
    factor_row = nsfa_load_df[col]
    top_neg = factor_row.sort_values().head(10)
    top_neg = pd.DataFrame(top_neg[top_neg<0])
    top_neg.index.name = '%s_REGIONS' % col
    top_neg.reset_index(inplace=True)
    topneg_df = topneg_df.merge(top_neg, left_index=True, right_index=True, how='outer')
    top_pos = factor_row.sort_values(ascending=False).head(10)
    top_pos = pd.DataFrame(top_pos[top_pos>0])
    top_pos.index.name = '%s_REGIONS' % col
    top_pos.reset_index(inplace=True)
    toppos_df = toppos_df.merge(top_pos, left_index=True, right_index=True, how='outer')
toppos_df.index = ['POS_%s' % (i+1) for i in toppos_df.index]
topneg_df.index = ['NEG_%s' % (i+1) for i in topneg_df.index]
topregions_df = pd.concat((toppos_df,topneg_df))
topregions_df.to_csv(topregions_output_file)


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
    columns = ['Age@AV45','Age@AV1451','Gender','APOE2_BIN','APOE4_BIN','Edu.(Yrs)',
               'Diag@AV45','Diag@AV1451','AV45_1_2_Diff','AV45_1_3_Diff',
               'AV45_NONTP_wcereb_BIN1.11',
               'AV45_NONTP_2_wcereb_BIN1.11',
               'AV45_NONTP_3_wcereb_BIN1.11']
    columns += ['UCB_FS_HC/ICV_AV45_1','UCB_FS_HC/ICV_AV1451_1','UCB_FS_HC/ICV_slope']
    columns += ['ADAS_AV45_1','ADAS_AV1451_1','ADASslope_postAV45']
    columns += ['AVLT_AV45_1','AVLT_AV1451_1','AVLT_slope_postAV45']
    columns += ['MMSEslope_postAV45','MMSE_AV45_1','MMSE_AV1451_1']
    columns += ['WMH_percentOfICV_AV45_1','WMH_percentOfICV_slope']
    columns += ['UW_MEM_AV45_1','UW_MEM_AV1451_1','UW_MEM_slope']
    columns += ['UW_EF_AV45_1','UW_EF_AV1451_1','UW_EF_slope']
    columns += ['CSF_ABETA_closest_AV45_1','CSF_ABETA_closest_AV45_1_BIN_192',
                'CSF_ABETA_closest_AV1451_1','CSF_ABETA_closest_AV1451_1_BIN_192',
                'CSF_ABETA_slope']
    columns += ['CSF_TAU_closest_AV45_1','CSF_TAU_closest_AV45_1_BIN_93',
                'CSF_TAU_closest_AV1451_1','CSF_TAU_closest_AV1451_1_BIN_93',
                'CSF_TAU_slope']
    columns += ['CSF_PTAU_closest_AV45_1','CSF_PTAU_closest_AV45_1_BIN_23',
                'CSF_PTAU_closest_AV1451_1','CSF_PTAU_closest_AV1451_1_BIN_23',
                'CSF_PTAU_slope']
    columns += ['GD_AV45_1','GD_slope']
    columns += ['CDR_GLOBAL_AV45_1',
               'CDR_MEMORY_AV45_1',
               'CDR_ORIENT_AV45_1',
               'CDR_JUDGE_AV45_1',
               'CDR_COMMUN_AV45_1',
               'CDR_HOME_AV45_1',
               'CDR_CARE_AV45_1',
               'CDR_GLOBAL_AV1451_1',
               'CDR_MEMORY_AV1451_1',
               'CDR_ORIENT_AV1451_1',
               'CDR_JUDGE_AV1451_1',
               'CDR_COMMUN_AV1451_1',
               'CDR_HOME_AV1451_1',
               'CDR_CARE_AV1451_1']
    columns += ['AV1451_Braak12_CerebGray_BL',
                'AV1451_Braak34_CerebGray_BL',
                'AV1451_Braak56_CerebGray_BL',
                'AV1451_Braak1_CerebGray_BL',
                'AV1451_Braak2_CerebGray_BL',
                'AV1451_Braak3_CerebGray_BL',
                'AV1451_Braak4_CerebGray_BL',
                'AV1451_Braak5_CerebGray_BL',
                'AV1451_Braak6_CerebGray_BL']
    other_df = master_df[columns]

print columns

# Apply threshold
result_df['positive_prior'] = (result_df['CORTICAL_SUMMARY_prior'] >= threshold).astype(int)
result_df['positive_post'] = (result_df['CORTICAL_SUMMARY_post'] >= threshold).astype(int)
result_df['ad_prior'] = (result_df['diag_prior'] == 'AD').astype(int)
result_df['ad_post'] = (result_df['diag_post'] == 'AD').astype(int)

# Combine result_df, igmm_prob_df, nsfa_act_df, scores_df, and other_df
if len(igmm_prob_df.index) > 0:
    combined_df = igmm_prob_df.merge(nsfa_act_df,left_index=True,right_index=True)
else:
    combined_df = nsfa_act_df
combined_df = combined_df.merge(result_df,left_index=True,right_index=True)
combined_df = combined_df.merge(other_df,left_index=True,right_index=True)
combined_df = combined_df.merge(scores_df,left_index=True,right_index=True)
if dod:
    combined_df.index.name = 'SCRNO'
else:
    combined_df.index.name = 'RID'

combined_df.to_csv(output_file,index=True)

# Save loading patterns
# savePatternAsAparc(nsfa_load_df, lut_file, bilateral, nsfa_output_template)
# savePatternAsAparc(igmm_load_df, lut_file, bilateral, igmm_output_template)
