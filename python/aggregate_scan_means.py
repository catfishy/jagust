'''
To be run on the output of the florbetapir preprocessing pipeline
The following files are neceessary and need to be specified in the main function:

    output = '/Users/ahorng/Documents/UCBERKELEYAV45_06_04_15.csv' # the aggregated output file
    registry = "/Users/ahorng/Documents/registry_clean.csv" # the registry, with null characters removed so it can be read (used to get viscodes)
    meta_pet = "/Users/ahorng/Documents/PET_META_LIST_edited.csv" # the meta list of pet scans (with irrelevant scans removed, used to get scan dates)
        -> irrelevant scans are ones where Sequence != "AV45 Coreg, Avg, Std Img and Vox Siz, Uniform Resolution"
    bl_means = "/Users/ahorng/Documents/AV45_preprocess_output_06_11_15/AV45_BL_means_11-Jun-2015_1066.csv"
    v2_means = "/Users/ahorng/Documents/AV45_preprocess_output_06_11_15/AV45_V2_means_11-Jun-2015_607.csv"
    v3_means = "/Users/ahorng/Documents/AV45_preprocess_output_06_11_15/AV45_V3_means_11-Jun-2015_86.csv"
    bl_sizes = "/Users/ahorng/Documents/AV45_preprocess_output_06_11_15/AV45_BL_roisize_11-Jun-2015_1066.csv"
    v2_sizes = "/Users/ahorng/Documents/AV45_preprocess_output_06_11_15/AV45_V2_roisize_11-Jun-2015_607.csv"
    v3_sizes = "/Users/ahorng/Documents/AV45_preprocess_output_06_11_15/AV45_V3_roisize_11-Jun-2015_86.csv"

'''
import os
import csv
from itertools import chain, izip, repeat
import numpy as np
from collections import defaultdict
from datetime import datetime, timedelta
import copy
import codecs
from glob import glob
import pandas as pd

from utils import *

ALL_TP = ['BL','Scan2','Scan3','Scan4','Scan5']

ALL_REGION_KEYS = ['BRAIN_STEM','BRAIN_STEM_SIZE','3RD_VENTRICLE','3RD_VENTRICLE_SIZE','4TH_VENTRICLE','4TH_VENTRICLE_SIZE','5TH_VENTRICLE','5TH_VENTRICLE_SIZE','CC_ANTERIOR','CC_ANTERIOR_SIZE',
                  'CC_CENTRAL','CC_CENTRAL_SIZE','CC_MID_ANTERIOR','CC_MID_ANTERIOR_SIZE','CC_MID_POSTERIOR','CC_MID_POSTERIOR_SIZE','CC_POSTERIOR','CC_POSTERIOR_SIZE',
                  'CSF','CSF_SIZE','CTX_LH_BANKSSTS','CTX_LH_BANKSSTS_SIZE','CTX_LH_CAUDALANTERIORCINGULATE','CTX_LH_CAUDALANTERIORCINGULATE_SIZE',
                  'CTX_LH_CAUDALMIDDLEFRONTAL','CTX_LH_CAUDALMIDDLEFRONTAL_SIZE','CTX_LH_CUNEUS','CTX_LH_CUNEUS_SIZE','CTX_LH_ENTORHINAL','CTX_LH_ENTORHINAL_SIZE',
                  'CTX_LH_FRONTALPOLE','CTX_LH_FRONTALPOLE_SIZE','CTX_LH_FUSIFORM','CTX_LH_FUSIFORM_SIZE','CTX_LH_INFERIORPARIETAL','CTX_LH_INFERIORPARIETAL_SIZE',
                  'CTX_LH_INFERIORTEMPORAL','CTX_LH_INFERIORTEMPORAL_SIZE','CTX_LH_INSULA','CTX_LH_INSULA_SIZE','CTX_LH_ISTHMUSCINGULATE','CTX_LH_ISTHMUSCINGULATE_SIZE',
                  'CTX_LH_LATERALOCCIPITAL','CTX_LH_LATERALOCCIPITAL_SIZE','CTX_LH_LATERALORBITOFRONTAL','CTX_LH_LATERALORBITOFRONTAL_SIZE','CTX_LH_LINGUAL','CTX_LH_LINGUAL_SIZE',
                  'CTX_LH_MEDIALORBITOFRONTAL','CTX_LH_MEDIALORBITOFRONTAL_SIZE','CTX_LH_MIDDLETEMPORAL','CTX_LH_MIDDLETEMPORAL_SIZE','CTX_LH_PARACENTRAL','CTX_LH_PARACENTRAL_SIZE',
                  'CTX_LH_PARAHIPPOCAMPAL','CTX_LH_PARAHIPPOCAMPAL_SIZE','CTX_LH_PARSOPERCULARIS','CTX_LH_PARSOPERCULARIS_SIZE','CTX_LH_PARSORBITALIS','CTX_LH_PARSORBITALIS_SIZE',
                  'CTX_LH_PARSTRIANGULARIS','CTX_LH_PARSTRIANGULARIS_SIZE','CTX_LH_PERICALCARINE','CTX_LH_PERICALCARINE_SIZE','CTX_LH_POSTCENTRAL','CTX_LH_POSTCENTRAL_SIZE',
                  'CTX_LH_POSTERIORCINGULATE','CTX_LH_POSTERIORCINGULATE_SIZE','CTX_LH_PRECENTRAL','CTX_LH_PRECENTRAL_SIZE','CTX_LH_PRECUNEUS','CTX_LH_PRECUNEUS_SIZE',
                  'CTX_LH_ROSTRALANTERIORCINGULATE','CTX_LH_ROSTRALANTERIORCINGULATE_SIZE','CTX_LH_ROSTRALMIDDLEFRONTAL','CTX_LH_ROSTRALMIDDLEFRONTAL_SIZE',
                  'CTX_LH_SUPERIORFRONTAL','CTX_LH_SUPERIORFRONTAL_SIZE','CTX_LH_SUPERIORPARIETAL','CTX_LH_SUPERIORPARIETAL_SIZE','CTX_LH_SUPERIORTEMPORAL','CTX_LH_SUPERIORTEMPORAL_SIZE',
                  'CTX_LH_SUPRAMARGINAL','CTX_LH_SUPRAMARGINAL_SIZE','CTX_LH_TEMPORALPOLE','CTX_LH_TEMPORALPOLE_SIZE','CTX_LH_TRANSVERSETEMPORAL','CTX_LH_TRANSVERSETEMPORAL_SIZE',
                  'CTX_LH_UNKNOWN','CTX_LH_UNKNOWN_SIZE','CTX_RH_BANKSSTS','CTX_RH_BANKSSTS_SIZE','CTX_RH_CAUDALANTERIORCINGULATE','CTX_RH_CAUDALANTERIORCINGULATE_SIZE',
                  'CTX_RH_CAUDALMIDDLEFRONTAL','CTX_RH_CAUDALMIDDLEFRONTAL_SIZE','CTX_RH_CUNEUS','CTX_RH_CUNEUS_SIZE','CTX_RH_ENTORHINAL','CTX_RH_ENTORHINAL_SIZE',
                  'CTX_RH_FRONTALPOLE','CTX_RH_FRONTALPOLE_SIZE','CTX_RH_FUSIFORM','CTX_RH_FUSIFORM_SIZE','CTX_RH_INFERIORPARIETAL','CTX_RH_INFERIORPARIETAL_SIZE',
                  'CTX_RH_INFERIORTEMPORAL','CTX_RH_INFERIORTEMPORAL_SIZE','CTX_RH_INSULA','CTX_RH_INSULA_SIZE','CTX_RH_ISTHMUSCINGULATE','CTX_RH_ISTHMUSCINGULATE_SIZE',
                  'CTX_RH_LATERALOCCIPITAL','CTX_RH_LATERALOCCIPITAL_SIZE','CTX_RH_LATERALORBITOFRONTAL','CTX_RH_LATERALORBITOFRONTAL_SIZE','CTX_RH_LINGUAL','CTX_RH_LINGUAL_SIZE',
                  'CTX_RH_MEDIALORBITOFRONTAL','CTX_RH_MEDIALORBITOFRONTAL_SIZE','CTX_RH_MIDDLETEMPORAL','CTX_RH_MIDDLETEMPORAL_SIZE','CTX_RH_PARACENTRAL','CTX_RH_PARACENTRAL_SIZE',
                  'CTX_RH_PARAHIPPOCAMPAL','CTX_RH_PARAHIPPOCAMPAL_SIZE','CTX_RH_PARSOPERCULARIS','CTX_RH_PARSOPERCULARIS_SIZE','CTX_RH_PARSORBITALIS','CTX_RH_PARSORBITALIS_SIZE',
                  'CTX_RH_PARSTRIANGULARIS','CTX_RH_PARSTRIANGULARIS_SIZE','CTX_RH_PERICALCARINE','CTX_RH_PERICALCARINE_SIZE','CTX_RH_POSTCENTRAL','CTX_RH_POSTCENTRAL_SIZE',
                  'CTX_RH_POSTERIORCINGULATE','CTX_RH_POSTERIORCINGULATE_SIZE','CTX_RH_PRECENTRAL','CTX_RH_PRECENTRAL_SIZE','CTX_RH_PRECUNEUS','CTX_RH_PRECUNEUS_SIZE',
                  'CTX_RH_ROSTRALANTERIORCINGULATE','CTX_RH_ROSTRALANTERIORCINGULATE_SIZE','CTX_RH_ROSTRALMIDDLEFRONTAL','CTX_RH_ROSTRALMIDDLEFRONTAL_SIZE',
                  'CTX_RH_SUPERIORFRONTAL','CTX_RH_SUPERIORFRONTAL_SIZE','CTX_RH_SUPERIORPARIETAL','CTX_RH_SUPERIORPARIETAL_SIZE','CTX_RH_SUPERIORTEMPORAL','CTX_RH_SUPERIORTEMPORAL_SIZE',
                  'CTX_RH_SUPRAMARGINAL','CTX_RH_SUPRAMARGINAL_SIZE','CTX_RH_TEMPORALPOLE','CTX_RH_TEMPORALPOLE_SIZE','CTX_RH_TRANSVERSETEMPORAL','CTX_RH_TRANSVERSETEMPORAL_SIZE',
                  'CTX_RH_UNKNOWN','CTX_RH_UNKNOWN_SIZE','LEFT_ACCUMBENS_AREA','LEFT_ACCUMBENS_AREA_SIZE','LEFT_AMYGDALA','LEFT_AMYGDALA_SIZE','LEFT_CAUDATE','LEFT_CAUDATE_SIZE',
                  'LEFT_CEREBELLUM_CORTEX','LEFT_CEREBELLUM_CORTEX_SIZE','LEFT_CEREBELLUM_WHITE_MATTER','LEFT_CEREBELLUM_WHITE_MATTER_SIZE','LEFT_CEREBRAL_WHITE_MATTER','LEFT_CEREBRAL_WHITE_MATTER_SIZE',
                  'LEFT_CHOROID_PLEXUS','LEFT_CHOROID_PLEXUS_SIZE','LEFT_HIPPOCAMPUS','LEFT_HIPPOCAMPUS_SIZE','LEFT_INF_LAT_VENT','LEFT_INF_LAT_VENT_SIZE','LEFT_LATERAL_VENTRICLE','LEFT_LATERAL_VENTRICLE_SIZE',
                  'LEFT_PALLIDUM','LEFT_PALLIDUM_SIZE','LEFT_PUTAMEN','LEFT_PUTAMEN_SIZE','LEFT_THALAMUS_PROPER','LEFT_THALAMUS_PROPER_SIZE','LEFT_VENTRALDC','LEFT_VENTRALDC_SIZE',
                  'LEFT_VESSEL','LEFT_VESSEL_SIZE','NON_WM_HYPOINTENSITIES','NON_WM_HYPOINTENSITIES_SIZE','OPTIC_CHIASM','OPTIC_CHIASM_SIZE','RIGHT_ACCUMBENS_AREA','RIGHT_ACCUMBENS_AREA_SIZE',
                  'RIGHT_AMYGDALA','RIGHT_AMYGDALA_SIZE','RIGHT_CAUDATE','RIGHT_CAUDATE_SIZE','RIGHT_CEREBELLUM_CORTEX','RIGHT_CEREBELLUM_CORTEX_SIZE','RIGHT_CEREBELLUM_WHITE_MATTER','RIGHT_CEREBELLUM_WHITE_MATTER_SIZE',
                  'RIGHT_CEREBRAL_WHITE_MATTER','RIGHT_CEREBRAL_WHITE_MATTER_SIZE','RIGHT_CHOROID_PLEXUS','RIGHT_CHOROID_PLEXUS_SIZE','RIGHT_HIPPOCAMPUS','RIGHT_HIPPOCAMPUS_SIZE','RIGHT_INF_LAT_VENT','RIGHT_INF_LAT_VENT_SIZE',
                  'RIGHT_LATERAL_VENTRICLE','RIGHT_LATERAL_VENTRICLE_SIZE','RIGHT_PALLIDUM','RIGHT_PALLIDUM_SIZE','RIGHT_PUTAMEN','RIGHT_PUTAMEN_SIZE','RIGHT_THALAMUS_PROPER','RIGHT_THALAMUS_PROPER_SIZE',
                  'RIGHT_VENTRALDC','RIGHT_VENTRALDC_SIZE','RIGHT_VESSEL','RIGHT_VESSEL_SIZE','WM_HYPOINTENSITIES','WM_HYPOINTENSITIES_SIZE']

ADNI_FIELDNAMES = ['RID','VISCODE','VISCODE2','EXAMDATE','CEREBELLUMGREYMATTER','WHOLECEREBELLUM',
                   'ERODED_SUBCORTICALWM','COMPOSITE_REF','FRONTAL','CINGULATE','PARIETAL','TEMPORAL',
                   'SUMMARYSUVR_WHOLECEREBNORM','SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF',
                   'SUMMARYSUVR_COMPOSITE_REFNORM','SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF'] + ALL_REGION_KEYS
ADNI_FIELDNAMES_EXTRA = ['RID','VISCODE','VISCODE2','EXAMDATE','CEREBELLUMGREYMATTER','CEREBELLUMWHITEMATTER','WHOLECEREBELLUM',
                   'ERODED_SUBCORTICALWM','COMPOSITE','COMPOSITE_REF','FRONTAL','CINGULATE','PARIETAL','TEMPORAL',
                   'SUMMARYSUVR_WHOLECEREBNORM','SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF',
                   'LEFT_PUTAMEN','RIGHT_PUTAMEN','LEFT_CAUDATE','RIGHT_CAUDATE','LEFT_PALLIDUM','RIGHT_PALLIDUM',
                   'SUMMARYSUVR_COMPOSITE_REFNORM','SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF'] + ALL_REGION_KEYS
DOD_FIELDNAMES = ['SCRNO','VISCODE','EXAMDATE','CEREBELLUMGREYMATTER','WHOLECEREBELLUM',
                  'FRONTAL','CINGULATE','PARIETAL','TEMPORAL',
                  'SUMMARYSUVR_WHOLECEREBNORM','SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF',
                  'SUMMARYSUVR_COMPOSITE_REFNORM','SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF'] + ALL_REGION_KEYS
DOD_FIELDNAMES_EXTRA = ['SCRNO','VISCODE','EXAMDATE','CEREBELLUMGREYMATTER','BRAIN_STEM','WHOLECEREBELLUM',
                  'FRONTAL','FRONTAL_SIZE','CINGULATE','CINGULATE_SIZE','PARIETAL','PARIETAL_SIZE','TEMPORAL','TEMPORAL_SIZE',
                  'COMPOSITE','COMPOSITE_REF','ERODED_SUBCORTICALWM',
                  'LEFT_PUTAMEN','RIGHT_PUTAMEN','LEFT_CAUDATE','RIGHT_CAUDATE','LEFT_PALLIDUM','RIGHT_PALLIDUM',
                  'LEFT_PUTAMEN_SIZE','RIGHT_PUTAMEN_SIZE','LEFT_CAUDATE_SIZE','RIGHT_CAUDATE_SIZE','LEFT_PALLIDUM_SIZE','RIGHT_PALLIDUM_SIZE',
                  'SUMMARYSUVR_WHOLECEREBNORM','SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF',
                  'SUMMARYSUVR_COMPOSITE_REFNORM','SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF'] + ALL_REGION_KEYS
DEP_FIELDNAMES = ADNI_FIELDNAMES
DEP_FIELDNAMES_EXTRA = ADNI_FIELDNAMES_EXTRA
TAU_FIELDNAMES = ['RID','VISCODE','VISCODE2','EXAMDATE','ERODED_SUBCORTICALWM','ERODED_SUBCORTICALWM_SIZE'] + ALL_REGION_KEYS
TAU_FIELDNAMES_EXTRA = ['RID','VISCODE','VISCODE2','EXAMDATE','BRAAK1', 'BRAAK1_SIZE', 'BRAAK2', 'BRAAK2_SIZE', 'BRAAK3', 'BRAAK3_SIZE',
                        'BRAAK4', 'BRAAK4_SIZE', 'BRAAK5', 'BRAAK5_SIZE', 'BRAAK6', 'BRAAK6_SIZE','ERODED_SUBCORTICALWM','ERODED_SUBCORTICALWM_SIZE'] + ALL_REGION_KEYS
DOD_TAU_FIELDNAMES = ['SCRNO','VISCODE','EXAMDATE','ERODED_SUBCORTICALWM','ERODED_SUBCORTICALWM_SIZE'] + ALL_REGION_KEYS
DOD_TAU_FIELDNAMES_EXTRA = ['SCRNO','VISCODE','EXAMDATE','BRAAK1', 'BRAAK1_SIZE', 'BRAAK2', 'BRAAK2_SIZE', 'BRAAK3', 'BRAAK3_SIZE',
                            'BRAAK4', 'BRAAK4_SIZE', 'BRAAK5', 'BRAAK5_SIZE', 'BRAAK6', 'BRAAK6_SIZE', 'ERODED_SUBCORTICALWM','ERODED_SUBCORTICALWM_SIZE'] + ALL_REGION_KEYS


def DFWeightedMean(df, keys):
    size_keys = ['%s_SIZE' % _ for _ in keys]
    sizes = df.loc[:,size_keys].sum(axis=1)
    ratios = df.loc[:,size_keys].divide(sizes, axis=0)
    ratios.columns = [_.replace('_SIZE','') for _ in ratios.columns]
    means = df.loc[:,keys].multiply(ratios).sum(axis=1)
    return (means, sizes)

def additionalTauCalculations(df, lut_table, keys=None):
    df = df.copy()
    cerebg = [translateColumn(_, lut_table) for _ in CEREBG]
    braak1 = [translateColumn(_, lut_table) for _ in BRAAK1]
    braak2 = [translateColumn(_, lut_table) for _ in BRAAK2]
    braak3 = [translateColumn(_, lut_table) for _ in BRAAK3]
    braak4 = [translateColumn(_, lut_table) for _ in BRAAK4]
    braak5 = [translateColumn(_, lut_table) for _ in BRAAK5]
    braak6 = [translateColumn(_, lut_table) for _ in BRAAK6]

    # calculate cereb gray
    means, sizes = DFWeightedMean(df, cerebg)
    df.loc[:,'CEREBELLUMGREYMATTER_SIZE'] = sizes
    df.loc[:,'CEREBELLUMGREYMATTER'] = means

    # calculate braak stages
    means, sizes = DFWeightedMean(df, braak1)
    df.loc[:,'BRAAK1_SIZE'] = sizes
    df.loc[:,'BRAAK1'] = means
    means, sizes = DFWeightedMean(df, braak2)
    df.loc[:,'BRAAK2_SIZE'] = sizes
    df.loc[:,'BRAAK2'] = means
    means, sizes = DFWeightedMean(df, braak3)
    df.loc[:,'BRAAK3_SIZE'] = sizes
    df.loc[:,'BRAAK3'] = means
    means, sizes = DFWeightedMean(df, braak4)
    df.loc[:,'BRAAK4_SIZE'] = sizes
    df.loc[:,'BRAAK4'] = means
    means, sizes = DFWeightedMean(df, braak5)
    df.loc[:,'BRAAK5_SIZE'] = sizes
    df.loc[:,'BRAAK5'] = means
    means, sizes = DFWeightedMean(df, braak6)
    df.loc[:,'BRAAK6_SIZE'] = sizes
    df.loc[:,'BRAAK6'] = means

    if keys:
        df = df.loc[:,keys]

    return df


def additionalAV45Calculations(df, lut_table, keys=None):
    '''
    Do additional calculations
    If keys given, filter/sort by the list of keys before outputting
    '''
    df = df.copy()
    cerebg = [translateColumn(_, lut_table) for _ in CEREBG]
    cerebw = [translateColumn(_, lut_table) for _ in CEREBW]
    cerebl = [translateColumn(_, lut_table) for _ in CEREBL]
    cerebr = [translateColumn(_, lut_table) for _ in CEREBR]
    frontal = [translateColumn(_, lut_table) for _ in FRONTAL]
    parietal = [translateColumn(_, lut_table) for _ in PARIETAL]
    temporal = [translateColumn(_, lut_table) for _ in TEMPORAL]
    cingulate = [translateColumn(_, lut_table) for _ in CINGULATE]
    wcereb = [translateColumn(_, lut_table) for _ in WHOLECEREBELLUM]
    compref = ['ERODED_SUBCORTICALWM', 'BRAIN_STEM', 'WHOLECEREBELLUM']

    # calculate composite
    frontal_means, frontal_sizes = DFWeightedMean(df, frontal)
    parietal_means, parietal_sizes = DFWeightedMean(df, parietal)
    cingulate_means, cingulate_sizes = DFWeightedMean(df, cingulate)
    temporal_means, temporal_sizes = DFWeightedMean(df, temporal)
    df.loc[:,'FRONTAL'] = frontal_means
    df.loc[:,'FRONTAL_SIZE'] = frontal_sizes
    df.loc[:,'PARIETAL'] = parietal_means
    df.loc[:,'PARIETAL_SIZE'] = parietal_sizes
    df.loc[:,'CINGULATE'] = cingulate_means
    df.loc[:,'CINGULATE_SIZE'] = cingulate_sizes
    df.loc[:,'TEMPORAL'] = temporal_means
    df.loc[:,'TEMPORAL_SIZE'] = temporal_sizes
    df.loc[:,'COMPOSITE_SIZE'] = frontal_sizes + parietal_sizes + cingulate_sizes + temporal_sizes
    df.loc[:,'COMPOSITE'] = pd.DataFrame([frontal_means,parietal_means,cingulate_means,temporal_means]).mean()

    # calculate cereb white
    means, sizes = DFWeightedMean(df, cerebw)
    df.loc[:,'CEREBELLUMWHITEMATTER_SIZE'] = sizes
    df.loc[:,'CEREBELLUMWHITEMATTER'] = means

    # calculate cereb gray
    means, sizes = DFWeightedMean(df, cerebg)
    df.loc[:,'CEREBELLUMGREYMATTER_SIZE'] = sizes
    df.loc[:,'CEREBELLUMGREYMATTER'] = means

    # calculate whole cerebellum
    left_means, left_sizes = DFWeightedMean(df, cerebl)
    right_means, right_sizes = DFWeightedMean(df, cerebr)
    df.loc[:,'WHOLECEREBELLUM_SIZE'] = left_sizes + right_sizes
    df.loc[:,'WHOLECEREBELLUM'] = pd.DataFrame([left_means,right_means]).mean()

    # calculate composite ref
    df.loc[:,'COMPOSITE_REF_SIZE'] = df.loc[:,['%s_SIZE' % _ for _ in compref]].sum(axis=1)
    df.loc[:,'COMPOSITE_REF'] = df.loc[:,compref].mean(axis=1)

    # SUVR'S
    df.loc[:,'SUMMARYSUVR_WHOLECEREBNORM'] = df.loc[:,'COMPOSITE'].divide(df.loc[:,'WHOLECEREBELLUM'], axis=0)
    df.loc[:,'SUMMARYSUVR_COMPOSITE_REFNORM'] = df.loc[:,'COMPOSITE'].divide(df.loc[:,'COMPOSITE_REF'], axis=0)

    # thresholds
    df.loc[:,'SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF'] = (df.loc[:,'SUMMARYSUVR_WHOLECEREBNORM'] >= 1.11).astype(int)
    df.loc[:,'SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF'] = (df.loc[:,'SUMMARYSUVR_COMPOSITE_REFNORM'] >= 0.79).astype(int)

    if keys:
        df = df.loc[:,keys]

    return df

def translateColumn(orig, lut_table):
    try:
        val = int(orig)
        if val == 0:
            return 'RID'
        return lut_table[val].replace('-','_').upper()
    except Exception as e:
        print e
        pass
    return orig

def getVisitCode(rid, date, registry, cutoff=60):
    regs = sorted(registry.get(rid,[]), key=lambda x: abs(x['EXAMDATE']-date).days)
    vc = vc2 = ''
    if len(regs) > 0:
        metadata = regs[0]
        if abs(date-metadata['EXAMDATE']).days <= cutoff:
            vc = metadata.get('VISCODE','')
            vc2 = metadata.get('VISCODE2','')
        else:
            print "CLOSEST VISCODE AT %s > %s (%s, %s)" % (abs(date-metadata['EXAMDATE']).days, cutoff, rid, date)
    return (vc, vc2)

def aggregateAllRegionFiles(mean_files, size_files, lut_table, meta_file, dod=False):
    '''
    if dod:
        - RID -> SCRNO
        - VISCODE2 is removed
    '''
    tp_df = []
    try:
        visits = importScanMeta(meta_file, with_viscode=True)
    except:
        visits = None
    print mean_files
    print size_files
    for i, tp in enumerate(ALL_TP):
        mean_file = mean_files[i]
        size_file = size_files[i]
        try:
            mean_df = pd.read_csv(mean_file)
            size_df = pd.read_csv(size_file)
        except Exception as e:
            print "Can't read %s, %s" % (mean_file, size_file)
            mean_df = pd.DataFrame(columns=['RID', 'EXAMDATE'])
            size_df = pd.DataFrame(columns=['RID', 'EXAMDATE'])
        mean_df.columns = [translateColumn(_, lut_table) for _ in mean_df.columns]
        size_df.columns = [translateColumn(_, lut_table) for _ in size_df.columns]
        # add in dates + viscodes
        if visits is not None:
            for row_idx in mean_df.index:
                rid = mean_df.ix[row_idx,'RID']
                subj_visits = visits[rid]
                try:
                    tp_date = sorted(subj_visits.keys())[i]
                except:
                    msg = 'No Date: %s, %s' % (rid,tp)
                    raise Exception(msg)
                vc = subj_visits[tp_date].get('VISCODE',np.nan)
                vc2 = subj_visits[tp_date].get('VISCODE2',np.nan)
                mean_df.ix[row_idx,'EXAMDATE'] = tp_date
                mean_df.ix[row_idx,'VISCODE'] = vc
                mean_df.ix[row_idx,'VISCODE2'] = vc2
        else:
            mean_df.ix['EXAMDATE'] = pd.NaT
            mean_df.ix['VISCODE'] = np.nan
            mean_df.ix['VISCODE2'] = np.nan
        merged = mean_df.merge(size_df, on='RID', suffixes=['','_SIZE'])
        tp_df.append(merged)

    # concat visits
    all_rows = pd.concat(tp_df, axis=0)
    all_rows.reset_index(inplace=True, drop=True)
    all_rows.sort_values(by=['RID','EXAMDATE'], inplace=True)

    # Order columns
    key_columns = ['RID', 'VISCODE', 'VISCODE2', 'EXAMDATE']
    other_columns = sorted(list(set(all_rows.columns) - set(key_columns)))
    column_order = key_columns + other_columns
    all_rows = all_rows.loc[:,column_order]

    # DOD modifications
    if dod:
        all_rows.rename(columns={'RID': 'SCRNO'}, inplace=True)
        all_rows.drop('VISCODE2', axis=1, inplace=True)

    print all_rows.shape
    return all_rows


def mergeRegularWithAllRegions(regular_output, allregions_output, output_file, dod=False):
    if dod:
        id_key = 'SCRNO'
    else:
        id_key = 'RID'

    regular_df = pd.read_csv(regular_output)
    allregions_df = pd.read_csv(allregions_output)
    regular_df.set_index([id_key,'EXAMDATE'],inplace=True)
    allregions_df.set_index([id_key,'EXAMDATE'],inplace=True)
    # specify columns to merge
    merge_columns = [_ for _ in allregions_df.columns if _ not in regular_df.columns]
    merge_columns = sorted(merge_columns)
    merge_df = allregions_df[merge_columns]
    all_df = regular_df.merge(merge_df, left_index=True, right_index=True, how='outer')
    # move update_stamp to the end (if necessary)
    if 'update_stamp' in all_df.columns:
        stamps = all_df['update_stamp']
        all_df.drop('update_stamp', axis=1, inplace=True)
        all_df.insert(len(all_df.columns), 'update_stamp',stamps)
    all_df.to_csv(output_file) # , float_format='%.8f'

def findPreprocessOutputFiles(folder_name, nontp=False, allregions=False):
    '''
    Assumes preprocess outputs include all freesurfer regions (no pre-aggregation)
    '''
    addon = "_nontp" if nontp else "_tp"
    mean_keys = ["*%s%s_means*" % (tp,addon) for tp in ALL_TP]
    size_keys = ["*%s%s_roisize*" % (tp,addon) for tp in ALL_TP]
    mean_files = []
    size_files = []
    for mk in mean_keys:
        fp = os.path.join(folder_name,mk)
        try:
            mean_files.append(glob(fp)[0])
        except:
            raise Exception("Couldn't find preprocess output %s" % fp)
    for sk in size_keys:
        fp = os.path.join(folder_name,sk)
        try:
            size_files.append(glob(fp)[0])
        except:
            raise Exception("Couldn't find preprocess output %s" % fp)
    return (mean_files,size_files)

def ADNINamingConventions(input_df):
    rename_dict = {'BRAIN_STEM': 'BRAINSTEM',
                   'BRAIN_STEM_SIZE' : 'BRAINSTEM_SIZE',
                   '3RD_VENTRICLE': 'VENTRICLE_3RD',
                   '3RD_VENTRICLE_SIZE': 'VENTRICLE_3RD_SIZE',
                   '4TH_VENTRICLE': 'VENTRICLE_4TH',
                   '4TH_VENTRICLE_SIZE': 'VENTRICLE_4TH_SIZE',
                   '5TH_VENTRICLE': 'VENTRICLE_5TH',
                   '5TH_VENTRICLE_SIZE': 'VENTRICLE_5TH_SIZE'}
    output_df = input_df.rename(columns=rename_dict)
    return output_df

def CreateLONIVersions(allregions_df, fieldnames, lut_table, output_folder, allregions_filename, regular_filename, merged_filename, dod=False):
    allregions_output = os.path.join(output_folder,'LONI_%s' % allregions_filename)
    regular_output = os.path.join(output_folder,'LONI_%s' % regular_filename)
    merged_output = os.path.join(output_folder,'LONI_%s' % merged_filename)
    full_df = additionalAV45Calculations(allregions_df, lut_table, keys=fieldnames)
    ADNINamingConventions(allregions_df).to_csv(allregions_output,index=False,float_format='%.4f')
    ADNINamingConventions(full_df).to_csv(regular_output,index=False,float_format='%.4f')
    #mergeRegularWithAllRegions(regular_output, allregions_output, merged_output, dod=dod)

if __name__ == "__main__":
    # freesurfer region lookup
    lut_file = "../FreeSurferColorLUT.txt"
    lut_table = importFreesurferLookup(lut_file)

    registry_file = "../docs/ADNI/REGISTRY.csv"
    dod_registry_file = "../docs/DOD/REGISTRY.csv"
    dep_registry_file = "../docs/DEP/REGISTRY.csv"

    meta_pet = "../docs/ADNI/AV45META.csv"
    meta_tau = '../docs/ADNI/TAUMETA.csv'

    dod_meta_pet = "../docs/DOD/AV45META.csv"
    dod_meta_tau = "../docs/DOD/TAUMETA.csv"

    dep_meta_pet = "../docs/DEP/AV45META.csv"
    dep_meta_tau = "../docs/DEP/TAUMETA.csv"

    # # registry imports
    # adni_registry = importRegistry(registry_file)
    # dod_registry = importDODRegistry(dod_registry_file)

    # # pet date imports
    # adni_av45_pet_dates = importScanMeta(meta_pet, with_viscode=True)
    # adni_tau_pet_dates = importScanMeta(meta_tau, with_viscode=True)
    # dod_av45_pet_dates = importScanMeta(dod_meta_pet, with_viscode=True)
    # dod_tau_pet_dates = importScanMeta(dod_meta_tau, with_viscode=True)

    #timestamp = datetime.now().strftime('%m_%d_%y')
    timestamp = '08_10_16'

    # preprocess output folders
    preprocess_folder = '../docs/preprocess_output/%s/' % timestamp
    adni_av45_preprocess_folder = '%sadni_av45' % preprocess_folder
    dod_av45_preprocess_folder = '%sdod_adni' % preprocess_folder
    dep_av45_preprocess_folder = '%sdep_adni_av45' % preprocess_folder
    adni_tau_preprocess_folder = '%sadni_av1451' % preprocess_folder
    dod_tau_preprocess_folder = '%sdod_adni_av1451' % preprocess_folder
    dep_tau_preprocess_folder = '%sdep_adni_av1451' % preprocess_folder

    # create output folder
    output_folder = '../output/%s' % timestamp
    mkdir_p(output_folder)

    # ADNI AV45 NONTP
    regular_filename = 'UCBERKELEYAV45_%s_regular_nontp.csv' % (timestamp)
    allregions_filename = 'UCBERKELEYAV45_%s_allregions_nontp.csv' % (timestamp)
    merged_filename = 'UCBERKELEYAV45_%s_merged_nontp.csv' % (timestamp)
    regular_output = os.path.join(output_folder, regular_filename)
    allregions_output = os.path.join(output_folder, allregions_filename)
    merged_output = os.path.join(output_folder, merged_filename)
    mean_files, size_files = findPreprocessOutputFiles(adni_av45_preprocess_folder, nontp=True)
    df = aggregateAllRegionFiles(mean_files, size_files, lut_table, meta_pet)
    df.to_csv(allregions_output,index=False,float_format='%.4f',date_format='%Y-%m-%d')
    full_df = additionalAV45Calculations(df, lut_table, keys=ADNI_FIELDNAMES_EXTRA)
    full_df.to_csv(regular_output,index=False,float_format='%.4f',date_format='%Y-%m-%d')
    CreateLONIVersions(df, ADNI_FIELDNAMES, lut_table, output_folder, allregions_filename, regular_filename, merged_filename)

    # ADNI AV45 TP
    regular_output = os.path.join(output_folder, 'UCBERKELEYAV45_%s_regular_tp.csv' % (timestamp))
    allregions_output = os.path.join(output_folder, 'UCBERKELEYAV45_%s_allregions_tp.csv' % (timestamp))
    merged_output = os.path.join(output_folder, 'UCBERKELEYAV45_%s_merged_tp.csv' % (timestamp))
    mean_files, size_files = findPreprocessOutputFiles(adni_av45_preprocess_folder, nontp=False)
    df = aggregateAllRegionFiles(mean_files, size_files, lut_table, meta_pet)
    df.to_csv(allregions_output,index=False,float_format='%.4f',date_format='%Y-%m-%d')
    full_df = additionalAV45Calculations(df, lut_table, keys=ADNI_FIELDNAMES_EXTRA)
    full_df.to_csv(regular_output,index=False,float_format='%.4f',date_format='%Y-%m-%d')

    # ADNI TAU TP
    regular_filename = 'UCBERKELEYAV1451_%s_regular_tp.csv' % (timestamp)
    allregions_filename = 'UCBERKELEYAV1451_%s_allregions_tp.csv' % (timestamp)
    merged_filename = 'UCBERKELEYAV1451_%s_merged_tp.csv' % (timestamp)
    regular_output = os.path.join(output_folder, regular_filename)
    allregions_output = os.path.join(output_folder, allregions_filename)
    merged_output = os.path.join(output_folder, merged_filename)
    mean_files, size_files = findPreprocessOutputFiles(adni_tau_preprocess_folder, nontp=False)
    df = aggregateAllRegionFiles(mean_files, size_files, lut_table, meta_tau)
    df.to_csv(allregions_output,index=False,float_format='%.4f',date_format='%Y-%m-%d')
    full_df = additionalTauCalculations(df, lut_table, keys=TAU_FIELDNAMES_EXTRA)
    full_df.to_csv(regular_output,index=False,float_format='%.4f',date_format='%Y-%m-%d')
    CreateLONIVersions(df, TAU_FIELDNAMES, lut_table, output_folder, allregions_filename, regular_filename, merged_filename)

    # DOD AV45 NONTP
    regular_filename = 'UCBERKELEYAV45_DOD_%s_regular_nontp.csv' % (timestamp)
    allregions_filename = 'UCBERKELEYAV45_DOD_%s_allregions_nontp.csv' % (timestamp)
    merged_filename = 'UCBERKELEYAV45_DOD_%s_merged_nontp.csv' % (timestamp)
    regular_output = os.path.join(output_folder, regular_filename)
    allregions_output = os.path.join(output_folder, allregions_filename)
    merged_output = os.path.join(output_folder, merged_filename)
    mean_files, size_files = findPreprocessOutputFiles(dod_av45_preprocess_folder, nontp=True)
    df = aggregateAllRegionFiles(mean_files, size_files, lut_table, dod_meta_pet, dod=True)
    df.to_csv(allregions_output,index=False,float_format='%.4f',date_format='%Y-%m-%d')
    full_df = additionalAV45Calculations(df, lut_table, keys=DOD_FIELDNAMES_EXTRA)
    full_df.to_csv(regular_output,index=False,float_format='%.4f',date_format='%Y-%m-%d')
    CreateLONIVersions(df, DOD_FIELDNAMES, lut_table, output_folder, allregions_filename, regular_filename, merged_filename, dod=True)

    # DOD TAU TP
    regular_filename = 'UCBERKELEYAV1451_DOD_%s_regular_tp.csv' % (timestamp)
    allregions_filename = 'UCBERKELEYAV1451_DOD_%s_allregions_tp.csv' % (timestamp)
    merged_filename = 'UCBERKELEYAV1451_DOD_%s_merged_tp.csv' % (timestamp)
    regular_output = os.path.join(output_folder, regular_filename)
    allregions_output = os.path.join(output_folder, allregions_filename)
    merged_output = os.path.join(output_folder, merged_filename)
    mean_files, size_files = findPreprocessOutputFiles(dod_tau_preprocess_folder, nontp=False)
    df = aggregateAllRegionFiles(mean_files, size_files, lut_table, dod_meta_tau, dod=True)
    df.to_csv(allregions_output,index=False,float_format='%.4f',date_format='%Y-%m-%d')
    full_df = additionalTauCalculations(df, lut_table, keys=DOD_TAU_FIELDNAMES_EXTRA)
    full_df.to_csv(regular_output,index=False,float_format='%.4f',date_format='%Y-%m-%d')
    CreateLONIVersions(df, DOD_TAU_FIELDNAMES, lut_table, output_folder, allregions_filename, regular_filename, merged_filename, dod=True)

    # DEP AV45 TP
    regular_filename = 'UCBERKELEYAV45_DEP_%s_regular_tp.csv' % (timestamp)
    allregions_filename = 'UCBERKELEYAV45_DEP_%s_allregions_tp.csv' % (timestamp)
    merged_filename = 'UCBERKELEYAV45_DEP_%s_merged_tp.csv' % (timestamp)
    regular_output = os.path.join(output_folder, regular_filename)
    allregions_output = os.path.join(output_folder, allregions_filename)
    merged_output = os.path.join(output_folder, merged_filename)
    mean_files, size_files = findPreprocessOutputFiles(dep_av45_preprocess_folder, nontp=False)
    df = aggregateAllRegionFiles(mean_files, size_files, lut_table, dep_meta_pet, dod=False)
    df.to_csv(allregions_output,index=False,float_format='%.4f',date_format='%Y-%m-%d')
    full_df = additionalAV45Calculations(df, lut_table, keys=DEP_FIELDNAMES_EXTRA)
    full_df.to_csv(regular_output,index=False,float_format='%.4f',date_format='%Y-%m-%d')
    CreateLONIVersions(df, DEP_FIELDNAMES, lut_table, output_folder, allregions_filename, regular_filename, merged_filename, dod=False)
