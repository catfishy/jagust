import re
from collections import defaultdict
from datetime import datetime, timedelta
from scipy import stats, optimize
import numpy as np

from utils import *


def syncAPOEData(old_headers, old_lines, apoe_file, registry_file):
    apoe_by_subj = importAPOE(apoe_file)

    to_add_headers = ['APOE4BIN']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after='Notes')
    new_lines = []

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        apoe4 = 0
        if subj_row['apgen1'] == 4 or subj_row['apgen2'] == 4:
            apoe4 = 1
        return {'APOE4BIN': apoe4}

    for old_l in old_lines:
        new_data = updateLine(old_l, apoe_by_subj, extraction_fn, 
                              pid_key='SCRNO', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    return (new_headers, new_lines)

def syncDemogData(old_headers, old_lines, demog_file, registry_file):
    demog_by_subj = importDemog(demog_file)
    
    to_add_headers = ['Sex', 'Age', 'Edu']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after='Notes')
    new_lines = []

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        if int(subj_row['gender']) == 1:
            sex = 'M'
        elif int(subj_row['gender']) == 2:
            sex = 'F'
        else:
            raise Exception("Bad gender %s" % subj_row['gender'])
        age = float(subj_row['age'])
        edu = int(subj_row['edu'])
        return {'Sex': sex,
                'Age': age,
                'Edu': edu}

    for old_l in old_lines:
        new_data = updateLine(old_l, demog_by_subj, extraction_fn, 
                              pid_key='SCRNO', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    return (new_headers, new_lines)

def syncCAPSData(old_headers, old_lines, curr_file, lifetime_file):
    caps_by_subj = importCAPS(curr_file, lifetime_file)

    to_add_headers = ['CAPS_CURRSCORE', 'CAPS_LIFETIME_SCORE']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    new_lines = []

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        return {'CAPS_CURRSCORE': subj_row['curr'],
                'CAPS_LIFETIME_SCORE': subj_row['life']}

    for old_l in old_lines:
        new_data = updateLine(old_l, caps_by_subj, extraction_fn,
                              pid_key='SCRNO', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)
    return (new_headers, new_lines)


def syncWMHData(old_headers, old_lines, wmh_file):
    wmh_by_subj = importWMH(wmh_file)

    to_add_headers = ['WMH_WHITMATHYP', 'WMH_percentOfICV']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    new_lines = []

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        point = sorted(subj_row, key=lambda x: x['EXAMDATE'])[0] # take first point
        return {'WMH_percentOfICV': point['wmh_percent'],
                'WMH_WHITMATHYP': point['wmh']}

    for old_l in old_lines:
        new_data = updateLine(old_l, wmh_by_subj, extraction_fn,
                              pid_key='SCRNO', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)
    return (new_headers, new_lines)


def syncMRIData(old_headers, old_lines, mri_file):
    mri_by_subj = importDODMRI(mri_file)
    to_add_headers = ["MRIDATE"]
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    new_lines = []

    def extraction_fn(subj, subj_row, old_line, patients_pets):
        first_mri = subj_row[0]
        return {"MRIDATE": first_mri}

    for old_l in old_lines:
        new_data = updateLine(old_l, mri_by_subj, extraction_fn,
                              pid_key='SCRNO', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    return (new_headers, new_lines)


def syncGDData(old_headers, old_lines, gd_file):
    gd_by_subj = importGD(gd_file)
    to_add_headers = ['GDtotal']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    new_lines = []

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        return {'GDtotal': subj_row['gdtotal']}

    for old_l in old_lines:
        new_data = updateLine(old_l, gd_by_subj, extraction_fn,
                              pid_key='SCRNO', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)
    return (new_headers, new_lines)


def syncAntidepData(old_headers, old_lines, backmeds_file, registry_file):
    antidep_by_subj = importDODAntidep(backmeds_file, registry_file)

    to_add_headers = ['ANTIDEP_USE','WHICH_ANTIDEP','SSRI']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    new_lines = []

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        bl = sorted(subj_row, key=lambda x: x['EXAMDATE'])[0]
        return {'ANTIDEP_USE': 1 if bl['ANTIDEP_USE'] else 0,
                'WHICH_ANTIDEP': bl['WHICH_ANTIDEP'],
                'SSRI': 1 if bl['SSRI'] else 0}

    for i, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, antidep_by_subj, extraction_fn, 
                              pid_key='SCRNO', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    return (new_headers, new_lines)


def syncCSFData(old_headers, old_lines, csf_file, registry_file):
    registry = importDODRegistry(registry_file)
    csf_by_subj = importCSF([csf_file], registry)

    to_add_headers = ['CSF_abeta', 'CSF_tau', 'CSF_ptau']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    new_lines = []

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        point = sorted(subj_row, key=lambda x: x['EXAMDATE'])[0]
        return {'CSF_abeta': point['abeta'],
                'CSF_tau': point['tau'],
                'CSF_ptau': point['ptau']}

    for i, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, csf_by_subj, extraction_fn, 
                              pid_key='SCRNO', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    return (new_headers, new_lines)

def syncADASData(old_headers, old_lines, adas_file, registry_file):
    registry = importDODRegistry(registry_file)
    adas_by_subj = importADASCog(None, adas_file, registry=registry)
    
    to_add_headers = ['ADAS_TOTSCORE']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    new_lines = []

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        point = sorted(subj_row, key=lambda x: x['EXAMDATE'])[0]
        return {'ADAS_TOTSCORE': point['TOTSCORE']}

    for i, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, adas_by_subj, extraction_fn, 
                              pid_key='SCRNO', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    return (new_headers, new_lines)

def syncAVLTData(old_headers, old_lines, avlt_file, registry_file):
    registry = importDODRegistry(registry_file)
    avlt_by_subj = importAVLT(avlt_file, registry=registry)

    to_add_headers = ['AVLT_total_6mths_Examdate']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    new_lines = []

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        sorted_points = sorted(subj_row, key=lambda x: x['EXAMDATE'])

        bl_av45, av45_2 = getAV45Dates(old_l, patient_pets=None)
        closest_1, closest_2, closest_3 = getClosestToScans(sorted_points, bl_av45, av45_2, None, 
                                                            day_limit=int(365/2), date_key='EXAMDATE')
        data = {}

        # Closest to AV45 value
        data['AVLT_total_6mths_Examdate'] = closest_1['TOTS'] if closest_1 else None

        # Long values
        long_values = extractLongitudinalFields(sorted_points, 'TOTS', 'AVLT.')
        long_dates = extractLongitudinalFields(sorted_points, 'EXAMDATE', 'AVLT_postAV45.')
        long_dates_diff = {k:(v-bl_av45).days for k,v in long_dates.iteritems()}
        long_dates_filtered = {k:(v/365.0) if v >= -90.0 else None for k,v in long_dates_diff.iteritems()}
        data.update(long_values)
        data.update(long_dates_filtered)

        # Slope values
        slope_pts = [(t,v) for t,v in zip(long_dates_filtered,long_values) if t is not None]
        print slope_pts
        data['AVLT_slope'] = slope(slope_pts)

        return data

    for i, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, avlt_by_subj, extraction_fn, 
                              pid_key='SCRNO', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    return (new_headers, new_lines)

def syncADNIMasterData(old_headers, old_lines, adnimaster_file):
    master_by_subj = importMaster(adnimaster_file)

    new_headers = old_headers
    new_lines = []

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        '''
        fill in the following keys:

        Sex -> Gender
        Age
        Edu -> Edu.(Yrs)
        APOE4BIN
        ADAS_TOTSCORE
        AVLT_total_6mths_Examdate
        CAPS_CURRSCORE
        CAPS_LIFETIME_SCORE
        WMH_WHITMATHYP
        WMH_percentOfICV
        GDtotal
        MRIDATE
        CSF_abeta
        CSF_tau
        CSF_ptau
        '''
        gender = subj_row['Gender']
        try:
            if int(subj_row['Gender']) == 1:
                gender = 'M'
        except:
            pass
        try:
            if int(subj_row['Gender']) == 2:
                gender = 'F'
        except:
            pass
        new_data = {'Sex': gender,
                    'Age': subj_row['Age@AV45'],
                    'Edu': subj_row['Edu.(Yrs)'],
                    'APOE4BIN': subj_row['APOE4_BIN'],
                    'ADAS_TOTSCORE': subj_row['ADAS_3MTH_AV45'],
                    'AVLT_total_6mths_Examdate': subj_row['AVLT_3MTHS_AV45'],
                    'CAPS_CURRSCORE': '',
                    'CAPS_LIFETIME_SCORE': '',
                    'WMH_WHITMATHYP': subj_row['WMH_WHITMATHYP.1'],
                    'WMH_percentOfICV': subj_row['WMH_percentOfICV.1'],
                    'GDtotal': subj_row['GD_AV45_6MTHS'],
                    'MRIDATE': '',
                    'CSF_abeta': subj_row['CSF_ABETA.1'],
                    'CSF_tau': subj_row['CSF_TAU.1'],
                    'CSF_ptau': subj_row['CSF_PTAU.1']
                    }
        return new_data

    for i, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, master_by_subj, extraction_fn, 
                              pid_key='SCRNO', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    return (new_headers, new_lines)

def syncDiagData(master_df, diag_file):
    diag_df = importADNIDiagnosis(diag_file, as_df=True)
    timepoints = 3

    # header order 
    headers = ['Diag_%s' % (_+1,) for _ in range(timepoints)]
    headers += ["Diag_Date_%s" % (_+1,) for _ in range(timepoints)]

    # extraction fn
    def extraction_fn(scrno, subj_rows):
        if len(subj_rows.index) > timepoints:
            raise Exception("Raise # Diag Timepoints to %s" % len(subj_rows.index))
        subj_rows.sort_values(by='EXAMDATE', inplace=True)
        subj_rows['SCRNO'] = scrno
        diag_long = pivotSubjectGroupValues(subj_rows, 'SCRNO','diag','Diag_')
        date_long = pivotSubjectGroupValues(subj_rows, 'SCRNO','EXAMDATE','Diag_Date_')
        merge_df = diag_long.merge(date_long,left_index=True,right_index=True)
        agg_row = merge_df.reset_index().iloc[0]
        return agg_row

    parsed_df = parseSubjectGroups(diag_df, extraction_fn)
    parsed_df.set_index('SCRNO',inplace=True)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after='Notes')
    return master_df


def syncStudyData(old_headers, old_lines, elig_file):
    elig_by_subj = importDODEligibility(elig_file)

    to_add_headers = ['PTGroup', 'Study', 'GroupNum']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after='Notes')
    new_lines = []

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        cohort = subj_row['cohort']
        new_data = {'PTGroup': '',
                    'Study': 1,
                    'GroupNum': ''}
        if cohort == 1:
            new_data['PTGroup'] = 'PTSD'
            new_data['GroupNum'] = 1
        elif cohort == 2:
            new_data['PTGroup'] = 'TBI'
            new_data['GroupNum'] = 2
        elif cohort == 3:
            new_data['PTGroup'] = 'C'
            new_data['GroupNum'] = 3
        elif cohort == 4:
            new_data['PTGroup'] = 'TBI+PTSD'
            new_data['GroupNum'] = 4
        elif cohort == '':
            # decide what to do here
            pass
        return new_data

    for i, old_l in enumerate(old_lines):
        new_data = {'PTGroup': '',
                    'Study': '',
                    'GroupNum': ''}
        subj = int(old_l['SCRNO'])
        if subj < 6000:
            new_data = {'PTGroup': 'ADNI C',
                        'Study': 0,
                        'GroupNum': 5}
        else:
            new_data = updateLine(old_l, elig_by_subj, extraction_fn, 
                          pid_key='SCRNO', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    return (new_headers, new_lines)

def syncAV1451Data(master_df, av1451_file):
    av1451_df = importAV1451(av1451_file, as_df=True)

    # create header order
    timepoints = 1
    row_indices = range(1,timepoints+1)
    headers = []
    headers += ['AV1451_%s_DATE' % i for i in row_indices]
    headers += ['AV1451_%s_Braak12_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_Braak34_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_Braak56_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_LEFT_CAUDATE_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_RIGHT_CAUDATE_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_LEFT_CHOROID_PLEXUS_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_RIGHT_CHOROID_PLEXUS_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_CEREBELLUMGREYMATTER' % i for i in row_indices]
    headers += ['AV1451_%s_BRAIN_STEM' % i for i in row_indices]
    headers += ['AV1451_%s_WHOLECEREBELLUM' % i for i in row_indices]

    # create extraction function
    def extraction_fn(scrno, subj_rows):
        new_data = {'SCRNO': scrno}
        subj_rows.sort_values(by='EXAMDATE', inplace=True)
        if len(subj_rows.index) > timepoints:
            raise Exception("Raise # AV1451 Timepoints to > %s" % timepoints)

        for i, (idx, row) in enumerate(subj_rows.iterrows()):
            row_i = i + 1

            new_data['AV1451_%s_DATE' % row_i] = row['EXAMDATE']
            new_data['AV1451_%s_BRAIN_STEM' % row_i] = float(row['BRAIN_STEM'])

            # get whole cerebellum
            cerebellum_parts = ['RIGHT_CEREBELLUM_CORTEX', 'RIGHT_CEREBELLUM_WHITE_MATTER', 'LEFT_CEREBELLUM_CORTEX', 'LEFT_CEREBELLUM_WHITE_MATTER']
            weights_values = [(float(row["%s_SIZE" % _]), float(row[_])) for _ in cerebellum_parts]
            wcereb = weightedMean(weights_values)
            new_data['AV1451_%s_WHOLECEREBELLUM' % row_i] = wcereb

            # get cerebellar gray
            cerebg_parts = ['RIGHT_CEREBELLUM_CORTEX', 'LEFT_CEREBELLUM_CORTEX']
            weights_values = [(float(row["%s_SIZE" % _]), float(row[_])) for _ in cerebg_parts]
            cerebg = weightedMean(weights_values)
            new_data['AV1451_%s_CEREBELLUMGREYMATTER' % row_i] = cerebg

            braak1 = (row['BRAAK1_SIZE'],row['BRAAK1'])
            braak2 = (row['BRAAK2_SIZE'],row['BRAAK2'])
            braak3 = (row['BRAAK3_SIZE'],row['BRAAK3'])
            braak4 = (row['BRAAK4_SIZE'],row['BRAAK4'])
            braak5 = (row['BRAAK5_SIZE'],row['BRAAK5'])
            braak6 = (row['BRAAK6_SIZE'],row['BRAAK6'])
            braak12 = weightedMean([braak1,braak2]) / cerebg
            braak34 = weightedMean([braak3,braak4]) / cerebg
            braak56 = weightedMean([braak5,braak6]) / cerebg

            # region suvrs
            new_data['AV1451_%s_Braak12_CerebGray' % row_i] = braak12
            new_data['AV1451_%s_Braak34_CerebGray' % row_i] = braak34
            new_data['AV1451_%s_Braak56_CerebGray' % row_i] = braak56
            new_data['AV1451_%s_LEFT_CAUDATE_CerebGray' % row_i] = float(row['LEFT_CAUDATE'])/cerebg
            new_data['AV1451_%s_RIGHT_CAUDATE_CerebGray' % row_i] = float(row['RIGHT_CAUDATE'])/cerebg
            new_data['AV1451_%s_LEFT_PUTAMEN_CerebGray' % row_i] = float(row['LEFT_PUTAMEN'])/cerebg
            new_data['AV1451_%s_RIGHT_PUTAMEN_CerebGray' % row_i] = float(row['RIGHT_PUTAMEN'])/cerebg
            new_data['AV1451_%s_LEFT_CHOROID_PLEXUS_CerebGray' % row_i] = float(row['LEFT_CHOROID_PLEXUS'])/cerebg
            new_data['AV1451_%s_RIGHT_CHOROID_PLEXUS_CerebGray' % row_i] = float(row['RIGHT_CHOROID_PLEXUS'])/cerebg 

        return pd.Series(new_data)

    parsed_df = parseSubjectGroups(av1451_df, extraction_fn)
    parsed_df.set_index('SCRNO',inplace=True)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df


def syncAV45Data(master_df, av45_file, diags_df):
    av45_df = importAV45(av45_file, as_df=True)
    adni_controls_df = set(diags_df[diags_df.diag.str.match('N')].index)
    # create header order
    row_indices = range(1,3)
    headers = []
    for r in row_indices:
        # Dates
        headers.append('AV45_%s_EXAMDATE' % r)
        # SUVRs
        headers.append('AV45_%s_Left-Putamen' % r)
        headers.append('AV45_%s_Right-Putamen' % r)
        headers.append('AV45_%s_Left-Caudate' % r)
        headers.append('AV45_%s_Right-Caudate' % r)
        headers.append('AV45_%s_Left-Pallidum' % r)
        headers.append('AV45_%s_Right-Pallidum' % r)
        headers.append('AV45_%s_Left_BG_avg' % r)
        headers.append('AV45_%s_Right_BG_avg' % r)
        headers.append('AV45_%s_BG_avg' % r)
        headers.append('AV45_%s_comp/wcerb' % r)
        headers.append('AV45_%s_wcerb_bin' % r)
        headers.append('AV45_%s_comp/brainstem' % r)
        headers.append('AV45_%s_comp/bigref' % r)
        headers.append('AV45_%s_comp/cerbg' % r)
        headers.append('AV45_%s_comp/wm70' % r)
        # Assymetry values
        headers.append('AV45_%s_frontal_asymmetry_negvalue_means_R<L' % r)
        headers.append('AV45_%s_cingulate_asymmetry_negvalue_means_R<L' % r)
        headers.append('AV45_%s_parietal_asymmetry_negvalue_means_R<L' % r)
        headers.append('AV45_%s_temporal_asymmetry_negvalue_means_R<L' % r)
        #headers.append('AV45_%s_asym_summary_absvals_negvalue_means_R<L' % r)
        headers.append('AV45_%s_frontal_MR_asymmetry' % r)
        headers.append('AV45_%s_cingulate_MR_asymmetry' % r)
        headers.append('AV45_%s_parietal_MR_asymmetry' % r)
        headers.append('AV45_%s_temporal__MR_asymmetry' % r)
        headers.append('AV45_%s_ventrical_MR_asymmetry' % r)
    # parse AV45 results
    parsed_df = parseSubjectGroups(av45_df, parseAV45SubjectRows)
    parsed_df.set_index('SCRNO',inplace=True)
    # filter old subjects
    valid_adni_subj = set(master_df.index) &  adni_controls_df
    old_subj = valid_adni_subj | set(master_df[master_df.index >= 6000].index)
    master_df = master_df.loc[list(old_subj),:]
    # update
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df

def parseAV45SubjectRows(scrno, subj_rows):
    # reference dicts
    left_frontal_keys = ['CTX_LH_CAUDALMIDDLEFRONTAL', 'CTX_LH_LATERALORBITOFRONTAL', 'CTX_LH_MEDIALORBITOFRONTAL', 'CTX_LH_PARSOPERCULARIS', 
                    'CTX_LH_PARSORBITALIS', 'CTX_LH_PARSTRIANGULARIS', 'CTX_LH_ROSTRALMIDDLEFRONTAL', 'CTX_LH_SUPERIORFRONTAL', 'CTX_LH_FRONTALPOLE']
    right_frontal_keys = ['CTX_RH_CAUDALMIDDLEFRONTAL', 'CTX_RH_LATERALORBITOFRONTAL', 'CTX_RH_MEDIALORBITOFRONTAL', 'CTX_RH_PARSOPERCULARIS', 
                     'CTX_RH_PARSORBITALIS', 'CTX_RH_PARSTRIANGULARIS', 'CTX_RH_ROSTRALMIDDLEFRONTAL', 'CTX_RH_SUPERIORFRONTAL', 'CTX_RH_FRONTALPOLE']
    left_cingulate_keys = ['CTX_LH_CAUDALANTERIORCINGULATE', 'CTX_LH_ISTHMUSCINGULATE', 'CTX_LH_POSTERIORCINGULATE', 'CTX_LH_ROSTRALANTERIORCINGULATE']
    right_cingulate_keys = ['CTX_RH_CAUDALANTERIORCINGULATE', 'CTX_RH_ISTHMUSCINGULATE', 'CTX_RH_POSTERIORCINGULATE', 'CTX_RH_ROSTRALANTERIORCINGULATE']
    left_parietal_keys = ['CTX_LH_INFERIORPARIETAL', 'CTX_LH_PRECUNEUS', 'CTX_LH_SUPERIORPARIETAL', 'CTX_LH_SUPRAMARGINAL']
    right_parietal_keys = ['CTX_RH_INFERIORPARIETAL', 'CTX_RH_PRECUNEUS', 'CTX_RH_SUPERIORPARIETAL', 'CTX_RH_SUPRAMARGINAL']
    left_temporal_keys = ['CTX_LH_MIDDLETEMPORAL', 'CTX_LH_SUPERIORTEMPORAL']
    right_temporal_keys = ['CTX_RH_MIDDLETEMPORAL', 'CTX_RH_SUPERIORTEMPORAL']
    left_bg_keys = ['LEFT_PUTAMEN', 'LEFT_CAUDATE', 'LEFT_PALLIDUM']
    right_bg_keys = ['RIGHT_PUTAMEN', 'RIGHT_CAUDATE', 'RIGHT_PALLIDUM']
    left_ventrical_keys = []
    right_ventrical_keys = []

    data = {'SCRNO': scrno}
    subj_rows.sort_values(by='EXAMDATE', inplace=True)

    for i, (idx, point) in enumerate(subj_rows.iterrows()):
        examdate = point['EXAMDATE']
        # extract necessary values
        wcereb = float(point['WHOLECEREBELLUM'])
        cerebg = float(point['CEREBELLUMGREYMATTER'])
        compositeroi = float(point['COMPOSITE'])
        bigref = float(point['COMPOSITE_REF'])
        wm70 = float(point['ERODED_SUBCORTICALWM'])
        brainstem = float(point['BRAIN_STEM'])
        cingulate = float(point['CINGULATE'])
        frontal = float(point['FRONTAL'])
        parietal = float(point['PARIETAL'])
        temporal = float(point['TEMPORAL'])
        leftputamen = float(point['LEFT_PUTAMEN'])
        rightputamen = float(point['RIGHT_PUTAMEN'])
        leftcaudate = float(point['LEFT_CAUDATE'])
        rightcaudate = float(point['RIGHT_CAUDATE'])
        leftpallidum = float(point['LEFT_PALLIDUM'])
        rightpallidum = float(point['RIGHT_PALLIDUM'])
        left_frontal = np.average([float(point[k]) for k in left_frontal_keys], weights=[float(point["%s_SIZE" % k]) for k in left_frontal_keys])
        right_frontal = np.average([float(point[k]) for k in right_frontal_keys], weights=[float(point["%s_SIZE" % k]) for k in right_frontal_keys])
        left_cingulate = np.average([float(point[k]) for k in left_cingulate_keys], weights=[float(point["%s_SIZE" % k]) for k in left_cingulate_keys])
        right_cingulate = np.average([float(point[k]) for k in right_cingulate_keys], weights=[float(point["%s_SIZE" % k]) for k in right_cingulate_keys])
        left_parietal = np.average([float(point[k]) for k in left_parietal_keys], weights=[float(point["%s_SIZE" % k]) for k in left_parietal_keys])
        right_parietal = np.average([float(point[k]) for k in right_parietal_keys], weights=[float(point["%s_SIZE" % k]) for k in right_parietal_keys])
        left_temporal = np.average([float(point[k]) for k in left_temporal_keys], weights=[float(point["%s_SIZE" % k]) for k in left_temporal_keys])
        right_temporal = np.average([float(point[k]) for k in right_temporal_keys], weights=[float(point["%s_SIZE" % k]) for k in right_temporal_keys])
        left_bg = np.average([float(point[k]) for k in left_bg_keys])
        right_bg = np.average([float(point[k]) for k in right_bg_keys])
        left_frontal_size = np.sum([float(point["%s_SIZE" % k]) for k in left_frontal_keys])
        right_frontal_size = np.sum([float(point["%s_SIZE" % k]) for k in right_frontal_keys])
        left_cingulate_size = np.sum([float(point["%s_SIZE" % k]) for k in left_cingulate_keys])
        right_cingulate_size = np.sum([float(point["%s_SIZE" % k]) for k in right_cingulate_keys])
        left_parietal_size = np.sum([float(point["%s_SIZE" % k]) for k in left_parietal_keys])
        right_parietal_size = np.sum([float(point["%s_SIZE" % k]) for k in right_parietal_keys])
        left_temporal_size = np.sum([float(point["%s_SIZE" % k]) for k in left_temporal_keys])
        right_temporal_size = np.sum([float(point["%s_SIZE" % k]) for k in right_temporal_keys])
        left_ventrical_size = np.sum([float(point["%s_SIZE" % k]) for k in left_ventrical_keys])
        right_ventrical_size = np.sum([float(point["%s_SIZE" % k]) for k in right_ventrical_keys])

        # Dates
        data['AV45_%s_EXAMDATE' % (i+1)] = examdate
        # SUVR values
        data['AV45_%s_Left-Putamen' % (i+1)] = leftputamen/wm70
        data['AV45_%s_Right-Putamen' % (i+1)] = rightputamen/wm70
        data['AV45_%s_Left-Caudate' % (i+1)] = leftcaudate/wm70
        data['AV45_%s_Right-Caudate' % (i+1)] = rightcaudate/wm70
        data['AV45_%s_Left-Pallidum' % (i+1)] = leftpallidum/wm70
        data['AV45_%s_Right-Pallidum' % (i+1)] = rightpallidum/wm70
        data['AV45_%s_Left_BG_avg' % (i+1)] = left_bg/wm70
        data['AV45_%s_Right_BG_avg' % (i+1)] = right_bg/wm70
        data['AV45_%s_BG_avg' % (i+1)] = np.mean([left_bg,right_bg])/wm70
        data['AV45_%s_comp/wcerb' % (i+1)] = compositeroi / wcereb
        data['AV45_%s_wcerb_bin' % (i+1)] = 1 if (compositeroi/wcereb) >= 1.11 else 0 # WHATS THE REAL THRESHOLD
        data['AV45_%s_comp/brainstem' % (i+1)] = compositeroi / brainstem
        data['AV45_%s_comp/bigref' % (i+1)] = compositeroi / bigref
        data['AV45_%s_comp/cerbg' % (i+1)] = compositeroi / cerebg
        data['AV45_%s_comp/wm70' % (i+1)] = compositeroi / wm70
        # Asymmetry values
        data['AV45_%s_frontal_asymmetry_negvalue_means_R<L' % (i+1)] = asymIndex(left_frontal, right_frontal)
        data['AV45_%s_cingulate_asymmetry_negvalue_means_R<L' % (i+1)] = asymIndex(left_cingulate, right_cingulate)
        data['AV45_%s_parietal_asymmetry_negvalue_means_R<L' % (i+1)] = asymIndex(left_parietal, right_parietal)
        data['AV45_%s_temporal_asymmetry_negvalue_means_R<L' % (i+1)] = asymIndex(left_temporal, right_temporal)
        left_summary = np.mean([left_frontal, left_cingulate, left_parietal, left_temporal])
        right_summary = np.mean([right_frontal, right_cingulate, right_parietal, right_temporal])
        '''
        data['AV45_%s_asym_summary_absvals_negvalue_means_R<L' % (i+1)] = np.mean([abs(asymIndex(left_frontal, right_frontal)),
                                                                                   abs(asymIndex(left_cingulate, right_cingulate)),
                                                                                   abs(asymIndex(left_parietal, right_parietal)),
                                                                                   abs(asymIndex(left_temporal, right_temporal))])
        '''
        data['AV45_%s_frontal_MR_asymmetry' % (i+1)] = asymIndex(left_frontal_size, right_frontal_size)
        data['AV45_%s_cingulate_MR_asymmetry' % (i+1)] = asymIndex(left_cingulate_size, right_cingulate_size)
        data['AV45_%s_parietal_MR_asymmetry' % (i+1)] = asymIndex(left_parietal_size, right_parietal_size)
        data['AV45_%s_temporal__MR_asymmetry' % (i+1)] = asymIndex(left_temporal_size, right_temporal_size)
        data['AV45_%s_ventrical_MR_asymmetry' % (i+1)] = asymIndex(left_ventrical_size, right_ventrical_size)
    
    return pd.Series(data)



def asymIndex(left, right):
    return (left-right)/np.mean([left,right])

def getAV45Dates(old_l, patient_pets=None):
    if patient_pets is None:
        patient_pets = []
    # Get AV45 Scan dates
    if old_l.get('AV45_1_EXAMDATE','') != '':
        bl_av45 = datetime.strptime(old_l['AV45_1_EXAMDATE'], '%m/%d/%y')
    elif len(patient_pets) >= 1:
        bl_av45 = patient_pets[0]
    else:
        bl_av45 = None
    if old_l.get('AV45_2_EXAMDATE','') != '':
        av45_2 = datetime.strptime(old_l['AV45_2_EXAMDATE'], '%m/%d/%y')
    elif len(patient_pets) >= 2:
        av45_2 = patient_pets[1]
    else:
        av45_2 = None
    return (bl_av45, av45_2)


def runPipeline():
    # syncing pipeline
    try:
        master_df = pd.read_csv(master_file)
        master_df.set_index('SCRNO',inplace=True)
        print len(master_df.index)
    except:
        master_df = pd.DataFrame(columns=['Notes'])
        new_lines = []

    # get diagnoses
    diags_df = extractDiagnosesFromMasterData(importMaster(av45_master_file, as_df=True))
    registry = importDODRegistry(registry_file)

    print "\nSYNCING AV45\n"
    master_df = syncAV45Data(master_df, av45_file, diags_df) # adds new patients
    print "\nSYNCING AV1451\n"
    master_df = syncAV1451Data(master_df, av1451_file) # adds new patients
    print "\nSYNCING DIAG\n"
    master_df = syncDiagData(master_df, diag_file)


    print master_df
    return

    print "\nSYNCING DIAG\n"
    new_headers, new_lines = syncDiagData(new_headers, new_lines, diag_file)
    print "\nSYNCING ANTIDEP\n"
    new_headers, new_lines = syncAntidepData(new_headers, new_lines, backmeds_file, registry_file)
    print "\nSYNCING APOE\n"
    new_headers, new_lines = syncAPOEData(new_headers, new_lines, apoe_file, registry_file)
    print "\nSYNCING ADAS\n"
    new_headers, new_lines = syncADASData(new_headers, new_lines, adas_file, registry_file)
    print "\nSYNCING AVLT\n"
    new_headers, new_lines = syncAVLTData(new_headers, new_lines, avlt_file, registry_file)
    print "\nSYNCING DEMOG\n"
    new_headers, new_lines = syncDemogData(new_headers, new_lines, demog_file, registry_file)
    print "\nSYNCING CAPS\n"
    new_headers, new_lines = syncCAPSData(new_headers, new_lines, caps_curr_file, caps_lifetime_file)
    print "\nSYNCING WMH\n"
    new_headers, new_lines = syncWMHData(new_headers, new_lines, wmh_file)
    print "\nSYNCING GD\n"
    new_headers, new_lines = syncGDData(new_headers, new_lines, gd_file)
    print "\nSYNCING MRI\n"
    new_headers, new_lines = syncMRIData(new_headers, new_lines, mri_file)
    print "\nSYNCING STUDY\n"
    new_headers, new_lines = syncStudyData(new_headers, new_lines, elig_file)
    print "\nSYNCING CSF\n"
    new_headers, new_lines = syncCSFData(new_headers, new_lines, csf_file, registry_file)
    print "\nSYNCING ADNI MASTER\n"
    new_headers, new_lines = syncADNIMasterData(new_headers, new_lines, av45_master_file)

    print "\nDUMPING CSV\n"
    dumpDFtoCSV(master_df,output_file,decimal_places=2)

if __name__ == '__main__':
    now = datetime.now()

    # IO files
    master_file = "" # not used
    output_file = "../DOD_DATA_%s.csv" % (now.strftime("%m_%d_%y"))

    # ADNI master file
    av45_master_file = '../FDG_AV45_COGdata/FDG_AV45_COGdata_01_19_16.csv'
    # AV45 file
    av45_file = '../output/04_05_16/UCBERKELEYAV45_DOD_04_05_16_regular_nontp.csv'
    # AV1451 file
    av1451_file = '../output/04_05_16/UCBERKELEYAV1451_DOD_04_05_16_regular_tp.csv'

    # Registry file
    registry_file = "../docs/DOD/REGISTRY.csv"
    # AVLT file
    avlt_file = "../docs/DOD/NEUROBAT.csv"
    # ADAS file
    adas_file = "../docs/DOD/ADAS.csv"
    # Demog file
    demog_file = "../docs/DOD/PTDEMOG.csv"
    # CSF file
    csf_file = "../docs/DOD/UPENNBIOMK.csv"
    # APOE file
    apoe_file = "../docs/DOD/APOERES.csv"
    # CAPS files
    caps_curr_file = "../docs/DOD/CAPSCURR.csv"
    caps_lifetime_file = "../docs/DOD/CAPSLIFE.csv"
    # WMH file
    wmh_file = "../docs/DOD/WMH_VOLUMES.csv"
    # GD file
    gd_file = "../docs/DOD/GDSCALE.csv"
    # MRi file
    mri_file = "../docs/DOD/MRIMETA.csv"
    # Elig file
    elig_file = "../docs/DOD/VAELG.csv"
    # diag file
    diag_file = "../docs/DOD/DXSUM.csv"
    # antidep file
    backmeds_file = '../docs/DOD/BACKMEDS.csv'

    runPipeline()
