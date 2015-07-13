import re
from collections import defaultdict
from datetime import datetime, timedelta
from scipy import stats, optimize
import numpy as np

from utils import *


def syncAV45Data(old_headers, old_lines, av45_file, registry_file, dump_to=None):
    registry = importDODRegistry(registry_file)
    av45_by_subj = importAV45(av45_file, registry=registry):

    new_headers = None
    new_lines = []
    old_subjects = set()
    for old_l in old_lines:
        # get subject ID
        try:
            subj = int(old_l['PID'])
        except Exception as e:
            continue
        old_subjects.add(subj)
        if subj not in av45_by_subj:
            print "No AV45 data for %s" % subj
            new_lines.append(old_l)
            continue
        updated_headers, new_data = parseAV45Entries(old_headers, av45_by_subj[subj])
        if new_headers if None:
            new_headers = updated_headers
        new_data = convertToCSVDataType(new_data, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # find new subjects - in av45 file but not in old master csv file
    new_subjects = list(set(av45_by_subj.keys()) - old_subjects)
    print "New AV45 subjects: %s" % new_subjects
    for ns in new_subjects:
        av45_data = av45_by_subj[ns]
        old_l = {k: '' for k in new_headers}
        old_l['PID'] = str(ns)
        updated_headers, new_data = parseAV45Entries(old_headers, av45_data)
        if new_headers if None:
            new_headers = updated_headers
        new_data = convertToCSVDataType(new_data, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def parseAV45Entries(old_headers, subj_rows):
    subj_rows = sorted(subj_rows, key=lambda x: x['EXAMDATE'])

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

    # Rearrange the headers
    row_indices = range(1,3)
    to_add_headers = []
    for r in row_indices:
        # Dates
        to_add_headers.append('AV45_%s_EXAMDATE' % r)
        # SUVRs
        to_add_headers.append('AV45_%s_Left-Putamen' % r)
        to_add_headers.append('AV45_%s_Right-Putamen' % r)
        to_add_headers.append('AV45_%s_Left-Caudate' % r)
        to_add_headers.append('AV45_%s_Right-Caudate' % r)
        to_add_headers.append('AV45_%s_Left-Pallidum' % r)
        to_add_headers.append('AV45_%s_Right-Pallidum' % r)
        to_add_headers.append('AV45_%s_Left_BG_avg' % r)
        to_add_headers.append('AV45_%s_Right_BG_avg' % r)
        to_add_headers.append('AV45_%s_BG_avg' % r)
        to_add_headers.append('AV45_%s_comp/wcerb' % r)
        to_add_headers.append('AV45_%s_wcerb_bin' % r)
        to_add_headers.append('AV45_%s_comp/brainstem' % r)
        to_add_headers.append('AV45_%s_comp/bigref' % r)
        to_add_headers.append('AV45_%s_comp/cerbg' % r)
        to_add_headers.append('AV45_%s_comp/wm70' % r)
        # Assymetry values
        to_add_headers.append('AV45_%s_frontal_asymmetry_negvalue_means_R<L' % r)
        to_add_headers.append('AV45_%s_cingulate_asymmetry_negvalue_means_R<L' % r)
        to_add_headers.append('AV45_%s_parietal_asymmetry_negvalue_means_R<L' % r)
        to_add_headers.append('AV45_%s_temporal_asymmetry_negvalue_means_R<L' % r)
        to_add_headers.append('AV45_%s_asym_summary_absvals_negvalue_means_R<L' % r)
        to_add_headers.append('AV45_%s_frontal_MR_asymmetry' % r)
        to_add_headers.append('AV45_%s_cingulate_MR_asymmetry' % r)
        to_add_headers.append('AV45_%s_parietal_MR_asymmetry' % r)
        to_add_headers.append('AV45_%s_temporal__MR_asymmetry' % r)
        to_add_headers.append('AV45_%s_ventrical_MR_asymmetry' % r)

    data = {}
    for i, point in enumerate(subj_rows):
        examdate = point['EXAMDATE']
        # extract necessary values
        cerebw = float(point['CEREBELLUMWHITEMATTER'])
        wcereb = float(point['WHOLECEREBELLUM'])
        cerebg = float(point['CEREBELLUMGREYMATTER'])
        compositeroi = float(point['COMPOSITE'])
        bigref = float(point['COMPOSITE_REF'])
        wm70 = float(point['ERODED_SUBCORTICALWM'])
        brainstem = float(point['BRAINSTEM'])
        cingulate = float(point['CINGULATE'])
        frontal = float(point['FRONTAL'])
        parietal = float(point['PARIETAL'])
        temporal = float(point['TEMPORAL'])
        leftputamen = float(point['LEFT-PUTAMEN'])
        rightputamen = float(point['RIGHT-PUTAMEN'])
        leftcaudate = float(point['LEFT-CAUDATE'])
        rightcaudate = float(point['RIGHT-CAUDATE'])
        leftpallidum = float(point['LEFT-PALLIDUM'])
        rightpallidum = float(point['RIGHT-PALLIDUM'])
        left_frontal = np.mean([float(point[k]) for k in left_frontal_keys])
        right_frontal = np.mean([float(point[k]) for k in right_frontal_keys])
        left_cingulate = np.mean([float(point[k]) for k in left_cingulate_keys])
        right_cingulate = np.mean([float(point[k]) for k in right_cingulate_keys])
        left_parietal = np.mean([float(point[k]) for k in left_parietal_keys])
        right_parietal = np.mean([float(point[k]) for k in right_parietal_keys])
        left_temporal = np.mean([float(point[k]) for k in left_temporal_keys])
        right_temporal = np.mean([float(point[k]) for k in right_temporal_keys])
     
        # Dates
        data['AV45_%s_EXAMDATE' % (i+1)] = examdate
        # SUVR values
        data['AV45_%s_Left-Putamen' % (i+1)] = leftputamen
        data['AV45_%s_Right-Putamen' % (i+1)] = rightputamen
        data['AV45_%s_Left-Caudate' % (i+1)] = leftcaudate
        data['AV45_%s_Right-Caudate' % (i+1)] = rightcaudate
        data['AV45_%s_Left-Pallidum' % (i+1)] = leftpallidum
        data['AV45_%s_Right-Pallidum' % (i+1)] = rightpallidum
        data['AV45_%s_Left_BG_avg' % (i+1)] = ''
        data['AV45_%s_Right_BG_avg' % (i+1)] = ''
        data['AV45_%s_BG_avg' % (i+1)] = ''
        data['AV45_%s_comp/wcerb' % (i+1)] = compositeroi / wcereb
        data['AV45_%s_wcerb_bin' % (i+1)] = 1 if data['comp/wcerb'] >= 1.11 else 0 # WHATS THE REAL THRESHOLD
        data['AV45_%s_comp/brainstem' % (i+1)] = compositeroi / brainstem
        data['AV45_%s_comp/bigref' % (i+1)] = compositeroi / bigref
        data['AV45_%s_comp/cerbg' % (i+1)] = compositeroi / cerebg
        data['AV45_%s_comp/wm70' % (i+1)] = compositeroi / wm70
        # Assymetry values
        data['AV45_%s_frontal_asymmetry_negvalue_means_R<L' % (i+1)] = right_frontal - left_frontal
        data['AV45_%s_cingulate_asymmetry_negvalue_means_R<L' % (i+1)] = right_cingulate - left_cingulate
        data['AV45_%s_parietal_asymmetry_negvalue_means_R<L' % (i+1)] = right_parietal - left_parietal
        data['AV45_%s_temporal_asymmetry_negvalue_means_R<L' % (i+1)] = right_temporal - left_temporal
        data['AV45_%s_asym_summary_absvals_negvalue_means_R<L' % (i+1)] = ''
        data['AV45_%s_frontal_MR_asymmetry' % (i+1)] = ''
        data['AV45_%s_cingulate_MR_asymmetry' % (i+1)] = ''
        data['AV45_%s_parietal_MR_asymmetry' % (i+1)] = ''
        data['AV45_%s_temporal__MR_asymmetry' % (i+1)] = ''
        data['AV45_%s_ventrical_MR_asymmetry' % (i+1)] = ''

    return (new_headers, data)


if __name__ == '__main__':
    # IO files
    master_file = "../DOD_DATA.csv"
    output_file = "../DOD_DATA_synced.csv"

    # LOOKUP files
    registry_file = "../docs/DOD_REGISTRY.csv"

    # AV45 files
    av45_file = "../DOD/AV45_DOD_LONI_07.13.15_extra.csv"


    # syncing pipeline
    try:
        new_headers, new_lines = parseCSV(master_file)
    except:
        new_headers = []
        new_lines = []

    print "\nSYNCING AV45\n"
    new_headers, new_lines = syncAV45Data(new_headers, new_lines, av45_file, registry_file, dump_to=None) # adds new patients





    '''
    print "\nSYNCING DIAGNOSES\n"
    new_headers, new_lines = syncDiagnosisData(new_headers, new_lines, diagnosis_file, registry_file, arm_file, pet_meta_file, dump_to=None) # refreshes av45 dates
    print "\nSYNCING TBMSYN\n"
    new_headers, new_lines = syncTBMSynData(new_headers, new_lines, tbm_file, registry_file, dump_to=None)
    print "\nSYNCING MMSE\n"
    new_headers, new_lines = syncMMSEData(new_headers, new_lines, mmse_file, registry_file, dump_to=None)
    print "\nSYNCING ADASCOG\n"
    new_headers, new_lines = syncADASCogData(new_headers, new_lines, adni1_adas_file, adnigo2_adas_file, registry_file, dump_to=None)
    print "\nSYNCING AVLT\n"
    new_headers, new_lines = syncAVLTData(new_headers, new_lines, neuro_battery_file, registry_file, dump_to=None)
    print "\nSYNCING CSF\n"
    new_headers, new_lines = syncCSFData(new_headers, new_lines, csf_files, registry_file, dump_to=None)
    print "\nSYNCING WMH\n"
    new_headers, new_lines = syncWMHData(new_headers, new_lines, wmh_file, registry_file, dump_to=output_file)
    '''