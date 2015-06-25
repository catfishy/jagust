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

import csv
from itertools import chain, izip, repeat
import numpy as np
from collections import defaultdict
from datetime import datetime, timedelta
import copy
import codecs

def readHeaderAndLines(csv_file, limit=None):
    bl_lines = []
    bl_header = None
    reader = csv.reader(open(csv_file, 'rU'))
    for i, l in enumerate(reader):
        if i == 0:
            bl_header = l
            continue
        elif limit is not None and i > limit:
            break
        if len(bl_header) != len(l):
            raise Exception("%s, %s" % (bl_header, l))
        bl_lines.append(l)
    return (bl_header, bl_lines)

def convertHeaderCodes(header):
    lookup = {0: 'RID',
              1025: 'CTX_LH_PRECUNEUS',
              1026: 'CTX_LH_ROSTRALANTERIORCINGULATE',
              1027: 'CTX_LH_ROSTRALMIDDLEFRONTAL',
              1028: 'CTX_LH_SUPERIORFRONTAL',
              1029: 'CTX_LH_SUPERIORPARIETAL',
              1030: 'CTX_LH_SUPERIORTEMPORAL',
              1031: 'CTX_LH_SUPRAMARGINAL',
              8: 'LEFT_CEREBELLUM_CORTEX',
              16: 'BRAINSTEM',
              47: 'RIGHT_CEREBELLUM_CORTEX',
              1032: 'CTX_LH_FRONTALPOLE',
              2026: 'CTX_RH_ROSTRALANTERIORCINGULATE',
              1003: 'CTX_LH_CAUDALMIDDLEFRONTAL',
              5000: 'WHOLECEREBELLUM',
              5001: 'LEFT_UNSEGMENTEDWHITEMATTER',
              5002: 'RIGHT_UNSEGMENTEDWHITEMATTER',
              5003: 'CEREBELLUMGREYMATTER',
              4000: 'ERODED_SUBCORTICALWM',
              2032: 'CTX_RH_FRONTALPOLE',
              3000: 'FRONTAL',
              3001: 'CINGULATE',
              3002: 'PARIETAL',
              3003: 'TEMPORAL',
              3004: 'COMPOSITE',
              2002: 'CTX_RH_CAUDALANTERIORCINGULATE',
              2003: 'CTX_RH_CAUDALMIDDLEFRONTAL',
              2008: 'CTX_RH_INFERIORPARIETAL',
              2010: 'CTX_RH_ISTHMUSCINGULATE',
              2012: 'CTX_RH_LATERALORBITOFRONTAL',
              2014: 'CTX_RH_MEDIALORBITOFRONTAL',
              2015: 'CTX_RH_MIDDLETEMPORAL',
              2018: 'CTX_RH_PARSOPERCULARIS',
              2019: 'CTX_RH_PARSORBITALIS',
              2020: 'CTX_RH_PARSTRIANGULARIS',
              2023: 'CTX_RH_POSTERIORCINGULATE',
              2025: 'CTX_RH_PRECUNEUS',
              1002: 'CTX_LH_CAUDALANTERIORCINGULATE',
              2027: 'CTX_RH_ROSTRALMIDDLEFRONTAL',
              2028: 'CTX_RH_SUPERIORFRONTAL',
              2029: 'CTX_RH_SUPERIORPARIETAL',
              2030: 'CTX_RH_SUPERIORTEMPORAL',
              2031: 'CTX_RH_SUPRAMARGINAL',
              1008: 'CTX_LH_INFERIORPARIETAL',
              1010: 'CTX_LH_ISTHMUSCINGULATE',
              1012: 'CTX_LH_LATERALORBITOFRONTAL',
              1014: 'CTX_LH_MEDIALORBITOFRONTAL',
              1015: 'CTX_LH_MIDDLETEMPORAL',
              1018: 'CTX_LH_PARSOPERCULARIS',
              1019: 'CTX_LH_PARSORBITALIS',
              1020: 'CTX_LH_PARSTRIANGULARIS',
              1023: 'CTX_LH_POSTERIORCINGULATE'}
    converted = [str(lookup.get(int(h),h)) for h in header]
    return converted

def combineMeansAndSize(mean_header, size_header, mean_row, size_row):
    omit = ['LEFT_CEREBELLUM_CORTEX',
            'RIGHT_CEREBELLUM_CORTEX',
            'LEFT_UNSEGMENTEDWHITEMATTER',
            'RIGHT_UNSEGMENTEDWHITEMATTER']
    omit_sizes = ['RID', 'FRONTAL', 'CINGULATE', 'PARIETAL', 'BRAINSTEM', 'TEMPORAL', 'COMPOSITE', 
                  'ERODED_SUBCORTICALWM', 'CEREBELLUMGREYMATTER', 'WHOLECEREBELLUM', 'SUMMARYSUVR_WHOLECEREBNORM',
                  'SUMMARYSUVR_COMPOSITE_REFNORM', 'COMPOSITE_REF', 'VISCODE', 'VISCODE2', 'EXAMDATE', 
                  'SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF', 'SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF', 'update_stamp']
    mean_values = dict(zip(mean_header, mean_row))
    size_values = dict(zip(size_header, size_row))
    rid = int(mean_values['RID'])

    header_list, mean_values, size_values = additionalCalculations(mean_header, mean_values, size_values)
    all_headers = []
    all_values = []
    for h in header_list:
        if h in omit:
            continue
        all_headers.append(str(h))
        all_values.append(mean_values[h])
        if h not in omit_sizes:
            all_headers.append(str(h) + '_SIZE')
            all_values.append(size_values[h])
    return (rid, all_headers, all_values)


def arrangePETRegistry(headers, lines):
    registry = defaultdict(dict)
    for l in lines:
        data = dict(zip(headers,l))
        if data['EXAMDATE'] == '' or data['VISCODE'] in set(['sc', 'scmri']) or data['VISCODE2'] in set(['sc', 'scmri']):
            continue
        subj = int(data['RID'])
        date = datetime.strptime(data['EXAMDATE'],'%Y-%m-%d')
        registry[subj][date] = {'VISCODE': data['VISCODE'],
                                'VISCODE2': data['VISCODE2'],
                                'update_stamp': data['update_stamp']}
    return registry

def arrangePETDates(headers, lines):
    dates = defaultdict(list)
    for l in lines:
        data = dict(zip(headers,l))
        subj = int(data['Subject'].split("_")[-1])
        date = datetime.strptime(data['Scan Date'],'%m/%d/%y')
        dates[subj].append(date)
    for k in dates.keys():
        dates[k] = sorted(dates[k])
        if len(dates[k]) > 3:
            print "%s, %s" % (k, len(dates[k]))
    return dates

def aggregatePreprocessingOutput(total_output, bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, meta_pet, registry):
    #et_meta_header, pet_meta_lines = readHeaderAndLines(meta_pet)
    registry = arrangePETRegistry(*readHeaderAndLines(registry))
    pet_dates = arrangePETDates(*readHeaderAndLines(meta_pet))

    num_bl = int(bl_means.split('_')[-1].replace('.csv',''))
    num_v2 = int(v2_means.split('_')[-1].replace('.csv',''))
    num_v3 = int(v3_means.split('_')[-1].replace('.csv',''))
    bl_header, bl_lines = readHeaderAndLines(bl_means, limit=num_bl)
    v2_header, v2_lines = readHeaderAndLines(v2_means, limit=num_v2)
    v3_header, v3_lines = readHeaderAndLines(v3_means, limit=num_v3)
    bl_size_header, bl_size_lines = readHeaderAndLines(bl_sizes, limit=num_bl)
    v2_size_header, v2_size_lines = readHeaderAndLines(v2_sizes, limit=num_v2)
    v3_size_header, v3_size_lines = readHeaderAndLines(v3_sizes, limit=num_v3)
    print "%s baseline scans" % len(bl_lines)
    print "%s visit 2 scans" % len(v2_lines)
    print "%s visit 3 scans" % len(v3_lines)
    mean_header = convertHeaderCodes(bl_header) # assuming headers are equivalent across files
    size_header = convertHeaderCodes(bl_size_header)
    # aggregate and write
    fieldnames = ['RID','VISCODE','VISCODE2','EXAMDATE','CEREBELLUMGREYMATTER','BRAINSTEM','WHOLECEREBELLUM',
                  'ERODED_SUBCORTICALWM','COMPOSITE_REF','FRONTAL','CINGULATE','PARIETAL','TEMPORAL',
                  'SUMMARYSUVR_WHOLECEREBNORM','SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF',
                  'SUMMARYSUVR_COMPOSITE_REFNORM','SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF','CTX_LH_CAUDALMIDDLEFRONTAL',
                  'CTX_LH_CAUDALMIDDLEFRONTAL_SIZE','CTX_LH_LATERALORBITOFRONTAL','CTX_LH_LATERALORBITOFRONTAL_SIZE',
                  'CTX_LH_MEDIALORBITOFRONTAL','CTX_LH_MEDIALORBITOFRONTAL_SIZE','CTX_LH_PARSOPERCULARIS',
                  'CTX_LH_PARSOPERCULARIS_SIZE','CTX_LH_PARSORBITALIS','CTX_LH_PARSORBITALIS_SIZE','CTX_LH_PARSTRIANGULARIS',
                  'CTX_LH_PARSTRIANGULARIS_SIZE','CTX_LH_ROSTRALMIDDLEFRONTAL','CTX_LH_ROSTRALMIDDLEFRONTAL_SIZE',
                  'CTX_LH_SUPERIORFRONTAL','CTX_LH_SUPERIORFRONTAL_SIZE','CTX_LH_FRONTALPOLE','CTX_LH_FRONTALPOLE_SIZE',
                  'CTX_RH_CAUDALMIDDLEFRONTAL','CTX_RH_CAUDALMIDDLEFRONTAL_SIZE','CTX_RH_LATERALORBITOFRONTAL',
                  'CTX_RH_LATERALORBITOFRONTAL_SIZE','CTX_RH_MEDIALORBITOFRONTAL','CTX_RH_MEDIALORBITOFRONTAL_SIZE',
                  'CTX_RH_PARSOPERCULARIS','CTX_RH_PARSOPERCULARIS_SIZE','CTX_RH_PARSORBITALIS','CTX_RH_PARSORBITALIS_SIZE',
                  'CTX_RH_PARSTRIANGULARIS','CTX_RH_PARSTRIANGULARIS_SIZE','CTX_RH_ROSTRALMIDDLEFRONTAL',
                  'CTX_RH_ROSTRALMIDDLEFRONTAL_SIZE','CTX_RH_SUPERIORFRONTAL','CTX_RH_SUPERIORFRONTAL_SIZE',
                  'CTX_RH_FRONTALPOLE','CTX_RH_FRONTALPOLE_SIZE','CTX_LH_CAUDALANTERIORCINGULATE','CTX_LH_CAUDALANTERIORCINGULATE_SIZE',
                  'CTX_LH_ISTHMUSCINGULATE','CTX_LH_ISTHMUSCINGULATE_SIZE','CTX_LH_POSTERIORCINGULATE','CTX_LH_POSTERIORCINGULATE_SIZE',
                  'CTX_LH_ROSTRALANTERIORCINGULATE','CTX_LH_ROSTRALANTERIORCINGULATE_SIZE','CTX_RH_CAUDALANTERIORCINGULATE',
                  'CTX_RH_CAUDALANTERIORCINGULATE_SIZE','CTX_RH_ISTHMUSCINGULATE','CTX_RH_ISTHMUSCINGULATE_SIZE','CTX_RH_POSTERIORCINGULATE',
                  'CTX_RH_POSTERIORCINGULATE_SIZE','CTX_RH_ROSTRALANTERIORCINGULATE','CTX_RH_ROSTRALANTERIORCINGULATE_SIZE',
                  'CTX_LH_INFERIORPARIETAL','CTX_LH_INFERIORPARIETAL_SIZE','CTX_LH_PRECUNEUS','CTX_LH_PRECUNEUS_SIZE','CTX_LH_SUPERIORPARIETAL',
                  'CTX_LH_SUPERIORPARIETAL_SIZE','CTX_LH_SUPRAMARGINAL','CTX_LH_SUPRAMARGINAL_SIZE','CTX_RH_INFERIORPARIETAL',
                  'CTX_RH_INFERIORPARIETAL_SIZE','CTX_RH_PRECUNEUS','CTX_RH_PRECUNEUS_SIZE','CTX_RH_SUPERIORPARIETAL',
                  'CTX_RH_SUPERIORPARIETAL_SIZE','CTX_RH_SUPRAMARGINAL','CTX_RH_SUPRAMARGINAL_SIZE','CTX_LH_MIDDLETEMPORAL',
                  'CTX_LH_MIDDLETEMPORAL_SIZE','CTX_LH_SUPERIORTEMPORAL','CTX_LH_SUPERIORTEMPORAL_SIZE','CTX_RH_MIDDLETEMPORAL',
                  'CTX_RH_MIDDLETEMPORAL_SIZE','CTX_RH_SUPERIORTEMPORAL','CTX_RH_SUPERIORTEMPORAL_SIZE','update_stamp']
    writer = csv.DictWriter(open(output, 'w'), fieldnames)
    writer.writeheader()
    count = 0
    for vis, (mean_line, size_line) in chain(izip(repeat('BL'), zip(bl_lines, bl_size_lines)), 
                                                 izip(repeat('V2'), zip(v2_lines, v2_size_lines)), 
                                                 izip(repeat('V3'), zip(v3_lines, v3_size_lines))):
        rid, all_header, all_values = combineMeansAndSize(copy.copy(mean_header), copy.copy(size_header), mean_line, size_line)

        # add on metadata
        subj_meta_list = pet_dates[rid]
        date = None
        if vis == 'BL':
            date = subj_meta_list[0]
        elif vis == 'V2':
            date = subj_meta_list[1]
        elif vis == 'V3':
            date = subj_meta_list[2]
        if date is None:
            raise Exception("Date not found: %s on %s" % (rid, vis))

        registry_list = registry[rid]
        cap_date = date+timedelta(days=30)
        possible_dates = [_ for _ in registry_list.keys() if _ <= cap_date]
        if len(possible_dates) == 0:
            raise Exception("No possible dates for %s (%s): %s" % (rid, date, registry_list.keys()))
        registry_date = sorted(possible_dates, key=lambda x: cap_date-x)[0]
        metadata = registry_list[registry_date]
        metadata['EXAMDATE'] = date.strftime('%Y-%m-%d')

        if (date-registry_date).days > 30:
            print "%s, %s: %s, %s -> %s (%s days)" % (rid, vis, metadata, date, registry_date, (date-registry_date).days)

        # insert viscode, viscode2, update_stamp
        data = dict(zip(all_header, all_values))
        data.update(metadata)
        for k in data.keys():
            if k not in fieldnames:
                data.pop(k,None)
        writer.writerow(data)

def additionalCalculations(headers, mean_values, size_values):
    '''
        SUMMARYSUVR_WHOLECEREBNORM
        SUMMARYSUVR_COMPOSITE_REFNORM
        COMPOSITE_REF
        SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF
        SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF



        VISCODE
        VISCODE2
        EXAMDATE

        update_stamp
    '''
    composite = float(mean_values['COMPOSITE'])
    wholecereb = float(mean_values['WHOLECEREBELLUM'])
    compref_components = [mean_values['ERODED_SUBCORTICALWM'], mean_values['BRAINSTEM'], mean_values['WHOLECEREBELLUM']]
    composite_ref = np.mean([float(_) for _ in compref_components])

    headers.append('SUMMARYSUVR_WHOLECEREBNORM')
    mean_values['SUMMARYSUVR_WHOLECEREBNORM'] = composite / wholecereb
    headers.append('SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF')
    mean_values['SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF'] = 1 if mean_values['SUMMARYSUVR_WHOLECEREBNORM'] >= 1.11 else 0
    headers.append('SUMMARYSUVR_COMPOSITE_REFNORM')
    mean_values['SUMMARYSUVR_COMPOSITE_REFNORM'] = composite / composite_ref
    headers.append('SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF')
    mean_values['SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF'] = 1 if mean_values['SUMMARYSUVR_COMPOSITE_REFNORM'] >= 0.79 else 0

    headers.append('COMPOSITE_REF')
    mean_values['COMPOSITE_REF'] = composite_ref

    return (headers, mean_values, size_values)


if __name__ == "__main__":
    output = '../UCBERKELEYAV45_06_22_15.csv'
    registry = "../docs/registry_clean.csv"
    meta_pet = "../docs/PET_META_LIST_edited.csv"
    bl_means = "../docs/AV45_preprocess_output_06_22_15/AV45_BL_means_24-Jun-2015_1089.csv"
    v2_means = "../docs/AV45_preprocess_output_06_22_15/AV45_V2_means_24-Jun-2015_619.csv"
    v3_means = "../docs/AV45_preprocess_output_06_22_15/AV45_V3_means_24-Jun-2015_92.csv"
    bl_sizes = "../docs/AV45_preprocess_output_06_22_15/AV45_BL_roisize_24-Jun-2015_1089.csv"
    v2_sizes = "../docs/AV45_preprocess_output_06_22_15/AV45_V2_roisize_24-Jun-2015_619.csv"
    v3_sizes = "../docs/AV45_preprocess_output_06_22_15/AV45_V3_roisize_24-Jun-2015_92.csv"
    aggregatePreprocessingOutput(output, bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, meta_pet, registry)


'''
## output keys ###

RID
VISCODE
VISCODE2
EXAMDATE


CEREBELLUMGREYMATTER
BRAINSTEM
WHOLECEREBELLUM
ERODED_SUBCORTICALWM
COMPOSITE_REF
FRONTAL
CINGULATE
PARIETAL
TEMPORAL
SUMMARYSUVR_WHOLECEREBNORM
SUMMARYSUVR_COMPOSITE_REFNORM


CTX_LH_CAUDALMIDDLEFRONTAL
CTX_LH_LATERALORBITOFRONTAL
CTX_LH_MEDIALORBITOFRONTAL
CTX_LH_PARSOPERCULARIS
CTX_LH_PARSORBITALIS
CTX_LH_PARSTRIANGULARIS
CTX_LH_ROSTRALMIDDLEFRONTAL
CTX_LH_SUPERIORFRONTAL
CTX_LH_FRONTALPOLE
CTX_RH_CAUDALMIDDLEFRONTAL
CTX_RH_LATERALORBITOFRONTAL
CTX_RH_MEDIALORBITOFRONTAL
CTX_RH_PARSOPERCULARIS
CTX_RH_PARSORBITALIS
CTX_RH_PARSTRIANGULARIS
CTX_RH_ROSTRALMIDDLEFRONTAL
CTX_RH_SUPERIORFRONTAL
CTX_RH_FRONTALPOLE
CTX_LH_CAUDALANTERIORCINGULATE
CTX_LH_ISTHMUSCINGULATE
CTX_LH_POSTERIORCINGULATE
CTX_LH_ROSTRALANTERIORCINGULATE
CTX_RH_CAUDALANTERIORCINGULATE
CTX_RH_ISTHMUSCINGULATE
CTX_RH_POSTERIORCINGULATE
CTX_RH_ROSTRALANTERIORCINGULATE
CTX_LH_INFERIORPARIETALE
CTX_LH_PRECUNEUS
CTX_LH_SUPERIORPARIETAL
CTX_LH_SUPRAMARGINAL
CTX_RH_INFERIORPARIETAL
CTX_RH_PRECUNEUS
CTX_RH_SUPERIORPARIETAL
CTX_RH_SUPRAMARGINAL
CTX_LH_MIDDLETEMPORAL
CTX_LH_SUPERIORTEMPORAL
CTX_RH_MIDDLETEMPORAL
CTX_RH_SUPERIORTEMPORAL
update_stamp
'''