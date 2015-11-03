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

ADNI_FIELDNAMES = ['RID','VISCODE','VISCODE2','EXAMDATE','CEREBELLUMGREYMATTER','CEREBELLUMWHITEMATTER','BRAINSTEM','WHOLECEREBELLUM',
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
ADNI_FIELDNAMES_EXTRA = ['RID','VISCODE','VISCODE2','EXAMDATE','CEREBELLUMGREYMATTER','CEREBELLUMWHITEMATTER','BRAINSTEM','WHOLECEREBELLUM',
                   'ERODED_SUBCORTICALWM','COMPOSITE','COMPOSITE_REF','FRONTAL','CINGULATE','PARIETAL','TEMPORAL',
                   'SUMMARYSUVR_WHOLECEREBNORM','SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF',
                   'LEFT-PUTAMEN','RIGHT-PUTAMEN','LEFT-CAUDATE','RIGHT-CAUDATE','LEFT-PALLIDUM','RIGHT-PALLIDUM',
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
DOD_FIELDNAMES = ['SCRNO','RID','VISCODE','EXAMDATE','CEREBELLUMGREYMATTER','BRAINSTEM','WHOLECEREBELLUM',
                  'FRONTAL','CINGULATE','PARIETAL','TEMPORAL',
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
                  'CTX_RH_MIDDLETEMPORAL_SIZE','CTX_RH_SUPERIORTEMPORAL','CTX_RH_SUPERIORTEMPORAL_SIZE']
DOD_FIELDNAMES_EXTRA = ['SCRNO','RID','VISCODE','EXAMDATE','CEREBELLUMGREYMATTER','BRAINSTEM','WHOLECEREBELLUM',
                  'FRONTAL','FRONTAL_SIZE','CINGULATE','CINGULATE_SIZE','PARIETAL','PARIETAL_SIZE','TEMPORAL','TEMPORAL_SIZE',
                  'COMPOSITE','COMPOSITE_REF','ERODED_SUBCORTICALWM',
                  'LEFT-PUTAMEN','RIGHT-PUTAMEN','LEFT-CAUDATE','RIGHT-CAUDATE','LEFT-PALLIDUM','RIGHT-PALLIDUM',
                  'LEFT-PUTAMEN_SIZE','RIGHT-PUTAMEN_SIZE','LEFT-CAUDATE_SIZE','RIGHT-CAUDATE_SIZE','LEFT-PALLIDUM_SIZE','RIGHT-PALLIDUM_SIZE',
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
                  'CTX_RH_MIDDLETEMPORAL_SIZE','CTX_RH_SUPERIORTEMPORAL','CTX_RH_SUPERIORTEMPORAL_SIZE']
ADNI_OMIT = ['LEFT_CEREBELLUM_CORTEX',
             'RIGHT_CEREBELLUM_CORTEX',
             'LEFT_UNSEGMENTEDWHITEMATTER',
             'RIGHT_UNSEGMENTEDWHITEMATTER']
ADNI_EXTRA_OMIT = ADNI_OMIT
DOD_OMIT = ADNI_OMIT
DOD_EXTRA_OMIT = ADNI_OMIT
ADNI_OMIT_SIZES = ['SCRNO', 'RID', 'PID', 'FRONTAL', 'CINGULATE', 'PARIETAL', 'BRAINSTEM', 'TEMPORAL', 'COMPOSITE', 
                   'ERODED_SUBCORTICALWM', 'CEREBELLUMGREYMATTER', 'WHOLECEREBELLUM', 'SUMMARYSUVR_WHOLECEREBNORM',
                   'SUMMARYSUVR_COMPOSITE_REFNORM', 'COMPOSITE_REF', 'VISCODE', 'VISCODE2', 'EXAMDATE', 
                   'SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF', 'SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF', 'update_stamp',
                   'CEREBELLUMWHITEMATTER', 'LEFT-PALLIDUM', 'LEFT-CAUDATE', 'LEFT-PUTAMEN', 
                   'RIGHT-PALLIDUM', 'RIGHT-CAUDATE', 'RIGHT-PUTAMEN']
ADNI_EXTRA_OMIT_SIZES = ADNI_OMIT_SIZES
DOD_OMIT_SIZES = ADNI_OMIT_SIZES
DOD_EXTRA_OMIT_SIZES = ['SCRNO', 'RID', 'PID', 'COMPOSITE', 'ERODED_SUBCORTICALWM', 'CEREBELLUMGREYMATTER',
                        'WHOLECEREBELLUM', 'SUMMARYSUVR_WHOLECEREBNORM',
                        'SUMMARYSUVR_COMPOSITE_REFNORM', 'COMPOSITE_REF', 'VISCODE', 'VISCODE2', 'EXAMDATE', 
                        'SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF', 'SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF', 'update_stamp',
                        'CEREBELLUMWHITEMATTER']



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
        l = map(float, l)
        if len(bl_header) != len(l):
            raise Exception("%s, %s" % (bl_header, l))
        bl_lines.append(l)
    return (bl_header, bl_lines)

def convertHeaderCodes(header):
    lookup = {0: 'RID',
            12: 'LEFT-PUTAMEN',
            51: 'RIGHT-PUTAMEN',
            11: 'LEFT-CAUDATE',
            50: 'RIGHT-CAUDATE',
            13: 'LEFT-PALLIDUM',
            52: 'RIGHT-PALLIDUM',
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
    if isinstance(header,list):
        converted = [str(lookup.get(int(h),h)) for h in header]
    elif isinstance(header,str) or isinstance(header,int) or isinstance(header,float):
        converted = str(lookup.get(int(header),header))
    else:
        raise Exception("Invalid header data type passed into conversion")
    return converted

def combineMeansAndSize(agg_type, mean_header, size_header, mean_values, size_values):
    if agg_type == 'adni':
        omit = ADNI_OMIT
        omit_sizes = ADNI_OMIT_SIZES
    elif agg_type == 'adni_extra':
        omit = ADNI_EXTRA_OMIT
        omit_sizes = ADNI_EXTRA_OMIT_SIZES
    elif agg_type == 'dod':
        omit = DOD_OMIT
        omit_sizes = DOD_OMIT_SIZES
    elif agg_type == 'dod_extra':
        omit = DOD_EXTRA_OMIT
        omit_sizes = DOD_EXTRA_OMIT_SIZES
    else:
        raise Exception("Bad agg_type")
    mean_values = {convertHeaderCodes(k):v for k,v in mean_values.iteritems()}
    size_values = {convertHeaderCodes(k):v for k,v in size_values.iteritems()}
    rid = int(float(mean_values['RID']))

    header_list, mean_values, size_values = additionalCalculations(mean_header, mean_values, size_values, agg_type)
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

def aggregatePreprocessingOutput(output, bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, meta_pet, registry, agg_type):
    assert agg_type in set(['adni', 'adni_extra', 'dod', 'dod_extra'])

    if agg_type == 'adni':
        pet_dates = importPetMETA(meta_pet)
        fieldnames = ADNI_FIELDNAMES
    elif agg_type == 'adni_extra':
        pet_dates = importPetMETA(meta_pet)
        fieldnames = ADNI_FIELDNAMES_EXTRA
    elif agg_type == 'dod':
        fieldnames = DOD_FIELDNAMES
    elif agg_type == 'dod_extra':
        fieldnames = DOD_FIELDNAMES_EXTRA

    num_bl = int(bl_means.split('_')[-1].replace('.csv',''))
    num_v2 = int(v2_means.split('_')[-1].replace('.csv',''))
    bl_header, bl_lines = parseCSV(bl_means)
    bl_lines = bl_lines[:num_bl]
    v2_header, v2_lines = parseCSV(v2_means)
    v2_lines = v2_lines[:num_v2]
    bl_size_header, bl_size_lines = parseCSV(bl_sizes)
    bl_size_lines = bl_size_lines[:num_bl]
    v2_size_header, v2_size_lines = parseCSV(v2_sizes)
    v2_size_lines = v2_size_lines[:num_v2]

    print "%s baseline scans" % len(bl_lines)
    print "%s visit 2 scans" % len(v2_lines)

    if agg_type in set(['adni', 'adni_extra']) and v3_means is not None and v3_sizes is not None:
        num_v3 = int(v3_means.split('_')[-1].replace('.csv',''))
        v3_header, v3_lines = parseCSV(v3_means)
        v3_lines = v3_lines[:num_v3]
        v3_size_header, v3_size_lines = parseCSV(v3_sizes)
        v3_size_lines = v3_size_lines[:num_v3]
        print "%s visit 3 scans" % len(v3_lines)
        total_iter_chain = chain(izip(repeat('BL'), zip(bl_lines, bl_size_lines)), 
                                 izip(repeat('V2'), zip(v2_lines, v2_size_lines)), 
                                 izip(repeat('V3'), zip(v3_lines, v3_size_lines)))
    else:
        total_iter_chain = chain(izip(repeat('BL'), zip(bl_lines, bl_size_lines)), 
                                 izip(repeat('V2'), zip(v2_lines, v2_size_lines)))
    
    # header list of converted here
    # headers for the row dictionary are converted in combineMeansAndSize
    mean_header = convertHeaderCodes(bl_header) # assuming headers are equivalent across files
    size_header = convertHeaderCodes(bl_size_header)

    # aggregate and write
    writer = csv.DictWriter(open(output, 'w'), fieldnames)
    writer.writeheader()
    count = 0
    for vis, (mean_line, size_line) in total_iter_chain:
        rid, all_header, all_values = combineMeansAndSize(agg_type, copy.copy(mean_header), copy.copy(size_header), mean_line, size_line)

        # add on metadata
        if agg_type in set(['adni', 'adni_extra']):
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

            subj_reg = registry.get(rid,[])
            if len(subj_reg) == 0:
                raise Exception("No possible dates for %s (%s)" % (rid, date))
            sorted_subj_reg = sorted(subj_reg, key=lambda x: abs(x['EXAMDATE']-date).days)
            metadata = sorted_subj_reg[0]

            if (date-metadata['EXAMDATE']).days > 60:
                print "%s, %s: %s, %s -> %s days" % (rid, vis, metadata, date, (date-metadata['EXAMDATE']).days)
        elif agg_type in set(['dod', 'dod_extra']):
            subj_registry = registry[rid]
            metadata = None
            if vis == 'BL':
                metadata = subj_registry[0]
            elif vis == 'V2':
                metadata = subj_registry[1]
            elif vis == 'V3':
                metadata = subj_registry[2]
            if metadata is None:
                raise Exception("Date not found: %s on %s" % (rid, vis))

        # insert viscode, viscode2, update_stamp
        data = dict(zip(all_header, all_values))
        data.update(metadata)
        for k in data.keys():
            if k not in fieldnames:
                data.pop(k, None)
        data = convertToCSVDataType(data, decimal_places=8)
        writer.writerow(data)

def additionalCalculations(headers, mean_values, size_values, agg_type):
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
    cerebWM = (float(mean_values['LEFT_UNSEGMENTEDWHITEMATTER']) + float(mean_values['RIGHT_UNSEGMENTEDWHITEMATTER']))/2.0 #- float(mean_values['CEREBELLUMGREYMATTER'])
    compref_components = [mean_values['ERODED_SUBCORTICALWM'], mean_values['BRAINSTEM'], mean_values['WHOLECEREBELLUM']]
    composite_ref = np.mean([float(_) for _ in compref_components])

    headers.append('CEREBELLUMWHITEMATTER')
    mean_values['CEREBELLUMWHITEMATTER'] = cerebWM
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

    if agg_type in set(['dod', 'dod_extra']):
        # convert RID -> SCRNO, then actual RID is added in with metadata
        headers.append('SCRNO')
        mean_values['SCRNO'] = int(mean_values['RID'])
    else:
        '''
        headers.append('PID')
        headers.remove('RID')
        mean_values['PID'] = int(mean_values['RID'])
        '''
        pass
    return (headers, mean_values, size_values)

def aggregateAllRegionFiles(bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, lut_table, pet_dates):
    # read into dataframes
    bl_means_df = pd.read_csv(bl_means)
    v2_means_df = pd.read_csv(v2_means)
    v3_means_df = pd.read_csv(v3_means)
    bl_sizes_df = pd.read_csv(bl_sizes)
    v2_sizes_df = pd.read_csv(v2_sizes)
    v3_sizes_df = pd.read_csv(v3_sizes)
    # convert headers
    bl_means_df.columns = ['RID' if int(_)==0 else lut_table[int(_)].replace('-','_').upper() for _ in bl_means_df.columns]
    v2_means_df.columns = ['RID' if int(_)==0 else lut_table[int(_)].replace('-','_').upper() for _ in v2_means_df.columns]
    v3_means_df.columns = ['RID' if int(_)==0 else lut_table[int(_)].replace('-','_').upper() for _ in v3_means_df.columns]
    bl_sizes_df.columns = ['RID' if int(_)==0 else lut_table[int(_)].replace('-','_').upper() for _ in bl_sizes_df.columns]
    v2_sizes_df.columns = ['RID' if int(_)==0 else lut_table[int(_)].replace('-','_').upper() for _ in v2_sizes_df.columns]
    v3_sizes_df.columns = ['RID' if int(_)==0 else lut_table[int(_)].replace('-','_').upper() for _ in v3_sizes_df.columns]
    # add exam dates
    bl_means_df['EXAMDATE'] = ''
    v2_means_df['EXAMDATE'] = ''
    v3_means_df['EXAMDATE'] = ''

    # fill in exam dates
    for i in bl_means_df.index:
        bl_means_df.ix[i,'EXAMDATE'] = pet_dates[bl_means_df.ix[i,'RID']][0].strftime('%m/%d/%Y')
    for i in v2_means_df.index:
        v2_means_df.ix[i,'EXAMDATE'] = pet_dates[v2_means_df.ix[i,'RID']][1].strftime('%m/%d/%Y')
    for i in v3_means_df.index:
        v3_means_df.ix[i,'EXAMDATE'] = pet_dates[v3_means_df.ix[i,'RID']][2].strftime('%m/%d/%Y')

    # merge sizes
    bl_means_df = bl_means_df.merge(bl_sizes_df, on='RID', suffixes=['','_SIZE'])
    v2_means_df = v2_means_df.merge(v2_sizes_df, on='RID', suffixes=['','_SIZE'])
    v3_means_df = v3_means_df.merge(v3_sizes_df, on='RID', suffixes=['','_SIZE'])

    # concat visits
    all_rows = pd.concat((bl_means_df, v2_means_df, v3_means_df), axis=0)
    all_rows.set_index(['RID', 'EXAMDATE'], inplace=True)
    all_rows = all_rows.sort_index()

    print all_rows.columns
    print all_rows.shape

    return all_rows



def findPreprocessOutputFiles(folder_name, nontp=False, allregions=False):
    bl_means = v2_means = v3_means = bl_sizes = v2_sizes = v3_sizes = None
    addon = "_nontp" if nontp else ""
    for filename in os.listdir(folder_name):
        if allregions and 'allregions' not in filename:
            continue
        if "BL%s_means" % addon in filename:
            bl_means = os.path.join(folder_name, filename)
        elif "V2%s_means" % addon in filename:
            v2_means = os.path.join(folder_name, filename)
        elif "V3%s_means" % addon in filename:
            v3_means = os.path.join(folder_name, filename)
        elif "BL%s_roisize" % addon in filename:
            bl_sizes = os.path.join(folder_name, filename)
        elif "V2%s_roisize" % addon in filename:
            v2_sizes = os.path.join(folder_name, filename)
        elif "V3%s_roisize" % addon in filename:
            v3_sizes = os.path.join(folder_name, filename)
    return (bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes)

if __name__ == "__main__":
    lut_file = "../FreeSurferColorLUT.txt"
    meta_pet = "../docs/ADNI/PET_META_LIST.csv"
    lut_table = importFreesurferLookup(lut_file)

    # for adni av45 nontp all regions
    output = '../output/UCBERKELEYAV45_11_02_15_allregions_nontp.csv'
    preprocess_folder =  '../docs/AV45_allregions_preprocess_output_11_02_15'
    bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes = findPreprocessOutputFiles(preprocess_folder, nontp=True, allregions=True)
    pet_dates = importPetMETA(meta_pet)
    df = aggregateAllRegionFiles(bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, lut_table, pet_dates)
    df.to_csv(output)

    # # for adni av45 all regions
    # output = '../output/UCBERKELEYAV45_11_02_15_allregions.csv'
    # preprocess_folder =  '../docs/AV45_preprocess_output_11_02_15'
    # bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes = findPreprocessOutputFiles(preprocess_folder, nontp=False, allregions=True)
    # pet_dates = importPetMETA(meta_pet)
    # df = aggregateAllRegionFiles(bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, lut_table, pet_dates)
    # df.to_csv(output)

    # # for adni av45 nontp
    # output = '../output/UCBERKELEYAV45_09_25_15_extra_nontp.csv'
    # preprocess_folder =  '../docs/AV45_preprocess_output_09_25_15'
    # registry = importRegistry("../docs/ADNI/REGISTRY.csv") 
    # agg_type = 'adni_extra'
    # bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes = findPreprocessOutputFiles(preprocess_folder, nontp=True)
    # aggregatePreprocessingOutput(output, bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, 
    #                              meta_pet, registry, agg_type)
    
    # # for adni av45
    # output = '../output/UCBERKELEYAV45_09_25_15_extra.csv'
    # preprocess_folder =  '../docs/AV45_preprocess_output_09_25_15'
    # registry = importRegistry("../docs/ADNI/REGISTRY.csv") 
    # agg_type = 'adni_extra'
    # bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes = findPreprocessOutputFiles(preprocess_folder, nontp=False)
    # aggregatePreprocessingOutput(output, bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, 
    #                              meta_pet, registry, agg_type)
    

    # # for adni dod (add in adni controls)
    # output = '../output/AV45_DOD_LONI_10.22.15_extra_withcontrols.csv'
    # temp_output = '../output/temp.csv'
    # preprocess_folder =  '../docs/AV45_DOD_preprocess_output_10_22_15'
    # adni_preprocess_folder = '../docs/AV45_preprocess_output_09_25_15'
    # registry = importDODRegistry("../docs/DOD/DOD_REGISTRY.csv")
    # registry_adni = importRegistry("../docs/ADNI/REGISTRY.csv")
    # meta_pet = None
    # meta_pet_adni = None
    # agg_type = 'dod_extra'
    # bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes = findPreprocessOutputFiles(preprocess_folder)
    # bl_means_adni, v2_means_adni, v3_means_adni, bl_sizes_adni, v2_sizes_adni, v3_sizes_adni = findPreprocessOutputFiles(adni_preprocess_folder, nontp=True)
    # aggregatePreprocessingOutput(output, bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, 
    #                              meta_pet, registry, agg_type)
    # aggregatePreprocessingOutput(temp_output, bl_means_adni, v2_means_adni, v3_means_adni, bl_sizes_adni, v2_sizes_adni, v3_sizes_adni, 
    #                              meta_pet_adni, registry_adni, agg_type)
    # print 'appending'
    # appendCSV(output, temp_output)
    # print 'removing'
    # removeFile(temp_output)


    # # # for adni dod (don't add in adni controls, for uploading)
    # output = '../output/AV45_DOD_LONI_10.22.15_extra.csv'
    # preprocess_folder =  '../docs/AV45_DOD_preprocess_output_10_22_15'
    # registry = importDODRegistry("../docs/DOD/DOD_REGISTRY.csv")
    # meta_pet = None
    # meta_pet_adni = None
    # agg_type = 'dod'
    # bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes = findPreprocessOutputFiles(preprocess_folder)
    # aggregatePreprocessingOutput(output, bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, 
    #                              meta_pet, registry, agg_type)

