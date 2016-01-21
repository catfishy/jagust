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

ADNI_FIELDNAMES = ['RID','VISCODE','VISCODE2','EXAMDATE','CEREBELLUMGREYMATTER','BRAIN_STEM','WHOLECEREBELLUM',
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
ADNI_FIELDNAMES_EXTRA = ['RID','VISCODE','VISCODE2','EXAMDATE','CEREBELLUMGREYMATTER','CEREBELLUMWHITEMATTER','BRAIN_STEM','WHOLECEREBELLUM',
                   'ERODED_SUBCORTICALWM','COMPOSITE','COMPOSITE_REF','FRONTAL','CINGULATE','PARIETAL','TEMPORAL',
                   'SUMMARYSUVR_WHOLECEREBNORM','SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF',
                   'LEFT_PUTAMEN','RIGHT_PUTAMEN','LEFT_CAUDATE','RIGHT_CAUDATE','LEFT_PALLIDUM','RIGHT_PALLIDUM',
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
DOD_FIELDNAMES = ['SCRNO','VISCODE','EXAMDATE','CEREBELLUMGREYMATTER','BRAIN_STEM','WHOLECEREBELLUM',
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
DOD_FIELDNAMES_EXTRA = ['SCRNO','VISCODE','EXAMDATE','CEREBELLUMGREYMATTER','BRAIN_STEM','WHOLECEREBELLUM',
                  'FRONTAL','FRONTAL_SIZE','CINGULATE','CINGULATE_SIZE','PARIETAL','PARIETAL_SIZE','TEMPORAL','TEMPORAL_SIZE',
                  'COMPOSITE','COMPOSITE_REF','ERODED_SUBCORTICALWM',
                  'LEFT_PUTAMEN','RIGHT_PUTAMEN','LEFT_CAUDATE','RIGHT_CAUDATE','LEFT_PALLIDUM','RIGHT_PALLIDUM',
                  'LEFT_PUTAMEN_SIZE','RIGHT_PUTAMEN_SIZE','LEFT_CAUDATE_SIZE','RIGHT_CAUDATE_SIZE','LEFT_PALLIDUM_SIZE','RIGHT_PALLIDUM_SIZE',
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
TAU_FIELDNAMES = ['RID','VISCODE','VISCODE2','EXAMDATE',
                  'BRAAK1', 'BRAAK1_SIZE', 'BRAAK2', 'BRAAK2_SIZE', 'BRAAK3', 'BRAAK3_SIZE',
                  'BRAAK4', 'BRAAK4_SIZE', 'BRAAK5', 'BRAAK5_SIZE', 'BRAAK6', 'BRAAK6_SIZE',
                  'CEREBELLUMGREYMATTER', 'CEREBELLUMGREYMATTER_SIZE']
DOD_TAU_FIELDNAMES = ['SCRNO','VISCODE','EXAMDATE',
                      'BRAAK1', 'BRAAK1_SIZE', 'BRAAK2', 'BRAAK2_SIZE', 'BRAAK3', 'BRAAK3_SIZE',
                      'BRAAK4', 'BRAAK4_SIZE', 'BRAAK5', 'BRAAK5_SIZE', 'BRAAK6', 'BRAAK6_SIZE',
                      'CEREBELLUMGREYMATTER', 'CEREBELLUMGREYMATTER_SIZE']



'''
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

            if (date-metadata['EXAMDATE']).days > 120:
                # don't accept the visit codes
                print "%s, %s: %s, %s -> %s days" % (rid, vis, metadata, date, (date-metadata['EXAMDATE']).days)
                metadata['VISCODE'] = ''
                metadata['VISCODE2'] = ''
                
            # overwrite with actual date
            metadata['EXAMDATE'] = date

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
        pass
    return (headers, mean_values, size_values)

'''

def DFWeightedMean(df, keys):
    size_keys = ['%s_SIZE' % _ for _ in keys]
    sizes = df.loc[:,size_keys].sum(axis=1)
    ratios = df.loc[:,size_keys].divide(sizes, axis=0)
    ratios.columns = [_.replace('_SIZE','') for _ in ratios.columns]
    means = df.loc[:,keys].multiply(ratios).sum(axis=1)
    return (means, sizes)

def additionalTauCalculations(df, lut_table, keys=None):
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

def aggregateAllRegionFiles(bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, lut_table, pet_dates, registry, dod=False):
    '''
    if dod:
        - RID -> SCRNO
        - VISCODE2 is removed
    '''
    # read into dataframes
    try:
        bl_means_df = pd.read_csv(bl_means)
        bl_sizes_df = pd.read_csv(bl_sizes)
    except Exception as e:
        print "Can't Read BL input: %s" % e
        bl_means_df = pd.DataFrame(columns=['RID', 'EXAMDATE'])
        bl_sizes_df = pd.DataFrame(columns=['RID', 'EXAMDATE'])
    try:
        v2_means_df = pd.read_csv(v2_means)
        v2_sizes_df = pd.read_csv(v2_sizes)
    except Exception as e:
        print "Can't Read V2 input: %s" % e
        v2_means_df = pd.DataFrame(columns=['RID', 'EXAMDATE'])
        v2_sizes_df = pd.DataFrame(columns=['RID', 'EXAMDATE'])
    try:
        v3_means_df = pd.read_csv(v3_means)
        v3_sizes_df = pd.read_csv(v3_sizes)
    except Exception as e:
        print "Can't Read V3 input: %s" % e
        v3_means_df = pd.DataFrame(columns=['RID', 'EXAMDATE'])
        v3_sizes_df = pd.DataFrame(columns=['RID', 'EXAMDATE'])

    # convert headers
    bl_means_df.columns = [translateColumn(_, lut_table) for _ in bl_means_df.columns]
    v2_means_df.columns = [translateColumn(_, lut_table) for _ in v2_means_df.columns]
    v3_means_df.columns = [translateColumn(_, lut_table) for _ in v3_means_df.columns]
    bl_sizes_df.columns = [translateColumn(_, lut_table) for _ in bl_sizes_df.columns]
    v2_sizes_df.columns = [translateColumn(_, lut_table) for _ in v2_sizes_df.columns]
    v3_sizes_df.columns = [translateColumn(_, lut_table) for _ in v3_sizes_df.columns]
    # add exam dates
    bl_means_df['EXAMDATE'] = ''
    v2_means_df['EXAMDATE'] = ''
    v3_means_df['EXAMDATE'] = ''

    # fill in exam dates + visit codes
    try:
        for i in bl_means_df.index:
            rid = bl_means_df.ix[i,'RID']
            date = pet_dates[rid][0]
            vc, vc2 = getVisitCode(rid, date, registry, cutoff=90)
            bl_means_df.ix[i,'EXAMDATE'] = date
            bl_means_df.ix[i,'VISCODE'] = vc
            bl_means_df.ix[i,'VISCODE2'] = vc2
        for i in v2_means_df.index:
            rid = v2_means_df.ix[i,'RID']
            date = pet_dates[rid][1]
            vc, vc2 = getVisitCode(rid, date, registry, cutoff=90)
            v2_means_df.ix[i,'EXAMDATE'] = date
            v2_means_df.ix[i,'VISCODE'] = vc
            v2_means_df.ix[i,'VISCODE2'] = vc2
        for i in v3_means_df.index:
            rid = v3_means_df.ix[i,'RID']
            date = pet_dates[rid][2]
            vc, vc2 = getVisitCode(rid, date, registry, cutoff=90)
            v3_means_df.ix[i,'EXAMDATE'] = date
            v3_means_df.ix[i,'VISCODE'] = vc
            v3_means_df.ix[i,'VISCODE2'] = vc2
    except Exception as e:
        print "PET DATE Problem: %s" % rid
        raise e

    # merge sizes
    bl_means_df = bl_means_df.merge(bl_sizes_df, on='RID', suffixes=['','_SIZE'])
    v2_means_df = v2_means_df.merge(v2_sizes_df, on='RID', suffixes=['','_SIZE'])
    v3_means_df = v3_means_df.merge(v3_sizes_df, on='RID', suffixes=['','_SIZE'])

    # concat visits
    all_rows = pd.concat((bl_means_df, v2_means_df, v3_means_df), axis=0)
    all_rows.reset_index(inplace=True, drop=True)
    all_rows.sort(['RID','EXAMDATE'], inplace=True)
    all_rows.loc[:,'EXAMDATE'] = all_rows.loc[:,'EXAMDATE'].apply(lambda x: x.strftime('%Y-%m-%d'))

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
    bl_means = v2_means = v3_means = bl_sizes = v2_sizes = v3_sizes = None
    addon = "_nontp" if nontp else ""
    for filename in os.listdir(folder_name):
        '''
        if allregions and 'allregions' not in filename:
            continue
        '''
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
    # freesurfer region lookup
    lut_file = "../FreeSurferColorLUT.txt"
    lut_table = importFreesurferLookup(lut_file)

    registry_file = "../docs/ADNI/REGISTRY.csv"
    dod_registry_file = "../docs/DOD/DOD_REGISTRY.csv"
    meta_pet = "../docs/ADNI/AV45META.csv"
    meta_tau = '../docs/ADNI/TAUMETA.csv'
    dod_meta_pet = "../docs/DOD/AV45META.csv"
    dod_meta_tau = "../docs/DOD/TAUMETA.csv"

    # registry imports
    adni_registry = importRegistry(meta_pet) 
    dod_registry = importDODRegistry(dod_registry_file)
    
    # pet date imports
    adni_av45_pet_dates = importScanMeta(meta_pet)
    adni_tau_pet_dates = importScanMeta(meta_tau)
    dod_av45_pet_dates = importScanMeta(dod_meta_pet)
    dod_tau_pet_dates = importScanMeta(dod_meta_tau)
    
    #timestamp = datetime.now().strftime('%m_%d_%y')
    timestamp = '01_03_16'

    # preprocess output folders
    adni_av45_preprocess_folder = '../docs/AV45_preprocess_output_01_03_16'
    dod_av45_preprocess_folder = '../docs/AV45_DOD_preprocess_output_01_03_16'
    adni_tau_preprocess_folder = '../docs/TAU_preprocess_output_01_03_16'
    dod_tau_preprocess_folder = '../docs/TAU_DOD_preprocess_output_01_03_16'

    # create output folder
    output_folder = '../output/%s' % timestamp
    mkdir_p(output_folder)

    # FOR LONI UPLOAD (ADNI AV45)
    regular_output = os.path.join(output_folder, 'LONI_UCBERKELEYAV45_%s_regular_nontp.csv' % (timestamp))
    allregions_output = os.path.join(output_folder, 'LONI_UCBERKELEYAV45_%s_allregions_nontp.csv' % (timestamp))
    merged_output = os.path.join(output_folder, 'LONI_UCBERKELEYAV45_%s_merged_nontp.csv' % (timestamp))
    bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes = findPreprocessOutputFiles(adni_av45_preprocess_folder, nontp=True)
    df = aggregateAllRegionFiles(bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, lut_table, adni_av45_pet_dates, adni_registry)
    df.to_csv(allregions_output,index=False,float_format='%.4f')
    full_df = additionalAV45Calculations(df, lut_table, keys=ADNI_FIELDNAMES)
    full_df.to_csv(regular_output,index=False,float_format='%.4f')
    mergeRegularWithAllRegions(regular_output, allregions_output, merged_output)

    # ADNI AV45 NONTP
    regular_output = os.path.join(output_folder, 'UCBERKELEYAV45_%s_regular_nontp.csv' % (timestamp))
    allregions_output = os.path.join(output_folder, 'UCBERKELEYAV45_%s_allregions_nontp.csv' % (timestamp))
    merged_output = os.path.join(output_folder, 'UCBERKELEYAV45_%s_merged_nontp.csv' % (timestamp))
    bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes = findPreprocessOutputFiles(adni_av45_preprocess_folder, nontp=True)
    df = aggregateAllRegionFiles(bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, lut_table, adni_av45_pet_dates, adni_registry)
    df.to_csv(allregions_output,index=False,float_format='%.4f')
    full_df = additionalAV45Calculations(df, lut_table, keys=ADNI_FIELDNAMES_EXTRA)
    full_df.to_csv(regular_output,index=False,float_format='%.4f')
    mergeRegularWithAllRegions(regular_output, allregions_output, merged_output)

    # ADNI AV45 TP
    regular_output = os.path.join(output_folder, 'UCBERKELEYAV45_%s_regular_tp.csv' % (timestamp))
    allregions_output = os.path.join(output_folder, 'UCBERKELEYAV45_%s_allregions_tp.csv' % (timestamp))
    merged_output = os.path.join(output_folder, 'UCBERKELEYAV45_%s_merged_tp.csv' % (timestamp))
    bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes = findPreprocessOutputFiles(adni_av45_preprocess_folder, nontp=False)
    df = aggregateAllRegionFiles(bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, lut_table, adni_av45_pet_dates, adni_registry)
    df.to_csv(allregions_output,index=False,float_format='%.4f')
    full_df = additionalAV45Calculations(df, lut_table, keys=ADNI_FIELDNAMES_EXTRA)
    full_df.to_csv(regular_output,index=False,float_format='%.4f')
    mergeRegularWithAllRegions(regular_output, allregions_output, merged_output)

    # ADNI TAU TP
    regular_output = os.path.join(output_folder, 'UCBERKELEYTAU_%s_regular_tp.csv' % (timestamp))
    allregions_output = os.path.join(output_folder, 'UCBERKELEYTAU_%s_allregions_tp.csv' % (timestamp))
    merged_output = os.path.join(output_folder, 'UCBERKELEYTAU_%s_merged_tp.csv' % (timestamp))
    bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes = findPreprocessOutputFiles(adni_tau_preprocess_folder, nontp=False)
    df = aggregateAllRegionFiles(bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, lut_table, adni_tau_pet_dates, adni_registry)
    df.to_csv(allregions_output,index=False,float_format='%.4f')
    full_df = additionalTauCalculations(df, lut_table, keys=TAU_FIELDNAMES)
    full_df.to_csv(regular_output,index=False,float_format='%.4f')
    mergeRegularWithAllRegions(regular_output, allregions_output, merged_output)

    # DOD AV45 NONTP
    regular_output = os.path.join(output_folder, 'UCBERKELEYAV45_DOD_%s_regular_nontp.csv' % (timestamp))
    allregions_output = os.path.join(output_folder, 'UCBERKELEYAV45_DOD_%s_allregions_nontp.csv' % (timestamp))
    merged_output = os.path.join(output_folder, 'UCBERKELEYAV45_DOD_%s_merged_nontp.csv' % (timestamp))
    bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes = findPreprocessOutputFiles(dod_av45_preprocess_folder, nontp=True)
    df = aggregateAllRegionFiles(bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, lut_table, dod_av45_pet_dates, dod_registry, dod=True)
    df.to_csv(allregions_output,index=False,float_format='%.4f')
    full_df = additionalAV45Calculations(df, lut_table, keys=DOD_FIELDNAMES_EXTRA)
    full_df.to_csv(regular_output,index=False,float_format='%.4f')
    mergeRegularWithAllRegions(regular_output, allregions_output, merged_output, dod=True)

    # DOD TAU TP
    regular_output = os.path.join(output_folder, 'UCBERKELEYTAU_DOD_%s_regular_tp.csv' % (timestamp))
    allregions_output = os.path.join(output_folder, 'UCBERKELEYTAU_DOD_%s_allregions_tp.csv' % (timestamp))
    merged_output = os.path.join(output_folder, 'UCBERKELEYTAU_DOD_%s_merged_tp.csv' % (timestamp))
    bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes = findPreprocessOutputFiles(dod_tau_preprocess_folder, nontp=False)
    df = aggregateAllRegionFiles(bl_means, v2_means, v3_means, bl_sizes, v2_sizes, v3_sizes, lut_table, dod_tau_pet_dates, dod_registry, dod=True)
    df.to_csv(allregions_output,index=False,float_format='%.4f')
    full_df = additionalTauCalculations(df, lut_table, keys=DOD_TAU_FIELDNAMES)
    full_df.to_csv(regular_output,index=False,float_format='%.4f')
    mergeRegularWithAllRegions(regular_output, allregions_output, merged_output, dod=True)

    ###############
    ### OLD CODE ##
    ###############

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

