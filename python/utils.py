import sys
import csv
import os
import errno
from collections import defaultdict, Counter
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import itertools
import numpy as np
import pandas as pd
from ggplot import *

import scipy.io as sio
from scipy.spatial.distance import euclidean
from scipy.stats import linregress, norm
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d

from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.mixture import GMM
import statsmodels.api as sm

#plt.style.use('ggplot')
#pd.options.display.mpl_style = 'default'


# DEFINITIONS FOR THE COMPOSITE REGION OF INTEREST
FRONTAL=[1003,1012,1014,1018,1019,1020,1027,1028,1032,2003,2012,2014,2018,2019,2020,2027,2028,2032]
PARIETAL=[1008,1025,1029,1031,2008,2025,2029,2031]
TEMPORAL=[1015,1030,2015,2030]
CINGULATE=[1002,1010,1023,1026,2002,2010,2023,2026]
WHOLECEREBELLUM = [7,8,46,47]
CEREBG=[8,47]
CEREBW=[7,46]
CEREBL=[7,8]
CEREBR=[46,47]
BIGREF = [2,41,7,8,46,47,16]
WHITEMATTER = [2,41]
SUMMARY = FRONTAL+PARIETAL+TEMPORAL+CINGULATE

# DEFINITIONS FOR BRAAK STAGES
BRAAK1 = [1006,2006]
BRAAK2 = [17,53]
BRAAK3 = [1016,1007,1013,18,2016,2007,2013,54]
BRAAK4 = [1015,10,1002,1026,1023,1010,1035,1009,1033,2015,49,2002,2026,2023,2010,2035,2009,2033]
BRAAK5 = [1028,1012,1014,1032,1003,1027,1018,1019,1020,11,12,1011,1031,1008,1030,13,1029,1025,1001,26,1034,2028,2012,2014,2032,2003,2027,2018,2019,2020,50,51,2011,2031,2008,2030,52,2029,2025,2001,58,2034]
BRAAK6 = [1021,1022,1005,1024,1017,2021,2022,2005,2024,2017]

# DEFINITIONS FOR LOBES (FOR ANALYZING SPATIAL PATTERNS)
FRONTAL_LOBE = [1028,2028,1003,2003,1027,2027,1032,2032,1012,2012,1014,2014,1020,2020,1019,2019,1018,2018]
PARIETAL_LOBE = [1029,2029,1008,2008,1025,2025,1031,2031]
CINGULATE_LOBE = [1002,2002,1023,2023,1026,2026,1010,2010]
TEMPORAL_LOBE = [1015,2015,1030,2030,1009,2009,1034,2034,1033,2033,1001,2001]
OCCIPITAL_LOBE = [1011,2011]
MEDIAL_OCCIPITAL = [1005,2005,1021,2021,1013,2013,1007,2007]
SENSORY_LOBE = [1024,2024,1022,2022,1017,2017]
BASAL_GANGLIA = [11,50,12,51,13,52,26,58]
THALAMUS = [10,49,28,60]
LIMBIC = [17,53,1016,2016,1006,2006,1035,2035,18,54]
LOBES = {'FRONTAL': FRONTAL_LOBE,
         'PARIETAL': PARIETAL_LOBE,
         'CINGULATE': CINGULATE_LOBE,
         'TEMPORAL': TEMPORAL_LOBE,
         'OCCIPITAL': OCCIPITAL_LOBE,
         'MEDIALOCCIPITAL': MEDIAL_OCCIPITAL,
         'SENSORY': SENSORY_LOBE,
         'BASALGANGLIA': BASAL_GANGLIA,
         'THALAMUS': THALAMUS,
         'LIMBIC': LIMBIC,
         'CEREBELLUM_GRAY': CEREBG,
         'CEREBELLUM_WHITE': CEREBW,
         'CEREBRAL_WHITE': WHITEMATTER}


# BLACKLISTED REGIONS
REGION_BLACKLIST = [0,30,62,80,81,82,77,251,252,253,254,255,1000,2000,1004,2004,85,24,14,15,72,4,43,75,76]

'''
import pandas as pd; from utils import fitGMM_1D
df = pd.read_csv('../pvcsummary_bl_patterns_and_change.csv')
fitGMM_1D(df['ctx-lh-middletemporal_change'],2,graph=True)
'''

'''
import pandas as pd; from utils import scatterWithConfidence
df = pd.read_json('../pvcsummary_bl_vs_change.json')
keys = [_.replace('_bl','') for _ in df.columns if '_bl' in _ and 'ctx' in _]
for k in keys:
    # xkey = '%s_bl' % k
    xkey = 'cortical_summary_bl'
    ykey = '%s_change' % k
    label = k
    scatterWithConfidence(df, label, xkey, ykey, '../plots/%s_pvc_bl_vs_change.png' % k)
'''
'''
import pandas as pd; from utils import fitLine; import numpy as np; import matplotlib.pyplot as plt
df = pd.read_json('../pvcsummary_bl_vs_change.json')
df = df[df['cortical_summary_bl']<1.36]
df = df[df['cortical_summary_bl']>0.5]
keys = [_.replace('_bl','') for _ in df.columns if '_bl' in _ and 'ctx' in _]
xnew = np.linspace(df['cortical_summary_bl'].min(),df['cortical_summary_bl'].max(), 200)
plt.figure(1)
groups = {'high':[],'mid':[],'low':[]}
testpts = [0.50131866629999999,
           0.5500176993,
           0.60029965890000003,
           0.65183440459999997,
           0.70266162240000007,
           0.75155736560000008,
           0.80735474860000001,
           0.84944201800000008,
           0.90135566740000006,
           0.95089355310000001,
           1.0093468031999999,
           1.0530114480999999,
           1.1023878562,
           1.1575373425,
           1.2002348844999999,
           1.2509682801999999,
           1.3011045959,
           1.3513167401000001]
test_ranks = {_:[] for _ in testpts}
halfway_line = plt.axvline(x=0.897, color='k', label='Threshold')
for k in keys:
    xkey = 'cortical_summary_bl'
    ykey = '%s_change' % k
    label = k
    xfit,yfit = fitLine(df, xkey, ykey, order=0)
    # test = 0.57305707080000001
    # ypos = [x for x,y in zip(xfit,yfit) if y >= 0][0]
    # positivity[label] = ypos
    # if ytest > 0.0:
    #     groups['high'].append(label)
    #     color = 'r'
    # elif ytest > -0.00005:
    #     groups['mid'].append(label)
    #     color = 'k'
    # else:
    #     groups['low'].append(label)
    #     color = 'b'
    color = 'k'
    ytest = None
    for i, cur_test in enumerate(testpts):
        cur_ytest = [y for x,y in zip(xfit,yfit) if x == cur_test][0]
        if i == 1:
            ytest = cur_ytest
        test_ranks[cur_test].append((cur_ytest,k))
    plt.plot(xfit,yfit,label=label)


rows = {}
for i in testpts:
    print "BL STATUS: %s" % i
    sortedvals = sorted(test_ranks[i],key=lambda x: x[0], reverse=True)
    new_row = dict(enumerate([b for a,b in sortedvals]))
    new_row_vals = dict(enumerate([a for a,b in sortedvals]))
    rows['%.2f' % i] = new_row
    rows['%.2f_val' % i] = new_row_vals

output = '../regionalchange_order0.csv'
out_df = pd.DataFrame(rows)
out_df.to_csv(output)

'''

def updateLine(old_line, new_data, extraction_fn, 
               pid_key='RID', pet_meta=None, decimal_places=4):
    try:
        subj = int(old_line[pid_key])
    except Exception as e:
        print "No subject column %s found" % pid_key
        return {}

    subj_row = new_data.get(subj,None)
    if subj_row is None:
        # print "No subj row found for %s" % (subj)
        return {}

    patient_pets = None
    if pet_meta is not None:
        patient_pets = sorted(pet_meta.get(subj,[]))

    new_data = extraction_fn(subj, subj_row, old_line, patient_pets) # patient_pets is passed in as context
    new_data = convertToCSVDataType(new_data, decimal_places=decimal_places)
    return new_data




def regCoeffZTest(B1, B2, se1, se2):
    pooledSE = np.sqrt(se1 + se2)
    z = (B1 - B2) / pooledSE
    pval = 2*(1 - norm.cdf(abs(z)))
    return (z, pval)


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def isnan(arg):
    try:
        if np.isnan(arg):
            return True
    except Exception as e:
        pass
    return False

def fitLine(df, xkey, ykey, order=0):
    curdf = df[[xkey, ykey]]
    curdf = curdf.sort(columns=xkey)
    x = curdf[xkey].tolist()
    y = curdf[ykey].tolist()
    
    # use gaussian filter
    gf = gaussian_filter1d(y, sigma=50, order=order, mode='reflect')
    gauss_x, gauss_y = (x, gf)

    # # use lowess
    # lowess = sm.nonparametric.lowess(gauss_y,gauss_x,frac=0.35)
    # lowess_x, lowess_y = (lowess[:,0], lowess[:,1])

    return (gauss_x, gauss_y)

def scatterWithConfidence(df, label, xkey, ykey, output):
    p = ggplot(aes(x=xkey, y=ykey), data=df) + \
                 geom_point(color='lightblue') + \
                 stat_smooth(color='black', se=True) + \
                 ggtitle(label)
    ggsave(p, output)


def fitGMM_1D(data, components, graph=False):
    wrapped = [[_] for _ in data]
    g = GMM(n_components=components,
            covariance_type='full',
            params='wmc',
            init_params='wmc')
    g.fit(wrapped)

    if graph:
        plt.figure()
        ax = plt.subplot(1,1,1)

        # assume one dimensional data
        means = [_[0] for _ in g.means_]
        covars = [_[0] for _ in g.covars_]
        weights = g.weights_
        stds = [np.sqrt(_[0]) for _ in covars]
        data.plot(kind='kde', legend=True)
        x = np.linspace(data.min(), data.max(), 100)
        for m,s,w in zip(means,stds,weights):
            f = norm.pdf(x,m,s)*w
            plt.plot(x,f,label="Norm(%.3f, %.3f), Weight: %.3f" % (m,s,w))
        # find halfway line
        pdiffs = []
        for a in np.linspace(data.min(), data.max(), 1000):
            pone, ptwo = tuple(g.predict_proba([a])[0])
            if a > 0 and pone > 0.4 and ptwo > 0.4 and pone < 0.6 and ptwo < 0.6:
                pdiffs.append((a,abs(pone-ptwo)))
        halfway = sorted(pdiffs,key=lambda x: x[1])[0]
        print halfway
        halfway_line = plt.axvline(x=halfway[0], color='r', label='50/50 Split: %.3f' % halfway[0])
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
        plt.show()

    return (g.weights_, g.means_, g.covars_)


def apk(actual, predicted, k=10):
    """
    Computes the average precision at k.
    This function computes the average prescision at k between two lists of
    items.
    Parameters
    ----------
    actual : list
             A list of elements that are to be predicted (order doesn't matter)
    predicted : list
                A list of predicted elements (order does matter)
    k : int, optional
        The maximum number of predicted elements
    Returns
    -------
    score : double
            The average precision at k over the input lists
    """
    if len(predicted)>k:
        predicted = predicted[:k]

    score = 0.0
    num_hits = 0.0

    for i,p in enumerate(predicted):
        if p in actual and p not in predicted[:i]:
            num_hits += 1.0
            score += num_hits / (i+1.0)

    if not actual:
        return 0.0

    return score / min(len(actual), k)

def mapk(actual, predicted, k=10):
    """
    Computes the mean average precision at k.
    This function computes the mean average prescision at k between two lists
    of lists of items.
    Parameters
    ----------
    actual : list
             A list of lists of elements that are to be predicted 
             (order doesn't matter in the lists)
    predicted : list
                A list of lists of predicted elements
                (order matters in the lists)
    k : int, optional
        The maximum number of predicted elements
    Returns
    -------
    score : double
            The mean average precision at k over the input lists
    """
    return np.mean([apk(a,p,k) for a,p in zip(actual, predicted)])

def dcg(relevances, rank=20): 
    relevances = np.asarray(relevances)[:rank] 
    n_relevances = len(relevances) 
    if n_relevances == 0: 
        return 0. 
    discounts = np.log2(np.arange(n_relevances) + 2) 
    return np.sum(relevances / discounts) 

def ndcg(relevances, rank=20): 
    best_dcg = dcg(sorted(relevances, reverse=True), rank) 
    if best_dcg == 0: 
        return 0. 
    return dcg(relevances, rank) / best_dcg

def importRoussetCSV(rousset_csv, translate_threshold=None):
    '''
    If include_threshold, translates to the new PVC threshold
    '''
    df = pd.read_csv(rousset_csv)

    slope_points = []
    by_subj = {}
    for subject, subj_rows in df.groupby('subject'):
        tp_results = {}
        rid = int(subject.split('-')[-1])
        for timepoint, tp_rows in subj_rows.groupby('timepoint'):
            pvc_df = tp_rows[['name','pvcval','nonpvcval','groupsize']]
            pvc_dict = pvc_df.set_index('name').to_dict(orient='index')
            pvc_suvr = float(pvc_dict['WHOLECEREB']['pvcval'])/float(pvc_dict['COMPOSITE']['pvcval'])
            nonpvc_suvr = float(pvc_dict['WHOLECEREB']['nonpvcval'])/float(pvc_dict['COMPOSITE']['nonpvcval'])
            slope_points.append((nonpvc_suvr,pvc_suvr))
            tp_results[timepoint] = pvc_dict
        by_subj[rid] = tp_results

    if translate_threshold is not None:
        # determine threshold
        slope, intercept, r, p, stderr = linregress([_[0] for _ in slope_points], [_[1] for _ in slope_points])
        p = np.poly1d([slope, intercept])
        new_threshold = round(p(translate_threshold),2)
        return (by_subj, new_threshold)
    else:
        return (by_subj, None)


def importRoussetResults(rousset_mat):
    data = loadMATFile(rousset_mat)
    result_list = [parseResult(_) for _ in data['result_list'][0]]
    subj_list = data['subj_list'][0]
    group_names = list(data['group_names'][0])
    sorted_data = dict(zip(subj_list,result_list))
    return (group_names, sorted_data)

def importRawRoussetResults(rousset_mat):
    data = loadMATFile(rousset_mat)
    result_list = [parseResult(_) for _ in data['result_list'][0]]
    subj_list = data['subj_list'][0]
    sorted_data = dict(zip(subj_list,result_list))
    return sorted_data


def parseAllRegionOutput(all_region_file, lut_file):
    blacklist = ['RID', 'EXAMDATE']
    lut_table = importFreesurferLookup(lut_file)
    lut_reverse = {v.upper().replace('-','_'):k for k,v in lut_table.iteritems()}
    df = pd.read_csv(all_region_file)
    df = df.fillna(0.0)
    index_lookup = {k:[lut_reverse[k]] for k in df.columns if k not in blacklist and 'SIZE' not in k}

    wcereb_names = [lut_table[_].upper().replace('-','_') for _ in WHOLECEREBELLUM]
    wcereb_sizes = ['%s_SIZE' % _ for _ in wcereb_names]
    bigref_names = [lut_table[_].upper().replace('-','_') for _ in BIGREF]
    bigref_sizes = ['%s_SIZE' % _ for _ in bigref_names]
    summary_names = [lut_table[_].upper().replace('-','_') for _ in SUMMARY]
    summary_sizes = ['%s_SIZE' % _ for _ in summary_names]
    wm_names = [lut_table[_].upper().replace('-','_') for _ in WHITEMATTER]
    wm_sizes = ['%s_SIZE' % _ for _ in wm_names]

    # calculate aggregate regions
    for i in df.index:
        wcereb_pts = [(df.loc[i,size],df.loc[i,name]) for name,size in zip(wcereb_names, wcereb_sizes)]
        df.loc[i,'whole_cerebellum'] = weightedMean(wcereb_pts)
        index_lookup['whole_cerebellum'] = WHOLECEREBELLUM
        bigref_pts = [(df.loc[i,size],df.loc[i,name]) for name,size in zip(bigref_names, bigref_sizes)]
        df.loc[i,'composite_ref'] = weightedMean(bigref_pts)
        index_lookup['composite_ref'] = BIGREF
        summary_pts = [(df.loc[i,size],df.loc[i,name]) for name,size in zip(summary_names, summary_sizes)]
        df.loc[i,'cortical_summary'] = weightedMean(summary_pts)
        index_lookup['cortical_summary'] = SUMMARY
        wm_pts = [(df.loc[i,size],df.loc[i,name]) for name,size in zip(wm_names, wm_sizes)]
        df.loc[i,'white_matter'] = weightedMean(wm_pts)
        index_lookup['white_matter'] = WHITEMATTER
        df.loc[i,'EXAMDATE'] = parseDate(df.loc[i,'EXAMDATE'])

    data_bl = {}
    data_scan2 = {}
    data_scan3 = {}
    grouped = df.groupby(by=['RID'])
    for rid, rows in grouped:
        rows = rows.sort(columns=['EXAMDATE']).reset_index(drop=True)
        if len(rows)>=1:
            data_bl[rid] = {k:v for k,v in dict(rows.iloc[0]).iteritems() if k not in blacklist and 'SIZE' not in k}
        if len(rows)>=2:
            data_scan2[rid] = {k:v for k,v in dict(rows.iloc[1]).iteritems() if k not in blacklist and 'SIZE' not in k}
        if len(rows)>=3:
            data_scan3[rid] = {k:v for k,v in dict(rows.iloc[2]).iteritems() if k not in blacklist and 'SIZE' not in k}
    return data_bl, data_scan2, data_scan3, index_lookup


def parseRawAV45Output(av45_file, registry_file, lut_file):
    lut_table = importFreesurferLookup(lut_file)
    registry =  importRegistry(registry_file)
    data = importAV45(av45_file, registry=registry)

    # save column names
    index_lookup = {lut_table[code].upper().replace('-','_'):[code] for code in SUMMARY}

    # parse data
    data_bl = {}
    data_scan2 = {}
    data_scan3 = {}
    for rid, rows in data.iteritems():
        rows = sorted(rows, key=lambda x: x['EXAMDATE'])
        if len(rows) > 0:
            data_row = rows[0]
            data_bl[rid] = {k: float(data_row[k] or 0) for k in index_lookup}
            data_bl[rid]['whole_cerebellum'] = float(data_row['WHOLECEREBELLUM'])
        if len(rows) > 1:
            data_row = rows[1]
            data_scan2[rid] = {k: float(data_row[k] or 0) for k in index_lookup}
            data_scan2[rid]['whole_cerebellum'] = float(data_row['WHOLECEREBELLUM'])
        if len(rows) > 2:
            data_row = rows[2]
            data_scan3[rid] = {k: float(data_row[k] or 0) for k in index_lookup}
            data_scan3[rid]['whole_cerebellum'] = float(data_row['WHOLECEREBELLUM'])

    return data_bl, data_scan2, data_scan3, index_lookup

def bilateralTranslations(lut_file):
    lut_table = importFreesurferLookup(lut_file, flip=True)
    bilateral_dict = defaultdict(list)
    for k,v in lut_table.iteritems():
        bilateral_key = k.lower().replace('-','_').replace('lh.','').replace('rh.','').replace('lh_','').replace('rh_','').replace('right_','').replace('left_','').upper()
        bilateral_dict[bilateral_key].append(v)

    # add lobe names and keys
    bilateral_dict = dict(bilateral_dict)
    bilateral_dict.update(LOBES)
    
    return bilateral_dict

def parseRawRousset_inner(data, translations=None):
    names = [unwrap(_) for _ in data['names']]
    sizes = [unwrap(_) for _ in data['sizes']]
    pvcvals = [unwrap(_) for _ in data['pvcvals']]
    indices = [unwrap(_) for _ in data['indices']]
    indices = [[_] if (not (isinstance(_, list) or isinstance(_, np.ndarray))) else _ for _ in indices]

    val_lookup = dict(zip(names, pvcvals))
    index_lookup = dict(zip(names, indices))
    size_lookup = dict(zip(names, sizes))

    og_names = set(index_lookup.keys())
    used_names = set()
    if translations:
        used_indices = set([b for a in indices for b in a])
        for k,v in translations.iteritems():
            relevant_names = [name for name, idx in index_lookup.iteritems() if len(set(idx) & set(v))>0]
            if len(relevant_names) == 0:
                continue
            regionsizes = [(size_lookup[name],val_lookup[name]) for name in relevant_names]
            val_lookup[k] = weightedMean(regionsizes)
            index_lookup[k] = [_ for name in relevant_names for _ in index_lookup[name]]
            size_lookup[k] = sum([size_lookup[name] for name in relevant_names])
            for name in relevant_names:
                used_names.add(name)

    # only keep translations and unused names
    keep_keys = (og_names - used_names) | set(translations.keys())
    index_lookup = {k:index_lookup[k] for k in index_lookup.keys() if k in keep_keys}
    val_lookup = {k:val_lookup[k] for k in val_lookup.keys() if k in keep_keys}
    size_lookup = {k:size_lookup[k] for k in size_lookup.keys() if k in keep_keys}

    # calculate weighted sum fields
    names_and_keys = [('whole_cerebellum', WHOLECEREBELLUM),
                      ('composite_ref', BIGREF),
                      ('cortical_summary', SUMMARY)]
    names_and_keys += LOBES.items()
    for name, keys in names_and_keys:
        regionsizes = [(sz,pv) for i,sz,pv in zip(indices, sizes, pvcvals) if len(set(i) & set(keys)) > 0]
        val_lookup[name] = weightedMean(regionsizes)
        index_lookup[name] = keys
    

    return val_lookup, index_lookup

def parseRawRousset(bl_file, scan2_file, scan3_file, translations=None):
    data_bl_raw = importRawRoussetResults(bl_file)
    data_scan2_raw = importRawRoussetResults(scan2_file)
    data_scan3_raw = importRawRoussetResults(scan3_file)
    data_bl = {}
    data_scan2 = {}
    data_scan3 = {}
    index_lookup = None
    for rid, data in data_bl_raw.iteritems():
        val_lookup, index_lookup = parseRawRousset_inner(data, translations=translations)
        data_bl[rid] = val_lookup
    for rid, data in data_scan2_raw.iteritems():
        val_lookup, index_lookup = parseRawRousset_inner(data, translations=translations)
        data_scan2[rid] = val_lookup
    for rid, data in data_scan3_raw.iteritems():
        val_lookup, index_lookup = parseRawRousset_inner(data, translations=translations)
        data_scan3[rid] = val_lookup
    return data_bl, data_scan2, data_scan3, index_lookup

def removeBlacklistedGroups(data_bl, data_scan2, data_scan3, index_lookup, suvr=True, ref_key='whole_cerebellum'):
    assert ref_key in ['whole_cerebellum','composite_ref']

    regions_to_remove = [k for k,v in index_lookup.iteritems() if len(list(set(REGION_BLACKLIST) & set(list(v)))) > 0]

    index_lookup = {k:v for k,v in index_lookup.iteritems() if k not in regions_to_remove}
    for rid,val in data_bl.iteritems():
        ref_value = 1.0
        if suvr:
            ref_value = val.get(ref_key,1.0)
        new_val = {k: v/ref_value for k,v in val.iteritems() if k not in regions_to_remove}
        data_bl[rid] = new_val
    for rid,val in data_scan2.iteritems():
        ref_value = 1.0
        if suvr:
            ref_value = val.get(ref_key,1.0)
        new_val = {k:v/ref_value for k,v in val.iteritems() if k not in regions_to_remove}
        data_scan2[rid] = new_val
    for rid,val in data_scan3.iteritems():
        ref_value = 1.0
        if suvr:
            ref_value = val.get(ref_key,1.0)
        new_val = {k:v/ref_value for k,v in val.iteritems() if k not in regions_to_remove}
        data_scan3[rid] = new_val
    return data_bl, data_scan2, data_scan3, index_lookup

def unwrap(arg):
    try:
        while len(arg) == 1 and type(arg) in set([list, tuple, np.ndarray]):
            arg = arg[0]
    except Exception as e:
        pass
    return arg

def zipFields(args):
    '''
    For structs imported from mat files,
    where field values are in a list and the field names are listed under dtype
    '''
    fields = args.dtype.names
    values = []
    for i,v in enumerate(args):
        try:
            v = unwrap(v)
            if len(v) > 1:
                v = zipFields(v)
        except Exception as e:
            pass
        values.append(v)
    return dict(zip(fields,values))

def parseResult(arg):
    arg = unwrap(arg)
    data = zipFields(arg)
    return data

def gap(data, nrefs=20, ks=range(10,70), use_pca=True):
    """
    Compute the Gap statistic for an nxm dataset in data.
    Either give a precomputed set of reference distributions in refs as an (n,m,k) scipy array,
    or state the number k of reference distributions in nrefs for automatic generation with a
    uniformed distribution within the bounding box of data.
    Give the list of k-values for which you want to compute the statistic in ks.
    """
    shape = data.shape

    if use_pca:
        pca_model = PCA(n_components=0.97, copy=True, whiten=False)
        pca_data = pca_model.fit_transform(data)
        pca_shape = pca_data.shape
        tops = pca_data.max(axis=0)
        bots = pca_data.min(axis=0)
        dists = np.matrix(np.diag(tops-bots))
        pca_rands = np.random.random_sample(size=(pca_shape[0],pca_shape[1],nrefs))
        rands = np.zeros(shape=(shape[0],shape[1],nrefs))
        for i in range(nrefs):
            # convert back to original space
            pca_samples = pca_rands[:,:,i]*dists+bots
            samples = pca_model.inverse_transform(pca_samples)
            rands[:,:,i] = samples
    else:
        shape = data.shape
        tops = data.max(axis=0)
        bots = data.min(axis=0)
        dists = np.matrix(np.diag(tops-bots))
        rands = np.random.random_sample(size=(shape[0],shape[1],nrefs))
        for i in range(nrefs):
            rands[:,:,i] = rands[:,:,i]*dists+bots

    gaps = np.zeros((len(ks),))
    sds = np.zeros((len(ks),))
    for i,k in enumerate(ks):
        print "k: %s" % k
        model = KMeans(n_clusters=k, n_jobs=1, copy_x=True, verbose=True)
        labels = model.fit_predict(data)
        centers = model.cluster_centers_

        print "fit data"
        disp = sum([euclidean(data[m,:],centers[labels[m],:]) for m in range(shape[0])])

        refdisps = np.zeros((rands.shape[2],))
        for j in range(rands.shape[2]):
            model = KMeans(n_clusters=k, n_jobs=1, copy_x=True, verbose=True)
            labels = model.fit_predict(rands[:,:,j])
            centers = model.cluster_centers_
            refdisps[j] = sum([euclidean(rands[m,:,j],centers[labels[m],:]) for m in range(shape[0])])
        print "fit random"
        
        gaps[i] = np.mean(np.log(refdisps))-np.log(disp)
        sds[i] = np.sqrt(1.0 + 1.0/nrefs) * np.std(np.log(refdisps))
        print gaps[i]
        print sds[i]

    valid_ks = []
    for i, g in enumerate(gaps[:-1]):
        if g >= gaps[i+1] - sds[i+1]:
            diff = g - (gaps[i+1] - sds[i+1])
            valid_ks.append((ks[i], diff))

    return gaps, valid_ks

def extractDiagnosesFromMasterData(master_data):
    diags = {}
    for rid, row in master_data.iteritems():
        diag = row['Init_Diagnosis'].strip()
        diags[rid] = diag
    return diags

def parseCSV(file_path, delimiter=',', use_second_line=False):
    if use_second_line:
        data = pd.read_csv(file_path, sep=delimiter, low_memory=False, header=[0,1])
        data.columns = data.columns.get_level_values(1)
    else:
        data = pd.read_csv(file_path, sep=delimiter, low_memory=False)
    lines = [dict(data.iloc[i].replace(np.nan, '')) for i in range(len(data))]
    headers = list(data.columns)
    return (headers, lines)

def createVisitIDLookup(lookupfile):
    headers, lines = parseCSV(lookupfile)
    data = defaultdict(dict)
    for line in lines:
        subj = int(line['RID'])
        vc = line['VISCODE']
        vc2 = line['VISCODE2']
        data[subj][vc] = vc2
    return data

def convertVisitName(vc_name):
    vc_name = ' - '.join([_.strip() for _ in vc_name.split('-')])
    if 'ADNI1/GO' in vc_name:
        vc_name = vc_name.replace('ADNI1/GO','').strip()
    if 'ADNIGO' in vc_name:
        vc_name = vc_name.replace('ADNIGO','').strip()
    if 'ADNI1' in vc_name:
        vc_name = vc_name.replace('ADNI1','').strip()
    if 'ADNI2' in vc_name:
        vc_name = vc_name.replace('ADNI2','').strip()
    if 'Cont Pt' in vc_name:
        vc_name = vc_name.replace('Cont Pt','Continuing Pt').strip()
    lookup = {'Baseline - New Pt': 'v03',
              'Month 24': 'm24',
              'Baseline': 'bl',
              'Year 4 TelCheck': 'v42',
              'Unscheduled': 'uns1',
              'Screening - New Pt': 'v01',
              'Year 5 Visit': 'v51',
              'Month 54': 'm54',
              'No Visit Defined': 'nv',
              'Year 2 TelCheck': 'v22',
              'Month 18': 'm18',
              'Screening MRI - New Pt': 'v02',
              'Month 78': 'm78',
              'Month 12': 'm12',
              'Month 72': 'm72',
              'Year 3 TelCheck': 'v32',
              'Month 30': 'm30',
              'Month 36': 'm36',
              'Year 5 TelCheck': 'v52',
              'Month 6 - New Pt': 'v05',
              'Year 4 Visit': 'v41',
              'Year 1 TelCheck': 'v12',
              'Month 3 MRI': 'm03',
              'Year 2 Visit': 'v21',
              'Year 3 Visit': 'v31',
              'Initial TelCheck - Continuing Pt': 'v07',
              'Month 3 MRI - New Pt': 'v04',
              'Month 48': 'm48',
              'Initial Visit - Continuing Pt': 'v06',
              'Screen Fail': 'f',
              'Screening MRI': 'scmri',
              'Month 6': 'm06',
              'Year 1 Visit': 'v11',
              'Month 42': 'm42',
              'Month 66': 'm66',
              'Month 60': 'm60',
              'Screening': 'sc'}
    return lookup.get(vc_name,None)

def convertToCSVDataType(new_data, decimal_places=2):
    '''
    Converts datetimes to mm/dd/yy
    Rounds floats to strings with 2 decimal places
    Converts Nones to ''
    '''
    for k in new_data.keys():
        if isinstance(new_data[k], datetime):
            new_data[k] = datetime.strftime(new_data[k], '%m/%d/%y')
        elif isinstance(new_data[k], float):
            if decimal_places is None:
                new_data[k] = str(new_data[k])
            else:
                new_data[k] = str(round(new_data[k],decimal_places))
        elif new_data[k] is None:
            new_data[k] = ''
    return new_data


def importRegistry(registry_file, include_all=False):
    df = pd.read_csv(registry_file,low_memory=False)
    df = df.loc[:,['RID', 'VISCODE', 'VISCODE2', 'EXAMDATE']]
    df['VISCODE2'] = df['VISCODE2'].fillna('')
    df['VISCODE'] = df['VISCODE'].fillna('')
    df.loc[:,'VISCODE'] = df.loc[:,'VISCODE'].apply(lambda x: x.lower().strip())
    df.loc[:,'VISCODE2'] = df.loc[:,'VISCODE2'].apply(lambda x: x.lower().strip())
    if not include_all:
        df = df[df.VISCODE!='sc'] # filter out screening visits
    df.loc[:,'EXAMDATE'] = parseAllDates(df.loc[:,'EXAMDATE'])
    df.dropna(axis=0,subset=['EXAMDATE'],inplace=True)
    grouped = df.groupby(by=['RID'])
    by_subj = {}
    for key, val in grouped:
        by_subj[key] = sorted([dict(_) for i,_ in val.iterrows()], key=lambda x: x['EXAMDATE'])
    return by_subj

def importDODRegistry(dod_registry_file):
    headers, lines = parseCSV(dod_registry_file)
    registry = defaultdict(list)
    for data in lines:
        subj = int(data['SCRNO'])
        rid = int(data['RID'])
        viscode = data['VISCODE'].strip().lower()
        if viscode.startswith('sc'):
            continue
        date_string = data['EXAMDATE']
        if date_string == '':
            continue
        date = parseDate(date_string)
        registry[subj].append({'VISCODE': viscode,
                               'EXAMDATE': date,
                               'update_stamp': data['update_stamp'],
                               'RID': rid})
    registry = dict(registry)
    for k in registry.keys():
        new_val = sorted(registry[k], key=lambda x: x['EXAMDATE'])
        registry[k] = new_val
    return registry


def importCAPS(caps_curr_file, caps_lifetime_file):
    curr_headers, curr_lines = parseCSV(caps_curr_file)
    life_headers, life_lines = parseCSV(caps_lifetime_file)
    data = defaultdict(dict)
    for line in curr_lines:
        subj = int(line['SCRNO'])
        try:
            score = int(line['CAPSSCORE'])
        except Exception as e:
            score = ''
        data[subj]['curr'] = score
    for line in life_lines:
        subj = int(line['SCRNO'])
        try:
            score = int(line['CAPSSCORE'])
        except Exception as e:
            score = ''
        data[subj]['life'] = score
    return dict(data)

def importADNIDiagnosis(diag_file, registry=None):
    diag_headers, diag_lines = parseCSV(diag_file)
    diag_by_subj = defaultdict(list)
    for line in diag_lines:
        if 'SCRNO' in line:
            subj = int(line['SCRNO'])
        else:
            subj = int(line['RID'])
        viscode = line.get('VISCODE','').strip().lower()
        viscode2 = line.get('VISCODE2','').strip().lower()
        
        examdate = line['EXAMDATE']
        change = line['DXCHANGE']

        current = line.get('DXCURREN','')
        conv = line.get('DXCONV','')
        conv_type = line.get('DXCONTYP','')
        rev_type = line.get('DXREV','')

        if change != '':
            change = str(int(change))
        if current != '':
            current = str(int(current))
        if conv != '':
            conv = str(int(conv))
        if conv_type != '':
            conv_type = str(int(conv_type)).replace('-4','').strip()
        if rev_type != '':
            rev_type = str(int(rev_type)).replace('-4','').strip()

        if examdate:
            examdate = parseDate(examdate)
        else:
            subj_listings = registry[subj]
            for listing in subj_listings:
                if listing['VISCODE'] == viscode and listing['VISCODE2'] == viscode2:
                    examdate = listing['EXAMDATE']
                    break
            if not examdate:
                print "Could not find exam date for %s (%s, %s)" % (subj, viscode, viscode2)
                continue

        '''
            DXCHANGE: (ADNI2 DIAG CHANGE): 
                1=Stable: NL to NL; 
                2=Stable: MCI to MCI; 
                3=Stable: Dementia to Dementia; 
                4=Conversion: NL to MCI; 
                5=Conversion: MCI to Dementia; 
                6=Conversion: NL to Dementia; 
                7=Reversion: MCI to NL; 
                8=Reversion: Dementia to MCI; 
                9=Reversion: Dementia to NL
            
            DXCURREN: (CURRENT DIAGNOSIS): 
                1=NL;
                2=MCI;
                3=AD
            DXCONV: 
                1=Yes - Conversion; 
                2=Yes - Reversion; 
                0=No
            DXCONTYP: 
                1=Normal Control to MCI; 
                2=Normal Control to AD; 
                3=MCI to AD
            DXREV: 
                1=MCI to Normal Control; 
                2=AD to MCI; 
                3=AD to Normal Control
        '''

        # if adni 1 coding, convert to adni2 coding
        if change == '' and current != '':
            if conv == '1':
                # conversion
                if conv_type == '1':
                    change = '4'
                elif conv_type == '2':
                    change = '6'
                elif conv_type == '3':
                    change = '5'
            elif conv == '2':
                # reversion
                if rev_type == '1':
                    change = '7'
                elif rev_type == '2':
                    change = '8'
                elif rev_type == '3':
                    change = '9'
            elif conv == '0':
                # stable
                if current == '1':
                    change = '1'
                elif current == '2':
                    change = '2'
                elif current == '3':
                    change = '3'

        if change == '':
            print "Couldn't convert to adni2 coding: %s, %s, %s" % (subj, viscode, viscode2)
            continue

        diag_by_subj[subj].append({'VISCODE': viscode,
                                   'VISCODE2': viscode2,
                                   'EXAMDATE': examdate,
                                   'change': int(change)})
    diag_by_subj = dict(diag_by_subj)
    for k,v in diag_by_subj.iteritems():
        diag_by_subj[k] = sorted(v, key=lambda x: x['EXAMDATE'])
    return diag_by_subj

def importFreesurferLookup(lut_file, flip=False):
    infile = open(lut_file, 'r')
    data = {}
    for line in infile.readlines():
        parts = [_.strip() for _ in line.split(' ') if _.strip() != '']
        if len(parts) < 2:
            continue
        try:
            idx = int(parts[0])
            if flip:
                data[parts[1]] = idx
            else:
                data[idx] = parts[1]
        except Exception as e:
            continue
    return data

def importAPOE(apoe_file):
    headers, lines = parseCSV(apoe_file)
    apoe_by_subj = {}
    for line in lines:
        subj = int(line['SCRNO'])
        apgen1 = int(line['APGEN1'])
        apgen2 = int(line['APGEN2'])
        apoe_by_subj[subj] = {'apgen1': apgen1, 'apgen2': apgen2}
    return apoe_by_subj

def importFDG(fdg_file):
    headers, lines = parseCSV(fdg_file)
    by_subj = defaultdict(list)
    for line in lines:
        subj = int(line['RID'])
        vc = line['VISCODE']
        vc2 = line['VISCODE2']
        date = parseDate(line['EXAMDATE'])
        roiname = line['ROINAME']
        mean = float(line['MEAN'])
        by_subj[subj].append({'vc': vc,
                              'vc2': vc2,
                              'EXAMDATE': date,
                              'ROI': roiname,
                              'MEAN': mean})
    return dict(by_subj)


def importExtractedFDG(fdg_file, registry):
    headers, lines = parseCSV(fdg_file, delimiter='\t')
    fdg_rows = []
    ts = datetime.strftime(datetime.now(),"%Y-%m-%d %H:%M:%S.0")
    for line in lines:
        subj = line['PTID']
        subj_id = int(subj.split('_')[-1])
        date_str = line['DATE'].strip().split(' ')[0].strip()
        examdate = parseDate(date_str)
        subject_registry = registry.get(subj_id,[])
        if len(subject_registry) == 0:
            raise Exception("NO SUBJECT REGISTRY")
        # find viscode in registry
        vc_name = line['VISCODE'].strip()
        vc = convertVisitName(vc_name)
        vc2 = None
        if vc is not None:
            # find vc2 from registry
            for regrow in subject_registry:
                if regrow['VISCODE'].strip() == vc:
                    vc2 = regrow['VISCODE2']
                    break
        if vc is None or vc2 is None:
            # find by examdate
            subj_sort = sorted(subject_registry, key=lambda x: abs(examdate-x['EXAMDATE']).days)
            closest = subj_sort[0]
            if abs(closest['EXAMDATE']-examdate).days < 90:
                vc = closest['VISCODE']
                vc2 = closest['VISCODE2']
            else:
                print "VISCODE NAME COULD NOT BE CONVERTED AND CLOSEST DATE OUT OF RANGE: %s, %s, %s" % (vc_name,vc,vc2)
        if vc is not None and vc2 is None:
            # copy over
            print "COPYING VC TO VC2"
            vc2 = vc
        if vc is None or vc2 is None:
            print "COULD NOT FIND VISCODES FOR %s, %s, %s" % (subj_id, examdate, line['VISCODE'])
            vc = ''
            vc2 = ''
        data = {'RID': subj_id,
                'VISCODE': vc,
                'VISCODE2': vc2,
                'UID': line['UID'].strip(),
                'ROINAME': line['ROINAME'].strip(),
                'ROILAT': line['ROILAT'].strip(),
                'EXAMDATE': examdate,
                'MEAN': float(line['MEAN']),
                'MEDIAN': float(line['MEDIAN']),
                'MODE': float(line['MODE']),
                'MIN': float(line['MIN']),
                'MAX': float(line['MAX']),
                'STDEV': float(line['STDEV']),
                'NANVOX': int(line['NANVOX']),
                'TOTVOX': int(line['TOTVOX']),
                'update_stamp': ts}
        fdg_rows.append(data)
    return fdg_rows


def importDemog(demog_file):
    headers, lines = parseCSV(demog_file)
    demog_by_subj = {}
    for line in lines:
        if 'SCRNO' in line:
            subj = int(line['SCRNO'])
        else:
            subj = int(line['RID'])
        gender = line['PTGENDER']
        gender = int(gender) if gender != '' else gender
        if 'PTAGE' in line:
            age = float(line['PTAGE'])
        else:
            age = None
        married = line['PTMARRY']
        married = int(married) if married != '' else married
        edu = line['PTEDUCAT']
        edu = int(edu) if edu != '' else edu
        dob = None
        if 'PTDOB' in line:
            dob = parseDate(line['PTDOB'])
        elif 'PTDOBMM' in line and 'PTDOBYY' in line:
            year = line['PTDOBYY']
            month = line['PTDOBMM']
            if year != '' and month != '':
                dob = datetime(year=int(year), month=int(month), day=1)
        demog_by_subj[subj] = {'gender': gender,
                               'age': age,
                               'married': married,
                               'edu': edu,
                               'dob': dob}
    return demog_by_subj

def importMMSE(mmse_file, registry=None):


    mmse_headers, mmse_lines = parseCSV(mmse_file)
    
    # restructure mmse lines by subject
    mmse_by_subject = defaultdict(list)
    for line in mmse_lines:
        subj = int(line['RID'])
        date_string = line['EXAMDATE']
        if not date_string and registry is not None:
            # get from registry
            viscodes = [line['VISCODE'].strip().lower(), line['VISCODE2'].strip().lower()]
            date = findVisitDate(registry, subj, viscodes)
            if date is None:
                print "Could not find visit in registry: %s, %s" % (subj, viscodes)
        else:
            date = parseDate(date_string)
        if date is not None:
            mmse_by_subject[subj].append((date,line))
    mmse_by_subject = dict(mmse_by_subject)
    for k,v in mmse_by_subject.iteritems():
        mmse_by_subject[k] = sorted(v, key=lambda x: x[0])
    return mmse_by_subject


def importAVLT(avlt_file, registry=None):
    headers, lines = parseCSV(avlt_file)
    # restructure by subject
    avlt_by_subj = defaultdict(list)
    for line in lines:
        if 'SCRNO' in line:
            subj = int(line['SCRNO'])
        else:
            subj = int(line['RID'])

        date_string = line.get('EXAMDATE')
        viscode = line['VISCODE'].strip().lower() if 'VISCODE' in line else ''
        viscode2 = line['VISCODE2'].strip().lower() if 'VISCODE2' in line else ''
        if not date_string and registry is not None:
            # get from registry
            viscodes = [viscode, viscode2]
            date = findVisitDate(registry, subj, viscodes)
            if date is None:
                print "Could not find visit in registry: %s, %s" % (subj, viscodes)
                continue
        else:
            date = parseDate(date_string)

        tots = [line['AVTOT%s' % _ ] for _ in range(1,6)]
        try:
            score_sum = 0.0
            for score_str in tots:
                new_score = float(score_str)
                if new_score < 0:
                    raise Exception("Invalid score part")
                score_sum += new_score
            test_score = score_sum
        except Exception as e:
            continue
        avlt_by_subj[subj].append({'VISCODE': viscode,
                                   'VISCODE2': viscode2,
                                   'EXAMDATE': date,
                                   'TOTS': test_score})
    return dict(avlt_by_subj)

def importADASCog(adni1_file, adnigo2_file, registry=None):
    try:
        adni1_headers, adni1_lines = parseCSV(adni1_file)
    except Exception as e:
        adni1_lines = []
    try:
        adni2_headers, adni2_lines = parseCSV(adnigo2_file)
    except Exception as e:
        adni2_lines = []
    # restructure by subject
    adas_by_subj = defaultdict(list)
    for line in adni1_lines:
        subj = int(line['RID'])
        viscode = line['VISCODE'].strip().lower()
        examdate = parseDate(line['EXAMDATE'])
        totscore = float(line['TOTAL11'])
        if totscore < 0:
            totscore = ''
        adas_by_subj[subj].append({'VISCODE': viscode,
                                   'VISCODE2': None,
                                   'EXAMDATE': examdate,
                                   'TOTSCORE': totscore})
    for line in adni2_lines:
        if 'SCRNO' in line:
            subj = int(line['SCRNO'])
        else:
            subj = int(line['RID'])
        viscode = line['VISCODE'].strip().lower()
        try:
            viscode2 = line['VISCODE2'].strip().lower()
        except Exception as e:
            viscode2 = ''
        raw_totscore = line['TOTSCORE']
        if raw_totscore == '':
            print "%s (%s, %s) has missing totscore" % (subj, viscode, viscode2)
            continue
        totscore = int(raw_totscore)
        if totscore < 0:
            totscore = ''
        examdate = None
        if registry is not None:
            examdate = findVisitDate(registry, subj, [viscode, viscode2])
        if examdate is None:
            print "No exam date for %s (%s, %s)" % (subj, viscode, viscode2)
            continue
        adas_by_subj[subj].append({'VISCODE': viscode,
                                   'VISCODE2': viscode2,
                                   'EXAMDATE': examdate,
                                   'TOTSCORE': totscore})
    return dict(adas_by_subj)

def importTBMSyn(tbm_file, registry=None):
    df = pd.read_csv(tbm_file, low_memory=False)
    df['EXAMDATE'] = df['EXAMDATE'].fillna('')
    df['EXAMDATEBL'] = df['EXAMDATEBL'].fillna('')

    for i in df.index:
        # parse bl examdate
        bl_date = df.loc[i,'EXAMDATEBL']
        df.loc[i,'EXAMDATEBL'] = parseDate(bl_date)
        # parse examdate
        raw_date = df.loc[i,'EXAMDATE']
        if raw_date == '':
            # try to find date from viscode
            date = findVisitDate(registry, df.loc[i,'RID'], [df.loc[i,'VISCODE']]) if registry is not None else None
            if date is None:
                date = np.nan
        else:
            date = parseDate(raw_date)
        df.loc[i,'EXAMDATE'] = date
        
    df.dropna(axis=0,subset=['EXAMDATE'],inplace=True)
    grouped = df.groupby(by=['RID'])
    by_subj = {}
    for key, val in grouped:
        by_subj[key] = [dict(_) for i,_ in val.iterrows()]
    return by_subj

def importWMH(wmh_file):
    headers, lines = parseCSV(wmh_file)
    data = defaultdict(list)
    for line in lines:
        # get subject ID
        try:
            if 'SCRNO' in line:
                subj = int(line['SCRNO'])
            else:
                subj = int(line['RID'])
        except Exception as e:
            continue
        vc = line['VISCODE'].strip().lower()
        try:
            vc2 = line['VISCODE2'].strip().lower()
        except:
            vc2 = '' 
        examdate = parseDate(line['EXAMDATE'])
        wmh = float(line['WHITMATHYP'])
        icv = float(line['ICV'])
        wmh_percent = wmh/icv
        data[subj].append({'vc': vc,
                           'vc2': vc2,
                           'EXAMDATE': examdate,
                           'wmh': wmh,
                           'icv': icv,
                           'wmh_percent': wmh_percent})
    return dict(data)

def importUW(uw_file, registry=None):
    df = pd.read_csv(uw_file, low_memory=False)
    df['EXAMDATE'] = df['EXAMDATE'].fillna('')
    for i in df.index:
        raw_date = df.loc[i,'EXAMDATE']
        if raw_date == '':
            # try to find date from viscode
            date = findVisitDate(registry, df.loc[i,'RID'], [df.loc[i,'VISCODE']]) if registry is not None else None
            if date is None:
                date = np.nan
        else:
            date = parseDate(raw_date)
        df.loc[i,'EXAMDATE'] = date
    df.dropna(axis=0,subset=['EXAMDATE'],inplace=True)
    
    grouped = df.groupby(by=['RID'])
    by_subj = {}
    for key, val in grouped:
        by_subj[key] = [dict(_) for i,_ in val.iterrows()]
    return by_subj

def importGD_DOD(gd_file):
    headers, lines = parseCSV(gd_file)
    data = defaultdict(dict)
    for line in lines:
        subj = int(line['SCRNO'])
        gdtotal = line['GDTOTAL']
        try:
            gdtotal = int(gdtotal)
        except:
            pass
        data[subj]['gdtotal']  = gdtotal
    return dict(data)

def importGD(gd_file, registry=None):
    headers, lines = parseCSV(gd_file)
    data = defaultdict(list)
    for line in lines:
        subj = int(line['RID'])
        vc = line['VISCODE'].strip()
        vc2 = line['VISCODE2'].strip()
        try:
            examdate = parseDate(line['EXAMDATE'])
        except:
            examdate = findVisitDate(registry, subj, [vc, vc2]) if registry is not None else None
        if examdate is None:
            print "Couldn't find GD examdate for %s, %s, %s" % (subj, vc, vc2)
            continue
        gdtotal = line['GDTOTAL']
        data[subj].append({'vc': vc,
                           'vc2': vc2,
                           'EXAMDATE': examdate,
                           'GDTOTAL': gdtotal})
    data = dict(data)
    return data

def importCSF(csf_files, registry=None):
    data = defaultdict(dict)
    for csffile in csf_files:
        print "Importing CSF file %s" % csffile
        headers, lines = parseCSV(csffile)
        for i, line in enumerate(lines):
            if 'SCRNO' in line: # for the dod patients
                subj = int(line['SCRNO'])
            else:
                subj = int(line['RID'])
            vc = line['VISCODE'].strip().lower()
            rundate = parseDate(line['RUNDATE'])
            try:
                vc2 = line['VISCODE2'].strip().lower()
            except:
                vc2 = ''
            try:
                examdate = parseDate(line['EXAMDATE'])
            except:
                examdate = findVisitDate(registry, subj, [vc, vc2]) if registry is not None else None
            try:
                abeta = float(line['ABETA'])
            except:
                abeta = None
            try:
                tau = float(line['TAU'])
            except:
                tau = None
            try:
                ptau = float(line['PTAU'])
            except:
                ptau = None
            visitkey = (vc,vc2)
            newdata = {'rundate': rundate,
                       'EXAMDATE': examdate,
                       'file': csffile,
                       'abeta': abeta,
                       'tau': tau,
                       'ptau': ptau,
                       'vc': vc,
                       'vc2' : vc2}
            if examdate in data[subj]:
                data[subj][examdate].append(newdata)
            else:
                data[subj][examdate] = [newdata]
    # take the timepoint with the most recent run date for each visit
    data = dict(data)
    for k in data.keys():
        userdata = data[k]
        flattened = []
        for examdate, inners in userdata.iteritems():
            chosen = sorted(inners, key=lambda x: x['rundate'])[-1]
            flattened.append(chosen)
        data[k] = sorted(flattened, key=lambda x: x['EXAMDATE'])
    return data

def importScanMeta(meta_file, with_viscode=False):
    df = pd.read_csv(meta_file)
    # get exam date column
    subj_col = None
    date_col = None
    viscode_cols = ['VISCODE']
    if 'SCRNO' in df.columns:
        subj_col = 'SCRNO'
    elif 'RID' in df.columns:
        subj_col = 'RID'
    if 'SCANDATE' in df.columns:
        date_col = 'SCANDATE'
    elif 'EXAMDATE' in df.columns:
        date_col = 'EXAMDATE'
    if 'VISCODE2' in df.columns:
        viscode_cols.append('VISCODE2')
    if subj_col is None or date_col is None:
        raise Exception("No subj/date column in %s" % meta_file)
    df.dropna(axis=0, subset=[subj_col, date_col], inplace=True)
    df.loc[:,date_col] = df.loc[:,date_col].apply(parseDate)
    scans = {}
    for pid, rows in df.groupby(subj_col):
        # IF TWO ENTRIES HAVE THE SAME VISIT CODE, TAKE THE ONE WITH THE LATER DATE
        rows = rows[[date_col] + viscode_cols]
        rows = rows.sort_values(by=date_col)
        dedup_rows = rows[rows.shift(-1) != rows]
        dedup_rows.dropna(inplace=True)
        if with_viscode:
            dedup_rows.set_index(date_col,inplace=True)
            scans[pid] = dedup_rows.to_dict(orient='index')
        else:
            scans[pid] = sorted(list(dedup_rows[date_col]))
    return scans

def importPetMETA(pet_meta_file, tracer='AV45'):
    headers, lines = parseCSV(pet_meta_file)
    pets = defaultdict(list)
    for row in lines:
        # filter
        if row['Sequence'].strip() != '%s Coreg, Avg, Std Img and Vox Siz, Uniform Resolution' % tracer:
            continue
        subj = int(row['Subject'].split('_')[-1].strip())
        new_date = parseDate(row['Scan Date'])
        pets[subj].append(new_date)
    pets = dict(pets)
    for k in pets.keys():
        pets[k] = sorted(list(set(pets[k])))
    return pets

def importARM(arm_file):
    '''
        1=NL - (ADNI1 1.5T only)
        2=LMCI - (ADNI1 1.5T only)
        3=AD - (ADNI1 1.5T only)
        4=NL - (ADNI1 PET+1.5T)
        5=LMCI - (ADNI1 PET+1.5T)
        6=AD - (ADNI1 PET+1.5T)
        7=NL - (ADNI1 3T+1.5T)
        8=LMCI - (ADNI1 3T+1.5T)
        9=AD - (ADNI1 3T+1.5T)
        10=EMCI; 
        11=SMC - (Significant Memory Concern)
    '''
    translation = {1: 'N',
                   2: 'LMCI',
                   3: 'AD',
                   4: 'N',
                   5: 'LMCI',
                   6: 'AD',
                   7: 'N',
                   8: 'LMCI',
                   9: 'AD',
                   10: 'EMCI',
                   11: 'SMC'}
    headers, lines = parseCSV(arm_file)
    arms = defaultdict(list)
    for data in lines:
        subj = int(data['RID'])
        status = data['ARM']
        if status == '':
            continue
        userdate = parseDate(data['USERDATE'])
        # convert status
        status_str = translation[status]
        arms[subj].append({'USERDATE': userdate,
                           'STATUS': status_str})
    return dict(arms)

def importDODEligibility(elig_file):
    headers, lines = parseCSV(elig_file)
    data = defaultdict(dict)
    for line in lines:
        subj = int(line['SCRNO'])
        try:
            cohort = int(line['COHORT'])
        except Exception as e:
            cohort = ''
        if cohort == '' and data[subj].get('cohort','') != '':
            # don't overwrite good value with bad one
            continue
        data[subj]['cohort'] = cohort
    return dict(data)

def importMedicalHistory(mhist_file):
    mhist_df = pd.read_csv(mhist_file)
    mhist_df.set_index('RID',inplace=True)
    mhist_df = mhist_df[['MHDESC']]
    mhist_df.loc[:,'MHDESC'] = mhist_df.MHDESC.apply(lambda x: x.lower())
    mhist_df.loc[:,'smoking'] = mhist_df.MHDESC.str.contains('smok')
    mhist_df.loc[:,'diabetes'] = mhist_df.MHDESC.str.contains('diabete')
    mhist_df.loc[:,'hyperlipidemia'] = mhist_df.MHDESC.str.contains('hyperlipidemia')
    smoking = set(mhist_df[mhist_df.smoking].index)
    diabetes = set(mhist_df[mhist_df.diabetes].index)
    hyperlipidemia = set(mhist_df[mhist_df.hyperlipidemia].index)
    by_rid = {}
    for rid in set(mhist_df.index):
        by_rid[rid] = {'smoking': rid in smoking,
                       'diabetes': rid in diabetes,
                       'hyperlipidemia': rid in hyperlipidemia}
    return by_rid

def importFAQ(faq_file, registry):
    faq_df = pd.read_csv(faq_file)
    faq_df = faq_df[['RID','EXAMDATE','VISCODE','VISCODE2','FAQTOTAL']]
    # fill in examdates
    def fillInDate(row, registry):
        if isnan(row['EXAMDATE']):
            vcs = [_ for _ in [row['VISCODE'],row['VISCODE2']] if not isnan(_)]
            found_date = findVisitDate(registry, row['RID'], vcs)
            if found_date is not None:
                return found_date
            else:
                return np.nan
        else:
            return parseDate(row['EXAMDATE'])
    faq_df['EXAMDATE'] = faq_df.apply(lambda x: fillInDate(x,registry), axis=1)
    faq_df = faq_df[['RID','EXAMDATE','FAQTOTAL']]
    faq_df.dropna(inplace=True)
    by_rid = {}
    for rid, ridrows in faq_df.groupby('RID'):
        rows = ridrows[['EXAMDATE','FAQTOTAL']].to_dict(orient='records')
        sorted_rows = sorted(rows, key=lambda x: x['EXAMDATE'])
        by_rid[rid] = sorted_rows
    return by_rid

def importNPI(npi_file):
    npi_df = pd.read_csv(npi_file)
    npi_df = npi_df[['RID','EXAMDATE','NPITOTAL']]
    npi_df.dropna(inplace=True)
    npi_df['EXAMDATE'] = npi_df.EXAMDATE.apply(lambda x: parseDate(x))
    by_rid = {}
    for rid, ridrows in npi_df.groupby('RID'):
        rows = ridrows[['EXAMDATE','NPITOTAL']].to_dict(orient='records')
        sorted_rows = sorted(rows, key=lambda x: x['EXAMDATE'])
        by_rid[rid] = sorted_rows
    return by_rid

def importDODMRI(mri_file):
    headers, lines = parseCSV(mri_file)
    data = defaultdict(list)
    for line in lines:
        subj = int(line['SCRNO'])
        try:
            conducted = int(line['MMCONDCT'])
        except Exception as e:
            continue
        if conducted == 0:
            continue
        date = parseDate(line['EXAMDATE'])
        data[subj].append(date)

    # sort
    data = dict(data)
    for k in data.keys():
        data[k] = sorted(data[k])
    return data

def importMRI(mri_file, magstrength_filter=None, filter_mprage=True):
    bad_sequences = set([])
    headers, lines = parseCSV(mri_file)
    data = defaultdict(list)
    for i, line in enumerate(lines):
        # get subject ID
        try:
            if 'SubjectID' in line:
                subj_id_whole = line['SubjectID']
            elif 'SUBJECT' in line:
                subj_id_whole = line['SUBJECT']
            else:
                raise Exception("No subj column found")
            subj = int(subj_id_whole.split('_')[-1])
        except Exception as e:
            print e
            continue

        if 'Sequence' in line:
            seq = line['Sequence'].strip().lower()
        elif 'SEQUENCE' in line:
            seq = line['SEQUENCE'].strip().lower()
        if 'accelerat' in seq:
            bad_sequences.add(seq)
            continue
        seq = seq.replace('adni','').strip()
        seq = seq.replace('mp rage', 'mprage')
        seq = seq.replace('mp-rage', 'mprage')
        seq = seq.replace('mp- rage', 'mprage')
        
        
        if filter_mprage and not ('mpr' in seq or 'spgr' in seq or 'n3m' in seq):
            bad_sequences.add(seq)
            continue
        
        if 'ScanDate' in line:
            new_date = line['ScanDate']
        elif 'SCANDATE' in line:
            new_date = line['SCANDATE']

        date = None
        try:
            date = parseDate(new_date)
        except Exception as e:
            pass
        if date is None:
            continue

        if 'Visit' in line:
            vc = line['Visit']
        else:
            vc = None

        MagStrength = line.get('MagStrength','')
        if magstrength_filter is not None:
            if MagStrength.strip().lower() != magstrength_filter.strip().lower():
                continue
        # find existing
        data[subj].append({'EXAMDATE' : date, 'vc': vc, 'strength': MagStrength})
    print bad_sequences
    return dict(data)


def importMaster(master_file):
    headers, lines = parseCSV(master_file, use_second_line=True)
    data = {}
    for i, line in enumerate(lines):
        # get subject ID
        try:
            subj = int(line['RID'])
        except Exception as e:
            print line.keys()
            continue
        data[subj] = line
    return data

def importBSI(bsi_file, include_failed=False):
    headers, lines = parseCSV(bsi_file)
    data = defaultdict(list)
    failed = 0
    bl_examdates = {}
    for i, line in enumerate(lines):
        subj = int(line['RID'])
        if not include_failed and int(line['QC_PASS']) == 0:
            failed += 1
            continue
        elif line['MRSEQUENCE'] == 'Acc':
            continue
        elif line['KMNDBCBBSI'] == '':
            # a baseline scan
            bl_examdate = parseDate(line['EXAMDATE'])
            bl_examdates[subj] = bl_examdate
            continue
        vc = line['VISCODE'].strip().lower()
        vc2 = line['VISCODE2'].strip().lower()
        examdate = parseDate(line['EXAMDATE'])
        dbcb_bsi = line['DBCBBSI']
        kmndbcb_bsi = line['KMNDBCBBSI']
        v_bsi = line['VBSI']
        h_bsi_r = line['HBSI_R']
        h_bsi_l = line['HBSI_L']
        data[subj].append({'VISCODE': vc,
                           'VISCODE2': vc2,
                           'EXAMDATE': examdate,
                           'BL_EXAMDATE': bl_examdates.get(subj,None),
                           'VBSI': v_bsi,
                           'HBSI_R': h_bsi_r,
                           'HBSI_L': h_bsi_l,
                           'WB_BSI': dbcb_bsi,
                           'WB_KNBSI': kmndbcb_bsi})
    print "BSI Failed: %s" % failed

    # fill in baseline times
    data = dict(data)
    for subj in data.keys():
        if subj in bl_examdates:
            lines = data[subj]
            for l in lines:
                l['BL_EXAMDATE'] = bl_examdates[subj]
            data[subj] = lines

    return dict(data)


def importAV1451(av1451_file):
    df = pd.read_csv(av1451_file)
    if 'SCRNO' in df.columns:
        subj_col = 'SCRNO'
    else:
        subj_col = 'RID'
    df.loc[:,'EXAMDATE'] = df.loc[:,'EXAMDATE'].apply(parseDate)
    subj_dict = {}
    for subj, rows in df.groupby(subj_col):
        subj_dict[int(subj)] = rows.set_index(subj_col).to_dict(orient='records')
    return subj_dict


def importAV45(av45_file, registry=None):
    df = pd.read_csv(av45_file)
    df.loc[:,'TP_SPECIFIC'] = True
    for i in df.index:
        raw_date = df.loc[i,'EXAMDATE']
        date = None
        if raw_date:
            date = parseDate(raw_date)
        else:
            print "Can't find examdate: %s" % i
            date = np.nan
        df.loc[i,'EXAMDATE'] = date
    # remove rows with no date
    df.dropna(axis=0,subset=['EXAMDATE'],inplace=True)

    # set id
    if 'SCRNO' in df.columns:
        group_col = 'SCRNO'
    elif 'RID' in df.columns:
        group_col = 'RID'
    elif 'PID' in df.columns:
        group_col = 'PID'
    else:
        raise Exception("Can't find subject column")
    grouped = df.groupby(by=[group_col])
    
    # convert to dictionary
    av45_by_subj = {}
    for key, val in grouped:
        av45_by_subj[key] = [dict(_) for i,_ in val.iterrows()]
    return av45_by_subj

def getMagStrength(mprage_df, imageuid):
    try:
        return mprage_df.loc[imageuid,'MagStrength']
    except Exception as e:
        print e
        return np.nan

def importFSVolumes(vol_file):
    columns = ['Subject', 'Scan', 'Date', 'MRI', 'EstimatedTotalIntraCranialVol', 'Left-Hippocampus', 'Right-Hippocampus']
    df = pd.read_csv(vol_file)
    df = df[columns]
    df.loc[:,'Date'] = df.loc[:,'Date'].apply(parseDate)
    df.loc[:,'Subject'] = df.loc[:,'Subject'].apply(lambda x: int(x.split('-')[-1]))
    df['HCV'] = df['Left-Hippocampus'] + df['Right-Hippocampus']
    df['HC/ICV'] = df['HCV'] / df['EstimatedTotalIntraCranialVol']
    df.dropna(inplace=True)

    by_subj = {subj: rows.sort_values(by='Date').to_dict(orient='record') for subj, rows in df.groupby('Subject')}
    return by_subj


def importUCSFFreesurfer(in_file, mprage_file, version='', include_failed=False, as_df=False):
    mprage_df = pd.read_csv(mprage_file)
    mprage_df.set_index('ImageUID',inplace=True)

    # pass qc
    df = pd.read_csv(in_file)
    if not include_failed:
        df = df[df['OVERALLQC'].isin(['Pass','Partial'])]
        # at least tempqc has to pass
        before = len(df.index)
        tempqc_passed = df.loc[:,'TEMPQC'] == 'Pass'
        failed = df[~tempqc_passed]
        df = df[tempqc_passed]
        failed_counter = Counter(list(failed.index))
        print "%s: FAILED TEMP QC: %s/%s" % (in_file, before-len(df.index), before)

    # filter by image type 
    if 'IMAGETYPE' in df.columns:
        df = df[df['IMAGETYPE']=='Non-Accelerated T1']

    # get relevant columns
    # right hippocampus volume: ST88SV
    # left hippocampus volume: ST29SV
    # intracranial volume: ST10CV
    key_columns = ['RID','VISCODE']
    if 'VISCODE2' in df.columns:
        key_columns.append('VISCODE2')
    filtered = df[key_columns].copy()
    parsed_dates = parseAllDates(df['EXAMDATE'].tolist())
    filtered.loc[:,'EXAMDATE'] = parsed_dates
    filtered.loc[:,'LEFT_HCV'] = df.loc[:,'ST29SV']
    filtered.loc[:,'RIGHT_HCV'] = df.loc[:,'ST88SV']
    filtered.loc[:,'HCV'] = df[['ST88SV','ST29SV']].sum(axis=1)
    filtered.loc[:,'ICV'] = df.loc[:,'ST10CV']
    filtered.loc[:,'FLDSTRENG'] = df.loc[:,'IMAGEUID'].apply(lambda x: getMagStrength(mprage_df, x))
    filtered.loc[:,'version'] = version
    filtered = filtered.dropna(subset=['EXAMDATE'])
    print 'FIELD STRENGTH DOMAIN: %s' % (set(filtered['FLDSTRENG']),)

    if as_df:
        return filtered
    else:
        data = {}
        for rid, rows in filtered.groupby('RID'):
            data[rid] = rows.to_dict(orient='records') 
        return data


def slope(points):
    '''
    Expects list of (x,y) points
    '''
    if len(points) >= 2:
        x = [_[0] for _ in points]
        y = [_[1] for _ in points]
        slope, intercept, r, p, stderr = linregress(x,y)
        return slope
    return None

def weightedMean(points):
    '''
    Expects list of (weight,value) points
    '''
    weights = np.array([_[0] for _ in points])
    values = np.array([_[1] for _ in points])
    weights = weights / float(sum(weights))
    return sum(weights*values)

def parseDate(date_str, registry=None):
    if date_str == '':
        raise Exception("Blank date string")
    date_str = date_str.replace('--','1')
    formats = ['%Y-%m-%d', '%m/%d/%y', '%m/%d/%Y']
    for f in formats:
        try:
            date = datetime.strptime(date_str,f)
            return date
        except:
            pass
    raise Exception("Cannot parse date: %s" % date_str)

def parseAllDates(date_str_list):
    parsed = []
    for d in date_str_list:
        try:
            new_date = parseDate(d)
        except:
            new_date = np.nan
        parsed.append(new_date)
    return parsed

def dumpCSV(file_path, headers, lines):
    print "DUMPING OUTPUT TO %s" % (file_path)
    writer = csv.DictWriter(open(file_path,'w'), fieldnames=headers)
    writer.writeheader()
    for l in lines:
        filtered_line = {}
        for k in headers:
            filtered_line[k] = l[k] if k in l else ''
        writer.writerow(filtered_line)

def rearrangeHeaders(new_headers, to_add, after=None):
    '''
    if after is None, then stick in the end of the headers
    '''
    for ta in to_add:
        if ta in new_headers:
            new_headers.remove(ta)
    if after is None:
        new_headers.extend(to_add)
    else:
        idx = new_headers.index(after) + 1
        new_headers = new_headers[:idx] + to_add + new_headers[idx:]
    return new_headers

def findVisitDate(registry, subj, viscodes):
    vcs = [_.lower().strip() for _ in viscodes]
    vcs = set([_ for _ in vcs if _ != ''])
    subj_registry = registry.get(subj,[])
    date = None
    for reg_row in subj_registry:
        row_visits = [reg_row.get('VISCODE','').lower().strip(), reg_row.get('VISCODE2','').lower().strip()]
        row_visits = set([_ for _ in row_visits if _ != ''])
        if len(vcs & row_visits) > 0:
            date = reg_row['EXAMDATE']
            break
    return date

def appendCSV(csv1, csv2):
    '''
    appends csv2 rows to csv1 (have to have same headers)
    '''
    h1, r1 = parseCSV(csv1)
    h2, r2 = parseCSV(csv2)
    all_rows = r1 + r2
    all_headers = h1
    dumpCSV(csv1, all_headers, all_rows)

def removeFile(filename):
    try:
        os.remove(filename)
    except OSError:
        print "File to remove does not exist"

def saveFakeAparcInput(out_file, values, index_lookup):
    sio.savemat(out_file, {'regions': [{'inds': index_lookup[k], 'val': v} for k,v in values.iteritems()]})

def printROIGrouping(names, groupings, lut_file):
    data = importFreesurferLookup(lut_file)
    for name, group in zip(names,groupings):
        print "\t%s" % name
        for g in group:
            print "\t\t%s: %s" % (g, data.get(g, 'Unknown'))

def createROIGrouping(names, indices, output_file):
    createMATFile('roigroups', [{'name': n, 'ind': d} for n,d in zip(names, indices)], output_file)

def createMATFile(name, data, output_file):
    sio.savemat(output_file, {name: data})

def loadMATFile(input_file):
    data = sio.loadmat(input_file)
    return data

def calculateCSVDifference(file1, file2, index='RID'):
    headers1, rows1 = parseCSV(file1)
    headers2, rows2 = parseCSV(file2)
    common_headers = list(set(headers1) & set(headers2))
    differences = {h:[] for h in common_headers}
    # index the rows
    index1 = {}
    for r in rows1:
        key = r.get(index)
        if key:
            index1[key.lower().strip()] = r
    index2 = {}
    for r in rows2:
        key = r.get(index)
        if key:
            index2[key.lower().strip()] = r
    common_keys = list(set(index1.keys()) & set(index2.keys()))
    for key, header in itertools.product(common_keys, common_headers):
        cur_diff = float(index1[key][header]) - float(index2[key][header])
        differences[header].append(cur_diff)
    diff_distr = {}
    for k,values in differences.iteritems():
        diff_distr[k] = (np.mean(values), np.std(values), np.max(values), np.min(values))
    return diff_distr


if __name__ == "__main__":
    pass
