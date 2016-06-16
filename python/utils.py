import sys
import csv
import os
import errno
from collections import defaultdict, Counter
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import itertools
import numpy as np
import re
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
import jellyfish

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
MEDIAL_TEMPORAL = [17,53,1006,2006,1016,2016,18,54,1035,2035]
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
         'MEDIAL_TEMPORAL': MEDIAL_TEMPORAL,
         'OCCIPITAL': OCCIPITAL_LOBE,
         'MEDIAL_OCCIPITAL': MEDIAL_OCCIPITAL,
         'SENSORY': SENSORY_LOBE,
         'BASALGANGLIA': BASAL_GANGLIA,
         'THALAMUS': THALAMUS,
         'LIMBIC': LIMBIC,
         'CEREBELLUM_GRAY': CEREBG,
         'CEREBELLUM_WHITE': CEREBW,
         'CEREBRAL_WHITE': WHITEMATTER}


# BLACKLISTED REGIONS
REGION_BLACKLIST = [0,30,62,80,81,82,77,251,252,253,254,255,1000,2000,1004,2004,85,24,14,15,72,4,43,75,76]


########## CSV UPDATING INFRASTRUCTURE ##############

def updateDataFrame(old_frame, new_frame, headers=None, after=None, restrict=True, wipe=True):
    '''
    Merges one DF into the other
    1. If header order specified, reorder
    2. Remove any columns by same name in original frame
    3. Merge together (after column, if specified)
    4. If restrict, only allow indices that existed in old_frame
    5. If wipe, wipe all values in the old_frame under headers
    '''
    # Reorder new columns
    if isinstance(headers,list):
        missing_headers = [_ for _ in headers if _ not in new_frame.columns]
        new_frame = pd.concat((new_frame, pd.DataFrame(columns=missing_headers)))
        new_frame = new_frame[headers]
    else:
        headers = list(new_frame.columns)
    # Wipe old values
    if wipe:
        to_wipe = [_ for _ in headers if _ in old_frame]
        for tw in to_wipe:
            old_frame[tw] = np.nan
    # Remove duplicate columns
    new_columns = set(new_frame.columns)
    old_columns = [_ for _ in old_frame.columns if _ not in new_columns]
    old_frame = old_frame[old_columns]
    # Merge
    if after is not None:
        if after not in old_columns:
            raise Exception("After column name not found")
        after_idx = old_columns.index(after)
        merged_columns = old_columns[:after_idx+1] + list(new_frame.columns) + old_columns[after_idx+1:]
    else:
        merged_columns = list(old_frame.columns) + list(new_frame.columns)
    how = 'outer' if not restrict else 'left'
    merged_frame = old_frame.merge(new_frame,left_index=True,right_index=True,how=how)
    merged_frame = merged_frame[merged_columns]
    return merged_frame

def parseSubjectGroups(df, extraction_fn):
    '''
    Takes in df indexed by subject
    Groups the subjects, then feeds each group through extraction_fn
    Which must spit out a Dataframe (one row)
    They get concatenated together

    extraction_fn(group_key, group_rows):
        return pd.DataFrame
    '''
    if not df.index.name:
        df.index.name = 'SID'
    index_name = df.index.name
    reset_df = df.reset_index()
    parsed_df = pd.concat((extraction_fn(i,_) for i,_ in reset_df.groupby(index_name)),axis=0)
    return parsed_df

def groupLongPivot(group_rows, index_column, val_column, key_prefix):
    df = group_rows[[index_column,val_column]]
    df['key'] = ['%s%s' % (key_prefix,i+1) for i in range(len(group_rows.index))]
    out_df = df.pivot_table(values=val_column,
                            index=index_column,
                            columns='key',
                            aggfunc=lambda x: x.iloc[0])

    return out_df

def groupSlope(group_rows, date_key, val_key, cutoff_date=None):
    subset = group_rows[[date_key, val_key]]
    subset.sort_values(by=date_key,inplace=True)
    if cutoff_date is None:
        bl_row = group_rows.iloc[0]
        cutoff_date = bl_row[date_key]
    subset = subset[subset[date_key] >= cutoff_date]
    subset.dropna(inplace=True)
    if len(subset.index) >= 2:
        subset[date_key] = subset[date_key].apply(lambda x: yrDiff(x,cutoff_date))
        grp_slope = slope(zip(subset[date_key],subset[val_key]))
    else:
        grp_slope = np.nan
    return grp_slope

def groupClosest(group_rows, date_key, val_key, comp_dates, day_limit=365):
    df = group_rows.copy()
    results = []
    for i,cd in enumerate(comp_dates):
        if isnan(cd):
            results.append(np.nan)
        else:
            closest_key = 'closest_%s' % i
            df[closest_key] = df[date_key].apply(lambda x: abs(x-cd).days)
            closest_row = df.sort_values(by=closest_key).iloc[0]
            closest_val = closest_row[val_key] if closest_row[closest_key] <= day_limit else np.nan
            results.append(closest_val)
    return results

def yrDiff(date1,date2):
    if isnan(date1) or isnan(date2):
        return np.nan
    else:
        return (date1-date2).days/365.25

def isnan(arg):
    try:
        if pd.isnull(arg):
            return True
    except Exception as e:
        pass
    return False

def getClosestToScans(points, bl_scan, scan_2, scan_3, day_limit=550, date_key='EXAMDATE'):
    used = set()
    # match up with av45 scans
    scan_bl_closest = scan_2_closest = scan_3_closest = None
    if bl_scan is not None:
        eligible = [p for p in points if p[date_key] not in used and (abs(bl_scan - p[date_key]).days <= day_limit)]
        cand = sorted(eligible, key=lambda x: abs(bl_scan-x[date_key]))
        if len(cand) > 0:
            scan_bl_closest = cand[0]
            used.add(scan_bl_closest[date_key])
    if scan_2 is not None:
        eligible = [p for p in points if p[date_key] not in used and (abs(scan_2 - p[date_key]).days <= day_limit)]
        cand = sorted(eligible, key=lambda x: abs(scan_2-x[date_key]))
        if len(cand) > 0:
            scan_2_closest = cand[0]
            used.add(scan_2_closest[date_key])
    if scan_3 is not None:
        eligible = [p for p in points if p[date_key] not in used and (abs(scan_3 - p[date_key]).days <= day_limit)]
        cand = sorted(eligible, key=lambda x: abs(scan_3-x[date_key]))
        if len(cand) > 0:
            scan_3_closest = cand[0]
            used.add(scan_3_closest[date_key])
    return scan_bl_closest, scan_2_closest, scan_3_closest

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


def df_count(df, time_cols, value_cols):
    '''
    Count the number of valid timepoints that would
    go into slope calculations

    Return a Series
    '''
    if len(time_cols) != len(value_cols):
        raise Exception("Time and Value columns not same length")
    timepoints = len(time_cols)

    def calc_count(row):
        times = [row[tc] for tc in time_cols]
        values = [row[vc] for vc in value_cols]
        points = [(t,v) for  t,v in zip(times,values) if not isnan(t) and not isnan(v)]
        return len(points)

    counts = df.apply(calc_count, axis=1)
    return counts


def df_slope(df, time_cols, value_cols, take_diff=False, exact=False):
    '''
    Calculate slope for each row in the dataframe
    Given the time columns and the value columns

    If take_diff, find the year diff between the
    time columns before calculating slope

    If exact, return None for a row if it doesn't
    have values for all timepoints

    Return a Series
    '''
    if len(time_cols) != len(value_cols):
        raise Exception("Time and Value columns not same length")
    timepoints = len(time_cols)

    def calc_slope(row):
        times = [row[tc] for tc in time_cols]
        values = [row[vc] for vc in value_cols]
        points = [(t,v) for  t,v in zip(times,values) if not isnan(t) and not isnan(v)]
        if exact and len(points) != timepoints:
            return None
        if take_diff:
            first_time = points[0][0]
            points = [(yrDiff(t,first_time),v) for t,v in points]
        return slope(points)

    slopes = df.apply(calc_slope, axis=1)
    return slopes

def df_mean(df, weight_cols, value_cols):
    '''
    Calculate weighted means for each row
    given the weight cols and value cols
    '''
    if len(weight_cols) != len(value_cols):
        raise Exception("Weight and Value cols not same length")

    def calc_mean(row):
        weights = [row[wc] for wc in weight_cols]
        values = [row[vc] for vc in value_cols]
        if any(map(isnan,weights)) or any(map(isnan,values)):
            return None
        return weightedMean(zip(weights,values))

    means = df.apply(calc_mean, axis=1)
    return means


def parseDate(date_str, ignore_error=False):
    if date_str == '' or isnan(date_str):
        if ignore_error:
            return None
        else:
            raise Exception("Blank date string")
    if isinstance(date_str,datetime):
        return date_str
    date_str = date_str.replace('--','1')
    formats = ['%Y-%m-%d', '%m/%d/%y', '%m/%d/%Y']
    for f in formats:
        try:
            date = datetime.strptime(date_str,f)
            return date
        except:
            pass
    if ignore_error:
        return None
    else:
        raise Exception("Cannot parse date: %s" % date_str)

def parseOrFindDate(df, date_key, registry=None):
    if date_key not in df.columns:
        df[date_key] = None
    df_null = df[df[date_key].isnull()]
    if registry is not None:
        df_null.loc[:,date_key] = df_null.apply(lambda x: findVisitDate(registry, x.name, [x.get('VISCODE',''),x.get('VISCODE2','')]), axis=1)
    df_nonnull = df[~df[date_key].isnull()]
    df_nonnull.loc[:,date_key] = df_nonnull.loc[:,date_key].apply(parseDate)
    df = pd.concat((df_nonnull,df_null))
    return df

def findVisitDate(registry, subj, viscodes):
    vcs = [_.lower().strip() for _ in viscodes if not isnan(_)]
    vcs = set([_ for _ in vcs if _ != ''])
    date = None
    if isinstance(registry, pd.DataFrame):
        if subj in registry.index:
            subj_registry = registry.loc[subj]
            rows = subj_registry.iterrows()
        else:
            return date
    elif isinstance(registry, dict):
        subj_registry = registry.get(subj,[])
        rows = enumerate(subj_registry)
    else:
        raise Exception("Invalid registry data type")
    for (vc,vc2), reg_row in rows:
        row_visits = [vc.lower().strip(), vc2.lower().strip()]
        row_visits = set([_ for _ in row_visits if _ != ''])
        if len(vcs & row_visits) > 0:
            date = reg_row['EXAMDATE']
            break
    return date

def convertToSubjDict(df, sort_by=None, extract=None, singular=False):
    by_subj = {}
    index_name = df.index.name
    for group_key, group_rows in df.reset_index().groupby(index_name):
        if sort_by in group_rows.columns:
            group_rows.sort_values(by=sort_by, inplace=True)
        if extract and extract in group_rows.columns:
            result = [_.get(extract) for _ in group_rows.to_dict(orient='records')]
        else:
            result = group_rows.to_dict(orient='records')
        if singular:
            by_subj[group_key] = result[0]
        else:
            by_subj[group_key] = result
    return by_subj

def dumpDFtoCSV(df,output_path,decimal_places=2):
    assert isinstance(decimal_places,int)
    float_format='%.' + str(decimal_places) + 'f'
    df.to_csv(output_path,
              sep=',',
              na_rep='',
              float_format=float_format,
              index=True,
              date_format='%m/%d/%y')

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

def extractLongitudinalFields(datapoints, key, prefix):
    data = {}
    for i, point in enumerate(datapoints):
        data['%s%s' % (prefix,i+1)] = point.get(key)
    return data

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

def appendCSV(csv1, csv2):
    '''
    appends csv2 rows to csv1 (have to have same headers)
    '''
    h1, r1 = parseCSV(csv1)
    h2, r2 = parseCSV(csv2)
    all_rows = r1 + r2
    all_headers = h1
    dumpCSV(csv1, all_headers, all_rows)


######################################################

################ SHELL UTILS ###################

def removeFile(filename):
    try:
        os.remove(filename)
    except OSError:
        print "File to remove does not exist"

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

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

######################################################

################# STATS/GRAPH UTILS ###################

def asymIndex(left, right):
    return (left-right)/np.mean([left,right])

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

def regCoeffZTest(B1, B2, se1, se2):
    pooledSE = np.sqrt(se1 + se2)
    z = (B1 - B2) / pooledSE
    pval = 2*(1 - norm.cdf(abs(z)))
    return (z, pval)


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


######################################################



def importRoussetCSV(rousset_csv, translate_threshold=1.11, as_df=False):
    '''
    If include_threshold, translates to the new PVC threshold
    '''
    df = pd.read_csv(rousset_csv)

    slope_points = []
    if as_df:
        by_subj = []
    else:
        by_subj = {}
    for subject, subj_rows in df.groupby('subject'):
        tp_results = {}
        if isinstance(subject,str):
            rid = int(subject.split('-')[-1])
        else:
            rid = int(subject)
        for timepoint, tp_rows in subj_rows.groupby('timepoint'):
            pvc_df = tp_rows[['name','pvcval','nonpvcval','groupsize']]
            pvc_dict = pvc_df.set_index('name').to_dict(orient='index')
            pvc_suvr = float(pvc_dict['COMPOSITE']['pvcval'])/float(pvc_dict['WHOLECEREB']['pvcval'])
            nonpvc_suvr = float(pvc_dict['COMPOSITE']['nonpvcval'])/float(pvc_dict['WHOLECEREB']['nonpvcval'])
            slope_points.append((nonpvc_suvr,pvc_suvr))
            tp_results[timepoint] = pvc_dict
        if as_df:
            for tp, values in tp_results.iteritems():
                pvc_values = {k:v['pvcval'] for k,v in values.iteritems()}
                size_values = {'%s_SIZE' % k : v['groupsize'] for k,v in values.iteritems()}
                pvc_values.update(size_values)
                pvc_values['RID'] = rid
                pvc_values['TP'] = tp
                by_subj.append(pvc_values)
        else:
            by_subj[rid] = tp_results

    if as_df:
        by_subj = pd.DataFrame(by_subj).set_index('RID')

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



def extractDiagnosesFromMasterData(master_df, as_df=True):
    df = master_df[['Init_Diagnosis']]
    df.columns = ['diag']
    df.dropna(inplace=True)
    df['diag'] = df['diag'].apply(lambda x: x.strip())
    if as_df:
        return df
    else:
        return convertToSubjDict(df,extract='diag', singular=True)


def parseCSV(file_path, delimiter=',', use_second_line=False, as_df=False):
    if use_second_line:
        data = pd.read_csv(file_path, sep=delimiter, low_memory=False, header=[0,1])
        data.columns = data.columns.get_level_values(1)
    else:
        data = pd.read_csv(file_path, sep=delimiter, low_memory=False)
    if as_df:
        return data
    else:
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


def importRegistry(registry_file, include_all=False, as_df=False):
    df = pd.read_csv(registry_file,low_memory=False)
    df = df[['RID', 'VISCODE', 'VISCODE2', 'EXAMDATE']]
    df['VISCODE2'] = df['VISCODE2'].fillna('')
    df['VISCODE'] = df['VISCODE'].fillna('')
    df.loc[:,'VISCODE'] = df.loc[:,'VISCODE'].apply(lambda x: x.lower().strip())
    df.loc[:,'VISCODE2'] = df.loc[:,'VISCODE2'].apply(lambda x: x.lower().strip())
    df['VISCODE'] = df['VISCODE'].apply(lambda x: str(x))
    df['VISCODE2'] = df['VISCODE2'].apply(lambda x: str(x))
    if not include_all:
        df = df[~df.VISCODE.str.contains('sc')] # filter out screening visits

    df.set_index(['RID','VISCODE','VISCODE2'],inplace=True)
    df = parseOrFindDate(df, 'EXAMDATE')
    df.dropna(subset=['EXAMDATE'],inplace=True)

    if as_df:
        df.sortlevel(inplace=True)
        return df
    else:
        return convertToSubjDict(df, sort_by='EXAMDATE')

def importDODRegistry(dod_registry_file, as_df=False):
    df = pd.read_csv(dod_registry_file)

    all_columns = ['SCRNO','RID','VISCODE','EXAMDATE']
    df = df[all_columns]
    df.dropna(inplace=True)
    df.set_index('SCRNO',inplace=True)

    #df = df[~df['VISCODE'].str.startswith('sc')]
    df.loc[:,'EXAMDATE'] = df.loc[:,'EXAMDATE'].apply(parseDate)

    if as_df:
        return df
    else:
        return convertToSubjDict(df, sort_by='EXAMDATE')

def importDODAntidep(backmeds_file, registry, as_df=False):
    '''
    FOR DOD BACKMEDS.CSV FILE ONLY
    '''

    SSRI = set(['citalopram',
                'celexa',
                'escitalopram',
                'lexapro',
                'fluoxetine',
                'prozac',
                'fluvoxamine',
                'luvox',
                'paroxetine',
                'paxil',
                'sertraline',
                'zoloft'])
    SNRI = set(['duloxetine',
                'cymbalta',
                'desvenlafaxine',
                'pristiq',
                'venlafaxine',
                'effexor',
                'milnacipran',
                'savella',
                'levomilnacipran',
                'fetzima'])

    def findSSRI(med_str):
        med_str = re.sub('[^a-zA-Z]+', ' ', med_str)
        meds = [_.strip().lower() for _ in med_str.split()]
        min_distances = []
        for m in meds:
            min_dist_ssri = min([jellyfish.damerau_levenshtein_distance(unicode(_), unicode(m)) for _ in SSRI])
            min_dist_snri = min([jellyfish.damerau_levenshtein_distance(unicode(_), unicode(m)) for _ in SNRI])
            if min_dist_ssri <= min_dist_snri:
                min_distances.append(min_dist_ssri)
        if len(min_distances) == 0:
            return False
        closest = min(min_distances)
        if closest <= 2:
            return True
        return False

    backmeds_df = pd.read_csv(backmeds_file)
    backmeds_df = backmeds_df[backmeds_df.VISCODE.str.lower() != 'sc']
    backmeds_df['ANTIDEP_USE'] = backmeds_df['KEYMED'].str.contains('6')
    backmeds_df['WHICH_ANTIDEP'] = backmeds_df['ANTDEPDES']
    backmeds_df.loc[backmeds_df.WHICH_ANTIDEP == '-4','WHICH_ANTIDEP'] = ''
    backmeds_df.loc[:,'SSRI'] = backmeds_df.loc[:,'WHICH_ANTIDEP'].apply(findSSRI)
    backmeds_df['WHICH_ANTIDEP'] = backmeds_df['WHICH_ANTIDEP'].apply(lambda x: x.replace(',',';'))

    # keep relevant columns
    backmeds_df = backmeds_df[['SCRNO','VISCODE','ANTIDEP_USE','WHICH_ANTIDEP','SSRI']]
    backmeds_df.set_index('SCRNO',inplace=True)

    backmeds_df = parseOrFindDate(backmeds_df, 'EXAMDATE', registry=registry)
    backmeds_df.dropna(subset=['EXAMDATE'],inplace=True)

    if as_df:
        backmeds_df['EXAMDATE'] = pd.to_datetime(backmeds_df['EXAMDATE'])
        return backmeds_df
    else:
        return convertToSubjDict(backmeds_df, sort_by='EXAMDATE')


def importCAPS(caps_curr_file, caps_lifetime_file, as_df=False):
    curr_df = pd.read_csv(caps_curr_file)
    curr_df = curr_df[['SCRNO','VISCODE','EXAMDATE','CAPSSCORE']]
    curr_df.rename(columns={'CAPSSCORE': 'curr', 'EXAMDATE':'EXAMDATE_curr'}, inplace=True)
    curr_df = curr_df.groupby(['SCRNO','VISCODE']).aggregate(lambda x: x.iloc[0]).reset_index()

    life_df = pd.read_csv(caps_lifetime_file)
    life_df = life_df[['SCRNO','VISCODE','EXAMDATE','CAPSSCORE']]
    life_df.rename(columns={'CAPSSCORE': 'life', 'EXAMDATE':'EXAMDATE_life'}, inplace=True)
    life_df = life_df.groupby(['SCRNO','VISCODE']).aggregate(lambda x: x.iloc[0]).reset_index()

    df = curr_df.merge(life_df,left_on=['SCRNO','VISCODE'],right_on=['SCRNO','VISCODE'],how='outer')
    df.set_index('SCRNO',inplace=True)
    if as_df:
        return df
    else:
        return convertToSubjDict(df)

def deriveDiagCode(diag_row):
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
    diag_row = diag_row.fillna('')
    change = diag_row.get('DXCHANGE','')
    change = int(change) if change != '' else None
    current = diag_row.get('DXCURREN','')
    current = int(current) if current != '' else None
    conv = diag_row.get('DXCONV','')
    conv = int(conv) if conv != '' else None
    conv_type = diag_row.get('DXCONTYP','')
    conv_type = int(conv_type) if conv_type != '' and int(conv_type) >= 0 else None
    rev_type = diag_row.get('DXREV','')
    rev_type = int(rev_type) if rev_type != '' and int(rev_type) >= 0 else None

    if change is None:
        if conv == 1:
            # conversion
            if conv_type == 1:
                change = 4
            elif conv_type == 2:
                change = 6
            elif conv_type == 3:
                change = 5
        elif conv == 2:
            # reversion
            if rev_type == 1:
                change = 7
            elif rev_type == 2:
                change = 8
            elif rev_type == 3:
                change = 9
        elif conv == 0:
            # stable
            if current == 1:
                change = 1
            elif current == 2:
                change = 2
            elif current == 3:
                change = 3

    return change


def importADNIDiagnosis(diag_file, registry=None, as_df=False):
    df = pd.read_csv(diag_file)
    if 'SCRNO' in df.columns:
        df.set_index('SCRNO',inplace=True)
    else:
        df.set_index('RID',inplace=True)

    all_columns = ['VISCODE','VISCODE2','EXAMDATE','DXCHANGE','DXCURREN','DXCONV','DXCONTYP','DXREV']
    all_columns = [_ for _ in all_columns if _ in df.columns]
    df = df[all_columns]

    df = parseOrFindDate(df, 'EXAMDATE', registry=registry)
    df['change'] = df.apply(lambda x: deriveDiagCode(x), axis=1)
    df.dropna(subset=['change','EXAMDATE'],inplace=True)

    conversion_dict = {1: 'N',
                       2: 'MCI',
                       3: 'AD',
                       4: 'MCI',
                       5: 'AD',
                       6: 'AD',
                       7: 'N',
                       8: 'MCI',
                       9: 'N'}
    df['diag'] = df['change'].apply(lambda x: conversion_dict[x])

    if as_df:
        df['EXAMDATE'] = pd.to_datetime(df['EXAMDATE'])
        return df
    else:
        return convertToSubjDict(df, sort_by='EXAMDATE')


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

def importAPOE(apoe_file, as_df=False):
    df = pd.read_csv(apoe_file)
    df.set_index('SCRNO',inplace=True)
    df = df[['APGEN1','APGEN2']]
    if as_df:
        return df
    else:
        return convertToSubjDict(df)
    return apoe_by_subj

def importFDG(fdg_file, registry, as_df=False):
    df = pd.read_csv(fdg_file)
    df.set_index('RID',inplace=True)
    columns = ['VISCODE','VISCODE2','EXAMDATE','ROINAME','MEAN']
    df = df[columns]
    df = parseOrFindDate(df, 'EXAMDATE', registry=registry)
    df.dropna(subset=['EXAMDATE'],inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df)

# UNUSED
'''
def importExtractedFDG(fdg_file, registry, as_df=False):
    headers, lines = parseCSV(fdg_file, delimiter='\t')
    fdg_rows = []

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
'''

def importDemog(demog_file, as_df=False):
    df = pd.read_csv(demog_file)
    if 'SCRNO' in df.columns:
        df.set_index('SCRNO',inplace=True)
    else:
        df.set_index('RID',inplace=True)

    columns = ['PTGENDER','PTAGE','PTMARRY','PTEDUCAT','PTDOB','PTDOBMM','PTDOBYY']
    columns = [_ for _ in columns if _ in df.columns]
    df = df[columns]

    def createDateFromParts(year, month):
        try:
            return datetime(year=int(year), month=int(month), day=1)
        except Exception as e:
            return np.nan

    if 'PTDOBMM' in columns and 'PTDOBYY' in columns:
        df['PTDOB'] = df.apply(lambda x: createDateFromParts(x['PTDOBYY'],x['PTDOBMM']), axis=1)
        df.drop(['PTDOBMM','PTDOBYY'],axis=1,inplace=True)
    elif 'PTDOB' in columns:
        df['PTDOB'] = df.PTDOB.apply(lambda x: parseDate(x, ignore_error=True))

    df.dropna(subset=['PTDOB'],inplace=True)
    df['PTDOB'] = pd.to_datetime(df['PTDOB'])
    index_name = df.index.name
    df.reset_index(inplace=True)
    df = df.groupby(index_name).max()

    if as_df:
        return df
    else:
        return convertToSubjDict(df)


def importMMSE(mmse_file, registry=None, as_df=False):
    df = pd.read_csv(mmse_file)
    if 'SCRNO' in df.columns:
        df.set_index('SCRNO',inplace=True)
    else:
        df.set_index('RID',inplace=True)

    columns = ['EXAMDATE','VISCODE','VISCODE2','MMSCORE']
    columns = [_ for _ in columns if _ in df.columns]
    df = df[columns]
    df = parseOrFindDate(df, 'EXAMDATE', registry=registry)
    df.dropna(subset=['EXAMDATE'],inplace=True)
    if as_df:
        return df
    else:
        return convertToSubjDict(df, sort_by='EXAMDATE')


def importAVLT(avlt_file, registry=None, as_df=False):
    df = pd.read_csv(avlt_file,low_memory=False)
    if 'SCRNO' in df.columns:
        df.set_index('SCRNO',inplace=True)
    else:
        df.set_index('RID',inplace=True)

    tp_columns = ['EXAMDATE', 'VISCODE', 'VISCODE2']
    avlt_columns = ['AVTOT%s' % _ for _ in range(1,6)]
    all_columns = [_ for _ in tp_columns+avlt_columns if _ in df]
    df = df[all_columns]
    df = parseOrFindDate(df, 'EXAMDATE', registry=registry)
    df = df[(df[avlt_columns] > 0).apply(all,axis=1)]
    df.dropna(inplace=True)
    df['TOTS'] = df[avlt_columns].sum(axis=1)

    if as_df:
        return df
    else:
        return convertToSubjDict(df,sort_by='EXAMDATE')

def importADASCog(adni1_file, adnigo2_file, registry, as_df=False):
    if adni1_file and os.path.isfile(adni1_file):
        df1 = pd.read_csv(adni1_file)
        df1.set_index('RID',inplace=True)
        all_columns = ['VISCODE','EXAMDATE','TOTAL11']
        df1 = df1[all_columns]
        df1['VISCODE2'] = None
        df1['VISCODE'] = df1['VISCODE'].apply(lambda x: x.lower().strip())
        df1['TOTSCORE'] = df1['TOTAL11']

    else:
        df1 = pd.DataFrame()
    if adnigo2_file and os.path.isfile(adnigo2_file):
        df2 = pd.read_csv(adnigo2_file)
        if 'SCRNO' in df2:
            df2.set_index('SCRNO',inplace=True)
        else:
            df2.set_index('RID',inplace=True)
        all_columns = ['VISCODE','VISCODE2','TOTSCORE']
        all_columns = [_ for _ in all_columns if _ in df2.columns]
        df2 = df2[all_columns]
    else:
        df2 = pd.DataFrame()

    df = pd.concat((df1,df2))
    df = parseOrFindDate(df, 'EXAMDATE', registry=registry)
    df.dropna(subset=['EXAMDATE','TOTSCORE'], inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df,sort_by='EXAMDATE')


def importTBMSyn(tbm_file, as_df=False):
    df = pd.read_csv(tbm_file, low_memory=False)
    all_columns = ['RID','VISCODE','VISCODE2','EXAMDATE','EXAMDATEBL','TBMSYNSCOR']
    df = df[all_columns]
    df.set_index('RID',inplace=True)
    df.dropna(inplace=True)
    df.loc[:,'EXAMDATE'] = df.loc[:,'EXAMDATE'].apply(parseDate)
    df.loc[:,'EXAMDATEBL'] = df.loc[:,'EXAMDATEBL'].apply(parseDate)

    if as_df:
        return df
    else:
        return convertToSubjDict(df,sort_by='EXAMDATE')

def importWMH(wmh_file, as_df=False):
    df = pd.read_csv(wmh_file)
    if 'SCRNO' in df:
        df.set_index('SCRNO',inplace=True)
    else:
        df.set_index('RID',inplace=True)

    all_columns = ['VISCODE','VISCODE2','EXAMDATE','WHITMATHYP','ICV']
    all_columns = [_ for _ in all_columns if _ in df.columns]
    df = df[all_columns]

    df.loc[:,'EXAMDATE'] = df.loc[:,'EXAMDATE'].apply(parseDate)
    df.loc[:,'wmh_percent'] = df['WHITMATHYP'] / df['ICV']
    df['wmh'] = df['WHITMATHYP']
    df['icv'] = df['ICV']

    if as_df:
        return df
    else:
        return convertToSubjDict(df,sort_by='EXAMDATE')

def importUW(uw_file, registry=None, as_df=False):
    df = pd.read_csv(uw_file, low_memory=False)
    df.set_index('RID',inplace=True)
    all_columns = ['VISCODE2','EXAMDATE','ADNI_MEM','ADNI_EF']
    df = df[all_columns]

    df = parseOrFindDate(df, 'EXAMDATE', registry=registry)
    df.dropna(subset=['EXAMDATE'], inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df,sort_by='EXAMDATE')

def importGD(gd_file, registry=None, as_df=False):
    df = pd.read_csv(gd_file)
    if 'SCRNO' in df:
        df.set_index('SCRNO',inplace=True)
    else:
        df.set_index('RID',inplace=True)

    all_columns = ['VISCODE','VISCODE2','GDTOTAL','EXAMDATE']
    all_columns = [_ for _ in all_columns if _ in df.columns]
    df = df[all_columns]

    df = parseOrFindDate(df, 'EXAMDATE', registry=registry)
    df.dropna(subset=['EXAMDATE'], inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df,sort_by='EXAMDATE')


def importCSF(csf_files, registry=None, as_df=False):
    dfs = []
    id_names = set()
    for csffile in csf_files:
        df = pd.read_csv(csffile)
        if 'SCRNO' in df:
            df.set_index('SCRNO',inplace=True)
            id_names.add('SCRNO')
        else:
            df.set_index('RID',inplace=True)
            id_names.add('RID')
        df.loc[:,'RUNDATE'] = df.loc[:,'RUNDATE'].apply(parseDate)
        df = parseOrFindDate(df, 'EXAMDATE', registry=registry)
        df.dropna(subset=['EXAMDATE','RUNDATE'],inplace=True)
        dfs.append(df)

    if len(id_names) != 1:
        raise Exception('Multiple ID columns used in CSF files')

    id_name = list(id_names)[0]
    df_all = pd.concat(dfs)
    df_all.index.name = id_name
    df_all.reset_index(inplace=True)

    idx = df_all['RUNDATE'] == df_all.groupby(by=[id_name,'EXAMDATE'])['RUNDATE'].transform(max)
    df_all = df_all[idx]
    df_all.set_index(id_name,inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df,sort_by='EXAMDATE')

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

    df.dropna(subset=[subj_col, date_col], inplace=True)
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

def importPetMETA(pet_meta_file, tracer='AV45', as_df=False):
    df = pd.read_csv(pet_meta_file)
    seq_str = '%s Coreg, Avg, Std Img and Vox Siz, Uniform Resolution' % tracer
    df.loc[:,'Sequence'] = df.loc[:,'Sequence'].apply(lambda x: x.strip())
    df = df[df['Sequence'].str.match(seq_str)]
    df['Subject'] = df.loc[:,'Subject'].apply(lambda x: int(x.split('_')[-1].strip()))
    df['Scan Date'] = df.loc[:,'Scan Date'].apply(parseDate)

    if as_df:
        df = df[['Subject','Scan Date']]
        def add_tp(rows):
            rows.sort_values('Scan Date', inplace=True)
            rows['TP'] = range(1,len(rows.index)+1)
            return rows
        df = df.groupby('Subject').apply(add_tp)
        df = df.pivot('Subject','TP','Scan Date')
        tps = [1,2,3,4]
        for tp in tps:
            if tp not in df.columns:
                df[tp] = np.nan
        return df
    else:
        scans = {pid: sorted(list(set(rows['Scan Date']))) for pid, rows in df.groupby('Subject')}
        return scans

def importARM(arm_file, as_df=False):
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

    df = pd.read_csv(arm_file)
    df.set_index('RID',inplace=True)
    all_columns = ['USERDATE','ARM']
    df = df[all_columns]
    df.loc[:,'USERDATE'] = df.loc[:,'USERDATE'].apply(parseDate)
    df.dropna(subset=['ARM'],inplace=True)
    df['STATUS'] = df.loc[:,'ARM'].apply(lambda x: translation[x])

    if as_df:
        return df
    else:
        return convertToSubjDict(df,sort_by='USERDATE')


def importDODEligibility(elig_file, as_df=False):
    df = pd.read_csv(elig_file)
    df.set_index('SCRNO',inplace=True)
    df = df[['COHORT']]
    df.dropna(inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df)


def importMedicalHistory(mhist_file, as_df=False):
    mhist_df = pd.read_csv(mhist_file)
    mhist_df.loc[:,'MHDESC'] = mhist_df.MHDESC.apply(lambda x: x.lower())
    mhist_df.loc[:,'smoking'] = mhist_df.MHDESC.str.contains('smok')
    mhist_df.loc[:,'diabetes'] = mhist_df.MHDESC.str.contains('diabete')
    mhist_df.loc[:,'hyperlipidemia'] = mhist_df.MHDESC.str.contains('hyperlipidemia')
    mhist_df = mhist_df[['RID','smoking','diabetes','hyperlipidemia']]
    df = mhist_df.groupby('RID').aggregate(any)
    if as_df:
        return df
    else:
        return convertToSubjDict(df)

def importASI(asi_file, registry, as_df=False):
    df = pd.read_csv(asi_file)
    if 'SCRNO' in df.columns:
        df.set_index('SCRNO',inplace=True)
    else:
        df.set_index('RID',inplace=True)

    columns = ['EXAMDATE','VISCODE','VISCODE2','ASI1A','ASI2B']
    columns = [_ for _ in columns if _ in df.columns]
    df = df[columns]

    df = parseOrFindDate(df, 'EXAMDATE', registry=registry)
    df.dropna(subset=['EXAMDATE'], inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df, sort_by='EXAMDATE')


def importCDR(cdr_file, registry, as_df=False):
    df = pd.read_csv(cdr_file)
    if 'SCRNO' in df.columns:
        df.set_index('SCRNO',inplace=True)
    else:
        df.set_index('RID',inplace=True)

    columns = ['EXAMDATE','VISCODE','VISCODE2']
    columns += ['CDSOB','CDMEMORY','CDORIENT','CDJUDGE',
                'CDCOMMUN','CDHOME','CDCARE','CDGLOBAL']
    columns = [_ for _ in columns if _ in df.columns]
    df = df[columns]

    df = parseOrFindDate(df, 'EXAMDATE', registry=registry)
    df.dropna(subset=['EXAMDATE'], inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df, sort_by='EXAMDATE')


def importFAQ(faq_file, registry, as_df=False):
    df = pd.read_csv(faq_file)
    if 'SCRNO' in df.columns:
        df.set_index('SCRNO',inplace=True)
    else:
        df.set_index('RID',inplace=True)

    all_columns = ['EXAMDATE','VISCODE','VISCODE2','FAQTOTAL']
    all_columns = [_ for _ in all_columns if _ in df.columns]
    df = df[all_columns]

    df = parseOrFindDate(df, 'EXAMDATE', registry=registry)
    df = df[['EXAMDATE','FAQTOTAL']]
    df.dropna(inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df, sort_by='EXAMDATE')


def importNPI(npi_file, as_df=False):
    df = pd.read_csv(npi_file)
    df = df[['RID','EXAMDATE','NPITOTAL']]
    df.dropna(inplace=True)
    df['EXAMDATE'] = df.loc[:,'EXAMDATE'].apply(parseDate)
    df.set_index('RID',inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df, sort_by='EXAMDATE')


def importDODMRI(mri_file, as_df=False):
    df = pd.read_csv(mri_file)
    df.set_index('SCRNO',inplace=True)
    columns = ['MMCONDCT','EXAMDATE']
    df = df[columns]
    df = df[df.MMCONDCT == 1]
    df['EXAMDATE'] = df.EXAMDATE.apply(lambda x: parseDate(x,ignore_error=True))
    df.dropna(inplace=True)
    df = df[['EXAMDATE']]
    if as_df:
        return df
    else:
        return convertToSubjDict(df, sort_by='EXAMDATE', extract='EXAMDATE')

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


def importMaster(master_file, as_df=False):
    df = pd.read_csv(master_file, low_memory=False, header=[0,1])
    df.columns = df.columns.get_level_values(1)
    df.set_index('RID',inplace=True)
    if as_df:
        return df
    else:
        df.fillna('',inplace=True)
        return convertToSubjDict(df)

def importBSI(bsi_file, include_failed=False, as_df=False):
    df = pd.read_csv(bsi_file)
    df.set_index('RID',inplace=True)
    columns = ['QC_PASS','MRSEQUENCE','KMNDBCBBSI','VISCODE', 'VISCODE2', 'DBCBBSI', 'VBSI', 'HBSI_R', 'HBSI_L']
    df = df[columns]
    df.rename(columns={'DBCBBSI': 'WB_BSI', 'KMNDBCBBSI': 'WB_KNBSI'}, inplace=True)

    df = df[~df['MRSEQUENCE'].str.match('Acc')]
    if not include_failed:
        df = df[df['QC_PASS'] != 0]
    df['EXAMDATE'] = df.EXAMDATE.apply(parseDate)
    df_bl = df[df['WB_KNBSI'].str.match('')]
    df_sub = df[~df['WB_KNBSI'].str.match('')]

    def getBaselineDate(subj, df_bl):
        try:
            return df_bl.loc[subj,'EXAMDATE']
        except:
            return None

    df_bl = df_bl[['EXAMDATE']].drop_duplicates()
    df_sub['BL_EXAMDATE'] = df_sub.apply(lambda x: getBaselineDate(x.name,df_bl), axis=1)
    if as_df:
        return df_sub
    else:
        return convertToSubjDict(df_sub,sort_by='EXAMDATE')


def importAV1451(av1451_file, as_df=False):
    df = pd.read_csv(av1451_file)
    if 'SCRNO' in df:
        df.set_index('SCRNO',inplace=True)
    else:
        df.set_index('RID',inplace=True)
    df.loc[:,'EXAMDATE'] = df.loc[:,'EXAMDATE'].apply(parseDate)
    if as_df:
        return df
    else:
        return convertToSubjDict(df,sort_by='EXAMDATE')

def importAV45(av45_file, as_df=False):
    df = pd.read_csv(av45_file)

    if 'SCRNO' in df:
        df.set_index('SCRNO',inplace=True)
    elif 'RID' in df:
        df.set_index('RID',inplace=True)
    elif 'PID' in df:
        df.set_index('PID',inplace=True)
    else:
        raise Exception("Can't find subject column")

    df['TP_SPECIFIC'] = True
    df.loc[:,'EXAMDATE'] = df.loc[:,'EXAMDATE'].apply(lambda x: parseDate(x,ignore_error=True))
    df.dropna(subset=['EXAMDATE'],inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df,sort_by='EXAMDATE')


def importFSVolumes(vol_file, as_df=False):
    df = pd.read_csv(vol_file)
    df.loc[:,'EXAMDATE'] = df.loc[:,'Date'].apply(parseDate)
    df.loc[:,'RID'] = df.loc[:,'Subject'].apply(lambda x: int(x.split('-')[-1]))
    df['HCV'] = df['Left-Hippocampus'] + df['Right-Hippocampus']
    df['ICV'] = df['EstimatedTotalIntraCranialVol']
    df['HC/ICV'] = df['HCV'] / df['ICV']

    df = df[['RID','EXAMDATE','HCV','ICV','HC/ICV']]
    df.dropna(inplace=True)
    df.set_index('RID',inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df,sort_by='EXAMDATE')


def importUCSFFreesurfer(in_file, mprage_file, version='', include_failed=False, as_df=False):
    mprage_df = pd.read_csv(mprage_file)
    mprage_df.set_index('ImageUID',inplace=True)

    def getMagStrength(mprage_df, imageuid):
        try:
            return mprage_df.loc[imageuid,'MagStrength']
        except Exception as e:
            print e
            return np.nan

    # pass qc
    df = pd.read_csv(in_file)
    if not include_failed:
        df = df[df['OVERALLQC'].isin(['Pass','Partial'])]
        df = df[df['TEMPQC'].str.match('Pass')]

    # filter by image type
    if 'IMAGETYPE' in df.columns:
        df = df[df['IMAGETYPE']=='Non-Accelerated T1']

    # get relevant columns
    # right hippocampus volume: ST88SV
    # left hippocampus volume: ST29SV
    # intracranial volume: ST10CV
    columns = ['RID','VISCODE','VISCODE2','EXAMDATE','ST29SV','ST88SV','ST10CV','IMAGEUID']
    columns = [_ for _ in columns if _ in df.columns]
    df = df[columns]
    df['HCV'] = df[['ST88SV','ST29SV']].sum(axis=1)
    df['IMAGEUID'] = df['IMAGEUID'].apply(lambda x: getMagStrength(mprage_df, x))
    df['version'] = version
    df.rename(columns={'ST29SV': 'LEFT_HCV', 'ST88SV': 'RIGHT_HCV', 'ST10CV': 'ICV', 'IMAGEUID': 'FLDSTRENG'},inplace=True)

    df.set_index('RID',inplace=True)
    df = parseOrFindDate(df, 'EXAMDATE', registry=None)
    df.dropna(subset=['EXAMDATE'],inplace=True)

    if as_df:
        return df
    else:
        return convertToSubjDict(df,sort_by='EXAMDATE')


# NOT USED
'''
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
'''
