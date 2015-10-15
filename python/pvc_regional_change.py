'''
Analysis of regional longitudinal data from processing streams:
- Non timepoint specific
- Timepoint specific
- PVC (Rousset, group 4, timepoint specific)
'''


from utils import *
import matplotlib.pyplot as plt
import copy
import numpy as np
import pylab
from scipy.stats import norm, mannwhitneyu, linregress, shapiro, ttest_ind


STABLE = [89,
          272,
          685,
          923,
          4018,
          4032,
          4041,
          4043,
          4050,
          4076,
          4084,
          4148,
          4158,
          4164,
          4173,
          4200,
          4208,
          4213,
          4218,
          4224,
          4369,
          4382,
          4393,
          4410,
          4427,
          4453,
          4505,
          4506,
          4516,
          4545,
          4643,
          4739,
          4762,
          4795,
          4878,
          5023,
          4429]

def extractRegionalValuesRousset(data, pvcval=True):
    if pvcval:
        key = 'pvcval'
    else:
        key = 'nonpvcval'
    data = data['group4']
    wholecereb = float(data['wholecereb'][key])
    bigref = float(data['bigref'][key])
    cingulate_uptake = data['cingulate'][key]/bigref
    cingulate_vol = int(data['cingulate']['size'])
    parietal_uptake = data['parietal'][key]/bigref
    parietal_vol = int(data['parietal']['size'])
    temporal_uptake = data['temporal'][key]/bigref
    temporal_vol = int(data['temporal']['size'])
    frontal_uptake = data['frontal'][key]/bigref
    frontal_vol = int(data['frontal']['size'])
    composite_uptake = data['composite'][key]/bigref
    composite_vol = int(data['composite']['size'])
    uptakes = {'cingulate': cingulate_uptake,
               'parietal': parietal_uptake,
               'temporal': temporal_uptake,
               'frontal': frontal_uptake,
               'composite': composite_uptake}
    sizes = {'cingulate': cingulate_vol,
               'parietal': parietal_vol,
               'temporal': temporal_vol,
               'frontal': frontal_vol,
               'composite': composite_vol}
    return uptakes, sizes

def extractRoussetResiduals(data):
    residuals = {}
    for k,v in data.iteritems():
        if v is None:
            continue
        residuals[k] = v['residuals']['rmse']
    return residuals

def extractRegionalValuesAV45(row, lut):
    wholecereb = float(row['WHOLECEREBELLUM'])
    bigref = float(row['COMPOSITE_REF'])
    cingulate_uptake = float(row['CINGULATE']) / bigref
    parietal_uptake = float(row['PARIETAL']) / bigref
    temporal_uptake = float(row['TEMPORAL']) / bigref
    frontal_uptake = float(row['FRONTAL']) / bigref
    composite_uptake = float(row['COMPOSITE'])/ bigref
    frontal=[1003,1012,1014,1018,1019,1020,1027,1028,1032,2003,2012,2014,2018,2019,2020,2027,2028,2032]
    cingulate=[1002,1010,1023,1026,2002,2010,2023,2026]
    parietal=[1008,1025,1029,1031,2008,2025,2029,2031]
    temporal=[1015,1030,2015,2030]
    cingulate_size_keys = ["%s_SIZE" % lut[_].replace('-','_').upper().strip() for _ in cingulate]
    parietal_size_keys = ["%s_SIZE" % lut[_].replace('-','_').upper().strip() for _ in parietal]
    temporal_size_keys = ["%s_SIZE" % lut[_].replace('-','_').upper().strip() for _ in temporal]
    frontal_size_keys = ["%s_SIZE" % lut[_].replace('-','_').upper().strip() for _ in frontal]
    cingulate_vol = sum([float(row[_]) for _ in cingulate_size_keys])
    parietal_vol = sum([float(row[_]) for _ in parietal_size_keys])
    temporal_vol = sum([float(row[_]) for _ in temporal_size_keys])
    frontal_vol = sum([float(row[_]) for _ in frontal_size_keys])
    composite_vol = sum([cingulate_vol, parietal_vol, temporal_vol, frontal_vol])
    uptakes = {'cingulate': cingulate_uptake,
               'parietal': parietal_uptake,
               'temporal': temporal_uptake,
               'frontal': frontal_uptake,
               'composite': composite_uptake}
    sizes = {'cingulate': cingulate_vol,
               'parietal': parietal_vol,
               'temporal': temporal_vol,
               'frontal': frontal_vol,
               'composite': composite_vol}
    return uptakes, sizes


def parseRoussetOutputs(bl_file, scan2_file, scan3_file, pvcval=True, grouping='group4'):
    group_bl, data_bl_raw = importRoussetResults(bl_file)
    group_scan2, data_scan2_raw = importRoussetResults(scan2_file)
    group_scan3, data_scan3_raw = importRoussetResults(scan3_file)
    data_bl = {}
    data_scan2 = {}
    data_scan3 = {}
    for k,v in data_bl_raw.iteritems():
        data_bl[k] = extractRegionalValuesRousset(v, pvcval=pvcval)
    for k,v in data_scan2_raw.iteritems():
        data_scan2[k] = extractRegionalValuesRousset(v, pvcval=pvcval)
    for k,v in data_scan3_raw.iteritems():
        data_scan3[k] = extractRegionalValuesRousset(v, pvcval=pvcval)
    return data_bl, data_scan2, data_scan3

def parseRoussetResiduals(bl_file, scan2_file, scan3_file):
    group_bl, data_bl_raw = importRoussetResults(bl_file)
    group_scan2, data_scan2_raw = importRoussetResults(scan2_file)
    group_scan3, data_scan3_raw = importRoussetResults(scan3_file)
    data_bl = {}
    data_scan2 = {}
    data_scan3 = {}
    for k,v in data_bl_raw.iteritems():
        data_bl[k] = extractRoussetResiduals(v)
    for k,v in data_scan2_raw.iteritems():
        data_scan2[k] = extractRoussetResiduals(v)
    for k,v in data_scan3_raw.iteritems():
        data_scan3[k] = extractRoussetResiduals(v)
    return data_bl, data_scan2, data_scan3

def parseAV45Output(av45_file, registry_file, lut_file):
    lut_table = importFreesurferLookup(lut_file)
    registry =  importRegistry(registry_file)
    data = importAV45(av45_file,registry=registry)
    data_bl = {}
    data_scan2 = {}
    data_scan3 = {}
    for rid, rows in data.iteritems():
        rows = sorted(rows, key=lambda x: x['EXAMDATE'])
        if len(rows) > 0:
            data_bl[rid] = extractRegionalValuesAV45(rows[0], lut_table)
        if len(rows) > 1:
            data_scan2[rid] = extractRegionalValuesAV45(rows[1], lut_table)
        if len(rows) > 2:
            data_scan3[rid] = extractRegionalValuesAV45(rows[2], lut_table)
    return data_bl, data_scan2, data_scan3

def findRegionalAnnualizedChange(rids, data_bl, data_scan2, data_scan3, master_data, annualize=True):
    # find annualized
    yrs_diff = {}
    vol_frontal_diff = {}
    vol_parietal_diff = {}
    vol_temporal_diff = {}
    vol_cingulate_diff = {}
    vol_composite_diff = {}
    uptake_frontal_diff = {}
    uptake_parietal_diff = {}
    uptake_temporal_diff = {}
    uptake_cingulate_diff = {}
    uptake_composite_diff = {}

    for rid, data in data_scan2.iteritems():
        if rid not in rids:
            continue
        yrs = float(master_data[rid]['AV45_1_2_Diff (Yrs)'])
        if not annualize:
            yrs = 1.0
        yrs_diff[rid] = yrs
        scan2_uptakes, scan2_sizes = data
        bl_uptakes, bl_sizes = data_bl[rid]
        vol_frontal_diff[rid] = (scan2_sizes['frontal'] - bl_sizes['frontal']) / yrs
        vol_parietal_diff[rid] = (scan2_sizes['parietal'] - bl_sizes['parietal']) / yrs
        vol_temporal_diff[rid] = (scan2_sizes['temporal'] - bl_sizes['temporal']) / yrs
        vol_cingulate_diff[rid] = (scan2_sizes['cingulate'] - bl_sizes['cingulate']) / yrs
        vol_composite_diff[rid] = (scan2_sizes['composite'] - bl_sizes['composite']) / yrs
        uptake_frontal_diff[rid] = (scan2_uptakes['frontal'] - bl_uptakes['frontal']) / yrs
        uptake_parietal_diff[rid] = (scan2_uptakes['parietal'] - bl_uptakes['parietal']) / yrs
        uptake_temporal_diff[rid] = (scan2_uptakes['temporal'] - bl_uptakes['temporal']) / yrs
        uptake_cingulate_diff[rid] = (scan2_uptakes['cingulate'] - bl_uptakes['cingulate']) / yrs
        uptake_composite_diff[rid] = (scan2_uptakes['composite'] - bl_uptakes['composite']) / yrs
    for rid, data in data_scan3.iteritems():
        if rid not in rids:
            continue
        yrs = float(master_data[rid]['AV45_1_3_Diff (yrs)'])
        if not annualize:
            yrs = 1.0
        yrs_diff[rid] = yrs
        scan3_uptakes, scan3_sizes = data
        bl_uptakes, bl_sizes = data_bl[rid]
        vol_frontal_diff[rid] = (scan3_sizes['frontal'] - bl_sizes['frontal']) / yrs
        vol_parietal_diff[rid] = (scan3_sizes['parietal'] - bl_sizes['parietal']) / yrs
        vol_temporal_diff[rid] = (scan3_sizes['temporal'] - bl_sizes['temporal']) / yrs
        vol_cingulate_diff[rid] = (scan3_sizes['cingulate'] - bl_sizes['cingulate']) / yrs
        vol_composite_diff[rid] = (scan3_sizes['composite'] - bl_sizes['composite']) / yrs
        uptake_frontal_diff[rid] = (scan3_uptakes['frontal'] - bl_uptakes['frontal']) / yrs
        uptake_parietal_diff[rid] = (scan3_uptakes['parietal'] - bl_uptakes['parietal']) / yrs
        uptake_temporal_diff[rid] = (scan3_uptakes['temporal'] - bl_uptakes['temporal']) / yrs
        uptake_cingulate_diff[rid] = (scan3_uptakes['cingulate'] - bl_uptakes['cingulate']) / yrs
        uptake_composite_diff[rid] = (scan3_uptakes['composite'] - bl_uptakes['composite']) / yrs

    vol_diff = {'frontal': vol_frontal_diff,
                'parietal': vol_parietal_diff,
                'temporal': vol_temporal_diff,
                'cingulate': vol_cingulate_diff,
                'composite': vol_composite_diff}
    uptake_diff = {'frontal': uptake_frontal_diff,
                   'parietal': uptake_parietal_diff,
                   'temporal': uptake_temporal_diff,
                   'cingulate': uptake_cingulate_diff,
                   'composite': uptake_composite_diff}

    return (vol_diff, uptake_diff, yrs_diff)

def percentPlausibleNoNoise(vol_diff, uptake_diff, yrs):
    # calculate percent with decreasing volume
    frontal_decreasing = len([rid for rid, val in vol_diff['frontal'].iteritems() if val <= 0.0]) / float(len(vol_diff['frontal']))
    parietal_decreasing = len([rid for rid, val in vol_diff['parietal'].iteritems() if val <= 0.0]) / float(len(vol_diff['parietal']))
    temporal_decreasing = len([rid for rid, val in vol_diff['temporal'].iteritems() if val <= 0.0]) / float(len(vol_diff['temporal']))
    cingulate_decreasing = len([rid for rid, val in vol_diff['cingulate'].iteritems() if val <= 0.0]) / float(len(vol_diff['cingulate']))
    composite_decreasing = len([rid for rid, val in vol_diff['composite'].iteritems() if val <= 0.0]) / float(len(vol_diff['composite']))
    # calculate percent increasing uptake
    frontal_increasing = len([rid for rid, val in uptake_diff['frontal'].iteritems() if val >= 0.0]) / float(len(uptake_diff['frontal']))
    parietal_increasing = len([rid for rid, val in uptake_diff['parietal'].iteritems() if val >= 0.0]) / float(len(uptake_diff['parietal']))
    temporal_increasing = len([rid for rid, val in uptake_diff['temporal'].iteritems() if val >= 0.0]) / float(len(uptake_diff['temporal']))
    cingulate_increasing = len([rid for rid, val in uptake_diff['cingulate'].iteritems() if val >= 0.0]) / float(len(uptake_diff['cingulate']))
    composite_increasing = len([rid for rid, val in uptake_diff['composite'].iteritems() if val >= 0.0]) / float(len(uptake_diff['composite']))
    print "Decreasing Volume"
    print "frontal: %s" % frontal_decreasing
    print "parietal: %s" % parietal_decreasing
    print "temporal: %s" % temporal_decreasing
    print "cingulate: %s" % cingulate_decreasing
    print "composite: %s" % composite_decreasing
    print "Increasing Uptake"
    print "frontal: %s" % frontal_increasing
    print "parietal: %s" % parietal_increasing
    print "temporal: %s" % temporal_increasing
    print "cingulate: %s" % cingulate_increasing
    print "composite: %s" % composite_increasing

def percentPlausible(vol_diff, uptake_diff, norm_fits, yrs):
    '''
    Determine percentage of dataset that is plausible given two measures:
    - volume only decreases
    - amyloid load only increases

    Also determine noise intervals for the proportions, given the scan-rescan errors for each region

    '''
    # calculate percent with decreasing volume
    frontal_decreasing = len([rid for rid, val in vol_diff['frontal'].iteritems() if val <= 0.0]) / float(len(vol_diff['frontal']))
    parietal_decreasing = len([rid for rid, val in vol_diff['parietal'].iteritems() if val <= 0.0]) / float(len(vol_diff['parietal']))
    temporal_decreasing = len([rid for rid, val in vol_diff['temporal'].iteritems() if val <= 0.0]) / float(len(vol_diff['temporal']))
    cingulate_decreasing = len([rid for rid, val in vol_diff['cingulate'].iteritems() if val <= 0.0]) / float(len(vol_diff['cingulate']))
    composite_decreasing = len([rid for rid, val in vol_diff['composite'].iteritems() if val <= 0.0]) / float(len(vol_diff['composite']))

    # calculate percent with increasing uptake
    # find sampling distr of proportions
    frontal_distr = {rid: (val-norm_fits['frontal'][0], norm_fits['frontal'][1]) for rid, val in uptake_diff['frontal'].iteritems()}
    parietal_distr = {rid: (val-norm_fits['parietal'][0], norm_fits['parietal'][1]) for rid, val in uptake_diff['parietal'].iteritems()}
    temporal_distr = {rid: (val-norm_fits['temporal'][0], norm_fits['temporal'][1]) for rid, val in uptake_diff['temporal'].iteritems()}
    cingulate_distr = {rid: (val-norm_fits['cingulate'][0], norm_fits['cingulate'][1]) for rid, val in uptake_diff['cingulate'].iteritems()}
    composite_distr = {rid: (val-norm_fits['composite'][0], norm_fits['composite'][1]) for rid, val in uptake_diff['composite'].iteritems()}
    frontal_cdf = {rid: 1.0 - norm.cdf(0.0, mean, std) for rid, (mean, std) in frontal_distr.iteritems()}
    parietal_cdf = {rid: 1.0 - norm.cdf(0.0, mean, std) for rid, (mean, std) in parietal_distr.iteritems()}
    temporal_cdf = {rid: 1.0 - norm.cdf(0.0, mean, std) for rid, (mean, std) in temporal_distr.iteritems()}
    cingulate_cdf = {rid: 1.0 - norm.cdf(0.0, mean, std) for rid, (mean, std) in cingulate_distr.iteritems()}
    composite_cdf = {rid: 1.0 - norm.cdf(0.0, mean, std) for rid, (mean, std) in composite_distr.iteritems()}

    frontal_n = float(len(uptake_diff['frontal']))
    frontal_mean = sum(frontal_cdf.values()) / float(len(uptake_diff['frontal']))
    frontal_std = np.sqrt(sum([(1.0-_)*_ for _ in frontal_cdf.values()]) / (frontal_n**2))
    parietal_n = float(len(uptake_diff['parietal']))
    parietal_mean = np.sqrt(sum(parietal_cdf.values()) / float(len(uptake_diff['parietal'])))
    parietal_std = sum([(1.0-_)*_ for _ in parietal_cdf.values()]) / (parietal_n**2)
    temporal_n = float(len(uptake_diff['temporal']))
    temporal_mean = sum(temporal_cdf.values()) / float(len(uptake_diff['temporal']))
    temporal_std = np.sqrt(sum([(1.0-_)*_ for _ in temporal_cdf.values()]) / (temporal_n**2))
    cingulate_n = float(len(uptake_diff['cingulate']))
    cingulate_mean = sum(cingulate_cdf.values()) / float(len(uptake_diff['cingulate']))
    cingulate_std = np.sqrt(sum([(1.0-_)*_ for _ in cingulate_cdf.values()]) / (cingulate_n**2))
    composite_n = float(len(uptake_diff['composite']))
    composite_mean = sum(composite_cdf.values()) / float(len(uptake_diff['composite']))
    composite_std = np.sqrt(sum([(1.0-_)*_ for _ in composite_cdf.values()]) / (composite_n**2))

    print "Decreasing Volume"
    print "frontal: %s" % frontal_decreasing
    print "parietal: %s" % parietal_decreasing
    print "temporal: %s" % temporal_decreasing
    print "cingulate: %s" % cingulate_decreasing
    print "composite: %s" % composite_decreasing
    print "Increasing Uptake"
    print "frontal: %s +/- %s" % (frontal_mean, frontal_std)
    print "parietal: %s +/- %s" % (parietal_mean, parietal_std)
    print "temporal: %s +/- %s" % (temporal_mean, temporal_std)
    print "cingulate: %s +/- %s" % (cingulate_mean, cingulate_std)
    print "composite: %s +/- %s" % (composite_mean, composite_std)

def plot_regional_volumechange_vs_uptakechange(colors, vol_diff, uptake_diff):
    frontal_x = []
    frontal_y = []
    frontal_c = []

    for rid in vol_diff['frontal']:
        frontal_x.append(vol_diff['frontal'][rid])
        frontal_y.append(uptake_diff['frontal'][rid])
        frontal_c.append(colors[rid])
    parietal_x = []
    parietal_y = []
    parietal_c = []
    for rid in vol_diff['parietal']:
        parietal_x.append(vol_diff['parietal'][rid])
        parietal_y.append(uptake_diff['parietal'][rid])
        parietal_c.append(colors[rid])
    temporal_x = []
    temporal_y = []
    temporal_c = []
    for rid in vol_diff['temporal']:
        temporal_x.append(vol_diff['temporal'][rid])
        temporal_y.append(uptake_diff['temporal'][rid])
        temporal_c.append(colors[rid])
    cingulate_x = []
    cingulate_y = []
    cingulate_c = []
    for rid in vol_diff['cingulate']:
        cingulate_x.append(vol_diff['cingulate'][rid])
        cingulate_y.append(uptake_diff['cingulate'][rid])
        cingulate_c.append(colors[rid])

    plt.figure(1)
    plt.scatter(frontal_x, frontal_y, c=frontal_c)
    plt.xlabel('Annualized Frontal Size Change')
    plt.ylabel('Annualized Frontal Uptake Change (post PVC)')
    plt.figure(2)
    plt.scatter(parietal_x, parietal_y, c=parietal_c)
    plt.xlabel('Annualized Parietal Size Change')
    plt.ylabel('Annualized Parietal Uptake Change (post PVC)')
    plt.figure(3)
    plt.scatter(temporal_x, temporal_y, c=temporal_c)
    plt.xlabel('Annualized Temporal Size Change')
    plt.ylabel('Annualized Temporal Uptake Change (post PVC)')
    plt.figure(4)
    plt.scatter(cingulate_x, cingulate_y, c=cingulate_c)
    plt.xlabel('Annualized Cingulate Size Change')
    plt.ylabel('Annualized Cingulate Uptake Change (post PVC)')
    plt.show()

def dumpRegionalDataset(output_file, data_bl, data_scan2, data_scan3, master_data):
    lines = []
    for rid, data in data_bl.iteritems():
        uptakes, sizes = data
        all_data = {'RID': rid,
                    'Timepoint': 'BL',
                    'Yrs_from_BL': 0.0}
        all_data.update({"%s_suvr" % k :v for k,v in uptakes.iteritems()})
        all_data.update({"%s_volume" % k :v for k,v in sizes.iteritems()})
        lines.append(all_data)
    for rid, data in data_scan2.iteritems():
        yrs = float(master_data[rid]['AV45_1_2_Diff (Yrs)'])
        uptakes, sizes = data
        all_data = {'RID': rid,
                    'Timepoint': 'Scan2',
                    'Yrs_from_BL': yrs}
        all_data.update({"%s_suvr" % k :v for k,v in uptakes.iteritems()})
        all_data.update({"%s_volume" % k :v for k,v in sizes.iteritems()})
        lines.append(all_data)
    for rid, data in data_scan3.iteritems():
        yrs = float(master_data[rid]['AV45_1_3_Diff (yrs)'])
        uptakes, sizes = data
        all_data = {'RID': rid,
                    'Timepoint': 'Scan3',
                    'Yrs_from_BL': yrs}
        all_data.update({"%s_suvr" % k :v for k,v in uptakes.iteritems()})
        all_data.update({"%s_volume" % k :v for k,v in sizes.iteritems()})
        lines.append(all_data)
    df = pd.DataFrame(lines)
    df.to_csv(output_file,index=False)

def fitNormalToUptakeChange(uptake_diff):
    '''
    remove positive data points, and mirror negative datapoints above the origin to get test-retest error distr
    '''
    normfits = {}
    for k,v in uptake_diff.iteritems():
        values = v.values()
        '''
        neg_values = [_ for _ in values if _ <= 0]
        mirrored_values = [-_ for _ in neg_values if _ < 0]
        all_values = neg_values + mirrored_values
        '''
        all_values = values
        (mu, std) = norm.fit(all_values)
        normfits[k] = (mu, std)
        
    return normfits

def diagGroupEffects(data, diags):
    keys = ['frontal', 'parietal', 'temporal', 'cingulate', 'composite']
    normals = list(set([k for k,v in diags.iteritems() if v in set(['N', 'SMC'])]) & set(data.keys()))
    ads = list(set([k for k,v in diags.iteritems() if v in set(['AD'])]) & set(data.keys()))
    results = {}
    for k in keys:
        normal = [data[_][0][k] for _ in normals]
        ad = [data[_][0][k] for _ in ads]
        u, pvalue = mannwhitneyu(normal, ad, use_continuity=True)
        u_max = len(normal) * len(ad)
        rank_biserial = 1.0 - (2*u/u_max)
        key_result = {'u': u,
                      'p': pvalue,
                      'rank_biserial': rank_biserial}
        results[k] = key_result

    for k,v in results.iteritems():
        print k
        print v


def stratifySubjects(bl_uptake, uptake_diff, norm_fits, yrs, diags, threshold=1.11, graph=False):
    '''
    Put subjects into groups:
    -> amyloid decreasing, stable, increasing
    -> diag group (N, EMCI, LMCI, AD)
    -> positivity threshold 
    '''
    print "THRESHOLD: %s" % threshold


    frontal_distr = {rid: (val-norm_fits['frontal'][0], norm_fits['frontal'][1]) for rid, val in uptake_diff['frontal'].iteritems()}
    parietal_distr = {rid: (val-norm_fits['parietal'][0], norm_fits['parietal'][1]) for rid, val in uptake_diff['parietal'].iteritems()}
    temporal_distr = {rid: (val-norm_fits['temporal'][0], norm_fits['temporal'][1]) for rid, val in uptake_diff['temporal'].iteritems()}
    cingulate_distr = {rid: (val-norm_fits['cingulate'][0], norm_fits['cingulate'][1]) for rid, val in uptake_diff['cingulate'].iteritems()}
    composite_distr = {rid: (val-norm_fits['composite'][0], norm_fits['composite'][1]) for rid, val in uptake_diff['composite'].iteritems()}
    frontal_cdf = {rid: 1.0 - norm.cdf(0.0, mean, std) for rid, (mean, std) in frontal_distr.iteritems()}
    parietal_cdf = {rid: 1.0 - norm.cdf(0.0, mean, std) for rid, (mean, std) in parietal_distr.iteritems()}
    temporal_cdf = {rid: 1.0 - norm.cdf(0.0, mean, std) for rid, (mean, std) in temporal_distr.iteritems()}
    cingulate_cdf = {rid: 1.0 - norm.cdf(0.0, mean, std) for rid, (mean, std) in cingulate_distr.iteritems()}
    composite_cdf = {rid: 1.0 - norm.cdf(0.0, mean, std) for rid, (mean, std) in composite_distr.iteritems()}

    #all_region_cdf_mean = {rid: np.mean([frontal_cdf[rid], parietal_cdf[rid], temporal_cdf[rid], cingulate_cdf[rid]]) for rid in frontal_cdf.keys()}
    all_region_cdf_mean = {rid: composite_cdf[rid] for rid in composite_cdf.keys()}

    # split by threshold
    high = [k for k,v in bl_uptake.iteritems() if v[0]['composite'] >= threshold]
    low = [k for k,v in bl_uptake.iteritems() if v[0]['composite'] < threshold]

    # split by longitudinal status
    increasing = [k for k,v in all_region_cdf_mean.iteritems() if v > 0.5]
    stable = [k for k,v in all_region_cdf_mean.iteritems() if v <= 0.5 and v >= 0.0]
    decreasing = [k for k,v in all_region_cdf_mean.iteritems() if v < 0.0]

    # graph histogram of slopes
    if graph:
        slope_values = uptake_diff['composite'].values()
        noise_distr_mean = norm_fits['composite'][0]
        noise_distr_std = norm_fits['composite'][1]
        pylab.figure()
        weights = np.ones_like(slope_values)/float(len(slope_values))
        n, bins, patches = pylab.hist(slope_values,40,weights=weights,histtype='stepfilled',label=['Composite Annualized Slope'])
        #pylab.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        # add a line showing the noise distribution
        y = np.array(pylab.normpdf(bins, noise_distr_mean, noise_distr_std))
        norm_weights = np.ones_like(y)/float(len(y))
        y *= norm_weights
        l = pylab.plot(bins, y, 'k--', linewidth=1.5)
        # first thirds
        first_third = norm.ppf(0.33, noise_distr_mean, noise_distr_std)
        second_third = norm.ppf(0.66, noise_distr_mean, noise_distr_std)
        half = norm.ppf(0.5, noise_distr_mean, noise_distr_std)
        y_max = max(y)
        pylab.plot([first_third, first_third],[0,y_max])
        pylab.plot([second_third, second_third],[0,y_max])
        pylab.plot([half, half],[0,y_max])
        pylab.legend()
        pylab.show()


    # mix longitudinal status and threshold
    increasing_high = list(set(increasing) & set(high))
    increasing_low = list(set(increasing) & set(low))
    stable_high = list(set(stable) & set(high))
    stable_low = list(set(stable) & set(low))
    decreasing_high = list(set(decreasing) & set(high))
    decreasing_low = list(set(decreasing) & set(low))

    # split by diagnosis
    increasing_high_n = [k for k in increasing_high if diags[k] in set(['N'])]
    stable_high_n = [k for k in stable_high if diags[k] in set(['N'])]
    decreasing_high_n = [k for k in decreasing_high if diags[k] in set(['N'])]
    increasing_low_n = [k for k in increasing_low if diags[k] in set(['N'])]
    stable_low_n = [k for k in stable_low if diags[k] in set(['N'])]
    decreasing_low_n = [k for k in decreasing_low if diags[k] in set(['N'])]

    increasing_high_emci = [k for k in increasing_high if diags[k] in set(['EMCI'])]
    stable_high_emci = [k for k in stable_high if diags[k] in set(['EMCI'])]
    decreasing_high_emci = [k for k in decreasing_high if diags[k] in set(['EMCI'])]
    increasing_low_emci = [k for k in increasing_low if diags[k] in set(['EMCI'])]
    stable_low_emci = [k for k in stable_low if diags[k] in set(['EMCI'])]
    decreasing_low_emci = [k for k in decreasing_low if diags[k] in set(['EMCI'])]

    increasing_high_lmci = [k for k in increasing_high if diags[k] in set(['LMCI'])]
    stable_high_lmci = [k for k in stable_high if diags[k] in set(['LMCI'])]
    decreasing_high_lmci = [k for k in decreasing_high if diags[k] in set(['LMCI'])]
    increasing_low_lmci = [k for k in increasing_low if diags[k] in set(['LMCI'])]
    stable_low_lmci = [k for k in stable_low if diags[k] in set(['LMCI'])]
    decreasing_low_lmci = [k for k in decreasing_low if diags[k] in set(['LMCI'])]

    increasing_high_ad = [k for k in increasing_high if diags[k] in set(['AD'])]
    stable_high_ad = [k for k in stable_high if diags[k] in set(['AD'])]
    decreasing_high_ad = [k for k in decreasing_high if diags[k] in set(['AD'])]
    increasing_low_ad = [k for k in increasing_low if diags[k] in set(['AD'])]
    stable_low_ad = [k for k in stable_low if diags[k] in set(['AD'])]
    decreasing_low_ad = [k for k in decreasing_low if diags[k] in set(['AD'])]

    print 'N'
    print "Decrease: %s/%s" % (len(decreasing_low_n),len(decreasing_high_n))
    print "Stable: %s/%s" % (len(stable_low_n),len(stable_high_n))
    print "Increase: %s/%s" % (len(increasing_low_n),len(increasing_high_n))
    print 'EMCI'
    print "Decrease: %s/%s" % (len(decreasing_low_emci),len(decreasing_high_emci))
    print "Stable: %s/%s" % (len(stable_low_emci),len(stable_high_emci))
    print "Increase: %s/%s" % (len(increasing_low_emci),len(increasing_high_emci))
    print 'LMCI'
    print "Decrease: %s/%s" % (len(decreasing_low_lmci),len(decreasing_high_lmci))
    print "Stable: %s/%s" % (len(stable_low_lmci),len(stable_high_lmci))
    print "Increase: %s/%s" % (len(increasing_low_lmci),len(increasing_high_lmci))
    print 'AD'
    print "Decrease: %s/%s" % (len(decreasing_low_ad),len(decreasing_high_ad))
    print "Stable: %s/%s" % (len(stable_low_ad),len(stable_high_ad))
    print "Increase: %s/%s" % (len(increasing_low_ad),len(increasing_high_ad))

    data = {'increasing_low' : {'N': increasing_low_n,
                            'EMCI': increasing_low_emci,
                            'LMCI': increasing_low_lmci,
                            'AD': increasing_low_ad},
            'stable_low' : {'N': stable_low_n,
                        'EMCI': stable_low_emci,
                        'LMCI': stable_low_lmci,
                        'AD': stable_low_ad},
            'decreasing_low': {'N': decreasing_low_n,
                           'EMCI': decreasing_low_emci,
                           'LMCI': decreasing_low_lmci,
                           'AD': decreasing_low_ad},
            'increasing_high' : {'N': increasing_high_n,
                            'EMCI': increasing_high_emci,
                            'LMCI': increasing_high_lmci,
                            'AD': increasing_high_ad},
            'stable_high' : {'N': stable_high_n,
                        'EMCI': stable_high_emci,
                        'LMCI': stable_high_lmci,
                        'AD': stable_high_ad},
            'decreasing_high': {'N': decreasing_high_n,
                           'EMCI': decreasing_high_emci,
                           'LMCI': decreasing_high_lmci,
                           'AD': decreasing_high_ad}}
    return data

def effectTests(pts1, pts2, name_one, name_two):
    W_1, pvalue_1 = shapiro(pts1)
    W_2, pvalue_2 = shapiro(pts2)
    print "(Shapiro) %s: %s\n%s: %s" % (name_one, pvalue_1, name_two, pvalue_2)
    t, pval = ttest_ind(pts1, pts2, equal_var=False)
    print "(TTEST) t: %s\npval: %s" % (t, pval)
    u, pvalue = mannwhitneyu(pts1, pts2, use_continuity=True)
    u_max = len(pts1) * len(pts2)
    rank_biserial = 1.0 - (2*u/u_max)
    print "(MANNWHITNEY): pval: %s\n rankbiserial: %s" % (pvalue, rank_biserial)

def calculatePVCvsNonPVCSlope(bl_file, scan2_file, scan3_file)
    # for PVC
    data_bl_pvc, data_scan2_pvc, data_scan3_pvc = parseRoussetOutputs(bl_file,scan2_file,scan3_file, pvcval=True)
    data_bl_nonpvc, data_scan2_nonpvc, data_scan3_nonpvc = parseRoussetOutputs(bl_file,scan2_file,scan3_file, pvcval=False)
    points = []
    for k,v in data_bl_pvc.iteritems():
        pvcval = v[0]['composite']
        nonpvcval = data_bl_nonpvc[k][0]['composite']
        res = res_bl[k]['group4']
        points.append((nonpvcval, pvcval, res))
    for k,v in data_scan2_pvc.iteritems():
        pvcval = v[0]['composite']
        nonpvcval = data_scan2_nonpvc[k][0]['composite']
        res = res_scan2[k]['group4']
        points.append((nonpvcval, pvcval, res))
    for k,v in data_scan3_pvc.iteritems():
        pvcval = v[0]['composite']
        nonpvcval = data_scan3_nonpvc[k][0]['composite']
        res = res_scan3[k]['group4']
        points.append((nonpvcval, pvcval, res))
    slope, intercept, r, p, stderr = linregress([_[0] for _ in points], [_[1] for _ in points])
    p = np.poly1d([slope, intercept])
    fit_threshold = p(fit_threshold)
    print "NEW THRESHOLD: %s" % fit_threshold



if __name__ == "__main__":
    # main files
    registry_file = "../docs/registry_clean.csv"
    master_file = '../FDG_AV45_COGdata_09_25_15.csv'

    # tp-specific file
    av45_file = "../output/UCBERKELEYAV45_09_25_15_extra.csv"

    # non tp-specific file
    av45_file_nontp = "../output/UCBERKELEYAV45_09_25_15_extra_nontp.csv"

    # Rousset output files
    rousset_matfile_bl_manual = '../output/Rousset_BL/rousset_outputs_manual.mat'
    rousset_matfile_scan2_manual = '../output/Rousset_Scan2/rousset_outputs_manual.mat'
    rousset_matfile_scan3_manual = '../output/Rousset_Scan3/rousset_outputs_manual.mat'
    rousset_matfile_bl_agg = '../output/Rousset_BL/rousset_outputs_agg.mat'
    rousset_matfile_scan2_agg = '../output/Rousset_Scan2/rousset_outputs_agg.mat'
    rousset_matfile_scan3_agg = '../output/Rousset_Scan3/rousset_outputs_agg.mat'
    rousset_agghigh_raw_bl = '../raw_agghigh_output_BL.mat'
    rousset_agghigh_raw_scan2 = '../raw_agghigh_output_Scan2.mat'
    rousset_agghigh_raw_scan3 = '../raw_agghigh_output_Scan3.mat'

    # get cluster values and rates of change
    '''
    raw_bl, raw_scan2, raw_scan3, index_lookup = parseRawRousset(rousset_agghigh_raw_bl, 
                                                                 rousset_agghigh_raw_scan2, 
                                                                 rousset_agghigh_raw_scan3)
    for k,v in raw_bl.iteritems():
        print k 
        print v
    sys.exit(1)
    '''

    # freesurfer region lookup
    lut_file = "../FreeSurferColorLUT.txt"
    master_data = importMaster(master_file)
    diags = {}
    for rid, row in master_data.iteritems():
        diag = row['Init_Diagnosis'].strip()
        diags[rid] = diag


    # grabbing residuals
    res_bl, res_scan2, res_scan3 = parseRoussetResiduals(rousset_matfile_bl_manual,rousset_matfile_scan2_manual,rousset_matfile_scan3_manual)
    res_bl_agg, res_scan2_agg, res_scan3_agg = parseRoussetResiduals(rousset_matfile_bl_agg,rousset_matfile_scan2_agg,rousset_matfile_scan3_agg)    
    group4_residuals = defaultdict(list)
    for k,v in res_bl.iteritems():
        group4_residuals[k].append(v['group4'])
    for k,v in res_scan2.iteritems():
        group4_residuals[k].append(v['group4'])
    for k,v in res_scan3.iteritems():
        group4_residuals[k].append(v['group4'])
    group4_residuals = dict(group4_residuals)
    for k,v in group4_residuals.iteritems():
        group4_residuals[k] = np.mean(v)
    res_values = sorted(group4_residuals.iteritems(), key=lambda x: x[1])
    group4_quantiles = {}
    for i,(k,v) in enumerate(res_values):
        qtile = float(i) / float(len(res_values))
        group4_quantiles[k] = qtile

    # split residuals into quartiles
    values_only = [_[1] for _ in res_values]
    keys_only = [_[0] for _ in res_values]
    quarter = len(res_values) / 4
    quartile_0 = [_ for _ in keys_only[:quarter]]
    quartile_25 = [_ for _ in keys_only[quarter:(quarter*2)]]
    quartile_50 = [_ for _ in keys_only[(quarter*2):(quarter*3)]]
    quartile_75 = [_ for _ in keys_only[(quarter*3):]]
    quartile_0_vals = [_ for _ in values_only[:quarter]]
    quartile_25_vals = [_ for _ in values_only[quarter:(quarter*2)]]
    quartile_50_vals = [_ for _ in values_only[(quarter*2):(quarter*3)]]
    quartile_75_vals = [_ for _ in values_only[(quarter*3):]]

    # graph residuals
    '''
    pylab.figure()
    n, bins, patches = pylab.hist([values_only],70,histtype='bar',label=['All Residuals'])
    pylab.legend()
    pylab.figure()
    n, bins, patches = pylab.hist([quartile_0_vals],70,histtype='step',label=['25 quartile'])
    n, bins, patches = pylab.hist([quartile_25_vals],70,histtype='step',label=['50 quartile'])
    n, bins, patches = pylab.hist([quartile_50_vals],70,histtype='step',label=['75 quartile'])
    n, bins, patches = pylab.hist([quartile_75_vals],70,histtype='step',label=['100 quartile'])
    pylab.legend()
    pylab.show()
    sys.exit(1)
    '''


    # composite threshold
    #fit_threshold = 1.11 # for whole cereb
    fit_threshold = 0.79 # for big ref (non PVC)

    # for PVC
    data_bl_pvc, data_scan2_pvc, data_scan3_pvc = parseRoussetOutputs(rousset_matfile_bl_manual,rousset_matfile_scan2_manual,rousset_matfile_scan3_manual, pvcval=True)
    data_bl_nonpvc, data_scan2_nonpvc, data_scan3_nonpvc = parseRoussetOutputs(rousset_matfile_bl_manual,rousset_matfile_scan2_manual,rousset_matfile_scan3_manual, pvcval=False)
    points = []
    for k,v in data_bl_pvc.iteritems():
        pvcval = v[0]['composite']
        nonpvcval = data_bl_nonpvc[k][0]['composite']
        res = res_bl[k]['group4']
        points.append((nonpvcval, pvcval, res))
    for k,v in data_scan2_pvc.iteritems():
        pvcval = v[0]['composite']
        nonpvcval = data_scan2_nonpvc[k][0]['composite']
        res = res_scan2[k]['group4']
        points.append((nonpvcval, pvcval, res))
    for k,v in data_scan3_pvc.iteritems():
        pvcval = v[0]['composite']
        nonpvcval = data_scan3_nonpvc[k][0]['composite']
        res = res_scan3[k]['group4']
        points.append((nonpvcval, pvcval, res))
    slope, intercept, r, p, stderr = linregress([_[0] for _ in points], [_[1] for _ in points])
    p = np.poly1d([slope, intercept])
    fit_threshold = p(fit_threshold)
    print "NEW THRESHOLD: %s" % fit_threshold
    
    # for TP
    data_bl_tp, data_scan2_tp, data_scan3_tp = parseAV45Output(av45_file, registry_file, lut_file)

    # for nonTP
    data_bl_nontp, data_scan2_nontp, data_scan3_nontp = parseAV45Output(av45_file_nontp, registry_file, lut_file)

    # find annualize change
    all_rids = list(set(data_scan2_nontp.keys() + data_scan3_nontp.keys()))
    
    print "Nontp vs TP, 2tp"
    vol_diff_1, uptake_diff_1, yrs_diff_1 = findRegionalAnnualizedChange(all_rids, data_bl_nontp, data_scan2_nontp, {}, master_data, annualize=True)
    vol_diff_2, uptake_diff_2, yrs_diff_2 = findRegionalAnnualizedChange(all_rids, data_bl_tp, data_scan2_tp, {}, master_data, annualize=True)
    effectTests(uptake_diff_1['composite'].values(), uptake_diff_2['composite'].values(), 'NONTP', 'TP')
    print "Nontp vs TP, 3tp"
    vol_diff_1, uptake_diff_1, yrs_diff_1 = findRegionalAnnualizedChange(all_rids, data_bl_nontp, {}, data_scan3_nontp, master_data, annualize=True)
    vol_diff_2, uptake_diff_2, yrs_diff_2 = findRegionalAnnualizedChange(all_rids, data_bl_tp, {}, data_scan3_tp, master_data, annualize=True)
    effectTests(uptake_diff_1['composite'].values(), uptake_diff_2['composite'].values(), 'NONTP', 'TP')
    print "Nontp vs TP, alltp"
    vol_diff_1, uptake_diff_1, yrs_diff_1 = findRegionalAnnualizedChange(all_rids, data_bl_nontp, data_scan2_nontp, data_scan3_nontp, master_data, annualize=True)
    vol_diff_2, uptake_diff_2, yrs_diff_2 = findRegionalAnnualizedChange(all_rids, data_bl_tp, data_scan2_tp, data_scan3_tp, master_data, annualize=True)
    effectTests(uptake_diff_1['composite'].values(), uptake_diff_2['composite'].values(), 'NONTP', 'TP')
    print "PVC vs Nontp, 2tp"
    vol_diff_1, uptake_diff_1, yrs_diff_1 = findRegionalAnnualizedChange(all_rids, data_bl_pvc, data_scan2_pvc, {}, master_data, annualize=True)
    vol_diff_2, uptake_diff_2, yrs_diff_2 = findRegionalAnnualizedChange(all_rids, data_bl_nontp, data_scan2_nontp, {}, master_data, annualize=True)
    effectTests(uptake_diff_1['composite'].values(), uptake_diff_2['composite'].values(), 'PVC', 'NONTP')
    print "PVC vs Nontp, 3tp"
    vol_diff_1, uptake_diff_1, yrs_diff_1 = findRegionalAnnualizedChange(all_rids, data_bl_pvc, {}, data_scan3_pvc, master_data, annualize=True)
    vol_diff_2, uptake_diff_2, yrs_diff_2 = findRegionalAnnualizedChange(all_rids, data_bl_nontp, {}, data_scan3_nontp, master_data, annualize=True)
    effectTests(uptake_diff_1['composite'].values(), uptake_diff_2['composite'].values(), 'PVC', 'NONTP')
    print "PVC vs Nontp, alltp"
    vol_diff_1, uptake_diff_1, yrs_diff_1 = findRegionalAnnualizedChange(all_rids, data_bl_pvc, data_scan2_pvc, data_scan3_pvc, master_data, annualize=True)
    vol_diff_2, uptake_diff_2, yrs_diff_2 = findRegionalAnnualizedChange(all_rids, data_bl_nontp, data_scan2_nontp, data_scan3_nontp, master_data, annualize=True)
    effectTests(uptake_diff_1['composite'].values(), uptake_diff_2['composite'].values(), 'PVC', 'NONTP')
    print "PVC vs TP, 2tp"
    vol_diff_1, uptake_diff_1, yrs_diff_1 = findRegionalAnnualizedChange(all_rids, data_bl_pvc, data_scan2_pvc, {}, master_data, annualize=True)
    vol_diff_2, uptake_diff_2, yrs_diff_2 = findRegionalAnnualizedChange(all_rids, data_bl_tp, data_scan2_tp, {}, master_data, annualize=True)
    effectTests(uptake_diff_1['composite'].values(), uptake_diff_2['composite'].values(), 'PVC', 'TP')
    print "PVC vs TP, 3tp"
    vol_diff_1, uptake_diff_1, yrs_diff_1 = findRegionalAnnualizedChange(all_rids, data_bl_pvc, {}, data_scan3_pvc, master_data, annualize=True)
    vol_diff_2, uptake_diff_2, yrs_diff_2 = findRegionalAnnualizedChange(all_rids, data_bl_tp, {}, data_scan3_tp, master_data, annualize=True)
    effectTests(uptake_diff_1['composite'].values(), uptake_diff_2['composite'].values(), 'PVC', 'TP')
    print "PVC vs TP, alltp"
    vol_diff_1, uptake_diff_1, yrs_diff_1 = findRegionalAnnualizedChange(all_rids, data_bl_pvc, data_scan2_pvc, data_scan3_pvc, master_data, annualize=True)
    vol_diff_2, uptake_diff_2, yrs_diff_2 = findRegionalAnnualizedChange(all_rids, data_bl_tp, data_scan2_tp, data_scan3_tp, master_data, annualize=True)
    effectTests(uptake_diff_1['composite'].values(), uptake_diff_2['composite'].values(), 'PVC', 'TP')


    sys.exit(1)

    # graph pvc vs nonpvcval, colored by residual quartile
    '''
    quarter = len(points)/4 
    x = []
    y = []
    c = []
    for i, (xs,ys,cs) in enumerate(points):
        if i < quarter:
            c.append('y')
        elif i < quarter*2:
            c.append('g')
        elif i < quarter*3:
            c.append('r')
        else:
            c.append('b')
        x.append(xs)
        y.append(ys)
    yhat = p(x)
    ydiff = (y-yhat)**2
    pylab.figure()
    n, bins, patches = pylab.hist([_ for i,_ in enumerate(ydiff) if c[i] == 'y'],30,normed=1,histtype='step',label=['25 quartile'])
    n, bins, patches = pylab.hist([_ for i,_ in enumerate(ydiff) if c[i] == 'g'],30,normed=1,histtype='step',label=['50 quartile'])
    n, bins, patches = pylab.hist([_ for i,_ in enumerate(ydiff) if c[i] == 'r'],30,normed=1,histtype='step',label=['75 quartile'])
    n, bins, patches = pylab.hist([_ for i,_ in enumerate(ydiff) if c[i] == 'b'],30,normed=1,histtype='step',label=['100 quartile'])
    pylab.legend()
    pylab.show()
    sys.exit(1)
    '''

    # dump the dataset
    '''
    output_file = "../nontp_regional_change.csv"
    dumpRegionalDataset(output_file, data_bl, data_scan2, data_scan3, master_data)
    '''

    # # find scan-rescan annualized noise
    # vol_diff_stable, uptake_diff_stable, yrs_diff_stable = findRegionalAnnualizedChange(STABLE, data_bl, data_scan2, data_scan3, master_data, annualize=True)
    # norm_fits = fitNormalToUptakeChange(uptake_diff_stable)
    # for k,v in norm_fits.iteritems():
    #     print k
    #     print "(%.5f, %.5f)" % v

    # graph residual quantiles versus uptake change
    '''
    points = []
    uptake_diff_composite = uptake_diff['composite']
    for k,v in uptake_diff_composite.iteritems():
        res = group4_quantiles[k]
        points.append((res, v))
    x = [_[0] for _ in points]
    y = [_[1] for _ in points]
    plt.figure()
    plt.scatter(x,y)
    plt.show()
    sys.exit(1)
    '''

    # stratify into subject groups
    subj_groups = stratifySubjects(data_bl, uptake_diff, norm_fits, yrs_diff, diags, threshold=fit_threshold, graph=False)
    print subj_groups

    # find percent plausible
    #percentPlausible(vol_diff, uptake_diff, norm_fits, yrs_diff)
    percentPlausibleNoNoise(vol_diff, uptake_diff, yrs_diff)

    # calculate diagnosis group effects
    diagGroupEffects(data_bl, diags)

    # plot volume change versus uptake change (marking high amyloid at baseline)
    '''
    colors = {}
    for rid in all_rids:
        if rid in subj_groups['increasing']['N']:
            c = 'k'
        elif rid in subj_groups['increasing']['EMCI']:
            c = 'y'
        elif rid in subj_groups['increasing']['LMCI']:
            c = 'g'
        elif rid in subj_groups['increasing']['AD']:
            c= 'r'
        else:
            c= 'w'
        colors[rid] = c
    plot_regional_volumechange_vs_uptakechange(colors, vol_diff, uptake_diff)
    '''

