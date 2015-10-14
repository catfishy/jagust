import sys
import csv
import os
import scipy.io as sio
from scipy.spatial.distance import euclidean
from sklearn.cluster import KMeans, AgglomerativeClustering
from collections import defaultdict
from datetime import datetime, timedelta
import itertools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA

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

def parseRawRousset(bl_file, scan2_file, scan3_file):
    wcereb_regions = set([7,8,46,47])
    data_bl_raw = importRawRoussetResults(bl_file)
    data_scan2_raw = importRawRoussetResults(scan2_file)
    data_scan3_raw = importRawRoussetResults(scan3_file)
    data_bl = {}
    data_scan2 = {}
    data_scan3 = {}
    index_lookup = None
    for rid, data in data_bl_raw.iteritems():
        names = [unwrap(_) for _ in data['names']]
        sizes = [unwrap(_) for _ in data['sizes']]
        pvcvals = [unwrap(_) for _ in data['pvcvals']]
        indices = [unwrap(_) for _ in data['indices']]
        indices = [[_] if (not (isinstance(_, list) or isinstance(_, np.ndarray))) else _ for _ in indices]
        vals = dict(zip(names, pvcvals))
        # add whole cereb
        regionsizes = [(sz,pv) for i,sz,pv in zip(indices, sizes, pvcvals) if len(set(i) & wcereb_regions) > 0]
        regionvalues = [_[1] for _ in regionsizes]
        regionweights = np.array([_[0] for _ in regionsizes])
        regionweights = regionweights / float(sum(regionweights))
        vals['whole_cerebellum'] = sum(regionweights*regionvalues)
        if index_lookup is None:
            index_lookup = dict(zip(names, indices))
        data_bl[rid] = vals
    for rid, data in data_scan2_raw.iteritems():
        names = [unwrap(_) for _ in data['names']]
        sizes = [unwrap(_) for _ in data['sizes']]
        pvcvals = [unwrap(_) for _ in data['pvcvals']]
        indices = [unwrap(_) for _ in data['indices']]
        indices = [[_] if (not (isinstance(_, list) or isinstance(_, np.ndarray))) else _ for _ in indices]
        vals = dict(zip(names, pvcvals))
        # add whole cereb
        regionsizes = [(sz,pv) for i,sz,pv in zip(indices, sizes, pvcvals) if len(set(i) & wcereb_regions) > 0]
        regionvalues = [_[1] for _ in regionsizes]
        regionweights = np.array([_[0] for _ in regionsizes])
        regionweights = regionweights / float(sum(regionweights))
        vals['whole_cerebellum'] = sum(regionweights*regionvalues)
        data_scan2[rid] = vals
    for rid, data in data_scan2_raw.iteritems():
        names = [unwrap(_) for _ in data['names']]
        sizes = [unwrap(_) for _ in data['sizes']]
        pvcvals = [unwrap(_) for _ in data['pvcvals']]
        indices = [unwrap(_) for _ in data['indices']]
        indices = [[_] if (not (isinstance(_, list) or isinstance(_, np.ndarray))) else _ for _ in indices]
        vals = dict(zip(names, pvcvals))
        # add whole cereb
        regionsizes = [(sz,pv) for i,sz,pv in zip(indices, sizes, pvcvals) if len(set(i) & wcereb_regions) > 0]
        regionvalues = [_[1] for _ in regionsizes]
        regionweights = np.array([_[0] for _ in regionsizes])
        regionweights = regionweights / float(sum(regionweights))
        vals['whole_cerebellum'] = sum(regionweights*regionvalues)
        data_scan3[rid] = vals
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


def parseCSV(file_path, delimiter=','):
    data = pd.read_csv(file_path, sep=delimiter, low_memory=False)
    data.rename(columns=lambda x: x.strip(), inplace=True)
    lines = [dict(data.iloc[i].replace(np.nan, '')) for i in range(len(data))]
    headers = list(data.columns.values)
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


def updateLine(old_line, new_data, extraction_fn, 
               pid_key='RID', pet_meta=None, decimal_places=4):
    try:
        subj = int(old_line[pid_key])
    except Exception as e:
        print "No subject column %s found" % pid_key
        return {}

    subj_row = new_data.get(subj,None)
    if subj_row is None:
        print "No subj row found for %s" % (subj)
        return {}

    patient_pets = None
    if pet_meta is not None:
        patient_pets = sorted(pet_meta.get(subj,[]))

    new_data = extraction_fn(subj, subj_row, old_line, patient_pets) # patient_pets is passed in as context
    new_data = convertToCSVDataType(new_data, decimal_places=decimal_places)
    return new_data


def importRegistry(registry_file, include_all=False):
    headers, lines = parseCSV(registry_file)
    registry = defaultdict(list)
    for data in lines:
        if not include_all and data['EXAMDATE'] == '':
            continue
        subj = int(data['RID'])
        date_str = data['EXAMDATE'].strip().lower()
        date = None
        try:
            date = datetime.strptime(date_str,'%Y-%m-%d')
        except Exception as e:
            pass
        try:
            date = datetime.strptime(date_str,'%m/%d/%y')
        except Exception as e:
            pass
        if not include_all and date is None:
            continue
        registry[subj].append({'VISCODE': data['VISCODE'].strip().lower(),
                               'VISCODE2': data['VISCODE2'].strip().lower(),
                               'EXAMDATE': date,
                               'update_stamp': data['update_stamp']})
    registry = dict(registry)
    for k in registry.keys():
        new_val = sorted(registry[k], key=lambda x: x['EXAMDATE'])
        registry[k] = new_val
    return registry

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
        date = datetime.strptime(date_string,'%Y-%m-%d')
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
        subj = int(line['RID'])
        viscode = line['VISCODE'].strip().lower()
        viscode2 = line['VISCODE2'].strip().lower()
        #if viscode == 'sc' or viscode2 == 'sc':
        #    continue
        examdate = line['EXAMDATE']
        change = line['DXCHANGE']
        current = line['DXCURREN']
        conv = line['DXCONV']
        conv_type = line['DXCONTYP']
        rev_type = line['DXREV']

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
            examdate = datetime.strptime(examdate,'%Y-%m-%d')
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

def importFreesurferLookup(lut_file):
    infile = open(lut_file, 'r')
    data = {}
    for line in infile.readlines():
        parts = [_.strip() for _ in line.split(' ') if _.strip() != '']
        if len(parts) < 2:
            continue
        try:
            idx = int(parts[0])
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
        date = datetime.strptime(line['EXAMDATE'],'%m/%d/%y')
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
        examdate = datetime.strptime(date_str,'%Y-%m-%d')
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
            try:
                dob = datetime.strptime(line['PTDOB'].replace('--','01'),'%m/%d/%Y')
            except:
                dob = datetime.strptime(line['PTDOB'].replace('--','01'),'%m/%d/%y')
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
            vs = line['VISCODE'].strip().lower()
            vs2 = line['VISCODE2'].strip().lower()
            subject_registry = registry.get(subj,[])
            date = None
            for v in subject_registry:
                item_date = v['EXAMDATE']
                if vs2 == v['VISCODE2']:
                    date = item_date
                    break
            if date is None:
                print "Could not find visit in registry: %s, %s, %s" % (subj, vs, vs2)
        else:
            date = datetime.strptime(date_string,'%m/%d/%y')
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
        viscode = line['VISCODE'].strip().lower()
        if 'VISCODE2' in line:
            viscode2 = line['VISCODE2'].strip().lower()
        else:
            viscode2 = ''
        examdate = line.get('EXAMDATE','')
        if examdate != '':
            examdate = datetime.strptime(examdate,'%Y-%m-%d')
        elif registry is not None:
            examdate = findVisitDate(registry, subj, viscode, viscode2)
            if not examdate:
                print "Could not find exam date for %s (%s, %s)" % (subj, viscode, viscode2)
                continue
        tots = [line['AVTOT%s' % _ ]for _ in range(1,6)]

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
                                   'EXAMDATE': examdate,
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
        examdate = datetime.strptime(line['EXAMDATE'],'%Y-%m-%d')
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
            examdate = findVisitDate(registry, subj, viscode, viscode2)
        if examdate is None:
            print "No exam date for %s (%s, %s)" % (subj, viscode, viscode2)
            continue
        adas_by_subj[subj].append({'VISCODE': viscode,
                                   'VISCODE2': viscode2,
                                   'EXAMDATE': examdate,
                                   'TOTSCORE': totscore})
    return dict(adas_by_subj)

def importTBMSyn(tbm_file):
    headers, lines = parseCSV(tbm_file)
    data = defaultdict(list)
    for i, line in enumerate(lines):  
        if line['ACCELERATED'] == 'Yes':
            continue      
        subj = int(line['RID'])
        vc = line['VISCODE'].strip().lower()
        vc2 = line['VISCODE2'].strip().lower()
        vcbl = line['VISCODE2BL'].strip().lower()
        bl_examdate = datetime.strptime(line['EXAMDATEBL'],'%Y-%m-%d')
        score = float(line['TBMSYNSCOR'])
        examdate = datetime.strptime(line['EXAMDATE'],'%Y-%m-%d')
        data[subj].append({'VISCODE': vc,
                           'VISCODE2': vc2,
                           'VISCODEBL': vcbl,
                           'EXAMDATE': examdate,
                           'BL_EXAMDATE': bl_examdate,
                           'SCORE': score})
    return data

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
        try:
            examdate = datetime.strptime(line['EXAMDATE'],'%m/%d/%y')
        except:
            examdate = datetime.strptime(line['EXAMDATE'],'%Y-%m-%d')
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


def importGD_DOD(gd_file):
    headers, lines = parseCSV(gd_file)
    data = defaultdict(dict)
    for line in lines:
        subj = int(line['SCRNO'])
        gdtotal = int(line['GDTOTAL'])
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
            examdate = datetime.strptime(line['EXAMDATE'],'%Y-%m-%d')
        except:
            examdate = findVisitDate(registry, subj, vc, vc2) if registry is not None else None
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
            try:
                rundate = datetime.strptime(line['RUNDATE'],'%m/%d/%y')
            except:
                rundate = datetime.strptime(line['RUNDATE'],'%Y-%m-%d')
            try:
                vc2 = line['VISCODE2'].strip().lower()
            except:
                vc2 = ''
            try:
                examdate = datetime.strptime(line['EXAMDATE'],'%m/%d/%y')
            except:
                examdate = findVisitDate(registry, subj, vc, vc2) if registry is not None else None
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
        print flattened
        data[k] = sorted(flattened, key=lambda x: x['EXAMDATE'])
    return data

def importPetMETA(pet_meta_file):
    headers, lines = parseCSV(pet_meta_file)
    pets = defaultdict(list)
    for row in lines:
        subj = int(row['Subject'].split('_')[-1].strip())
        new_date = datetime.strptime(row['Scan Date'], '%m/%d/%y')
        pets[subj].append(new_date)
    pets = dict(pets)
    for k in pets.keys():
        pets[k] = sorted(pets[k])
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
        userdate = datetime.strptime(data['USERDATE'],'%Y-%m-%d')
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
        try:
            date = datetime.strptime(line['EXAMDATE'], '%Y-%m-%d')
        except:
            date = datetime.strptime(line['EXAMDATE'], '%m/%d/%y')
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
            date = datetime.strptime(new_date,'%Y-%m-%d')
        except Exception as e:
            pass
        try:
            date = datetime.strptime(new_date,'%m/%d/%y')
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
    headers, lines = parseCSV(master_file)
    data = {}
    for i, line in enumerate(lines):
        # get subject ID
        try:
            subj = int(line['RID'])
        except Exception as e:
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
            bl_examdate = datetime.strptime(line['EXAMDATE'],'%Y-%m-%d')
            bl_examdates[subj] = bl_examdate
            continue
        vc = line['VISCODE'].strip().lower()
        vc2 = line['VISCODE2'].strip().lower()
        examdate = datetime.strptime(line['EXAMDATE'],'%Y-%m-%d')
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


def importAV45(av45_file, av45_nontp_file=None ,registry=None):
    av45_headers, av45_lines = parseCSV(av45_file)

    if av45_nontp_file:
        av45_nontp_headers, av45_nontp_lines = parseCSV(av45_nontp_file)
    else:
        av45_nontp_headers = []
        av45_nontp_lines = []

    av45_by_subj = defaultdict(list)
    for line in av45_lines:
        if 'SCRNO' in line:
            subj = int(line.pop('SCRNO',None))
        elif 'RID' in line:
            subj = int(line.pop('RID',None))
        elif 'PID' in line:
            subj = int(line.pop('PID',None))
        else:
            raise Exception("Can't find subject column in AV45 file")

        viscode = line['VISCODE'].strip().lower()
        if 'VISCODE2' in line:
            viscode2 = line['VISCODE2'].strip().lower()
        else:
            viscode2 = ''
        examdate = line.get('EXAMDATE',None)
        if examdate:
            try:
                examdate = datetime.strptime(examdate,'%Y-%m-%d')
            except:
                examdate = datetime.strptime(examdate,'%m/%d/%y')
        elif registry is not None:
            examdate = findVisitDate(registry, subj, viscode, viscode2)
        if not examdate:
            print "Could not find exam date for %s (%s, %s)" % (subj, viscode, viscode2)
            #continue
        line['EXAMDATE'] = examdate
        line['TP_SPECIFIC'] = True
        av45_by_subj[subj].append(line)

    # add in nontp data if given
    for line in av45_nontp_lines:
        if 'SCRNO' in line:
            subj = int(line.pop('SCRNO',None))
        elif 'RID' in line:
            subj = int(line.pop('RID',None))
        elif 'PID' in line:
            subj = int(line.pop('PID',None))
        else:
            raise Exception("Can't find subject column in AV45 file")

        if subj in av45_by_subj:
            continue

        viscode = line['VISCODE'].strip().lower()
        if 'VISCODE2' in line:
            viscode2 = line['VISCODE2'].strip().lower()
        else:
            viscode2 = ''
        examdate = line.get('EXAMDATE',None)
        if examdate:
            try:
                examdate = datetime.strptime(examdate,'%Y-%m-%d')
            except:
                examdate = datetime.strptime(examdate,'%m/%d/%y')
        elif registry is not None:
            examdate = findVisitDate(registry, subj, viscode, viscode2)
        if not examdate:
            print "Could not find exam date for %s (%s, %s)" % (subj, viscode, viscode2)
            #continue
        line['EXAMDATE'] = examdate
        line['TP_SPECIFIC'] = False
        av45_by_subj[subj].append(line)

    return dict(av45_by_subj)

def importCrossSectionFreesurfer(crossfree_file, include_failed=False):
    headers, lines = parseCSV(crossfree_file)
    data = defaultdict(list)
    failed = 0
    for i, line in enumerate(lines):
        if not include_failed and line['OVERALLQC'] == 'Fail':
            failed += 1
            continue
        elif not include_failed and line['OVERALLQC'] == 'Partial' and line['VENTQC'] == 'Fail':
            failed += 1
            continue
        subj = int(line['RID'])
        vc = line['VISCODE'].strip().lower()
        vc2 = line['VISCODE2'].strip().lower()
        examdate = datetime.strptime(line['EXAMDATE'],'%Y-%m-%d')
        inner_data = {k: v for k,v in line.iteritems() if k.startswith('ST') and k not in set(['STATUS'])}
        data[subj].append({'VISCODE': vc,
                           'VISCODE2': vc2,
                           'EXAMDATE': examdate,
                           'inner_data': inner_data})
    print "CROSS FREESURFER failed: %s" % failed
    return dict(data)


def importLongitudinalFreesurfer(longfree_file, include_failed = False):
    headers, lines = parseCSV(longfree_file)
    data = defaultdict(list)
    failed = 0
    for i, line in enumerate(lines):
        if not include_failed and line['OVERALLQC'] == 'Fail':
            failed += 1
            continue
        elif not include_failed and line['OVERALLQC'] == 'Partial' and line['VENTQC'] == 'Fail':
            failed += 1
            continue
        subj = int(line['RID'])
        vc = line['VISCODE'].strip().lower()
        vc2 = line['VISCODE2'].strip().lower()
        examdate = datetime.strptime(line['EXAMDATE'],'%Y-%m-%d')
        inner_data = {k: v for k,v in line.iteritems() if k.startswith('ST') and k not in set(['STATUS'])}
        data[subj].append({'VISCODE': vc,
                           'VISCODE2': vc2,
                           'EXAMDATE': examdate,
                           'inner_data': inner_data})
    print "LONG FREESURFER failed: %s" % failed
    return dict(data)


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

def findVisitDate(registry, subj, vc, vc2):
    subj_listings = registry.get(subj,[])
    examdate = None
    for listing in subj_listings:
        if vc2 != '' and listing['VISCODE'] == vc and listing['VISCODE2'] == vc2:
            examdate = listing['EXAMDATE']
            break
        elif vc2 == '' and listing['VISCODE'] == vc:
            examdate = listing['EXAMDATE']
            break
    return examdate

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
    # Calculate DOD raw output differences
    '''
    lut_file = "../FreeSurferColorLUT.txt"
    data = importFreesurferLookup(lut_file)
    file1 = '../docs/AV45_DOD_preprocess_output_12_18_14/AV45_BL_means_extraregions_18-Dec-2014_99.csv'
    file2 = '../docs/AV45_DOD_preprocess_output_08_05_15/AV45_BL_means_extraregions_05-Aug-2015_121.csv'
    diff_distr = calculateCSVDifference(file1, file2, index='0')
    plt.figure(1)
    keys = diff_distr.keys()
    x = range(len(keys))
    avgs = [diff_distr[k][0] for k in keys]
    stds = [diff_distr[k][1] for k in keys]
    maxs = [diff_distr[k][2] for k in keys]
    mins = [diff_distr[k][3] for k in keys]
    plt.errorbar(x, avgs, stds)
    plt.plot(x,maxs)
    plt.plot(x,mins)
    x1,x2,y1,y2 = plt.axis()
    plt.axis((0,len(x),y1,y2))
    locs, labels = plt.xticks()
    keys = [data.get(int(k),k) for k in keys]
    plt.xticks(x, keys, rotation='vertical')
    plt.show()
    sys.exit(1)
    '''

    # Scatter plot of composite SUVR values
    file1 = '../docs/DOD/AV45_BL_means_Dec2014_99.csv'
    file2 = '../DOD_DATA_08_05_15.csv'
    headers1, rows1 = parseCSV(file1)
    headers2, rows2 = parseCSV(file2)
    data1 = {r['PID']:r['comp/wcerb'] for r in rows1}
    data2 = {r['PID']:r['AV45_1_comp/wcerb'] for r in rows2}
    points = []
    for pid in list(set(data1.keys()) & set(data2.keys())):
        if pid < 6000:
            continue
        xval = float(data1[pid])
        yval = float(data2[pid])
        points.append((xval,yval))
    x = [_[0] for _ in points]
    y = [_[1] for _ in points]
    plt.figure(1)
    plt.scatter(x,y)
    plt.plot([0.8,1.7],[0.8,1.7])
    
    plt.xlabel('Old Value')
    plt.ylabel('New Value')
    plt.show()
    # Import FDG raw data
    '''
    fdg_extract_file = "../docs/FDG_preprocess_output_01_20_15/Tue-Jan-20-15-38-15-2015_FDGmetaroi.csv"
    registry_file = "../docs/registry_clean.csv"
    registry = importRegistry(registry_file)
    fdg_rows = importExtractedFDG(fdg_extract_file, registry=registry)
    fdg_rows = [convertToCSVDataType(_, decimal_places=5) for _ in fdg_rows]
    fdg_output = '../UCBERKELEY_FDG_07_29_15.csv'
    headers = ['RID','VISCODE','VISCODE2','UID','EXAMDATE','ROINAME','ROILAT','MEAN','MEDIAN','MODE','MIN','MAX','STDEV','NANVOX','TOTVOX','update_stamp']
    dumpCSV(fdg_output, headers, fdg_rows)
    '''

    # Get ventricle sections
    '''
    for k,v in data.iteritems():
        if 'ventricle' in v.lower():
            print "%s: %s" % (k,v)
    sys.exit(1)
    '''

    # Suzanne Other ROI
    '''
    to_lookup = [5,14,15,18,24,26,28,30,31,44,54,58,60,62,63,72,77,80,85,
                 251,252,253,254,255,1000,1001,1004,1005,1007,1009,1010,
                 1013,1016,1017,1021,1022,1024,1026,1030,1033,1034,1035,
                 2000,2001,2004,2005,2007,2009,2010,2013,2016,2017,2021,
                 2022,2024,2026,2030,2033,2034,2035]
    '''

    # All in aparc+aseg
    to_lookup = [0,2,4,5,7,8,10,11,12,13,14,15,16,17,18,24,26,28,30,31,41,43,
                 44,46,47,49,50,51,52,53,54,58,60,62,63,72,77,80,85,251,252,
                 253,254,255,1000,1001,1002,1003,1004,1005,1006,1007,1008,1009,
                 1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,
                 1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,
                 1034,1035,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,
                 2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,
                 2022,2023,2024,2025,2026,2027,2028,2029,2030,2031,2032,2033,
                 2034,2035]


    # DOD Output
    '''
    to_lookup = [3000,3001,3002,3003,3004,4000,8,47,5003,5001,5002,5000,16,12,
                 51,11,50,13,52,1003,1012,1014,1018,1019,1020,1027,1028,1032,
                 2003,2012,2014,2018,2019,2020,2027,2028,2032,1002,1010,1023,
                 1026,2002,2010,2023,2026,1008,1025,1029,1031,2008,2025,2029,
                 2031,1015,1030,2015,2030]
    '''

    '''
    lut_file = "../FreeSurferColorLUT.txt"
    data = importFreesurferLookup(lut_file)
    for idx in to_lookup:
        print "%s: %s" % (idx, data.get(idx, 'Unknown'))
    '''
