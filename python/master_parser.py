'''

Pull in the master FDG_AV45_COG DATA

remove null char:
tr -d '\000' < file1 > file2

'''
import re
from collections import defaultdict
from datetime import datetime, timedelta
from scipy import optimize
from scipy.stats import linregress
import numpy as np
import itertools

from utils import *

'''
def syncTemplate(old_headers, old_lines, input_file, dump_to=None):
    data = importTemplate(input_file)

    to_add_headers = []
    after = None
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=after)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        pass

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, data, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)
'''



def syncRoussetResults(old_headers, old_lines, rousset_matfile, timepoint, dump_to=None):
    assert timepoint in set(['BL', 'Scan2', 'Scan3'])

    group_names, sorted_data = importRoussetResults(rousset_matfile)

    valid_timepoints = ['BL', 'Scan2', 'Scan3']
    valid_groupings = ['summary']
    to_add_headers = []
    for vg in valid_groupings:
        to_add_headers += ['AV45_PVC_%s_CorticalSummary/Wholecereb_%s' % (vg, tp) for tp in valid_timepoints]
        to_add_headers += ['AV45_PVC_%s_CorticalSummary_slope_2points' % (vg), 'AV45_PVC_%s_CorticalSummary_slope_3points' % (vg)]
        to_add_headers += ['AV45_PVC_%s_CorticalSummary/Bigref_%s' % (vg, tp) for tp in valid_timepoints]
        to_add_headers += ['AV45_PVC_%s_CorticalSummary/Bigref_slope_2points' % (vg), 'AV45_PVC_%s_CorticalSummary/Bigref_slope_3points' % (vg)]
        to_add_headers += ['AV45_PVC_%s_HemiWM/Wholecereb_%s' % (vg, tp) for tp in valid_timepoints]
        to_add_headers += ['AV45_PVC_%s_HemiWM_slope_2points' % (vg), 'AV45_PVC_%s_HemiWM_slope_3points' % (vg)]
        to_add_headers += ['AV45_PVC_%s_Wholecereb_%s' % (vg, tp) for tp in valid_timepoints]
        to_add_headers += ['AV45_PVC_%s_CerebGM_%s' % (vg, tp) for tp in valid_timepoints]
        to_add_headers += ['AV45_PVC_%s_CerebWM_%s' % (vg, tp) for tp in valid_timepoints]
        if vg == 'summary':
            # SOME HARDCODED THRESHOLDING
            to_add_headers += ['AV45_PVC_summary_accumulator.0081', 'AV45_PVC_summary_wcereb1.27_BL']
    after = old_headers[max(i for i,_ in enumerate(old_headers) if _.startswith('AV45_') and 'PVC' not in _)] # last element that contains 'AV45'
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=after)
    timepoint_headers = [t for t in to_add_headers if timepoint in t]
    wipeKeys(old_lines, timepoint_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        #valid_groupings = ['summary']
        new_subj_data = {}
        for k in subj_row.keys():
            if k in valid_groupings:
                v = subj_row[k]
                if v is None:
                    continue
                wholecereb = float(v['wholecereb']['pvcval'])
                bigref = float(v['bigref']['pvcval'])
                for region, values in v.iteritems():
                    pvcval = float(values.get('pvcval',0.0))
                    if region == 'composite':
                        new_subj_data['AV45_PVC_%s_CorticalSummary/Wholecereb' % k] = pvcval/wholecereb
                        new_subj_data['AV45_PVC_%s_CorticalSummary/Bigref' % k] = pvcval/bigref
                    elif region == 'cerebGM':
                        new_subj_data['AV45_PVC_%s_CerebGM' % k] = pvcval
                    elif region == 'hemiWM':
                        new_subj_data['AV45_PVC_%s_HemiWM/Wholecereb' % k] = pvcval/wholecereb
                    elif region == 'wholecereb':
                        new_subj_data['AV45_PVC_%s_Wholecereb' % k] = pvcval
                    elif region == 'cerebWM':
                        new_subj_data['AV45_PVC_%s_CerebWM' % k] = pvcval
        return new_subj_data

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, sorted_data, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        new_data = {'%s_%s' % (k, timepoint): v for k,v in new_data.iteritems()}
        old_l.update(new_data)
        if timepoint == 'BL':
            bl_val = old_l.get('AV45_PVC_summary_CorticalSummary/Wholecereb_BL',0.0)
            try:
                old_l['AV45_PVC_summary_wcereb1.27_BL'] = 1 if float(bl_val) >= 1.27 else 0
            except ValueError as e:
                old_l['AV45_PVC_summary_wcereb1.27_BL'] = ''
        elif timepoint == 'Scan2':
            # try to calculate slope with 2 points
            for vg in valid_groupings:
                try:
                    diff1 = float(old_l.get('AV45_1_2_Diff'))
                    x = [0.0, diff1]

                    val1 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Wholecereb_BL' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Wholecereb_Scan2' % vg))
                    slope, intercept, r, p, stderr = linregress(x, [val1, val2])
                    old_l['AV45_PVC_%s_CorticalSummary_slope_2points' % (vg)] = slope
                    if vg == 'summary':
                        old_l['AV45_PVC_summary_accumulator.0081'] = 1 if slope >= 0.0081 else 0

                    val1 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Bigref_BL' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Bigref_Scan2' % vg))
                    slope, intercept, r, p, stderr = linregress(x, [val1, val2])
                    old_l['AV45_PVC_%s_CorticalSummary/Bigref_slope_2points' % (vg)] = slope

                    val1 = float(old_l.get('AV45_PVC_%s_HemiWM/Wholecereb_BL' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_HemiWM/Wholecereb_Scan2' % vg))
                    slope, intercept, r, p, stderr = linregress(x, [val1, val2])
                    old_l['AV45_PVC_%s_HemiWM_slope_2points' % (vg)] = slope
                except Exception as e:
                    old_l['AV45_PVC_%s_CorticalSummary_slope_2points' % (vg)] = ''
                    old_l['AV45_PVC_%s_CorticalSummary/Bigref_slope_2points' % (vg)] = ''
                    old_l['AV45_PVC_%s_HemiWM_slope_2points' % (vg)] = ''
        elif timepoint == 'Scan3':
            # try to calculate slope with 3 points
            for vg in valid_groupings:
                try:
                    diff1 = float(old_l.get('AV45_1_2_Diff'))
                    diff2 = float(old_l.get('AV45_1_3_Diff'))
                    x = [0.0, diff1, diff2]

                    val1 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Wholecereb_BL' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Wholecereb_Scan2' % vg))
                    val3 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Wholecereb_Scan3' % vg))
                    slope, intercept, r, p, stderr = linregress(x, [val1, val2, val3])
                    old_l['AV45_PVC_%s_CorticalSummary_slope_3points' % (vg)] = slope
                    if vg == 'summary':
                        old_l['AV45_PVC_summary_accumulator.0081'] = 1 if slope >= 0.0081 else 0

                    val1 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Bigref_BL' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Bigref_Scan2' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Bigref_Scan3' % vg))
                    slope, intercept, r, p, stderr = linregress(x, [val1, val2, val3])
                    old_l['AV45_PVC_%s_CorticalSummary/Bigref_slope_3points' % (vg)] = slope

                    val1 = float(old_l.get('AV45_PVC_%s_HemiWM/Wholecereb_BL' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_HemiWM/Wholecereb_Scan2' % vg))
                    val3 = float(old_l.get('AV45_PVC_%s_HemiWM/Wholecereb_Scan3' % vg))
                    slope, intercept, r, p, stderr = linregress(x, [val1, val2, val3])
                    old_l['AV45_PVC_%s_HemiWM_slope_3points' % (vg)] = slope
                except Exception as e:
                    old_l['AV45_PVC_%s_CorticalSummary_slope_3points' % (vg)] = ''
                    old_l['AV45_PVC_%s_CorticalSummary/Bigref_slope_3points' % (vg)] = ''
                    old_l['AV45_PVC_%s_HemiWM_slope_3points' % (vg)] = ''
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def syncFAQData(old_headers, old_lines, faq_file, registry_file, dump_to=None):
    tmpts = 12
    registry = importRegistry(registry_file)
    faq_by_subj = importFAQ(faq_file, registry)

    to_add_headers = []
    to_add_headers += ['FAQTOTAL.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FAQTOTAL_timePostAV45.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FAQTOTAL_AV45_6MTHS', 'FAQTOTAL_AV45_DATE',
                       'FAQTOTAL_AV45_2_6MTHS', 'FAQTOTAL_AV45_2_DATE',
                       'FAQTOTAL_AV45_3_6MTHS', 'FAQTOTAL_AV45_3_DATE',
                       'FAQTOTAL_slope']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        if len(subj_row) > tmpts:
            raise Exception("INCREASE NUMBER OF FAQ TIMEPOINTS TO > %s" % (len(subj_row)))
        new_subj_data = {}
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        
        # get closest to av45s
        six_months = 31 * 6
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToAV45(subj_row, bl_av45, av45_2, av45_3, day_limit=six_months)
        new_subj_data['FAQTOTAL_AV45_6MTHS'] = av45_bl_closest['FAQTOTAL'] if av45_bl_closest else ''
        new_subj_data['FAQTOTAL_AV45_DATE'] = av45_bl_closest['EXAMDATE'] if av45_bl_closest else ''
        new_subj_data['FAQTOTAL_AV45_2_6MTHS'] = av45_2_closest['FAQTOTAL'] if av45_2_closest else ''
        new_subj_data['FAQTOTAL_AV45_2_DATE'] = av45_2_closest['EXAMDATE'] if av45_2_closest else ''
        new_subj_data['FAQTOTAL_AV45_3_6MTHS'] = av45_3_closest['FAQTOTAL'] if av45_3_closest else ''
        new_subj_data['FAQTOTAL_AV45_3_DATE'] = av45_3_closest['EXAMDATE'] if av45_3_closest else ''
        
        post_av45_points = []
        for i in range(tmpts):
            if i < len(subj_row):
                cur_score = subj_row[i]
                date = cur_score['EXAMDATE']
                tot_score = cur_score['FAQTOTAL']
                if bl_av45 is not None:
                    timediff = (date - bl_av45).days / 365.0
                    if timediff > (-90/365.0):
                        post_av45_points.append((timediff, tot_score))
                else:
                    timediff = ''
                new_subj_data['FAQTOTAL.%s' % (i+1)] = tot_score
                new_subj_data['FAQTOTAL_timePostAV45.%s' % (i+1)] = timediff
        # get slope
        faq_slope = slope(post_av45_points)
        new_subj_data['FAQTOTAL_slope'] = faq_slope if faq_slope is not None else ''
        return new_subj_data


    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, faq_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)
    return (new_headers, new_lines)

def syncNPIData(old_headers, old_lines, npi_file, dump_to=None):
    tmpts = 5
    npi_by_subj = importNPI(npi_file)

    to_add_headers = []
    to_add_headers += ['NPITOTAL.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['NPITOTAL_timePostAV45.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['NPITOTAL_AV45_6MTHS', 'NPITOTAL_AV45_DATE',
                       'NPITOTAL_AV45_2_6MTHS', 'NPITOTAL_AV45_2_DATE',
                       'NPITOTAL_AV45_3_6MTHS', 'NPITOTAL_AV45_3_DATE',
                       'NPITOTAL_slope']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        if len(subj_row) > tmpts:
            raise Exception("INCREASE NUMBER OF NPI TIMEPOINTS TO > %s" % (len(subj_row)))

        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        new_subj_data = {}

        # get closest to av45s
        six_months = 31 * 6
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToAV45(subj_row, bl_av45, av45_2, av45_3, day_limit=six_months)
        new_subj_data['NPITOTAL_AV45_6MTHS'] = av45_bl_closest['NPITOTAL'] if av45_bl_closest else ''
        new_subj_data['NPITOTAL_AV45_DATE'] = av45_bl_closest['EXAMDATE'] if av45_bl_closest else ''
        new_subj_data['NPITOTAL_AV45_2_6MTHS'] = av45_2_closest['NPITOTAL'] if av45_2_closest else ''
        new_subj_data['NPITOTAL_AV45_2_DATE'] = av45_2_closest['EXAMDATE'] if av45_2_closest else ''
        new_subj_data['NPITOTAL_AV45_3_6MTHS'] = av45_3_closest['NPITOTAL'] if av45_3_closest else ''
        new_subj_data['NPITOTAL_AV45_3_DATE'] = av45_3_closest['EXAMDATE'] if av45_3_closest else ''
        
        post_av45_points = []
        for i in range(tmpts):
            if i < len(subj_row):
                cur_score = subj_row[i]
                date = cur_score['EXAMDATE']
                tot_score = cur_score['NPITOTAL']
                if bl_av45 is not None:
                    timediff = (date - bl_av45).days / 365.0
                    if timediff > (-90/365.0):
                        post_av45_points.append((timediff, tot_score))
                else:
                    timediff = ''
                new_subj_data['NPITOTAL.%s' % (i+1)] = tot_score
                new_subj_data['NPITOTAL_timePostAV45.%s' % (i+1)] = timediff
        # get slope
        npi_slope = slope(post_av45_points)
        new_subj_data['NPITOTAL_slope'] = npi_slope if npi_slope is not None else ''
        return new_subj_data

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, npi_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)
    return (new_headers, new_lines)

def syncMHISTData(old_headers, old_lines, mhist_file, dump_to=None):
    mhist_by_subj = importMedicalHistory(mhist_file)

    to_add_headers = ['SMOKING','DIABETES','HYPERLIPIDEMIA']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after='Handedness')
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        new_subj_data = {'SMOKING': 1 if subj_row['smoking'] else 0,
                         'DIABETES': 1 if subj_row['diabetes'] else 0,
                         'HYPERLIPIDEMIA': 1 if subj_row['hyperlipidemia'] else 0}
        return new_subj_data

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, mhist_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)
    return (new_headers, new_lines)

def syncGDData(old_headers, old_lines, gd_file, registry_file, dump_to=None):
    tmpts = 10
    registry = importRegistry(registry_file, include_all=True)
    gd_by_subj = importGD(gd_file, registry=registry)

    to_add_headers = []
    to_add_headers += ['GD_TOTAL.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['GD_timePostAV45.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['GD_AV45_6MTHS', 'GD_AV45_DATE',
                       'GD_AV45_2_6MTHS', 'GD_AV45_2_DATE',
                       'GD_AV45_3_6MTHS', 'GD_AV45_3_DATE',
                       'GD_slope']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        scores = sorted(subj_row, key=lambda x: x['EXAMDATE'])
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        new_subj_data = {}
        
        six_months = 31 * 6
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToAV45(scores, bl_av45, av45_2, av45_3, day_limit=six_months)
        if av45_bl_closest:
            new_subj_data['GD_AV45_6MTHS'] = av45_bl_closest['GDTOTAL']
            new_subj_data['GD_AV45_DATE'] = av45_bl_closest['EXAMDATE']
        else:
            new_subj_data['GD_AV45_6MTHS'] = ''
            new_subj_data['GD_AV45_DATE'] = ''
        if av45_2_closest:
            new_subj_data['GD_AV45_2_6MTHS'] = av45_2_closest['GDTOTAL']
            new_subj_data['GD_AV45_2_DATE'] = av45_2_closest['EXAMDATE']
        else:
            new_subj_data['GD_AV45_2_6MTHS'] = ''
            new_subj_data['GD_AV45_2_DATE'] = ''
        if av45_3_closest:
            new_subj_data['GD_AV45_3_6MTHS'] = av45_3_closest['GDTOTAL']
            new_subj_data['GD_AV45_3_DATE'] = av45_3_closest['EXAMDATE']
        else:
            new_subj_data['GD_AV45_3_6MTHS'] = ''
            new_subj_data['GD_AV45_3_DATE'] = ''
        
        post_av45_points = []
        for i in range(tmpts):
            if i < len(scores):
                cur_score = scores[i]
                date = cur_score['EXAMDATE']
                tot_score = cur_score['GDTOTAL']
                try:
                    tot_score = int(tot_score)
                except:
                    tot_score = ''
                if bl_av45 is not None:
                    timediff = (date - bl_av45).days / 365.0
                    if tot_score != '' and timediff > (-90/365.0):
                        post_av45_points.append((timediff, tot_score))
                else:
                    timediff = ''
                new_subj_data['GD_TOTAL.%s' % (i+1)] = tot_score
                new_subj_data['GD_timePostAV45.%s' % (i+1)] = timediff
        # get slope
        gd_slope = slope(post_av45_points)
        new_subj_data['GD_slope'] = gd_slope if gd_slope is not None else ''
        return new_subj_data

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, gd_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def syncADASCogData(old_headers, old_lines, adni1_adas_file, adnigo2_adas_file, registry_file, dump_to=None):
    tmpts=11
    registry = importRegistry(registry_file)
    adas_by_subj = importADASCog(adni1_adas_file, adnigo2_adas_file, registry=registry)

    to_add_headers = []
    to_add_headers += ['ADAScog.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIME_ADAS.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIMEreltoAV45_ADAS.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIMEpostAV45_ADAS.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['ADAS_3MTH_AV45','ADAS_3MTHS_AV45DATE',
                       'ADAS_post_AV45_followuptime','ADASslope_postAV45',
                       'ADAS_AV45_2_3MTHS','ADAS_AV45_2_DATE',
                       'ADAS_AV45_3_3MTHS','ADAS_AV45_3_DATE']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        tests = sorted(subj_row, key=lambda x : x['EXAMDATE'])
        
        dates = [t['EXAMDATE'] for t in tests]
        if len(set(dates)) != len(dates):
            print "%s has Dup dates, skipping: %s" % (subj, dates)
            return {}

        first_scan_date = dates[0]
        new_subj_data = {}
        all_slope_points = []
        post_slope_points = []
        max_followup = None
        for i in range(11):
            if i < len(tests):
                test_results = tests[i]
                test_date = test_results['EXAMDATE']
                test_score = test_results['TOTSCORE']
                if test_score != '':
                    test_score = float(test_score)
                    
                diff_from_first = (test_date-first_scan_date).days / 365.0
                if test_score != '':
                    all_slope_points.append((diff_from_first, test_score))
                if bl_av45 is not None:
                    timediff = (test_date-bl_av45).days / 365.0
                    new_subj_data['TIMEreltoAV45_ADAS.%s' % (i+1)] = timediff
                    if timediff > -93.0/365.0:
                        max_followup = timediff
                        new_subj_data['TIMEpostAV45_ADAS.%s' % (i+1)] = timediff
                        if test_score != '':
                            post_slope_points.append((timediff, test_score))
                    if abs(timediff) <= (93.0/365.0) and 'ADAS_3MTH_AV45' not in new_subj_data:
                        new_subj_data['ADAS_3MTH_AV45'] = test_score
                        new_subj_data['ADAS_3MTHS_AV45DATE'] = test_date
                if av45_2 is not None:
                    timediff = (test_date-av45_2).days / 365.0
                    if abs(timediff) <= (93.0/365.0) and 'ADAS_AV45_2_3MTHS' not in new_subj_data:
                        new_subj_data['ADAS_AV45_2_3MTHS'] = test_score
                        new_subj_data['ADAS_AV45_2_DATE'] = test_date
                if av45_3 is not None:
                    timediff = (test_date-av45_3).days / 365.0
                    if abs(timediff) <= (93.0/365.0) and 'ADAS_AV45_3_3MTHS' not in new_subj_data:
                        new_subj_data['ADAS_AV45_3_3MTHS'] = test_score
                        new_subj_data['ADAS_AV45_3_DATE'] = test_date

                new_subj_data['ADAScog_DATE%s' % (i+1)] = test_date
                new_subj_data['ADAScog.%s' % (i+1)] = test_score
                new_subj_data['TIME_ADAS.%s' % (i+1)] = diff_from_first
        new_subj_data['ADAS_post_AV45_followuptime'] = max(max_followup, 0.0) if max_followup is not None else ''

        # get slopes
        adas_slope = slope(post_slope_points)
        new_subj_data['ADASslope_postAV45'] = adas_slope if adas_slope is not None else ''
        return new_subj_data

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, adas_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def syncMMSEData(old_headers, old_lines, mmse_file, registry_file, dump_to=None):
    tmpts = 11
    registry = importRegistry(registry_file, include_all=True)
    mmse_by_subj = importMMSE(mmse_file, registry=registry)

    # add these new headers, if not present
    to_add_headers = []
    to_add_headers += ['MMSCORE.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIMEpostAV45_MMSE.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['MMSE_post_AV45_followuptime',
                       'MMSEslope_postAV45',
                       'MMSE_AV45_3MTHS','MMSE_AV45_DATE',
                       'MMSE_AV45_2_3MTHS','MMSE_AV45_2_DATE',
                       'MMSE_AV45_3_3MTHS','MMSE_AV45_3_DATE']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        tests = sorted(subj_row, key=lambda x: x[0])
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)

        new_subj_data = {}
        max_followup = None
        post_av45_points = []
        for i in range(tmpts):
            if i < len(tests):
                test_date, test_results = tests[i]
                test_score = test_results['MMSCORE']
                try:
                    test_score = int(test_score)
                    if test_score < 0:
                        test_score = ''
                except:
                    test_score = ''

                if bl_av45 is not None:
                    timediff = (test_date-bl_av45).days / 365.0
                    if abs(timediff) <= (93.0/365.0) and 'MMSE_AV45_3MTHS' not in new_subj_data:
                        new_subj_data['MMSE_AV45_3MTHS'] = test_score
                        new_subj_data['MMSE_AV45_DATE'] = test_date
                    if timediff > (-93.0/365.0):
                        max_followup = timediff
                        new_subj_data['TIMEpostAV45_MMSE.%s' % (i+1)] = timediff
                        if test_score != '':
                            post_av45_points.append((timediff, test_score))
                if av45_2 is not None:
                    timediff = (test_date-av45_2).days / 365.0
                    if abs(timediff) <= (93.0/365.0) and 'MMSE_AV45_2_3MTHS' not in new_subj_data:
                        new_subj_data['MMSE_AV45_2_3MTHS'] = test_score
                        new_subj_data['MMSE_AV45_2_DATE'] = test_date
                if av45_3 is not None:
                    timediff = (test_date-av45_3).days / 365.0
                    if abs(timediff) <= (93.0/365.0) and 'MMSE_AV45_3_3MTHS' not in new_subj_data:
                        new_subj_data['MMSE_AV45_3_3MTHS'] = test_score
                        new_subj_data['MMSE_AV45_3_DATE'] = test_date

                new_subj_data['MMSE_DATE%s' % (i+1)] = test_date
                new_subj_data['MMSCORE.%s' % (i+1)] = test_score
        new_subj_data['MMSE_post_AV45_followuptime'] = max(max_followup, 0.0) if max_followup is not None else ''

        # get post av45 slope
        mmse_slope = slope(post_av45_points)
        new_subj_data['MMSEslope_postAV45'] = mmse_slope if mmse_slope is not None else ''
        return new_subj_data

    new_lines = []
    for i, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, mmse_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return new_headers, new_lines

def syncAVLTData(old_headers, old_lines, neuro_battery_file, registry_file, dump_to=None):
    tmpts = 11
    registry = importRegistry(registry_file, include_all=True)
    avlt_by_subj = importAVLT(neuro_battery_file, registry=registry)

    # no change in headers
    to_add_headers = []
    to_add_headers += ['AVLT.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIME_AVLT.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIMEreltoAV45_AVLT.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIMEpostAV45_AVLT.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['AVLT_3MTHS_AV45','AVLT_3MTHSAV45_Date',
                       'AVLT_AV45_2_3MTHS','AVLT_AV45_2_DATE',
                       'AVLT_post_AV45_followuptime',
                       'AVLT_slope_all','AVLTslope_postAV45']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        tests = sorted(subj_row, key=lambda x : x['EXAMDATE'])
        
        dates = [t['EXAMDATE'] for t in tests]
        if len(set(dates)) != len(dates):
            print "%s has Dup dates, skipping: %s" % (subj, dates)
            return {}

        new_subj_data = {}
        all_slope_points = []
        post_slope_points = []
        max_followup = ''
        first_scan_date = dates[0]
        for i in range(tmpts):
            if i < len(tests):
                test_results = tests[i]
                test_date = test_results['EXAMDATE']
                test_score = test_results['TOTS']
                diff_from_first = (test_date-first_scan_date).days / 365.0
                diff_from_bl = (test_date - bl_av45).days / 365.0 if bl_av45 is not None else ''
                diff_from_scan2 = (test_date - av45_2).days / 365.0 if av45_2 is not None else ''
                if test_score != '':
                    all_slope_points.append((diff_from_first,test_score))
                if diff_from_bl != '':
                    max_followup = diff_from_bl
                    if test_score != '' and diff_from_bl >= (-93.0/365.0):
                        post_slope_points.append((diff_from_first,test_score))
                        new_subj_data['TIMEpostAV45_AVLT.%s' % (i+1)] = diff_from_bl
                    if abs(diff_from_bl) <= (-93.0/365.0) and 'AVLT_3MTHS_AV45' not in new_subj_data:
                        new_subj_data['AVLT_3MTHS_AV45'] = test_score
                        new_subj_data['AVLT_3MTHSAV45_Date'] = test_date
                if diff_from_scan2 != '':
                    if abs(diff_from_scan2) <= (-93.0/365.0) and 'AVLT_AV45_2_3MTHS' not in new_subj_data:
                        new_subj_data['AVLT_AV45_2_3MTHS'] = test_score
                        new_subj_data['AVLT_AV45_2_DATE'] = test_date
                new_subj_data['AVLT_DATE.%s' % (i+1)] = test_date
                new_subj_data['AVLT.%s' % (i+1)] = test_score
                new_subj_data['TIME_AVLT.%s' % (i+1)] = diff_from_first
                new_subj_data['TIMEreltoAV45_AVLT.%s' % (i+1)] = diff_from_bl
        new_subj_data['AVLT_post_AV45_followuptime'] = max(max_followup, 0.0) if max_followup is not None else ''

        # get slopes
        all_slope = slope(all_slope_points)
        post_slope = slope(post_slope_points)
        new_subj_data['AVLT_slope_all'] = all_slope if all_slope is not None else ''
        new_subj_data['AVLTslope_postAV45'] = post_slope if post_slope is not None else ''
        return new_subj_data

    new_lines = []
    new_values = 0
    total = 0
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, avlt_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)   
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def syncDiagnosisData(old_headers, old_lines, diag_file, registry_file, demog_file, arm_file, pet_meta_file, dump_to=None):
    registry = importRegistry(registry_file)
    diag_by_subj = importADNIDiagnosis(diag_file, registry=registry)
    demogs = importDemog(demog_file)
    arm = importARM(arm_file)
    pet_meta = importPetMETA(pet_meta_file)

    pivot_date = datetime(day=1, month=6, year=2015)
    pivot_date_closest_diag_key = 'Closest_DX_Jun15'
    pivot_date_closest_date_key = 'DX_Jun15_closestdate'

    to_add_headers = ['AV45_Date','AV45_2_Date','AV45_3_Date',
                      'BD MM-YY','Age@AV45','Age@AV45_2','Age@AV45_3',
                      'AV45_1_2_Diff','AV45_1_3_Diff',
                      'MCItoADConv(fromav45)','MCItoADConvDate','MCItoADconv_','AV45_MCItoAD_ConvTime','Baseline_MCItoAD_ConvTime',
                      'Diag@AV45_long','Diag@AV45_2_long','Diag@AV45_3_long','FollowupTimetoDX','Baseline','Init_Diagnosis',
                      pivot_date_closest_diag_key,pivot_date_closest_date_key]
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after='APOE4_NUM')
    wipeKeys(old_lines, to_add_headers)

    # find initial diags per subject
    for subj in diag_by_subj.keys():
        init_diag = None
        if subj in arm:
            sorted_arm = sorted(arm[subj], key=lambda x: x['USERDATE'])
            init_diags = list(set([_['STATUS'] for _ in sorted_arm]))
            init_diag = init_diags[-1]
        for point in diag_by_subj[subj]:
            point['init_diag'] = init_diag

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        subj_diags = sorted(subj_row, key=lambda x: x['EXAMDATE'])

        init_diag = subj_diags[0]['init_diag']
        if init_diag is None:
            raise Exception("NO INITIAL DIAGNOSIS FOUND FOR %s" % subj)

        # get mci type
        mci_type = 'LMCI'
        if init_diag in set(['LMCI', 'EMCI', 'SMC']):
            mci_type = init_diag

        # convert codes to diagnoses base on init diag
        for idx in range(len(subj_diags)):
            if idx == 0:
                subj_diags[idx]['diag'] = init_diag
                continue
            change = subj_diags[idx]['change']
            if change in set([1,2,3]): # stable
                subj_diags[idx]['diag'] = subj_diags[idx-1]['diag']
            elif change in set([5,6]): # conversion to AD
                subj_diags[idx]['diag'] = 'AD'
            elif change in set([4]): # conversion to MCI
                subj_diags[idx]['diag'] = mci_type
            elif change in set([8]): # reversion to MCI
                subj_diags[idx]['diag'] = mci_type
            elif change in set([7,9]): # reversion to normal
                subj_diags[idx]['diag'] = 'N'

        # Get AV45 Scan dates
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        patient_dob = demogs.get(subj, {}).get('dob',None)
        av45_age = ''
        av45_age2 = ''
        av45_age3 = ''
        if patient_dob is not None:
            if bl_av45:
                av45_age = (bl_av45-patient_dob).days / 365.25
            if av45_2:
                av45_age2 = (av45_2-patient_dob).days / 365.25
            if av45_3:
                av45_age3 = (av45_3-patient_dob).days / 365.25
        
        # Date differences between scans
        av45_1_2_diff = ((av45_2 - bl_av45).days/365.0) if (bl_av45 is not None and av45_2 is not None) else ''
        av45_1_3_diff = ((av45_3 - bl_av45).days/365.0) if (bl_av45 is not None and av45_3 is not None) else ''

        # find very first visit date in registry (excluding scmri)
        subj_registry = [_ for _ in registry[subj] if _['VISCODE'] != 'scmri' and _['VISCODE2'] != 'scmri']
        sorted_registry = sorted(subj_registry, key=lambda x: x['EXAMDATE'])
        baseline_registry = sorted_registry[0]

        # find pivot data point
        sorted_by_pivot = sorted(subj_diags, key=lambda x: abs(x['EXAMDATE'] - pivot_date))
        closest_to_pivot = sorted_by_pivot[0]

        new_data = {'AV45_Date': bl_av45,
                    'AV45_2_Date': av45_2,
                    'AV45_3_Date': av45_3,
                    'BD MM-YY': patient_dob,
                    'Age@AV45': av45_age,
                    'Age@AV45_2' : av45_age2,
                    'Age@AV45_3' : av45_age3,
                    'AV45_1_2_Diff': av45_1_2_diff,
                    'AV45_1_3_Diff': av45_1_3_diff,
                    'MCItoADConv(fromav45)': '0',
                    'MCItoADConvDate': '',
                    'MCItoADconv_': '',
                    'AV45_MCItoAD_ConvTime': '',
                    'Baseline_MCItoAD_ConvTime': '',
                    'Diag@AV45_long': '',
                    'Diag@AV45_2_long': '',
                    'Diag@AV45_3_long': '',
                    'FollowupTimetoDX': (pivot_date - baseline_registry['EXAMDATE']).days / 365.0,
                    'Baseline': baseline_registry['EXAMDATE'],
                    'Init_Diagnosis': init_diag,
                    pivot_date_closest_diag_key: closest_to_pivot['diag'],
                    pivot_date_closest_date_key: closest_to_pivot['EXAMDATE']}
        
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToAV45(subj_diags, bl_av45, av45_2, av45_3)
        if av45_bl_closest:
            new_data['Diag@AV45_long'] = av45_bl_closest['diag']
        if av45_2_closest:
            new_data['Diag@AV45_2_long'] = av45_2_closest['diag']
        if av45_3_closest:
            new_data['Diag@AV45_3_long'] = av45_3_closest['diag'] 

        for i, diag_row in enumerate(subj_diags):
            # parse change based on init diag
            if diag_row['change'] == 5:
                new_data['MCItoADConv(fromav45)'] = '1'
                new_data['MCItoADConvDate'] = new_data['MCItoADconv_'] = diag_row['EXAMDATE']
                new_data['Baseline_MCItoAD_ConvTime'] = (diag_row['EXAMDATE'] - new_data['Baseline']).days / 365.0
                if bl_av45 is not None:
                    new_data['AV45_MCItoAD_ConvTime'] = (diag_row['EXAMDATE'] - bl_av45).days / 365.0
        return new_data

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, diag_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=pet_meta)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def syncAV45Data(old_headers, old_lines, av45_file, av45_nontp_file, registry_file, dump_to=None):
    '''
    This function does not follow the pattern of the other sync functions because
    the header/data update is accomplished in a nested function, and it also allows
    for the possibility of introducing new subjects via the dataset being synced in
    '''
    registry = importRegistry(registry_file)
    av45_by_subj = importAV45(av45_file, av45_nontp_file, registry=registry)
    print "SUBJECTS WITH AV45: %s" % len(av45_by_subj)
    new_headers = None
    new_lines = []
    old_subjects = set()
    clear_av45 = set()
    for old_l in old_lines:
        # get subject ID
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            continue 
        old_subjects.add(subj)

        updated_headers, new_data = parseAV45Entries(old_headers, av45_by_subj.get(subj,[]))
        if new_headers is None:
            new_headers = updated_headers
        new_data = convertToCSVDataType(new_data, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # find new subjects - in av45 file but not in old master csv file
    print len(set(av45_by_subj.keys()))
    print len(old_subjects)
    new_subjects = list(set(av45_by_subj.keys()) - old_subjects)
    print new_subjects
    for ns in new_subjects:
        new_subj_row = {k: '' for k in new_headers}
        new_subj_row['RID'] = str(ns)
        updated_headers, new_columns = parseAV45Entries(old_headers, av45_by_subj[ns])
        new_subj_row.update(new_columns)
        new_lines.append(new_subj_row)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def parseAV45Entries(old_headers, subj_av45):
    subj_av45 = sorted(subj_av45, key=lambda x: x['EXAMDATE'])
    exam_times = [_['EXAMDATE'] for _ in subj_av45]
    exam_timedeltas = [(_-exam_times[0]).days / 365.0 for _ in exam_times]

    wm70_composite_keys = ['AV45_WM70/composite', 'AV45_2_WM70/composite', 'AV45_3_WM70/composite']
    wm70_cerebg_keys = ['AV45_WM70/cerebg', 'AV45_2_WM70/cerebg', 'AV45_3_WM70/cerebg']
    wm70_wcereb_keys = ['AV45_WM70/wcereb', 'AV45_2_WM70/wcereb', 'AV45_3_WM70/wcereb']
    unilateral_keys = ['AV45_LeftPutamen/WM70','AV45_RightPutamen/WM70',
                       'AV45_LeftCaudate/WM70','AV45_RightCaudate/WM70',
                       'AV45_LeftPallidum/WM70','AV45_RightPallidum/WM70']
    bigref_keys = ['AV45_BigRef','AV45_2_BigRef','AV45_3_BigRef'] # assuming composite ROI
    wm70_keys = ['AV45_WM70','AV45_2_WM70','AV45_3_WM70'] # assuming composite ROI
    cerebg_keys = ['AV45_cerebg','AV45_2_cerebg','AV45_3_cerebg'] # assuming composite ROI
    wcereb_keys = ['AV45_wcereb','AV45_2_wcereb','AV45_3_wcereb'] # assuming composite ROI
    brainstem_keys = ['AV45_brainstem','AV45_2_brainstem','AV45_3_brainstem'] # assuming composite ROI
    wmratio_keys = ['AV45_WMratio','AV45_2_WMratio','AV45_3_WMratio']
    frontal_bigref_keys = ['AV45_Frontal/BigRef','AV45_2_Frontal/BigRef','AV45_3_Frontal/BigRef']
    cingulate_bigref_keys = ['AV45_Cingulate/BigRef','AV45_2_Cingulate/BigRef','AV45_3_Cingulate/BigRef']
    parietal_bigref_keys = ['AV45_Parietal/BigRef','AV45_2_Parietal/BigRef','AV45_3_Parietal/BigRef']
    temporal_bigref_keys = ['AV45_Temporal/BigRef','AV45_2_Temporal/BigRef','AV45_3_Temporal/BigRef']

    # generate additional keys and arrange into header list
    all_tp_keys = ['AV45_TP_SPECIFIC_BL', 'AV45_TP_SPECIFIC_SCAN2', 'AV45_TP_SPECIFIC_SCAN3']
    all_wm70_composite_keys = wm70_composite_keys + \
                         ["%s_pchange" % _ for _ in wm70_composite_keys[1:]] + \
                         ["%s_pchange_ABS" % _ for _ in wm70_composite_keys[1:]] + \
                         ["%s_diff" % _ for _ in wm70_composite_keys[1:]] + \
                         ["%s_diff_ABS" % _ for _ in wm70_composite_keys[1:]] + \
                         ['AV45_WM70/composite_Slope_2pts', 'AV45_WM70/composite_Slope_3pts']
    all_wm70_cerebg_keys = wm70_cerebg_keys + \
                           ["%s_pchange" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ["%s_pchange_ABS" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ["%s_diff" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ["%s_diff_ABS" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ['AV45_WM70/cerebg_Slope_2pts', 'AV45_WM70/cerebg_Slope_3pts']
    all_wm70_wcereb_keys = wm70_wcereb_keys + \
                           ["%s_pchange" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ["%s_pchange_ABS" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ["%s_diff" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ["%s_diff_ABS" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ['AV45_WM70/wcereb_Slope_2pts', 'AV45_WM70/wcereb_Slope_3pts']
    all_unilateral_keys = unilateral_keys
    all_bigref_keys = bigref_keys + \
                      ["%s_pchange" % _ for _ in bigref_keys[1:]] + \
                      ["%s_pchange_ABS" % _ for _ in bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in bigref_keys[1:]] + \
                      ["%s_diff_ABS" % _ for _ in bigref_keys[1:]] + \
                      ["%s_BIN.79" % _ for _ in bigref_keys] + \
                      ['AV45_BigRef_Slope_2pts', 'AV45_BigRef_Slope_3pts']
    all_wm70_keys = wm70_keys + \
                    ["%s_pchange" % _ for _ in wm70_keys[1:]] + \
                    ["%s_pchange_ABS" % _ for _ in wm70_keys[1:]] + \
                    ["%s_diff" % _ for _ in wm70_keys[1:]] + \
                    ["%s_diff_ABS" % _ for _ in wm70_keys[1:]] + \
                    ["%s_BIN.62" % _ for _ in wm70_keys] + \
                    ['AV45_WM70_Slope_2pts', 'AV45_WM70_Slope_3pts']
    all_cerebg_keys = cerebg_keys + \
                    ["%s_pchange" % _ for _ in cerebg_keys[1:]] + \
                    ["%s_pchange_ABS" % _ for _ in cerebg_keys[1:]] + \
                    ["%s_diff" % _ for _ in cerebg_keys[1:]] + \
                    ["%s_diff_ABS" % _ for _ in cerebg_keys[1:]] + \
                    ["%s_BIN1.26" % _ for _ in cerebg_keys] + \
                    ['AV45_cerebg_Slope_2pts', 'AV45_cerebg_Slope_3pts']
    all_wcereb_keys = wcereb_keys + \
                    ["%s_pchange" % _ for _ in wcereb_keys[1:]] + \
                    ["%s_pchange_ABS" % _ for _ in wcereb_keys[1:]] + \
                    ["%s_diff" % _ for _ in wcereb_keys[1:]] + \
                    ["%s_diff_ABS" % _ for _ in wcereb_keys[1:]] + \
                    ["%s_BIN1.11" % _ for _ in wcereb_keys] + \
                    ['AV45_wcereb_Slope_2pts', 'AV45_wcereb_Slope_3pts']
    all_brainstem_keys = brainstem_keys + \
                    ["%s_pchange" % _ for _ in brainstem_keys[1:]] + \
                    ["%s_pchange_ABS" % _ for _ in brainstem_keys[1:]] + \
                    ["%s_diff" % _ for _ in brainstem_keys[1:]] + \
                    ["%s_diff_ABS" % _ for _ in brainstem_keys[1:]] + \
                    ["%s_BIN.79" % _ for _ in brainstem_keys] + \
                    ['AV45_brainstem_Slope_2pts', 'AV45_brainstem_Slope_3pts']
    all_wmratio_keys = wmratio_keys + \
                      ["%s_pchange" % _ for _ in wmratio_keys[1:]] + \
                      ["%s_pchange_ABS" % _ for _ in wmratio_keys[1:]] + \
                      ["%s_diff" % _ for _ in wmratio_keys[1:]] + \
                      ["%s_diff_ABS" % _ for _ in wmratio_keys[1:]] + \
                      ['AV45_WMratio_Slope_2pts', 'AV45_WMratio_Slope_3pts']
    all_frontal_bigref_keys = frontal_bigref_keys + \
                      ["%s_pchange" % _ for _ in frontal_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in frontal_bigref_keys[1:]] + \
                      ['AV45_Frontal/BigRef_Slope_2pts', 'AV45_Frontal/BigRef_Slope_3pts']
    all_cingulate_bigref_keys = cingulate_bigref_keys + \
                      ["%s_pchange" % _ for _ in cingulate_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in cingulate_bigref_keys[1:]] + \
                      ['AV45_Cingulate/BigRef_Slope_2pts', 'AV45_Cingulate/BigRef_Slope_3pts']
    all_parietal_bigref_keys = parietal_bigref_keys + \
                      ["%s_pchange" % _ for _ in parietal_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in parietal_bigref_keys[1:]] + \
                      ['AV45_Parietal/BigRef_Slope_2pts', 'AV45_Parietal/BigRef_Slope_3pts']
    all_temporal_bigref_keys = temporal_bigref_keys + \
                      ["%s_pchange" % _ for _ in temporal_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in temporal_bigref_keys[1:]] + \
                      ['AV45_Temporal/BigRef_Slope_2pts', 'AV45_Temporal/BigRef_Slope_3pts']
    all_av45_key_lists = [all_tp_keys,
                          all_wm70_composite_keys,
                          all_wm70_cerebg_keys,
                          all_wm70_wcereb_keys,
                          all_unilateral_keys,
                          all_bigref_keys,
                          all_wm70_keys,
                          all_cerebg_keys,
                          all_wcereb_keys,
                          all_brainstem_keys,
                          all_wmratio_keys,
                          all_frontal_bigref_keys,
                          all_cingulate_bigref_keys,
                          all_parietal_bigref_keys,
                          all_temporal_bigref_keys]
    pchange_diff_lists = [all_wm70_composite_keys,
                          all_wm70_cerebg_keys,
                          all_wm70_wcereb_keys,
                          all_bigref_keys,
                          all_wm70_keys,
                          all_cerebg_keys,
                          all_wcereb_keys,
                          all_brainstem_keys,
                          all_wmratio_keys,
                          all_frontal_bigref_keys,
                          all_cingulate_bigref_keys,
                          all_parietal_bigref_keys,
                          all_temporal_bigref_keys]
    abs_lists = [all_wm70_composite_keys,
                 all_wm70_cerebg_keys,
                 all_wm70_wcereb_keys,
                 all_bigref_keys,
                 all_wm70_keys,
                 all_cerebg_keys,
                 all_wcereb_keys,
                 all_brainstem_keys,
                 all_wmratio_keys]
    all_av45_keys = [_ for l in all_av45_key_lists for _ in l]
    new_headers = rearrangeHeaders(old_headers, all_av45_keys, after='AV45_1_3_Diff')

    data = {k:'' for k in all_av45_keys}
    wm70_0 = None
    # fill in values
    for i, point in enumerate(subj_av45):
        # extract necessary values
        cerebw = float(point['CEREBELLUMWHITEMATTER'])
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
        tp_spec = 1 if point['TP_SPECIFIC'] else 0

        # fill in basic keys
        data[wm70_composite_keys[i]] = wm70/compositeroi
        data[wm70_cerebg_keys[i]] = wm70/cerebg
        data[wm70_wcereb_keys[i]] = wm70/wcereb
        data[bigref_keys[i]] = compositeroi/bigref
        data[wm70_keys[i]] = compositeroi/wm70
        data[cerebg_keys[i]] = compositeroi/cerebg
        data[wcereb_keys[i]] = compositeroi/wcereb
        data[brainstem_keys[i]] = compositeroi/brainstem
        data[frontal_bigref_keys[i]] = frontal/bigref
        data[cingulate_bigref_keys[i]] = cingulate/bigref
        data[parietal_bigref_keys[i]] = parietal/bigref
        data[temporal_bigref_keys[i]] = temporal/bigref
        data["%s_BIN.79" % bigref_keys[i]] = 1 if (compositeroi/bigref) > 0.79 else 0
        data["%s_BIN.62" % wm70_keys[i]] = 1 if (compositeroi/wm70) > 0.62 else 0
        data["%s_BIN1.26" % cerebg_keys[i]] = 1 if (compositeroi/cerebg) > 1.26 else 0
        data["%s_BIN1.11" % wcereb_keys[i]] = 1 if (compositeroi/wcereb) > 1.11 else 0
        data["%s_BIN.79" % brainstem_keys[i]] = 1 if (compositeroi/brainstem) > 0.79 else 0


        # fill in derivative keys
        if i == 0:
            # fill in unilateral
            data['AV45_LeftPutamen/WM70'] = leftputamen/wm70
            data['AV45_RightPutamen/WM70'] = rightputamen/wm70
            data['AV45_LeftCaudate/WM70'] = leftcaudate/wm70
            data['AV45_RightCaudate/WM70'] = rightcaudate/wm70
            data['AV45_LeftPallidum/WM70'] = leftpallidum/wm70
            data['AV45_RightPallidum/WM70'] = rightpallidum/wm70
            data[wmratio_keys[i]] = (compositeroi/wcereb)
            data['AV45_TP_SPECIFIC_BL'] = tp_spec
            wm70_0 = wm70
        elif i == 1 or i == 2:
            data[wmratio_keys[i]] = (compositeroi/wcereb) / (wm70/wm70_0)
            for pl in pchange_diff_lists:
                diff = (data[pl[i]] - data[pl[0]]) / exam_timedeltas[i] # annualized
                data["%s_pchange" % pl[i]] = diff / data[pl[0]]
                data["%s_diff" % pl[i]] = diff
            for al in abs_lists:
                data["%s_pchange_ABS" % al[i]] = abs(data["%s_pchange" % al[i]])
                data["%s_diff_ABS" % al[i]] = abs(data["%s_diff" % al[i]])
            if i == 1:
                data['AV45_TP_SPECIFIC_SCAN2'] = tp_spec
                times = exam_timedeltas[:2]
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in bigref_keys[:2]])
                data['AV45_BigRef_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in wm70_keys[:2]])
                data['AV45_WM70_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in cerebg_keys[:2]])
                data['AV45_cerebg_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in wcereb_keys[:2]])
                data['AV45_wcereb_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in brainstem_keys[:2]])
                data['AV45_brainstem_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in wm70_composite_keys[:2]])
                data['AV45_WM70/composite_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in wm70_cerebg_keys[:2]])
                data['AV45_WM70/cerebg_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in wm70_wcereb_keys[:2]])
                data['AV45_WM70/wcereb_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in wmratio_keys[:2]])
                data['AV45_WMratio_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in frontal_bigref_keys[:2]])
                data['AV45_Frontal/BigRef_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in cingulate_bigref_keys[:2]])
                data['AV45_Cingulate/BigRef_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in parietal_bigref_keys[:2]])
                data['AV45_Parietal/BigRef_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in temporal_bigref_keys[:2]])
                data['AV45_Temporal/BigRef_Slope_2pts'] = slope
            if i == 2:
                data['AV45_TP_SPECIFIC_SCAN3'] = tp_spec
                times = exam_timedeltas[:3]
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in bigref_keys[:3]])
                data['AV45_BigRef_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in wm70_keys[:3]])
                data['AV45_WM70_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in cerebg_keys[:3]])
                data['AV45_cerebg_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in wcereb_keys[:3]])
                data['AV45_wcereb_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in brainstem_keys[:3]])
                data['AV45_brainstem_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in wm70_composite_keys[:3]])
                data['AV45_WM70/composite_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in wm70_cerebg_keys[:3]])
                data['AV45_WM70/cerebg_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in wm70_wcereb_keys[:3]])
                data['AV45_WM70/wcereb_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in wmratio_keys[:3]])
                data['AV45_WMratio_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in frontal_bigref_keys[:3]])
                data['AV45_Frontal/BigRef_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in cingulate_bigref_keys[:3]])
                data['AV45_Cingulate/BigRef_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in parietal_bigref_keys[:3]])
                data['AV45_Parietal/BigRef_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = linregress(times, [data[_] for _ in temporal_bigref_keys[:3]])
                data['AV45_Temporal/BigRef_Slope_3pts'] = slope
    return (new_headers, data)

def syncFDGData(old_headers, old_lines, fdg_file, registry_file, dump_to=None):
    fdg_by_subj = importFDG(fdg_file)

    to_add_headers = []
    to_add_headers += ['FDG_pons_Vis%s' % (i+1) for i in range(10)]
    to_add_headers += ['FDGVis%s_ReltoAV45' % (i+1) for i in range(10)]
    to_add_headers += ['FDGVis%s_Date' % (i+1) for i in range(10)]
    to_add_headers += ['FDG_postAV45_slope','FDG_postAV45_followuptime']
    to_add_headers += ['FDG_PONS_AV45_6MTHS', 'FDG_PONS_AV45_DATE',
                       'FDG_PONS_AV45_2_6MTHS', 'FDG_PONS_AV45_2_DATE',
                       'FDG_PONS_AV45_3_6MTHS', 'FDG_PONS_AV45_3_DATE']
    new_headers= rearrangeHeaders(old_headers, to_add_headers, after='FDG_Bin_Baseline')
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        # collapse subj rows per visit
        av45_bl, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        new_data = {}

        examdates = list(set([_['EXAMDATE'] for _ in subj_row]))
        examdates = sorted(examdates)
        subj_fdg = []
        for ed in examdates:
            rel_rows = [_ for _ in subj_row if _['EXAMDATE'] == ed]
            means = [_['MEAN'] for _ in rel_rows]
            pons_avg = np.mean(means)
            subj_fdg.append({'EXAMDATE': ed, 'PONS': pons_avg})

        # get closest to av45s
        six_months = 31 * 6
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToAV45(subj_fdg, av45_bl, av45_2, av45_3, day_limit=six_months)
        new_data['FDG_PONS_AV45_6MTHS'] = av45_bl_closest['PONS'] if av45_bl_closest else ''
        new_data['FDG_PONS_AV45_DATE'] = av45_bl_closest['EXAMDATE'] if av45_bl_closest else ''
        new_data['FDG_PONS_AV45_2_6MTHS'] = av45_2_closest['PONS'] if av45_2_closest else ''
        new_data['FDG_PONS_AV45_2_DATE'] = av45_2_closest['EXAMDATE'] if av45_2_closest else ''
        new_data['FDG_PONS_AV45_3_6MTHS'] = av45_3_closest['PONS'] if av45_3_closest else ''
        new_data['FDG_PONS_AV45_3_DATE'] = av45_3_closest['EXAMDATE'] if av45_3_closest else ''

        slope_points = []
        for i in range(10):
            if i < len(subj_fdg):
                datapoint = subj_fdg[i]
                new_data['FDG_pons_Vis%s' % (i+1)] = datapoint['PONS']
                new_data['FDGVis%s_Date' % (i+1)] = datapoint['EXAMDATE']
                if av45_bl:
                    rel_time = (datapoint['EXAMDATE'] - av45_bl).days / 365.0
                    new_data['FDGVis%s_ReltoAV45' % (i+1)] = rel_time
                    if rel_time >= -(30.0/365.0):
                        slope_points.append((rel_time, datapoint['PONS']))
                else:
                    new_data['FDGVis%s_ReltoAV45' % (i+1)] = ''
            else:
                new_data['FDG_pons_Vis%s' % (i+1)] = ''
                new_data['FDGVis%s_ReltoAV45' % (i+1)] = ''
        if len(slope_points) >= 2:
            raw_dates = [_[0] for _ in slope_points]
            raw_scores = [_[1] for _ in slope_points]
            slope, intercept, r, p, stderr = linregress(raw_dates, raw_scores)
            new_data['FDG_postAV45_slope'] = slope
            new_data['FDG_postAV45_followuptime'] = raw_dates[-1]
        else:
            new_data['FDG_postAV45_slope'] = ''
            new_data['FDG_postAV45_followuptime'] = ''

        return new_data
    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, fdg_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def syncTBMSynData(old_headers, old_lines, tbm_file, registry_file, dump_to=None):
    tbm_by_subj = importTBMSyn(tbm_file)

    # add new headers as needed
    to_add_headers = ['TBMSyn_DATE.%s' % (i+1) for i in range(10)]
    to_add_headers.extend(['TBMSyn_SCORE.%s' % (i+1) for i in range(10)])
    to_add_headers.extend(['TBMSyn_postAV45.%s' % (i+1) for i in range(10)])
    to_add_headers.append('TBMSyn_BL_DATE')
    to_add_headers.append('TBMSyn_count')
    #to_add_headers.append('TBMSyn_SLOPE')
    #to_add_headers.append('TBMSyn_slope_followuptime')
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        subj_tbm = sorted(subj_row, key=lambda x: x['EXAMDATE'])
        av45_bl, av45_2, av45_3 = getAV45Dates(old_l)
        bl_examdate = subj_tbm[0]['EXAMDATEBL']
        if av45_bl is None:
            print "No Baseline AV45 date for %s" % subj
            return {}
        # check that TBMSyn baseline is within range of AV45 baseline
        if abs(av45_bl - bl_examdate).days > 210:
            print "TBM BL and AV45 BL %s days apart for subj %s" % (abs(av45_bl - bl_examdate).days, subj)
            return {}

        new_data = {}
        slope_points = []
        for i in range(10):
            if i < len(subj_tbm):
                datapoint = subj_tbm[i]
                new_data['TBMSyn_DATE.%s' % (i+1)] = datapoint['EXAMDATE']
                new_data['TBMSyn_SCORE.%s' % (i+1)] = datapoint['TBMSYNSCOR']
                timediff = (datapoint['EXAMDATE']-datapoint['EXAMDATEBL']).days / 365.0
                slope_points.append((timediff,datapoint['TBMSYNSCOR']))
                new_data['TBMSyn_postAV45.%s' % (i+1)] = timediff
            else:
                new_data['TBMSyn_DATE.%s' % (i+1)] = ''
                new_data['TBMSyn_SCORE.%s' % (i+1)] = ''
                new_data['TBMSyn_postAV45.%s' % (i+1)] = ''
        '''
        # calculate slope
        if len(slope_points) >= 2:
            raw_dates = [_[0] for _ in slope_points]
            raw_scores = [_[1] for _ in slope_points]
            diff_dates = [0.0] + raw_dates
            diff_scores = [0.0] + raw_scores
            #slope, intercept, r, p, stderr = linregress(diff_dates, diff_scores)
            output = optimize.fmin(lambda b, x, y: ((b*x-y)**2).sum(), x0=0.1, args=(diff_dates, diff_scores))
            new_data['TBMSyn_SLOPE'] = output[0]
            new_data['TBMSyn_slope_followuptime'] = raw_dates[-1]
        else:
            new_data['TBMSyn_SLOPE'] = ''
            new_data['TBMSyn_slope_followuptime'] = ''
        '''
        new_data['TBMSyn_count'] = len(slope_points)
        new_data['TBMSyn_BL_DATE'] = bl_examdate
        return new_data

    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, tbm_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def syncCSFData(old_headers, old_lines, csf_files, registry_file, dump_to=None):
    tmpts = 7
    registry = importRegistry(registry_file)
    csf_by_subj = importCSF(csf_files, registry)
    
    # add new headers as needed
    to_add_headers = []
    to_add_headers += ['CSF_ABETA.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['CSF_ABETApostAV45.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['CSF_ABETA_slope', 'CSF_ABETA_closest_AV45', 
                       'CSF_ABETA_closest_AV45_2', 'CSF_ABETA_closest_AV45_3',
                       'CSF_ABETA_closest_AV45_BIN_192',
                       'CSF_ABETA_closest_AV45_2_BIN_192',
                       'CSF_ABETA_closest_AV45_3_BIN_192']
    to_add_headers += ['CSF_TAU.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['CSF_TAUpostAV45.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['CSF_TAU_slope', 'CSF_TAU_closest_AV45', 
                       'CSF_TAU_closest_AV45_2', 'CSF_TAU_closest_AV45_3',
                       'CSF_TAU_closest_AV45_BIN_93',
                       'CSF_TAU_closest_AV45_2_BIN_93',
                       'CSF_TAU_closest_AV45_3_BIN_93']
    to_add_headers += ['CSF_PTAU.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['CSF_PTAUpostAV45.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['CSF_PTAU_slope', 'CSF_PTAU_closest_AV45', 
                       'CSF_PTAU_closest_AV45_2', 'CSF_PTAU_closest_AV45_3',
                       'CSF_PTAU_closest_AV45_BIN_23',
                       'CSF_PTAU_closest_AV45_2_BIN_23',
                       'CSF_PTAU_closest_AV45_3_BIN_23']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    # get dates for each subject csf
    for subj in csf_by_subj.keys():
        subj_csf = csf_by_subj[subj]
        subj_csf_filtered = []
        # find the examdate for each point, then sort
        subj_listings = registry[subj]

        for i in range(len(subj_csf)):
            csf_row = subj_csf[i]
            viscode = csf_row['vc']
            viscode2 = csf_row['vc2']
            # check if exam date already in there
            if csf_row['EXAMDATE'] is not None:
                subj_csf_filtered.append(csf_row)
                continue
            # if not, get date from registry listings
            examdate = None
            for listing in subj_listings:
                if listing['VISCODE'] == viscode or listing['VISCODE2'] == viscode2:
                    examdate = listing['EXAMDATE']
                    break
            if examdate is None:
                print "No exam date for %s (%s, %s)" % (subj, viscode, viscode2)
                continue
            else:
                csf_row['EXAMDATE'] = examdate
                subj_csf_filtered.append(csf_row)
        subj_csf = sorted(subj_csf_filtered, key=lambda x: x['EXAMDATE'])
        csf_by_subj[subj] = subj_csf

    def extraction_fn(subj, subj_csf, old_l, patient_pets):
        # find av45 dates
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)

        new_data = {}
        ptau_slope_points = []
        abeta_slope_points = []
        tau_slope_points = []
        for i in range(tmpts):
            if i < len(subj_csf):
                datapoint = subj_csf[i]
                examdate = datapoint['EXAMDATE']
                timediff = ((examdate-bl_av45).days / 365.0) if bl_av45 else ''
                if timediff <= -(90.0/365.0):
                    timediff = ''
                abeta_val = datapoint['abeta']
                ptau_val = datapoint['ptau']
                tau_val = datapoint['tau']
                abeta_val = float(abeta_val) if abeta_val is not None else ''
                ptau_val = float(ptau_val) if ptau_val is not None else ''
                tau_val = float(tau_val) if tau_val is not None else ''
                new_data['CSF_ABETA.%s' % (i+1)] = abeta_val
                new_data['CSF_PTAU.%s' % (i+1)] = ptau_val
                new_data['CSF_TAU.%s' % (i+1)] = tau_val
                new_data['CSF_ABETApostAV45.%s' % (i+1)] = timediff
                new_data['CSF_PTAUpostAV45.%s' % (i+1)] = timediff
                new_data['CSF_TAUpostAV45.%s' % (i+1)] = timediff
                if timediff:
                    if abeta_val != '':
                        abeta_slope_points.append((timediff, abeta_val))
                    if tau_val != '':
                        tau_slope_points.append((timediff, tau_val))
                    if ptau_val != '':
                        ptau_slope_points.append((timediff, ptau_val))

        # match up with av45 scans
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToAV45(subj_csf, bl_av45, av45_2, av45_3, day_limit=93)
        if av45_bl_closest:
            tau_val = av45_bl_closest['tau']
            ptau_val = av45_bl_closest['ptau']
            abeta_val = av45_bl_closest['abeta']
            abeta_val = float(abeta_val) if abeta_val is not None else ''
            ptau_val = float(ptau_val) if ptau_val is not None else ''
            tau_val = float(tau_val) if tau_val is not None else ''
            new_data['CSF_TAU_closest_AV45'] = tau_val
            new_data['CSF_PTAU_closest_AV45'] = ptau_val
            new_data['CSF_ABETA_closest_AV45'] = abeta_val
            new_data['CSF_TAU_closest_AV45_BIN_93'] = 1 if tau_val >= 93 else 0
            new_data['CSF_PTAU_closest_AV45_BIN_23'] = 1 if ptau_val >= 23 else 0
            new_data['CSF_ABETA_closest_AV45_BIN_192'] = 1 if abeta_val <= 192 else 0
        if av45_2_closest:
            tau_val = av45_2_closest['tau']
            ptau_val = av45_2_closest['ptau']
            abeta_val = av45_2_closest['abeta']
            abeta_val = float(abeta_val) if abeta_val is not None else ''
            ptau_val = float(ptau_val) if ptau_val is not None else ''
            tau_val = float(tau_val) if tau_val is not None else ''
            new_data['CSF_TAU_closest_AV45_2'] = tau_val
            new_data['CSF_PTAU_closest_AV45_2'] = ptau_val
            new_data['CSF_ABETA_closest_AV45_2'] = abeta_val
            new_data['CSF_TAU_closest_AV45_2_BIN_93'] = 1 if tau_val >= 93 else 0
            new_data['CSF_PTAU_closest_AV45_2_BIN_23'] = 1 if ptau_val >= 23 else 0
            new_data['CSF_ABETA_closest_AV45_2_BIN_192'] = 1 if abeta_val <= 192 else 0
        if av45_3_closest:
            tau_val = av45_3_closest['tau']
            ptau_val = av45_3_closest['ptau']
            abeta_val = av45_3_closest['abeta']
            abeta_val = float(abeta_val) if abeta_val is not None else ''
            ptau_val = float(ptau_val) if ptau_val is not None else ''
            tau_val = float(tau_val) if tau_val is not None else ''
            new_data['CSF_TAU_closest_AV45_3'] = tau_val
            new_data['CSF_PTAU_closest_AV45_3'] = ptau_val
            new_data['CSF_ABETA_closest_AV45_3'] = abeta_val
            new_data['CSF_TAU_closest_AV45_3_BIN_93'] = 1 if tau_val >= 93 else 0
            new_data['CSF_PTAU_closest_AV45_3_BIN_23'] = 1 if ptau_val >= 23 else 0
            new_data['CSF_ABETA_closest_AV45_3_BIN_192'] = 1 if abeta_val <= 192 else 0

        # calculate slope
        if len(ptau_slope_points) >= 2:
            ptau_slope = slope(ptau_slope_points)
            new_data['CSF_PTAU_slope'] = ptau_slope if ptau_slope is not None else ''
        if len(tau_slope_points) >= 2:
            tau_slope = slope(tau_slope_points)
            new_data['CSF_TAU_slope'] = tau_slope if tau_slope is not None else ''
        if len(abeta_slope_points) >= 2:
            abeta_slope = slope(abeta_slope_points)
            new_data['CSF_ABETA_slope'] = abeta_slope if abeta_slope is not None else ''
        return new_data


    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, csf_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def syncUWData(old_headers, old_lines, uw_file, registry_file, dump_to=None):
    tmpts=12
    registry = importRegistry(registry_file)
    uw_by_subj = importUW(uw_file, registry=registry)

    to_add_headers = []
    to_add_headers += ['UW_MEM_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['UW_MEM_postAV45_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['UW_MEM_postAV45_count', 'UW_MEM_slope', 'UW_MEM_BL_3months']
    to_add_headers += ['UW_EF_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['UW_EF_postAV45_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['UW_EF_postAV45_count', 'UW_EF_slope', 'UW_EF_BL_3months']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        subj_uw = sorted(subj_row, key=lambda x: x['EXAMDATE'])
        if len(subj_uw) > tmpts:
            raise Exception("INCREASE NUMBER OF UW TIMEPOINTS TO > %s" % (len(subj_uw)))
        
        new_data = {}
        mem_slope_points = []
        ef_slope_points = []
        for i in range(tmpts):
            if i < len(subj_uw):
                datapoint = subj_uw[i]
                examdate = datapoint['EXAMDATE']
                cur_mem = datapoint['ADNI_MEM']
                cur_ef = datapoint['ADNI_EF']
                timediff = ((examdate-bl_av45).days / 365.0) if bl_av45 else ''
                if timediff <= -(90.0/365.0):
                    timediff = ''
                if not np.isnan(cur_mem):
                    new_data['UW_MEM_%s' % (i+1)] = cur_mem
                    new_data['UW_MEM_postAV45_%s' % (i+1)] = timediff
                    if timediff != '':
                        mem_slope_points.append((timediff, cur_mem))
                        if abs(timediff) <= (90.0/365.0) and 'UW_MEM_BL_3months' not in new_data:
                            new_data['UW_MEM_BL_3months'] = cur_mem
                if not np.isnan(cur_ef):
                    new_data['UW_EF_%s' % (i+1)] = cur_ef
                    new_data['UW_EF_postAV45_%s' % (i+1)] = timediff
                    if timediff != '':
                        ef_slope_points.append((timediff, cur_ef))
                        if abs(timediff) <= (90.0/365.0) and 'UW_EF_BL_3months' not in new_data:
                            new_data['UW_EF_BL_3months'] = cur_ef
        # calc slope
        new_data['UW_MEM_postAV45_count'] = len(mem_slope_points)
        new_data['UW_MEM_slope'] = slope(mem_slope_points)
        new_data['UW_EF_postAV45_count'] = len(ef_slope_points)
        new_data['UW_EF_slope'] = slope(ef_slope_points)
        return new_data

    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, uw_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def syncUCSFFreesurferCrossData(old_headers, old_lines, ucsf_files, mprage_file, dump_to=None):
    tmpts = 12

    fs_by_subj_list = [importUCSFFreesurfer(filename, mprage_file, version=fsversion, include_failed=False) for fsversion,filename in ucsf_files]
    fs_by_subj = {}
    for k,v in itertools.chain(*[_.iteritems() for _ in fs_by_subj_list]):
        if k not in fs_by_subj:
            fs_by_subj[k] = v
        else:
            fs_by_subj[k] = fs_by_subj[k] + v

    to_add_headers = []
    to_add_headers += ['FSX_HC/ICV_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FSX_postAV45_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FSX_postAV45_count', 'FSX_HC/ICV_slope', 'FSX_HC/ICV_BL_3months']
    to_add_headers += ['FSX_MRI_STRENGTH_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FSX_FSVERSION_%s' % (i+1) for i in range(tmpts)]
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        '''
        If RID<2000, use only 1.5T scans
        If RID>2000, use only 3T scans
        '''
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        mristrengths = [_['FLDSTRENG'] for _ in subj_row]
        # filter by field strength
        if int(subj) < 2000:
            subj_row = [_ for _ in subj_row if float(_['FLDSTRENG']) == 1.5]
        else:
            subj_row = [_ for _ in subj_row if float(_['FLDSTRENG']) == 3.0]
        if len(subj_row) == 0:
            print "%s had all potential scans eliminated by field strength: %s" % (subj,mristrengths)
            return {}

        # take care of scans on the same date
        by_date = defaultdict(list)
        for scan in subj_row:
            by_date[scan['EXAMDATE']].append(scan)
        subj_fs = []
        for scandate in sorted(by_date.keys()):
            scans_on_date = by_date[scandate]
            if len(scans_on_date) == 1:
                subj_fs.append(scans_on_date[0])
            else:
                # average if possible, else take first one
                if len(set([_['FLDSTRENG'] for _ in scans_on_date])) == 1 and len(set([_['version'] for _ in scans_on_date])) == 1:
                    avg_row = {}
                    avg_row['EXAMDATE'] = scandate
                    avg_row['ICV'] = np.mean([float(_['ICV']) for _ in scans_on_date])
                    avg_row['HCV'] = np.mean([float(_['HCV']) for _ in scans_on_date])
                    avg_row['FLDSTRENG'] = scans_on_date[0]['FLDSTRENG']
                    avg_row['version'] = scans_on_date[0]['version']
                else:
                    avg_row = scans_on_date[0]
                subj_fs.append(avg_row)

        if len(subj_fs) > tmpts:
            raise Exception("INCREASE NUMBER OF FSX TIMEPOINTS TO > %s" % (len(subj_fs)))

        new_data = {}
        slope_points = []
        bl_icv = subj_fs[0]['ICV']
        for i in range(tmpts):
            if i < len(subj_fs):
                datapoint = subj_fs[i]
                examdate = datapoint['EXAMDATE']
                hc_icv = datapoint['HCV']/bl_icv
                try:
                    timediff = ((examdate-bl_av45).days / 365.0) if bl_av45 else ''
                    if timediff <= -(90.0/365.0):
                        timediff = ''
                    if timediff != '':
                        slope_points.append((timediff, hc_icv))
                except:
                    print "BAD EXAMDATE"
                    continue
                if timediff != '' and abs(timediff) <= (90.0/365.0) and 'FSX_HC/ICV_BL_3months' not in new_data:
                    new_data['FSX_HC/ICV_BL_3months'] = hc_icv
                new_data['FSX_HC/ICV_%s' % (i+1)] = hc_icv
                new_data['FSX_postAV45_%s' % (i+1)] = timediff
                new_data['FSX_FSVERSION_%s' % (i+1)] = datapoint['version']
                new_data['FSX_MRI_STRENGTH_%s' % (i+1)] = datapoint['FLDSTRENG']
        
        # calc slope
        new_data['FSX_postAV45_count'] = len(slope_points)
        new_data['FSX_HC/ICV_slope'] = slope(slope_points)
        
        return new_data

    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, fs_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def syncUCSFFreesurferLongData(old_headers, old_lines, ucsf_files, mprage_file, dump_to=None):
    tmpts = 12

    fs_by_subj_list = [importUCSFFreesurfer(filename, mprage_file, version=fsversion, include_failed=False) for fsversion,filename in ucsf_files]
    fs_by_subj = {}
    for k,v in itertools.chain(*[_.iteritems() for _ in fs_by_subj_list]):
        if k not in fs_by_subj:
            fs_by_subj[k] = v
        else:
            fs_by_subj[k] = fs_by_subj[k] + v

    to_add_headers = []
    to_add_headers += ['FSL_HC/ICV_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FSL_postAV45_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FSL_postAV45_count', 'FSL_HC/ICV_slope', 'FSL_HC/ICV_BL_3months']
    to_add_headers += ['FSL_MRI_STRENGTH_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FSL_FSVERSION_%s' % (i+1) for i in range(tmpts)]
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        mristrengths = [_['FLDSTRENG'] for _ in subj_row]
        # filter by field strength
        if int(subj) < 2000:
            subj_row = [_ for _ in subj_row if float(_['FLDSTRENG']) == 1.5]
        else:
            subj_row = [_ for _ in subj_row if float(_['FLDSTRENG']) == 3.0]
        if len(subj_row) == 0:
            print "%s had all potential scans eliminated by field strength: %s" % (subj,mristrengths)
            return {}

        # take care of scans on the same date
        by_date = defaultdict(list)
        for scan in subj_row:
            by_date[scan['EXAMDATE']].append(scan)
        subj_fs = []
        for scandate in sorted(by_date.keys()):
            scans_on_date = by_date[scandate]
            if len(scans_on_date) == 1:
                subj_fs.append(scans_on_date[0])
            else:
                if len(set([_['FLDSTRENG'] for _ in scans_on_date])) == 1 and len(set([_['version'] for _ in scans_on_date])) == 1:
                    # average
                    avg_row = {}
                    avg_row['EXAMDATE'] = scandate
                    avg_row['ICV'] = np.mean([float(_['ICV']) for _ in scans_on_date])
                    avg_row['HCV'] = np.mean([float(_['HCV']) for _ in scans_on_date])
                    avg_row['FLDSTRENG'] = scans_on_date[0]['FLDSTRENG']
                    avg_row['version'] = scans_on_date[0]['version']
                else:
                    avg_row = scans_on_date[0]
                subj_fs.append(avg_row)

        if len(subj_fs) > tmpts:
            raise Exception("INCREASE NUMBER OF FSL TIMEPOINTS TO > %s" % (len(subj_fs)))

        new_data = {}
        slope_points = []
        bl_icv = float(subj_fs[0]['ICV'])
        for i in range(tmpts):
            if i < len(subj_fs):
                datapoint = subj_fs[i]
                examdate = datapoint['EXAMDATE']
                hc_icv = datapoint['HCV']/bl_icv
                try:
                    timediff = ((examdate-bl_av45).days / 365.0) if bl_av45 else ''
                    if timediff <= -(90.0/365.0):
                        timediff = ''
                    if timediff != '':
                        slope_points.append((timediff, hc_icv))
                except:
                    print "BAD EXAMDATE"
                    continue
                if timediff != '' and abs(timediff) <= (90.0/365.0) and 'FSL_HC/ICV_BL_3months' not in new_data:
                    new_data['FSL_HC/ICV_BL_3months'] = hc_icv
                new_data['FSL_HC/ICV_%s' % (i+1)] = datapoint['HCV']/bl_icv
                new_data['FSL_postAV45_%s' % (i+1)] = timediff
                new_data['FSL_FSVERSION_%s' % (i+1)] = datapoint['version']
                new_data['FSL_MRI_STRENGTH_%s' % (i+1)] = datapoint['FLDSTRENG']
        # Calculate slope
        new_data['FSL_postAV45_count'] = len(slope_points)
        new_data['FSL_HC/ICV_slope'] = slope(slope_points)

        return new_data

    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, fs_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def syncWMHData(old_headers, old_lines, wmh_file, registry_file, dump_to=None):
    tmpts = 7
    wmh_by_subj = importWMH(wmh_file)
    registry = importRegistry(registry_file)

    to_add_headers = ['WMH_percentOfICV.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['WMH_postAV45.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['WMH_WHITMATHYP.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['WMH_percentOfICV_AV45_6MTHS', 'WMH_percentOfICV_AV45_DATE',
                       'WMH_percentOfICV_AV45_2_6MTHS', 'WMH_percentOfICV_AV45_2_DATE',
                       'WMH_percentOfICV_AV45_3_6MTHS', 'WMH_percentOfICV_AV45_3_DATE',
                       'WMH_percentOfICV_slope']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        subj_wmh = sorted(subj_row, key=lambda x: x['EXAMDATE'])
        new_data = {}

        # get closest to av45s
        six_months = 31 * 6
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToAV45(subj_wmh, bl_av45, av45_2, av45_3, day_limit=six_months)
        new_data['WMH_percentOfICV_AV45_6MTHS'] = av45_bl_closest['wmh_percent'] if av45_bl_closest else ''
        new_data['WMH_percentOfICV_AV45_DATE'] = av45_bl_closest['EXAMDATE'] if av45_bl_closest else ''
        new_data['WMH_percentOfICV_AV45_2_6MTHS'] = av45_2_closest['wmh_percent'] if av45_2_closest else ''
        new_data['WMH_percentOfICV_AV45_2_DATE'] = av45_2_closest['EXAMDATE'] if av45_2_closest else ''
        new_data['WMH_percentOfICV_AV45_3_6MTHS'] = av45_3_closest['wmh_percent'] if av45_3_closest else ''
        new_data['WMH_percentOfICV_AV45_3_DATE'] = av45_3_closest['EXAMDATE'] if av45_3_closest else ''
        
        slope_points = []
        for i in range(tmpts):
            if i < len(subj_wmh):
                datapoint = subj_wmh[i]
                examdate = datapoint['EXAMDATE']
                wmh_raw = datapoint['wmh']
                wmh_val = datapoint['wmh_percent']
                timediff = ((examdate-bl_av45).days / 365.0) if bl_av45 else ''
                if timediff <= -(90.0/365.0):
                    timediff = ''
                new_data['WMH_percentOfICV.%s' % (i+1)] = wmh_val
                new_data['WMH_postAV45.%s' % (i+1)] = timediff
                new_data['WMH_WHITMATHYP.%s' % (i+1)] = wmh_raw
                if timediff and wmh_val:
                    slope_points.append((timediff, wmh_val))
            else:
                new_data['WMH_percentOfICV.%s' % (i+1)] = ''
                new_data['WMH_postAV45.%s' % (i+1)] = ''
                new_data['WMH_WHITMATHYP.%s' % (i+1)] = ''
        # calculate slope
        new_data['WMH_percentOfICV_slope'] = slope(slope_points)
        return new_data

    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, wmh_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def getClosestToAV45(points, bl_av45, av45_2, av45_3, day_limit=550):
    used = set()
    # match up with av45 scans
    av45_bl_closest = av45_2_closest = av45_3_closest = None
    if bl_av45 is not None:
        eligible = [p for p in points if p['EXAMDATE'] not in used and (abs(bl_av45 - p['EXAMDATE']).days <= day_limit)]
        cand = sorted(eligible, key=lambda x: abs(bl_av45-x['EXAMDATE']))
        if len(cand) > 0:
            av45_bl_closest = cand[0]
            used.add(av45_bl_closest['EXAMDATE'])
    if av45_2 is not None:
        eligible = [p for p in points if p['EXAMDATE'] not in used and (abs(av45_2 - p['EXAMDATE']).days <= day_limit)]
        cand = sorted(eligible, key=lambda x: abs(av45_2-x['EXAMDATE']))
        if len(cand) > 0:
            av45_2_closest = cand[0]
            used.add(av45_2_closest['EXAMDATE'])
    if av45_3 is not None:
        eligible = [p for p in points if p['EXAMDATE'] not in used and (abs(av45_3 - p['EXAMDATE']).days <= day_limit)]
        cand = sorted(eligible, key=lambda x: abs(av45_3-x['EXAMDATE']))
        if len(cand) > 0:
            av45_3_closest = cand[0]
            used.add(av45_3_closest['EXAMDATE'])
    return av45_bl_closest, av45_2_closest, av45_3_closest


def eliminateColumns(headers, lines):
    to_remove = ['LastClinicalVisit',
                 'DIFF_LastClinicalVisit_AV45',
                 'ClinicalVisitClosestto_AV45',
                 'ClinicalVisitClosestto_AV45_2',
                 'MMSEclosest_1', 
                 'MMSEclosest_AV45',
                 'datediff_MMSE_AV45', 
                 'MMSE_3mths_AV45',
                 'PostAV45Followup',
                 'abeta_closest_AV45',
                 'abeta_bin192',
                 'abeta_closest_AV45_2',
                 'abeta_2_bin192',
                 'abeta_diff',
                 'abeta_slope',
                 'tau_closest_AV45',
                 'tau_bin93',
                 'tau_closest_AV45_2',
                 'ptau_closest_AV45',
                 'ptau_bin23',
                 'ptau_closest_2',
                 'WMH_CLOSEST_AV45',
                 'WMH_CLOSEST_',
                 'WMH_slope',
                 'ClosestCSF_AV45',
                 'diff_CSF_AV45',
                 'AV45CSF<3months',
                 'ClosestCSF_AV45_2',
                 'diff_CSF_AV45_2',
                 'CSF<3mths_AV45_2',
                 'LastCSF',
                 'LastCSF_AV45_1_TimeDiff',
                 'LastCSFAbeta',
                 'ADAS_slope_all',
                 'FDG_slope',
                 'TBMSyn_SLOPE',
                 'TBMSyn_slope_followuptime',
                 'Hippocampal_Volume_3months',
                 'Left_Volume',
                 'Right_Volume',
                 'ICV',
                 'Hippocampal_Volume_normalizd',
                 'HippVol_norm_BIN0.0047',
                 'ICV_pchange(yrs)',]
    to_remove += ['WHITMATHYP.%s' % (i+1) for i in range(5)]
    to_remove += ['ABETA.%s' % (i+1) for i in range(7)]
    to_remove += ['TAU.%s' % (i+1) for i in range(7)]
    to_remove += ['PTAU.%s' % (i+1) for i in range(7)]
    to_remove += ['EXAM_tau.%s' % (i+1) for i in range(7)]
    to_remove += ['ADAScog_DATE%s' % (i+1) for i in range(11)]
    to_remove += ['MMSE_DATE%s' % (i+1) for i in range(11)]
    to_remove += ['AVLT_DATE.%s' % (i+1) for i in range(11)]
    to_remove += ['TBMSyn_DATE.%s' % (i+1) for i in range(11)]

    for tm in to_remove:
        while tm in headers:
            headers.remove(tm)
    new_lines = []
    for l in lines:
        for tm in to_remove:
            l.pop(tm, None)
        new_lines.append(l)
    return headers, new_lines

def getAV45Dates(old_l, patient_pets=None):
    if patient_pets is None:
        patient_pets = []
    # Get AV45 Scan dates
    if old_l.get('AV45_Date','') != '':
        bl_av45 = datetime.strptime(old_l['AV45_Date'], '%m/%d/%y')
    elif len(patient_pets) >= 1:
        bl_av45 = patient_pets[0]
    else:
        bl_av45 = None
    if old_l.get('AV45_2_Date','') != '':
        av45_2 = datetime.strptime(old_l['AV45_2_Date'], '%m/%d/%y')
    elif len(patient_pets) >= 2:
        av45_2 = patient_pets[1]
    else:
        av45_2 = None
    if old_l.get('AV45_3_Date','') != '':
        av45_3 = datetime.strptime(old_l['AV45_3_Date'], '%m/%d/%y')
    elif len(patient_pets) >= 3:
        av45_3 = patient_pets[2]
    else:
        av45_3 = None
    return (bl_av45, av45_2, av45_3)

def wipeKeys(lines, keys):
    empty_dict = {k:'' for k in keys}
    new_lines = []
    for l in lines:
        l.update(empty_dict)
        new_lines.append(l)
    return l

def addCategories(output_file):
    colcats_df = pd.read_csv('column_categories.csv')
    colcats_df.set_index('COLUMN',inplace=True)
    df = pd.read_csv(output_file, low_memory=False)
    cats = []
    for col in df.columns:
        try:
            cat = colcats_df.loc[col.strip(),'CATEGORY']
        except Exception as e:
            cat = 'Unknown'
            print col
        cats.append(cat)
    df.columns = [cats, df.columns]
    df.to_csv(output_file, index=False)

def mergeSlopes(output_file):
    # merge 2/3 point slope columns into one column (replacing 2pt slopes where 3pt slopes are available)
    # add merged column in after 3 point column
    df = pd.read_csv(output_file, low_memory=False)
    twopoint_cols = [i for i in df.columns if '2pts' in i or '2points' in i]
    threepoint_cols = [i for i in df.columns if '3pts' in i or '3points' in i]
    for two_col in twopoint_cols:
        # find three point col
        new_col = two_col.replace('_2pts','').replace('_2points','')
        if new_col in df.columns:
            continue
        three_col = two_col.replace('2pts','3pts').replace('2points', '3points')
        if three_col not in threepoint_cols:
            print "Could not find 3 point slope column for %s" % (three_col,)
            continue
        three_col_ind = df.columns.get_loc(three_col)
        df.insert(three_col_ind+1, new_col, df.loc[:,two_col])
        df[new_col].update(df[three_col])
    df.to_csv(output_file, index=False)

def runPipeline():
    # syncing pipeline
    new_headers, new_lines = parseCSV(master_file, use_second_line=True)
    print "\nELIMINATING COLUMNS\n"
    new_headers, new_lines = eliminateColumns(new_headers, new_lines)
    print "\nSYNCING AV45\n"
    new_headers, new_lines = syncAV45Data(new_headers, new_lines, av45_file, av45_nontp_file, registry_file, dump_to=None) # adds new patients
    print "\nSYNCING DIAGNOSES\n"
    new_headers, new_lines = syncDiagnosisData(new_headers, new_lines, diagnosis_file, registry_file, demog_file, arm_file, pet_meta_file, dump_to=None) # refreshes av45 dates
    print "\nSYNCING FAQ\n"
    new_headers, new_lines = syncFAQData(new_headers, new_lines, faq_file, registry_file, dump_to=None)
    print "\nSYNCING NPI\n"
    new_headers, new_lines = syncNPIData(new_headers, new_lines, npi_file, dump_to=None)
    print "\nSYNCING MHIST\n"
    new_headers, new_lines = syncMHISTData(new_headers, new_lines, mhist_file, dump_to=None)
    print "\nSYNCING UW NEURO\n"
    new_headers, new_lines = syncUWData(new_headers, new_lines, uw_file, registry_file, dump_to=None) # refreshes av45 dates
    print "\nSYNCING UCSF LONG FreeSurfer\n"
    new_headers, new_lines = syncUCSFFreesurferLongData(new_headers, new_lines, ucsf_long_files, mprage_file, dump_to=None)
    print "\nSYCNING UCSF CROSS FreeSurfer"
    new_headers, new_lines = syncUCSFFreesurferCrossData(new_headers, new_lines, ucsf_cross_files, mprage_file, dump_to=None)
    print "\nSYNCING ROUSSET BL\n"
    new_headers, new_lines = syncRoussetResults(new_headers, new_lines, rousset_matfile_bl, 'BL', dump_to=None)
    print "\nSYNCING ROUSSET SCAN2\n"
    new_headers, new_lines = syncRoussetResults(new_headers, new_lines, rousset_matfile_scan2, 'Scan2', dump_to=None)
    print "\nSYNCING ROUSSET SCAN3\n"
    new_headers, new_lines = syncRoussetResults(new_headers, new_lines, rousset_matfile_scan3, 'Scan3', dump_to=None)
    print "\nSYNCING FDG\n"
    new_headers, new_lines = syncFDGData(new_headers, new_lines, fdg_file, registry_file, dump_to=None)
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
    print "\nSYNCING GD\n"
    new_headers, new_lines = syncGDData(new_headers, new_lines, gd_file, registry_file, dump_to=None)
    print "\nSYNCING WMH\n"
    new_headers, new_lines = syncWMHData(new_headers, new_lines, wmh_file, registry_file, dump_to=output_file)

    # postprocessing
    mergeSlopes(output_file)
    addCategories(output_file)

if __name__ == '__main__':
    now = datetime.now()

    # IO files
    master_file = "../FDG_AV45_COGdata_10_27_15.csv"
    output_file = "../FDG_AV45_COGdata_%s.csv" % (now.strftime("%m_%d_%y"))
    
    # ADNI docs
    registry_file = "../docs/ADNI/REGISTRY.csv"
    arm_file = "../docs/ADNI/ARM.csv"
    pet_meta_file = "../docs/ADNI/PET_META_LIST.csv"
    demog_file = "../docs/ADNI/PTDEMOG.csv"
    diagnosis_file = '../docs/ADNI/DXSUM_PDXCONV_ADNIALL.csv'
    csf_files = ['../docs/ADNI/UPENNBIOMK5_10_31_13.csv',
                 '../docs/ADNI/UPENNBIOMK6_07_02_13.csv',
                 '../docs/ADNI/UPENNBIOMK7.csv',
                 '../docs/ADNI/UPENNBIOMK8.csv']
    wmh_file = '../docs/ADNI/UCD_ADNI2_WMH_10_26_15.csv'
    fdg_file = '../docs/ADNI/UCBERKELEYFDG_07_30_15.csv'
    mprage_file = '../docs/ADNI/MPRAGEMETA.csv'
    mhist_file = '../docs/ADNI/RECMHIST.csv'
    faq_file = '../docs/ADNI/FAQ.csv'
    npi_file = '../docs/ADNI/NPI.csv'
    
    # Output files
    av45_file = None # overwrites nontp if specified
    av45_nontp_file = "../output/01_03_16/UCBERKELEYAV45_01_03_16_merged_nontp.csv"
    rousset_matfile_bl = '../output/Rousset_BL/rousset_output_BL.mat'
    rousset_matfile_scan2 = '../output/Rousset_Scan2/rousset_output_Scan2.mat'
    rousset_matfile_scan3 = '../output/Rousset_Scan3/rousset_output_Scan3.mat'

    # Cog files
    mmse_file = "../cog_tests/MMSE.csv"
    adni1_adas_file = '../cog_tests/ADASSCORES.csv'
    adnigo2_adas_file = '../cog_tests/ADAS_ADNIGO2.csv'
    neuro_battery_file = '../cog_tests/NEUROBAT.csv'
    gd_file = '../cog_tests/GDSCALE.csv'
    uw_file = '../cog_tests/UWNPSYCHSUM_10_26_15.csv'
    
    # MR files
    tbm_file = '../mr_docs/Mayo/MAYOADIRL_MRI_TBMSYN_05_07_15.csv'
    ucsf_long_files = [('4.4','../mr_docs/UCSF/longitudinal/UCSFFSL_11_02_15.csv'),
                       ('5.1','../mr_docs/UCSF/longitudinal/UCSFFSL51Y1_11_02_15.csv')]
    ucsf_cross_files = [('4.3','../mr_docs/UCSF/cross_section/UCSFFSX_11_02_15.csv'),
                        ('5.1','../mr_docs/UCSF/cross_section/UCSFFSX51_11_02_15.csv'),
                        ('5.1','../mr_docs/UCSF/cross_section/UCSFFSX51_ADNI1_3T_11_02_15.csv')]

    # Run pipeline
    runPipeline()
