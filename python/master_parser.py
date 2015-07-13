'''

Pull in the master FDG_AV45_COG DATA

remove null char:
tr -d '\000' < file1 > file2

'''
import re
from collections import defaultdict
from datetime import datetime, timedelta
from scipy import stats, optimize

from utils import *

def syncADASCogData(old_headers, old_lines, adni1_adas_file, adnigo2_adas_file, registry_file, dump_to=None):
    adni1_headers, adni1_lines = parseCSV(adni1_adas_file)
    adni2_headers, adni2_lines = parseCSV(adnigo2_adas_file)
    registry = importRegistry(registry_file)

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
        subj = int(line['RID'])
        viscode = line['VISCODE'].strip().lower()
        viscode2 = line['VISCODE2'].strip().lower()
        raw_totscore = line['TOTSCORE']
        if raw_totscore == '':
            print "%s (%s, %s) has missing totscore" % (subj, viscode, viscode2)
            continue
        totscore = int(raw_totscore)
        if totscore < 0:
            totscore = ''
        examdate = None
        subj_listings = registry[subj]
        for listing in subj_listings:
            if listing['VISCODE'] == viscode and listing['VISCODE2'] == viscode2:
                examdate = listing['date']
                break
        if examdate is None:
            print "No exam date for %s (%s, %s)" % (subj, viscode, viscode2)
            continue
        adas_by_subj[subj].append({'VISCODE': viscode,
                                   'VISCODE2': viscode2,
                                   'EXAMDATE': examdate,
                                   'TOTSCORE': totscore})
    adas_by_subj = dict(adas_by_subj)

    new_headers = old_headers
    new_lines = []
    new_values = 0
    total = 0

    to_add_headers = ['ADAS_AV45_3_3MTHS', 'ADAS_AV45_3_DATE']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after='ADAS_AV45_2_DATE')
    to_add_headers = ['ADAS_post_AV45_followuptime']
    new_headers = rearrangeHeaders(new_headers, to_add_headers, after='ADAS_3MTH_AV45')

    # remove date fields
    for i in range(11):
        key = 'ADAScog_DATE%s' % (i+1)
        if key in new_headers:
            new_headers.remove(key)

    for linenum, old_l in enumerate(old_lines):
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            continue

        unsorted_tests = adas_by_subj.get(subj,[])
        if len(unsorted_tests) == 0:
            print "No ADAS tests found for %s" % (subj)
            new_lines.append(old_l)
            continue

        tests = sorted(unsorted_tests, key=lambda x : x['EXAMDATE'])
        dates = [t['EXAMDATE'] for t in tests]
        if len(set(dates)) != len(dates):
            print "%s has Dup dates, skipping: %s" % (subj, dates)
            new_lines.append(old_l)
            continue
        first_scan_date = dates[0]


        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=None)

        new_subj_data = {}
        all_values = []
        post_values = []
        all_times = []
        post_times = []
        max_followup_counter = None
        for i in range(11):
            if i >= len(tests):
                test_date = ''
                test_score = ''
                diff_from_first = ''
            else:
                test_results = tests[i]
                test_date = test_results['EXAMDATE']
                test_score = test_results['TOTSCORE']
                if test_score != '':
                    test_score = float(test_score)
                diff_from_first = (test_date-first_scan_date).days / 365.0
                if test_score != '':
                    all_values.append(test_score)
                    all_times.append(diff_from_first)
            count = i+1
            new_subj_data['ADAScog_DATE%s' % count] = test_date
            new_subj_data['ADAScog.%s' % count] = test_score
            new_subj_data['TIME_ADAS.%s' % count] = diff_from_first
            if bl_av45 is not None and test_date != '':
                rel_time_days = (test_date - bl_av45).days
                rel_time = round(rel_time_days / 365.0, 2)
                new_subj_data['TIMEreltoAV45_ADAS.%s' % count] = rel_time
                if abs(rel_time_days) <= 93 and 'ADAS_3MTH_AV45' not in new_subj_data:
                    new_subj_data['ADAS_3MTH_AV45'] = test_score
                    new_subj_data['ADAS_3MTHS_AV45DATE'] = test_date
                if rel_time >= (-93.0/365.0):
                    if test_score != '':
                        post_values.append(test_score)
                        post_times.append(diff_from_first)
                    new_subj_data['TIMEpostAV45_ADAS.%s' % count] = rel_time
                    max_followup_counter = rel_time
                else:
                    new_subj_data['TIMEpostAV45_ADAS.%s' % count] = ''
            if av45_2 is not None and test_date != '':
                rel_time_days = (test_date - av45_2).days 
                if abs(rel_time_days) <= 93 and 'ADAS_AV45_2_3MTHS' not in new_subj_data:
                    new_subj_data['ADAS_AV45_2_3MTHS'] = test_score
                    new_subj_data['ADAS_AV45_2_DATE'] = test_date
            if av45_3 is not None and test_date != '':
                rel_time_days = (test_date - av45_3).days 
                if abs(rel_time_days) <= 93 and 'ADAS_AV45_3_3MTHS' not in new_subj_data:
                    new_subj_data['ADAS_AV45_3_3MTHS'] = test_score
                    new_subj_data['ADAS_AV45_3_DATE'] = test_date
        if max_followup_counter is not None:
            new_subj_data['ADAS_post_AV45_followuptime'] = max(max_followup_counter, 0.0)
        # fill in the blanks
        fill_in = ['ADAS_3MTH_AV45', 'ADAS_3MTHS_AV45DATE', 
                   'ADAS_AV45_2_3MTHS','ADAS_AV45_2_DATE',
                   'ADAS_AV45_3_3MTHS','ADAS_AV45_3_DATE',
                   'ADAS_post_AV45_followuptime']
        for f_key in fill_in:
            if f_key not in new_subj_data:
                new_subj_data[f_key] = ''
        
        # get slopes
        if len(all_values) >= 2:
            slope, intercept, r, p, stderr = stats.linregress(all_times, all_values)
            new_subj_data['ADAS_slope_all'] = slope
        else:
            new_subj_data['ADAS_slope_all'] = ''
        if len(post_values) >= 2:
            slope, intercept, r, p, stderr = stats.linregress(post_times, post_values)
            new_subj_data['ADASslope_postAV45'] = slope
        else:
            new_subj_data['ADASslope_postAV45'] = ''

        
        # do comparison:
        #print "SUBJ: %s" % subj
        #print "first: %s, second: %s" % (bl_av45, av45_2)
        changed = False
        for k in sorted(new_subj_data.keys()):
            old_value = old_l.get(k, None) # new columns potentially added
            new_value = new_subj_data[k]
            if k.startswith('ADAScog.') and old_value == '' and new_value != '':
                new_values += 1
                changed = True
            #print "\t%s: %s -> %s" % (k, old_value, new_value)
        if changed:
            total += 1
        
        new_subj_data = convertToCSVDataType(new_subj_data)
        old_l.update(new_subj_data)
        new_lines.append(old_l)

    print "Total subj changed: %s" % total
    print "Total new tests: %s" % new_values

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def syncMMSEData(old_headers, old_lines, mmse_file, registry_file, dump_to=None):
    registry = importRegistry(registry_file)
    mmse_by_subject = importMMSE(mmse_file, registry=registry)

    # add these new headers, if not present
    to_add_headers = ['MMSE_post_AV45_followuptime',
                      'MMSE_AV45_3MTHS',
                      'MMSE_AV45_DATE',
                      'MMSE_AV45_2_3MTHS',
                      'MMSE_AV45_2_DATE',
                      'MMSE_AV45_3_3MTHS',
                      'MMSE_AV45_3_DATE',
                      'MMSEslope_postAV45']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after='MMSCORE.11')

    # remove date fields
    for i in range(11):
        key = 'MMSE_DATE%s' % (i+1)
        if key in new_headers:
            new_headers.remove(key)

    new_lines = []
    new_values = 0
    total = 0
    for i, old_l in enumerate(old_lines):
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            continue

        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=None)

        # create new subject data
        tests = mmse_by_subject[subj]
        new_subj_data = {}
        max_followup_counter = None
        post_av45_points = []
        for i in range(11):
            if i >= len(tests):
                test_date = ''
                test_score = ''
            else:
                test_date, test_results = tests[i]
                if test_results['MMSCORE'] == '':
                    test_score = ''
                else:
                    test_score = int(test_results['MMSCORE'])
                    if test_score < 0:
                        test_score = ''

            count = i+1
            new_subj_data['MMSE_DATE%s' % count] = test_date
            new_subj_data['MMSCORE.%s' % count] = test_score

            # pair up with av45 scan date if possible
            if bl_av45 is not None and test_date != '':
                rel_time_days = (test_date - bl_av45).days 
                if abs(rel_time_days) <= 93 and 'MMSE_AV45_3MTHS' not in new_subj_data:
                    new_subj_data['MMSE_AV45_3MTHS'] = test_score
                    new_subj_data['MMSE_AV45_DATE'] = test_date
                # also try to put value into time post av45 columns
                if rel_time_days >= -93.0:
                    annualized_time = rel_time_days / 365.0
                    new_subj_data['TIMEpostAV45_MMSE.%s' % count] = annualized_time
                    max_followup_counter = annualized_time
                    if test_score != '':
                        post_av45_points.append((annualized_time, test_score))
            if ('TIMEpostAV45_MMSE.%s' % count) not in new_subj_data:
                new_subj_data['TIMEpostAV45_MMSE.%s' % count] = ''
            # pair up with subsequent av45 scans
            if av45_2 is not None and test_date != '':
                rel_time_days = (test_date - av45_2).days 
                if abs(rel_time_days) <= 93 and 'MMSE_AV45_2_3MTHS' not in new_subj_data:
                    new_subj_data['MMSE_AV45_2_3MTHS'] = test_score
                    new_subj_data['MMSE_AV45_2_DATE'] = test_date
            if av45_3 is not None and test_date != '':
                rel_time_days = (test_date - av45_3).days 
                if abs(rel_time_days) <= 93 and 'MMSE_AV45_3_3MTHS' not in new_subj_data:
                    new_subj_data['MMSE_AV45_3_3MTHS'] = test_score
                    new_subj_data['MMSE_AV45_3_DATE'] = test_date
        if max_followup_counter is not None:
            new_subj_data['MMSE_post_AV45_followuptime'] = max(max_followup_counter, 0.0)

        # get post av45 slope
        if len(post_av45_points) >= 2:
            post_times = [_[0] for _ in post_av45_points]
            post_values = [_[1] for _ in post_av45_points]
            slope, intercept, r, p, stderr = stats.linregress(post_times, post_values)
            new_subj_data['MMSEslope_postAV45'] = slope

        # fill in the blanks
        for f_key in to_add_headers:
            if f_key not in new_subj_data:
                new_subj_data[f_key] = ''

        # do comparison:
        #print "SUBJ: %s" % subj
        changed = False
        for k in sorted(new_subj_data.keys()):
            old_value = old_l.get(k,'')
            new_value = new_subj_data[k]
            if k.startswith('MMSCORE') and old_value == '' and new_value != '':
                new_values += 1
                changed = True
            #print "\t%s: %s -> %s" % (k, old_value, new_value)
        if changed:
            total += 1

        new_data = convertToCSVDataType(new_subj_data)
        old_l.update(new_data)
        new_lines.append(old_l)

    print "%s subj with new tests (%s new tests)" % (total, new_values)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return new_headers, new_lines

def syncAVLTData(old_headers, old_lines, neuro_battery_file, registry_file, dump_to=None):
    registry = importRegistry(registry_file)
    avlt_by_subj = importAVLT(neuro_battery_file, registry=registry)

    # no change in headers
    to_add_headers = ['AVLT_post_AV45_followuptime']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after='AVLT_3MTHS_AV45')

    # remove date fields
    for i in range(11):
        key = 'AVLT_DATE.%s' % (i+1)
        if key in new_headers:
            new_headers.remove(key)

    new_lines = []
    new_values = 0
    total = 0
    subj_new_counts = {}
    for linenum, old_l in enumerate(old_lines):
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            continue

        unsorted_tests = avlt_by_subj.get(subj,[])
        if len(unsorted_tests) == 0:
            print "No AVLT tests found for %s" % (subj)
            new_lines.append(old_l)
            continue

        tests = sorted(unsorted_tests, key=lambda x : x['EXAMDATE'])
        dates = [t['EXAMDATE'] for t in tests]
        if len(set(dates)) != len(dates):
            print "%s has Dup dates, skipping: %s" % (subj, dates)
            new_lines.append(old_l)
            continue
        first_scan_date = dates[0]

        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=None)

        new_subj_data = {}
        all_values = []
        post_values = []
        all_times = []
        post_times = []
        max_followup_counter = None
        for i in range(11):
            if i >= len(tests):
                test_date = ''
                test_date_string = ''
                test_score = ''
                diff_from_first = ''
            else:
                test_results = tests[i]
                test_date = test_results['EXAMDATE']
                test_date_string = test_date.strftime("%m/%d/%y")
                test_score = test_results['TOTS']
                diff_from_first = (test_date-first_scan_date).days / 365.0
                if test_score != '':
                    all_values.append(test_score) 
                    all_times.append(diff_from_first)
            count = i+1
            new_subj_data['AVLT_DATE.%s' % count] = test_date_string
            new_subj_data['AVLT.%s' % count] = test_score
            new_subj_data['TIME_AVLT.%s' % count] = diff_from_first
            if bl_av45 is not None and test_date != '':
                rel_time_days = (test_date - bl_av45).days
                rel_time = rel_time_days / 365.0
                new_subj_data['TIMEreltoAV45_AVLT.%s' % count] = rel_time
                if abs(rel_time_days) <= 93 and 'AVLT_3MTHS_AV45' not in new_subj_data:
                    new_subj_data['AVLT_3MTHS_AV45'] = test_score
                    new_subj_data['AVLT_3MTHSAV45_Date'] = test_date_string
                if rel_time >= (-93.0/365.0):
                    if test_score != '':
                        post_values.append(test_score)
                        post_times.append(diff_from_first)
                    new_subj_data['TIMEpostAV45_AVLT.%s' % count] = rel_time
                    max_followup_counter = rel_time
                else:
                    new_subj_data['TIMEpostAV45_AVLT.%s' % count] = ''
            if av45_2 is not None and test_date != '':
                rel_time_days = (test_date - av45_2).days 
                if abs(rel_time_days) <= 93 and 'AVLT_AV45_2_3MTHS' not in new_subj_data:
                    new_subj_data['AVLT_AV45_2_3MTHS'] = test_score
                    new_subj_data['AVLT_AV45_2_DATE'] = test_date_string
        if max_followup_counter is not None:
            new_subj_data['AVLT_post_AV45_followuptime'] = max(max_followup_counter, 0.0)

        # fill in the blanks
        fill_in = ['AVLT_3MTHS_AV45', 'AVLT_3MTHSAV45_Date', 
                   'AVLT_AV45_2_3MTHS','AVLT_AV45_2_DATE',
                   'AVLT_post_AV45_followuptime']
        for f_key in fill_in:
            if f_key not in new_subj_data:
                new_subj_data[f_key] = ''

        # get slopes
        if len(all_values) >= 2:
            try:
                slope, intercept, r, p, stderr = stats.linregress(all_times, all_values)
            except Exception as e:
                print "%s vs %s" % (all_times, all_values)
                raise e
            new_subj_data['AVLT_slope_all'] = slope
        else:
            new_subj_data['AVLT_slope_all'] = ''
        if len(post_values) >= 2:
            slope, intercept, r, p, stderr = stats.linregress(post_times, post_values)
            new_subj_data['AVLTslope_postAV45'] = slope
        else:
            new_subj_data['AVLTslope_postAV45'] = ''

        # do comparison:
        #print "SUBJ: %s" % subj
        #print "first: %s, second: %s" % (bl_av45, av45_2)
        #print "first date: %s" % (first_scan_date)
        changed = False
        new_for_subj = 0
        for k in sorted(new_subj_data.keys()):
            old_value = old_l.get(k,None)
            new_value = new_subj_data[k]
            if k.startswith('AVLT.') and not bool(old_value) and bool(new_value):
                new_values += 1
                new_for_subj += 1
                changed = True
            #print "\t%s: %s -> %s" % (k, old_value, new_value)
        if changed:
            total += 1
        subj_new_counts[subj] = new_for_subj
        
        new_subj_data = convertToCSVDataType(new_subj_data)
        old_l.update(new_subj_data)
        new_lines.append(old_l)

    print "Total subj changed: %s" % total
    print "Total new tests: %s" % new_values
    #print "New test counts per subj: %s" % sorted(subj_new_counts.items(), key=lambda x: x[1], reverse=True)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def syncDiagnosisData(old_headers, old_lines, diag_file, registry_file, arm_file, pet_meta_file, dump_to=None):
    diag_headers, diag_lines = parseCSV(diag_file)
    registry = importRegistry(registry_file)
    arm = importARM(arm_file)
    pet_meta = importPetMETA(pet_meta_file)

    # restructure by subject
    diag_by_subj = defaultdict(list)
    for line in diag_lines:
        subj = int(line['RID'])
        viscode = line['VISCODE'].strip().lower()
        viscode2 = line['VISCODE2'].strip().lower()
        #if viscode == 'sc' or viscode2 == 'sc':
        #    continue
        examdate = line['EXAMDATE']
        change = line['DXCHANGE'].strip()
        current = line['DXCURREN'].strip()
        conv = line['DXCONV'].strip()
        conv_type = line['DXCONTYP'].replace('-4','').strip()
        rev_type = line['DXREV'].replace('-4','').strip()
        if examdate:
            examdate = datetime.strptime(examdate,'%Y-%m-%d')
        else:
            subj_listings = registry[subj]
            for listing in subj_listings:
                if listing['VISCODE'] == viscode and listing['VISCODE2'] == viscode2:
                    examdate = listing['date']
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

    new_headers = old_headers 
    new_lines = []
    new_values = 0
    new_conversions = 0
    total = 0

    pivot_date = datetime(day=14, month=10, year=2014)
    pivot_date_closest_diag_key = 'Closest_DX_Jun15'
    pivot_date_closest_date_key = 'DX_Jun15_closestdate'
    pivot_diag_regex = re.compile('^Closest_DX_')
    pivot_date_regex = re.compile('^DX_.*_closestdate$')

    # replace old pivot date keys
    for i in range(len(new_headers)):
        if re.search(pivot_diag_regex, new_headers[i]):
            new_headers[i] = pivot_date_closest_diag_key
        elif re.search(pivot_date_regex, new_headers[i]):
            new_headers[i] = pivot_date_closest_date_key

    # add 'Diag@AV45_3_long' header if not there
    if 'Diag@AV45_3_long' not in new_headers:
        new_headers.insert(new_headers.index('Diag@AV45_2_long')+1, 'Diag@AV45_3_long')

    for linenum, old_l in enumerate(old_lines):
        # get subject ID
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            continue

        # get subject diagnoses
        subj_diags = diag_by_subj.get(subj,None)
        if subj_diags is None:
            print "No diag data found for %s" % subj
            new_lines.append(old_l)
            continue
        else:
            subj_diags = sorted(subj_diags, key=lambda x: x['EXAMDATE'])

        # check arm for init diag
        init_diag = None
        if subj in arm:
            sorted_arm = sorted(arm[subj], key=lambda x: x['USERDATE'])
            init_diags = list(set([_['STATUS'] for _ in sorted_arm]))
            init_diag = init_diags[-1]

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
        patient_pets = sorted(pet_meta.get(subj,[]))
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        
        # Date differences between scans
        av45_1_2_diff = ((av45_2 - bl_av45).days/365.0) if (bl_av45 is not None and av45_2 is not None) else ''
        av45_1_3_diff = ((av45_3 - bl_av45).days/365.0) if (bl_av45 is not None and av45_3 is not None) else ''


        # find very first visit date in registry (excluding scmri)
        subj_registry = [_ for _ in registry[subj] if _['VISCODE'] != 'scmri' and _['VISCODE2'] != 'scmri']
        sorted_registry = sorted(subj_registry, key=lambda x: x['date'])
        baseline_registry = sorted_registry[0]

        # find pivot data point
        sorted_by_pivot = sorted(subj_diags, key=lambda x: abs(x['EXAMDATE'] - pivot_date))
        closest_to_pivot = sorted_by_pivot[0]

        '''
        Baseline: date
        Init_Diagnosis: N/EMCI/LMCI/AD
        Diag@AV45_long: N/EMCI/LMCI/AD
        Diag@AV45_2_long: N/EMCI/LMCI/AD
        Closest_DX_xxx: N/EMCI/LMCI/AD
        DX_xxx_closestdate: date
        MCItoADConv(fromav45): 0/1
        MCItoADConvDate: date
        MCItoADconv_: date
        AV45_MCItoAD_ConvTime: yrs
        Baseline_MCItoAD_ConvTime: yrs
        FollowupTimetoDX: yrs
        '''

        new_data = {'AV45_Date': bl_av45,
                    'AV45_2_Date': av45_2,
                    'AV45_3_Date': av45_3,
                    'AV45_1_2_Diff (Yrs)': av45_1_2_diff,
                    'AV45_1_3_Diff (yrs)': av45_1_3_diff,
                    'MCItoADConv(fromav45)': '0',
                    'MCItoADConvDate': '',
                    'MCItoADconv_': '',
                    'AV45_MCItoAD_ConvTime': '',
                    'Baseline_MCItoAD_ConvTime': '',
                    'Diag@AV45_long': '',
                    'Diag@AV45_2_long': '',
                    'Diag@AV45_3_long': '',
                    'FollowupTimetoDX': (pivot_date - baseline_registry['date']).days / 365.0,
                    'Baseline': baseline_registry['date'],
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
                
                    
        # convert datetimes to strings
        new_data = convertToCSVDataType(new_data)
        
        # compare
        '''
        print "SUBJ: %s" % subj
        print "ARM: %s" % arm.get(subj,[])
        print "First AV: %s" % bl_av45
        print "Second AV: %s" % av45_2
        print "Diags: %s" % ([(str(_['EXAMDATE']),_['change'],_['diag']) for _ in subj_diags])
        print 'Reg Dates: %s' % sorted_registry
        '''
        changed = False
        ignore_change = set([pivot_date_closest_diag_key, 
                             pivot_date_closest_date_key, 
                             'FollowupTimetoDX', 
                             'MCItoADConvDate', 
                             'MCItoADconv_',
                             'AV45_MCItoAD_ConvTime',
                             'MCItoADConv(fromav45)',
                             'Baseline_MCItoAD_ConvTime',
                             'AV45_Date',
                             'AV45_2_Date'])
        
        for k in sorted(new_data.keys()):
            old_value = old_l.get(k)
            new_value = new_data.get(k)
            if k not in ignore_change and old_value != new_value:
                changed = True
            if k == 'MCItoADConv(fromav45)' and old_value == '0' and new_value == '1':
                new_conversions += 1

        if changed:
            for k in sorted(new_data.keys()):
                if k not in ignore_change:
                    old_value = old_l.get(k)
                    new_value = new_data.get(k)
                    #print "\t%s: %s -> %s" % (k, old_value, new_value)
            new_values += 1


        old_l.update(new_data)
        new_lines.append(old_l)

    print "Changed: %s" % new_values
    print "Total: %s" % len(new_lines)
    print "New conversions: %s" % new_conversions


    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def syncAV45Data(old_headers, old_lines, av45_file, registry_file, dump_to=None):
    '''
    This function adds new subject lines to the output
    '''
    registry = importRegistry(registry_file)
    av45_by_subj = importAV45(av45_file, registry=registry):

    new_headers = None
    new_lines = []
    old_subjects = set([])
    for old_l in old_lines:
        # get subject ID
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            continue
        old_subjects.add(subj)
        if subj not in av45_by_subj:
            print "No AV45 data for %s" % subj
            new_lines.append(old_l)
            continue

        updated_headers, new_data = parseAV45Entries(old_headers, av45_by_subj[subj])
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
    all_wm70_composite_keys = wm70_composite_keys + \
                         ["%s_pchange" % _ for _ in wm70_composite_keys[1:]] + \
                         ["%s_pchange_ABS" % _ for _ in wm70_composite_keys[1:]] + \
                         ["%s_diff" % _ for _ in wm70_composite_keys[1:]] + \
                         ["%s_diff_ABS" % _ for _ in wm70_composite_keys[1:]]
    all_wm70_cerebg_keys = wm70_cerebg_keys + \
                           ["%s_pchange" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ["%s_pchange_ABS" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ["%s_diff" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ["%s_diff_ABS" % _ for _ in wm70_cerebg_keys[1:]]
    all_wm70_wcereb_keys = wm70_wcereb_keys + \
                           ["%s_pchange" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ["%s_pchange_ABS" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ["%s_diff" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ["%s_diff_ABS" % _ for _ in wm70_wcereb_keys[1:]]
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
                      ["%s_diff_ABS" % _ for _ in wmratio_keys[1:]]
    all_frontal_bigref_keys = frontal_bigref_keys + \
                      ["%s_pchange" % _ for _ in frontal_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in frontal_bigref_keys[1:]]
    all_cingulate_bigref_keys = cingulate_bigref_keys + \
                      ["%s_pchange" % _ for _ in cingulate_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in cingulate_bigref_keys[1:]]
    all_parietal_bigref_keys = parietal_bigref_keys + \
                      ["%s_pchange" % _ for _ in parietal_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in parietal_bigref_keys[1:]]
    all_temporal_bigref_keys = temporal_bigref_keys + \
                      ["%s_pchange" % _ for _ in temporal_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in temporal_bigref_keys[1:]]
    all_av45_key_lists = [all_wm70_composite_keys,
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
    new_headers = rearrangeHeaders(old_headers, all_av45_keys, after = 'LastCSFAbeta')
    
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
                times = exam_timedeltas[:2]
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in bigref_keys[:2]])
                data['AV45_BigRef_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wm70_keys[:2]])
                data['AV45_WM70_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in cerebg_keys[:2]])
                data['AV45_cerebg_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wcereb_keys[:2]])
                data['AV45_wcereb_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in brainstem_keys[:2]])
                data['AV45_brainstem_Slope_2pts'] = slope
            if i == 2:
                times = exam_timedeltas[:3]
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in bigref_keys[:3]])
                data['AV45_BigRef_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wm70_keys[:3]])
                data['AV45_WM70_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in cerebg_keys[:3]])
                data['AV45_cerebg_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wcereb_keys[:3]])
                data['AV45_wcereb_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in brainstem_keys[:3]])
                data['AV45_brainstem_Slope_3pts'] = slope

    return (new_headers, data)


def syncTBMSynData(old_headers, old_lines, tbm_file, registry_file, dump_to=None):
    tbm_by_subj = importTBMSyn(tbm_file)

    # add new headers as needed
    tbm_columns = ['TBMSyn_DATE.%s' % (i+1) for i in range(10)]
    tbm_columns.extend(['TBMSyn_SCORE.%s' % (i+1) for i in range(10)])
    tbm_columns.append('TBMSyn_count')
    tbm_columns.append('TBMSyn_SLOPE')
    new_headers = rearrangeHeaders(old_headers, tbm_columns, after=None)

    # remove date fields
    for i in range(11):
        key = 'TBMSyn_DATE.%s' % (i+1)
        if key in new_headers:
            new_headers.remove(key)
    
    new_lines = []
    for old_l in old_lines:
        # get subject ID
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            continue

        if subj not in tbm_by_subj:
            #print "No TBM data for %s" % subj
            new_lines.append(old_l)
            continue

        subj_tbm = sorted(tbm_by_subj[subj], key=lambda x: x['EXAMDATE'])
        av45_bl, av45_2, av45_3 = getAV45Dates(old_l)
        bl_examdate = subj_tbm[0]['BL_EXAMDATE']
        if av45_bl is None:
            print "No Baseline AV45 date for %s" % subj
            new_lines.append(old_l)
            continue
        # check that TBMSyn baseline is within range of AV45 baseline
        if abs(av45_bl - bl_examdate).days > 210:
            print "TBM BL and AV45 BL %s days apart for subj %s" % (abs(av45_bl - bl_examdate).days, subj)
            new_lines.append(old_l)
            continue

        new_data = {}
        slope_points = []
        for i in range(10):
            if i < len(subj_tbm):
                datapoint = subj_tbm[i]
                new_data['TBMSyn_DATE.%s' % (i+1)] = datapoint['EXAMDATE']
                new_data['TBMSyn_SCORE.%s' % (i+1)] = datapoint['SCORE']
                timediff = (datapoint['EXAMDATE']-datapoint['BL_EXAMDATE']).days / 365.0
                slope_points.append((timediff,datapoint['SCORE']))
            else:
                new_data['TBMSyn_DATE.%s' % (i+1)] = ''
                new_data['TBMSyn_SCORE.%s' % (i+1)] = ''

        # calculate slope
        if len(slope_points) >= 2:
            raw_dates = [_[0] for _ in slope_points]
            raw_scores = [_[1] for _ in slope_points]
            diff_dates = [0.0] + raw_dates
            diff_scores = [0.0] + raw_scores
            #slope, intercept, r, p, stderr = stats.linregress(diff_dates, diff_scores)
            output = optimize.fmin(lambda b, x, y: ((b*x-y)**2).sum(), x0=0.1, args=(diff_dates, diff_scores))
            new_data['TBMSyn_SLOPE'] = output[0]
        else:
            new_data['TBMSyn_SLOPE'] = ''
        new_data['TBMSyn_count'] = len(slope_points)

        new_data = convertToCSVDataType(new_data, decimal_places=5) # TBMSyn needs more decimal places
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def syncCSFData(old_headers, old_lines, csf_files, registry_file, dump_to=None):
    csf_by_subj = importCSF(csf_files)
    registry = importRegistry(registry_file)

    # add new headers as needed
    csf_headers = []
    csf_headers += ['CSF_ABETA.%s' % (i+1) for i in range(7)]
    csf_headers += ['CSF_ABETApostAV45.%s' % (i+1) for i in range(7)]
    csf_headers += ['CSF_ABETA_slope', 'CSF_ABETA_closest_AV45', 
                    'CSF_ABETA_closest_AV45_2', 'CSF_ABETA_closest_AV45_3',
                    'CSF_ABETA_closest_AV45_BIN_192',
                    'CSF_ABETA_closest_AV45_2_BIN_192',
                    'CSF_ABETA_closest_AV45_3_BIN_192']

    csf_headers += ['CSF_TAU.%s' % (i+1) for i in range(7)]
    csf_headers += ['CSF_TAUpostAV45.%s' % (i+1) for i in range(7)]
    csf_headers += ['CSF_TAU_slope', 'CSF_TAU_closest_AV45', 
                    'CSF_TAU_closest_AV45_2', 'CSF_TAU_closest_AV45_3',
                    'CSF_TAU_closest_AV45_BIN_93',
                    'CSF_TAU_closest_AV45_2_BIN_93',
                    'CSF_TAU_closest_AV45_3_BIN_93']

    csf_headers += ['CSF_PTAU.%s' % (i+1) for i in range(7)]
    csf_headers += ['CSF_PTAUpostAV45.%s' % (i+1) for i in range(7)]
    csf_headers += ['CSF_PTAU_slope', 'CSF_PTAU_closest_AV45', 
                    'CSF_PTAU_closest_AV45_2', 'CSF_PTAU_closest_AV45_3',
                    'CSF_PTAU_closest_AV45_BIN_23',
                    'CSF_PTAU_closest_AV45_2_BIN_23',
                    'CSF_PTAU_closest_AV45_3_BIN_23']

    new_headers = rearrangeHeaders(old_headers, csf_headers, after=None)

    lengths = set()
    new_lines = []
    for old_l in old_lines:
        # get subject ID
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            continue
        if subj not in csf_by_subj:
            #print "No TBM data for %s" % subj
            new_lines.append(old_l)
            continue

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
                    examdate = listing['date']
                    break
            if examdate is None:
                print "No exam date for %s (%s, %s)" % (subj, viscode, viscode2)
                continue
            else:
                csf_row['EXAMDATE'] = examdate
                subj_csf_filtered.append(csf_row)
        subj_csf = sorted(subj_csf_filtered, key=lambda x: x['EXAMDATE'])
        lengths.add(len(subj_csf))
        # find av45 dates
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=None)

        new_data = {}
        ptau_slope_points = []
        abeta_slope_points = []
        tau_slope_points = []
        for i in range(7):
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
            else:
                new_data['CSF_ABETA.%s' % (i+1)] = ''
                new_data['CSF_PTAU.%s' % (i+1)] = ''
                new_data['CSF_TAU.%s' % (i+1)] = ''
                new_data['CSF_ABETApostAV45.%s' % (i+1)] = ''
                new_data['CSF_PTAUpostAV45.%s' % (i+1)] = ''
                new_data['CSF_TAUpostAV45.%s' % (i+1)] = ''

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
            raw_dates = [_[0] for _ in ptau_slope_points]
            raw_scores = [_[1] for _ in ptau_slope_points]
            slope, intercept, r, p, stderr = stats.linregress(raw_dates, raw_scores)
            new_data['CSF_PTAU_slope'] = slope
        if len(tau_slope_points) >= 2:
            raw_dates = [_[0] for _ in tau_slope_points]
            raw_scores = [_[1] for _ in tau_slope_points]
            slope, intercept, r, p, stderr = stats.linregress(raw_dates, raw_scores)
            new_data['CSF_TAU_slope'] = slope
        if len(abeta_slope_points) >= 2:
            raw_dates = [_[0] for _ in abeta_slope_points]
            raw_scores = [_[1] for _ in abeta_slope_points]
            slope, intercept, r, p, stderr = stats.linregress(raw_dates, raw_scores)
            new_data['CSF_ABETA_slope'] = slope

        new_data = convertToCSVDataType(new_data, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def syncWMHData(old_headers, old_lines, wmh_file, registry_file, dump_to=None):
    wmh_by_subj = importWMH(wmh_file)
    registry = importRegistry(registry_file)

    to_add_headers = ['WMH_percentOfICV.%s' % (i+1) for i in range(7)]
    to_add_headers += ['WMH_DATE.%s' % (i+1) for i in range(7)]
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    
    new_lines = []
    lengths = set()
    for old_l in old_lines:
        # get subject ID
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            continue

        if subj not in wmh_by_subj:
            new_lines.append(old_l)
            continue

        subj_wmh = sorted(wmh_by_subj[subj], key=lambda x: x['EXAMDATE'])
        # av45_bl, av45_2, av45_3 = getAV45Dates(old_l)
        lengths.add(len(subj_wmh))

        new_data = {}
        for i in range(7):
            if i < len(subj_wmh):
                datapoint = subj_wmh[i]
                new_data['WMH_percentOfICV.%s' % (i+1)] = datapoint['wmh_percent']
                new_data['WMH_DATE.%s' % (i+1)] = datapoint['EXAMDATE']
            else:
                new_data['WMH_percentOfICV.%s' % (i+1)] = ''
                new_data['WMH_DATE.%s' % (i+1)] = ''

        new_data = convertToCSVDataType(new_data, decimal_places=5) # TBMSyn needs more decimal places
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
                 'MMSEclosest_1', 'MMSEclosest_AV45',
                 'datediff_MMSE_AV45', 'MMSE_3mths_AV45',
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
                 'ptau_closest_2']
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
    if old_l['AV45_Date'] != '':
        bl_av45 = datetime.strptime(old_l['AV45_Date'], '%m/%d/%y')
    elif len(patient_pets) >= 1:
        bl_av45 = patient_pets[0]
    else:
        bl_av45 = None
    if old_l['AV45_2_Date'] != '':
        av45_2 = datetime.strptime(old_l['AV45_2_Date'], '%m/%d/%y')
    elif len(patient_pets) >= 2:
        av45_2 = patient_pets[1]
    else:
        av45_2 = None
    if old_l['AV45_3_Date'] != '':
        av45_3 = datetime.strptime(old_l['AV45_3_Date'], '%m/%d/%y')
    elif len(patient_pets) >= 3:
        av45_3 = patient_pets[2]
    else:
        av45_3 = None
    return (bl_av45, av45_2, av45_3)


if __name__ == '__main__':

    '''
    REMEMBER TO DELETE ALL OLD AV45 FIELDS BEFORE RUNNING
    '''

    # IO files
    master_file = "../FDG_AV45_COGdata.csv"
    output_file = "../FDG_AV45_COGdata_synced.csv"
    # LOOKUP files
    registry_file = "../docs/registry_clean.csv"
    arm_file = "../docs/ARM.csv"
    pet_meta_file = "../docs/PET_META_LIST_edited.csv"
    # AV45 File
    av45_file = "../output/UCBERKELEYAV45_07_08_15_extra.csv"
    # MMSE files
    mmse_file = "../cog_tests/MMSE.csv"
    # ADAS-COG files
    adni1_adas_file = '../cog_tests/ADASSCORES.csv'
    adnigo2_adas_file = '../cog_tests/ADAS_ADNIGO2.csv'
    # Neuropsychological Battery file
    neuro_battery_file = '../cog_tests/NEUROBAT.csv'
    # Diagnosis file
    diagnosis_file = '../docs/DXSUM_PDXCONV_ADNIALL.csv'
    # TBMsyn file
    tbm_file = '../mr_docs/Mayo/MAYOADIRL_MRI_TBMSYN_05_07_15.csv'
    # CSF files
    csf_files = ['../docs/UPENNBIOMK5_firstset.csv','../docs/UPENNBIOMK6_secondset.csv','../docs/UPENNBIOMK7_thirdset.csv','../docs/UPENNBIOMK8_fourthset.csv']
    # WMH file
    wmh_file = '../docs/UCD_ADNI2_WMH_03_12_15.csv'

    # syncing pipeline
    new_headers, new_lines = parseCSV(master_file)
    print "\nELIMINATING COLUMNS\n"
    new_headers, new_lines = eliminateColumns(new_headers, new_lines)
    print "\nSYNCING AV45\n"
    new_headers, new_lines = syncAV45Data(new_headers, new_lines, av45_file, registry_file, dump_to=None) # adds new patients
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
for graphing new counts:
[['Test', 'Subj w/ New Test', 'Total New Tests'], ['MMSE', '598', '703'], ['ADAS-cog', '341', '360'], ['AVLT', '318', '333']]


for ARM:
    1=NL - 1.5T only;
    2=MCI - 1.5T only;
    3=AD - 1.5T only;
    4=NL - PET+1.5T;
    5=MCI - PET+1.5T;
    6=AD - PET+1.5T;
    7=NL - 3T+1.5T;
    8=MCI - 3T+1.5T;
    9=AD - 3T+1.5T;

'''
