'''

Pull in the master FDG_AV45_COG DATA

remove null char:
tr -d '\000' < file1 > file2

'''
import csv
import re
from collections import defaultdict
from datetime import datetime, timedelta
from scipy import stats

def parseCSV(file_path):
    reader = csv.DictReader(open(file_path,'rU'))
    lines = [l for l in reader]
    headers = reader.fieldnames
    return (headers, lines)

def importRegistry(registry_file):
    headers, lines = parseCSV(registry_file)
    registry = defaultdict(list)
    for data in lines:
        if data['EXAMDATE'] == '':
            continue
        subj = int(data['RID'])
        date = None
        try:
            date = datetime.strptime(data['EXAMDATE'],'%Y-%m-%d')
        except Exception as e:
            pass
        try:
            date = datetime.strptime(data['EXAMDATE'],'%m/%d/%y')
        except Exception as e:
            pass
        if date is None:
            continue
        registry[subj].append({'VISCODE': data['VISCODE'].strip().lower(),
                               'VISCODE2': data['VISCODE2'].strip().lower(),
                               'date': date,
                               'update_stamp': data['update_stamp']})
    return registry

def importPetMETA(pet_meta_file):
    headers, lines = parseCSV(pet_meta_file)
    pets = defaultdict(list)
    for row in lines:
        subj = int(row['Subject'].split('_')[-1].strip())
        new_date = datetime.strptime(row['Scan Date'], '%m/%d/%y')
        pets[subj].append(new_date)
    return dict(pets)

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
        status = data['ARM'].strip()
        if status == '':
            continue
        status = int(status)
        userdate = datetime.strptime(data['USERDATE'],'%Y-%m-%d')
        # convert status
        status_str = translation[status]
        arms[subj].append({'USERDATE': userdate,
                           'STATUS': status_str})
    return dict(arms)

def dumpCSV(file_path, headers, lines):
    writer = csv.DictWriter(open(file_path,'w'), fieldnames=headers)
    writer.writeheader()
    for l in lines:
        filtered_line = {k:v for k,v in l.iteritems() if k in headers}
        writer.writerow(filtered_line)


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

    # add some new headers
    if 'ADAS_AV45_3_DATE' not in new_headers or 'ADAS_AV45_3_3MTHS' not in new_headers:
        if 'ADAS_AV45_3_3MTHS' in new_headers:
            # remove it
            new_headers.remove('ADAS_AV45_3_3MTHS')
        elif 'ADAS_AV45_3_DATE' in new_headers:
            # remove it
            new_headers.remove('ADAS_AV45_3_DATE')
        idx = new_headers.index('ADAS_AV45_2_DATE')
        new_headers.insert(idx+1, 'ADAS_AV45_3_DATE')
        new_headers.insert(idx+1, 'ADAS_AV45_3_3MTHS')
    if 'ADAS_post_AV45_followuptime' not in new_headers:
        idx = new_headers.index('ADAS_3MTH_AV45')
        new_headers.insert(idx, 'ADAS_post_AV45_followuptime')


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


        # Get AV45 Scan dates
        if old_l['AV45_Date'] != '':
            bl_av45 = datetime.strptime(old_l['AV45_Date'], '%m/%d/%y')
        else:
            bl_av45 = None
        if old_l['AV45_2_Date'] != '':
            av45_2 = datetime.strptime(old_l['AV45_2_Date'], '%m/%d/%y')
        else:
            av45_2 = None
        if old_l['AV45_3_Date'] != '':
            av45_3 = datetime.strptime(old_l['AV45_3_Date'], '%m/%d/%y')
        else:
            av45_3 = None


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
                test_score = test_results['TOTSCORE']
                if test_score != '':
                    test_score = round(float(test_score),2) # CAST TO INT
                diff_from_first = round((test_date-first_scan_date).days / 365.0, 2)
                if test_score != '':
                    all_values.append(test_score)
                    all_times.append(diff_from_first)
            count = i+1
            new_subj_data['ADAScog_DATE%s' % count] = test_date_string
            new_subj_data['ADAScog.%s' % count] = test_score
            new_subj_data['TIME_ADAS.%s' % count] = diff_from_first
            if bl_av45 is not None and test_date != '':
                rel_time_days = (test_date - bl_av45).days
                rel_time = round(rel_time_days / 365.0, 2)
                new_subj_data['TIMEreltoAV45_ADAS.%s' % count] = rel_time
                if abs(rel_time_days) <= 93 and 'ADAS_3MTH_AV45' not in new_subj_data:
                    new_subj_data['ADAS_3MTH_AV45'] = test_score
                    new_subj_data['ADAS_3MTHS_AV45DATE'] = test_date_string
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
                    new_subj_data['ADAS_AV45_2_DATE'] = test_date_string
            if av45_3 is not None and test_date != '':
                rel_time_days = (test_date - av45_3).days 
                if abs(rel_time_days) <= 93 and 'ADAS_AV45_3_3MTHS' not in new_subj_data:
                    new_subj_data['ADAS_AV45_3_3MTHS'] = test_score
                    new_subj_data['ADAS_AV45_3_DATE'] = test_date_string
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
            new_subj_data['ADAS_slope_all'] = round(slope,2)
        else:
            new_subj_data['ADAS_slope_all'] = ''
        if len(post_values) >= 2:
            slope, intercept, r, p, stderr = stats.linregress(post_times, post_values)
            new_subj_data['ADASslope_postAV45'] = round(slope,2)
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
        
        old_l.update(new_subj_data)
        new_lines.append(old_l)

    print "Total subj changed: %s" % total
    print "Total new tests: %s" % new_values

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def syncMMSEData(old_headers, old_lines, mmse_file, registry_file, dump_to=None):
    mmse_headers, mmse_lines = parseCSV(mmse_file)
    registry = importRegistry(registry_file)
    
    # restructure mmse lines by subject
    mmse_by_subject = defaultdict(list)
    for line in mmse_lines:
        subj = int(line['RID'])
        date_string = line['EXAMDATE']
        if not date_string:
            # get from registry
            vs = line['VISCODE'].strip().lower()
            vs2 = line['VISCODE2'].strip().lower()
            subject_registry = registry[subj]
            date = None
            for v in subject_registry:
                item_date = v['date']
                if vs2 == v['VISCODE2'] and item_date != '' and item_date is not None:
                    date = item_date
                    break
            if date == None:
                print "Could not find visit in registry: %s, %s, %s" % (subj, vs, vs2)
                for _ in subject_registry:
                    print _
        else:
            date = datetime.strptime(date_string,'%m/%d/%y')
        if date is not None:
            mmse_by_subject[subj].append((date,line))
    
    mmse_by_subject = dict(mmse_by_subject)
    for k,v in mmse_by_subject.iteritems():
        mmse_by_subject[k] = sorted(v, key=lambda x: x[0])

    # Fill in MMSE columns for each subject
    new_mmse_columns = {}
    for subj, tests in mmse_by_subject.iteritems():
        new_subj_data = {}
        for i in range(11):
            if i >= len(tests):
                test_date = ''
                test_score = ''
            else:
                test_date, test_results = tests[i]
                test_date = test_date.strftime("%m/%d/%y")
                if test_results['MMSCORE'] == '':
                    test_score = ''
                else:
                    test_score = int(test_results['MMSCORE'])
                    if test_score < 0:
                        test_score = ''

            count = i+1
            new_subj_data['MMSE_DATE%s' % count] = test_date
            new_subj_data['MMSCORE.%s' % count] = test_score
        # save
        new_mmse_columns[subj] = new_subj_data

    # overwrite old lines
    new_headers = old_headers # no change in headers
    new_lines = []
    new_values = 0
    total = 0
    for i, old_l in enumerate(old_lines):
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            continue
        #print "%s: %s" % (subj, new_mmse_columns[subj])

        # add av45-relative columns
        bl_date_string = old_l['AV45_Date']
        date_used = False
        if bl_date_string != '':
            date_used = True
            bl_av45 = datetime.strptime(old_l['AV45_Date'], '%m/%d/%y')
            closest_mmse_date, closest_results = sorted(mmse_by_subject[subj], key=lambda x: abs(x[0] - bl_av45))[0]
        else:
            closest_mmse_date, closest_results = sorted(mmse_by_subject[subj], key=lambda x: x[0])[0]
        test_score = closest_results['MMSCORE']
        if test_score != '':
            test_score = int(test_score)
        new_mmse_columns[subj]['MMSEclosest_1'] = closest_mmse_date.strftime('%m/%d/%y')
        new_mmse_columns[subj]['MMSEclosest_AV45'] = test_score
        if date_used:
            date_diff = abs(bl_av45 - closest_mmse_date).days
            new_mmse_columns[subj]['datediff_MMSE_AV45'] = abs(bl_av45 - closest_mmse_date).days / 365.25
            new_mmse_columns[subj]['MMSE_3mths_AV45'] = test_score if date_diff < 93 else ''
        else:
            new_mmse_columns[subj]['datediff_MMSE_AV45'] = ''
            new_mmse_columns[subj]['MMSE_3mths_AV45'] = ''

        # do comparison:
        #print "SUBJ: %s" % subj
        changed = False
        for k in sorted(new_mmse_columns[subj].keys()):
            old_value = old_l[k]
            new_value = new_mmse_columns[subj][k]
            if k.startswith('MMSCORE') and old_value == '' and new_value != '':
                new_values += 1
                changed = True
            #print "\t%s: %s -> %s" % (k, old_value, new_value)
        if changed:
            total += 1

        old_l.update(new_mmse_columns[subj])
        new_lines.append(old_l)

    print "%s subj with new tests (%s new tests)" % (total, new_values)
    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return new_headers, new_lines

def syncAVLTData(old_headers, old_lines, neuro_battery_file, registry_file, dump_to=None):
    neurobat_headers, neurobat_lines = parseCSV(neuro_battery_file)
    registry = importRegistry(registry_file)

    # restructure by subject
    avlt_by_subj = defaultdict(list)
    for line in neurobat_lines:
        subj = int(line['RID'])
        viscode = line['VISCODE'].strip().lower()
        viscode2 = line['VISCODE2'].strip().lower()
        examdate = line['EXAMDATE']
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
            #print "Invalid scores found for %s (%s, %s): %s" % (subj, viscode, viscode2, tots)
            continue

        avlt_by_subj[subj].append({'VISCODE': viscode,
                                   'VISCODE2': viscode2,
                                   'EXAMDATE': examdate,
                                   'TOTS': test_score})

    new_headers = old_headers # no change in headers

    # add potential new headers
    if 'AVLT_post_AV45_followuptime' not in new_headers:
        idx = new_headers.index('AVLT_3MTHS_AV45')
        new_headers.insert(idx, 'AVLT_post_AV45_followuptime')

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
        # Get AV45 Scan dates
        if old_l['AV45_Date'] != '':
            bl_av45 = datetime.strptime(old_l['AV45_Date'], '%m/%d/%y')
        else:
            bl_av45 = None
        if old_l['AV45_2_Date'] != '':
            av45_2 = datetime.strptime(old_l['AV45_2_Date'], '%m/%d/%y')
        else:
            av45_2 = None


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
                diff_from_first = round((test_date-first_scan_date).days / 365.0, 2)
                if test_score != '':
                    all_values.append(test_score) 
                    all_times.append(diff_from_first)
            count = i+1
            new_subj_data['AVLT_DATE.%s' % count] = test_date_string
            new_subj_data['AVLT.%s' % count] = test_score
            new_subj_data['TIME_AVLT.%s' % count] = diff_from_first
            if bl_av45 is not None and test_date != '':
                rel_time_days = (test_date - bl_av45).days
                rel_time = round(rel_time_days / 365.0, 2)
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
            new_subj_data['AVLT_slope_all'] = round(slope,2)
        else:
            new_subj_data['AVLT_slope_all'] = ''
        if len(post_values) >= 2:
            slope, intercept, r, p, stderr = stats.linregress(post_times, post_values)
            new_subj_data['AVLTslope_postAV45'] = round(slope,2)
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
        if old_l['AV45_Date'] != '':
            bl_av45 = datetime.strptime(old_l['AV45_Date'], '%m/%d/%y')
        elif len(patient_pets) >= 1:
            bl_av45 = patient_pets[0]
        else:
            bl_av45 = ''
        if old_l['AV45_2_Date'] != '':
            av45_2 = datetime.strptime(old_l['AV45_2_Date'], '%m/%d/%y')
        elif len(patient_pets) >= 2:
            av45_2 = patient_pets[1]
        else:
            av45_2 = ''
        if old_l['AV45_3_Date'] != '':
            av45_3 = datetime.strptime(old_l['AV45_3_Date'], '%m/%d/%y')
        elif len(patient_pets) >= 3:
            av45_3 = patient_pets[2]
        else:
            av45_3 = ''


        
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
        
        if bl_av45 != '':
            av45_bl_closest = sorted(subj_diags, key=lambda x: abs(bl_av45-x['EXAMDATE']))[0]
            if abs(bl_av45 - av45_bl_closest['EXAMDATE']).days < 550:
                new_data['Diag@AV45_long'] = av45_bl_closest['diag']
        if av45_2 != '':
            av45_2_closest = sorted(subj_diags, key=lambda x: abs(av45_2-x['EXAMDATE']))[0]
            if abs(av45_2 - av45_2_closest['EXAMDATE']).days < 550:
                new_data['Diag@AV45_2_long'] = av45_2_closest['diag']
        if av45_3 != '':
            av45_3_closest = sorted(subj_diags, key=lambda x: abs(av45_3-x['EXAMDATE']))[0]
            if abs(av45_3 - av45_3_closest['EXAMDATE']).days < 550:
                new_data['Diag@AV45_3_long'] = av45_3_closest['diag']     


        for i, diag_row in enumerate(subj_diags):
            # parse change based on init diag
            if diag_row['change'] == 5:
                new_data['MCItoADConv(fromav45)'] = '1'
                new_data['MCItoADConvDate'] = new_data['MCItoADconv_'] = diag_row['EXAMDATE']
                new_data['Baseline_MCItoAD_ConvTime'] = (diag_row['EXAMDATE'] - new_data['Baseline']).days / 365.0
                if bl_av45 != '':
                    new_data['AV45_MCItoAD_ConvTime'] = (diag_row['EXAMDATE'] - bl_av45).days / 365.0
                
                    
        # convert datetimes to strings
        for k in new_data.keys():
            if isinstance(new_data[k], datetime):
                new_data[k] = datetime.strftime(new_data[k], '%m/%d/%y')
            elif isinstance(new_data[k], float):
                new_data[k] = str(round(new_data[k],2))
            elif new_data[k] is None:
                new_data[k] = ''
        
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
    av45_headers, av45_lines = parseCSV(av45_file)
    registry = importRegistry(registry_file)

    av45_by_subj = defaultdict(list)
    for line in av45_lines:
        subj = int(line.pop('RID',None))
        viscode = line['VISCODE'].strip().lower()
        viscode2 = line['VISCODE2'].strip().lower()
        examdate = line.get('EXAMDATE',None)
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
        av45_by_subj[subj].append(line)

    new_headers = old_headers
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

        new_data = parseAV45Entries(av45_by_subj[subj])
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
        new_columns = parseAV45Entries(av45_by_subj[ns])
        new_subj_row.update(new_columns)
        new_lines.append(new_subj_row)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def parseAV45Entries(subj_av45):
    pass


def syncTBMSynData(old_headers, old_lines, tbm_file, registry_file, dump_to=None):
    tbm_headers, tbm_lines = parseCSV(tbm_file)
    registry = importRegistry(registry_file)

    tbm_by_subj = defaultdict(list)
    for line in tbm_lines:
        subj = int(line.pop('RID',None))
        viscode = line['VISCODE'].strip().lower()
        viscode2 = line['VISCODE2'].strip().lower()
        examdate = line.get('EXAMDATE',None)
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
        av45_by_subj[subj].append(line)

    new_headers = old_headers
    new_lines = []
    for old_l in old_lines:
        # get subject ID
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            continue

        if subj not in tbm_by_subj:
            print "No TBM data for %s" % subj
            new_lines.append(old_l)
            continue

        new_data = {}
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)


def eliminateColumns(headers, lines):
    to_remove = ['LastClinicalVisit',
                 'DIFF_LastClinicalVisit_AV45',
                 'ClinicalVisitClosestto_AV45',
                 'ClinicalVisitClosestto_AV45_2']
    for tm in to_remove:
        if tm in headers:
            headers.remove(tm)
    new_lines = []
    for l in lines:
        for tm in to_remove:
            l.pop(tm, None)
        new_lines.append(l)
    return headers, new_lines


if __name__ == '__main__':
    # Input/output/lookup files
    master_file = "../FDG_AV45_COGdata.csv"
    output_file = "../FDG_AV45_COGdata_synced.csv"
    registry_file = "../docs/registry_clean.csv"
    arm_file = "../docs/ARM.csv"
    pet_meta_file = "../docs/PET_META_LIST_edited.csv"
    # AV45 File
    av45_file = "../output/UCBERKELEYAV45_06_25_15.csv"
    # MMSE files
    mmse_file = "../cog_tests/MMSE.csv"
    # ADAS-COG files
    adni1_adas_file = '../cog_tests/ADASSCORES.csv'
    adnigo2_adas_file = '../cog_tests/ADAS_ADNIGO2.csv'
    # Neuropsychological Battery file
    neuro_battery_file = '../cog_tests/NEUROBAT.csv'
    # Diagnosis file
    diagnosis_file = '../docs/DXSUM_PDXCONV_ADNIALL.csv'


    # syncing pipeline
    new_headers, new_lines = parseCSV(master_file)
    new_headers, new_lines = eliminateColumns(new_headers, new_lines)
    new_headers, new_lines = syncAV45Data(new_headers, new_lines, av45_file, registry_file, dump_to=None)
    new_headers, new_lines = syncDiagnosisData(new_headers, new_lines, diagnosis_file, registry_file, arm_file, pet_meta_file, dump_to=None)
    new_headers, new_lines = syncMMSEData(new_headers, new_lines, mmse_file, registry_file, dump_to=None)
    new_headers, new_lines = syncADASCogData(new_headers, new_lines, adni1_adas_file, adnigo2_adas_file, registry_file, dump_to=None)
    new_headers, new_lines = syncAVLTData(new_headers, new_lines, neuro_battery_file, registry_file, dump_to=output_file)
    


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
