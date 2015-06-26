'''

Pull in the master FDG_AV45_COG DATA

remove null char:
tr -d '\000' < file1 > file2

'''
import csv
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
        date = datetime.strptime(data['EXAMDATE'],'%Y-%m-%d')
        registry[subj].append({'VISCODE': data['VISCODE'].strip().lower(),
                               'VISCODE2': data['VISCODE2'].strip().lower(),
                               'date': date,
                               'update_stamp': data['update_stamp']})
    return registry

def dumpCSV(file_path, headers, lines):
    writer = csv.DictWriter(open(file_path,'w'), fieldnames=headers)
    writer.writeheader()
    for l in lines:
        writer.writerow(l)


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

    new_headers = old_headers # no change in headers
    new_lines = []
    new_values = 0
    total = 0
    for linenum, old_l in enumerate(old_lines):
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            #print "Line %s: Can't convert '%s' to subj" % (linenum+1, old_l['RID'])
            new_lines.append(old_l)
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


        new_subj_data = {}
        all_values = []
        post_values = []
        all_times = []
        post_times = []
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
                test_score = round(float(test_results['TOTSCORE']),2) # CAST TO INT
                diff_from_first = round((test_date-first_scan_date).days / 365.0, 2)
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
                    post_values.append(test_score)
                    post_times.append(diff_from_first)
                    new_subj_data['TIMEpostAV45_ADAS.%s' % count] = rel_time
                else:
                    new_subj_data['TIMEpostAV45_ADAS.%s' % count] = ''
            if av45_2 is not None and test_date != '':
                rel_time_days = (test_date - av45_2).days 
                if abs(rel_time_days) <= 93 and 'ADAS_AV45_2_3MTHS' not in new_subj_data:
                    new_subj_data['ADAS_AV45_2_3MTHS'] = test_score
                    new_subj_data['ADAS_AV45_2_DATE'] = test_date_string
        # fill in the blanks
        fill_in = ['ADAS_3MTH_AV45', 'ADAS_3MTHS_AV45DATE', 'ADAS_AV45_2_3MTHS','ADAS_AV45_2_DATE']
        for f_key in fill_in:
            if f_key not in new_subj_data:
                new_subj_data[f_key] = ''
        
        # get slopes
        if len(all_values) >= 2:
            slope, intercept, r, p, stderr = stats.linregress(all_times, all_values)
            '''
            print "For %s vs %s" % (all_values, all_times)
            print old_l['ADAS_slope_all']
            print slope
            '''
            new_subj_data['ADAS_slope_all'] = round(slope,2)
        else:
            new_subj_data['ADAS_slope_all'] = ''
        if len(post_values) >= 2:
            slope, intercept, r, p, stderr = stats.linregress(post_times, post_values)
            '''
            print "For %s vs %s" % (post_values, post_times)
            print old_l['ADASslope_postAV45']
            print slope
            '''
            new_subj_data['ADASslope_postAV45'] = round(slope,2)
        else:
            new_subj_data['ADASslope_postAV45'] = ''

        
        # do comparison:
        print "SUBJ: %s" % subj
        print "first: %s, second: %s" % (bl_av45, av45_2)
        changed = False
        for k in sorted(new_subj_data.keys()):
            old_value = old_l[k]
            new_value = new_subj_data[k]
            if k.startswith('ADAScog.') and old_value == '' and new_value != '':
                new_values += 1
                changed = True
            print "\t%s: %s -> %s" % (k, old_value, new_value)
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
            #print "Line %s: Can't convert '%s' to subj" % (i+1, old_l['RID'])
            new_lines.append(old_l)
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
            print "Invalid scores found for %s (%s, %s): %s" % (subj, viscode, viscode2, tots)
            continue

        avlt_by_subj[subj].append({'VISCODE': viscode,
                                   'VISCODE2': viscode2,
                                   'EXAMDATE': examdate,
                                   'TOTS': test_score})

    new_headers = old_headers # no change in headers
    new_lines = []
    new_values = 0
    total = 0
    subj_new_counts = {}
    for linenum, old_l in enumerate(old_lines):
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            #print "Line %s: Can't convert '%s' to subj" % (linenum+1, old_l['RID'])
            new_lines.append(old_l)
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
                else:
                    new_subj_data['TIMEpostAV45_AVLT.%s' % count] = ''
            if av45_2 is not None and test_date != '':
                rel_time_days = (test_date - av45_2).days 
                if abs(rel_time_days) <= 93 and 'AVLT_AV45_2_3MTHS' not in new_subj_data:
                    new_subj_data['AVLT_AV45_2_3MTHS'] = test_score
                    new_subj_data['AVLT_AV45_2_DATE'] = test_date_string
        # fill in the blanks
        fill_in = ['AVLT_3MTHS_AV45', 'AVLT_3MTHSAV45_Date', 'AVLT_AV45_2_3MTHS','AVLT_AV45_2_DATE']
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
            '''
            print "For %s vs %s" % (all_values, all_times)
            print old_l['AVLT_slope_all']
            print slope
            '''
            new_subj_data['AVLT_slope_all'] = round(slope,2)
        else:
            new_subj_data['AVLT_slope_all'] = ''
        if len(post_values) >= 2:
            slope, intercept, r, p, stderr = stats.linregress(post_times, post_values)
            '''
            print "For %s vs %s" % (post_values, post_times)
            print old_l['AVLTslope_postAV45']
            print slope
            '''
            new_subj_data['AVLTslope_postAV45'] = round(slope,2)
        else:
            new_subj_data['AVLTslope_postAV45'] = ''

        # do comparison:
        print "SUBJ: %s" % subj
        print "first: %s, second: %s" % (bl_av45, av45_2)
        print "first date: %s" % (first_scan_date)
        changed = False
        new_for_subj = 0
        for k in sorted(new_subj_data.keys()):
            old_value = old_l[k]
            new_value = new_subj_data[k]
            if k.startswith('AVLT.') and not bool(old_value) and bool(new_value):
                new_values += 1
                new_for_subj += 1
                changed = True
            print "\t%s: %s -> %s" % (k, old_value, new_value)
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


if __name__ == '__main__':
    # Input/output/lookup files
    master_file = "../FDG_AV45_COGdata.csv"
    output_file = "../FDG_AV45_COGdata_synced.csv"
    registry_file = "../docs/registry_clean.csv"
    # MMSE files
    mmse_file = "../cog_tests/MMSE.csv"
    # ADAS-COG files
    adni1_adas_file = '../cog_tests/ADASSCORES.csv'
    adnigo2_adas_file = '../cog_tests/ADAS_ADNIGO2.csv'
    # Neuropsychological Battery file
    neuro_battery_file = '../cog_tests/NEUROBAT.csv'

    # syncing pipeline
    old_headers, old_lines = parseCSV(master_file)
    new_headers, new_lines = syncMMSEData(old_headers, old_lines, mmse_file, registry_file, dump_to=None)
    new_headers, new_lines = syncADASCogData(new_headers, new_lines, adni1_adas_file, adnigo2_adas_file, registry_file, dump_to=None)
    new_headers, new_lines = syncAVLTData(new_headers, new_lines, neuro_battery_file, registry_file, dump_to=output_file)

'''
for graphing new counts:
[['Test', 'Subj w/ New Test', 'Total New Tests'], ['MMSE', '598', '703'], ['ADAS-cog', '341', '360'], ['AVLT', '318', '333']]


ADAS keys:

'ADAScog_DATE1'
'ADAScog_DATE2'
'ADAScog_DATE3'
'ADAScog_DATE4'
'ADAScog_DATE5'
'ADAScog_DATE6'
'ADAScog_DATE7'
'ADAScog_DATE8'
'ADAScog_DATE9'
'ADAScog_DATE10'
'ADAScog_DATE11'
'ADAScog.1'
'ADAScog.2'
'ADAScog.3'
'ADAScog.4'
'ADAScog.5'
'ADAScog.6'
'ADAScog.7'
'ADAScog.8'
'ADAScog.9'
'ADAScog.10'
'ADAScog.11'
'TIME_ADAS.1'
'TIME_ADAS.2'
'TIME_ADAS.3'
'TIME_ADAS.4'
'TIME_ADAS.5'
'TIME_ADAS.6'
'TIME_ADAS.7'
'TIME_ADAS.8'
'TIME_ADAS.9'
'TIME_ADAS.10'
'TIME_ADAS.11'
'TIMEreltoAV45_ADAS.1'
'TIMEreltoAV45_ADAS.2'
'TIMEreltoAV45_ADAS.3'
'TIMEreltoAV45_ADAS.4'
'TIMEreltoAV45_ADAS.5'
'TIMEreltoAV45_ADAS.6'
'TIMEreltoAV45_ADAS.7'
'TIMEreltoAV45_ADAS.8'
'TIMEreltoAV45_ADAS.9'
'TIMEreltoAV45_ADAS.10'
'TIMEreltoAV45_ADAS.11'
'TIMEpostAV45_ADAS.1'
'TIMEpostAV45_ADAS.2'
'TIMEpostAV45_ADAS.3'
'TIMEpostAV45_ADAS.4'
'TIMEpostAV45_ADAS.5'
'TIMEpostAV45_ADAS.6'
'TIMEpostAV45_ADAS.7'
'TIMEpostAV45_ADAS.8'
'TIMEpostAV45_ADAS.9'
'TIMEpostAV45_ADAS.10'
'TIMEpostAV45_ADAS.11'
'ADAS_3MTH_AV45'
'ADAS_3MTHS_AV45DATE'
'ADAS_AV45_2_3MTHS'
'ADAS_AV45_2_DATE'
'ADAS_slope_all'
'ADASslope_postAV45'

AVLT keys

'AVLT_DATE.1'
'AVLT_DATE.2'
'AVLT_DATE.3'
'AVLT_DATE.4'
'AVLT_DATE.5'
'AVLT_DATE.6'
'AVLT_DATE.7'
'AVLT_DATE.8'
'AVLT_DATE.9'
'AVLT_DATE.10'
'AVLT_DATE.11'
'AVLT.1'
'AVLT.2'
'AVLT.3'
'AVLT.4'
'AVLT.5'
'AVLT.6'
'AVLT.7'
'AVLT.8'
'AVLT.9'
'AVLT.10'
'AVLT.11'
'TIME_AVLT.1'
'TIME_AVLT.2'
'TIME_AVLT.3'
'TIME_AVLT.4'
'TIME_AVLT.5'
'TIME_AVLT.6'
'TIME_AVLT.7'
'TIME_AVLT.8'
'TIME_AVLT.9'
'TIME_AVLT.10'
'TIME_AVLT.11'
'TIMEreltoAV45_AVLT.1'
'TIMEreltoAV45_AVLT.2'
'TIMEreltoAV45_AVLT.3'
'TIMEreltoAV45_AVLT.4'
'TIMEreltoAV45_AVLT.5'
'TIMEreltoAV45_AVLT.6'
'TIMEreltoAV45_AVLT.7'
'TIMEreltoAV45_AVLT.8'
'TIMEreltoAV45_AVLT.9'
'TIMEreltoAV45_AVLT.10'
'TIMEreltoAV45_AVLT.11'
'TIMEpostAV45_AVLT.1'
'TIMEpostAV45_AVLT.2'
'TIMEpostAV45_AVLT.3'
'TIMEpostAV45_AVLT.4'
'TIMEpostAV45_AVLT.5'
'TIMEpostAV45_AVLT.6'
'TIMEpostAV45_AVLT.7'
'TIMEpostAV45_AVLT.8'
'TIMEpostAV45_AVLT.9'
'TIMEpostAV45_AVLT.10'
'TIMEpostAV45_AVLT.11'
'AVLT_3MTHS_AV45'
'AVLT_3MTHSAV45_Date'
'AVLT_AV45_2_3MTHS'
'AVLT_AV45_2_DATE'
'PostAV45Followup'
'AVLT_slope_all'
'AVLTslope_postAV45'

'''