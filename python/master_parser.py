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
    writer = csv.DictWriter(open(file_path,'rU'), fieldnames=headers)
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

    '''
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
    '''

    new_headers = old_headers # no change in headers
    new_lines = []
    new_values = 0
    total = 0
    for i, old_l in enumerate(old_lines):
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            print "Line %s: Can't convert '%s' to subj" % (i+1, old_l['RID'])
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
                test_score = float(test_results['TOTSCORE']) # CAST TO INT
                diff_from_first = (test_date-first_scan_date).days / 365.0
                all_values.append(test_score)
                all_times.append(diff_from_first)
            count = i+1
            new_subj_data['ADAScog_DATE%s' % count] = test_date_string
            new_subj_data['ADAScog.%s' % count] = test_score
            new_subj_data['TIME_ADAS.%s' % count] = diff_from_first
            if bl_av45 is not None and test_date != '':
                rel_time_days = (test_date - bl_av45).days 
                rel_time = rel_time_days / 365.0
                new_subj_data['TIMEreltoAV45_ADAS.%s' % count] = rel_time
                if rel_time >= 0:
                    post_values.append(test_score)
                    post_times.append(diff_from_first)
                    new_subj_data['TIMEpostAV45_ADAS.%s' % count] = rel_time
                    if rel_time_days <= 93 and 'ADAS_3MTH_AV45' not in new_subj_data:
                        new_subj_data['ADAS_3MTH_AV45'] = test_score
                        new_subj_data['ADAS_3MTHS_AV45DATE'] = test_date_string
            if av45_2 is not None and test_date != '' and test_score != '':
                rel_time_days = (test_date - av45_2).days 
                if rel_time_days >= 0 and rel_time_days <= 93 and 'ADAS_AV45_2_3MTHS' not in new_subj_data:
                    new_subj_data['ADAS_AV45_2_3MTHS'] = test_score
                    new_subj_data['ADAS_AV45_2_DATE'] = test_date_string
        # fill in the blanks
        fill_in = ['ADAS_3MTH_AV45', 'ADAS_3MTHS_AV45DATE', 'ADAS_AV45_2_3MTHS','ADAS_AV45_2_DATE']
        for f_key in fill_in:
            if f_key not in new_subj_data:
                new_subj_data[f_key] = ''
        
        # get slopes
        if len(all_values) >= 2:
            all_times = [float("{0:.2f}".format(_)) for _ in all_times]
            slope, intercept, r, p, stderr = stats.linregress(all_times, all_values)
            print "For %s vs %s" % (all_values, all_times)
            print old_l['ADAS_slope_all']
            print slope
            new_subj_data['ADAS_slope_all'] = slope
        else:
            new_subj_data['ADAS_slope_all'] = ''
        if len(post_values) >= 2:
            post_times = [float("{0:.2f}".format(_)) for _ in post_times]
            slope, intercept, r, p, stderr = stats.linregress(post_times, post_values)
            print "For %s vs %s" % (post_values, post_times)
            print old_l['ADASslope_postAV45']
            print slope
            new_subj_data['ADASslope_postAV45'] = slope
        else:
            new_subj_data['ADASslope_postAV45'] = ''

        '''
        # do comparison:
        print "SUBJ: %s" % subj
        changed = False
        for k in sorted(new_mmse_columns[subj].keys()):
            old_value = old_l[k]
            new_value = new_mmse_columns[subj][k]
            if 'MMSCORE' in k and old_value == '' and new_value != '':
                new_values += 1
                changed = True
            print "\t%s: %s -> %s" % (k, old_value, new_value)
        if changed:
            total += 1
        '''
        old_l.update(new_subj_data)
        new_lines.append(old_l)

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
            print "Line %s: Can't convert '%s' to subj" % (i+1, old_l['RID'])
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
            if 'MMSCORE' in k and old_value == '' and new_value != '':
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

    # syncing pipeline
    old_headers, old_lines = parseCSV(master_file)
    new_headers, new_lines = syncMMSEData(old_headers, old_lines, mmse_file, registry_file, dump_to=None)
    new_headers, new_lines = syncADASCogData(new_headers, new_lines, adni1_adas_file, adnigo2_adas_file, registry_file, dump_to=None)


'''
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

CDR keys

'CDR_sum_boxes'


'''