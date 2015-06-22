'''

Pull in the master FDG_AV45_COG DATA
'''
import csv
from collections import defaultdict
from datetime import datetime, timedelta

def parseCSV(file_path):
    reader = csv.DictReader(open(file_path,'rU'))
    lines = [l for l in reader]
    headers = reader.fieldnames
    return (headers, lines)

def importRegistry(registry_file):
    headers, lines = parseCSV(registry_file)
    registry = defaultdict(dict)
    for data in lines:
        if data['EXAMDATE'] == '' or data['VISCODE'] in set(['sc', 'scmri']) or data['VISCODE2'] in set(['sc', 'scmri']):
            continue
        subj = int(data['RID'])
        date = datetime.strptime(data['EXAMDATE'],'%Y-%m-%d')
        registry[subj][date] = {'VISCODE': data['VISCODE'].strip().lower(),
                                'VISCODE2': data['VISCODE2'].strip().lower(),
                                'update_stamp': data['update_stamp']}
    return registry

def dumpCSV(file_path, headers, lines):
    writer = csv.DictWriter(open(file_path,'rU'), fieldnames=headers)
    for l in lines:
        writer.writerow(l)

def syncMMSEData(old_headers, old_lines, mmse_file, registry_file, dump_to=None):
    relevant_old_headers = ['MMSE_DATE1','MMSE_DATE2','MMSE_DATE3','MMSE_DATE4','MMSE_DATE5',
                            'MMSE_DATE6','MMSE_DATE7','MMSE_DATE8','MMSE_DATE9','MMSE_DATE10','MMSE_DATE11',
                            'MMSCORE.1','MMSCORE.2','MMSCORE.3','MMSCORE.4','MMSCORE.5',
                            'MMSCORE.6','MMSCORE.7','MMSCORE.8','MMSCORE.9','MMSCORE.10','MMSCORE.11',
                            'MMSEclosest_1','datediff_MMSE_AV45','MMSEclosest_AV45','MMSE_3mths_AV45']

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
            for k,v in subject_registry.iteritems():
                if vs == v['VISCODE'] and vs2 == v['VISCODE2']:
                    date = k
            if date == None:
                raise Exception("Could not find visit in registry: %s, %s, %s" % (subj, vs, vs2))
        else:
            date = datetime.strptime(date_string,'%Y-%m-%d')
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
                test_score = int(test_results['MMSCORE'])
            count = i+1
            new_subj_data['MMSE_DATE%s' % count] = test_date
            new_subj_data['MMSCORE.%s' % count] = test_score
        # find date of closest av45 scan, and corresponding cog test
        pass
        # save
        new_mmse_columns[subj] = new_subj_data

    # overwrite old lines
    new_headers = old_headers # no change in headers
    new_lines = []
    for old_l in old_lines:
        subj = int(old_l['RID'])
        old_l.update(new_mmse_columns[subj])
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return new_headers, new_lines

if __name__ == '__main__':
    master_file = "../FDG_AV45_COGdata.csv"
    output_file = "../FDG_AV45_COGdata_synced.csv"
    mmse_file = "../cog_tests/MMSE.csv"
    registry_file = "../docs/registry_clean.csv"

    # syncing pipeline
    old_headers, old_lines = parseCSV(master_file)
    new_headers, new_lines = syncMMSEData(old_headers, old_lines, mmse_file, dump_to=None)

    #dumpCSV(output_file, new_headers, new_lines)

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
'AVLT_AV45_2_DATE'
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