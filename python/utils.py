import csv
from collections import defaultdict
from datetime import datetime, timedelta

def parseCSV(file_path):
    reader = csv.DictReader(open(file_path,'rU'))
    lines = [l for l in reader]
    headers = reader.fieldnames
    return (headers, lines)

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
            new_data[k] = str(round(new_data[k],decimal_places))
        elif new_data[k] is None:
            new_data[k] = ''
    return new_data

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

def importDODRegistry(dod_registry_file):
    headers, lines = parseCSV(dod_registry_file)
    registry = defaultdict(list)
    for data in lines:
        subj = int(data['SCRNO'])
        viscode = data['VISCODE'].strip().lower()
        if viscode.startswith('sc'):
            continue
        date_string = data['EXAMDATE']
        if date_string == '':
            continue
        date = datetime.strptime(date_string,'%Y-%m-%d')
        registry[subj].append({'VISCODE': viscode,
                               'EXAMDATE': date,
                               'update_stamp': data['update_stamp']})
    registry = dict(registry)
    for k in registry.keys():
        new_val = sorted(registry[k], key=lambda x: x['EXAMDATE'])
        registry[k] = new_val
    return registry


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


def importCSF(csf_files):
    data = defaultdict(dict)
    for csffile in csf_files:
        print csffile
        headers, lines = parseCSV(csffile)
        for i, line in enumerate(lines):
            print line
            subj = int(line['RID'])
            vc = line['VISCODE'].strip().lower()
            try:
                vc2 = line['VISCODE2'].strip().lower()
            except:
                vc2 = ''
            rundate = datetime.strptime(line['RUNDATE'],'%m/%d/%y')
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
                       'abeta': abeta,
                       'tau': tau,
                       'ptau': ptau}
            if visitkey in data[subj]:
                data[subj][visitkey].append(newdata)
            else:
                data[subj][visitkey] = [newdata]
    # take the timepoint with the most recent run date for each visit
    data = dict(data)
    for k in data.keys():
        userdata = data[k]
        flattened = []
        for (vc, vc2), inners in userdata.iteritems():
            chosen = sorted(inners, key=lambda x: x['rundate'])[-1]
            chosen.update({'vc': vc, 'vc2': vc2})
            flattened.append(chosen)
        data[k] = flattened
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


def importMRI(mri_file):
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
        '''
        if line['Sequence'] != 'MPRAGE':
            continue
        '''
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
        
        '''
        if not ('mpr' in seq or 'spgr' in seq or 'n3m' in seq):
            bad_sequences.add(seq)
            continue
        '''
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

        data[subj].append({'EXAMDATE' : date, 'vc': vc})
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
    writer = csv.DictWriter(open(file_path,'w'), fieldnames=headers)
    writer.writeheader()
    for l in lines:
        filtered_line = {}
        for k in headers:
            filtered_line[k] = l[k] if k in l else ''
        writer.writerow(filtered_line)