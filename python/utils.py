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
                               'EXAMDATE': date,
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
    return dict(diag_by_subj)


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
            subj = int(line['RID'])
        except Exception as e:
            continue
        vc = line['VISCODE'].strip().lower()
        vc2 = line['VISCODE2'].strip().lower()
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
                examdate = findVisitDate(registry, subj, vc, vc2) if registry is not None else Nones
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


def importAV45(av45_file, registry=None):
    av45_headers, av45_lines = parseCSV(av45_file)
    av45_by_subj = defaultdict(list)
    for line in av45_lines:
        if 'RID' in line:
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
            continue
        line['EXAMDATE'] = examdate
        av45_by_subj[subj].append(line)
    return dict(av45_by_subj)

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


if __name__ == "__main__":
    lut_file = "../FreeSurferColorLUT.txt"
    data = importFreesurferLookup(lut_file)
    to_lookup = [5,14,15,18,24,26,28,30,31,44,54,58,60,62,63,72,77,80,85,
                 251,252,253,254,255,1000,1001,1004,1005,1007,1009,1010,
                 1013,1016,1017,1021,1022,1024,1026,1030,1033,1034,1035,
                 2000,2001,2004,2005,2007,2009,2010,2013,2016,2017,2021,
                 2022,2024,2026,2030,2033,2034,2035]
    for idx in to_lookup:
        print "%s: %s" % (idx, data[idx])


