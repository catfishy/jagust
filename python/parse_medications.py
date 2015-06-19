'''
it's a list of peoples' medications, current and past.
the question is whether the use of certain medications at baseline is associated with subsequent 
cognitive decline (and potentially amyloid positivity).  

recent epidemiological studies showing that use of benzodiazepines & antihistamines
are associated with future onset of dementia.  

we have to deal with both misspellings and use of different names for the same medication (i.e. alprazolam/xanax).  

reduce this huge amount of medication data into a few variables representing whether people 
are on a benzo at their baseline visit 
(probably using the baseline amyloid scan date which is in that big excel file in the 
shared ADNI dropbox folder called FDG_AV45_COGdata), and some way to categorize dosage/frequency info.  

'''


import csv
from collections import defaultdict
from datetime import datetime
import re

KEYMED_TRANSLATE = {1:'Aricept',
                    2:'Cognex',
                    3:'Exelon',
                    4:'Namenda',
                    5:'Razadyne',
                    6:'Anti-depressant medication',
                    7:'Other behavioral medication',
                    0:'None of the above'}

REASON_TRANSLATE = {1:'Adverse Event',
                    2:'Therapeutic Use',
                    3:'Prophylaxis/Non-therapeutic Use'}

def parseBackMeds(backmed_file):
    reader = csv.DictReader(open(backmed_file, 'rU'))
    subj_dict = defaultdict(list)
    for l in reader:
        subj = int(l['RID'])
        new_data = {}
        new_data['VISCODE'] = l['VISCODE']
        new_data['VISCODE2'] = l['VISCODE2']
        new_data['INPUTDATE'] = datetime.strptime(l['USERDATE'],'%Y-%m-%d')
        new_data['KEYMED'] = map(int, l['KEYMED'].split(':')) # At this visit, is participant on any of the following medication?, 1=Aricept;2=Cognex;3=Exelon;4=Namenda;5=Razadyne;6=Anti-depressant medication;7=Other behavioral medication;0=None of the above
        new_data['PHASE'] = l['Phase']

        #key = (subj, l['VISCODE'].strip(), l['VISCODE2'].strip())
        key = subj
        subj_dict[key].append(new_data)
    return subj_dict

def parseFuzzyDate(date_string):
    '''
    Assumes --/--/-- format, month/day/year
    '''
    date_string = date_string.replace('-4','').strip()
    if not date_string:
        return None
    date_string = re.sub('--','/1/', date_string)
    date_string = re.sub('[/]+','/', date_string)
    date_string = re.sub('^/','', date_string)
    parts = date_string.split('/')
    if len(parts) != 3:
        raise Exception("Invalid parts: %s" % date_string)
    for i,p in enumerate(parts):
        p = p.replace('-','').strip()
        if p == '':
            p = '1'
        parts[i] = p
    year = parts[-1]
    if len(year) == 4:
        date = datetime.strptime('/'.join(parts), '%m/%d/%Y')
    elif len(year) == 2:
        date = datetime.strptime('/'.join(parts), '%m/%d/%y')
    else:
        raise Exception("Invalid year: %s" % date_string)
    return date

def parseRecMeds(recmed_file):
    reader = csv.DictReader(open(recmed_file, 'rU'))
    subj_dict = defaultdict(list)
    for i,l in enumerate(reader):
        subj = int(l['RID'])
        new_data = {}
        new_data['VISCODE'] = l['VISCODE']
        new_data['VISCODE2'] = l['VISCODE2']
        new_data['INPUTDATE'] = datetime.strptime(l['USERDATE'],'%Y-%m-%d')
        new_data['PHASE'] = l['Phase']
        #new_data['RECNO'] = l['RECNO']
        new_data['EXAMDATE'] = parseFuzzyDate(l['EXAMDATE'])
        new_data['CMMED'] = l['CMMED'].replace('-4','').lower() # Concurrent Non-study Medication (including key background medication), <display>select id,medname from adnigo.medlist where lower(medname) like lower('%$S1.qe%') order by lower(medname)</display><process target_id_field="CMMEDID">select id,medname from adnigo.medlist where lower(medname)=lower('$S1.qe')</process>
        new_data['CMDOSE'] = l['CMDOSE'] # Dose
        new_data['CMFREQNC'] = l['CMFREQNC'].replace('-4','').lower() # Frequency, <display>select freqcode,freqcode||' - '||frequency from adnigo.freqlist where lower(freqcode||' - '||frequency) like lower('%$S1.qe%') order by lower(freqcode)</display><process target_id_field="CMFREQID">select freqcode,freqcode||' - '||frequency from adnigo.freqlist where freqcode||' - '||frequency='$S1.qe'</process>
        new_data['CMROUTE'] = re.sub('[^a-z]+', '', l['CMROUTE'].replace('-4','').lower()) # Route, <display>select routecode,routecode||' - '||route from adnigo.routelist where lower(routecode||' - '||route) like lower('%$S1.qe%') order by lower(routecode)</display><process target_id_field="CMROUTEID">select routecode,routecode||' - '||route from adnigo.routelist where routecode||' - '||route='$S1.qe'</process>
        new_data['CMREASON'] = l['CMREASON'] # Reason prescribed, 1=Adverse Event;2=Therapeutic Use;3=Prophylaxis/Non-therapeutic Use
        new_data['CMBGN'] = parseFuzzyDate(l['CMBGN']) # Date Began
        new_data['CMEND'] = parseFuzzyDate(l['CMEND']) # date ended
        new_data['CMCONT'] = True if l['CMCONT'] == '1' else False # Continuing, 1=True 0=False
        new_data['CMCOMM'] = l['CMCOMM'].replace('-4','').lower() # comments

        if not new_data['CMMED']:
            continue

        #key = (subj, l['VISCODE'].strip(), l['VISCODE2'].strip())
        key = subj
        subj_dict[key].append(new_data)
    return subj_dict

def aggregateBackAndConcurrent(recmed_dict, backmed_dict):
    aggregate = defaultdict(dict)

    for subj,v in recmed_dict.items():
        for record in v:
            visit_key = (record['VISCODE'], record['VISCODE2'], record['PHASE'])
            record_data  = {'DOSE': record['CMDOSE'],
                            'FREQ': record['CMFREQNC'],
                            'END': record['CMEND'],
                            'REASON': record['CMREASON'], # CONVERT
                            'CONT': record['CMCONT'],
                            'BGN': record['CMBGN'],
                            'MED': record['CMMED'],
                            'INPUTDATE': record['INPUTDATE'],
                            'ROUTE': record['CMROUTE'],
                            'COMMENT': record['CMCOMM']}
            if visit_key in aggregate[subj]:
                aggregate[subj][visit_key].append(record_data)
            else:
                aggregate[subj][visit_key] = [record_data]
    for subj,v in backmed_dict.items():
        for record in v:
            visit_key = (record['VISCODE'], record['VISCODE2'], record['PHASE'])
            meds = [KEYMED_TRANSLATE[_] for _ in record['KEYMED']]
            record_data  = {'DOSE': None,
                            'FREQ': None,
                            'END': False,
                            'REASON': None,
                            'CONT': True,
                            'BGN': None,
                            'MED': meds,
                            'INPUTDATE': record['INPUTDATE'],
                            'ROUTE': None,
                            'COMMENT': None}
            if visit_key in aggregate[subj]:
                aggregate[subj][visit_key].append(record_data)
            else:
                aggregate[subj][visit_key] = [record_data]


    for k,v in aggregate.items():
        print "Subj: %s" % (k)
        for visit, records in v.items():
            print visit
            for r in sorted(records, key=lambda x: x['INPUTDATE']):
                print "\t%s: %s" % (r['INPUTDATE'],r['MED'])
        print "\n"


if __name__ == "__main__":
    backmed_file = "../docs/BACKMEDS.csv" # key background medications
    recmed_file = "../docs/RECCMEDS.csv" # concurrent medications
    
    recmed_dict = parseRecMeds(recmed_file)
    backmed_dict = parseBackMeds(backmed_file)

    aggregateBackAndConcurrent(recmed_dict, backmed_dict)


'''
a visit:

VISCODE2: sc
VISCODE: v01
PHASE: ADNI2
EXAMDATE: None (optional - omit for now, fill in later)


each visit will have multiple recorded events...

for concurrent:

CMDOSE: 3.25
CMFREQNC: 3
CMEND: None
CMREASON: 2
CMCONT: True
CMBGN: 2008-01-01 00:00:00
CMMED: zopiclone
input_date: 2013-08-30 00:00:00
CMROUTE: 
CMCOMM: subject has been told she cannot take zopiclone within 2 days of visits and has agreed to comply.


for back:

KEYMED: [0]
input_date: 2013-08-29 00:00:00



'''





