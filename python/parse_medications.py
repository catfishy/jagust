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

def parseBackMeds(backmed_file):
    reader = csv.DictReader(open(backmed_file, 'rU'))
    subj_dict = defaultdict(list)
    for l in reader:
        subj = int(l['RID'])
        new_data = {}
        new_data['VISCODE'] = l['VISCODE']
        new_data['VISCODE2'] = ['VISCODE2']
        new_data['input_date'] = datetime.strptime(l['USERDATE'],'%Y-%m-%d')
        new_data['KEYMED'] = map(int, l['KEYMED'].split(':')) # At this visit, is participant on any of the following medication?, 1=Aricept;2=Cognex;3=Exelon;4=Namenda;5=Razadyne;6=Anti-depressant medication;7=Other behavioral medication;0=None of the above
        new_data['PHASE'] = l['Phase']
        subj_dict[subj].append(new_data)
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
        new_data['VISCODE2'] = ['VISCODE2']
        new_data['input_date'] = datetime.strptime(l['USERDATE'],'%Y-%m-%d')
        new_data['PHASE'] = l['Phase']
        #ew_data['RECNO'] = l['RECNO']
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

        subj_dict[subj].append(new_data)
    return subj_dict

if __name__ == "__main__":
    backmed_file = "/Users/ahorng/Documents/jagust/docs/BACKMEDS.csv" # key background medications
    recmed_file = "/Users/ahorng/Documents/jagust/docs/RECCMEDS.csv" # concurrent medications
    
    recmed_dict = parseRecMeds(recmed_file)
    backmed_dict = parseBackMeds(backmed_file)


