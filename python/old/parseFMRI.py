import pandas as pd

from utils import *

fmri_file = '../FMRI_idaSearch_9_10_2015.csv'
master_file = '../FDG_AV45_COGdata_09_10_15.csv'
output_file = '../FMRI_RestingState_meta.csv'
diag_file = '../docs/DXSUM_PDXCONV_ADNIALL.csv'
fmri_frame = pd.read_csv(fmri_file, low_memory=False)
master_frame = pd.read_csv(master_file, low_memory=False)

arm_file = "../docs/ARM.csv"
arm = importARM(arm_file)

by_subj = fmri_frame.groupby('Subject ID')

new_headers = ['Subject ID', 'Diagnosis @ Enrollment', 'Scans', 'Sex', 'Research Group', 'Imaging Protocol', 'Study Dates', 'Ages', 'Descriptions']
new_rows = []
for subject, rows in by_subj:
    print subject
    rid = int(subject.split('_')[-1])

    masterrow = master_frame.loc[master_frame['RID'] == rid]
    if masterrow.shape[0] == 0:
        initdiag = ''
    else:
        initdiag = masterrow['Init_Diagnosis'].iloc[0]

    if initdiag == '' and rid in arm:
        sorted_arm = sorted(arm[rid], key=lambda x: x['USERDATE'])
        init_diags = list(set([_['STATUS'] for _ in sorted_arm]))
        initdiag = init_diags[-1]

    if initdiag == '':
        print "Could not find diag for %s" % subject

    subj_data = {}
    subj_data['Diagnosis @ Enrollment'] = initdiag

    rows = rows.drop_duplicates(subset=['Study Date']).sort(['Age'])

    subj_data['Subject ID'] = subject
    subj_data['Imaging Protocol'] = rows['Imaging Protocol'].unique()[0]
    subj_data['Sex'] = rows['Sex'].unique()[0]
    subj_data['Research Group'] = rows['Research Group'].unique()[0]
    subj_data['Study Dates'] = rows['Study Date'].tolist()
    subj_data['Ages'] = rows['Age'].tolist()
    subj_data['Scans'] = rows.shape[0]
    subj_data['Descriptions'] = rows['Description'].tolist()
    new_rows.append(subj_data)

dumpCSV(output_file, new_headers, new_rows)