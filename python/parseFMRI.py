import pandas as pd

from utils import *

fmri_file = '../FMRI_idaSearch_9_10_2015.csv'
master_file = '../FDG_AV45_COGdata_09_10_15.csv'
output_file = '../FMRI_RestingState_meta.csv'

fmri_frame = pd.read_csv(fmri_file, low_memory=False)
master_frame = pd.read_csv(master_file, low_memory=False)

by_subj = fmri_frame.groupby('Subject ID')

new_headers = ['Subject ID', 'Diagnosis @ Enrollment', 'Scans', 'Sex', 'Research Group', 'Imaging Protocol', 'Study Dates', 'Ages', 'Descriptions']
new_rows = []
for subject, rows in by_subj:
	print subject
	rid = int(subject.split('_')[-1])
	masterrow = master_frame.loc[master_frame['RID'] == rid]
	if masterrow.shape[0] == 0:
		print "Could not find Diag for %s" % subject
		initdiag = ''
	else:
		initdiag = masterrow['Init_Diagnosis'].iloc[0]
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