'''
Compares MRI and PET metadata files
'''

import csv
from collections import defaultdict
from datetime import datetime, timedelta

def compareDates(mri_file, pet_file):
    mri_reader = csv.DictReader(open(mri_file, 'rU'))
    patient_mris = defaultdict(list)
    for row in mri_reader:
        new_date = row['ScanDate']
        new_date = datetime.strptime(new_date, '%Y-%m-%d')
        patient_mris[row['SubjectID']].append(new_date)

    pet_reader = csv.DictReader(open(pet_file, 'rU'))
    patient_pets = defaultdict(list)
    for row in pet_reader:
        new_date = row['Scan Date']
        new_date = datetime.strptime(new_date, '%m/%d/%y')
        patient_pets[row['Subject']].append(new_date)

    datapoints = {}
    for k in patient_pets.keys():
        points = 0
        pets = sorted(list(set(patient_pets[k])))
        mris = sorted(list(set(patient_mris[k])))
        for p in pets:
            # check if there is an mri within 3 months of scan, if there is, pop it
            found = False
            for i,m in enumerate(mris):
                if abs(m-p) < timedelta(days=90):
                    points += 1
                    found = True
                    mris.pop(i)
                    break
            '''
            if not found:
                print "No MRI found: %s, %s" % (p, mris)
            '''
        datapoints[k] = points

    by_count = defaultdict(int)
    for k,v in datapoints.iteritems():
        by_count[v] += 1

    for k,v in by_count.items():
        print "%s, %s" % (k,v)



if __name__ == "__main__":
    pet_meta = "/Users/ahorng/Documents/PET_META_LIST_edited.csv"
    mri_meta = "/Users/ahorng/Documents/MPRAGEMETA.csv"
    compareDates(mri_meta, pet_meta)
