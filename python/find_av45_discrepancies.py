'''
compares our av45 output csv (from aggregate_scan_means) to the 
pet meta list, to determine which scans we are missing/dropping out

* todo: see if the pet scans exist in our local db
'''
import csv
from collections import defaultdict
from datetime import datetime

def find_pet_discrepancies(local_csv, pet_meta_csv):
    meta_reader = csv.DictReader(open(pet_meta_csv, 'rU'))
    meta_pets = defaultdict(list)
    for row in meta_reader:
        subj = int(row['Subject'].split('_')[-1])
        new_date = datetime.strptime(row['Scan Date'], '%m/%d/%y')
        meta_pets[subj].append(new_date)

    local_reader = csv.DictReader(open(local_csv, 'rU'))
    local_pets = defaultdict(list)
    for row in local_reader:
        subj = int(row['RID'])
        new_date = datetime.strptime(row['EXAMDATE'], '%Y-%m-%d')
        local_pets[subj].append(new_date)

    unequal_subj = []
    for k in meta_pets.keys():
        meta_vals = sorted(meta_pets[k])
        local_vals = sorted(local_pets[k])
        if meta_vals != local_vals:
            unequal_subj.append(k)
            print "%s Unequal:" % k
            print "\tMETA: %s" % (meta_vals,)
            print "\tLOCAL: %s" % (local_vals,)

    print unequal_subj



if __name__ == "__main__":
    local_csv = "../UCBERKELEYAV45_06_22_15.csv"
    pet_meta = "../docs/PET_META_LIST_edited.csv"
    find_pet_discrepancies(local_csv, pet_meta)