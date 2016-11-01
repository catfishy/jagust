'''
Finding time specific PET-MRI stats

Looks for PET visits that also included an MRI within a given n-day window
'''


import csv
from collections import defaultdict
from datetime import datetime, timedelta
from tabulate import tabulate

def loadScanMeta(pet_meta, mri_meta):
    pet_reader = csv.DictReader(open(pet_meta, 'rU'))
    mri_reader = csv.DictReader(open(mri_meta, 'rU'))

    pet_dict = defaultdict(list)
    for l in pet_reader:
        subj = int(l['Subject'].split('_')[-1])
        date = datetime.strptime(l['Scan Date'],'%m/%d/%y')
        pet_dict[subj].append(date)

    mri_dict = defaultdict(list)
    for l in mri_reader:
        # only look at processed
        if l['Orig/Proc'] != 'Processed':
            continue

        subj = int(l['SubjectID'].split('_')[-1])
        date = datetime.strptime(l['ScanDate'],'%Y-%m-%d')
        mri_dict[subj].append(date)

    return (pet_dict, mri_dict)

def findTimeSpecificPairs(pet_dict, mri_dict, window=90):
    '''
    print "%s pet scan subjects" % (len(pet_dict.keys()))
    print "%s mri scan subjects" % (len(mri_dict.keys()))
    print "%s day window" % window
    '''
    by_datapoints = defaultdict(int)

    by_timepoint_pet = defaultdict(int)
    by_timepoint_mri = defaultdict(int)

    # pet counts
    pet_counts = defaultdict(int)
    for subj in pet_dict.keys():
        pet_counts[len(pet_dict[subj])] += 1
    #print pet_counts

    for subj in pet_dict.keys():
        sorted_pets = sorted(list(set(pet_dict[subj])))
        '''
        if len(sorted_pets) >= 2:
            time_between = [sorted_pets[i+1]-sorted_pets[i] for i,p in enumerate(sorted_pets[1:])]
            print time_between
        '''
        sorted_mris = sorted(list(set(mri_dict[subj])))
        matches = 0
        if len(sorted_pets) == 0:
            continue
        for tp, p in enumerate(sorted_pets):
            by_timepoint_pet[tp] += 1
        for tp, p in enumerate(sorted_pets):
            matched = False
            for i,m in enumerate(sorted_mris):
                if abs(p-m) <= timedelta(days=window):
                    matches += 1
                    matched = True
                    sorted_mris = sorted_mris[i+1:]
                    by_timepoint_mri[tp] += 1
                    break
            
            if not matched:
                break
            
        by_datapoints[matches] += 1

    by_timepoint_pet = dict(by_timepoint_pet)
    by_timepoint_mri = dict(by_timepoint_mri)

    by_datapoints = dict(by_datapoints)
    totals = sum(by_datapoints.values())
    bl_totals = sum([by_datapoints[1],by_datapoints[2],by_datapoints[3]])
    for k in by_datapoints.keys():
        if k == 0:
            by_datapoints[k] = "%s" % by_datapoints[k]
        else:
            by_datapoints[k] = "%s/%s" % (by_datapoints[k],pet_counts[k])

    return by_datapoints, totals, bl_totals, by_timepoint_pet, by_timepoint_mri


if __name__ == "__main__":
    pet_meta = "/Users/andyhorng/Documents/PET_META_LIST_edited.csv"
    mri_meta = "/Users/andyhorng/Documents/MPRAGEMETA.csv"
    pet_dict, mri_dict = loadScanMeta(pet_meta, mri_meta)

    windows = [10, 30, 60, 90, 120, 180, 400, 3000]
    table = []
    table.append(['Window (Days)', 'Subj w/ 0 mris', 'Subj w/ 1 mri', 'Subj w/ 2 mris', 'Subj w/ 3 mris', 'Total Count', 'Total with BL'])
    count_table = []
    count_table.append(['Window (Days)', 'MRI @ BL', 'PET @ BL', 'MRI @ V2', 'PET @ V2', 'MRI @ V3', 'PET @ V3'])
    for w in windows:
        res, totals, bl_totals, pet, mri = findTimeSpecificPairs(pet_dict, mri_dict, window=w)
        table.append([w, res[0], res[1], res[2], res[3], totals, bl_totals])
        count_table.append([w, mri[0], pet[0], mri[1], pet[1], mri[2], pet[2]])
    print tabulate(table)
    print tabulate(count_table)
    '''
    to_check = [5150,5154,5165,5209,5250,5241,5290,5288,914,5166,5263,5259,5271,5283,5253,5287,5227,5282,4095,5280,5273,5275,4862]

    for subj in to_check:
        print subj
        print list(set(pet_dict[subj]))
        print list(set(mri_dict[subj]))
    '''
    