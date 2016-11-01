from utils import *

mri_meta_file = "../docs/MPRAGEMETA.csv"


def findBothStrengthMRIS():
    ones = importMRI(mri_meta_file, magstrength_filter='1.5')
    threes = importMRI(mri_meta_file, magstrength_filter='3.0')

    all_subj = list(set(ones.keys()) | set(threes.keys()))

    overlaps = []
    for s in all_subj:
        onescans = {_['vc']: _ for _ in ones.get(s,[])}
        threescans = {_['vc']: _ for _ in threes.get(s,[])}
        for k in list(set(onescans.keys()) & set(threescans.keys())):
            overlaps.append({'RID': s, 'Visit': k, 'Date_1.5': onescans[k]['EXAMDATE'], 'Date_3.0': threescans[k]['EXAMDATE']})
    dumpCSV('mri_1.5_and_3.0.csv', ['RID', 'Visit', 'Date_1.5', 'Date_3.0'], overlaps)
    return overlaps

findBothStrengthMRIS()