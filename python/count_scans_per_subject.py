from utils import *
from collections import defaultdict


if __name__ == '__main__':
    input_file = "../docs/idaSearch_8_10_2015.csv"
    headers, rows = parseCSV(input_file, delimiter=',')
    search_counter = defaultdict(list)
    for r in rows:
        subj = r['Subject ID']
        search_counter[subj].append(r)

    input_file = "../docs/PET_META_LIST.csv"
    headers, rows = parseCSV(input_file, delimiter=',')
    meta_counter = defaultdict(list)
    for r in rows:
        if "AV45 Coreg, Avg, Std Img and Vox Siz, Uniform Resolution" not in r['Sequence']:
            continue
        subj = r['Subject']
        meta_counter[subj].append(r)

    search_counts = {0:0, 1:0, 2:0, 3:0, 4:0}
    meta_counts = {0:0, 1:0, 2:0, 3:0, 4:0}
    for k,v in search_counter.iteritems():
        search_counts[len(v)] += 1
    for k,v in meta_counter.iteritems():
        meta_counts[len(v)] += 1

    print search_counts
    print meta_counts