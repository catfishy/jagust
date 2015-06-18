import csv
from collections import defaultdict


def main(file_path):
    reader = csv.reader(open(file_path,'r'))
    count = 0
    labels = None
    by_subj = defaultdict(list)
    for count, r in enumerate(reader):
        if count == 0:
            labels = r
            continue
        data = dict(zip(labels, r))
        vis_code = data['VISCODE2']
        rid = data['RID']
        by_subj[rid].append(vis_code)
    
    old_counts = defaultdict(int)
    for k,v in by_subj.items():
        for i,m in enumerate(v):
            old_counts[i+1] += 1

    for k,v in old_counts.items():
        print k
        print v


if __name__ == "__main__":
    file_path = "/Users/ahorng/Documents/UCBERKELEYAV45_01_26_15.csv"
    main(file_path)