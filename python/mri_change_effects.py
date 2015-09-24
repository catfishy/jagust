from utils import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


if __name__ == "__main__":
    nontp_file = '../output/UCBERKELEYAV45_09_21_15_extra_nontp.csv'
    tp_file = '../output/UCBERKELEYAV45_09_21_15_extra.csv'
    registry_file = "../docs/registry_clean.csv"
    master_file = '../FDG_AV45_COGdata_09_23_15.csv'

    registry = importRegistry(registry_file)
    tp_data = importAV45(tp_file, av45_nontp_file=None ,registry=registry)
    nontp_data = importAV45(nontp_file, av45_nontp_file=None, registry=registry)
    master_data = importMaster(master_file) 

    size_changes = {}
    uptake_diffs = {}
    pvc_diffs = {}
    tp_points = []
    pvc_points = []
    for rid, rows in nontp_data.iteritems():
        if len(rows) < 2:
            continue
        rows = sorted(rows, key=lambda x: x['EXAMDATE'])

        # get tp rows
        tp_rows = tp_data.get(rid,[])
        if len(tp_rows) < 2:
            continue
        tp_rows = sorted(tp_rows, key=lambda x: x['EXAMDATE'])

        # get master row
        master_row = master_data.get(rid)
        if not master_row:
            continue

        tp_first = tp_rows[0]
        tp_second = tp_rows[1]
        nontp_first = rows[0]
        nontp_second = rows[1]
        size_keys = sorted([_ for _ in tp_rows[0].keys() if '_SIZE' in _])
        
        # calculate size changes
        size_change = [abs(int(tp_second[_]) - int(tp_first[_])) for _ in size_keys]
        size_changes[rid] = np.mean(size_change)
        
        # calculate uptake difference between tp and nontp
        if float(tp_first['SUMMARYSUVR_WHOLECEREBNORM']) >= 1.11:
            tp_change = float(tp_second['SUMMARYSUVR_WHOLECEREBNORM']) - float(tp_first['SUMMARYSUVR_WHOLECEREBNORM'])
            nontp_change = float(nontp_second['SUMMARYSUVR_WHOLECEREBNORM']) - float(nontp_first['SUMMARYSUVR_WHOLECEREBNORM'])
            change_diff = abs(tp_change - nontp_change)
            uptake_diffs[rid] = change_diff
            tp_points.append((size_changes[rid], uptake_diffs[rid]))

        # calculate uptake difference between pvc and nonpvc
        try:
            if float(master_row['AV45_wcereb']) >= 1.11:
                first_pvc = float(master_row['AV45_PVC_agghigh_CorticalSummary/Wholecereb_BL'])
                second_pvc = float(master_row['AV45_PVC_agghigh_CorticalSummary/Wholecereb_Scan2'])
                pvc_change = second_pvc - first_pvc
                pvc_change_diff = abs(pvc_change - nontp_change)
                pvc_diffs[rid] = pvc_change_diff
                pvc_points.append((size_changes[rid], pvc_diffs[rid]))
        except Exception as e:
            pass

    
    # plot it
    plt.figure(1)
    x = [_[0] for _ in tp_points]
    y = [_[1] for _ in tp_points]
    rho, p = spearmanr(x,y)
    plt.scatter(x,y,c='r')
    plt.xlabel('Mean Size Difference')
    plt.ylabel('Cortical Summary Difference (TP vs Non-TP)')
    props = dict(boxstyle='round,pad=0.8', facecolor='wheat', alpha=0.7)
    text = "$\\rho=%s$\n$p=%s$" % (rho, p)
    plt.annotate(text, xy=(0.05,0.85), xycoords='axes fraction', bbox=props, fontsize=13)

    plt.figure(2)
    x = [_[0] for _ in pvc_points]
    y = [_[1] for _ in pvc_points]
    rho, p = spearmanr(x,y)
    plt.scatter(x,y,c='b')
    plt.xlabel('Mean Size Difference')
    plt.ylabel('Cortical Summary Difference (PVC vs Non-PVC)')
    props = dict(boxstyle='round,pad=0.8', facecolor='wheat', alpha=0.7)
    text = "$\\rho=%s$\n$p=%s$" % (rho, p)
    plt.annotate(text, xy=(0.05,0.85), xycoords='axes fraction', bbox=props, fontsize=13)

    plt.show()
