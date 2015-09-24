import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from collections import defaultdict
from scipy.stats import norm, mannwhitneyu, linregress
from sklearn.mixture import GMM, VBGMM, DPGMM
from utils import *
import itertools
from scipy import linalg
import matplotlib as mpl


def graph_tp_vs_nontp():
    registry = importRegistry(reg_file)
    tp_data = importAV45(av45_tp_file, registry=registry)
    nontp_data = importAV45(av45_nontp_file, registry=registry)

    # find subjects with more than 2 timepoints
    valid_tp = {k:v for k,v in tp_data.iteritems() if len(v) >= 2}
    valid_nontp = {k:v for k,v in nontp_data.iteritems() if len(v) >= 2}

    overlap = list(set(valid_tp.keys()) & set(valid_nontp.keys()))

    points = []
    for k in overlap:
        tp_rows = sorted(valid_tp[k], key=lambda x: x['EXAMDATE'])
        nontp_rows = sorted(valid_nontp[k], key=lambda x: x['EXAMDATE'])

        bl_tp = float(tp_rows[0]['SUMMARYSUVR_WHOLECEREBNORM'])
        bl_nontp = float(nontp_rows[0]['SUMMARYSUVR_WHOLECEREBNORM'])

        change_tp = float(tp_rows[1]['SUMMARYSUVR_WHOLECEREBNORM']) - bl_tp
        change_nontp = float(nontp_rows[1]['SUMMARYSUVR_WHOLECEREBNORM']) - bl_nontp

        '''
        points.append((bl_tp, change_tp, 'g'))
        points.append((bl_nontp, change_nontp, 'r'))
        '''
        difference = change_tp-change_nontp
        if abs(difference) > .10:
            print k
        points.append((bl_tp, difference, 'b'))

    x = [_[0] for _ in points]
    y = [_[1] for _ in points]
    c = [_[2] for _ in points]

    plt.plot()
    plt.scatter(x,y,c=c)
    plt.show()

def graph_pvc_vs_nonpvc():
    data = importMaster(big_csv_file)
    # find valid tp one and two
    points = []
    for k,v in data.iteritems():
        bl_spec = 0
        scan2_spec = 0
        try:
            bl_spec = int(v['AV45_TP_SPECIFIC_BL'])
            scan2_spec = int(v['AV45_TP_SPECIFIC_SCAN2'])
        except Exception as e:
            pass

        if bl_spec != 1 or scan2_spec != 1:
            continue
        if str(v['AV45_CorticalSummary/Wholecereb_PVC_BL']).strip() == '' or str(v['AV45_CorticalSummary/Wholecereb_PVC_Scan2']).strip() == '':
            continue
        pvc_difference = float(v['AV45_CorticalSummary/Wholecereb_PVC_Scan2']) - float(v['AV45_CorticalSummary/Wholecereb_PVC_BL'])
        pvc_bl = float(v['AV45_CorticalSummary/Wholecereb_PVC_BL'])
        nonpvc_difference = float(v['AV45_2_wcereb']) - float(v['AV45_wcereb'])
        nonpvc_bl = float(v['AV45_wcereb'])

        difference = pvc_difference-nonpvc_difference
        if abs(difference) > 0.5:
            print k
        '''
        points.append((pvc_bl, difference, 'g'))
        '''
        points.append((pvc_bl, pvc_difference, 'r'))
        points.append((nonpvc_bl, nonpvc_difference, 'b'))
        

    print "points: %s" % len(points)
    x = [_[0] for _ in points]
    y = [_[1] for _ in points]
    c = [_[2] for _ in points]

    plt.plot()
    plt.scatter(x,y,c=c)
    plt.show()



if __name__ == "__main__":
    reg_file = "../docs/registry_clean.csv"
    av45_tp_file = "../output/UCBERKELEYAV45_09_21_15_extra.csv"
    av45_nontp_file = "../output/UCBERKELEYAV45_09_21_15_extra_nontp.csv"

    big_csv_file = "../FDG_AV45_COGdata_09_21_15.csv"


    graph_pvc_vs_nonpvc()
