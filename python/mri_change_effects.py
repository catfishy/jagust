from utils import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

size_keys = ['CTX_LH_CAUDALMIDDLEFRONTAL_SIZE',
            'CTX_LH_LATERALORBITOFRONTAL_SIZE',
            'CTX_LH_MEDIALORBITOFRONTAL_SIZE',
            'CTX_LH_PARSOPERCULARIS_SIZE',
            'CTX_LH_PARSORBITALIS_SIZE',
            'CTX_LH_PARSTRIANGULARIS_SIZE',
            'CTX_LH_ROSTRALMIDDLEFRONTAL_SIZE',
            'CTX_LH_SUPERIORFRONTAL_SIZE',
            'CTX_LH_FRONTALPOLE_SIZE',
            'CTX_RH_CAUDALMIDDLEFRONTAL_SIZE',
            'CTX_RH_LATERALORBITOFRONTAL_SIZE',
            'CTX_RH_MEDIALORBITOFRONTAL_SIZE',
            'CTX_RH_PARSOPERCULARIS_SIZE',
            'CTX_RH_PARSORBITALIS_SIZE',
            'CTX_RH_PARSTRIANGULARIS_SIZE',
            'CTX_RH_ROSTRALMIDDLEFRONTAL_SIZE',
            'CTX_RH_SUPERIORFRONTAL_SIZE',
            'CTX_RH_FRONTALPOLE_SIZE',
            'CTX_LH_CAUDALANTERIORCINGULATE_SIZE',
            'CTX_LH_ISTHMUSCINGULATE_SIZE',
            'CTX_LH_POSTERIORCINGULATE_SIZE',
            'CTX_LH_ROSTRALANTERIORCINGULATE_SIZE',
            'CTX_RH_CAUDALANTERIORCINGULATE_SIZE',
            'CTX_RH_ISTHMUSCINGULATE_SIZE',
            'CTX_RH_POSTERIORCINGULATE_SIZE',
            'CTX_RH_ROSTRALANTERIORCINGULATE_SIZE',
            'CTX_LH_INFERIORPARIETAL_SIZE',
            'CTX_LH_PRECUNEUS_SIZE',
            'CTX_LH_SUPERIORPARIETAL_SIZE',
            'CTX_LH_SUPRAMARGINAL_SIZE',
            'CTX_RH_INFERIORPARIETAL_SIZE',
            'CTX_RH_PRECUNEUS_SIZE',
            'CTX_RH_SUPERIORPARIETAL_SIZE',
            'CTX_RH_SUPRAMARGINAL_SIZE',
            'CTX_LH_MIDDLETEMPORAL_SIZE',
            'CTX_LH_SUPERIORTEMPORAL_SIZE',
            'CTX_RH_MIDDLETEMPORAL_SIZE',
            'CTX_RH_SUPERIORTEMPORAL_SIZE',
            ]


if __name__ == "__main__":
    nontp_file = '../output/UCBERKELEYAV45_09_25_15_extra_nontp.csv'
    tp_file = '../output/UCBERKELEYAV45_09_25_15_extra.csv'
    registry_file = "../docs/registry_clean.csv"
    master_file = '../FDG_AV45_COGdata_09_25_15.csv'

    registry = importRegistry(registry_file)
    tp_data = importAV45(tp_file, av45_nontp_file=None ,registry=registry)
    nontp_data = importAV45(nontp_file, av45_nontp_file=None, registry=registry)
    master_data = importMaster(master_file) 

    size_changes = {}
    uptake_diffs = {}
    hemiwm_uptake_diffs = {}
    pvc_diffs = {}
    pvc_hemiwm_diffs = {}
    tp_points = []
    tp_hemiwm_points = []
    pvc_points = []
    pvc_hemiwm_points = []
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

        init_diag = master_row['Init_Diagnosis'].strip()

        colors = {'EMCI': 'r',
                  'LMCI': 'm',
                  'N' : 'g',
                  'AD': 'b',
                  'SMC': 'c'}

        color = colors[init_diag]

        tp_first = tp_rows[0]
        tp_second = tp_rows[1]
        nontp_first = rows[0]
        nontp_second = rows[1]

        # check exam dates
        if nontp_first['EXAMDATE'] != tp_first['EXAMDATE']:
            print "First EXAMdates don't match: %s" % rid
            continue
        if nontp_second['EXAMDATE'] != tp_second['EXAMDATE']:
            print "Second EXAMdates don't match: %s" % rid
            continue


        tp_change = float(tp_second['SUMMARYSUVR_WHOLECEREBNORM']) - float(tp_first['SUMMARYSUVR_WHOLECEREBNORM'])
        nontp_change = float(nontp_second['SUMMARYSUVR_WHOLECEREBNORM']) - float(nontp_first['SUMMARYSUVR_WHOLECEREBNORM'])
        tp_hemiwm_change = float(tp_second['ERODED_SUBCORTICALWM']) - float(tp_first['ERODED_SUBCORTICALWM'])
        nontp_hemiwm_change = float(nontp_second['ERODED_SUBCORTICALWM']) - float(nontp_first['ERODED_SUBCORTICALWM'])
       
        # calculate size changes
        size_change = np.mean([(int(tp_second[_]) - int(tp_first[_])) for _ in size_keys])
        size_changes[rid] = size_change
        
        # calculate uptake difference between tp and nontp
        if float(tp_first['SUMMARYSUVR_WHOLECEREBNORM']) >= 1.11:
            change_diff = tp_change - nontp_change
            hemiwm_change_diff = tp_hemiwm_change - nontp_hemiwm_change
            uptake_diffs[rid] = change_diff
            hemiwm_uptake_diffs[rid] = hemiwm_change_diff
            tp_points.append((size_changes[rid], uptake_diffs[rid], color))
            tp_hemiwm_points.append((size_changes[rid], hemiwm_uptake_diffs[rid], color))


        # calculate uptake difference between pvc and nonpvc
        try:
            if float(master_row['AV45_wcereb']) >= 1.11:
                first_pvc = float(master_row['AV45_PVC_agghigh_CorticalSummary/Wholecereb_BL'])
                second_pvc = float(master_row['AV45_PVC_agghigh_CorticalSummary/Wholecereb_Scan2'])
                first_hemiwm_pvc = float(master_row['AV45_PVC_agghigh_HemiWM/Wholecereb_BL'])
                second_hemiwm_pvc = float(master_row['AV45_PVC_agghigh_HemiWM/Wholecereb_Scan2'])
                pvc_change = second_pvc - first_pvc
                pvc_change_diff = pvc_change - nontp_change
                pvc_hemiwm_change = second_hemiwm_pvc - first_hemiwm_pvc
                pvc_hemiwm_change_diff = pvc_hemiwm_change - tp_hemiwm_change
                pvc_diffs[rid] = pvc_change_diff
                pvc_hemiwm_diffs[rid] = pvc_hemiwm_change_diff
                pvc_points.append((size_changes[rid], pvc_diffs[rid], color))
                pvc_hemiwm_points.append((size_changes[rid], pvc_hemiwm_diffs[rid], color))


                if pvc_hemiwm_diffs[rid] < -1.0:
                    print rid
                    print "size diff: %s" % size_changes[rid]
                    print "first pvc summary: %s" % first_pvc
                    print "second pvc summary: %s" % second_pvc
                    print "first pvc hemiwm: %s" % first_hemiwm_pvc
                    print "second pvc hemiwm: %s" % second_hemiwm_pvc
                    print 'first nonpvc hemiwm: %s' % tp_first['ERODED_SUBCORTICALWM']
                    print "second nonpvc hemiwm: %s" % tp_second['ERODED_SUBCORTICALWM']
                    print "size diff: %s, pvchemiwm diff: %s" % (size_changes[rid], pvc_hemiwm_diffs[rid])
                    print "size diff: %s, pvcsumm diff: %s" % (size_changes[rid], pvc_diffs[rid])

        except Exception as e:
            pass

    
    # plot cortical summary
    '''
    plt.figure(1)
    x = [_[0] for _ in tp_points]
    y = [_[1] for _ in tp_points]
    c = [_[2] for _ in tp_points]
    rho, p = spearmanr(x,y)
    plt.scatter(x,y,c=c)
    plt.xlabel('Mean Size Difference')
    plt.ylabel('Cortical Summary Difference (TP vs Non-TP)')
    props = dict(boxstyle='round,pad=0.8', facecolor='wheat', alpha=0.7)
    text = "$\\rho=%s$\n$p=%s$" % (rho, p)
    plt.annotate(text, xy=(0.05,0.85), xycoords='axes fraction', bbox=props, fontsize=13)

    plt.figure(2)
    x = [_[0] for _ in pvc_points]
    y = [_[1] for _ in pvc_points]
    c = [_[2] for _ in pvc_points]
    rho, p = spearmanr(x,y)
    plt.scatter(x,y,c=c)
    plt.xlabel('Mean Size Difference')
    plt.ylabel('Cortical Summary Difference (PVC vs Non-PVC)')
    props = dict(boxstyle='round,pad=0.8', facecolor='wheat', alpha=0.7)
    text = "$\\rho=%s$\n$p=%s$" % (rho, p)
    plt.annotate(text, xy=(0.05,0.85), xycoords='axes fraction', bbox=props, fontsize=13)
    plt.show()
    '''

    
    # plot hemiwm
    plt.figure(1)
    x = [_[0] for _ in tp_hemiwm_points]
    y = [_[1] for _ in tp_hemiwm_points]
    c = [_[2] for _ in tp_hemiwm_points]
    rho, p = spearmanr(x,y)
    plt.scatter(x,y,c=c)
    plt.xlabel('Mean Size Difference')
    plt.ylabel('Hemi WM Change Difference (TP vs Non-TP)')
    props = dict(boxstyle='round,pad=0.8', facecolor='wheat', alpha=0.7)
    text = "$\\rho=%s$\n$p=%s$" % (rho, p)
    plt.annotate(text, xy=(0.05,0.85), xycoords='axes fraction', bbox=props, fontsize=13)

    plt.figure(2)
    x = [_[0] for _ in pvc_hemiwm_points]
    y = [_[1] for _ in pvc_hemiwm_points]
    c = [_[2] for _ in pvc_hemiwm_points]
    rho, p = spearmanr(x,y)
    plt.scatter(x,y,c=c)
    plt.xlabel('Mean Size Difference')
    plt.ylabel('Hemi WM Change Difference (PVC vs Non-PVC)')
    props = dict(boxstyle='round,pad=0.8', facecolor='wheat', alpha=0.7)
    text = "$\\rho=%s$\n$p=%s$" % (rho, p)
    plt.annotate(text, xy=(0.05,0.85), xycoords='axes fraction', bbox=props, fontsize=13)

    plt.show()
    
