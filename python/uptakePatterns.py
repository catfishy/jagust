from utils import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from collections import defaultdict

UPTAKE_KEYS =  ['LEFT-PUTAMEN',
                'RIGHT-PUTAMEN',
                'LEFT-CAUDATE',
                'RIGHT-CAUDATE',
                'LEFT-PALLIDUM',
                'RIGHT-PALLIDUM',
                'CTX_LH_CAUDALMIDDLEFRONTAL',
                'CTX_LH_LATERALORBITOFRONTAL',
                'CTX_LH_MEDIALORBITOFRONTAL',
                'CTX_LH_PARSOPERCULARIS',
                'CTX_LH_PARSORBITALIS',
                'CTX_LH_PARSTRIANGULARIS',
                'CTX_LH_ROSTRALMIDDLEFRONTAL',
                'CTX_LH_SUPERIORFRONTAL',
                'CTX_LH_FRONTALPOLE',
                'CTX_RH_CAUDALMIDDLEFRONTAL',
                'CTX_RH_LATERALORBITOFRONTAL',
                'CTX_RH_MEDIALORBITOFRONTAL',
                'CTX_RH_PARSOPERCULARIS',
                'CTX_RH_PARSORBITALIS',
                'CTX_RH_PARSTRIANGULARIS',
                'CTX_RH_ROSTRALMIDDLEFRONTAL',
                'CTX_RH_SUPERIORFRONTAL',
                'CTX_RH_FRONTALPOLE',
                'CTX_LH_CAUDALANTERIORCINGULATE',
                'CTX_LH_ISTHMUSCINGULATE',
                'CTX_LH_POSTERIORCINGULATE',
                'CTX_LH_ROSTRALANTERIORCINGULATE',
                'CTX_RH_CAUDALANTERIORCINGULATE',
                'CTX_RH_ISTHMUSCINGULATE',
                'CTX_RH_POSTERIORCINGULATE',
                'CTX_RH_ROSTRALANTERIORCINGULATE',
                'CTX_LH_INFERIORPARIETAL',
                'CTX_LH_PRECUNEUS',
                'CTX_LH_SUPERIORPARIETAL',
                'CTX_LH_SUPRAMARGINAL',
                'CTX_RH_INFERIORPARIETAL',
                'CTX_RH_PRECUNEUS',
                'CTX_RH_SUPERIORPARIETAL',
                'CTX_RH_SUPRAMARGINAL',
                'CTX_LH_MIDDLETEMPORAL',
                'CTX_LH_SUPERIORTEMPORAL',
                'CTX_RH_MIDDLETEMPORAL',
                'CTX_RH_SUPERIORTEMPORAL']


def deriveUptakePattern(row):
    wcereb = float(row['WHOLECEREBELLUM']) 
    comp_ref = float(row['COMPOSITE_REF'])
    # normalize by either comp_ref or wcereb
    pattern = np.array([float(row[_])/comp_ref for _ in UPTAKE_KEYS])
    pattern = pattern/float(sum(pattern))
    return pattern


if __name__ == "__main__":
    nontp_file = '../output/UCBERKELEYAV45_09_25_15_extra_nontp.csv'
    tp_file = '../output/UCBERKELEYAV45_09_25_15_extra.csv'
    registry_file = "../docs/registry_clean.csv"
    master_file = '../FDG_AV45_COGdata_09_25_15.csv'

    registry = importRegistry(registry_file)
    tp_data = importAV45(tp_file, av45_nontp_file=None ,registry=registry)
    nontp_data = importAV45(nontp_file, av45_nontp_file=None, registry=registry)
    master_data = importMaster(master_file) 

    # TODO: USE PVC CLUSTER VALUES AS FEATURES
    # TODO: USE APARC+ASEG VALUES AS FEATURES (look into analyzeCentroid)


    # only using nontp data
    uptake_patterns = {}
    diagnoses = {}
    summary = {}
    for rid, rows in nontp_data.iteritems():
        master_row = master_data.get(rid)
        if not master_row:
            continue
        diagnoses[rid] = master_row['Init_Diagnosis'].strip()
        # filter out non-normals
        if diagnoses[rid] != 'N':
            continue
        for row in rows:
            # filter out positives
            if float(row['SUMMARYSUVR_WHOLECEREBNORM']) > 1.11:
                continue
            uptake_vector = deriveUptakePattern(row)
            date = row['EXAMDATE']
            uptake_patterns[(rid, date)] = uptake_vector
            summary[(rid, date)] = float(row['SUMMARYSUVR_WHOLECEREBNORM'])
    keys = []
    data = []
    for k,v in uptake_patterns.iteritems():
        keys.append(k)
        data.append(v)
    data = np.array(data)
    print data.shape

    '''
    ks=range(1,30)
    gaps, valid_ks = gap(data, nrefs=100, ks=ks, use_pca=False)
    print "GAPS"
    print gaps
    print sorted(valid_ks, key=lambda x: x[1], reverse=True)
    fig = plt.figure(1)
    plt.plot(ks, gaps)
    plt.show()
    '''

    clusters=9
    model = KMeans(n_clusters=clusters, n_jobs=1, copy_x=True, verbose=True)
    labels = model.fit_predict(data)
    centers = model.cluster_centers_
    by_label = defaultdict(list)
    for i,l in enumerate(labels):
        by_label[l].append(keys[i])
    by_label = dict(by_label)

    for c in centers:
        print c

    # TODO: MAP CLUSTER CENTERS TO MODEL aparc+aseg, and visualize

