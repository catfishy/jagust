from utils import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.mixture import GMM, VBGMM, DPGMM
from collections import defaultdict
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D

UPTAKE_KEYS =  ['CTX_LH_CAUDALMIDDLEFRONTAL',
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


def deriveCorticalSummaryUptakePattern(input_data, master_data, lut_dict):
    # only using nontp data
    uptake_patterns = {}
    diagnoses = {}
    summary = {}
    for rid, rows in input_data.iteritems():
        master_row = master_data.get(rid)
        if not master_row:
            continue
        diagnoses[rid] = master_row['Init_Diagnosis'].strip()
        '''
        # filter out non-normals
        if diagnoses[rid] != 'N':
            continue
        '''
        for row in rows:
            '''
            # filter out positives
            if float(row['SUMMARYSUVR_WHOLECEREBNORM']) > 1.11:
                continue
            '''
            #wcereb = float(row['WHOLECEREBELLUM']) 
            #comp_ref = float(row['COMPOSITE_REF'])
            # normalize by either comp_ref or wcereb
            try:
                pattern = np.array([float(row[_]) for _ in UPTAKE_KEYS])
                pattern = pattern/float(sum(pattern))
            except:
                continue
            date = row['EXAMDATE']
            uptake_patterns[(rid, date)] = pattern
            summary[(rid, date)] = float(row['SUMMARYSUVR_WHOLECEREBNORM'])
    keys = []
    data = []
    for k,v in uptake_patterns.iteritems():
        keys.append(k)
        data.append(v)
    data = np.array(data)
    print data.shape

    # reverse translate regions
    lut_by_name = {b.lower().strip().replace('-','_') : int(a) for a,b in lut_dict.iteritems()}

    regions = [[lut_by_name[_.lower().strip()]] for _ in UPTAKE_KEYS]

    return keys, regions, data

def deriveRoussetClusterUptakePattern(rousset_data, master_data):
    # parse rousset data
    subj_list = unwrap(rousset_data['subj_list']) # list
    result_list = unwrap(rousset_data['result_list'])
    subj_data = {}
    for subj, result in zip(subj_list, result_list):
        result = parseResult(result)
        # unwrap 
        for k,v in result.iteritems():
            try:
                v = [unwrap(_) for _ in v]
            except:
                pass
            result[k] = v
        subj_data[subj] = result

    # parse uptake pattern data
    regions = None
    keys = []
    data = []
    for rid, result in subj_data.iteritems():
        if regions is None:
            regions = result['indices']
        
        master_row = master_data.get(rid)
        if not master_row:
            print "no master row"
            continue

        # filter out non-normals
        '''
        diag = master_row['Init_Diagnosis'].strip()
        if diag != 'N':
            continue
        '''
        # filter out positives at baseline
        if float(master_row['AV45_wcereb']) >= 1.11:
            continue


        # get baseline scandate
        data.append(result['pvcvals'])
        keys.append(rid)
    
    data = np.array(data)
    print data.shape
    return keys, regions, data



if __name__ == "__main__":
    nontp_file = '../output/UCBERKELEYAV45_09_25_15_extra_nontp.csv'
    tp_file = '../output/UCBERKELEYAV45_09_25_15_extra.csv'
    registry_file = "../docs/registry_clean.csv"
    master_file = '../FDG_AV45_COGdata_09_25_15.csv'

    lut_file = "../FreeSurferColorLUT.txt"
    lut_table = importFreesurferLookup(lut_file)

    registry = importRegistry(registry_file)
    tp_data = importAV45(tp_file, av45_nontp_file=None ,registry=registry)
    nontp_data = importAV45(nontp_file, av45_nontp_file=None, registry=registry)
    master_data = importMaster(master_file) 

    rousset_mat = "../raw_agghigh_output.mat"
    rousset_data = loadMATFile(rousset_mat)

    diags = {}
    for rid, row in master_data.iteritems():
        diag = row['Init_Diagnosis'].strip()
        diags[rid] = diag

    # get dataset
    keys, regions, data = deriveCorticalSummaryUptakePattern(nontp_data, master_data, lut_table) # using preprocessed data
    #keys, regions, data = deriveRoussetClusterUptakePattern(rousset_data, master_data) # using rousset output

    # PCA
    '''
    pca_model = PCA(n_components=3, copy=True, whiten=False)
    obs_pca = pca_model.fit_transform(data)
    x = [_[0] for _ in obs_pca]
    y = [_[1] for _ in obs_pca]
    z = [_[2] for _ in obs_pca]
    c = ['r' if diags[k]=='AD' else 'b' for k,m in keys]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z)
    plt.show()
    sys.exit(1)
    '''

    
    ks=range(1,30)
    gaps, valid_ks = gap(data, nrefs=100, ks=ks, use_pca=False)
    print "GAPS"
    print gaps
    print sorted(valid_ks, key=lambda x: x[1], reverse=True)
    fig = plt.figure(1)
    plt.plot(ks, gaps)
    plt.show()
    sys.exit(1)
    

    # Run GMM 
    '''
    g = GMM(n_components=1, 
              covariance_type='full', 
              random_state=None, 
              thresh=None, 
              tol=0.00001, 
              min_covar=0.00001, 
              n_iter=200,
              params='wmc', 
              init_params='wmc')
    '''
    g = DPGMM(n_components=1, 
              covariance_type='full', 
              alpha=10,
              tol=0.00001, 
              min_covar=0.00001, 
              n_iter=200,
              params='wmc', 
              init_params='wmc')
    g.fit(data)
    print np.round(g.weights_, 2)
    print np.round(g.means_, 2)
    try:
        print np.round(g.precs_, 2) 
    except:
        print np.round(g.covars_, 2)
    print g.converged_
    sys.exit(1)

    # Run KMeans
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

