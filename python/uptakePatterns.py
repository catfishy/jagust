from utils import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.mixture import GMM, VBGMM, DPGMM
from collections import defaultdict
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

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


def deriveCorticalSummaryUptakePattern(input_data, master_data, lut_dict):
    # only using nontp data
    bl_patterns = {}
    scan2_patterns = {}
    scan3_patterns = {}
    diagnoses = {}
    summary = {}
    for rid, rows in input_data.iteritems():
        master_row = master_data.get(rid)
        if not master_row:
            continue
        diagnoses[rid] = master_row['Init_Diagnosis'].strip()
        # filter out AD
        if diagnoses[rid] in set(['AD', 'LMCI']):
            continue
        for i, row in enumerate(sorted(rows, key=lambda x: x['EXAMDATE'])):

            # filter out positives at baseline
            if i == 0 and float(row['SUMMARYSUVR_WHOLECEREBNORM']) > 1.11:
                break

            wcereb = float(row['WHOLECEREBELLUM']) 
            comp_ref = float(row['COMPOSITE_REF'])
            # normalize by either comp_ref or wcereb
            try:
                pattern = np.array([float(row[_])/comp_ref for _ in UPTAKE_KEYS])
                pattern = pattern/float(sum(pattern))
            except:
                continue
            if i == 0:
                bl_patterns[rid] = pattern
                summkey = 'BL'
            elif i == 1:
                scan2_patterns[rid] = pattern
                summkey = 'Scan2'
            elif i == 2:
                scan3_patterns[rid] = pattern
                summkey = 'Scan3'
            summary[(rid, summkey)] = float(row['SUMMARYSUVR_WHOLECEREBNORM'])
    
    lut_by_name = {b.lower().strip() : int(a) for a,b in lut_dict.iteritems()}
    regions = [lut_by_name[_.lower().strip().replace('_','-')] for _ in UPTAKE_KEYS]

    bl_data = pd.DataFrame.from_dict(bl_patterns, orient='index')
    scan2_data = pd.DataFrame.from_dict(scan2_patterns, orient='index')
    scan3_data = pd.DataFrame.from_dict(scan3_patterns, orient='index')
    bl_data.columns = regions
    scan2_data.columns = regions
    scan3_data.columns = regions
    return summary, bl_data, scan2_data, scan3_data

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
    registry_file = "../docs/registry_clean.csv"
    master_file = '../FDG_AV45_COGdata_09_25_15.csv'
    lut_file = "../FreeSurferColorLUT.txt"
    lut_table = importFreesurferLookup(lut_file)
    registry = importRegistry(registry_file)
    master_data = importMaster(master_file) 
    diags = extractDiagnosesFromMasterData(master_data)
    
    # load av45 data 
    nontp_file = '../output/UCBERKELEYAV45_09_25_15_extra_nontp.csv'
    tp_file = '../output/UCBERKELEYAV45_09_25_15_extra.csv'
    tp_data = importAV45(tp_file, av45_nontp_file=None ,registry=registry)
    nontp_data = importAV45(nontp_file, av45_nontp_file=None, registry=registry)
    # rousset_mat = "../raw_agghigh_output.mat"
    # rousset_data = loadMATFile(rousset_mat)

    

    # get dataset
    summary, data_BL, data_Scan2, data_Scan3 = deriveCorticalSummaryUptakePattern(nontp_data, master_data, lut_table) # using preprocessed data
    data_all = pd.concat([data_BL, data_Scan2, data_Scan3]).as_matrix()

    pca_model = PCA(n_components=0.95, copy=True, whiten=False)
    pca_model.fit(data_all)

    # get years between BL and Scan2, annualized pattern transitions, annualized summary uptake changes
    yr_diff = {}
    uptake_diff = {}
    pattern_diff = {}
    for rid in data_Scan2.index:
        yrs = float(master_data[rid]['AV45_1_2_Diff (Yrs)'])
        first_pattern = pca_model.transform(data_BL.loc[rid].as_matrix())[0]
        second_pattern = pca_model.transform(data_Scan2.loc[rid].as_matrix())[0]
        pattern_delta = (second_pattern - first_pattern) / yrs
        uptake_delta = (summary[(rid, 'Scan2')] - summary[(rid, 'BL')]) / yrs
        yr_diff[rid] = yrs
        uptake_diff[rid] = abs(uptake_delta)
        pattern_diff[rid] = pattern_delta
    for rid in data_Scan3.index:
        yrs = float(master_data[rid]['AV45_1_3_Diff (yrs)'])
        first_pattern = pca_model.transform(data_BL.loc[rid].as_matrix())[0]
        second_pattern = pca_model.transform(data_Scan3.loc[rid].as_matrix())[0]
        pattern_delta = (second_pattern - first_pattern) / yrs
        uptake_delta = (summary[(rid, 'Scan3')] - summary[(rid, 'BL')]) / yrs
        yr_diff[rid] = yrs
        uptake_diff[rid] = abs(uptake_delta)
        pattern_diff[rid] = pattern_delta

    points = []
    for rid, uptake_diff in uptake_diff.iteritems():
        pattern_diff_abs = np.linalg.norm(pattern_diff[rid], ord=None)
        points.append((uptake_diff, pattern_diff_abs))
    x = [_[0] for _ in points]
    y = [_[1] for _ in points]
    plt.figure(1)
    plt.scatter(x,y)
    plt.show()


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

    # Choose number of clusters/components
    '''
    ks=range(1,40)
    gaps, valid_ks = gap(data, nrefs=50, ks=ks, use_pca=False)
    print "GAPS"
    print gaps
    print sorted(valid_ks, key=lambda x: x[1], reverse=True)
    fig = plt.figure(1)
    plt.plot(ks, gaps)
    plt.show()
    sys.exit(1)
    '''
    clusters=15


    # Run GMM 
    '''
    g = GMM(n_components=clusters, 
              covariance_type='full', 
              random_state=None, 
              thresh=None, 
              tol=0.00001, 
              min_covar=0.00001, 
              n_iter=200,
              params='wmc', 
              init_params='wmc')
    g.fit(data_all)
    print np.round(g.weights_, 2)
    print np.round(g.means_, 2)
    try:
        print np.round(g.precs_, 2) 
    except:
        print np.round(g.covars_, 2)
    print g.converged_
    print g.aic(data_all)
    print g.bic(data_all)

    # convert data to component responsibilities
    bl_resp = pd.DataFrame(g.predict_proba(data_BL))
    scan2_resp = pd.DataFrame(g.predict_proba(data_Scan2))
    scan3_resp = pd.DataFrame(g.predict_proba(data_Scan3))
    print bl_resp
    print scan2_resp
    print scan3_resp
    '''


    # Run KMeans
    '''
    model = KMeans(n_clusters=clusters, n_jobs=1, copy_x=True, verbose=True)
    labels = model.fit_predict(data)
    centers = model.cluster_centers_
    by_label = defaultdict(list)
    for i,l in enumerate(labels):
        by_label[l].append(keys[i])
    by_label = dict(by_label)

    for c in centers:
        print c
    '''
    # TODO: MAP CLUSTER CENTERS TO MODEL aparc+aseg, and visualize

