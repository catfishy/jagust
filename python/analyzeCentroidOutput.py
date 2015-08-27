from scipy.spatial.distance import euclidean
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import kneighbors_graph
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools
from scipy.sparse import csr_matrix

from utils import *

'''
NORMALIZE THE CENTROID DIFFERENCES BY SOME STANDARD MEASURE OF BRAIN SIZE
'''
def createGroupingFile(clusters, filepath):
    # convert back to indices
    names = []
    all_groups = []
    for c,v in clusters.iteritems():
        names.append(c)
        all_groups.append(v)
    createROIGrouping(names, all_groups, filepath)
    printROIGrouping(names, all_groups)


def parseResult(data, lut_table):
    '''
    region values are ordered by sorted subjects (low->high)
    '''
    region_list = data['region_list']
    regions = data['regions']
    subj_list = data['subj_list']
    subjects = data['subjects']

    while len(subjects) <= 1:
        subjects = subjects[0]
    while len(regions) <= 1:
        regions = regions[0]
    while len(region_list) <= 1:
        region_list = region_list[0]
    while len(subj_list) <= 1:
        subj_list = subj_list[0]

    by_region = {}
    for region_code, region in zip(regions, region_list):
        while len(region) <= 1:
            region = region[0]
        key = lut_table.get(region_code, region_code)
        by_region[lut_table.get(region_code, region_code)] = region

    by_subj = {}
    for subj_id, subj_result in zip(subjects, subj_list):
        while len(subj_result) <= 1:
            subj_result = subj_result[0]
        subj_data = {}
        for i, row in enumerate(subj_result):
            lut_index = i+1
            row = row[0]
            if len(row) == 0:
                continue
            lut_name = lut_table.get(lut_index, lut_index)
            fields = row.dtype.names
            rowdata = row[0]
            extracted = {}
            for k,v in zip(fields, rowdata):
                while (isinstance(v, np.ndarray) or isinstance(v, np.void)) and len(v) <= 1:
                    v = v[0]
                extracted[k] = v
            subj_data[lut_name] = extracted
        by_subj[subj_id] = subj_data

    return by_region, by_subj, subjects

def createConnectivityGraph(region_list):
    '''
    Separate ventricles from white matter from other
    '''
    ventricles = ['Left-Lateral-Ventricle','4th-Ventricle','CSF','3rd-Ventricle','5th-Ventricle','Left-Inf-Lat-Vent','Right-Inf-Lat-Vent']
    whitematter = ['Optic-Chiasm','CC_Central','Left-Cerebral-White-Matter','Left-Cerebellum-White-Matter','CC_Anterior','Right-Cerebellum-White-Matter','ctx-rh-corpuscallosum','CC_Mid_Posterior','CC_Posterior','WM-hypointensities','CC_Mid_Anterior','Right-Cerebral-White-Matter']
    #ambiguous = ['Optic-Chiasm','ctx-rh-unknown','ctx-lh-unknown']
    ambiguous = []
    others = [_ for _ in region_list if _ not in ventricles+whitematter+ambiguous]

    graph = np.zeros((len(region_list),len(region_list)))

    for v1, v2 in itertools.combinations_with_replacement(ventricles, 2):
        v1_index = region_list.index(v1)
        v2_index = region_list.index(v2)
        graph[v1_index, v2_index] = 1
        graph[v2_index, v1_index] = 1
    for v1, v2 in itertools.combinations_with_replacement(whitematter, 2):
        v1_index = region_list.index(v1)
        v2_index = region_list.index(v2)
        graph[v1_index, v2_index] = 1
        graph[v2_index, v1_index] = 1
    for v1, v2 in itertools.combinations_with_replacement(others, 2):
        v1_index = region_list.index(v1)
        v2_index = region_list.index(v2)
        graph[v1_index, v2_index] = 1
        graph[v2_index, v1_index] = 1
    '''
    for a in ambiguous:
        a_index = region_list.index(a)
        graph[a_index,:] = 1
        graph[:,a_index] = 1
    '''
    sparse = csr_matrix(graph)
    return connectivity

if __name__ == '__main__':
    output_mat = "../aggOutput.mat"
    data = loadMATFile(output_mat)
    lut_file = "../FreeSurferColorLUT.txt"
    lut_table = importFreesurferLookup(lut_file)

    by_region, by_subj, subjects = parseResult(data, lut_table)

    # find centroid of centroids + normalize it somehow!!! (divide by original image x,y,z)
    centroids = {}
    for r in by_region.keys():
        subj_centroids = []
        for s in subjects:
            cur_centroid = by_subj[s].get(r,{}).get('centroid', None)
            if cur_centroid is not None:
                subj_centroids.append(cur_centroid)
        centroids[r] = np.mean(np.array(subj_centroids), axis=0)

    binding_distributions = []
    centroid_vectors = []
    segments = []
    variances = {}
    means = {}
    for k,v in by_region.iteritems():
        v = np.nan_to_num(v)
        binding_distributions.append(v)
        centroid_vectors.append(centroids[k])
        segments.append(k)
        variances[k] = np.std(v)
        means[k] = np.mean(v)

    print segments

    '''
    sorted_means = sorted(means.items(), key=lambda x: x[1], reverse=True)
    for k,v in sorted_means:
        print "%s: %s +/- %s" % (k,v, variances[k])
    '''

    obs = np.array(binding_distributions)

    # Append centroids
    '''
    appended = []
    for i,x in enumerate(obs):
        appended.append(np.append(x, centroid_vectors[i]))
    obs = np.array(appended)
    '''

    # Scale (not really necessary as scans are already prenormalized)
    '''
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    obs = scaler.fit_transform(obs)
    '''

    '''
    # CHOOSE NUMBER OF CLUSTERS
    ks=range(18,60)
    gaps_raw, valid_ks_raw = gap(obs, nrefs=40, ks=ks)

    print "RAW GAPS"
    print gaps_raw
    print sorted(valid_ks_raw, key=lambda x: x[1], reverse=True)

    fig = plt.figure()
    plt.plot(ks, gaps_raw)
    plt.show()
    sys.exit(1)
    '''
    k = [26,55]
    
    # RUN HIERARCHICAL CLUSTERING
    # RUN K_MEANS
    scores = []


    # Get connectivity graph
    connectivity = createConnectivityGraph(segments)

    '''
    # PCA
    pca_model = PCA(n_components=0.99, copy=True, whiten=False)
    obs_pca = pca_model.fit_transform(obs)
    for x in pca_model.explained_variance_ratio_:
        print x
    '''
    # no pca
    obs_pca = obs

    print obs_pca.shape
    for cur_k in k:
        '''
        # KMEANS
        model = KMeans(n_clusters=cur_k, max_iter=1000, tol=0.00001,
                       n_jobs=-1, copy_x=True, verbose=False)
        labels = model.fit_predict(obs_pca)
        try:
            score = silhouette_score(obs,labels,metric='euclidean')
        except Exception as e:
            print e
            score = np.nan
        scores.append((cur_k, 'kmeans', score, model))
        '''
        
        # AGGLOMERATIVE
        model = AgglomerativeClustering(n_clusters=cur_k, 
                                        affinity='euclidean', 
                                        linkage='ward',
                                        connectivity=connectivity)
        labels = model.fit_predict(obs_pca)
        try:
            score = silhouette_score(obs,labels,metric='euclidean')
        except Exception as e:
            print e
            score = np.nan
        scores.append((cur_k, 'agg', score, model))
        
    best = sorted(scores, key=lambda x: x[2], reverse=True)[:40]
    for c in best:
        print "Clusters: %s, Neighbors: %s, Score: %s" % (c[0],c[1],c[2])
        labels = c[3].fit_predict(obs_pca)
        by_cluster_label = {_:[] for _ in list(set(labels))}
        for cluster, segment in zip(labels, segments):
            segment_index = [k for k,v in lut_table.iteritems() if v == segment][0]
            #by_cluster_label[cluster].append(segment_index)
            by_cluster_label[cluster].append(segment)
        roinames = []
        varnames = []
        for cluster, indices in by_cluster_label.iteritems():
            cluster_varname = 'Cluster%s' % cluster
            print '%s=%s' % (cluster_varname,indices)
            roinames.append("%s" % cluster_varname) 
            varnames.append(cluster_varname)
        print "all_groups=[%s]" % (','.join(varnames),)
        print "names=%s" % roinames
        print "\n\n"
        
