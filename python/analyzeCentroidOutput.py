from scipy.spatial.distance import euclidean
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import kneighbors_graph
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from utils import *

'''
NORMALIZE THE CENTROID DIFFERENCES BY SOME STANDARD MEASURE OF BRAIN SIZE
'''

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

    # Scale
    scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
    obs = scaler.fit_transform(obs)

    # PCA
    pca_model = PCA(n_components=0.96, copy=True, whiten=False)
    obs_pca = pca_model.fit_transform(obs)
    for x in pca_model.explained_variance_ratio_:
        print x
    print obs_pca.shape

    for op in obs_pca:
        print op


    # RUN HIERARCHICAL CLUSTERING
    # RUN K_MEANS
    scores = []
    for s in range(25,75):
        # connectivity (use centroids)
        connectivity = kneighbors_graph(centroid_vectors, n_neighbors=s, include_self=False)
        for k in range(15o,65):
            '''
            model = KMeans(n_clusters=k, n_jobs=-1, copy_x=True)
            labels = model.fit_predict(obs_pca)
            '''
            model = AgglomerativeClustering(n_clusters=k, 
                                            affinity='euclidean', 
                                            linkage='ward',
                                            connectivity=connectivity)
            labels = model.fit_predict(obs_pca)
            try:
                score = silhouette_score(obs_pca,labels,metric='euclidean')
            except Exception as e:
                print e
                score = np.nan
            scores.append((k, s, score, model))

    best = sorted(scores, key=lambda x: x[2], reverse=True)[0]
    labels = best[3].fit_predict(obs_pca)
    by_cluster_label = {_:[] for _ in list(set(labels))}
    for cluster, segment in zip(labels, segments):
        by_cluster_label[cluster].append(segment)

    for k,v in by_cluster_label.iteritems():
        print "Cluster %s" % k
        print "Segments: %s" % (v,)

    fig = plt.figure()
    ax = Axes3D(fig)
    xs = [_[0] for _ in scores]
    ys = [_[1] for _ in scores]
    zs = [_[2] for _ in scores]
    # put 0s on the y-axis, and put the y axis on the z-axis
    ax.scatter(xs=xs, ys=ys, zs=zs, zdir='z')
    #plt.plot([_[0] for _ in scores], [_[1] for _ in scores])
    plt.show()