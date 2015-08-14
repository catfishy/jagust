from scipy.spatial.distance import euclidean
from scipy.cluster.vq import kmeans, whiten, vq, kmeans2
from sklearn.metrics import silhouette_score
import numpy as np

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

    binding_distributions = []
    segments = []
    for k,v in by_region.iteritems():
        binding_distributions.append(v)
        segments.append(k)

    #obs = whiten(np.array(binding_distributions))
    obs = np.array(binding_distributions)
    
    scores = []
    for k in range(10,32):
        centroid, labels = kmeans2(data=obs, k=k, iter=10, minit='points')
        print labels
        try:
            score = silhouette_score(obs,labels,metric='euclidean')
        except:
            score = np.nan
        scores.append((k,score))
    for s in scores:
        print s

    sys.exit(1)

    by_cluster_label = {_:[] for _ in list(set(labels))}
    for cluster, segment in zip(labels, segments):
        by_cluster_label[cluster].append(segment)

    for k,v in by_cluster_label.iteritems():
        print "Cluster %s" % k
        print "Segments: %s" % (v,)
