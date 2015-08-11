from scipy.spatial.distance import euclidean
from scipy.cluster.vq import kmeans, whiten, vq
from utils import *

'''
NORMALIZE THE CENTROID DIFFERENCES BY SOME STANDARD MEASURE OF BRAIN SIZE
'''


def parseResult(data, lut_table):
    results = data['results']
    while len(results) <= 1:
        results = results[0]
    to_return = {}
    for i, row in enumerate(results):
        row = row[0]
        lut_index = i+1
        lut_name = lut_table.get(lut_index, lut_index)
        if len(row) == 0:
            continue
        fields = row.dtype.names
        rowdata = row[0]
        extracted = {}
        for k,v in zip(fields, rowdata):
            while (isinstance(v, np.ndarray) or isinstance(v, np.void)) and len(v) <= 1:
                v = v[0]
            print type(v)
            extracted[k] = v
        to_return[lut_name] = extracted
    return to_return




if __name__ == '__main__':
    output_mat = "../4269.mat"
    data = loadMATFile(output_mat)
    lut_file = "../FreeSurferColorLUT.txt"
    lut_table = importFreesurferLookup(lut_file)

    translated = parseResult(data, lut_table)
    rows = [(k,v) for k,v in translated.iteritems()]
    '''
    for k,v in translated.iteritems():
        print k
        print v
    '''

    centroids = [v['centroid'] for k,v in rows]
    obs = whiten(np.array(centroids))

    codebook, distortion = kmeans(obs=obs, k_or_guess=30)
    code, dist = vq(obs, codebook)

    print codebook
    print code
    print dist

    by_cluster = {_:[] for _ in list(set(code))}
    for label, row in zip(code, rows):
        by_cluster[label].append(row)

    for k,v in by_cluster.iteritems():
        print "Cluster %s" % k
        print "Parcels: %s" % ([a for a,b in v],)
