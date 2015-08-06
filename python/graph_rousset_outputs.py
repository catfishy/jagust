import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

from utils import *



def parseResult(arg):
    while type(arg) != np.ndarray or len(arg) <= 1:
        arg = arg[0]
    to_return = []
    for i, a in enumerate(arg):
        fields = a.dtype.names
        data = a[0][0]
        extracted = {}
        for k,v in zip(fields,data):
            while isinstance(v, np.ndarray) and len(v) <= 1:
                v = v[0]
            extracted[k] = v
        to_return.append(extracted)
    return to_return



if __name__ == "__main__":
    ad = [4009, 5241, 4338, 5067, 4676, 4732, 4379, 4172, 5138, 4924, 4755, 4280, 5224, 4223, 4772, 5120, 4911, 4568, 4827, 4730, 5032, 5196, 4152, 5028, 4692, 4863, 5187, 4195, 4307, 4892, 4211, 4879, 4853, 4962, 5146, 4910, 4583, 5119, 4657]
    normal = [4119, 4499, 4427, 4393, 4269, 4148, 4643, 4222, 4084, 5040, 4387, 4609, 4516, 4367, 4762, 4337, 4496, 4139, 4580, 4357, 4020, 4177, 4488, 4224, 4382, 4599, 4424, 4586, 4279, 4585, 4173, 4428, 4097, 4090, 4921, 4835, 4449, 4598, 4637, 4555, 4441, 4255, 4060, 4576, 4644, 4612, 4290, 4093, 4578, 4587, 4100, 4474, 4120, 4385, 4508, 4616, 4151, 4010, 4291, 4277, 4433]

    rousset_mat = '../rousset_outputs_080515.mat'
    data = loadMATFile(rousset_mat)
    result_list = [parseResult(_) for _ in data['result_list'][0]]
    subj_list = data['subj_list'][0]
    groups = list(data['groups'][0])

    # sort by temporal_suvr_nonpvc
    sorted_data = sorted(zip(subj_list,result_list), key=lambda x: np.mean([_['temporal_suvr_nonpvc'] for _ in x[1]]))

    print [_[0] for _ in sorted_data]

    only_graph=[0,2,4]
    grouping=3

    plt.figure(1)
    plt.subplot(251)
    plt.title('Temporal SUVR')
    scatter_points = []
    for subj, subjdata in sorted_data:
        startval = np.mean([_['temporal_suvr_nonpvc'] for _ in subjdata])
        tograph = [_['temporal_suvr_pvc'] for _ in subjdata]
        color = 'b'
        if subj in normal:
            color = 'g'
        elif subj in ad:
            color = 'r'
        # line plot
        '''
        x = [0] + groups
        y = [startval] + tograph
        plt.plot([x[_] for _ in only_graph],[y[_] for _ in only_graph], color)
        plt.xticks(only_graph)
        '''
        # scatter plot
        scatter_points.append((startval, tograph[grouping], color))
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    x1,x2,y1,y2 = plt.axis()
    plt.plot([x1,x2],[x1,x2])
    plt.xlabel('NON-PVC')
    plt.ylabel('PVC')

    #plt.figure(2)
    '''
    plt.subplot(252)
    plt.title('Residuals')
    by_group = defaultdict(list)
    for subj, subjdata in sorted_data:
        tograph = [_['rmse_residual'] for _ in subjdata]
        for i,v in enumerate(tograph):
            by_group[i].append(v)
        color = 'b'
        if subj in normal:
            color = 'g'
        elif subj in ad:
            color = 'r'
        plt.plot(groups, tograph, color)
    avgs = [np.mean(by_group[i]) for i,g in enumerate(groups)]
    stds = [np.std(by_group[i]) for i,g in enumerate(groups)]
    print stds
    plt.errorbar(groups, avgs, stds)
    x1,x2,y1,y2 = plt.axis()
    plt.axis((0,5,y1,y2))
    '''

    #plt.figure(9)
    plt.subplot(252)
    plt.title('hemiWM')
    scatter_points = []
    for subj, subjdata in sorted_data:
        startval = np.mean([_['hemiWM_nonpvc'] for _ in subjdata])
        tograph = [_['hemiWM_pvc'] for _ in subjdata]
        color = 'b'
        if subj in normal:
            color = 'g'
        elif subj in ad:
            color = 'r'
        '''
        x = [0] + groups
        y = [startval] + tograph
        plt.plot([x[_] for _ in only_graph],[y[_] for _ in only_graph], color)
        plt.xticks(only_graph)
        '''
        # scatter plot
        scatter_points.append((startval, tograph[grouping], color))
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    x1,x2,y1,y2 = plt.axis()
    plt.plot([x1,x2],[x1,x2])
    plt.xlabel('NON-PVC')
    plt.ylabel('PVC')

    #plt.figure(3)
    plt.subplot(253)
    plt.title('Whole Cereb Ref')
    scatter_points = []
    for subj, subjdata in sorted_data:
        startval = np.mean([_['ref_nonpvc'] for _ in subjdata])
        tograph = [_['ref_pvc'] for _ in subjdata]
        color = 'b'
        if subj in normal:
            color = 'g'
        elif subj in ad:
            color = 'r'
        '''
        x = [0] + groups
        y = [startval] + tograph
        plt.plot([x[_] for _ in only_graph],[y[_] for _ in only_graph], color)
        plt.xticks(only_graph)
        '''
        # scatter plot
        scatter_points.append((startval, tograph[grouping], color))
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    x1,x2,y1,y2 = plt.axis()
    plt.plot([x1,x2],[x1,x2])
    plt.xlabel('NON-PVC')
    plt.ylabel('PVC')

    #plt.figure(4)
    plt.subplot(254)
    plt.title('Composite SUVR')
    scatter_points = []
    for subj, subjdata in sorted_data:
        rawcomposite = np.mean([np.mean([_['parietal_suvr_nonpvc'],_['frontal_suvr_nonpvc'],_['temporal_suvr_nonpvc'],_['cingulate_suvr_nonpvc']])*_['ref_nonpvc'] for _ in subjdata])
        startval = np.mean([np.mean([_['parietal_suvr_nonpvc'],_['frontal_suvr_nonpvc'],_['temporal_suvr_nonpvc'],_['cingulate_suvr_nonpvc']]) for _ in subjdata])

        print subj
        print "frontal raw avg: %s" % np.mean([_['frontal_suvr_nonpvc']*_['ref_nonpvc'] for _ in subjdata])
        print "cingulate raw avg: %s" % np.mean([_['cingulate_suvr_nonpvc']*_['ref_nonpvc'] for _ in subjdata])
        print "parietal raw avg: %s" % np.mean([_['parietal_suvr_nonpvc']*_['ref_nonpvc'] for _ in subjdata])
        print "temporal raw avg: %s" % np.mean([_['temporal_suvr_nonpvc']*_['ref_nonpvc'] for _ in subjdata])
        print "composite raw avg: %s" % (rawcomposite)
        print "cerebgm ref: %s" % np.mean([_['cerebGM_nonpvc'] for _ in subjdata])
        print "leftcereb: %s" % np.mean([_['left_wholecereb_nonpvc'] for _ in subjdata])
        print "rightcereb: %s" % np.mean([_['right_wholecereb_nonpvc'] for _ in subjdata])
        print "wcereb ref: %s" % np.mean([_['ref_nonpvc'] for _ in subjdata])
        
        print "cerebwm ref: %s" % np.mean([_['cerebWM_nonpvc'] for _ in subjdata])
        print "brainstem ref: %s" % np.mean([_['brainstem_nonpvc'] for _ in subjdata])
        

        tograph = [np.mean([_['parietal_suvr_pvc'],_['frontal_suvr_pvc'],_['temporal_suvr_pvc'],_['cingulate_suvr_pvc']]) for _ in subjdata]
        color = 'b'
        if subj in normal:
            color = 'g'
        elif subj in ad:
            color = 'r'
        '''
        x = [0] + groups
        y = [startval] + tograph
        plt.plot([x[_] for _ in only_graph],[y[_] for _ in only_graph], color)
        plt.xticks(only_graph)
        '''
        # scatter plot
        scatter_points.append((startval, tograph[grouping], color))
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    x1,x2,y1,y2 = plt.axis()
    plt.plot([x1,x2],[x1,x2])
    plt.xlabel('NON-PVC')
    plt.ylabel('PVC')

    #plt.figure(5)
    plt.subplot(255)
    plt.title('Occipital SUVR')
    scatter_points = []
    for subj, subjdata in sorted_data:
        startval = np.mean([_['occipital_suvr_nonpvc'] for _ in subjdata])
        tograph = [_['occipital_suvr_pvc'] for _ in subjdata]
        color = 'b'
        if subj in normal:
            color = 'g'
        elif subj in ad:
            color = 'r'
        '''
        x = [0] + groups
        y = [startval] + tograph
        plt.plot([x[_] for _ in only_graph],[y[_] for _ in only_graph], color)
        plt.xticks(only_graph)
        '''
        # scatter plot
        scatter_points.append((startval, tograph[grouping], color))
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    x1,x2,y1,y2 = plt.axis()
    plt.plot([x1,x2],[x1,x2])
    plt.xlabel('NON-PVC')
    plt.ylabel('PVC')

    #plt.figure(6)
    plt.subplot(256)
    plt.title('Cingulate SUVR')
    scatter_points = []
    for subj, subjdata in sorted_data:
        startval = np.mean([_['cingulate_suvr_nonpvc'] for _ in subjdata])
        tograph = [_['cingulate_suvr_pvc'] for _ in subjdata]
        color = 'b'
        if subj in normal:
            color = 'g'
        elif subj in ad:
            color = 'r'
        '''
        x = [0] + groups
        y = [startval] + tograph
        plt.plot([x[_] for _ in only_graph],[y[_] for _ in only_graph], color)
        plt.xticks(only_graph)
        '''
        # scatter plot
        scatter_points.append((startval, tograph[grouping], color))
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    x1,x2,y1,y2 = plt.axis()
    plt.plot([x1,x2],[x1,x2])
    plt.xlabel('NON-PVC')
    plt.ylabel('PVC')

    #plt.figure(7)
    plt.subplot(257)
    plt.title('Frontal SUVR')
    scatter_points = []
    for subj, subjdata in sorted_data:
        startval = np.mean([_['frontal_suvr_nonpvc'] for _ in subjdata])
        tograph = [_['frontal_suvr_pvc'] for _ in subjdata]
        color = 'b'
        if subj in normal:
            color = 'g'
        elif subj in ad:
            color = 'r'
        '''
        x = [0] + groups
        y = [startval] + tograph
        plt.plot([x[_] for _ in only_graph],[y[_] for _ in only_graph], color)
        plt.xticks(only_graph)
        '''
        # scatter plot
        scatter_points.append((startval, tograph[grouping], color))
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    x1,x2,y1,y2 = plt.axis()
    plt.plot([x1,x2],[x1,x2])
    plt.xlabel('NON-PVC')
    plt.ylabel('PVC')

    #plt.figure(8)
    plt.subplot(258)
    plt.title('Parietal SUVR')
    scatter_points = []
    for subj, subjdata in sorted_data:
        startval = np.mean([_['parietal_suvr_nonpvc'] for _ in subjdata])
        tograph = [_['parietal_suvr_pvc'] for _ in subjdata]
        color = 'b'
        if subj in normal:
            color = 'g'
        elif subj in ad:
            color = 'r'
        '''
        x = [0] + groups
        y = [startval] + tograph
        plt.plot([x[_] for _ in only_graph],[y[_] for _ in only_graph], color)
        plt.xticks(only_graph)
        '''
        # scatter plot
        scatter_points.append((startval, tograph[grouping], color))
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    x1,x2,y1,y2 = plt.axis()
    plt.plot([x1,x2],[x1,x2])
    plt.xlabel('NON-PVC')
    plt.ylabel('PVC')

    #plt.figure(9)
    plt.subplot(259)
    plt.title('Brainstem')
    scatter_points = []
    for subj, subjdata in sorted_data:
        startval = np.mean([_['brainstem_nonpvc'] for _ in subjdata])
        tograph = [_['brainstem_pvc'] for _ in subjdata]
        color = 'b'
        if subj in normal:
            color = 'g'
        elif subj in ad:
            color = 'r'
        '''
        x = [0] + groups
        y = [startval] + tograph
        plt.plot([x[_] for _ in only_graph],[y[_] for _ in only_graph], color)
        plt.xticks(only_graph)
        '''
        # scatter plot
        scatter_points.append((startval, tograph[grouping], color))
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    x1,x2,y1,y2 = plt.axis()
    plt.plot([x1,x2],[x1,x2])
    plt.xlabel('NON-PVC')
    plt.ylabel('PVC')

    #plt.figure(10)
    plt.subplot(2,5,10)
    plt.title('Basal Ganglia SUVR')
    scatter_points = []
    for subj, subjdata in sorted_data:
        startval = np.mean([_['basal_ganglia_suvr_nonpvc'] for _ in subjdata])
        tograph = [_['basal_ganglia_suvr_pvc'] for _ in subjdata]
        color = 'b'
        if subj in normal:
            color = 'g'
        elif subj in ad:
            color = 'r'
        '''
        x = [0] + groups
        y = [startval] + tograph
        plt.plot([x[_] for _ in only_graph],[y[_] for _ in only_graph], color)
        plt.xticks(only_graph)
        '''
        # scatter plot
        scatter_points.append((startval, tograph[grouping], color))
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    x1,x2,y1,y2 = plt.axis()
    plt.plot([x1,x2],[x1,x2])
    plt.xlabel('NON-PVC')
    plt.ylabel('PVC')

    plt.show()

    '''
    print result_list
    print subj_list
    print groups
    '''

