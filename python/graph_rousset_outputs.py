import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from collections import defaultdict
from scipy.stats import norm, mannwhitneyu

from utils import *

AD = [4009, 5241, 4338, 5067, 4676, 4732, 4379, 4172, 5138, 4924, 4755, 4280, 5224, 4223, 4772, 5120, 4911, 4568, 4827, 4730, 5032, 5196, 4152, 5028, 4692, 4863, 5187, 4195, 4307, 4892, 4211, 4879, 4853, 4962, 5146, 4910, 4583, 5119, 4657]
NORMAL = [4119, 4499, 4427, 4393, 4269, 4148, 4643, 4222, 4084, 5040, 4387, 4609, 4516, 4367, 4762, 4337, 4496, 4139, 4580, 4357, 4020, 4177, 4488, 4224, 4382, 4599, 4424, 4586, 4279, 4585, 4173, 4428, 4097, 4090, 4921, 4835, 4449, 4598, 4637, 4555, 4441, 4255, 4060, 4576, 4644, 4612, 4290, 4093, 4578, 4587, 4100, 4474, 4120, 4385, 4508, 4616, 4151, 4010, 4291, 4277, 4433]


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

def subplot_scatter(data, key, grouping, plot_raw=False):
    scatter_points = []
    for subj, subjdata in data:
        if plot_raw and 'suvr' in key:
            startval = np.mean([_['%s_nonpvc' % key]*_['ref_nonpvc'] for _ in subjdata])
            tograph = [_['%s_pvc' % key]*_['ref_pvc'] for _ in subjdata]
        else:
            startval = np.mean([_['%s_nonpvc' % key] for _ in subjdata])
            tograph = [_['%s_pvc' % key] for _ in subjdata]
        color = 'b'
        if subj in NORMAL:
            color = 'g'
        elif subj in AD:
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


def slope_graph(data, key, grouping, plot_raw=False):
    points = {'b': [], 'g': [], 'r': []}
    for subj, subjdata in data:
        if plot_raw and 'suvr' in key:
            startval = np.mean([_['%s_nonpvc' % key]*_['ref_nonpvc'] for _ in subjdata])
            tograph = [_['%s_pvc' % key]*_['ref_pvc'] for _ in subjdata]
        else:
            startval = np.mean([_['%s_nonpvc' % key] for _ in subjdata])
            tograph = [_['%s_pvc' % key] for _ in subjdata]
        color = 'b'
        if subj in NORMAL:
            color = 'g'
        elif subj in AD:
            color = 'r'
        # scatter plot
        slope = tograph[grouping] - startval
        points[color].append(slope)
    all_values = [v for _ in points.values() for v in _]
    min_v = min(all_values)
    max_v = max(all_values)

    u, pvalue = mannwhitneyu(points['r'], points['g'], use_continuity=True)
    u_max = len(points['g']) * len(points['r'])
    rank_biserial = 1.0 - (2*u/u_max)
    print key
    print "%s/%s" % (u,u_max)
    print pvalue
    print rank_biserial

    stats = {}
    for k,v in points.iteritems():
        if v:
            weights = np.ones_like(v)/len(v)
            n, bins, patches = plt.hist(v, 80, weights=weights, facecolor=k, alpha=0.30)
            (mu, sigma) = norm.fit(v)
            print mu
            print sigma
            x = np.linspace(min_v, max_v, 100)
            y = mlab.normpdf(x, mu, sigma)
            weights = np.ones_like(y)/len(y)
            y = np.multiply(y,weights)
            plt.plot(x, y, '%s--' % k, linewidth=3)
            if k == 'g': # Normal
                stats['N'] = (mu, sigma)
            elif k == 'r': # AD
                stats['AD'] = (mu, sigma)
    full_text = "$\mathrm{N:}\ \mu=%.3f,\ \sigma=%.3f$\n$\mathrm{AD:}\ \mu=%.3f,\ \sigma=%.3f$" % (stats['N'][0], stats['N'][1], stats['AD'][0], stats['AD'][1])
    mwu_text = "$U=%i$\n$\\rho=%.3e$\n$r_{RankBiserial}=%.3f$\n$n_N=%i$\n$n_{AD}=%i$" % (u, pvalue, rank_biserial, len(points['g']), len(points['r']))
    plt.xlim(min_v, max_v)
    props = dict(boxstyle='round,pad=0.8', facecolor='wheat', alpha=0.7)
    plt.annotate(full_text, xy=(0.05,0.87), xycoords='axes fraction', bbox=props, fontsize=13)
    plt.annotate(mwu_text, xy=(0.6,0.75), xycoords='axes fraction', bbox=props, fontsize=13)
    plt.xlabel('(PVC value) - (Non-PVC value)')
    plt.ylabel('Frequency')

def mean_graph(data, key, grouping, plot_raw=False):
    points = {'b': [], 'g': [], 'r': []}
    for subj, subjdata in data:
        if plot_raw and 'suvr' in key:
            startval = np.mean([_['%s_nonpvc' % key]*_['ref_nonpvc'] for _ in subjdata])
            tograph = [_['%s_pvc' % key]*_['ref_pvc'] for _ in subjdata]
        else:
            startval = np.mean([_['%s_nonpvc' % key] for _ in subjdata])
            tograph = [_['%s_pvc' % key] for _ in subjdata]
        color = 'b'
        if subj in NORMAL:
            color = 'g'
        elif subj in AD:
            color = 'r'
        # scatter plot
        val = tograph[grouping]
        points[color].append(val)
    all_values = [v for _ in points.values() for v in _]
    min_v = min(all_values)
    max_v = max(all_values)

    u, pvalue = mannwhitneyu(points['r'], points['g'], use_continuity=True)
    u_max = len(points['g']) * len(points['r'])
    rank_biserial = 1.0 - (2*u/u_max)
    print key
    print "%s/%s" % (u,u_max)
    print pvalue
    print rank_biserial

    stats = {}
    for k,v in points.iteritems():
        if v:
            weights = np.ones_like(v)/len(v)
            n, bins, patches = plt.hist(v, 80, weights=weights, facecolor=k, alpha=0.30)
            (mu, sigma) = norm.fit(v)
            x = np.linspace(min_v, max_v, 100)
            y = mlab.normpdf(x, mu, sigma)
            weights = np.ones_like(y)/len(y)
            y = np.multiply(y,weights)
            plt.plot(x, y, '%s--' % k, linewidth=3)
            if k == 'g': # Normal
                stats['N'] = (mu, sigma)
            elif k == 'r': # AD
                stats['AD'] = (mu, sigma)
    full_text = "$\mathrm{N:}\ \mu=%.3f,\ \sigma=%.3f$\n$\mathrm{AD:}\ \mu=%.3f,\ \sigma=%.3f$" % (stats['N'][0], stats['N'][1], stats['AD'][0], stats['AD'][1])
    mwu_text = "$U=%i$\n$\\rho=%.3e$\n$r_{RankBiserial}=%.3f$\n$n_N=%i$\n$n_{AD}=%i$" % (u, pvalue, rank_biserial, len(points['g']), len(points['r']))
    plt.xlim(min_v, max_v)
    props = dict(boxstyle='round,pad=0.8', facecolor='wheat', alpha=0.7)
    plt.annotate(full_text, xy=(0.05,0.87), xycoords='axes fraction', bbox=props, fontsize=13)
    plt.annotate(mwu_text, xy=(0.6,0.75), xycoords='axes fraction', bbox=props, fontsize=13)
    plt.xlabel('PVC')
    plt.ylabel('Frequency')

def subplot_scatter_between_group(data, key, groupA, groupB):
    scatter_points = []
    for subj, subjdata in data:
        groupA_val = subjdata[groupA]['%s_pvc' % key]
        groupB_val = subjdata[groupB]['%s_pvc' % key]
        color = 'b'
        if subj in NORMAL:
            color = 'g'
        elif subj in AD:
            color = 'r'
        # scatter plot
        scatter_points.append(([groupA_val], [groupB_val], color))
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    plt.xlabel('GROUPING %s' % (groupA+1,))
    plt.ylabel('GROUPING %s' % (groupB+1,))

def subplot_residuals(data):
    by_group = defaultdict(list)
    scatter_points = []
    for subj, subjdata in data:
        tograph = [_['rmse_residual'] for _ in subjdata]
        color = 'b'
        if subj in NORMAL:
            color = 'g'
            xmod = -0.07
        elif subj in AD:
            color = 'r'
            xmod = 0.07
        for i,v in enumerate(tograph):
            by_group[i].append(v)
            actual_x = i+1+xmod
            scatter_points.append((actual_x,v,color))
    avgs = [np.mean(by_group[i]) for i,g in enumerate(groups)]
    stds = [np.std(by_group[i]) for i,g in enumerate(groups)]
    plt.errorbar(groups, avgs, stds)
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    x1,x2,y1,y2 = plt.axis()
    plt.axis((0,5,y1,y2))
    plt.xlabel('GROUPINGS')
    plt.ylabel('PVC Residuals')

def addCompositeSUVR(result_list):
    # add composite
    for subjdata in result_list:
        for groupdata in subjdata:
            groupdata['composite_suvr_nonpvc'] = np.mean([groupdata['parietal_suvr_nonpvc'],groupdata['frontal_suvr_nonpvc'],groupdata['temporal_suvr_nonpvc'],groupdata['cingulate_suvr_nonpvc']])
            groupdata['composite_suvr_pvc'] = np.mean([groupdata['parietal_suvr_pvc'],groupdata['frontal_suvr_pvc'],groupdata['temporal_suvr_pvc'],groupdata['cingulate_suvr_pvc']])
            groupdata['composite_wm_suvr_nonpvc'] = np.mean([groupdata['parietal_suvr_nonpvc'],groupdata['frontal_suvr_nonpvc'],groupdata['temporal_suvr_nonpvc'],groupdata['cingulate_suvr_nonpvc']] * groupdata['ref_nonpvc'] / groupdata['hemiWM_nonpvc'])
            groupdata['composite_wm_suvr_pvc'] = np.mean([groupdata['parietal_suvr_pvc'],groupdata['frontal_suvr_pvc'],groupdata['temporal_suvr_pvc'],groupdata['cingulate_suvr_pvc']] * groupdata['ref_nonpvc'] / groupdata['hemiWM_nonpvc'])
    return result_list


if __name__ == "__main__":
    rousset_mat = '../rousset_outputs_080515.mat'
    data = loadMATFile(rousset_mat)
    result_list = [parseResult(_) for _ in data['result_list'][0]]
    subj_list = data['subj_list'][0]
    groups = list(data['groups'][0])

    # sort by temporal_suvr_nonpvc
    result_list = addCompositeSUVR(result_list)
    sorted_data = zip(subj_list,result_list)

    
    '''
    plot the residuals for each group
    '''
    '''
    plt.figure(1)
    plt.title('Residuals')
    subplot_residuals(sorted_data)
    plt.show()
    '''
    
    '''
    plot pvc value distribution, by group
    '''
    '''
    grouping=3
    subplot_height = 2
    subplot_length = 5
    subplot_indices = range(1,(subplot_height*subplot_length)+1)
    subplot_order = ['ref', 'brainstem', 'hemiWM', 'composite_wm_suvr', 'temporal_suvr', 'frontal_suvr', 'occipital_suvr', 'cingulate_suvr', 'composite_suvr', 'parietal_suvr']
    plt.figure(1)
    for subplot_i, subplot_key in zip(subplot_indices, subplot_order):
        ax = plt.subplot(subplot_height, subplot_length, subplot_i)
        plt.title(subplot_key.replace('_',' ').replace('suvr', 'SUVR'))
        mean_graph(sorted_data, subplot_key, grouping, plot_raw=False)
    plt.show()
    sys.exit(1)
    '''

    '''
    plot pvc-nonpvc slope, by group
    '''
    grouping=3
    subplot_height = 2
    subplot_length = 5
    subplot_indices = range(1,(subplot_height*subplot_length)+1)
    subplot_order = ['ref', 'brainstem', 'hemiWM', 'composite_wm_suvr', 'temporal_suvr', 'frontal_suvr', 'occipital_suvr', 'cingulate_suvr', 'composite_suvr', 'parietal_suvr']
    plt.figure(1)
    for subplot_i, subplot_key in zip(subplot_indices, subplot_order):
        ax = plt.subplot(subplot_height, subplot_length, subplot_i)
        plt.title(subplot_key.replace('_',' ').replace('suvr', 'SUVR'))
        slope_graph(sorted_data, subplot_key, grouping, plot_raw=False)
    plt.show()
    sys.exit(1)

    '''
    plot pvc vs nonpvc, by group
    '''
    # grouping to plot, by index of group (group num - 1)
    grouping=3
    subplot_height = 2
    subplot_length = 5
    subplot_indices = range(1,(subplot_height*subplot_length)+1)
    subplot_order = ['ref', 'brainstem', 'hemiWM', 'composite_wm_suvr', 'temporal_suvr', 'frontal_suvr', 'occipital_suvr', 'cingulate_suvr', 'composite_suvr', 'parietal_suvr']
    plt.figure(1)
    for subplot_i, subplot_key in zip(subplot_indices, subplot_order):
        ax = plt.subplot(subplot_height, subplot_length, subplot_i)
        plt.title(subplot_key.replace('_',' ').replace('suvr', 'SUVR'))
        subplot_scatter(sorted_data, subplot_key, grouping)
        if 'suvr' in subplot_key or subplot_key == 'ref':
            ax.set_xlim([0.2, 1.8])
            ax.set_ylim([0.2, 2.2])
        elif subplot_key in set(['brainstem', 'hemiWM']):
            ax.set_xlim([1.0, 3.0])
            ax.set_ylim([1.0, 4.0])
        x1,x2,y1,y2 = plt.axis()
        plt.plot([x1,x2],[x1,x2],color='b')
    plt.show()
    

    '''
    plot pvc vs nonpvc, by group, non-suvr values
    '''
    '''
    # grouping to plot, by index of group (group num - 1)
    grouping=3
    subplot_height = 2
    subplot_length = 5
    subplot_indices = range(1,(subplot_height*subplot_length)+1)
    subplot_order = ['ref', 'brainstem', 'hemiWM', 'basal_ganglia_suvr', 'temporal_suvr', 'frontal_suvr', 'occipital_suvr', 'cingulate_suvr', 'composite_suvr', 'parietal_suvr']
    plt.figure(1)
    for subplot_i, subplot_key in zip(subplot_indices, subplot_order):
        plt.subplot(subplot_height, subplot_length, subplot_i)
        plt.title(subplot_key.replace('_',' ').replace('suvr', '').strip())
        subplot_scatter(sorted_data, subplot_key, grouping, plot_raw=True)
    plt.show()
    '''


    '''
    plot group A vs group B, pvc
    '''
    '''
    groupA = 0
    groupB = 3
    subplot_height = 2
    subplot_length = 5
    subplot_indices = range(1,(subplot_height*subplot_length)+1)
    subplot_order = ['ref', 'brainstem', 'hemiWM', 'basal_ganglia_suvr', 'temporal_suvr', 'frontal_suvr', 'occipital_suvr', 'cingulate_suvr', 'composite_suvr', 'parietal_suvr']

    plt.figure(1)
    for subplot_i, subplot_key in zip(subplot_indices, subplot_order):
        plt.subplot(subplot_height, subplot_length, subplot_i)
        plt.title(subplot_key.replace('_',' ').replace('suvr', 'SUVR'))
        subplot_scatter_between_group(sorted_data, subplot_key, groupA, groupB)
    plt.show()
    '''

