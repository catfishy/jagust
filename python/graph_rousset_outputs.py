import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from collections import defaultdict
from scipy.stats import norm, mannwhitneyu

from utils import *


def zipFields(args):
    '''
    For structs imported from mat files,
    where field values are in a list and the field names are listed under dtype
    '''
    fields = args.dtype.names
    values = []
    for i,v in enumerate(args):
        try:
            v = unwrap(v)
            if len(v) > 1:
                v = zipFields(v)
        except Exception as e:
            pass
        values.append(v)
    return dict(zip(fields,values))

def parseResult(arg):
    arg = unwrap(arg)
    data = zipFields(arg)
    return data

'''
def parseResult(arg):
    while type(arg) != np.ndarray or len(arg) <= 1:
        print type(arg)
        print len(arg)
        if len(arg) > 1:
            print arg
        arg = arg[0]
    print arg
    to_return = []  
    for i, a in enumerate(arg):
        fields = a.dtype.names
        print a
        data = a[0][0]
        extracted = {}
        for k,v in zip(fields,data):
            while isinstance(v, np.ndarray) and len(v) <= 1:
                v = v[0]
            extracted[k] = v
        to_return.append(extracted)
    return to_return
'''

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

    for subj, subjdata in data.iteritems():
        groupdata = subjdata[grouping]
        roidata = groupdata[key]
        refdata = groupdata['wholecereb']

        pvcval = roidata['pvcval']
        non_pvcval = roidata['nonpvcval']

        # normalize
        if not plot_raw and key in set(['composite']):
            pvcval /= float(refdata['pvcval'])
            non_pvcval /= float(refdata['nonpvcval'])

        color = 'b'
        if subj in NORMAL:
            color = 'g'
        elif subj in AD:
            color = 'r'

        # scatter plot
        points[color].append(pvcval)
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
    groups = ['group2', 'group4', 'agglow', 'agghigh']
    for subj, subjdata in data.iteritems():
        for groupname, groupdata in subjdata.iteritems():
            color = 'b'
            if subj in NORMAL:
                color = 'g'
                xmod = -0.07
            elif subj in AD:
                color = 'r'
                xmod = 0.07
            else:
                continue
            index = groups.index(groupname)
            v = groupdata['residuals']['rmse']
            by_group[index].append(v)
            actual_x = index+1+xmod
            scatter_points.append((actual_x,v,color))
    avgs = [np.mean(by_group[i]) for i,g in enumerate(groups)]
    stds = [np.std(by_group[i]) for i,g in enumerate(groups)]
    plt.errorbar(range(1,len(groups)+1), avgs, stds)
    plt.scatter([_[0] for _ in scatter_points], [_[1] for _ in scatter_points],c=[_[2] for _ in scatter_points])
    x1,x2,y1,y2 = plt.axis()
    plt.axis((0,len(groups)+1,y1,y2))
    plt.xlabel('GROUPINGS')
    plt.ylabel('PVC Residuals')

def composite_hemiWM_scatter(data, grouping):
    pvc_points = []
    nonpvc_points = []

    for subj, subjdata in data.iteritems():
        groupdata = subjdata[grouping]
        composite_data = groupdata['composite']
        hemiwm_data = groupdata['hemiWM']
        refdata = groupdata['wholecereb']
        pvc_ref = float(refdata['pvcval'])
        nonpvc_ref = float(refdata['nonpvcval'])

        composite_pvcval = composite_data['pvcval'] / pvc_ref
        composite_nonpvcval = composite_data['nonpvcval'] / nonpvc_ref
        hemiwm_pvcval = hemiwm_data['pvcval'] / pvc_ref
        hemiwm_nonpvcval = hemiwm_data['nonpvcval'] / nonpvc_ref

        color = 'b'
        if subj in NORMAL:
            color = 'g'
        elif subj in AD:
            color = 'r'
        else:
            continue

        # scatter plot
        pvc_points.append((hemiwm_pvcval, composite_pvcval, color))
        nonpvc_points.append((hemiwm_nonpvcval, composite_nonpvcval, color))

    pvc_x = [_[0] for _ in pvc_points]
    pvc_y = [_[1] for _ in pvc_points]
    pvc_color = [_[2] for _ in pvc_points]
    nonpvc_x = [_[0] for _ in nonpvc_points]
    nonpvc_y = [_[1] for _ in nonpvc_points]
    nonpvc_color = [_[2] for _ in nonpvc_points]

    min_x = min(pvc_x+nonpvc_x)
    max_x = max(pvc_x+nonpvc_x)
    min_y = min(pvc_y+nonpvc_y)
    max_y = max(pvc_y+nonpvc_y)
    min_x -= (max_x-min_x)*0.1
    max_x += (max_x-min_x)*0.1
    min_y -= (max_y-min_y)*0.1
    max_y += (max_y-min_y)*0.1

    plt.figure(1)
    plt.scatter(pvc_x, pvc_y, c=pvc_color)
    x1,x2,y1,y2 = plt.axis([min_x, max_x, min_y, max_y])
    plt.plot([x1,x2],[1.11,1.11])
    plt.xlabel('hemiWM')
    plt.ylabel('composite')

    plt.figure(2)
    ax = plt.scatter(nonpvc_x, nonpvc_y, c=nonpvc_color)
    x1,x2,y1,y2 = plt.axis([min_x, max_x, min_y, max_y])
    plt.plot([x1,x2],[1.11,1.11])
    plt.xlabel('hemiWM')
    plt.ylabel('composite')

if __name__ == "__main__":

    master_file = "../FDG_AV45_COGdata_08_21_15.csv"
    subj_data = importMaster(master_file)
    diag_key = 'Init_Diagnosis'
    by_diag = defaultdict(list)
    for k,v in subj_data.iteritems():
        new_key = v.get('Init_Diagnosis')
        if new_key:
            by_diag[new_key].append(k)

    AD = by_diag['AD']
    NORMAL = by_diag['N']
    MCI = by_diag['LMCI'] + by_diag['EMCI']


    # load agglomerative grouping results
    rousset_mat = '../rousset_output_low_high.mat'
    data = loadMATFile(rousset_mat)
    agg_result_list = [parseResult(_) for _ in data['result_list'][0]]
    agg_subj_list = data['subj_list'][0]
    agg_groups = list(data['group_names'][0])
    agg_sorted_data = dict(zip(agg_subj_list,agg_result_list))

    # load regular grouping results
    rousset_mat = '../rousset_output_2_4.mat'
    data = loadMATFile(rousset_mat)
    result_list = [parseResult(_) for _ in data['result_list'][0]]
    subj_list = data['subj_list'][0]
    groups = list(data['group_names'][0])
    sorted_data = dict(zip(subj_list,result_list))

    # aggregate
    for k,v in sorted_data.iteritems():
        sorted_data[k].update(agg_sorted_data.get(k,{}))


    
    '''
    plot the residuals for each group
    '''
    '''
    plt.figure(1)
    plt.title('Residuals')
    subplot_residuals(sorted_data)
    plt.show()
    sys.exit(1)
    '''

    '''
    Plot scatter plot (composite vs hemiwm)
    '''
    grouping='agghigh'
    composite_hemiWM_scatter(sorted_data, grouping)
    plt.show()
    sys.exit(1)


    '''
    plot pvc value distribution, by group
    '''
    grouping='agghigh'
    subplot_height = 1
    subplot_length = 5
    subplot_indices = range(1,(subplot_height*subplot_length)+1)
    subplot_order = ['wholecereb',
                     'composite',
                     'hemiWM',
                     'cerebWM',
                     'cerebGM']
    plt.figure(1)
    for subplot_i, subplot_key in zip(subplot_indices, subplot_order):
        ax = plt.subplot(subplot_height, subplot_length, subplot_i)
        plt.title(subplot_key)
        mean_graph(sorted_data, subplot_key, grouping)
    plt.show()
    sys.exit(1)
    

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

