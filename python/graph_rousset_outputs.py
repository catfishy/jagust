import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from collections import defaultdict
from scipy.stats import norm, mannwhitneyu, linregress
from sklearn.mixture import GMM, VBGMM, DPGMM
from utils import *
import itertools
from scipy import linalg
import matplotlib as mpl



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

def subplot_residuals_raw(groupname, data):
    residual_vectors = []
    colors = []
    assert groupname in ['group2', 'group4', 'agglow', 'agghigh']
    for subj, subjdata in data.iteritems():
        groupdata = subjdata[groupname]
        if subj in NORMAL:
            color = 'g'
        elif subj in AD:
            color = 'r'
        else:
            color = 'b'
        v = groupdata['residuals']['raw']
        residual_vectors.append(v)
        colors.append(color)
    scatter_points = []
    for i, group_res in enumerate(zip(*residual_vectors)):
        for c, val in zip(colors, group_res):
            scatter_points.append((i+1, val, c))
    '''
    avgs = [np.mean(by_group[i]) for i,g in enumerate(groups)]
    stds = [np.std(by_group[i]) for i,g in enumerate(groups)]
    plt.errorbar(range(1,len(groups)+1), avgs, stds)
    '''
    x = [_[0] for _ in scatter_points]
    y = [_[1] for _ in scatter_points]
    y_mod = abs(max(y)-min(y))*0.1
    color = [_[2] for _ in scatter_points]
    plt.scatter(x, y, c=color)
    x1,x2,y1,y2 = plt.axis()
    plt.axis((0,max(x)+1,min(y)-y_mod,max(y)+y_mod))
    plt.xlabel('Groups')
    plt.ylabel('PVC Residuals')

def composite_hemiWM_GMM(data, grouping, threshold=1.11):
    g = GMM(n_components=2, 
                covariance_type='full', 
                random_state=None, 
                thresh=None, 
                tol=0.00001, 
                min_covar=0.00001, 
                n_iter=200,
                params='wmc', 
                init_params='wmc')
    pvc_points = []
    for subj, subjdata in data.iteritems():
        groupdata = subjdata[grouping]
        composite_data = groupdata['composite']
        hemiwm_data = groupdata['hemiWM']
        pvc_ref = float(groupdata['wholecereb']['pvcval'])
        composite_pvcval = composite_data['pvcval'] / pvc_ref
        hemiwm_pvcval = hemiwm_data['pvcval'] / pvc_ref
        # scatter plot
        pvc_points.append(np.array([hemiwm_pvcval, composite_pvcval]))

    g.fit(pvc_points)
    print np.round(g.weights_, 2)
    print np.round(g.means_, 2)
    try:
        print np.round(g.precs_, 2) 
    except:
        print np.round(g.covars_, 2)
    print g.converged_

    # predict
    probs = g.predict_proba(pvc_points)
    colors = []
    pvc_above = []
    pvc_below = []
    for i, p in enumerate(probs):
        if p[0] > 0.75:
            colors.append('r')
            pvc_above.append(pvc_points[i])
        elif p[1] > 0.75:
            colors.append('g')
            pvc_below.append(pvc_points[i])
        else:
            colors.append('b')

    # plot
    pvc_x = [_[0] for _ in pvc_points]
    pvc_y = [_[1] for _ in pvc_points]

    min_x = min(pvc_x)
    max_x = max(pvc_x)
    min_y = min(pvc_y)
    max_y = max(pvc_y)
    min_x -= (max_x-min_x)*0.1
    max_x += (max_x-min_x)*0.1
    min_y -= (max_y-min_y)*0.1
    max_y += (max_y-min_y)*0.1

    plt.plot()
    splot = plt.subplot(1, 1, 1)
    plt.title('%s post-PVC' % grouping)
    plt.scatter(pvc_x, pvc_y, c=colors)
    x1,x2,y1,y2 = plt.axis([min_x, max_x, min_y, max_y])
    plt.plot([x1,x2],[threshold,threshold])
    slope, intercept, r, p, stderr = linregress([_[0] for _ in pvc_above], [_[1] for _ in pvc_above])
    plt.plot([min_x, max_x], [(min_x*slope)+intercept, (max_x*slope)+intercept])
    slope, intercept, r, p, stderr = linregress([_[0] for _ in pvc_below], [_[1] for _ in pvc_below])
    plt.plot([min_x, max_x], [(min_x*slope)+intercept, (max_x*slope)+intercept])


    color_iter = itertools.cycle(['r', 'g', 'b', 'c', 'm'])
    for i, (mean, covar, color) in enumerate(zip(g.means_, g._get_covars(), color_iter)):
        v, w = linalg.eigh(covar)
        u = w[0] / linalg.norm(w[0])
        # Plot an ellipse to show the Gaussian component
        angle = np.arctan(u[1] / u[0])
        angle = 180 * angle / np.pi  # convert to degrees
        ell = mpl.patches.Ellipse(mean, v[0], v[1], 180 + angle, color=color)
        ell.set_clip_box(splot.bbox)
        ell.set_alpha(0.5)
        splot.add_artist(ell)
    
    plt.xlabel('Hemi WM / wcereb')
    plt.ylabel('Cortical Summary / wcereb')
    plt.show()



def composite_hemiWM_scatter(data, grouping, threshold=1.11):
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

        if subj in NORMAL:
            color = 'g'
        elif subj in AD:
            color = 'r'
        else:
            color = 'b'

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

    # linearly fit points above/below thresholds
    pre_pvc_above = [_ for _ in nonpvc_points if _[1] >= 1.11]
    pre_pvc_below = [_ for _ in nonpvc_points if _[1] < 1.11]
    post_pvc_above = [_ for _ in pvc_points if _[1] >= threshold]
    post_pvc_below = [_ for _ in pvc_points if _[1] < threshold]


    plt.subplot(1, 2, 1)
    plt.title('%s post-PVC' % grouping)
    plt.scatter(pvc_x, pvc_y, c=pvc_color)
    x1,x2,y1,y2 = plt.axis([min_x, max_x, min_y, max_y])
    plt.plot([x1,x2],[threshold,threshold])
    slope, intercept, r, p, stderr = linregress([_[0] for _ in post_pvc_above], [_[1] for _ in post_pvc_above])
    plt.plot([min_x, max_x], [(min_x*slope)+intercept, (max_x*slope)+intercept])
    slope, intercept, r, p, stderr = linregress([_[0] for _ in post_pvc_below], [_[1] for _ in post_pvc_below])
    plt.plot([min_x, max_x], [(min_x*slope)+intercept, (max_x*slope)+intercept])
    plt.xlabel('Hemi WM / wcereb')
    plt.ylabel('Cortical Summary / wcereb')

    plt.subplot(1, 2, 2)
    plt.title('%s pre-PVC' % grouping)
    ax = plt.scatter(nonpvc_x, nonpvc_y, c=nonpvc_color)
    x1,x2,y1,y2 = plt.axis([min_x, max_x, min_y, max_y])
    plt.plot([x1,x2],[1.11,1.11])
    slope, intercept, r, p, stderr = linregress([_[0] for _ in pre_pvc_above], [_[1] for _ in pre_pvc_above])
    plt.plot([min_x, max_x], [(min_x*slope)+intercept, (max_x*slope)+intercept])
    slope, intercept, r, p, stderr = linregress([_[0] for _ in pre_pvc_below], [_[1] for _ in pre_pvc_below])
    plt.plot([min_x, max_x], [(min_x*slope)+intercept, (max_x*slope)+intercept])
    plt.xlabel('Hemi WM / wcereb')
    plt.ylabel('Cortical Summary / wcereb')

def thresholdCorrection(grouping, data):
    scatter_points = []

    for subj, subjdata in data.iteritems():
        groupdata = subjdata[grouping]
        composite_data = groupdata['composite']
        hemiwm_data = groupdata['hemiWM']
        refdata = groupdata['wholecereb']
        pvc_ref = float(refdata['pvcval'])
        nonpvc_ref = float(refdata['nonpvcval'])
        composite_pvcval = composite_data['pvcval'] / pvc_ref
        composite_nonpvcval = composite_data['nonpvcval'] / nonpvc_ref

        if subj in NORMAL:
            color = 'g'
        elif subj in AD:
            color = 'r'
        else:
            color = 'b'

        # scatter plot
        scatter_points.append((composite_nonpvcval, composite_pvcval, color))
    x = [_[0] for _ in scatter_points]
    y = [_[1] for _ in scatter_points]
    color = [_[2] for _ in scatter_points]
    y_mod = abs(max(y)-min(y))*0.1
    x_mod = abs(max(x)-min(x))*0.1
    min_x,max_x,min_y,max_y = (min(x)-x_mod,max(x)+x_mod,min(y)-y_mod,max(y)+y_mod)

    # linear fit
    slope, intercept, r, p, stderr = linregress(x, y)
    low_fit_point = min_x*slope + intercept
    high_fit_point = max_x*slope + intercept
    threshold_fit_point = 1.11*slope + intercept

    print "%s Corrected Threshold: %s" % (grouping, threshold_fit_point)

    plt.scatter(x, y, c=color)
    x1,x2,y1,y2 = plt.axis()
    plt.axis((min_x,max_x,min_y,max_y))
    plt.xlabel('NonPVC')
    plt.ylabel('PVC')
    plt.plot([1.11, 1.11],[min_y,max_y])
    plt.plot([min_x, max_x],[low_fit_point, high_fit_point])
    plt.plot([min_x, max_x],[threshold_fit_point, threshold_fit_point])



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
    (agg_group_names, agg_sorted_data) = importRoussetResults(rousset_mat)
    # load regular grouping results
    rousset_mat = '../rousset_output_2_4.mat'
    (group_names, sorted_data) = importRoussetResults(rousset_mat)

    # aggregate
    for k,v in sorted_data.iteritems():
        sorted_data[k].update(agg_sorted_data.get(k,{}))

    groupings = ['group2', 'group4', 'agglow', 'agghigh', 'agglowtwo']
    thresholds = {'group2': 1.05225257559,
                  'group4': 1.10017792162,
                  'agglow': 1.26905730506,
                  'agghigh': 1.26780842565,
                  'agglowtwo': 1.27033990344}

    '''
    plot the nonpvc to pvc linear correction
    '''
    '''
    for i, groupname in enumerate(groupings):
        fig_num = i+1
        plt.figure(fig_num)
        thresholdCorrection(groupname, sorted_data)
    plt.show()
    sys.exit(1)
    '''
    
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
    plot the raw residuals 
    '''
    '''
    plt.figure(1)
    groupname = 'agghigh'
    subplot_residuals_raw(groupname, sorted_data)
    plt.show()
    sys.exit(1)
    '''

    '''
    GMM to composite/hemiWM scatter plot
    '''
    grouping='agghigh'
    composite_hemiWM_GMM(sorted_data, grouping, threshold=thresholds[grouping])
    sys.exit(1)

    '''
    Plot scatter plot (composite vs hemiwm)
    '''

    for i, grouping in enumerate(groupings):
        fig_num = i+1
        plt.figure(fig_num)
        composite_hemiWM_scatter(sorted_data, grouping, threshold=thresholds[grouping])
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

