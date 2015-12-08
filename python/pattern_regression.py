import pandas as pd
import numpy as np
import random
from collections import Counter, defaultdict
import itertools
import multiprocessing as mp

from scipy.stats import f_oneway, norm
import matplotlib.pyplot as plt
from matplotlib import cm
from pylab import get_cmap
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from sklearn.mixture import GMM, VBGMM, DPGMM
from sklearn.decomposition import PCA, KernelPCA
from sklearn.cross_validation import LeaveOneOut
from sklearn.mixture import GMM, VBGMM, DPGMM
from sklearn.tree import DecisionTreeRegressor
from sklearn.preprocessing import StandardScaler
import networkx as nx


from utils import saveFakeAparcInput, importFreesurferLookup, importMaster

pd.options.display.mpl_style = 'default'

# extract patterns and summary value, summary change
pattern_keys = ['BRAIN_STEM', 'CTX_LH_BANKSSTS', 'CTX_LH_CAUDALANTERIORCINGULATE', 'CTX_LH_CAUDALMIDDLEFRONTAL', 'CTX_LH_CUNEUS', 'CTX_LH_ENTORHINAL', 'CTX_LH_FRONTALPOLE', 'CTX_LH_FUSIFORM', 'CTX_LH_INFERIORPARIETAL', 'CTX_LH_INFERIORTEMPORAL', 'CTX_LH_INSULA', 'CTX_LH_ISTHMUSCINGULATE', 'CTX_LH_LATERALOCCIPITAL', 'CTX_LH_LATERALORBITOFRONTAL', 'CTX_LH_LINGUAL', 'CTX_LH_MEDIALORBITOFRONTAL', 'CTX_LH_MIDDLETEMPORAL', 'CTX_LH_PARACENTRAL', 'CTX_LH_PARAHIPPOCAMPAL', 'CTX_LH_PARSOPERCULARIS', 'CTX_LH_PARSORBITALIS', 'CTX_LH_PARSTRIANGULARIS', 'CTX_LH_PERICALCARINE', 'CTX_LH_POSTCENTRAL', 'CTX_LH_POSTERIORCINGULATE', 'CTX_LH_PRECENTRAL', 'CTX_LH_PRECUNEUS', 'CTX_LH_ROSTRALANTERIORCINGULATE', 'CTX_LH_ROSTRALMIDDLEFRONTAL', 'CTX_LH_SUPERIORFRONTAL', 'CTX_LH_SUPERIORPARIETAL', 'CTX_LH_SUPERIORTEMPORAL', 'CTX_LH_SUPRAMARGINAL', 'CTX_LH_TEMPORALPOLE', 'CTX_LH_TRANSVERSETEMPORAL', 'CTX_LH_UNKNOWN', 'CTX_RH_BANKSSTS', 'CTX_RH_CAUDALANTERIORCINGULATE', 'CTX_RH_CAUDALMIDDLEFRONTAL', 'CTX_RH_CUNEUS', 'CTX_RH_ENTORHINAL', 'CTX_RH_FRONTALPOLE', 'CTX_RH_FUSIFORM', 'CTX_RH_INFERIORPARIETAL', 'CTX_RH_INFERIORTEMPORAL', 'CTX_RH_INSULA', 'CTX_RH_ISTHMUSCINGULATE', 'CTX_RH_LATERALOCCIPITAL', 'CTX_RH_LATERALORBITOFRONTAL', 'CTX_RH_LINGUAL', 'CTX_RH_MEDIALORBITOFRONTAL', 'CTX_RH_MIDDLETEMPORAL', 'CTX_RH_PARACENTRAL', 'CTX_RH_PARAHIPPOCAMPAL', 'CTX_RH_PARSOPERCULARIS', 'CTX_RH_PARSORBITALIS', 'CTX_RH_PARSTRIANGULARIS', 'CTX_RH_PERICALCARINE', 'CTX_RH_POSTCENTRAL', 'CTX_RH_POSTERIORCINGULATE', 'CTX_RH_PRECENTRAL', 'CTX_RH_PRECUNEUS', 'CTX_RH_ROSTRALANTERIORCINGULATE', 'CTX_RH_ROSTRALMIDDLEFRONTAL', 'CTX_RH_SUPERIORFRONTAL', 'CTX_RH_SUPERIORPARIETAL', 'CTX_RH_SUPERIORTEMPORAL', 'CTX_RH_SUPRAMARGINAL', 'CTX_RH_TEMPORALPOLE', 'CTX_RH_TRANSVERSETEMPORAL', 'CTX_RH_UNKNOWN', 'LEFT_ACCUMBENS_AREA', 'LEFT_AMYGDALA', 'LEFT_CAUDATE', 'LEFT_CEREBELLUM_CORTEX', 'LEFT_CEREBELLUM_WHITE_MATTER', 'LEFT_CEREBRAL_WHITE_MATTER', 'LEFT_HIPPOCAMPUS', 'LEFT_PALLIDUM', 'LEFT_PUTAMEN', 'LEFT_THALAMUS_PROPER', 'LEFT_VENTRALDC', 'RIGHT_ACCUMBENS_AREA', 'RIGHT_AMYGDALA', 'RIGHT_CAUDATE', 'RIGHT_CEREBELLUM_CORTEX', 'RIGHT_CEREBELLUM_WHITE_MATTER', 'RIGHT_CEREBRAL_WHITE_MATTER', 'RIGHT_HIPPOCAMPUS', 'RIGHT_PALLIDUM', 'RIGHT_PUTAMEN', 'RIGHT_THALAMUS_PROPER', 'RIGHT_VENTRALDC']
prior_keys = ['%s_prior' % _ for _ in pattern_keys]
post_keys = ['%s_post' % _ for _ in pattern_keys]
rchange_keys = ['%s_change' % _ for _ in pattern_keys]
result_keys = ['CORTICAL_SUMMARY_post', 'CORTICAL_SUMMARY_prior', 'CORTICAL_SUMMARY_change', 'diag_prior', 'diag_post']


def importCogSlopes():
    master_file = '../FDG_AV45_COGdata_11_10_15.csv'
    master_data = importMaster(master_file)


def sample(data):
    return data[np.random.randint(0,len(data),(1,len(data)))[0]]

def mp_bootstrap_wrapper(arg):
    key, group1_vals, group2_vals = arg
    group1_boot = bootstrap_mean(group1_vals)
    group2_boot = bootstrap_mean(group2_vals)
    pvalue = bootstrap_t_test(group1_vals, group2_vals)
    return (key, pvalue, group1_boot, group2_boot)

def bootstrap_t_test(treatment, control, nboot=100000):
    treatment = np.array(treatment)
    control = np.array(control)
    treatment_len = len(treatment)
    tstat = abs(np.mean(treatment)-np.mean(control))
    # adjust mean
    treatment = treatment - np.mean(treatment) + np.mean(control)
    # concatenate and scramble
    Z = np.concatenate((treatment,control))
    np.random.shuffle(Z)
    tboot = np.zeros(nboot)
    idx = np.random.randint(0,len(Z),(nboot,len(Z)))
    for i in xrange(nboot):
        sboot = Z[idx[i]]
        tboot[i] = np.mean(sboot[:treatment_len]) - np.mean(sboot[treatment_len:])
    pvalue = np.sum(abs(tboot)>=tstat) / float(nboot)
    return pvalue

def bootstrap_mean(vals, nboot=1000):
    vals = np.array(vals)
    tboot = np.zeros(nboot)
    for i in xrange(nboot):
        sboot = sample(vals)
        tboot[i] = np.mean(sboot)
    return np.mean(tboot)

def chooseBestModel(patterns, components, alpha, covar_type):
    # Cluster
    models_scores = {}
    for i in range(30):
        g = DPGMM(n_components=components, 
                  covariance_type=covar_type,
                  alpha=alpha,
                  tol=1e-6,
                  n_iter=700,
                  params='wmc', 
                  init_params='wmc',
                  verbose=False)
        print "Fitting: %s" % i
        g.fit(patterns)
        print g.converged_
        y_ = g.predict(patterns)
        print Counter(y_)
        membership = np.zeros((len(patterns.index),components))
        for i,member in enumerate(y_):
            membership[i][member] = 1
        bound = g.lower_bound(patterns, membership)
        score = np.mean(g.score(patterns))
        models_scores[bound] = g
    best_score = max(models_scores.keys())
    best_model = models_scores[best_score]
    return best_model

def gmm_sweep_alpha(components, data, covar_type='full', alphas=None, graph=False):
    # model selection for alpha by looking at lower bound on marginal likelihood, or log prob of X under model
    if alphas is None:
        alphas = np.linspace(0.01,5,100)
    alpha_to_clusters = {}
    points = []
    for a in alphas:
        num_clusters = []
        var_clusters = []
        scores = []
        print a
        for i in range(3):
            g = DPGMM(n_components=components, 
                      covariance_type=covar_type,
                      alpha=a,
                      tol=1e-5,
                      n_iter=700,
                      params='wmc', 
                      init_params='wmc',
                      verbose=False)
            g.fit(data)
            print g.converged_
            y_ = g.predict(data)
            counter = dict(Counter(y_))
            score = np.mean(g.score(data))
            scores.append(score)
            num_clusters.append(len(counter.keys()))
            var_clusters.append(np.var(counter.values()))
        print scores
        points.append((a,np.mean(num_clusters),np.mean(var_clusters)))
        print points
        alpha_to_clusters[a] = np.mean(scores)
    if graph:
        plt.figure(1)
        alphas = [_[0] for _ in points]
        counts = [_[1] for _ in points]
        stds = [np.sqrt(_[2]) for _ in points]
        scores = [alpha_to_clusters[_] for _ in alphas]
        plt.plot(alphas, counts, label='counts')
        plt.plot(alphas, stds, label='stds')
        plt.plot(alphas, scores, label='scores')
        plt.legend()
        plt.show()
    return alpha_to_clusters

def compareTwoGroups(group1, group2, df, keys):
    group1_members = df[df['membership']==group1]
    group2_members = df[df['membership']==group2]
    print "GRP %s (%s) vs GRP %s (%s)" % (group1, len(group1_members), group2, len(group2_members))
    # bonferri correct
    hypos = len(keys)
    level = 0.05 / float(hypos)
    print level
    results = {}

    def generateArgs(group1_members, group2_members, keys):
        for i, key in enumerate(keys):
            group1_vals = group1_members[key].tolist()
            group2_vals = group2_members[key].tolist()
            yield (key, group1_vals, group2_vals)

    pool = mp.Pool(processes=4)
    for response in pool.imap(mp_bootstrap_wrapper, generateArgs(group1_members, group2_members, keys), 1):
        key, pvalue, group1_boot, group2_boot = response
        results[key] = (pvalue, group1_boot, group2_boot)
        if pvalue <= level:
            print "\t%s: %s (%s, %s)" % (key, pvalue, group1_boot, group2_boot)
    return results

def group_comparisons(df, groups, keys):
    data = {}
    for group1, group2 in itertools.combinations(groups,2):
        data[(group1,group2)] = compareTwoGroups(group1, group2, df, keys)
    return data

def parseRawInput(patterns_csv, ref_key):
    # read in data
    raw_df = pd.read_csv(patterns_csv)
    raw_df.set_index('rid',inplace=True,drop=True)
    # convert to SUVR units and calculate change
    priors = [_ for _ in raw_df.columns if '_prior' in _ and 'diag' not in _]
    prior_ref = raw_df.loc[:,'%s_prior' % ref_key]
    raw_df[priors] = raw_df[priors].divide(prior_ref,axis='index')
    posts = [_ for _ in raw_df.columns if '_post' in _ and 'diag' not in _]
    post_ref = raw_df.loc[:,'%s_post' % ref_key]
    raw_df[posts] = raw_df[posts].divide(post_ref,axis='index')
    yrs = raw_df['yrs']
    for prior_key in priors:
        post_key = prior_key.replace('_prior','_post')
        change_key = prior_key.replace('_prior','_change')
        raw_df[change_key] = (raw_df[post_key]-raw_df[prior_key]).divide(yrs,axis='index')

    # convert to patterns
    pattern_prior_df = raw_df[prior_keys].copy()
    pattern_prior_df.columns = [_.replace('_prior','') for _ in pattern_prior_df.columns]
    pattern_prior_df = pattern_prior_df.divide(pattern_prior_df.sum(axis=1),axis=0)
    pattern_prior_df['timepoint'] = 'prior'
    pattern_prior_df.set_index(raw_df.index, inplace=True)
    pattern_prior_df.set_index('timepoint', append=True, inplace=True)
    pattern_prior_df.dropna(inplace=True)

    pattern_post_df = raw_df[post_keys].copy()
    pattern_post_df.columns = [_.replace('_post','') for _ in pattern_post_df.columns]
    pattern_post_df = pattern_post_df.divide(pattern_post_df.sum(axis=1),axis=0)
    pattern_post_df['timepoint'] = 'post'
    pattern_post_df.set_index(raw_df.index, inplace=True)
    pattern_post_df.set_index('timepoint', append=True, inplace=True)
    pattern_post_df.dropna(inplace=True)

    #pattern_df = pd.concat((pattern_prior_df,pattern_post_df))
    
    uptake_prior_df = raw_df[prior_keys]
    uptake_prior_df.columns = [_.replace('_prior','') for _ in uptake_prior_df.columns]
    uptake_prior_df.dropna(inplace=True)

    uptake_post_df = raw_df[post_keys]
    uptake_post_df.columns = [_.replace('_post','') for _ in uptake_post_df.columns]
    uptake_post_df.dropna(inplace=True)

    result_df = raw_df[result_keys]
    rchange_df = raw_df[rchange_keys]

    return (pattern_prior_df, pattern_post_df, uptake_prior_df, uptake_post_df, result_df, rchange_df)

def loadResults(results_file):
    results = pd.read_csv(results_file)
    results.set_index('rid',inplace=True,drop=True)
    return results

def bigGroups(result_df, threshold=10):
    # determine big groups
    groups = np.array(list(set(result_df.membership_prior) | set(result_df.membership_post)))
    groups = groups[~np.isnan(groups)]
    big_groups = []
    for g in groups:
        prior_members = set(result_df[result_df.membership_prior==g].index)
        post_members = set(result_df[result_df.membership_post==g].index)
        allmembers = len(prior_members | post_members)
        print "%s: %s" % (g,allmembers)
        if allmembers >= threshold:
            big_groups.append(g)
    return big_groups

def fitNormalCdf(values, threshold):
    g_1 = GMM(n_components=1, 
              covariance_type='full',
              tol=1e-6,
              n_iter=700,
              params='wmc', 
              init_params='wmc')
    g_1.fit([[_] for _ in values])
    mu_1, sigma_1 = (g_1.means_[0][0],np.sqrt(g_1.covars_[0][0][0]))
    cdf_val = norm.cdf(threshold,loc=mu_1,scale=sigma_1)
    return cdf_val

def parseConversions(groups, result_df, threshold, master_cols):
    '''
    'APOE4_BIN'
    'UW_MEM_slope'
    'UW_EF_slope'
    'Gender'
    'Handedness'
    'WMH_slope'
    '''
    # group diagnostic/pattern/positivity conversions
    diags = ['N', 'SMC', 'MCI', 'AD']
    conversions = {}
    for g in groups:
        cur_conversions = {}
        conversions[g] = {}
        members_prior = result_df[result_df.membership_prior==g]
        members_prior_with_post = members_prior.dropna(axis=0,subset=['membership_post'])
        members_post = result_df[result_df.membership_post==g]
        count_prior = float(len(members_prior.index))
        count_prior_with_post = float(len(members_prior_with_post.index))
        count_post = float(len(members_post.index))
        cur_conversions['count_prior'] = int(count_prior)
        cur_conversions['count_post'] = int(count_post)
        members_prior_rid = list(members_prior.index)
        master_rows = result_df.loc[members_prior_rid,master_cols]
        master_avgs = dict(master_rows.mean())
        cortical_summary_values = np.array(list(members_prior['CORTICAL_SUMMARY_prior']) + list(members_post['CORTICAL_SUMMARY_post']))
        under_threshold = fitNormalCdf(cortical_summary_values, threshold)
        cur_conversions.update(master_avgs)
        cur_conversions['cortical_summary_change'] = members_prior['CORTICAL_SUMMARY_change'].mean()
        cur_conversions['cortical_summary'] = cortical_summary_values.mean()
        cur_conversions['below_threshold'] = under_threshold
        if int(g) == 67:
            print members_prior[['membership_prior','membership_post','diag_prior','diag_post']]
            print members_post[['membership_prior','membership_post','diag_prior','diag_post']]
        # diag counts:
        prior_counts = Counter(list(members_prior['diag_prior']))
        post_counts = Counter(list(members_post['diag_post']))
        diag_counts = {}
        for diag in diags:
            if diag == 'MCI':
                count = prior_counts.get("LMCI",0) + post_counts.get("LMCI",0) + prior_counts.get("EMCI",0) + post_counts.get("EMCI",0)
            else:
                count = prior_counts.get(diag,0) + post_counts.get(diag,0)
            diag_counts[diag] = count
        total_counts = float(sum(diag_counts.values()))
        diag_counts = {k:v/total_counts for k,v in diag_counts.iteritems()}
        cur_conversions.update(diag_counts)
        # diagnostic conversions
        for diag1, diag2 in itertools.product(diags,repeat=2):
            if diag1 == 'MCI':
                prior_members = set(members_prior[members_prior.diag_prior=='EMCI'].index) | set(members_prior[members_prior.diag_prior=='LMCI'].index)
            else:
                prior_members = set(members_prior[members_prior.diag_prior==diag1].index)
            if diag2 == 'MCI':
                post_members = set(members_prior[members_prior.diag_post=='EMCI'].index) | set(members_prior[members_prior.diag_post=='LMCI'].index)
            else:
                post_members = set(members_prior[members_prior.diag_post==diag2].index)
            cur_conversions['%s-%s' % (diag1, diag2)] = len(prior_members & post_members)/count_prior_with_post
        # positivity conversions
        cur_conversions['neg-neg'] = len(members_prior[(members_prior.CORTICAL_SUMMARY_prior<threshold) & (members_prior.CORTICAL_SUMMARY_post<threshold)].index)/count_prior_with_post
        cur_conversions['neg-pos'] = len(members_prior[(members_prior.CORTICAL_SUMMARY_prior<threshold) & (members_prior.CORTICAL_SUMMARY_post>=threshold)].index)/count_prior_with_post
        cur_conversions['pos-pos'] = len(members_prior[(members_prior.CORTICAL_SUMMARY_prior>=threshold) & (members_prior.CORTICAL_SUMMARY_post>=threshold)].index)/count_prior_with_post
        cur_conversions['pos-neg'] = len(members_prior[(members_prior.CORTICAL_SUMMARY_prior>=threshold) & (members_prior.CORTICAL_SUMMARY_post<threshold)].index)/count_prior_with_post
        # pattern conversions
        for g2 in groups:
            cur_conversions[g2] = len(members_prior[members_prior.membership_post==g2].index)/count_prior_with_post
        conversions[g] = cur_conversions
    conversions = pd.DataFrame(conversions).T
    conversions.index.name = 'pattern'
    return conversions

def saveAparcs(alpha, components, groups, pattern_members, uptake_members, change_members):
    # save group patterns and change as fake aparc inputs (for visualizing)
    lut_file = "../FreeSurferColorLUT.txt"
    lut_table = importFreesurferLookup(lut_file)
    index_lookup = {v.replace('-','_').upper():[k] for k,v in lut_table.iteritems()}
    aparc_input_template = "../output/fake_aparc_inputs/dpgmm_alpha%s_comp%s_group_%s_%s"
    for g in groups:
        out_file_pattern = aparc_input_template % (alpha, components, g, 'pattern')
        out_file_change = aparc_input_template % (alpha, components, g, 'change')
        out_file_uptake = aparc_input_template % (alpha, components, g, 'uptake')
        pattern_values = pattern_members[pattern_members['membership']==g][pattern_keys].mean().to_dict()
        uptake_values = uptake_members[uptake_members['membership']==g][pattern_keys].mean().to_dict()
        change_values = change_members[change_members['membership']==g][rchange_keys].mean().to_dict()
        change_values = {k.replace('_change',''):v for k,v in change_values.iteritems()}
        saveFakeAparcInput(out_file_pattern, pattern_values, index_lookup)
        saveFakeAparcInput(out_file_change, change_values, index_lookup)
        saveFakeAparcInput(out_file_uptake, uptake_values, index_lookup)

def mergeResults(result_df, pattern_prior_df, pattern_post_df, uptake_prior_df, uptake_post_df, rchange_df):
    # merge results with regional change and uptakes
    change_members = rchange_df.merge(result_df, left_index=True, right_index=True)
    change_members['membership'] = change_members['membership_prior']
    uptake_members_prior = uptake_prior_df.merge(result_df[['membership_prior','CORTICAL_SUMMARY_prior']], left_index=True, right_index=True)
    uptake_members_prior['membership'] = uptake_members_prior['membership_prior']
    uptake_members_prior['CORTICAL_SUMMARY'] = uptake_members_prior['CORTICAL_SUMMARY_prior']
    uptake_members_post = uptake_post_df.merge(result_df[['membership_post','CORTICAL_SUMMARY_post']], left_index=True, right_index=True)
    uptake_members_post['membership'] = uptake_members_post['membership_post']
    uptake_members_post['CORTICAL_SUMMARY'] = uptake_members_post['CORTICAL_SUMMARY_post']
    uptake_members = pd.concat((uptake_members_prior,uptake_members_post))
    pattern_prior_df = pattern_prior_df.xs('prior',level='timepoint')
    pattern_post_df = pattern_post_df.xs('post',level='timepoint')
    pattern_members_prior = pattern_prior_df.merge(result_df[['membership_prior','CORTICAL_SUMMARY_prior']], left_index=True, right_index=True)
    pattern_members_post = pattern_post_df.merge(result_df[['membership_post','CORTICAL_SUMMARY_post']], left_index=True, right_index=True)
    pattern_members_prior['membership'] = pattern_members_prior['membership_prior']
    pattern_members_prior['CORTICAL_SUMMARY'] = pattern_members_prior['CORTICAL_SUMMARY_prior']
    pattern_members_post['membership'] = pattern_members_post['membership_post']
    pattern_members_post['CORTICAL_SUMMARY'] = pattern_members_post['CORTICAL_SUMMARY_post']
    pattern_members = pd.concat((pattern_members_prior,pattern_members_post))
    return (uptake_members, pattern_members, change_members)

def scaleRawInput(pattern_df, scale_type='original'):
    assert scale_type in set(['original', 'pca', 'kernelpca'])
    if scale_type == 'original':
        # Scale in original space
        scaler = StandardScaler().fit(pattern_df)
        pattern_df_raw = pd.DataFrame(scaler.transform(pattern_df))
        pattern_df_raw.set_index(pattern_df.index, inplace=True)
        raw_keys = pattern_df_raw.columns
        return (pattern_df_raw[raw_keys], scaler)
    elif scale_type == 'pca':
        # convert to PCA space + whiten(scale)
        pca_model = PCA(n_components=len(pattern_df.columns), copy=True, whiten=True)
        pattern_df_pca = pd.DataFrame(pca_model.fit_transform(pattern_df))
        pattern_df_pca.set_index(pattern_df.index, inplace=True)
        pca_keys = pattern_df_pca.columns
        return (pattern_df_pca[pca_keys], scaler)
    elif scale_type == 'kernelpca':
        # convert to Kernel PCA Projection and scale
        kpca_model = KernelPCA(n_components=len(pattern_df.columns),kernel='rbf',fit_inverse_transform=False,gamma=1,alpha=1.0)
        X_kpca = kpca_model.fit_transform(pattern_df)
        pattern_df_kpca = pd.DataFrame(X_kpca)
        scaler = StandardScaler().fit(pattern_df_kpca)
        pattern_df_kpca = pd.DataFrame(scaler.transform(pattern_df_kpca))
        pattern_df_kpca.set_index(pattern_df.index, inplace=True)
        kpca_keys = pattern_df_kpca.columns
        return (pattern_df_kpca[kpca_keys], scaler)

def graphNetworkConversions(groups, conversions, iterations=50):
    G = nx.DiGraph()
    cmap = get_cmap("cool")
    for g1 in groups:
        for g2 in groups:
            if g1 == g2:
                if g1 not in G.nodes():
                    G.add_node(int(g1))
                continue
            conv_percent = conversions.loc[g1,g2]
            if conv_percent > 0.1:
                G.add_edge(int(g1),int(g2),weight=conv_percent)
    # edges sizes are % of pattern members converted to new pattern
    edgewidth = []
    edgelabels = {}
    for (u,v,d) in G.edges(data=True):
        edgewidth.append(d['weight']*20)
        edgelabels[(u,v)] = round(d['weight'],2)
    norm_cmap=ScalarMappable(norm=Normalize(min(edgewidth),max(edgewidth)),cmap=cmap)
    edgecolors = norm_cmap.to_rgba(edgewidth)
    # node sizes are number of prior+post members in group
    sizes = {}
    for n in G.nodes():
        sizes[n] = conversions.loc[n,'count_prior'] + conversions.loc[n,'count_post']
    print sizes
    nodesize = [sizes[_]*80 for _ in G]
    pos = nx.spring_layout(G,iterations=iterations)
    #pos = nx.circular_layout(G)
    plt.figure(figsize=(10,10))
    #nx.draw_networkx_edges(G,pos,width=edgewidth,edge_color='m',alpha=0.3)
    nodes = nx.draw_networkx_nodes(G,pos,node_size=nodesize,node_color='w',alpha=0.4)
    nodes.set_edgecolor('k')
    nx.draw_networkx_edges(G,pos,alpha=0.25,node_size=0,width=4,edge_color=edgecolors)
    #nx.draw_networkx_edge_labels(G,pos,edge_labels=edgelabels)
    nx.draw_networkx_labels(G,pos,font_size=14,font_family='sans-serif')
    plt.axis('off')
    plt.show()



if __name__ == '__main__':
    ref_key = 'COMPOSITE_REF'
    patterns_csv = '../datasets/pvc_allregions_uptake_change.csv'
    threshold = 0.91711
    pattern_prior_df, pattern_post_df, uptake_prior_df, uptake_post_df, result_df, rchange_df = parseRawInput(patterns_csv, ref_key)

    # import master
    master_cols = ['APOE4_BIN','UW_MEM_slope','UW_EF_slope','Gender','Handedness','WMH_percentOfICV_slope','CSF_TAU_slope','FSX_HC/ICV_slope'] # ADD HIPPOCAMPAL VOLUME CHANGE, CSF TAU
    master_csv = '../FDG_AV45_COGdata_12_03_15.csv'
    master_df = pd.read_csv(master_csv, low_memory=False, header=[0,1])
    master_df.columns = master_df.columns.get_level_values(1)
    master_df.set_index('RID', inplace=True)
    master_df = master_df.loc[:,master_cols]

    # # Scale inputs
    # patterns_only, scaler = scaleRawInput(pattern_prior_df, scale_type=scale_type)
    # post_patterns_only = pd.DataFrame(scaler.transform(pattern_post_df))
    # post_patterns_only.set_index(pattern_post_df.index, inplace=True)

    # # Sweep for alpha
    # scale_type = 'original'
    # components = 100
    # patterns_only, scaler = scaleRawInput(pattern_prior_df, scale_type=scale_type)
    # alpha_to_clusters = gmm_sweep_alpha(100, patterns_only, covar_type='spherical')
    # print sorted(alpha_to_clusters.items(), key=lambda x:x[1], reverse=True)

    # # Choose alpha + do model selection
    # alpha = 1.87495
    # best_model = chooseBestModel(patterns_only, components, alpha, 'spherical')

    # # Predict pattern groups and add to result df
    # y_df = pd.DataFrame(best_model.predict(patterns_only))
    # y_df.columns = ['group']
    # y_df.set_index(patterns_only.index, inplace=True)
    # y_df_post = pd.DataFrame(best_model.predict(post_patterns_only))
    # y_df_post.columns = ['group']
    # y_df_post.set_index(post_patterns_only.index, inplace=True)
    # for rid in result_df.index:
    #     result_df.loc[rid,'membership_prior'] = y_df.loc[(rid,'prior'),'group']
    #     if (rid, 'post') in y_df_post.index:
    #         result_df.loc[rid,'membership_post'] = y_df_post.loc[(rid,'post'),'group']

    # # save result
    # result_df.to_csv('../dpgmm_alpha1.87_spherical_result.csv')

    # sys.exit(1)

    # load result instead
    result_df = loadResults('../dpgmm_alpha1.87_spherical_result.csv')

    # generate conversion data
    result_df = result_df.merge(master_df,left_index=True,right_index=True)
    big_groups = bigGroups(result_df)
    conversions = parseConversions(big_groups, result_df, threshold, master_cols)
    #conversions.to_csv('../dpgmm_alpha%s_spherical_conversions.csv' % (round(alpha,2)))
    positive_patterns = list(conversions[conversions['cortical_summary']>=threshold].index)
    negative_patterns = list(conversions[conversions['neg-neg']>=0.8].index)
    transition_patterns = list(set(conversions.index) - (set(positive_patterns) | set(negative_patterns)))
    uptake_members, pattern_members, change_members = mergeResults(result_df, pattern_prior_df, pattern_post_df, uptake_prior_df, uptake_post_df, rchange_df)
    # saveAparcs(round(alpha,2), components, big_groups, pattern_members, uptake_members, change_members)
    # graphNetworkConversions(positive_patterns, conversions, iterations=50)

    # for explicitly comparing two patterns
    compareTwoGroups(34, 7, uptake_members, pattern_keys)
    compareTwoGroups(34, 7, pattern_members, pattern_keys)
    compareTwoGroups(34, 7, change_members, change_keys)

    # between group bootstrap hypothesis test for summary/regional change
    change_keys = rchange_keys
    change_pvalues = group_comparisons(change_members, positive_patterns, change_keys)
    flattened = []
    for k,v in change_pvalues.iteritems():
        row = {}
        grp1, grp2 = k
        row['GROUP1'] = grp1
        row['GROUP2'] = grp2
        for k_reg, v_reg in v.iteritems():
            p_reg, m1_reg, m2_reg = v_reg
            row['%s_%s' % ('PVALUE', k_reg)] = p_reg
            row['%s_%s' % ('MEAN1', k_reg)] = m1_reg
            row['%s_%s' % ('MEAN2', k_reg)] = m2_reg
        flattened.append(row)
    flattened_change_df = pd.DataFrame(flattened)
    flattened_change_df.columns = [_.replace('_change','') for _ in flattened_change_df.columns]
    flattened_change_df.loc[:,'TYPE'] = 'change'
    flattened_change_df.set_index(['GROUP1','GROUP2','TYPE'],inplace=True)

    # between group bootstrap hypothesis test for summary/regional uptake
    uptake_keys = pattern_keys
    uptake_pvalues = group_comparisons(uptake_members, positive_patterns, pattern_keys)
    flattened = []
    for k,v in uptake_pvalues.iteritems():
        row = {}
        grp1, grp2 = k
        row['GROUP1'] = grp1
        row['GROUP2'] = grp2
        for k_reg, v_reg in v.iteritems():
            p_reg, m1_reg, m2_reg = v_reg
            row['%s_%s' % ('PVALUE', k_reg)] = p_reg
            row['%s_%s' % ('MEAN1', k_reg)] = m1_reg
            row['%s_%s' % ('MEAN2', k_reg)] = m2_reg
        flattened.append(row)
    flattened_uptake_df = pd.DataFrame(flattened)
    flattened_uptake_df.loc[:,'TYPE'] = 'uptake'
    flattened_uptake_df.set_index(['GROUP1','GROUP2','TYPE'],inplace=True)

    # between group bootstrap hypothesis test for summary/regional patterns
    pattern_pvalues = group_comparisons(pattern_members, big_groups, pattern_keys)
    flattened = []
    for k,v in pattern_pvalues.iteritems():
        row = {}
        grp1, grp2 = k
        row['GROUP1'] = grp1
        row['GROUP2'] = grp2
        for k_reg, v_reg in v.iteritems():
            p_reg, m1_reg, m2_reg = v_reg
            row['%s_%s' % ('PVALUE', k_reg)] = p_reg
            row['%s_%s' % ('MEAN1', k_reg)] = m1_reg
            row['%s_%s' % ('MEAN2', k_reg)] = m2_reg
        flattened.append(row)
    flattened_pattern_df = pd.DataFrame(flattened)
    flattened_pattern_df.loc[:,'TYPE'] = 'pattern'
    flattened_pattern_df.set_index(['GROUP1','GROUP2','TYPE'],inplace=True)

    change_pvalues = flattened_change_df[[_ for _ in flattened_change_df.columns if 'PVALUE' in _]]
    uptake_pvalues = flattened_uptake_df[[_ for _ in flattened_uptake_df.columns if 'PVALUE' in _]]
    pattern_pvalues = flattened_pattern_df[[_ for _ in flattened_pattern_df.columns if 'PVALUE' in _]]


    # aggregate effects
    effects = []
    keys = set(change_pvalues.keys()) & set(uptake_pvalues.keys())
    for k in keys:
        group1, group2 = k
        row = {'group1': group1,
               'group2': group2}
        uptake_res = uptake_pvalues[k]
        change_res = change_pvalues[k]
        for ck in uptake_keys:
            uptake_pvalue = uptake_res[ck][0]
            group1_uptake = uptake_res[ck][1]
            group2_uptake = uptake_res[ck][2]
            cur_row = {'group1_%s_uptake' % ck: group1_uptake,
                       'group2_%s_uptake' % ck: group2_uptake,
                       'pvalue_%s_uptake' % ck: uptake_pvalue}
            row.update(cur_row)
        for ck in change_keys:
            change_pvalue = change_res[ck][0]
            group1_change = change_res[ck][1]
            group2_change = change_res[ck][2]
            cur_row = {'group1_%s_change' % ck: group1_change,
                       'group2_%s_change' % ck: group2_change,
                       'pvalue_%s_change' % ck: change_pvalue}
            row.update(cur_row)
        effects.append(row)

    effects_df = pd.DataFrame(effects)
    columns = ['group1', 'group2']
    for ck in uptake_keys:
        columns += ['group1_%s_uptake' % ck, 'group2_%s_uptake' % ck, 'pvalue_%s_uptake' % ck]

    for ck in change_keys:
        columns += ['group1_%s_change' % ck, 'group2_%s_change' % ck, 'pvalue_%s_change' % ck]

    effects_df = effects_df.ix[:,columns]


    # plot baseline value versus change scatter plot
    fig, ax = plt.subplots(1)
    cmap = cm.get_cmap('gist_rainbow')
    big_group_result_df = result_df[result_df['membership_prior'].isin(big_groups)]
    ax.scatter(big_group_result_df['CORTICAL_SUMMARY_prior'],big_group_result_df['CORTICAL_SUMMARY_change'],c=big_group_result_df['membership_prior'],cmap=cmap,s=50)
    plt.legend()
    plt.show()


    # plot against summary uptake values
    summary_prior_long = pd.melt(result_df, id_vars=['membership_prior'], value_vars='CORTICAL_SUMMARY_prior')
    summary_post_long = pd.melt(result_df, id_vars=['membership_post'], value_vars='CORTICAL_SUMMARY_post')
    plt.figure(1)
    for g in big_groups:
        prior_members = summary_prior_long[summary_prior_long.membership_prior==g][['variable','value']]
        post_members = summary_post_long[summary_post_long.membership_post==g][['variable','value']]
        members = pd.concat((prior_members,post_members))
        members['value'].plot(kind='kde', label=g, alpha=0.5)

    plt.legend()
    plt.show()

    # plot against summary change values
    summary_change_long = pd.melt(result_df, id_vars=['membership_prior'], value_vars='CORTICAL_SUMMARY_change')
    plt.figure(1)
    for g in big_groups:
        members = summary_change_long[summary_change_long.membership_prior==g]
        members['value'].plot(kind='kde', label=g, alpha=0.5)

    plt.legend()
    plt.show()

    # plot against wmh slope
    wmh_slope_long = pd.melt(result_df, id_vars=['membership_prior'], value_vars='WMH_slope')
    plt.figure(1)
    for g in positive_patterns:
        members = wmh_slope_long[wmh_slope_long.membership_prior==g]
        members['value'].plot(kind='kde', label=g, alpha=0.5)

    plt.legend()
    plt.show()

    # plot against uw ef
    ef_slope_long = pd.melt(result_df, id_vars=['membership_prior'], value_vars='UW_EF_slope')
    plt.figure(1)
    for g in positive_patterns:
        members = ef_slope_long[ef_slope_long.membership_prior==g]
        members['value'].plot(kind='kde', label=g, alpha=0.5)

    plt.legend()
    plt.show()

    # plot against uw mem
    mem_slope_long = pd.melt(result_df, id_vars=['membership_prior'], value_vars='UW_MEM_slope')
    plt.figure(1)
    for g in positive_patterns:
        members = mem_slope_long[mem_slope_long.membership_prior==g]
        members['value'].plot(kind='kde', label=g, alpha=0.5)

    plt.legend()
    plt.show()


    # plot against regional change
    rchange_members = rchange_df.merge(result_df, left_index=True, right_index=True)
    for i,rchange_key in enumerate(rchange_keys):
        long = pd.melt(rchange_members, id_vars='membership', value_vars=rchange_key)
        plt.figure(i+1)
        for g in groups:
            members = long[long['membership']==g]
            if len(members) < 3:
                continue
            members['value'].plot(kind='density', label='%s_%s' % (rchange_key, g))
        plt.legend()

    # plot against regional bl value
    uptake_members = uptake_prior_df.merge(result_df, left_index=True, right_index=True)
    for i,pattern_key in enumerate(pattern_keys):
        long = pd.melt(uptake_members, id_vars='membership', value_vars=pattern_key)
        plt.figure(i+1)
        for g in groups:
            members = long[long['membership']==g]
            if len(members) < 3:
                continue
            members['value'].plot(kind='density', label='%s_%s' % (pattern_key, g))
        plt.legend()





# # Leave one out X-validation
# loo = LeaveOneOut(n=len(pattern_df_pca.index))
# target = 'cortical_summary_bl'
# points = []
# for train_index, test_index in loo:
#     train_rows = pattern_df_pca.iloc[train_index]
#     test_rows = pattern_df_pca.iloc[test_index]
#     X_train = train_rows[pca_keys]
#     Y_train = train_rows[target].tolist()
#     X_test = test_rows[pca_keys]
#     Y_test = test_rows[target].tolist()
#     regressor = DecisionTreeRegressor(random_state=0)
#     regressor.fit(X_train,Y_train)
#     score = regressor.predict(X_test)
#     points.append((Y_test[0],score[0]))
#     print "%s -> %s" % (score, Y_test)
# plt.figure(1)
# plt.scatter([_[0] for _ in points],[_[1] for _ in points])
# plt.show()
