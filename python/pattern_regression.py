import pandas as pd
import numpy as np
import random
from collections import Counter, defaultdict
import itertools
import multiprocessing as mp
import cPickle


from scipy.stats import f_oneway, norm
import matplotlib.pyplot as plt
from matplotlib import cm
from pylab import get_cmap
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from pandas.stats.api import ols
import seaborn as sns
from sklearn.mixture import GMM, VBGMM, DPGMM
from sklearn.mixture.dpgmm import digamma
from sklearn.decomposition import PCA, KernelPCA
from sklearn.cross_validation import LeaveOneOut
from sklearn.mixture import GMM, VBGMM, DPGMM
from sklearn.tree import DecisionTreeRegressor
from sklearn.preprocessing import StandardScaler
from statsmodels.sandbox.stats.multicomp import multipletests
import networkx as nx

from utils import saveFakeAparcInput, importFreesurferLookup, importMaster, bilateralTranslations

pd.options.display.mpl_style = 'default'

bilateral=True

result_keys = ['CORTICAL_SUMMARY_post', 'CORTICAL_SUMMARY_prior', 'CORTICAL_SUMMARY_change', 'diag_prior', 'diag_post']
master_keys = ['APOE4_BIN','Age@AV45','Gender','Handedness',
               'UW_MEM_BL_3months','UW_MEM_slope',
               'UW_EF_BL_3months','UW_EF_slope',
               'WMH_percentOfICV_AV45_6MTHS','WMH_percentOfICV_slope',
               'CSF_TAU_closest_AV45','CSF_TAU_slope',
               'CSF_ABETA_closest_AV45','CSF_ABETA_slope',
               'FSX_HC/ICV_BL_3months','FSX_HC/ICV_slope',
               'FAQTOTAL_AV45_6MTHS','FAQTOTAL_slope',
               'NPITOTAL_AV45_6MTHS','NPITOTAL_slope',
               'MMSEslope_postAV45','MMSE_AV45_3MTHS',
               '','FDG_postAV45_slope']
summary_keys = ['CORTICAL_SUMMARY_post', 'CORTICAL_SUMMARY_prior', 'CORTICAL_SUMMARY_change']
pattern_keys = ['BRAIN_STEM', 'CTX_LH_BANKSSTS', 'CTX_LH_CAUDALANTERIORCINGULATE', 'CTX_LH_CAUDALMIDDLEFRONTAL', 'CTX_LH_CUNEUS', 'CTX_LH_ENTORHINAL', 'CTX_LH_FRONTALPOLE', 
                'CTX_LH_FUSIFORM', 'CTX_LH_INFERIORPARIETAL', 'CTX_LH_INFERIORTEMPORAL', 'CTX_LH_INSULA', 'CTX_LH_ISTHMUSCINGULATE', 'CTX_LH_LATERALOCCIPITAL', 
                'CTX_LH_LATERALORBITOFRONTAL', 'CTX_LH_LINGUAL', 'CTX_LH_MEDIALORBITOFRONTAL', 'CTX_LH_MIDDLETEMPORAL', 'CTX_LH_PARACENTRAL', 'CTX_LH_PARAHIPPOCAMPAL', 
                'CTX_LH_PARSOPERCULARIS', 'CTX_LH_PARSORBITALIS', 'CTX_LH_PARSTRIANGULARIS', 'CTX_LH_PERICALCARINE', 'CTX_LH_POSTCENTRAL', 'CTX_LH_POSTERIORCINGULATE', 
                'CTX_LH_PRECENTRAL', 'CTX_LH_PRECUNEUS', 'CTX_LH_ROSTRALANTERIORCINGULATE', 'CTX_LH_ROSTRALMIDDLEFRONTAL', 'CTX_LH_SUPERIORFRONTAL', 'CTX_LH_SUPERIORPARIETAL', 
                'CTX_LH_SUPERIORTEMPORAL', 'CTX_LH_SUPRAMARGINAL', 'CTX_LH_TEMPORALPOLE', 'CTX_LH_TRANSVERSETEMPORAL', 'CTX_RH_BANKSSTS', 'CTX_RH_CAUDALANTERIORCINGULATE', 
                'CTX_RH_CAUDALMIDDLEFRONTAL', 'CTX_RH_CUNEUS', 'CTX_RH_ENTORHINAL', 'CTX_RH_FRONTALPOLE', 'CTX_RH_FUSIFORM', 'CTX_RH_INFERIORPARIETAL', 
                'CTX_RH_INFERIORTEMPORAL', 'CTX_RH_INSULA', 'CTX_RH_ISTHMUSCINGULATE', 'CTX_RH_LATERALOCCIPITAL', 'CTX_RH_LATERALORBITOFRONTAL', 'CTX_RH_LINGUAL', 
                'CTX_RH_MEDIALORBITOFRONTAL', 'CTX_RH_MIDDLETEMPORAL', 'CTX_RH_PARACENTRAL', 'CTX_RH_PARAHIPPOCAMPAL', 'CTX_RH_PARSOPERCULARIS', 'CTX_RH_PARSORBITALIS', 
                'CTX_RH_PARSTRIANGULARIS', 'CTX_RH_PERICALCARINE', 'CTX_RH_POSTCENTRAL', 'CTX_RH_POSTERIORCINGULATE', 'CTX_RH_PRECENTRAL', 'CTX_RH_PRECUNEUS', 
                'CTX_RH_ROSTRALANTERIORCINGULATE', 'CTX_RH_ROSTRALMIDDLEFRONTAL', 'CTX_RH_SUPERIORFRONTAL', 'CTX_RH_SUPERIORPARIETAL', 'CTX_RH_SUPERIORTEMPORAL', 
                'CTX_RH_SUPRAMARGINAL', 'CTX_RH_TEMPORALPOLE', 'CTX_RH_TRANSVERSETEMPORAL', 'LEFT_ACCUMBENS_AREA', 'LEFT_AMYGDALA', 'LEFT_CAUDATE', 'LEFT_CEREBELLUM_CORTEX', 
                'LEFT_CEREBELLUM_WHITE_MATTER', 'LEFT_CEREBRAL_WHITE_MATTER', 'LEFT_HIPPOCAMPUS', 'LEFT_PALLIDUM', 'LEFT_PUTAMEN', 'LEFT_THALAMUS_PROPER', 
                'LEFT_VENTRALDC', 'RIGHT_ACCUMBENS_AREA', 'RIGHT_AMYGDALA', 'RIGHT_CAUDATE', 'RIGHT_CEREBELLUM_CORTEX', 'RIGHT_CEREBELLUM_WHITE_MATTER', 
                'RIGHT_CEREBRAL_WHITE_MATTER', 'RIGHT_HIPPOCAMPUS', 'RIGHT_PALLIDUM', 'RIGHT_PUTAMEN', 'RIGHT_THALAMUS_PROPER', 'RIGHT_VENTRALDC']
if bilateral: # truncate keys
    pattern_keys = list(set([_.replace('LH_','').replace('RH_','').replace('RIGHT_','').replace('LEFT_','') for _ in pattern_keys]))

prior_keys = ['%s_prior' % _ for _ in pattern_keys]
post_keys = ['%s_post' % _ for _ in pattern_keys]
change_keys = ['%s_change' % _ for _ in pattern_keys]


class DPGMM_SCALE(DPGMM):

    def __init__(self, n_components=1, covariance_type='diag', 
                 gamma_shape=1.0, gamma_inversescale=0.1,
                 random_state=None, tol=1e-3, verbose=0, n_iter=10, 
                 min_covar=None, params='wmc', init_params='wmc'):
        self.gamma_shape = float(gamma_shape)
        self.gamma_inversescale = float(gamma_inversescale)
        super(DPGMM_SCALE, self).__init__(n_components, covariance_type, alpha=gamma_shape/gamma_inversescale,
                                          random_state=random_state, tol=tol, min_covar=min_covar,
                                          n_iter=n_iter, params=params,
                                          init_params=init_params, verbose=verbose)

    def _update_concentration(self, z):
        """
            Update priors for the scaling parameter + update scale parameter,
            then update the concentration parameters for each cluster
        """
        # update scale gamma prior variational parameters
        new_gamma_shape = self.gamma_shape + self.n_components - 1
        dgamma = np.zeros(self.n_components-1)
        for i in range(self.n_components-1):
            dgamma[i] = digamma(self.gamma_[i,2]) - digamma(self.gamma_[i,1] + self.gamma_[i,2])
        new_gamma_inversescale = self.gamma_inversescale - np.sum(dgamma)
        del dgamma

        # update scale and concentration parameters
        self.alpha = new_gamma_shape / new_gamma_inversescale
        sz = np.sum(z, axis=0)
        self.gamma_.T[1] = 1. + sz
        self.gamma_.T[2].fill(0)
        for i in range(self.n_components - 2, -1, -1):
            self.gamma_[i, 2] = self.gamma_[i + 1, 2] + sz[i]
        self.gamma_.T[2] += self.alpha


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

def chooseBestModel(patterns, components, gamma_shape=7.5, gamma_inversescale=1.0, covar_type='spherical'):
    # Cluster
    models_scores = {}
    for i in range(30):
        g = DPGMM_SCALE(n_components=components,
                        gamma_shape=gamma_shape,
                        gamma_inversescale=gamma_inversescale,
                        tol=1e-6,
                        n_iter=700,
                        params='wmc',
                        init_params='wmc',
                        verbose=False)
        print "Fitting: %s" % i
        g.fit(patterns)
        print g.converged_
        print g.alpha
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

'''
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
        for i in range(10):
            g = DPGMM(n_components=components, 
                      covariance_type=covar_type,
                      alpha=a,
                      tol=1e-3,
                      n_iter=700,
                      params='wmc', 
                      init_params='wmc',
                      verbose=False)
            g.fit(data)
            print g.converged_
            y_ = g.predict(data)
            counter = dict(Counter(y_))
            score = np.mean(g.score(data))
            bic = g.bic(data)
            scores.append(bic)
            num_clusters.append(len(counter.keys()))
            var_clusters.append(np.var(counter.values()))
        print np.mean(scores)
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
'''

def compareTwoGroups(group1, group2, df, keys, adjust=True):
    group1_members = df[df['membership']==group1]
    group2_members = df[df['membership']==group2]
    print "GRP %s (%s) vs GRP %s (%s)" % (group1, len(group1_members), group2, len(group2_members))
    # bonferri correct
    hypos = len(keys)
    level = 0.05
    results = {}
    def generateArgs(group1_members, group2_members, keys):
        for i, key in enumerate(keys):
            group1_vals = group1_members[key].dropna().tolist()
            group2_vals = group2_members[key].dropna().tolist()
            if len(group1_vals) == 0:
                print "\tNO GROUP1 values for %s" % key
                continue
            if len(group2_vals) == 0:
                print "\tNO GROUP2 values for %s" % key
                continue
            yield (key, group1_vals, group2_vals)
    for arg in generateArgs(group1_members, group2_members, keys):
        key, pvalue, group1_boot, group2_boot = mp_bootstrap_wrapper(arg)
        results[key] = (pvalue, group1_boot, group2_boot)
        print "\t%s: %s (%s, %s)" % (key,pvalue,group1_boot,group2_boot)
    result_keys = results.keys()
    result_pvals = [results[_][0] for _ in result_keys]
    if adjust:
        reject, pvals_corrected, alphaSidak, alphaBonf = multipletests(result_pvals,alpha=level,method='hs')
    else:
        pvals_corrected = result_pvals
    for i,k in enumerate(result_keys):
        results[k] = (pvals_corrected[i],results[k][1],results[k][2])
        if pvals_corrected[i] <= level:
            print "\tSIG %s: %s (%s, %s)" % (k,pvals_corrected[i],results[k][1],results[k][2])
    return results

def group_comparisons(df, groups, keys, adjust=True):
    data = {}
    for group1, group2 in itertools.combinations(groups,2):
        data[(group1,group2)] = compareTwoGroups(group1, group2, df, keys, adjust=adjust)
    return data

def parseRawInput(patterns_csv, ref_key):
    # read in data
    raw_df = pd.read_csv(patterns_csv)
    raw_df.set_index('rid',inplace=True,drop=True)

    # convert to SUVR units and calculate change
    priors = [_ for _ in raw_df.columns if '_prior' in _ and 'diag' not in _]
    posts = [_ for _ in raw_df.columns if '_post' in _ and 'diag' not in _]
    prior_ref = raw_df.loc[:,'%s_prior' % ref_key]
    post_ref = raw_df.loc[:,'%s_post' % ref_key]

    raw_df[priors] = raw_df[priors].divide(prior_ref,axis='index')
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
    rchange_df = raw_df[change_keys]

    return (pattern_prior_df, pattern_post_df, uptake_prior_df, uptake_post_df, result_df, rchange_df)

def loadResults(results_file):
    results = pd.read_csv(results_file)
    results.set_index('rid',inplace=True,drop=True)
    return results

def bigGroups(result_df, threshold=7):
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

def parseConversions(groups, result_df, threshold, master_keys):
    # group diagnostic/pattern/positivity conversions
    diags = ['N', 'SMC', 'MCI', 'AD']
    conversions = {}
    for g in groups:
        cur_conversions = {}
        conversions[g] = {}
        members_prior = result_df[result_df.membership_prior==g]
        members_post = result_df[result_df.membership_post==g]
        # add membership counts
        cur_conversions['count_prior'] = len(members_prior.index)
        cur_conversions['count_post'] = len(members_post.index)
        # add master field avgs
        master_rows = result_df.loc[members_prior.index,master_keys]
        cur_conversions.update(dict(master_rows.mean()))
        # add cortical summary values
        cortical_summary_values = np.array(list(members_prior['CORTICAL_SUMMARY_prior']) + list(members_post['CORTICAL_SUMMARY_post']))
        cur_conversions['cortical_summary'] = cortical_summary_values.mean()
        cur_conversions['below_threshold'] = fitNormalCdf(cortical_summary_values, threshold)
        # add cortical summary change
        cur_conversions['cortical_summary_change'] = members_prior['CORTICAL_SUMMARY_change'].mean()
        # add diagnostic counts
        prior_counts = Counter(list(members_prior['diag_prior']))
        prior_counts['MCI']  = prior_counts.get('EMCI',0) + prior_counts.get('LMCI',0)
        post_counts = Counter(list(members_post['diag_post']))
        post_counts['MCI']  = post_counts.get('EMCI',0) + post_counts.get('LMCI',0)
        diag_counts = {diag: prior_counts.get(diag,0) + post_counts.get(diag,0) for diag in diags}
        total_counts = float(sum(diag_counts.values()))
        if total_counts > 0:
            diag_counts = {k:v/total_counts for k,v in diag_counts.items()}
        cur_conversions.update(diag_counts)
        # add diagnostic conversions
        diag_conversions = defaultdict(int)
        for rid, row in members_prior.iterrows():
            prior_diag = row['diag_prior']
            post_diag = row['diag_post']
            if pd.isnull(row['diag_prior']) or pd.isnull(row['diag_post']):
                continue
            if prior_diag in ['EMCI','LMCI']:
                prior_diag = 'MCI'
            if post_diag in ['EMCI','LMCI']:
                post_diag = 'MCI'
            diag_conversions['%s-%s' % (prior_diag,post_diag)] += 1
        total_counts = float(sum(diag_conversions.values()))
        if total_counts > 0:
            diag_conversions = {k:v/total_counts for k,v in diag_conversions.items()}
        cur_conversions.update(diag_conversions)
        # add positivity conversions
        pos_conversions = {}
        pos_conversions['neg-neg'] = len(members_prior[(members_prior.CORTICAL_SUMMARY_prior<threshold) & (members_prior.CORTICAL_SUMMARY_post<threshold)].index)
        pos_conversions['neg-pos'] = len(members_prior[(members_prior.CORTICAL_SUMMARY_prior<threshold) & (members_prior.CORTICAL_SUMMARY_post>=threshold)].index)
        pos_conversions['pos-pos'] = len(members_prior[(members_prior.CORTICAL_SUMMARY_prior>=threshold) & (members_prior.CORTICAL_SUMMARY_post>=threshold)].index)
        pos_conversions['pos-neg'] = len(members_prior[(members_prior.CORTICAL_SUMMARY_prior>=threshold) & (members_prior.CORTICAL_SUMMARY_post<threshold)].index)
        total_counts = float(sum(pos_conversions.values()))
        if total_counts > 0:
            pos_conversions = {k:v/total_counts for k,v in pos_conversions.items()}
        cur_conversions.update(pos_conversions)
        # pattern conversions
        pattern_conversions = {}
        for g2 in groups:
            pattern_conversions[g2] = len(members_prior[members_prior.membership_post==g2].index)
        total_counts = float(sum(pattern_conversions.values()))
        if total_counts > 0:
            pattern_conversions = {k:v/total_counts for k,v in pattern_conversions.items()}
        cur_conversions.update(pattern_conversions)
        # add to results
        conversions[g] = cur_conversions
    conversions = pd.DataFrame(conversions).T
    conversions.index.name = 'pattern'
    return conversions

def saveAparcs(alpha, components, groups, pattern_members, uptake_members, change_members, index_lookup):
    aparc_input_template = "../output/fake_aparc_inputs/dpgmm_alpha%s_comp%s_group_%s_%s"
    for g in groups:
        out_file_pattern = aparc_input_template % (alpha, components, g, 'pattern')
        out_file_change = aparc_input_template % (alpha, components, g, 'change')
        out_file_uptake = aparc_input_template % (alpha, components, g, 'uptake')
        pattern_values = pattern_members[pattern_members['membership']==g][pattern_keys].mean().to_dict()
        uptake_values = uptake_members[uptake_members['membership']==g][pattern_keys].mean().to_dict()
        change_values = change_members[change_members['membership']==g][change_keys].mean().to_dict()
        change_values = {k.replace('_change',''):v for k,v in change_values.iteritems()}
        saveFakeAparcInput(out_file_pattern, pattern_values, index_lookup)
        saveFakeAparcInput(out_file_change, change_values, index_lookup)
        saveFakeAparcInput(out_file_uptake, uptake_values, index_lookup)

def mergeResults(result_df, pattern_prior_df, pattern_post_df, uptake_prior_df, uptake_post_df, rchange_df, with_timepoint=False):
    # merge results with regional change and uptakes
    change_members = rchange_df.copy()
    change_members['membership'] = result_df['membership_prior']
    if with_timepoint:
        change_members.index.name='rid'
        change_members.loc[:,'timepoint'] = 'prior'
        change_members.set_index('timepoint',append=True,inplace=True)
    prior_members = result_df.copy()
    prior_members['membership'] = result_df['membership_prior']
    if with_timepoint:
        prior_members.index.name='rid'
        prior_members.loc[:,'timepoint'] = 'prior'
        prior_members.set_index('timepoint',append=True,inplace=True)
    uptake_members_prior = uptake_prior_df.merge(result_df[['membership_prior','CORTICAL_SUMMARY_prior']], left_index=True, right_index=True)
    uptake_members_prior['membership'] = uptake_members_prior['membership_prior']
    uptake_members_prior['CORTICAL_SUMMARY'] = uptake_members_prior['CORTICAL_SUMMARY_prior']
    uptake_members_post = uptake_post_df.merge(result_df[['membership_post','CORTICAL_SUMMARY_post']], left_index=True, right_index=True)
    uptake_members_post['membership'] = uptake_members_post['membership_post']
    uptake_members_post['CORTICAL_SUMMARY'] = uptake_members_post['CORTICAL_SUMMARY_post']
    if with_timepoint:
        uptake_members_prior.index.name='rid'
        uptake_members_post.index.name='rid'
        uptake_members_prior.loc[:,'timepoint'] = 'prior'
        uptake_members_post.loc[:,'timepoint'] = 'post'
        uptake_members_prior.set_index('timepoint',append=True,inplace=True)
        uptake_members_post.set_index('timepoint',append=True,inplace=True)
    uptake_members = pd.concat((uptake_members_prior,uptake_members_post))
    pattern_prior_df = pattern_prior_df.xs('prior',level='timepoint')
    pattern_post_df = pattern_post_df.xs('post',level='timepoint')
    pattern_members_prior = pattern_prior_df.merge(result_df[['membership_prior','CORTICAL_SUMMARY_prior']], left_index=True, right_index=True)
    pattern_members_post = pattern_post_df.merge(result_df[['membership_post','CORTICAL_SUMMARY_post']], left_index=True, right_index=True)
    pattern_members_prior['membership'] = pattern_members_prior['membership_prior']
    pattern_members_prior['CORTICAL_SUMMARY'] = pattern_members_prior['CORTICAL_SUMMARY_prior']
    pattern_members_post['membership'] = pattern_members_post['membership_post']
    pattern_members_post['CORTICAL_SUMMARY'] = pattern_members_post['CORTICAL_SUMMARY_post']
    if with_timepoint:
        pattern_members_prior.index.name='rid'
        pattern_members_post.index.name='rid'
        pattern_members_prior.loc[:,'timepoint'] = 'prior'
        pattern_members_post.loc[:,'timepoint'] = 'post'
        pattern_members_prior.set_index('timepoint',append=True,inplace=True)
        pattern_members_post.set_index('timepoint',append=True,inplace=True)
    pattern_members = pd.concat((pattern_members_prior,pattern_members_post))
    return (uptake_members, pattern_members, change_members, prior_members)

def scaleRawInput(pattern_df, scale_type='original', whiten=True):
    assert scale_type in set(['original', 'pca', 'kernelpca'])
    if scale_type == 'original':
        # Scale in original space
        scaler = StandardScaler().fit(pattern_df)
        pattern_df_raw = pd.DataFrame(scaler.transform(pattern_df))
        pattern_df_raw.set_index(pattern_df.index, inplace=True)
        return (pattern_df_raw, scaler)
    elif scale_type == 'pca':
        # convert to PCA space + whiten(scale)
        pca_model = PCA(n_components=len(pattern_df.columns), copy=True, whiten=whiten)
        pattern_df_pca = pd.DataFrame(pca_model.fit_transform(pattern_df))
        pattern_df_pca.set_index(pattern_df.index, inplace=True)
        return (pattern_df_pca, pca_model)
    elif scale_type == 'kernelpca':
        # convert to Kernel PCA Projection and scale
        kpca_model = KernelPCA(n_components=len(pattern_df.columns),kernel='rbf',fit_inverse_transform=False,gamma=1,alpha=1.0)
        X_kpca = kpca_model.fit_transform(pattern_df)
        pattern_df_kpca = pd.DataFrame(X_kpca)
        scaler = StandardScaler().fit(pattern_df_kpca)
        pattern_df_kpca = pd.DataFrame(scaler.transform(pattern_df_kpca))
        pattern_df_kpca.set_index(pattern_df.index, inplace=True)
        return (pattern_df_kpca, scaler)

def graphNetworkConversions(groups, conversions, iterations=50, threshold=0.0, alternate_nodes=None):
    G = nx.DiGraph()
    cmap = get_cmap("cool")
    if alternate_nodes:
        alternate_nodes = [int(_) for _ in alternate_nodes]
    for g1 in groups:
        for g2 in groups:
            if g1 == g2:
                if g1 not in G.nodes():
                    G.add_node(int(g1))
                continue
            conv_percent = conversions.loc[g1,g2]
            if conv_percent > threshold:
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
    nodecolor = []
    for n in G.nodes():
        sizes[n] = conversions.loc[n,'count_prior'] + conversions.loc[n,'count_post']
        if n in alternate_nodes:
            nodecolor.append('r')
        else:
            nodecolor.append('w')
    nodesize = [sizes[_]*80 for _ in G]
    pos = nx.spring_layout(G,iterations=iterations)
    #pos = nx.circular_layout(G)
    plt.figure(figsize=(10,10))
    #nx.draw_networkx_edges(G,pos,width=edgewidth,edge_color='m',alpha=0.3)
    nodes = nx.draw_networkx_nodes(G,pos,
                                   node_size=nodesize,
                                   node_color=nodecolor,
                                   alpha=0.4)
    nodes.set_edgecolor('k')
    nx.draw_networkx_edges(G,pos,alpha=0.25,node_size=0,width=4,edge_color=edgecolors)
    #nx.draw_networkx_edge_labels(G,pos,edge_labels=edgelabels)
    nx.draw_networkx_labels(G,pos,font_size=14,font_family='sans-serif')
    plt.axis('off')
    plt.show()

def plotValueDensity(result_df, groups, value_key):
    long = pd.melt(result_df, id_vars=['membership_prior'], value_vars=value_key)
    for i, (k, patterns) in enumerate(groups.iteritems()):
        plt.figure(i)
        plt.title("%s (%s)" % (value_key,k))
        for g in patterns:
            members = long[long.membership_prior==g]
            print len(members.index)
            members['value'].plot(kind='density', label=g, alpha=0.5)
        plt.legend()
    plt.show()

def plotValueBox(result_df, groups, value_key, save=False):
    long = pd.melt(result_df, id_vars=['membership_prior'], value_vars=value_key)
    for i, (groupk, patterns) in enumerate(groups.iteritems()):
        by_group = pd.DataFrame()
        for g in patterns:
            members = long[long.membership_prior==g][['value']]
            members['pattern'] = g
            by_group = pd.concat((by_group,members))
        dfg = by_group.groupby('pattern')
        labels = ['%s\n$n$=%d'%(k, len(v)) for k, v in dfg]
        counts = [len(v) for k,v in dfg]
        total = float(sum(counts))
        widths = [2*c/total for c in counts]
        plt.figure(i)
        bplot = sns.violinplot(x='pattern',y='value',data=by_group,palette='Set3')
        bplot.set_xticklabels(labels)
        plt.title("%s (%s)" % (value_key,groupk))
        plt.suptitle("")
        if save:
            plt.savefig('../boxplot_%s_%s.png' % (value_key.replace('/','_'), groupk), dpi=200)
        else:
            plt.show()
        plt.close()

def flattenGroupComparisonResults(change_pvalues):
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
    flattened_df = pd.DataFrame(flattened)
    flattened_df.set_index(['GROUP1','GROUP2'],inplace=True)
    return flattened_df

def regionRanks(groups, members, keys):
    '''
    remove white matter, ref regions
    '''
    blacklist = ['WHITE_MATTER','BRAIN_STEM','VENTRALDC']
    filtered_keys = [k for k in keys if not any([bl in k for bl in blacklist])]
    ranks = pd.DataFrame()
    for g in groups:
        ranked_regions = members[members.membership_prior==g].mean().loc[filtered_keys].sort(ascending=False,inplace=False)
        ranked_regions = pd.DataFrame(ranked_regions).reset_index()
        ranked_regions.columns = ['%s_RANK' % int(g),'%s_UPTAKE' % int(g)]
        ranks = pd.concat((ranks,ranked_regions.T))
    return ranks


def compareGroupPairs_worker(args):
    g1, g2, members, keys = args
    expanded = {}
    res = compareTwoGroups(g1, g2, members, keys, adjust=False)
    for k, (pval, m1, m2) in res.iteritems():
        expanded['%s_PVALUE' % k] = pval
        expanded['%s_MEAN1' % k] = m1
        expanded['%s_MEAN2' % k] = m2 
    expanded['GROUP1'] = g1
    expanded['GROUP2'] = g2
    return expanded

def compareGroupPairs(pairs, members, keys):
    key_order = ['GROUP1','GROUP2']
    key_order += [subkey for key in keys for subkey in ['%s_PVALUE' % key,'%s_MEAN1' % key,'%s_MEAN2' % key]]
    all_results = []
    def argGenerator(pairs, members, keys):
        for g1, g2 in pairs:
            yield (g1, g2, members, keys)
    pool = mp.Pool(processes=4)
    all_results = []
    for expanded in pool.imap(compareGroupPairs_worker, argGenerator(pairs,members,keys), 1):
        all_results.append(expanded)
    all_results_df = pd.DataFrame(all_results)
    all_results_df = all_results_df[key_order]
    return all_results_df


if __name__ == '__main__':
    # import master
    master_csv = '../FDG_AV45_COGdata_12_03_15.csv'
    master_df = pd.read_csv(master_csv, low_memory=False, header=[0,1])
    master_df.columns = master_df.columns.get_level_values(1)
    master_df.set_index('RID', inplace=True)
    master_df = master_df.loc[:,master_keys]

    # parse input
    ref_key = 'COMPOSITE_REF'
    patterns_csv = '../datasets/pvc_allregions_uptake_change_bilateral.csv'
    threshold = 0.91711
    membership_conf = 0.50
    components = 100

    pattern_prior_df, pattern_post_df, uptake_prior_df, uptake_post_df, result_df, rchange_df = parseRawInput(patterns_csv, ref_key)

    # Scale inputs
    scale_type = 'original'
    patterns_only, scaler = scaleRawInput(pattern_prior_df, scale_type=scale_type)
    post_patterns_only = pd.DataFrame(scaler.transform(pattern_post_df))
    post_patterns_only.set_index(pattern_post_df.index, inplace=True)

    # # Look at top pca components (top 10)
    # top = 10
    # columns = list(pattern_prior_df.columns)
    # pca_patterns_only, pca_model = scaleRawInput(patterns_only, scale_type='pca', whiten=False)
    # pca_post_patterns_only = pd.DataFrame(pca_model.transform(post_patterns_only))
    # pca_post_patterns_only.set_index(post_patterns_only.index, inplace=True)
    # pca_patterns_only = pca_patterns_only[range(top)]
    # pca_post_patterns_only = pca_post_patterns_only[range(top)]
    # top_components = []
    # for i in range(top):
    #     explained_ratio = pca_model.explained_variance_ratio_[i]
    #     component = pca_model.components_[i]
    #     comps = zip(columns,component)
    #     comps = sorted(comps,key=lambda x: abs(x[1]), reverse=True)
    #     comps = pd.DataFrame(comps, columns=['REGION','WEIGHT'])
    #     top_components.append(comps)

    # # Choose alpha + do model selection
    # best_model = chooseBestModel(patterns_only, components, gamma_shape=1.0, gamma_inversescale=0.1, covar_type='spherical')
    # alpha = best_model.alpha
    # with open('../dpgmm_alpha%s_spherical_bilateral_model.pkl' % round(alpha,2), 'wb') as fid:
    #     cPickle.dump(best_model, fid)   

    # Load model instead
    with open('../dpgmm_alpha13.93_spherical_bilateral_model.pkl', 'rb') as fid:
        best_model = cPickle.load(fid)
        alpha = best_model.alpha

    # # Predict pattern groups and add to result df
    # y_df = pd.DataFrame(best_model.predict(patterns_only))
    # y_df.columns = ['group']
    # y_df.set_index(patterns_only.index, inplace=True)
    # y_proba_df = pd.DataFrame(best_model.predict_proba(patterns_only)).set_index(patterns_only.index).xs('prior',level='timepoint')
    # probsums = y_proba_df.sum()
    # components_distr_max = y_proba_df[probsums[probsums >= 0.1].index].max(axis=1)
    # confident_assignments_prior = components_distr_max[components_distr_max>=membership_conf].index

    # y_df_post = pd.DataFrame(best_model.predict(post_patterns_only))
    # y_df_post.columns = ['group']
    # y_df_post.set_index(post_patterns_only.index, inplace=True)
    # y_proba_df_post = pd.DataFrame(best_model.predict_proba(post_patterns_only)).set_index(post_patterns_only.index).xs('post',level='timepoint')
    # probsums = y_proba_df_post.sum()
    # components_distr_max = y_proba_df_post[probsums[probsums > 0.1].index].max(axis=1)
    # confident_assignments_post = components_distr_max[components_distr_max>=membership_conf].index
    # for rid in result_df.index:
    #     result_df.loc[rid,'membership_prior'] = y_df.loc[(rid,'prior'),'group']
    #     for pca_col in pca_patterns_only.columns:
    #         result_df.loc[rid,'pca_%s_prior' % pca_col] = pca_patterns_only.loc[(rid, 'prior'),pca_col]
    #     if (rid, 'post') in y_df_post.index:
    #         result_df.loc[rid,'membership_post'] = y_df_post.loc[(rid,'post'),'group']
    #         for pca_col in pca_post_patterns_only.columns:
    #             result_df.loc[rid,'pca_%s_post' % pca_col] = pca_post_patterns_only.loc[(rid, 'post'),pca_col]
    #     if rid in confident_assignments_prior:
    #         result_df.loc[rid,'membership_prior_conf'] = result_df.loc[rid,'membership_prior']
    #     if rid in confident_assignments_post:
    #         result_df.loc[rid,'membership_post_conf'] = result_df.loc[rid,'membership_post']

    # # save
    # result_df.to_csv('../dpgmm_alpha%s_spherical_bilateral_result.csv' % round(alpha,2))

    # load results
    result_df = loadResults('../dpgmm_alpha%s_spherical_bilateral_result.csv' % round(alpha,2))

    # use group memberships filtered by confidence level
    result_df.drop(['membership_prior','membership_post'], axis=1, inplace=True)
    result_df.rename(columns={'membership_prior_conf':'membership_prior','membership_post_conf':'membership_post'},inplace=True)

    # generate conversion data
    result_df = result_df.merge(master_df,left_index=True,right_index=True)
    big_groups = bigGroups(result_df, threshold=3)
    conversions = parseConversions(big_groups, result_df, threshold, master_keys)
    # conversions.to_csv('../dpgmm_alpha%s_spherical_bilateral_conversions.csv' % (round(alpha,2)))
    positive_patterns = list(conversions[conversions['pos-pos']>=0.8].index)
    negative_patterns = list(conversions[conversions['neg-neg']>=0.9].index)
    transition_patterns = list(set(conversions.index) - (set(positive_patterns) | set(negative_patterns)))
    groups = {'positive': positive_patterns, 'negative': negative_patterns, 'transition': transition_patterns}
    uptake_members, pattern_members, change_members, prior_members = mergeResults(result_df, pattern_prior_df, pattern_post_df, uptake_prior_df, uptake_post_df, rchange_df)

    # create giant dataframe of merged results
    uptake_tp, pattern_tp, change_tp, prior_tp = mergeResults(result_df, pattern_prior_df, pattern_post_df, uptake_prior_df, uptake_post_df, rchange_df, with_timepoint=True)
    merged_members = prior_tp.copy().reset_index()
    merged_members = merged_members.merge(uptake_tp[pattern_keys].reset_index(), on=['rid','timepoint'], how='outer')
    merged_members = merged_members.merge(change_tp[change_keys].reset_index(), on=['rid','timepoint'], how='outer')
    all_keys = master_keys+summary_keys+pattern_keys+change_keys

    # # save fake aparcs
    # lut_file = "../FreeSurferColorLUT.txt"
    # if bilateral:
    #     index_lookup = bilateralTranslations(lut_file)
    # else:
    #     lut_table = importFreesurferLookup(lut_file)
    #     index_lookup = {v.replace('-','_').upper():[k] for k,v in lut_table.iteritems()}

    # saveAparcs(round(alpha,2), components, big_groups, pattern_members, uptake_members, change_members, index_lookup)

    # # pca component correlations
    # y_cols = ['CORTICAL_SUMMARY_prior',
    #           'CORTICAL_SUMMARY_change',
    #           'UW_MEM_BL_3months',
    #           'UW_EF_BL_3months',
    #           'CSF_TAU_closest_AV45',
    #           'CSF_ABETA_closest_AV45',
    #           'FSX_HC/ICV_BL_3months']
    # x_cols = ['pca_prior','CORTICAL_SUMMARY_prior']
    # props={'pad':10}
    # for y in y_cols:
    #     fig, axes = plt.subplots(nrows=1, ncols=2)
    #     # pca
    #     pca_x = ['pca_%s_prior' % i for i in range(5)]
    #     model = ols(x=result_df[pca_x],y=result_df[y])
    #     x = model.beta.x
    #     intercept = model.beta.intercept
    #     r2 = '{0:.3f}'.format(model.r2)
    #     equation = 'Y ~ %sx + %s' % ('{0:.3f}'.format(x),'{0:.3f}'.format(intercept))
    #     ax = axes[0]
    #     result_df.plot(ax=ax, kind='scatter',x='pca_0_prior',y=y)
    #     xvals = np.array(sorted(list(result_df['pca_0_prior'])))
    #     ax.plot(xvals, xvals*x + intercept, 'r')
    #     ax.set_title('Principal Component')
    #     ax.annotate('%s\nR2=%s' % (equation,r2), xy=(0.05,0.95), xycoords='axes fraction', bbox=props, verticalalignment='top', fontsize=15, color='k')
    #     # cortical summary
    #     model = ols(x=result_df['CORTICAL_SUMMARY_prior'],y=result_df[y])
    #     x = model.beta.x
    #     intercept = model.beta.intercept
    #     r2 = '{0:.3f}'.format(model.r2)
    #     equation = 'Y ~ %sx + %s' % ('{0:.3f}'.format(x),'{0:.3f}'.format(intercept))
    #     ax = axes[1]
    #     result_df.plot(ax=ax, kind='scatter',x='CORTICAL_SUMMARY_prior',y=y)
    #     xvals = np.array(sorted(list(result_df['CORTICAL_SUMMARY_prior'])))
    #     ax.plot(xvals, xvals*x + intercept, 'r')
    #     ax.set_title('Cortical Summary SUVR')
    #     ax.annotate('%s\nR2=%s' % (equation,r2), xy=(0.05,0.95), xycoords='axes fraction', bbox=props, verticalalignment='top', fontsize=15, color='k')

    # region ranking
    uptake_ranks = regionRanks(big_groups, uptake_members, pattern_keys)
    uptake_ranks.to_csv('../dpgmm_alpha%s_spherical_bilateral_uptake_ranks.csv' % round(alpha,2))
    change_ranks = regionRanks(big_groups, change_members, change_keys)
    change_ranks.to_csv('../dpgmm_alpha%s_spherical_bilateral_change_ranks.csv' % round(alpha,2))
    pattern_ranks = regionRanks(big_groups, pattern_members, pattern_keys)
    pattern_ranks.to_csv('../dpgmm_alpha%s_spherical_bilateral_pattern_ranks.csv' % round(alpha,2))

    # plot pattern flows
    graphNetworkConversions(transition_patterns+positive_patterns, conversions, iterations=50, threshold=0.1, alternate_nodes=transition_patterns)
    graphNetworkConversions(negative_patterns+transition_patterns, conversions, iterations=50, threshold=0.1, alternate_nodes=transition_patterns)
    graphNetworkConversions(negative_patterns+positive_patterns+transition_patterns, conversions, iterations=50, threshold=0.1, alternate_nodes=negative_patterns+transition_patterns)

    # plot group value distributions
    plotValueBox(result_df, groups, 'CORTICAL_SUMMARY_prior', save=True)
    plotValueBox(result_df, groups, 'CORTICAL_SUMMARY_change', save=True)
    plotValueBox(result_df, groups, 'WMH_percentOfICV_slope', save=True)
    plotValueBox(result_df, groups, 'UW_EF_slope', save=True)
    plotValueBox(result_df, groups, 'UW_MEM_slope', save=True)
    plotValueBox(result_df, groups, 'CSF_TAU_slope', save=True)
    plotValueBox(result_df, groups, 'CSF_ABETA_slope', save=True)
    plotValueBox(result_df, groups, 'CSF_ABETA_closest_AV45', save=True)
    plotValueBox(result_df, groups, 'FSX_HC/ICV_slope', save=True)

    # # derive cortical summary prior p values between all pairs
    # pvalue_df = compareGroupPairs(itertools.combinations(big_groups,2), merged_members, ['CORTICAL_SUMMARY_prior'])
    # pvalue_df.to_csv('../dpgmm_alpha%s_spherical_bilateral_cortical_summary_pvalues.csv' % (round(alpha,2)), index=False)
    
    # read cortical summary pvalues
    pvalue_df = pd.read_csv('../dpgmm_alpha%s_spherical_bilateral_cortical_summary_pvalues.csv' % (round(alpha,2)))
    pvalue_df.set_index(['GROUP1','GROUP2'],inplace=True)
    sig_same = pvalue_df[pvalue_df.CORTICAL_SUMMARY_prior_PVALUE>0.01]
    sig_same_pairs = [_[0] for _ in zip(sig_same.index)]

    # compare groups with non-significantly different cortical summary prior distributions
    pairs_comparisons_df = compareGroupPairs(sig_same_pairs, merged_members, keys)
    pairs_comparisons_df.set_index(['GROUP1','GROUP2'], inplace=True)
    pairs_comparisons_df = pairs_comparisons_df.T
    pairs_comparisons_df.to_csv('../dpgmm_alpha%s_spherical_bilateral_pair_comparisons.csv' % (round(alpha,2)))
    pvalue_indices = [_ for _ in pairs_comparisons_df.index if 'PVALUE' in _]
    # pvalues only
    pairs_comparisons_df_pvalues = pairs_comparisons_df.loc[pvalue_indices,:]
    pairs_comparisons_df_pvalues.to_csv('../dpgmm_alpha%s_spherical_bilateral_pair_pvalues.csv' % (round(alpha,2)), na_rep=-1)
    # means only
    reset_df = pairs_comparisons_df.T.reset_index().T
    mean1_indices = ['GROUP1'] + [_ for _ in pairs_comparisons_df.index if 'MEAN1' in _]
    mean2_indices = ['GROUP2'] + [_ for _ in pairs_comparisons_df.index if 'MEAN2' in _]
    clean_indices = [_.replace('1','').replace('_MEAN','') for _ in mean1_indices]
    mean1_cleaned = reset_df[mean1_indices]
    mean1_cleaned.columns = clean_indices
    mean1_cleaned.set_index('GROUP',inplace=True)
    mean2_cleaned = reset_df[mean2_indices]
    mean2_cleaned.columns = clean_indices
    mean2_cleaned.set_index('GROUP',inplace=True)
    all_cleaned = pd.concat((mean1_cleaned,mean2_cleaned))
    mean_averages = all_cleaned.groupby(level='GROUP').mean().T
    mean_averages.to_csv('../dpgmm_alpha%s_spherical_bilateral_pair_means.csv' % (round(alpha,2)), na_rep='NA')

    '''
    master_pvalues = group_comparisons(prior_members, transition_patterns, master_keys + summary_keys)
    flattened_master_df = flattenGroupComparisonResults(master_pvalues)
    flattened_master_df.to_csv('../dpgmm_alpha%s_spherical_bilateral_transition_master_pvalues.csv' % (round(alpha,2)))

    master_pvalues = group_comparisons(prior_members, positive_patterns, master_keys + summary_keys)
    flattened_master_df = flattenGroupComparisonResults(master_pvalues)
    flattened_master_df.to_csv('../dpgmm_alpha%s_spherical_bilateral_positive_master_pvalues.csv' % (round(alpha,2)))

    master_pvalues = group_comparisons(prior_members, negative_patterns, master_keys + summary_keys)
    flattened_master_df = flattenGroupComparisonResults(master_pvalues)
    flattened_master_df.to_csv('../dpgmm_alpha%s_spherical_bilateral_negative_master_pvalues.csv' % (round(alpha,2)))

    # between group bootstrap hypothesis tests
    change_pvalues = group_comparisons(change_members, transition_patterns, change_keys)
    flattened_change_df = flattenGroupComparisonResults(change_pvalues)
    flattened_change_df.to_csv('../dpgmm_alpha%s_spherical_bilateral_transition_change_pvalues.csv' % (round(alpha,2)))

    pattern_pvalues = group_comparisons(pattern_members, transition_patterns, pattern_keys)
    flattened_pattern_df = flattenGroupComparisonResults(pattern_pvalues)
    flattened_pattern_df.to_csv('../dpgmm_alpha%s_spherical_bilateral_transition_pattern_pvalues.csv' % (round(alpha,2)))

    uptake_pvalues = group_comparisons(uptake_members, transition_patterns, pattern_keys)
    flattened_uptake_df = flattenGroupComparisonResults(uptake_pvalues)
    flattened_uptake_df.to_csv('../dpgmm_alpha%s_spherical_bilateral_transition_uptake_pvalues.csv' % (round(alpha,2)))
    '''

    # plot baseline value versus change scatter plot
    fig, ax = plt.subplots(1)
    cmap = cm.get_cmap('gist_rainbow')
    big_group_result_df = result_df[result_df['membership_prior'].isin(big_groups)]
    ax.scatter(big_group_result_df['CORTICAL_SUMMARY_prior'],big_group_result_df['CORTICAL_SUMMARY_change'],c=big_group_result_df['membership_prior'],cmap=cmap,s=50)
    plt.legend()
    plt.show()

    # plot against regional change
    rchange_members = rchange_df.merge(result_df, left_index=True, right_index=True)
    for i,rchange_key in enumerate(change_keys):
        long = pd.melt(rchange_members, id_vars='membership', value_vars=rchange_key)
        plt.figure(i+1)
        for g in big_groups:
            members = long[long['membership']==g]
            members['value'].plot(kind='density', label='%s_%s' % (rchange_key, g))
        plt.legend()

    # plot against regional bl value
    uptake_members = uptake_prior_df.merge(result_df, left_index=True, right_index=True)
    for i,pattern_key in enumerate(pattern_keys):
        long = pd.melt(uptake_members, id_vars='membership', value_vars=pattern_key)
        plt.figure(i+1)
        for g in big_groups:
            members = long[long['membership']==g]
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
