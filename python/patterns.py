import pandas as pd
import numpy as np
import random
from collections import Counter, defaultdict
import itertools
import multiprocessing as mp
import cPickle

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, get_cmap

from scipy.stats import f_oneway, norm, chi2_contingency, linregress
from pandas.stats.api import ols
import seaborn as sns
import networkx as nx
import numpy.testing as npt
from statsmodels.sandbox.stats.multicomp import multipletests

from sklearn.mixture import GMM, VBGMM, DPGMM
from sklearn.mixture.dpgmm import digamma, _bound_state_log_lik, log_normalize
from sklearn.utils.validation import check_is_fitted
from sklearn.utils import check_array
from sklearn.decomposition import PCA, KernelPCA
from sklearn.cross_validation import LeaveOneOut
from sklearn.preprocessing import StandardScaler


from utils import saveFakeAparcInput, importFreesurferLookup, importMaster, bilateralTranslations, regCoeffZTest, slope

#plt.style.use('ggplot')
#pd.options.display.mpl_style = 'default'

bilateral=True
normals_only=False

patterns_csv = '../datasets/pvc_allregions_uptake_change_bilateral.csv'
patterns_split_csv = '../datasets/pvc_allregions_uptake_bilateral.csv'
lut_file = "../FreeSurferColorLUT.txt"
master_csv = '../FDG_AV45_COGdata/FDG_AV45_COGdata_01_26_16.csv'
membership_conf = 0.50
components = 1067
ref_key = 'WHOLE_CEREBELLUM'
threshold = 1.2813
# ref_key = 'COMPOSITE_REF'
# threshold = 0.91711

result_keys = ['CORTICAL_SUMMARY_post', 'CORTICAL_SUMMARY_prior', 'CORTICAL_SUMMARY_change', 'diag_prior', 'diag_post']
lobe_keys = ['FRONTAL','PARIETAL','CINGULATE','TEMPORAL','OCCIPITAL','MEDIALOCCIPITAL',
             'SENSORY','BASALGANGLIA','LIMBIC','CEREBGM','CEREBWM','CEREBRAL_WHITE','BRAIN_STEM']
master_keys = ['Age@AV45','Gender','APOE2_BIN','APOE4_BIN','Edu.(Yrs)','SMOKING','DIABETES',
               'UW_MEM_BL_3months','UW_MEM_slope',
               'UW_EF_BL_3months','UW_EF_slope',
               'WMH_percentOfICV_AV45_6MTHS','WMH_percentOfICV_slope',
               'CSF_TAU_closest_AV45','CSF_TAU_slope',
               'CSF_ABETA_closest_AV45','CSF_ABETA_slope',
               'FSX_HC/ICV_BL_3months','FSX_HC/ICV_slope',
               'FDG_PONS_AV45_6MTHS','FDG_postAV45_slope']
cat_keys = ['APOE4_BIN','APOE2_BIN','Gender','Handedness','SMOKING','DIABETES']
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
                'RIGHT_ACCUMBENS_AREA', 'RIGHT_AMYGDALA', 'RIGHT_CAUDATE', 'RIGHT_CEREBELLUM_CORTEX', 'RIGHT_CEREBELLUM_WHITE_MATTER', 
                'RIGHT_CEREBRAL_WHITE_MATTER', 'RIGHT_HIPPOCAMPUS', 'RIGHT_PALLIDUM', 'RIGHT_PUTAMEN', 'RIGHT_THALAMUS_PROPER']
if bilateral: # truncate keys
    pattern_keys = list(set([_.replace('LH_','').replace('RH_','').replace('RIGHT_','').replace('LEFT_','') for _ in pattern_keys]))

prior_keys = ['%s_prior' % _ for _ in pattern_keys]
post_keys = ['%s_post' % _ for _ in pattern_keys]
change_keys = ['%s_change' % _ for _ in pattern_keys]
lobe_change_keys = ['%s_change' % _ for _ in lobe_keys]


class DPGMM_SCALE(DPGMM):

    def __init__(self, n_components=1, covariance_type='diag', 
                 gamma_shape=1.0, gamma_inversescale=0.000001,
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


def sample(data):
    return data[np.random.randint(0,len(data),(1,len(data)))[0]]

def mp_bootstrap_wrapper(arg):
    key, group1_vals, group2_vals = arg
    group1_boot = bootstrap_mean(group1_vals)
    group2_boot = bootstrap_mean(group2_vals)
    if key in cat_keys:
        pvalue = chi2_test(group1_vals, group2_vals)
        print "chi2 %s: %s" % (key,pvalue)
    else:
        pvalue = bootstrap_t_test(group1_vals, group2_vals)
    return (key, pvalue, group1_boot, group2_boot)

def chi2_test(group1_vals, group2_vals):
    g1_counter = Counter(group1_vals)
    g2_counter = Counter(group2_vals)
    all_keys = list(set(g1_counter.keys()) | set(g2_counter.keys()))
    freqs = []
    freqs.append({g1_counter.get(k,0) for k in all_keys})
    freqs.append({g2_counter.get(k,0) for k in all_keys})
    freqs_df = pd.DataFrame(freqs)
    chi2, p, dof, expected = chi2_contingency(freqs_df,correction=True)
    return float(p)

def bootstrap_t_test(treatment, control, nboot=100000, stat_fn=np.mean):
    treatment = np.array(treatment)
    control = np.array(control)
    treatment_len = len(treatment)
    tstat = abs(stat_fn(treatment)-stat_fn(control))
    # adjust mean
    treatment = treatment - np.mean(treatment) + np.mean(control)
    # concatenate and scramble
    Z = np.concatenate((treatment,control))
    np.random.shuffle(Z)
    tboot = np.zeros(nboot)
    idx = np.random.randint(0,len(Z),(nboot,len(Z)))
    for i in xrange(nboot):
        sboot = Z[idx[i]]
        tboot[i] = stat_fn(sboot[:treatment_len]) - stat_fn(sboot[treatment_len:])
    pvalue = np.sum(abs(tboot)>=tstat) / float(nboot)
    return pvalue

def bootstrap_mean(vals, nboot=1000):
    vals = np.array(vals)
    tboot = np.zeros(nboot)
    for i in xrange(nboot):
        sboot = sample(vals)
        tboot[i] = np.mean(sboot)
    return np.mean(tboot)

def chooseBestModel(patterns, components, covar_type, gamma_shape=1.0, gamma_inversescale=0.01):
    # Cluster
    models_scores = {}
    for i in range(40):
        g = DPGMM_SCALE(n_components=components,
                        gamma_shape=gamma_shape,
                        gamma_inversescale=gamma_inversescale,
                        covariance_type=covar_type,
                        tol=1e-4,
                        n_iter=1000,
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


def parseRawDataset(data_csv, master_csv, tracer='AV45', ref_key='WHOLECEREB', bilateral=True):
    df = pd.read_csv(data_csv)
    df.loc[:,'subject'] = df.loc[:,'subject'].apply(lambda x: int(x.split('-')[-1]))
    value_df = pd.pivot_table(df, values='pvcval', index=['subject','timepoint'], columns='name')
    size_df = pd.pivot_table(df, values='groupsize', index=['subject','timepoint'], columns='name')
    value_df.columns = [_.upper().replace('-','_') for _ in value_df.columns]
    size_df.columns = ["%s_SIZE" % (_.upper().replace('-','_'),) for _ in size_df.columns]
    original_keys = value_df.columns
    pivot_df = value_df.merge(size_df,left_index=True,right_index=True)

    # drop subjects with no timepoint-specific baseline
    valid_subjects = set([subj for subj,tp in pivot_df.index.values if tp == 'BL'])
    valid_indices = [(subj,tp) for subj,tp in pivot_df.index.values if subj in valid_subjects]
    pivot_df = pivot_df.loc[valid_indices,:]

    # weighted average into bilateral regions
    if bilateral:
        joins = defaultdict(list)
        for o_k in original_keys:
            new_k = o_k.replace('LH_','').replace('RH_','').replace('RIGHT_','').replace('LEFT_','')
            joins[new_k].append(o_k)
        for k,v in joins.items():
            if len(v) > 1:
                weighted = pd.DataFrame([pivot_df[o_k]*pivot_df['%s_SIZE' % (o_k,)]  for o_k in v]).T.sum(axis=1)
                weight_sum = pivot_df[["%s_SIZE" % o_k for o_k in v]].sum(axis=1)
                weighted_avg = weighted/weight_sum
                weighted_avg = pd.DataFrame(weighted_avg, columns=[k])
                pivot_df = pivot_df.merge(weighted_avg, left_index=True, right_index=True)
        pivot_df = pivot_df[joins.keys()]

    # convert to SUVR
    if ref_key not in pivot_df.columns:
        raise Exception("Ref Key %s not found in dataset" % ref_key)
    pivot_df = pivot_df.divide(pivot_df[ref_key],axis='index')

    # Import master csv
    master_df = pd.read_csv(master_csv, low_memory=False, header=[0,1])
    master_df.columns = master_df.columns.get_level_values(1)
    master_df.set_index('RID',inplace=True)

    # Get Diagnoses
    if tracer == 'AV45':
        diag_keys = {'BL': 'Diag@AV45_long',
                     'Scan2': 'Diag@AV45_2_long',
                     'Scan3': 'Diag@AV45_3_long'}
    elif tracer == 'AV1451':
        diag_keys = {'BL': 'Diag@AV1451',
                     'Scan2': 'Diag@AV45_2',
                     'Scan3': 'Diag@AV45_3'}
    else:
        raise Exception("Diag keys for tracer %s not specified" % tracer)
    diags = []
    for subj, tp in pivot_df.index.values:
        diags.append({'subject': subj, 'timepoint': tp, 'diag': master_df.loc[subj,diag_keys[tp]]})
    diags_df = pd.DataFrame(diags).set_index(['subject','timepoint'])

    # separate patterns/uptakes/lobes
    if len(set(pattern_keys) - set(pivot_df.columns)) > 0:
        raise Exception("Some pattern keys are unavailable: %s" % (set(pattern_keys) - set(pivot_df.columns),))
    uptake_df = pivot_df[pattern_keys].copy()
    pattern_df = uptake_df.divide(uptake_df.sum(axis=1),axis=0)
    if len(set(lobe_keys) - set(pivot_df.columns)) > 0:
        raise Exception("Some lobe keys are unavailable: %s" % (set(lobe_keys) - set(pivot_df.columns),))
    lobes_df = pivot_df[lobe_keys].copy()

    # split by timepoint
    pattern_bl_df = pattern_df.loc[(slice(None),'BL'),:].reset_index(level=1,drop=True)
    pattern_scan2_df = pattern_df.loc[(slice(None),'Scan2'),:].reset_index(level=1,drop=True)
    pattern_scan3_df = pattern_df.loc[(slice(None),'Scan3'),:].reset_index(level=1,drop=True)
    uptake_bl_df = uptake_df.loc[(slice(None),'BL'),:].reset_index(level=1,drop=True)
    uptake_scan2_df = uptake_df.loc[(slice(None),'Scan2'),:].reset_index(level=1,drop=True)
    uptake_scan3_df = uptake_df.loc[(slice(None),'Scan3'),:].reset_index(level=1,drop=True)
    lobes_bl_df = lobes_df.loc[(slice(None),'BL'),:].reset_index(level=1,drop=True)
    lobes_scan2_df = lobes_df.loc[(slice(None),'Scan2'),:].reset_index(level=1,drop=True)
    lobes_scan3_df = lobes_df.loc[(slice(None),'Scan3'),:].reset_index(level=1,drop=True)

    # Get years between scans
    yrs = {}
    # for subj in pattern_scan3_df.index:
    #     yrs[subj] = master_df.loc[subj,'%s_1_3_Diff' % tracer]
    # for subj in list(set(pattern_scan2_df.index) - set(pattern_scan3_df.index)):
    #     yrs[subj] = master_df.loc[subj,'%s_1_2_Diff' % tracer]
    for subj in pattern_scan2_df.index:
        yrs[subj] = master_df.loc[subj,'%s_1_2_Diff' % tracer]
    yr_df = pd.DataFrame(yrs.items(), columns=['subject','yrs']).set_index('subject')

    # join patterns/lobes/uptakes into prior/post
    pattern_prior_df = pattern_bl_df.copy()
    lobes_prior_df = lobes_bl_df.copy()
    uptake_prior_df = uptake_bl_df.copy()

    # pattern_post_df = pd.concat((pattern_scan3_df,pattern_scan2_df.loc[list(set(pattern_scan2_df.index) - set(pattern_scan3_df.index)),:]), axis=0)
    # lobes_post_df = pd.concat((lobes_scan3_df,lobes_scan2_df.loc[list(set(lobes_scan2_df.index) - set(lobes_scan3_df.index)),:]), axis=0)
    # uptake_post_df = pd.concat((uptake_scan3_df,uptake_scan2_df.loc[list(set(uptake_scan2_df.index) - set(uptake_scan3_df.index)),:]), axis=0)
    pattern_post_df = pattern_scan2_df
    lobes_post_df = lobes_scan2_df
    uptake_post_df = uptake_scan2_df

    # calculate uptake change and lobe change
    uptake_change = uptake_post_df.subtract(uptake_prior_df.loc[uptake_post_df.index,:]).divide(yr_df.yrs, axis='index')
    lobes_change = lobes_post_df.subtract(lobes_prior_df.loc[lobes_post_df.index,:]).divide(yr_df.yrs, axis='index')

    # create result df
    results = []
    for subj, rows in pivot_df.groupby(level=0):
        try:
            summary_prior = rows.loc[(subj,'BL'),'COMPOSITE']
            diag_prior = diags_df.loc[(subj,'BL'),'diag']
        except Exception as e:
            print "%s has no BL row" % subj
            continue
        indices = set(rows.index.values)
        # if (subj,'Scan3') in indices:
        #     yr_diff = yr_df.loc[subj,'yrs']
        #     summary_post = rows.loc[(subj,'Scan3'),'COMPOSITE']
        #     diag_post = diags_df.loc[(subj,'Scan3'),'diag']
        if (subj,'Scan2') in indices:
            yr_diff = yr_df.loc[subj,'yrs']
            summary_post = rows.loc[(subj,'Scan2'),'COMPOSITE']
            diag_post = diags_df.loc[(subj,'Scan2'),'diag']
        else:
            yr_diff = None
            summary_post = None
            diag_post = None
        summary_change = None
        if summary_post is not None:
            summary_change = (summary_post-summary_prior)/yr_diff
        to_add = {'CORTICAL_SUMMARY_post': summary_post,
                  'CORTICAL_SUMMARY_prior': summary_prior,
                  'diag_prior': diag_prior,
                  'diag_post': diag_post,
                  'CORTICAL_SUMMARY_change': summary_change,
                  'subject': subj}
        results.append(to_add)
    result_df = pd.DataFrame(results).set_index('subject')

    # create return object
    data = {'pattern_bl_df': pattern_bl_df,
            'pattern_scan2_df': pattern_scan2_df,
            'pattern_scan3_df': pattern_scan3_df,
            'pattern_prior_df': pattern_prior_df,
            'pattern_post_df': pattern_post_df,
            'uptake_prior_df': uptake_prior_df,
            'uptake_post_df': uptake_post_df,
            'lobes_prior_df': lobes_prior_df,
            'lobes_post_df': lobes_post_df,
            'change_df': uptake_change,
            'lobes_change_df': lobes_change,
            'result_df': result_df}
    to_return = {}
    for k,v in data.iteritems():
        if 'prior' in k:
            v['timepoint'] = 'prior'
            v.set_index('timepoint', append=True, inplace=True)
        elif 'post' in k:
            v['timepoint'] = 'post'
            v.set_index('timepoint', append=True, inplace=True)
        to_return[k] = v
    return to_return


def smallGroups(result_df, threshold=50):
    # determine small groups
    groups = np.array(list(set(result_df.membership_prior) | set(result_df.membership_post)))
    groups = groups[~np.isnan(groups)]
    small_groups = []
    for g in groups:
        prior_members = set(result_df[result_df.membership_prior==g].index)
        post_members = set(result_df[result_df.membership_post==g].index)
        allmembers = len(prior_members | post_members)
        if allmembers <= threshold:
            small_groups.append(g)
    return small_groups

def bigGroups(result_df, threshold=10):
    # determine big groups
    groups = np.array(list(set(result_df.membership_prior) | set(result_df.membership_post)))
    groups = groups[~np.isnan(groups)]
    big_groups = []
    for g in groups:
        prior_members = set(result_df[result_df.membership_prior==g].index)
        post_members = set(result_df[result_df.membership_post==g].index)
        allmembers = len(prior_members | post_members)
        print "%s: %s" % (g, allmembers)
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
    return (mu_1, sigma_1, cdf_val)

def graphNormalFits(groups, result_df, threshold):
    x = np.linspace(0,3.5,1000)
    plt.figure(1)
    for g in groups:
        members_prior = result_df[result_df.membership_prior==g]
        members_post = result_df[result_df.membership_post==g]
        cortical_summary_values = np.array(list(members_prior['CORTICAL_SUMMARY_prior']) + list(members_post['CORTICAL_SUMMARY_post']))
        mu, sigma, cdf_val = fitNormalCdf(cortical_summary_values, threshold)
        y = norm.pdf(x, loc=mu, scale=sigma)
        plt.plot(x,y)
    plt.plot([threshold, threshold],[0,10])
    plt.show()



def parseConversions(groups, result_df, threshold, master_keys):
    # group diagnostic/pattern/positivity conversions
    diags = ['N', 'SMC', 'MCI', 'AD']
    conversions = {}
    nonzero_conversions = set()
    for g in groups:
        cur_conversions = {}
        conversions[g] = {}
        members_prior = result_df[result_df.membership_prior==g]
        members_post = result_df[result_df.membership_post==g]

        # add membership counts
        cur_conversions['count_prior'] = len(members_prior.index)
        cur_conversions['count_post'] = len(members_post.index)

        # add cortical summary values
        cortical_summary_values = np.array(list(members_prior['CORTICAL_SUMMARY_prior']) + list(members_post['CORTICAL_SUMMARY_post']))
        cur_conversions['cortical_summary'] = cortical_summary_values.mean()
        
        # calculate below threshold with fit
        # mu, sigma, cdf_val = fitNormalCdf(cortical_summary_values, threshold)
        # cur_conversions['below_threshold'] = cdf_val
        
        # calculate raw below threshold percent
        cs_below_threshold = cortical_summary_values[cortical_summary_values <= threshold]
        cur_conversions['below_threshold'] = float(len(cs_below_threshold)) / float(len(cortical_summary_values))

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
        nonzero_conversions = nonzero_conversions | set(diag_conversions.keys())
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
        
        # # add master field avgs
        # master_rows = result_df.loc[members_prior.index,master_keys]
        # master_counts = {"%s_N" % k: v for k,v in dict(master_rows.count()).iteritems()}
        # cur_conversions.update(master_counts)
        # cur_conversions.update(dict(master_rows.mean()))

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

    # CHOOSE HEADERS TO OUTPUT
    header_order = []
    header_order += ['count_prior', 'count_post']
    header_order += ['cortical_summary', 'below_threshold', 'cortical_summary_change']
    header_order += ['N', 'SMC', 'MCI', 'AD']
    header_order += ['neg-neg', 'neg-pos', 'pos-neg', 'pos-pos']
    header_order += sorted(list(nonzero_conversions))
    # for mk in master_keys:
    #     header_order.append(mk)
    #     header_order.append("%s_N" % mk)

    conversions = pd.DataFrame(conversions).T
    conversions.index.name = 'pattern'
    conversions = conversions[header_order]
    return conversions

def saveAparcs(alpha, components, groups, pattern_members, uptake_members, change_members, lut_file):
    if bilateral:
        index_lookup = bilateralTranslations(lut_file)
    else:
        lut_table = importFreesurferLookup(lut_file)
        index_lookup = {v.replace('-','_').upper():[k] for k,v in lut_table.iteritems()}
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

def saveLobeAparcs(alpha, groups, lobe_tp, lut_file):
    if bilateral:
        index_lookup = bilateralTranslations(lut_file)
    else:
        index_lookup = importFreesurferLookup(lut_file)

    all_lobe_members = lobePatterns(lobe_tp, groups, pattern=True)
    aparc_input_template = "../output/fake_aparc_inputs/dpgmm_alpha%s_group_%s_%s"
    for g in groups:
        lobe_members = all_lobe_members[g]
        out_file = aparc_input_template % (alpha,g,'lobes')
        lobe_values = lobe_members[lobe_keys].mean().to_dict()
        saveFakeAparcInput(out_file, lobe_values, index_lookup)


def mergeResults(result_df, pattern_prior_df, pattern_post_df, uptake_prior_df, uptake_post_df, lobe_prior_df, lobe_post_df, rchange_df):
    membership_prior = result_df[['membership_prior','CORTICAL_SUMMARY_prior']]
    membership_prior.columns = ['membership','CORTICAL_SUMMARY']
    membership_prior.loc[:,'timepoint'] = 'prior'
    membership_post = result_df[['membership_post','CORTICAL_SUMMARY_post']]
    membership_post.columns = ['membership','CORTICAL_SUMMARY']
    membership_post.loc[:,'timepoint'] = 'post'

    # changes (add 'prior' timepoint)
    change_members = rchange_df.copy().merge(membership_prior[['membership','timepoint']], left_index=True, right_index=True)
    change_members.index.name = 'rid'
    change_members.set_index('timepoint',append=True,inplace=True)

    # memberships
    membership_df = pd.concat((membership_prior,membership_post)).set_index('timepoint', append=True)

    # uptakes
    uptake_members = pd.concat((uptake_prior_df, uptake_post_df))
    uptake_members = uptake_members.merge(membership_df, how='left', left_index=True, right_index=True)

    # patterns
    pattern_members = pd.concat((pattern_prior_df, pattern_post_df))
    pattern_members = pattern_members.merge(membership_df, how='left', left_index=True, right_index=True)

    # lobes
    lobe_members = pd.concat((lobe_prior_df, lobe_post_df))
    lobe_members = lobe_members.merge(membership_df, how='left', left_index=True, right_index=True)

    return (uptake_members, pattern_members, lobe_members, change_members)


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

def plotValueScatter(result_df, groups, xys, fit_reg=False, test=True, show=False):  
    figures = []
    for i, (x,y) in enumerate(xys):
        members = []
        valid_groups = []
        for grp in groups:
            grp_members = result_df[result_df.membership_prior==grp].dropna(subset=[x,y])
            if len(grp_members) < 2:
                print "Not enough members in %s with values" % grp
            else:
                valid_groups.append(grp)
                members.append(grp_members)

        if len(valid_groups) == 0:
            print "No valid groups"
            continue

        group_result_df = result_df[result_df['membership_prior'].isin(valid_groups)]
        xmin = group_result_df[x].min()
        xmax = group_result_df[x].max()
        xranges = xmax-xmin
        xmin -= 0.05*xranges
        xmax += 0.05*xranges
        ymin = group_result_df[y].min()
        ymax = group_result_df[y].max()
        yranges = ymax-ymin
        ymin -= 0.05*yranges
        ymax += 0.05*yranges

        newplot = sns.lmplot(x=x,y=y,data=group_result_df,hue='membership_prior', fit_reg=fit_reg)
        newplot.set(ylim=(ymin,ymax), xlim=(xmin,xmax))

        if test:
            if len(members) != 2:
                print "can only test 2 groups"
            else:
                members1 = members[0]
                members2 = members[1]
                x1 = list(members1[x])
                x2 = list(members2[x])
                y1 = list(members1[y])
                y2 = list(members2[y])
                slope_z, slope_p, int_z, int_p = LinRegZTest(x1, y1, x2, y2)
                label = "SLOPE_P: %s, INT_P: %s" % (round(slope_p,5), round(int_p,5))
                if slope_p <= 0.05 or int_p <= 0.05:
                    print "%s vs %s: %s" % (x,y,label)
                newplot.fig.text(0.33, 1.02, label, fontsize=16)

        figures.append(newplot)
    if show:
        sns.plt.show()

def plotValueDensity(result_df, groups, value_key):
    long_df = pd.melt(result_df, id_vars=['membership_prior'], value_vars=value_key)
    plt.figure(1)
    plt.title("Cortical Summary SUVR Density, by Pattern Group (n>3)")
    groups_members = [long_df[long_df.membership_prior==g] for g in groups]
    groups_means = [_['value'].mean() for _ in groups_members]
    means_members = zip(groups_means, groups_members, groups)
    means_members_sorted = sorted(means_members, key=lambda x: x[0])
    for mean_val, members, g in means_members_sorted:
        sns.kdeplot(members)
        members['value'].plot(kind='kde', label="Pattern %s (n=%s, SUVR=%s)" % (int(g),int(len(members.index)),mean_val), alpha=0.8)
    plt.xlabel
    plt.legend()
    plt.show()


def plotValueDensity(result_df, groups, value_key):
    long_df = pd.melt(result_df, id_vars=['membership_prior'], value_vars=value_key)
    plt.figure(1)
    plt.title("Cortical Summary SUVR Density, by Pattern Group (n>3)")
    groups_members = [long_df[long_df.membership_prior==g] for g in groups]
    groups_means = [_['value'].mean() for _ in groups_members]
    means_members = zip(groups_means, groups_members, groups)
    means_members_sorted = sorted(means_members, key=lambda x: x[0])
    for mean_val, members, g in means_members_sorted:
        sns.kdeplot(members['value'], label="Pattern %s (n=%s, SUVR=%s)" % (int(g),int(len(members.index)),mean_val))
    plt.xlabel('Cortical Summary SUVR (partial volume corrected, whole cerebellum reference region)')
    plt.ylabel('Density')
    plt.legend()
    plt.show()


def plotValueBox(result_df, groups, value_key, save=False):
    long_df = pd.melt(result_df, id_vars=['membership_prior'], value_vars=value_key)
    by_group = pd.DataFrame()
    for g in groups:
        members = long_df[long_df.membership_prior==g][['value']]
        members['pattern'] = g
        by_group = pd.concat((by_group,members))
    dfg = by_group.groupby('pattern')
    labels = ['%s\n$n$=%d'%(k, len(v)) for k, v in dfg]
    counts = [len(v) for k,v in dfg]
    total = float(sum(counts))
    widths = [2*c/total for c in counts]
    bplot = sns.violinplot(x='pattern',y='value',data=by_group)
    bplot.set_xticklabels(labels, rotation=45)
    plt.title("%s" % value_key)
    plt.suptitle("")
    if save:
        plt.savefig('../boxplot_%s.png' % value_key.replace('/','_'), dpi=400)
        plt.close()
    else:
        plt.show()
    

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
    #blacklist = ['WHITE_MATTER','BRAIN_STEM','VENTRALDC']
    blacklist = []
    filtered_keys = [k for k in keys if not any([bl in k for bl in blacklist])]
    ranks = pd.DataFrame()
    for g in groups:
        ranked_regions = members[members.membership==g].mean().loc[filtered_keys].sort(ascending=False,inplace=False)
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

def diagnosisChi2Test(pairs, diags):
    diags_comparisons = []
    for g1, g2 in pairs:
        cont_table = diags.loc[[g1,g2]]
        chi2, p, dof, expected = chi2_contingency(cont_table,correction=True)
        diags_comparisons.append({'GROUP1': g1,
                                  'GROUP2': g2,
                                  'DIAGS_PVALUE': p})
    diags_comparisons_df = pd.DataFrame(diags_comparisons)
    return diags_comparisons_df

def generateFileRoot(alpha, model=False):
    print alpha
    file_root = '../dpgmm_alpha%s' % round(alpha,2)
    if bilateral:
        file_root += '_bilateral'
    if normals_only and not model:
        file_root += '_normals'
    return file_root

def loadResults(alpha, master_csv):
    master_df = pd.read_csv(master_csv, low_memory=False, header=[0,1])
    master_df.columns = master_df.columns.get_level_values(1)
    master_df.set_index('RID', inplace=True)
    master_df = master_df.loc[:,master_keys]
    result_df = pd.read_csv('%s_result.csv' % generateFileRoot(alpha, model=True))
    result_df.set_index('subject', inplace=True)
    result_df.sort_index(inplace=True)
    result_df.drop(['membership_prior','membership_post'], axis=1, inplace=True)
    result_df.rename(columns={'membership_prior_conf':'membership_prior','membership_post_conf':'membership_post'},inplace=True)
    result_df = result_df.merge(master_df,left_index=True,right_index=True,how='left')
    if normals_only:
        result_df = result_df[result_df.diag_prior.isin(['N','SMC'])]
    return result_df

def scaleInput(pattern_prior_df, pattern_post_df):
    scale_type = 'original'
    patterns_only, scaler = scaleRawInput(pattern_prior_df, scale_type=scale_type)
    post_patterns_only = pd.DataFrame(scaler.transform(pattern_post_df))
    post_patterns_only.set_index(pattern_post_df.index, inplace=True)
    return (patterns_only, post_patterns_only)

def trainDPGMM(pattern_prior_df, pattern_post_df, result_df, covar_type):
    # Scale inputs
    patterns_only, post_patterns_only = scaleInput(pattern_prior_df, pattern_post_df)

    # Choose alpha + do model selection
    components = len(patterns_only.index)
    best_model = chooseBestModel(patterns_only, components, covar_type, gamma_shape=1.0, gamma_inversescale=0.000000001)
    alpha = best_model.alpha
    with open("%s_%s_model.pkl" % (generateFileRoot(alpha, model=True),covar_type), 'wb') as fid:
        cPickle.dump(best_model, fid)   

def generateResults(model, pattern_prior_df, pattern_post_df, result_df, conf_filter=False):
    alpha = model.alpha

    # Scale inputs
    patterns_only, post_patterns_only = scaleInput(pattern_prior_df, pattern_post_df)

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

    # Predict pattern groups and add to result df
    y_df = pd.DataFrame(model.predict(patterns_only))
    y_df.columns = ['group']
    y_df.set_index(patterns_only.index, inplace=True)
    y_proba_df = pd.DataFrame(model.predict_proba(patterns_only)).set_index(patterns_only.index).xs('prior',level='timepoint')
    probsums = y_proba_df.sum()
    components_distr_max = y_proba_df[probsums[probsums >= 0.1].index].max(axis=1)
    if conf_filter:
        confident_assignments_prior = components_distr_max[components_distr_max>=membership_conf].index
    else:
        confident_assignments_prior = components_distr_max.index

    y_df_post = pd.DataFrame(model.predict(post_patterns_only))
    y_df_post.columns = ['group']
    y_df_post.set_index(post_patterns_only.index, inplace=True)
    y_proba_df_post = pd.DataFrame(model.predict_proba(post_patterns_only)).set_index(post_patterns_only.index).xs('post',level='timepoint')
    probsums = y_proba_df_post.sum()
    components_distr_max = y_proba_df_post[probsums[probsums > 0.1].index].max(axis=1)
    if conf_filter:
        confident_assignments_post = components_distr_max[components_distr_max>=membership_conf].index
    else:
        confident_assignments_post = components_distr_max.index
    for rid in result_df.index:
        result_df.loc[rid,'membership_prior'] = y_df.loc[(rid,'prior'),'group']
        # for pca_col in pca_patterns_only.columns:
        #     result_df.loc[rid,'pca_%s_prior' % pca_col] = pca_patterns_only.loc[(rid, 'prior'),pca_col]
        if (rid, 'post') in y_df_post.index:
            result_df.loc[rid,'membership_post'] = y_df_post.loc[(rid,'post'),'group']
            # for pca_col in pca_post_patterns_only.columns:
            #     result_df.loc[rid,'pca_%s_post' % pca_col] = pca_post_patterns_only.loc[(rid, 'post'),pca_col]
        if rid in confident_assignments_prior:
            result_df.loc[rid,'membership_prior_conf'] = result_df.loc[rid,'membership_prior']
        if rid in confident_assignments_post:
            result_df.loc[rid,'membership_post_conf'] = result_df.loc[rid,'membership_post']

    # save
    result_df.to_csv("%s_result.csv" % generateFileRoot(alpha, model=True))
    return result_df

def pcaCorrelations(result_df):
    # pca component correlations
    y_cols = ['CORTICAL_SUMMARY_prior',
              'CORTICAL_SUMMARY_change',
              'UW_MEM_BL_3months',
              'UW_EF_BL_3months',
              'CSF_TAU_closest_AV45',
              'CSF_ABETA_closest_AV45',
              'FSX_HC/ICV_BL_3months']
    x_cols = ['pca_prior','CORTICAL_SUMMARY_prior']
    props={'pad':10}
    for y in y_cols:
        fig, axes = plt.subplots(nrows=1, ncols=2)
        # pca
        pca_x = ['pca_%s_prior' % i for i in range(5)]
        model = ols(x=result_df[pca_x],y=result_df[y])
        x = model.beta.x
        intercept = model.beta.intercept
        r2 = '{0:.3f}'.format(model.r2)
        equation = 'Y ~ %sx + %s' % ('{0:.3f}'.format(x),'{0:.3f}'.format(intercept))
        ax = axes[0]
        result_df.plot(ax=ax, kind='scatter',x='pca_0_prior',y=y)
        xvals = np.array(sorted(list(result_df['pca_0_prior'])))
        ax.plot(xvals, xvals*x + intercept, 'r')
        ax.set_title('Principal Component')
        ax.annotate('%s\nR2=%s' % (equation,r2), xy=(0.05,0.95), xycoords='axes fraction', bbox=props, verticalalignment='top', fontsize=15, color='k')
        # cortical summary
        model = ols(x=result_df['CORTICAL_SUMMARY_prior'],y=result_df[y])
        x = model.beta.x
        intercept = model.beta.intercept
        r2 = '{0:.3f}'.format(model.r2)
        equation = 'Y ~ %sx + %s' % ('{0:.3f}'.format(x),'{0:.3f}'.format(intercept))
        ax = axes[1]
        result_df.plot(ax=ax, kind='scatter',x='CORTICAL_SUMMARY_prior',y=y)
        xvals = np.array(sorted(list(result_df['CORTICAL_SUMMARY_prior'])))
        ax.plot(xvals, xvals*x + intercept, 'r')
        ax.set_title('Cortical Summary SUVR')
        ax.annotate('%s\nR2=%s' % (equation,r2), xy=(0.05,0.95), xycoords='axes fraction', bbox=props, verticalalignment='top', fontsize=15, color='k')

def bootstrap_regression(x, y, nboot=1000):
    '''
    1. Estimate regression coefficients on original sample
    2. Calculate fitted value and residual for each observation
    3. Bootstrap sample the residuals
    4. Calculate bootstrapped Y values (Yfit + bootstrapped residual)
    5. Regress bootstrapped Y values on fixed X values to get bootstrap regression coefficients
    6. Use bootstrapped regression coefficients to get SE and CI
    '''
    # regression on original sample
    m_og, b_og, r_value_og, p_value_og, std_err_og = linregress(x, y)
    regfn_og = np.poly1d([m_og, b_og])
    # residuals for each observation
    residuals_og = []
    for x_i, y_i in zip(x,y):
        pass

    vals = np.array(vals)
    tboot = np.zeros(nboot)
    for i in xrange(nboot):
        sboot = sample(vals)
        tboot[i] = np.mean(sboot)
    return np.mean(tboot)


def LinRegZTest(x1, y1, x2, y2):
    m1, b1, r_value1, p_value1, std_err1 = linregress(x1, y1)
    m2, b2, r_value2, p_value2, std_err2 = linregress(x2, y2)
    slope_z, slope_pvalue = regCoeffZTest(m1, m2, std_err1, std_err2)
    int_z, int_pvalue = regCoeffZTest(b1, b2, std_err1, std_err2)
    return (slope_z, slope_pvalue, int_z, int_pvalue)

def lobePatterns(lobe_tp, groups, pattern=True):
    prior_pts = lobe_tp.xs('prior',level='timepoint')
    lobe_patterns = {}
    for i,g in enumerate(groups):
        members = prior_pts[prior_pts.membership_prior == g]
        members = members[lobe_keys]
        if pattern:
            members = members.divide(members.sum(axis=1),axis=0)
        members.reset_index(inplace=True)
        members['membership'] = g
        lobe_patterns[g] = members
    return lobe_patterns

def graphLobes(lobe_tp, groups, suvrs, pattern=True, save=True):
    y_limits = (0,3.0)
    if pattern:
        y_limits = (0,0.2)
    prior_pts = lobe_tp.xs('prior',level='timepoint')
    lobe_patterns = lobePatterns(lobe_tp, groups, pattern=pattern)
    for i,g in enumerate(groups):
        fig = plt.figure(i+1)
        #fig.set_size_inches(5.5,4)
        #fig.set_dpi(300)
        members = lobe_patterns[g]
        long_df = pd.melt(members, id_vars=['rid'], value_vars=lobe_keys)
        bplot = sns.violinplot(x='variable', y='value', data=long_df, size=4, aspect=2)
        bplot.set(ylim=y_limits)
        plt.xticks(rotation=80)
        title = 'Group %s, n=%s, mean SUVR: %s' % (int(g), len(members.index), suvrs[g])
        plt.title(title)
        plt.suptitle("")
        if save:
            plt.savefig('../pattern_%s_lobeplot.png' % int(g), bbox_inches='tight')
    plt.show()

def graphRegionPatterns(pattern_tp, groups, suvrs, save=True):
    y_limits = (0,0.07)
    prior_pts = pattern_tp.xs('prior',level='timepoint')
    for i,g in enumerate(groups):
        fig = plt.figure(i+1)
        members = prior_pts[prior_pts.membership_prior == g]
        long_df = pd.melt(members, id_vars=['rid'], value_vars=pattern_keys)
        bplot = sns.violinplot(x='variable', y='value', data=long_df, scale='width')
        bplot.set(ylim=y_limits)
        plt.xticks(rotation=80)
        title = 'Group %s, n=%s, mean SUVR: %s' % (int(g), len(members.index), suvrs[g])
        plt.title(title)
        plt.suptitle("")
        if save:
            plt.savefig('../pattern_%s_regionplot.png' % int(g), bbox_inches='tight')
    plt.show()

def testLobePatternDifferences(lobe_tp, groups, pattern=True):
    lobe_patterns = lobePatterns(lobe_tp, groups, pattern=pattern)
    full_results = {}
    for g1,g2 in itertools.combinations(groups,2):
        df = pd.concat((lobe_patterns[g1],lobe_patterns[g2]))
        res = compareTwoGroups(g1, g2, df, lobe_keys, adjust=False)
        full_results[(g1,g2)] = res
    return full_results


def saveAV45LMEDataset(model, master_csv, pattern_prior_df, result_df, confident=False):
    # Scale inputs
    bl_patterns_only, scaler = scaleRawInput(pattern_prior_df, scale_type='original')

    # get membership prob
    proba_df_bl = pd.DataFrame(model.predict_proba(bl_patterns_only)).set_index(bl_patterns_only.index)

    # groups
    groups_df_bl = pd.DataFrame(model.predict(bl_patterns_only)).set_index(bl_patterns_only.index)
    groups_df_bl.columns = ['group']

    # Drop patterns that don't have >1 members
    member_counts_bl = proba_df_bl.sum()
    valid_bl = set(member_counts_bl[member_counts_bl>1].index)
    proba_df_bl = proba_df_bl[list(valid_bl)]

    # get Max confidences (group assignments)
    proba_max_bl = proba_df_bl.max(axis=1)
    confident_bl = proba_max_bl[proba_max_bl>=0.5].index

    # Round off membership probs
    proba_df_bl = proba_df_bl.applymap(lambda x: round(x,3))

    if confident:
        proba_df_bl = proba_df_bl.loc[confident_bl,:]

    # Merge in patterns
    proba_df_bl = proba_df_bl.merge(groups_df_bl, left_index=True, right_index=True, how='left')
    proba_df_bl.index = proba_df_bl.index.droplevel(1)

    # Get master df
    master_df = pd.read_csv(master_csv, low_memory=False, header=[0,1])
    master_df.columns = master_df.columns.get_level_values(1)
    master_df.set_index('RID', inplace=True)

    # Add demographic data (baseline age, sex, e2 status, e4 status, education)
    demo_df_bl = master_df.loc[:,['Age@AV45','Gender','APOE2_BIN','APOE4_BIN','Edu.(Yrs)','Diag@AV45_long']].copy()
    demo_df_bl.Gender = demo_df_bl.Gender-1
    proba_df_bl = proba_df_bl.merge(demo_df_bl, left_index=True, right_index=True, how='left')

    # add slopes
    slopes = result_df[['CORTICAL_SUMMARY_change', 'CORTICAL_SUMMARY_prior']]
    full_df = proba_df_bl.merge(slopes, left_index=True, right_index=True, how='left')

    # Rename columns
    col_translations = {'Diag@AV45_long':'diag_prior',
                        'Age@AV45':'Age.AV45',
                        'CORTICAL_SUMMARY_change':'CORTICAL_SUMMARY_slope'}
    full_df.rename(columns=col_translations,inplace=True)
    full_df.index.name='RID'
    full_df.dropna(inplace=True)
    print full_df.shape

    # drops rows with no pattern membership
    member_sums = full_df[list(valid_bl)].sum(axis=1)
    notmember = set(member_sums[member_sums==0].index)
    full_df = full_df.drop(notmember) 
    all_output_file = '%s_AV45_ALL_longdata_slope.csv' % (generateFileRoot(model.alpha))
    full_df.to_csv(all_output_file)


def saveLMEDataset(model, master_csv, pattern_prior_df, result_df, dep_val_var, dep_time_var, categorize=False, confident=False):
    '''
    If categorize, encode pattern membership by one-hot, where membership is determined by argmax(membership_probs)
    '''
    # Scale inputs
    prior_patterns_only, scaler = scaleRawInput(pattern_prior_df, scale_type='original')

    # Get membership probabilities
    proba_df = pd.DataFrame(model.predict_proba(prior_patterns_only)).set_index(prior_patterns_only.index).xs('prior',level='timepoint')
    groups_df = pd.DataFrame(model.predict(prior_patterns_only)).set_index(prior_patterns_only.index).xs('prior',level='timepoint')
    groups_df.columns = ['group']
    # Remove components with insignificant membership
    probsums = proba_df.sum()
    valid_components = probsums[probsums >= 0.1].index
    proba_df = proba_df[valid_components]
    proba_max_df = proba_df.max(axis=1)
    confident_assignments = proba_max_df[proba_max_df>=0.5].index

    # Convert to categorical 
    if categorize:
        for c in proba_df.columns:
            proba_df[c] = (proba_df[c]==proba_max_df)
        proba_df = proba_df.applymap(lambda x: 1 if x else 0)
        cont_suffix = ''
    else:
        proba_df = proba_df.applymap(lambda x: round(x,3))
        cont_suffix = '_continuous'

    # Filter for confident group assignments
    if confident:
        proba_df = proba_df.loc[confident_assignments,:]

    # Drop patterns that don't have >1 members
    member_counts = proba_df.sum()
    valid_patterns = list(member_counts[member_counts>1].index)
    proba_df = proba_df[valid_patterns]

    # Merge in pattern factor
    proba_df = proba_df.merge(groups_df, left_index=True, right_index=True, how='left')

    # Get master df
    master_df = pd.read_csv(master_csv, low_memory=False, header=[0,1])
    master_df.columns = master_df.columns.get_level_values(1)
    master_df.set_index('RID', inplace=True)


    # Add demographic data (baseline age, sex, e2 status, e4 status, education)
    demo_df = master_df.loc[:,['Age@AV45','Gender','APOE2_BIN','APOE4_BIN','Edu.(Yrs)','Diag@AV45_long']].copy()
    demo_df.Gender = demo_df.Gender-1
    proba_df = proba_df.merge(demo_df, left_index=True, right_index=True, how='left')
    col_translations = {'Diag@AV45_long':'diag_prior',
                        'Diag@AV45_2_long':'diag_prior',
                        'Age@AV45':'Age.AV45',
                        'Age@AV45_2':'Age.AV45'}
    proba_df.rename(columns=col_translations,inplace=True)

    # add in baseline AV45
    proba_df = proba_df.merge(result_df[['CORTICAL_SUMMARY_prior']], left_index=True, right_index=True, how='left')


    # Add dependent variable and test times
    dep_value_keys = ['%s%s' % (dep_val_var,i) for i in range(30) if '%s%s' % (dep_val_var,i) in master_df.columns]
    dep_keypairs = {}
    for dv_key in dep_value_keys:
        dt_key = dv_key.replace(dep_val_var,dep_time_var)
        timepoint_index = int(dv_key.replace(dep_val_var,''))
        dep_keypairs[timepoint_index] = (dt_key, dv_key)
    dep_by_subject = defaultdict(list)
    for tp_index in sorted(dep_keypairs.keys()):
        dt_key, dv_key = dep_keypairs[tp_index]
        tp_df = master_df[[dt_key,dv_key]].dropna()
        for pid, row in tp_df.iterrows():
            if pid not in proba_df.index:
                continue
            tp_tuple = (float(row[dt_key]),float(row[dv_key]))
            dep_by_subject[pid].append(tp_tuple)
    
    # create longitudinal dataset
    full_df = pd.DataFrame()
    for pid, timepoints in dep_by_subject.iteritems():
        timepoints = sorted(timepoints, key=lambda x: x[0])
        base_time = timepoints[0][0]
        for time, val in timepoints:
            adjusted_time = time-base_time
            new_row = proba_df.loc[pid].copy()
            new_row[dep_val_var] = val
            new_row['YrsPostBL'] = round(adjusted_time,1)
            full_df = pd.concat((full_df, new_row),axis=1)

    # create slope dataset
    slope_df = pd.DataFrame()
    for pid, timepoints in dep_by_subject.iteritems():
        timepoints = sorted(timepoints, key=lambda x: x[0])
        base_time = timepoints[0][0]
        adj_timepoints = []
        for time, val in timepoints:
            adjusted_time = time-base_time
            adj_timepoints.append((round(adjusted_time,1),val))
        row_slope = slope(adj_timepoints)
        if row_slope is not None:
            new_row = proba_df.loc[pid].copy()
            new_row['%s_slope' % dep_val_var] = row_slope
            slope_df = pd.concat((slope_df,new_row),axis=1)

    full_df = full_df.T
    full_df.index.name='RID'
    full_df.dropna(inplace=True)
    slope_df = slope_df.T
    slope_df.index.name='RID'
    slope_df.dropna(inplace=True)

    # drops rows with no pattern membership
    member_sums = full_df[valid_patterns].sum(axis=1)
    notmember = set(member_sums[member_sums==0].index)
    full_df = full_df.drop(notmember) 
    member_sums = slope_df[valid_patterns].sum(axis=1)
    notmember = set(member_sums[member_sums==0].index)
    slope_df = slope_df.drop(notmember)

    # write out
    clean_varname = dep_val_var.replace('/','_').replace('.','_').strip()
    all_output_file = '%s_%s_ALL_longdata.csv' % (generateFileRoot(model.alpha),clean_varname)
    slope_output_file = '%s_%s_ALL_slopedata.csv' % (generateFileRoot(model.alpha),clean_varname)
    full_df.to_csv(all_output_file)
    slope_df.to_csv(slope_output_file)

def createContrasts(columns, time_key, out_file=None):
    new_df = pd.DataFrame(columns=columns)
    cols = [_ for _ in columns if 'X' in _ and time_key in _]
    for c in cols:
        new_df.loc[len(new_df)+1]=0
        new_df.loc[len(new_df),c]=1
    for c1,c2 in itertools.combinations(cols,2):
        new_df.loc[len(new_df)+1]=0
        new_df.loc[len(new_df),c1]=1
        new_df.loc[len(new_df),c2]=-1
    if out_file:
        new_df.to_csv(out_file,index=False)

def createLMContrast(columns,out_file=None):
    '''
    create contrasts between X* variables
    '''
    new_df = pd.DataFrame(columns=columns)
    cols = ["X%s" % i for i in range(100) if 'X%s' % i in columns]
    for c in cols:
        new_df.loc[len(new_df)+1]=0
        new_df.loc[len(new_df),c]=1
    for c1,c2 in itertools.combinations(cols,2):
        new_df.loc[len(new_df)+1]=0
        new_df.loc[len(new_df),c1]=1
        new_df.loc[len(new_df),c2]=-1
    if out_file:
        new_df.to_csv(out_file,index=False)

def printMeanLobePatterns(lobe_members,valid_patterns,cortical_summary,out_file):
    valid_means = []
    for i in valid_patterns:
        members = lobe_members[lobe_members.membership==i]
        members = members[lobe_keys]
        member_patterns = members.divide(members.sum(axis=1),axis=0)
        member_patterns_means = member_patterns.mean().to_dict()
        new_row = {'pattern': i,
                   'mean_suvr': cortical_summary[i]}
        new_row.update(member_patterns_means)
        valid_means.append(new_row)
    df = pd.DataFrame(valid_means)
    df.to_csv(out_file,index=False)


if __name__ == '__main__':
    data_csv = '../datasets/pvc_adni_av45/mostregions_output.csv'
    master_csv = '../FDG_AV45_COGdata/FDG_AV45_COGdata_03_23_16.csv'
    data = parseRawDataset(data_csv, master_csv, tracer='AV45', ref_key='WHOLECEREB', bilateral=True)

    pattern_bl_df = data['pattern_bl_df']
    pattern_scan2_df = data['pattern_scan2_df']
    pattern_scan3_df = data['pattern_scan3_df']
    pattern_prior_df = data['pattern_prior_df']
    pattern_post_df = data['pattern_post_df']
    uptake_prior_df = data['uptake_prior_df']
    uptake_post_df = data['uptake_post_df']
    lobes_prior_df = data['lobes_prior_df']
    lobes_post_df = data['lobes_post_df']
    lobes_change_df = data['lobes_change_df']
    result_df = data['result_df']
    rchange_df = data['change_df']
    pattern_col_order = list(pattern_prior_df.columns)

    # # Calculate
    # result_df = trainDPGMM(pattern_prior_df, pattern_post_df, result_df, 'spherical')
    # result_df = trainDPGMM(pattern_prior_df, pattern_post_df, result_df, 'tied')

    # Load
    # model_file = '../dpgmm_alpha14.36_bilateral_diag_model.pkl'
    # model_file = '../dpgmm_alpha13.13_bilateral_tied_model.pkl'
    model_file = '../dpgmm_alpha12.66_bilateral_spherical_model.pkl'

    best_model = cPickle.load(open(model_file, 'rb'))
    alpha = best_model.alpha
    result_df = generateResults(best_model, pattern_prior_df, pattern_post_df, result_df, conf_filter=False)

    # get groups
    big_groups = bigGroups(result_df, threshold=3)
    all_groups = bigGroups(result_df, threshold=2)

    # generate conversion data
    conversions = parseConversions(all_groups, result_df, threshold, master_keys)
    conversions.to_csv('%s_conversions.csv' % generateFileRoot(alpha))
    cortical_summary = conversions.cortical_summary.to_dict()
    positive_patterns = list(conversions[conversions['pos-pos']>=0.8].index)
    negative_patterns = list(conversions[conversions['neg-neg']>=0.9].index)
    transition_patterns = list(set(conversions.index) - (set(positive_patterns) | set(negative_patterns)))
    groups = {'positive': positive_patterns, 'negative': negative_patterns, 'transition': transition_patterns}

    # merge memberships into data
    uptake_members, pattern_members, lobe_members, change_members = mergeResults(result_df, pattern_prior_df, pattern_post_df, uptake_prior_df, uptake_post_df, lobes_prior_df, lobes_post_df, rchange_df)

    # # dump pattern lobe means
    # printMeanLobePatterns(lobe_members,all_groups,cortical_summary,'%s_pattern_means.csv' % generateFileRoot(alpha))

    # # save fake aparcs
    # saveAparcs(round(alpha,2), components, big_groups, pattern_members, uptake_members, change_members, lut_file)
    # saveLobeAparcs(round(alpha,2), big_groups, lobe_tp, lut_file)

    # dump LMM Dataset
    LME_dep_variables = [('UW_MEM_','UW_MEM_postAV45_'),
                         ('UW_EF_','UW_EF_postAV45_'),
                         ('AV45','')]
    for dep_val_var, dep_time_var in LME_dep_variables:
        if dep_val_var == 'AV45':
            saveAV45LMEDataset(best_model, master_csv, pattern_prior_df, result_df, confident=True)
        else:
            saveLMEDataset(best_model, master_csv, pattern_prior_df, result_df, dep_val_var, dep_time_var, categorize=False, confident=True)


    # lobe violin plots
    graphLobes(lobe_tp,big_groups,cortical_summary)

    # lobe pattern differences
    lobe_sigtests = testLobePatternDifferences(lobe_tp, groups, pattern=True)

    # plot densities
    plotValueDensity(result_df, big_groups, 'CORTICAL_SUMMARY_prior')


    # region ranking
    uptake_ranks = regionRanks(big_groups, uptake_members, pattern_keys)
    change_ranks = regionRanks(big_groups, change_members, change_keys)
    pattern_ranks = regionRanks(big_groups, pattern_members, pattern_keys)
    uptake_ranks.T.to_csv('%s_uptake_ranks.csv' % generateFileRoot(alpha))
    change_ranks.T.to_csv('%s_change_ranks.csv' % generateFileRoot(alpha))
    pattern_ranks.T.to_csv('%s_pattern_ranks.csv' % generateFileRoot(alpha))

    # plot pattern flows
    graphNetworkConversions(transition_patterns+positive_patterns, conversions, iterations=50, threshold=0.1, alternate_nodes=transition_patterns)
    graphNetworkConversions(negative_patterns+transition_patterns, conversions, iterations=50, threshold=0.1, alternate_nodes=transition_patterns)
    graphNetworkConversions(negative_patterns+positive_patterns+transition_patterns, conversions, iterations=50, threshold=0.1, alternate_nodes=negative_patterns+transition_patterns)

    # plot group value distributions
    boxplot_keys = ['CORTICAL_SUMMARY_prior', 'CORTICAL_SUMMARY_change']
    boxplot_keys += [k for k in master_keys if k not in cat_keys]
    for k in boxplot_keys:
        plotValueBox(result_df, big_groups, k, save=True)

    # # derive cortical summary prior p values between all pairs
    # # merge them all together
    # merged_members = uptake_members.copy()
    # merged_members = merged_members.merge(change_members, left_index=True, right_index=True)
    # merged_members = merged_members.mergee(lobe_members, left_index=True, right_index=True)
    # all_keys = master_keys+summary_keys+pattern_keys+change_keys
    # pvalue_df = compareGroupPairs(itertools.combinations(big_groups,2), merged_members, ['CORTICAL_SUMMARY_prior'])
    # pvalue_df.to_csv('%s_cortical_summary_pvalues.csv' % generateFileRoot(alpha), index=False)

    # read cortical summary pvalues
    pvalue_df = pd.read_csv('%s_cortical_summary_pvalues.csv' % generateFileRoot(alpha))
    pvalue_df.set_index(['GROUP1','GROUP2'],inplace=True)
    sig_same = pvalue_df[pvalue_df.CORTICAL_SUMMARY_prior_PVALUE>0.05]
    sig_same_pairs = [_[0] for _ in zip(sig_same.index)]

    # diagnosis chi2 contingency test
    diags = conversions[['AD','MCI','SMC','N']]
    diags = diags.multiply(conversions[['count_prior','count_post']].sum(axis=1),axis=0).astype(int)
    diags_comparisons_df = diagnosisChi2Test(sig_same_pairs, diags)

    # compare groups with non-significantly different cortical summary prior distributions
    pairs_comparisons_df = compareGroupPairs(sig_same_pairs, merged_members, all_keys)
    pairs_comparisons_df = pairs_comparisons_df.merge(diags_comparisons_df,on=['GROUP1','GROUP2'])
    pairs_comparisons_df.set_index(['GROUP1','GROUP2'], inplace=True)
    pairs_comparisons_df = pairs_comparisons_df.T
    pairs_comparisons_df.to_csv('%s_pair_comparisons.csv' % generateFileRoot(alpha))
    pvalue_indices = [_ for _ in pairs_comparisons_df.index if 'PVALUE' in _]
    # pvalues only
    pairs_comparisons_df_pvalues = pairs_comparisons_df.loc[pvalue_indices,:]
    pairs_comparisons_df_pvalues.to_csv('%s_pair_pvalues.csv' % generateFileRoot(alpha), na_rep=-1)
    # means only
    reset_df = pairs_comparisons_df.T.reset_index()
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
    mean_averages.to_csv('%s_pair_means.csv' % generateFileRoot(alpha), na_rep='NA')


    # SCATTER PLOTS! (FSX, WMH, FDG, EF, MEM)
    bl_changes = [('FSX_HC/ICV_BL_3months', 'FSX_HC/ICV_slope'),
                  ('UW_EF_BL_3months','UW_EF_slope'),
                  ('UW_MEM_BL_3months','UW_MEM_slope'),
                  ('WMH_percentOfICV_AV45_6MTHS','WMH_percentOfICV_slope'),
                  ('FDG_PONS_AV45_6MTHS','FDG_postAV45_slope'),
                  ('CSF_TAU_closest_AV45','CSF_TAU_slope'),
                  ('CSF_ABETA_closest_AV45','CSF_ABETA_slope'),]
    to_scatter = [('CORTICAL_SUMMARY_prior', 'CORTICAL_SUMMARY_change')]
    for bl, change in bl_changes:
        to_add = [('CORTICAL_SUMMARY_prior', bl),
                  ('CORTICAL_SUMMARY_prior', change),
                  ('CORTICAL_SUMMARY_change', change),
                  (bl, change)]
        to_scatter += to_add

    # [(2.0, 7.0), (2.0, 26.0), (7.0, 8.0), (7.0, 10.0), (7.0, 26.0), (8.0, 10.0), (13.0, 22.0), (14.0, 22.0)]
    plotValueScatter(result_df, [13,22], xys=to_scatter, fit_reg=True, test=True)
    plotValueScatter(result_df, [14,22], xys=to_scatter, fit_reg=True, test=True)
    plotValueScatter(result_df, [2,7], xys=to_scatter, fit_reg=True, test=True)
    plotValueScatter(result_df, [2,26], xys=to_scatter, fit_reg=True, test=True)
    plotValueScatter(result_df, [7,8], xys=to_scatter, fit_reg=True, test=True)
    plotValueScatter(result_df, [7,10], xys=to_scatter, fit_reg=True, test=True)
    plotValueScatter(result_df, [7,26], xys=to_scatter, fit_reg=True, test=True)
    plotValueScatter(result_df, [8,10], xys=to_scatter, fit_reg=True, test=True)

