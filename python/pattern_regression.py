import pandas as pd
import numpy as np
import random
from collections import Counter, defaultdict
from scipy.stats import f_oneway
from __future__ import division
import matplotlib.pyplot as plt
from matplotlib import cm
import itertools

from sklearn.mixture import GMM, VBGMM, DPGMM
from sklearn.decomposition import PCA, KernelPCA
from sklearn.cross_validation import LeaveOneOut
from sklearn.mixture import GMM, VBGMM, DPGMM
from sklearn.tree import DecisionTreeRegressor
from sklearn.preprocessing import StandardScaler

pd.options.display.mpl_style = 'default'

def sample(data):
    return data[np.random.randint(0,len(data),(1,len(data)))[0]]

def bootstrap_t_test(treatment, control, nboot = 10000):
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
    for i in xrange(nboot):
        sboot = sample(Z)
        tboot[i] = np.mean(sboot[:treatment_len]) - np.mean(sboot[treatment_len:])
    pvalue = np.sum(abs(tboot)>=tstat) / float(nboot)
    return pvalue

def bootstrap_mean(vals, nboot=10000):
    vals = np.array(vals)
    tboot = np.zeros(nboot)
    for i in xrange(nboot):
        sboot = sample(vals)
        tboot[i] = np.mean(sboot)
    return np.mean(tboot)


def gmm_sweep_alpha(components, data, covar_type='full'):
    # model selection for alpha by looking at lower bound on marginal likelihood, or log prob of X under model
    alphas = [0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
    alpha_to_clusters = {}
    for a in alphas:
        clusters = []
        scores = []
        print a
        for i in range(5):
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
            score = np.mean(g.score(data))
            scores.append(score)
            clusters.append(len(set(y_)))
        print scores
        print clusters
        alpha_to_clusters[a] = np.mean(scores)
    return alpha_to_clusters

def group_comparisons(df, groups, keys):
    data=defaultdict(dict)
    for group1, group2 in itertools.combinations(groups,2):
        group1_members = df[df['membership']==group1]
        group2_members = df[df['membership']==group2]
        print "GRP %s (%s) vs GRP %s (%s)" % (group1, len(group1_members), group2, len(group2_members))
        for i, key in enumerate(keys):
            group1_vals = group1_members[key].tolist()
            group2_vals = group2_members[key].tolist()
            group1_boot = bootstrap_mean(group1_vals)
            group2_boot = bootstrap_mean(group2_vals)
            pvalue = bootstrap_t_test(group1_vals, group2_vals, nboot = 10000)
            data[(group1,group2)][key] = (pvalue, group1_boot, group2_boot)
            # if pvalue <= 0.05:
            #     print "\t%s: %s" % (key, pvalue)
            print "\t%s: %s (%s, %s)" % (key, pvalue, group1_boot, group2_boot)
    return data

# extract patterns and summary value, summary change
pattern_keys = ['BRAIN_STEM', 'CTX_LH_BANKSSTS', 'CTX_LH_CAUDALANTERIORCINGULATE', 'CTX_LH_CAUDALMIDDLEFRONTAL', 'CTX_LH_CUNEUS', 'CTX_LH_ENTORHINAL', 'CTX_LH_FRONTALPOLE', 'CTX_LH_FUSIFORM', 'CTX_LH_INFERIORPARIETAL', 'CTX_LH_INFERIORTEMPORAL', 'CTX_LH_INSULA', 'CTX_LH_ISTHMUSCINGULATE', 'CTX_LH_LATERALOCCIPITAL', 'CTX_LH_LATERALORBITOFRONTAL', 'CTX_LH_LINGUAL', 'CTX_LH_MEDIALORBITOFRONTAL', 'CTX_LH_MIDDLETEMPORAL', 'CTX_LH_PARACENTRAL', 'CTX_LH_PARAHIPPOCAMPAL', 'CTX_LH_PARSOPERCULARIS', 'CTX_LH_PARSORBITALIS', 'CTX_LH_PARSTRIANGULARIS', 'CTX_LH_PERICALCARINE', 'CTX_LH_POSTCENTRAL', 'CTX_LH_POSTERIORCINGULATE', 'CTX_LH_PRECENTRAL', 'CTX_LH_PRECUNEUS', 'CTX_LH_ROSTRALANTERIORCINGULATE', 'CTX_LH_ROSTRALMIDDLEFRONTAL', 'CTX_LH_SUPERIORFRONTAL', 'CTX_LH_SUPERIORPARIETAL', 'CTX_LH_SUPERIORTEMPORAL', 'CTX_LH_SUPRAMARGINAL', 'CTX_LH_TEMPORALPOLE', 'CTX_LH_TRANSVERSETEMPORAL', 'CTX_LH_UNKNOWN', 'CTX_RH_BANKSSTS', 'CTX_RH_CAUDALANTERIORCINGULATE', 'CTX_RH_CAUDALMIDDLEFRONTAL', 'CTX_RH_CUNEUS', 'CTX_RH_ENTORHINAL', 'CTX_RH_FRONTALPOLE', 'CTX_RH_FUSIFORM', 'CTX_RH_INFERIORPARIETAL', 'CTX_RH_INFERIORTEMPORAL', 'CTX_RH_INSULA', 'CTX_RH_ISTHMUSCINGULATE', 'CTX_RH_LATERALOCCIPITAL', 'CTX_RH_LATERALORBITOFRONTAL', 'CTX_RH_LINGUAL', 'CTX_RH_MEDIALORBITOFRONTAL', 'CTX_RH_MIDDLETEMPORAL', 'CTX_RH_PARACENTRAL', 'CTX_RH_PARAHIPPOCAMPAL', 'CTX_RH_PARSOPERCULARIS', 'CTX_RH_PARSORBITALIS', 'CTX_RH_PARSTRIANGULARIS', 'CTX_RH_PERICALCARINE', 'CTX_RH_POSTCENTRAL', 'CTX_RH_POSTERIORCINGULATE', 'CTX_RH_PRECENTRAL', 'CTX_RH_PRECUNEUS', 'CTX_RH_ROSTRALANTERIORCINGULATE', 'CTX_RH_ROSTRALMIDDLEFRONTAL', 'CTX_RH_SUPERIORFRONTAL', 'CTX_RH_SUPERIORPARIETAL', 'CTX_RH_SUPERIORTEMPORAL', 'CTX_RH_SUPRAMARGINAL', 'CTX_RH_TEMPORALPOLE', 'CTX_RH_TRANSVERSETEMPORAL', 'CTX_RH_UNKNOWN', 'LEFT_ACCUMBENS_AREA', 'LEFT_AMYGDALA', 'LEFT_CAUDATE', 'LEFT_CEREBELLUM_CORTEX', 'LEFT_CEREBELLUM_WHITE_MATTER', 'LEFT_CEREBRAL_WHITE_MATTER', 'LEFT_HIPPOCAMPUS', 'LEFT_PALLIDUM', 'LEFT_PUTAMEN', 'LEFT_THALAMUS_PROPER', 'LEFT_VENTRALDC', 'RIGHT_ACCUMBENS_AREA', 'RIGHT_AMYGDALA', 'RIGHT_CAUDATE', 'RIGHT_CEREBELLUM_CORTEX', 'RIGHT_CEREBELLUM_WHITE_MATTER', 'RIGHT_CEREBRAL_WHITE_MATTER', 'RIGHT_HIPPOCAMPUS', 'RIGHT_PALLIDUM', 'RIGHT_PUTAMEN', 'RIGHT_THALAMUS_PROPER', 'RIGHT_VENTRALDC']
prior_keys = ['%s_prior' % _ for _ in pattern_keys]
post_keys = ['%s_post' % _ for _ in pattern_keys]
rchange_keys = ['%s_change' % _ for _ in pattern_keys]
result_keys = ['CORTICAL_SUMMARY_post', 'CORTICAL_SUMMARY_prior', 'CORTICAL_SUMMARY_change', 'diag_prior', 'diag_post']

# read in data
ref_key = 'COMPOSITE_REF'
patterns_csv = '../datasets/pvc_allregions_uptake_change.csv'
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

pattern_post_df = raw_df[post_keys].copy()
pattern_post_df.columns = [_.replace('_post','') for _ in pattern_post_df.columns]
pattern_post_df = pattern_post_df.divide(pattern_post_df.sum(axis=1),axis=0)
pattern_post_df['timepoint'] = 'post'
pattern_post_df.set_index(raw_df.index, inplace=True)
pattern_post_df.set_index('timepoint', append=True, inplace=True)


# concat prior and post patterns
pattern_df = pd.concat((pattern_prior_df,pattern_post_df))
uptake_prior_df = raw_df[prior_keys]
uptake_prior_df.columns = [_.replace('_prior','') for _ in uptake_prior_df.columns]
uptake_post_df = raw_df[post_keys]
uptake_post_df.columns = [_.replace('_post','') for _ in uptake_post_df.columns]
result_df = raw_df[result_keys]
rchange_df = raw_df[rchange_keys]

# Scale in original space
scaler = StandardScaler().fit(pattern_df)
pattern_df_raw = pd.DataFrame(scaler.transform(pattern_df))
pattern_df_raw.set_index(pattern_df.index, inplace=True)
raw_keys = pattern_df_raw.columns

# convert to PCA space + whiten(scale)
pca_model = PCA(n_components=len(pattern_df.columns), copy=True, whiten=True)
pattern_df_pca = pd.DataFrame(pca_model.fit_transform(pattern_df))
pattern_df_pca.set_index(pattern_df.index, inplace=True)
pca_keys = pattern_df_pca.columns

# convert to Kernel PCA Projection and scale
kpca_model = KernelPCA(n_components=None,kernel='rbf',fit_inverse_transform=False,gamma=1,alpha=1.0)
X_kpca = kpca_model.fit_transform(pattern_df)
pattern_df_kpca = pd.DataFrame(X_kpca)
scaler = StandardScaler().fit(pattern_df_kpca)
pattern_df_kpca = pd.DataFrame(scaler.transform(pattern_df_kpca))
pattern_df_kpca.set_index(pattern_df.index, inplace=True)
kpca_keys = pattern_df_kpca.columns

# GMM
raw_patterns_only = pattern_df_raw[raw_keys]
raw_alpha_to_clusters = gmm_sweep_alpha(30, raw_patterns_only, covar_type='diag')
print sorted(raw_alpha_to_clusters.items(), key=lambda x:x[1], reverse=True)

pca_patterns_only = pattern_df_pca[pca_keys]
pca_alpha_to_clusters = gmm_sweep_alpha(20, pca_patterns_only, covar_type='diag')
print sorted(pca_alpha_to_clusters.items(), key=lambda x:x[1], reverse=True)

kpca_patterns_only = pattern_df_kpca[kpca_keys]
kpca_alpha_to_clusters = gmm_sweep_alpha(50, kpca_patterns_only, covar_type='diag')
print sorted(kpca_alpha_to_clusters.items(), key=lambda x:x[1], reverse=True)

# Cluster
alpha = 30
components = 30
models_scores = {}
patterns_only = raw_patterns_only
for i in range(20):
    g = DPGMM(n_components=components, 
              covariance_type='diag',
              alpha=alpha,
              tol=1e-6,
              n_iter=700,
              params='wmc', 
              init_params='wmc',
              verbose=False)
    print "Fitting: %s" % i
    g.fit(patterns_only)
    print g.converged_
    y_ = g.predict(patterns_only)
    print Counter(y_)
    membership = np.zeros((len(patterns_only.index),components))
    for i,member in enumerate(y_):
        membership[i][member] = 1
    bound = g.lower_bound(patterns_only, membership)
    score = np.mean(g.score(patterns_only))
    models_scores[bound] = g

best_score = max(models_scores.keys())
best_model = models_scores[best_score]
y_ = best_model.predict(patterns_only)
probs = best_model.predict_proba(patterns_only)

# add membership to result df and determine significant clusters
y_df = pd.DataFrame(y_)
y_df.columns = ['group']
y_df.set_index(patterns_only.index, inplace=True)
for rid in result_df.index:
    result_df.loc[rid,'membership_prior'] = y_df.loc[(rid,'prior'),'group']
    result_df.loc[rid,'membership_post'] = y_df.loc[(rid,'post'),'group']

groups = list(set(y_))
big_groups = []
for g in groups:
    prior_members = set(result_df[result_df.membership_prior==g].index)
    post_members = set(result_df[result_df.membership_post==g].index)
    allmembers = len(prior_members | post_members)
    print "%s: %s" % (g,allmembers)
    if allmembers >= 10:
        big_groups.append(g)

# print group diagnoses
for g in big_groups:
    members_prior = result_df[result_df['membership_prior']==g]['diag_prior'].tolist()
    members_post = result_df[result_df['membership_post']==g]['diag_post'].tolist()
    diags = dict(Counter(members_prior+members_post))
    print "Group %s" % g
    print "\tN: %s" % diags.get('N',0)
    #print "\tSMC: %s" % diags.get('SMC',0)
    print "\tEMCI: %s" % diags.get('EMCI',0)
    print "\tLCMI: %s" % diags.get('LMCI',0)
    print "\tAD: %s" % diags.get('AD',0)


# between group bootstrap hypothesis test for summary/regional change
change_members = rchange_df.merge(result_df, left_index=True, right_index=True)
change_members['membership'] = change_members['membership_prior']
change_keys = rchange_keys
#change_keys = ['cortical_summary_change']
change_pvalues = group_comparisons(change_members, big_groups, change_keys)

# between group bootstrap hypothesis test for summary/regional uptake
# merge uptakes
uptake_members_prior = uptake_prior_df.merge(result_df[['membership_prior','CORTICAL_SUMMARY_prior']], left_index=True, right_index=True)
uptake_members_prior['membership'] = uptake_members_prior['membership_prior']
uptake_members_prior['CORTICAL_SUMMARY'] = uptake_members_prior['CORTICAL_SUMMARY_prior']
uptake_members_post = uptake_post_df.merge(result_df[['membership_post','CORTICAL_SUMMARY_post']], left_index=True, right_index=True)
uptake_members_post['membership'] = uptake_members_post['membership_post']
uptake_members_post['CORTICAL_SUMMARY'] = uptake_members_post['CORTICAL_SUMMARY_post']
uptake_members = pd.concat((uptake_members_prior,uptake_members_post))
#uptake_keys = ['CORTICAL_SUMMARY']
uptake_keys = pattern_keys
uptake_pvalues = group_comparisons(uptake_members, big_groups, pattern_keys)

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
