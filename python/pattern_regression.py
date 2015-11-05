import pandas as pd
import numpy as np
import random
from collections import Counter
from scipy.stats import f_oneway
from __future__ import division

from sklearn.mixture import GMM, VBGMM, DPGMM
from sklearn.decomposition import PCA, KernelPCA
from sklearn.cross_validation import LeaveOneOut
from sklearn.mixture import GMM, VBGMM, DPGMM
from sklearn.tree import DecisionTreeRegressor
from sklearn.preprocessing import StandardScaler

# read in data
pattern_csv = '../pvcsummary_bl_patterns_and_change.csv'
bl_csv = '../pvcsummary_bl_vs_change.json'
raw_df = pd.read_json(bl_csv)
raw_df.set_index('rid', inplace=True)
df = pd.read_csv(pattern_csv)
df.set_index('rid',inplace=True)


def sample(data):
    sample = [random.choice(data) for _ in xrange(len(data))]
    return sample

def bootstrap_t_test(treatment, control, nboot = 1000, direction = "less"):
    treatment_len = len(treatment)
    control_len = len(control)
    Z = treatment+control
    tstat = abs(np.mean(treatment)-np.mean(control))
    tboot = np.zeros(nboot)
    for i in xrange(nboot):
        sboot = sample(Z)
        #sboot = pd.DataFrame(np.array(sboot), columns=['treat', 'vals'])
        #tboot[i] = np.mean(sboot['vals'][sboot['treat'] == 1]) - np.mean(sboot['vals'][sboot['treat'] == 0]) - tstat
        tboot[i] = np.mean(sboot[:treatment_len]) - np.mean(sboot[treatment_len:])
    pvalue = len([_ for _ in tboot if abs(_) >= tstat])/float(nboot)
    return pvalue


# extract patterns and summary value, summary change
pattern_keys = ['ctx-lh-caudalanteriorcingulate_bl',
                'ctx-lh-caudalmiddlefrontal_bl',
                'ctx-lh-frontalpole_bl',
                'ctx-lh-inferiorparietal_bl',
                'ctx-lh-isthmuscingulate_bl',
                'ctx-lh-lateralorbitofrontal_bl',
                'ctx-lh-medialorbitofrontal_bl',
                'ctx-lh-middletemporal_bl',
                'ctx-lh-parsopercularis_bl',
                'ctx-lh-parsorbitalis_bl',
                'ctx-lh-parstriangularis_bl',
                'ctx-lh-posteriorcingulate_bl',
                'ctx-lh-precuneus_bl',
                'ctx-lh-rostralanteriorcingulate_bl',
                'ctx-lh-rostralmiddlefrontal_bl',
                'ctx-lh-superiorfrontal_bl',
                'ctx-lh-superiorparietal_bl',
                'ctx-lh-superiortemporal_bl',
                'ctx-lh-supramarginal_bl',
                'ctx-rh-caudalanteriorcingulate_bl',
                'ctx-rh-caudalmiddlefrontal_bl',
                'ctx-rh-frontalpole_bl',
                'ctx-rh-inferiorparietal_bl',
                'ctx-rh-isthmuscingulate_bl',
                'ctx-rh-lateralorbitofrontal_bl',
                'ctx-rh-medialorbitofrontal_bl',
                'ctx-rh-middletemporal_bl',
                'ctx-rh-parsopercularis_bl',
                'ctx-rh-parsorbitalis_bl',
                'ctx-rh-parstriangularis_bl',
                'ctx-rh-posteriorcingulate_bl',
                'ctx-rh-precuneus_bl',
                'ctx-rh-rostralanteriorcingulate_bl',
                'ctx-rh-rostralmiddlefrontal_bl',
                'ctx-rh-superiorfrontal_bl',
                'ctx-rh-superiorparietal_bl',
                'ctx-rh-superiortemporal_bl',
                'ctx-rh-supramarginal_bl']
rchange_keys = ['ctx-lh-caudalanteriorcingulate_change',
                'ctx-lh-caudalmiddlefrontal_change',
                'ctx-lh-frontalpole_change',
                'ctx-lh-inferiorparietal_change',
                'ctx-lh-isthmuscingulate_change',
                'ctx-lh-lateralorbitofrontal_change',
                'ctx-lh-medialorbitofrontal_change',
                'ctx-lh-middletemporal_change',
                'ctx-lh-parsopercularis_change',
                'ctx-lh-parsorbitalis_change',
                'ctx-lh-parstriangularis_change',
                'ctx-lh-posteriorcingulate_change',
                'ctx-lh-precuneus_change',
                'ctx-lh-rostralanteriorcingulate_change',
                'ctx-lh-rostralmiddlefrontal_change',
                'ctx-lh-superiorfrontal_change',
                'ctx-lh-superiorparietal_change',
                'ctx-lh-superiortemporal_change',
                'ctx-lh-supramarginal_change',
                'ctx-rh-caudalanteriorcingulate_change',
                'ctx-rh-caudalmiddlefrontal_change',
                'ctx-rh-frontalpole_change',
                'ctx-rh-inferiorparietal_change',
                'ctx-rh-isthmuscingulate_change',
                'ctx-rh-lateralorbitofrontal_change',
                'ctx-rh-medialorbitofrontal_change',
                'ctx-rh-middletemporal_change',
                'ctx-rh-parsopercularis_change',
                'ctx-rh-parsorbitalis_change',
                'ctx-rh-parstriangularis_change',
                'ctx-rh-posteriorcingulate_change',
                'ctx-rh-precuneus_change',
                'ctx-rh-rostralanteriorcingulate_change',
                'ctx-rh-rostralmiddlefrontal_change',
                'ctx-rh-superiorfrontal_change',
                'ctx-rh-superiorparietal_change',
                'ctx-rh-superiortemporal_change',
                'ctx-rh-supramarginal_change']
result_keys = ['cortical_summary_bl', 
               'cortical_summary_change',
               'diag']

pattern_df = df[pattern_keys]
result_df = df[result_keys]
rchange_df = df[rchange_keys]
uptake_df = raw_df[pattern_keys]

# convert to PCA space (99% variance)
pca_model = PCA(n_components=len(pattern_df.columns), copy=True, whiten=True)
pattern_df_pca = pd.DataFrame(pca_model.fit_transform(pattern_df))
for _ in pca_model.explained_variance_ratio_:
    print _

pattern_df_pca.set_index(pattern_df.index, inplace=True)
pca_keys = pattern_df_pca.columns
pattern_df_pca = pattern_df_pca.merge(result_df, left_index=True, right_index=True)


# Leave one out X-validation
loo = LeaveOneOut(n=len(pattern_df_pca.index))
target = 'cortical_summary_bl'
points = []
for train_index, test_index in loo:
    train_rows = pattern_df_pca.iloc[train_index]
    test_rows = pattern_df_pca.iloc[test_index]
    X_train = train_rows[pca_keys]
    Y_train = train_rows[target].tolist()
    X_test = test_rows[pca_keys]
    Y_test = test_rows[target].tolist()
    regressor = DecisionTreeRegressor(random_state=0)
    regressor.fit(X_train,Y_train)
    score = regressor.predict(X_test)
    points.append((Y_test[0],score[0]))
    print "%s -> %s" % (score, Y_test)
plt.figure(1)
plt.scatter([_[0] for _ in points],[_[1] for _ in points])
plt.show()

# convert to Kernel PCA Projection
kpca_model = KernelPCA(n_components=None,
                       kernel='rbf', 
                       fit_inverse_transform=True, 
                       gamma=1,
                       alpha=1.0)
X_kpca = kpca_model.fit_transform(pattern_df)
X_back = kpca_model.inverse_transform(X_kpca)
pattern_df_kpca = pd.DataFrame(X_kpca)
kpca_keys = pattern_df_kpca.columns
# scale
scaler = StandardScaler().fit(pattern_df_kpca)
pattern_df_kpca = pd.DataFrame(scaler.transform(pattern_df_kpca))
pattern_df_kpca['rid'] = result_df.index
pattern_df_kpca.set_index('rid', inplace=True)
pattern_df_kpca = pattern_df_kpca.merge(result_df, left_index=True, right_index=True)

# GMM
patterns_only = pattern_df_kpca[kpca_keys]

# model selection for alpha by looking at lower bound on marginal likelihood, or log prob of X under model
alphas = [0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1e3, 5e3, 1e4]
alpha_to_clusters = {}
components = 20
for a in alphas:
    clusters = []
    scores = []
    for i in range(5):
        g = DPGMM(n_components=components, 
                  covariance_type='diag',
                  alpha=a,
                  tol=1e-3,
                  n_iter=700,
                  params='wmc', 
                  init_params='wmc',
                  verbose=False)
        g.fit(patterns_only)
        print g.converged_
        y_ = g.predict(patterns_only)
        score = np.mean(g.score(patterns_only))
        scores.append(score)
        clusters.append(len(set(y_)))
    print scores
    print clusters
    alpha_to_clusters[a] = np.mean(scores)

print sorted(alpha_to_clusters.items(), key=lambda x:x[1], reverse=True)

# Cluster
alpha = 100
components = 20
models_scores = {}
for i in range(10):
    g = DPGMM(n_components=components, 
              covariance_type='diag',
              alpha=alpha,
              tol=1e-3,
              n_iter=700,
              params='wmc', 
              init_params='wmc',
              verbose=False)
    g.fit(patterns_only)
    print g.converged_
    y_ = g.predict(patterns_only)
    print Counter(y_)
    membership = np.zeros((len(patterns_only.index),components))
    for i,member in enumerate(y_):
        membership[i][member] = 1
    bound = g.lower_bound(patterns_only, membership)
    score = np.mean(g.score(patterns_only))
    models_scores[score] = g

best_score = max(models_scores.keys())
best_model = models_scores[best_score]
y_ = best_model.predict(patterns_only)
probs = best_model.predict_proba(patterns_only)
to_remove = []
for i, (y, yprob) in enumerate(zip(y_, probs)):
    print yprob
    print yprob[y]
    if yprob[y] < 0.99:
        to_remove.append(i)

for c in set(y_):
    print c
    print '\t%s' % len([_ for _ in y_ if _ == c])

groups = list(set(y_))
result_df['membership'] = y_


# anova for group difference in regional change
rchange_members = rchange_df.merge(result_df, left_index=True, right_index=True)
for i,rchange_key in enumerate(rchange_keys):
    groupvals = {}
    for g in groups:
        members = rchange_members[rchange_members['membership']==g]
        groupvals[g] = members[rchange_key].tolist()
    Fvalue,pvalue = f_oneway(*groupvals.values())
    print "%s: %s, %s" % (rchange_key, Fvalue, pvalue)


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
uptake_members = uptake_df.merge(result_df, left_index=True, right_index=True)
for i,pattern_key in enumerate(pattern_keys):
    long = pd.melt(uptake_members, id_vars='membership', value_vars=pattern_key)
    plt.figure(i+1)
    for g in groups:
        members = long[long['membership']==g]
        if len(members) < 3:
            continue
        members['value'].plot(kind='density', label='%s_%s' % (rchange_key, g))
    plt.legend()


# plot against summary uptake/change values
summary_bl_long = pd.melt(result_df, id_vars='membership',value_vars='cortical_summary_bl')
summary_change_long = pd.melt(result_df, id_vars='membership',value_vars='cortical_summary_change')
plt.figure()
for g in groups:
    members = summary_bl_long[summary_bl_long['membership']==g]
    if len(members) < 3:
        print "Cluster %s has only %s members" % (g, len(members))
        continue
    members['value'].plot(kind='density', label=g, alpha=0.5)

plt.legend()
plt.show()



print "%s: %s" % (g.n_components, g.aic(patterns_only))
