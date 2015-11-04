import pandas as pd
from sklearn.mixture import GMM, VBGMM, DPGMM
from sklearn.decomposition import PCA
from sklearn.cross_validation import LeaveOneOut

from sklearn.mixture import GMM, VBGMM, DPGMM
from sklearn.tree import DecisionTreeRegressor

# read in data
pattern_csv = '../pvcsummary_bl_patterns_and_change.csv'
df = pd.read_csv(pattern_csv)
df.set_index('rid',inplace=True)

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
result_keys = ['cortical_summary_bl', 
               'cortical_summary_change',
               'diag']
pattern_df = df[pattern_keys]
result_df = df[result_keys]

# convert to PCA space (99% variance)
pca_model = PCA(n_components=0.99, copy=True, whiten=True)
pattern_df_pca = pd.DataFrame(pca_model.fit_transform(pattern_df))
pattern_df_pca.set_index(pattern_df.index, inplace=True)
keys = pattern_df_pca.columns
pattern_df_pca = pattern_df_pca.merge(result_df, left_index=True, right_index=True)
target = 'cortical_summary_bl'

# Leave one out X-validation
loo = LeaveOneOut(n=len(pattern_df_pca.index))
points = []
for train_index, test_index in loo:
    train_rows = pattern_df_pca.iloc[train_index]
    test_rows = pattern_df_pca.iloc[test_index]
    X_train = train_rows[keys]
    Y_train = train_rows[target].tolist()
    X_test = test_rows[keys]
    Y_test = test_rows[target].tolist()
    regressor = DecisionTreeRegressor(random_state=0)
    regressor.fit(X_train,Y_train)
    score = regressor.predict(X_test)
    points.append((Y_test[0],score[0]))
    print "%s -> %s" % (score, Y_test)
plt.figure(1)
plt.scatter([_[0] for _ in points],[_[1] for _ in points])
plt.show()

# GMM
patterns_only = pattern_df_pca[keys]
g = DPGMM(n_components=10, 
          covariance_type='diag',
          alpha=5,
          tol=1e-10,
          random_state=None, 
          params='wmc', 
          init_params='wmc',
          verbose=True)
g.fit(patterns_only)
y_ = g.predict(patterns_only)
print y_
print "%s: %s" % (g.n_components, g.aic(patterns_only))
