import pandas as pd
import numpy as np
import random
from collections import Counter, defaultdict
from scipy.stats import f_oneway, norm
from scipy.misc import logsumexp
from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import itertools

from sklearn.mixture import GMM, VBGMM, DPGMM, log_multivariate_normal_density
from sklearn.decomposition import PCA, KernelPCA
from sklearn.cross_validation import LeaveOneOut
from sklearn.mixture import GMM, VBGMM, DPGMM
from sklearn.tree import DecisionTreeRegressor
from sklearn.preprocessing import StandardScaler

pd.options.display.mpl_style = 'default'

# read in data
uptake_csv = '../datasets/nontp_allregions_uptake.csv'
df = pd.read_csv(uptake_csv)
df.set_index('rid', inplace=True)
df.sort(inplace=True)

summary_bigref = pd.DataFrame(df['cortical_summary']/df['composite_ref'])
summary_bigref.set_index(df.index, inplace=True)
summary_bigref.columns = ['summary_bigref']
summary_wcereb = pd.DataFrame(df['cortical_summary']/df['whole_cerebellum'])
summary_wcereb.set_index(df.index,inplace=True)
summary_wcereb.columns = ['summary_wcereb']
wm_wcereb = pd.DataFrame(df['white_matter']/df['whole_cerebellum'])
wm_wcereb.set_index(df.index,inplace=True)
wm_wcereb.columns = ['wm_wcereb']

wm_vs_summary = wm_wcereb.merge(summary_wcereb, left_index=True, right_index=True)

# 1D GMM
g_1 = GMM(n_components=2, 
        covariance_type='full',
        tol=1e-6,
        n_iter=700,
        params='wmc', 
        init_params='wmc')
g_1.fit(summary_bigref)
# calculate prob, disregard weights
lpr = log_multivariate_normal_density(summary_bigref,g_1.means_,g_1.covars_,g_1.covariance_type)
logprob = logsumexp(lpr,axis=1)
responsibilities = np.exp(lpr - logprob[:, np.newaxis])
probs = pd.DataFrame(responsibilities)
probs.set_index(summary_bigref.index,inplace=True)
probs.columns = ['prob_0','prob_1']
probs.loc[:,'color'] = 'k'
probs.loc[probs.prob_0>=0.90, 'color'] = 'r'
probs.loc[probs.prob_1>=0.90, 'color'] = 'b'
# plot 1D GMM
delta= 0.0001
x = np.arange(0.5, 1.2, delta)
mu_1, sigma_1 = (g_1.means_[0][0],np.sqrt(g_1.covars_[0][0]))
mu_2, sigma_2 = (g_1.means_[1][0],np.sqrt(g_1.covars_[1][0]))
intervals_1 = norm.interval(0.99,loc=mu_1,scale=sigma_1)
intervals_2 = norm.interval(0.99,loc=mu_2,scale=sigma_2)
interval_1_x = np.arange(intervals_1[0][0],intervals_1[1][0],delta)
interval_2_x = np.arange(intervals_2[0][0],intervals_2[1][0],delta)
print intervals_1
print intervals_2
Z1 = mlab.normpdf(x,mu_1,sigma_1)
Z2 = mlab.normpdf(x,mu_2,sigma_2)
Z = (Z2-Z1)
diffpts = zip(x,Z)
diffpts = [(a,b) for a,b in diffpts if a > 0.6 and a < 1.0]
zeropt = sorted(diffpts, key=lambda x: abs(x[1]))[0][0]
min_interval = min(intervals_1[1][0],intervals_2[1][0])
max_interval = max(intervals_1[0][0],intervals_2[0][0])
summary_bigref.plot(kind='density')
plt.plot(x,Z1,label='gaussian one')
plt.plot(x,Z2,label='gaussian two')
plt.plot(x,Z,label='gaussian diff')
plt.axvline(x=zeropt, label='Threshold @ %s' % zeropt)
plt.axvline(x=min_interval, label='Threshold @ %s' % min_interval)
plt.axvline(x=max_interval, label='Threshold @ %s' % max_interval)
plt.fill_between(interval_1_x,mlab.normpdf(interval_1_x,mu_1,sigma_1),alpha=0.5)
plt.fill_between(interval_2_x,mlab.normpdf(interval_2_x,mu_2,sigma_2),alpha=0.5)
plt.legend()
plt.show()

# 2D GMM
g_2 = GMM(n_components=2, 
          covariance_type='full',
          tol=1e-6,
          n_iter=700,
          params='wmc', 
          init_params='wmc')
g_2.fit(wm_vs_summary)

# calculate prob with equal weights
lpr = log_multivariate_normal_density(wm_vs_summary,g_2.means_,g_2.covars_,g_2.covariance_type)
logprob = logsumexp(lpr,axis=1)
responsibilities = np.exp(lpr - logprob[:, np.newaxis])
probs = pd.DataFrame(responsibilities)
probs.set_index(wm_vs_summary.index,inplace=True)
probs.columns = ['prob_0','prob_1']
probs.loc[:,'color'] = 'k'
probs.loc[probs.prob_0>0.55, 'color'] = 'r'
probs.loc[probs.prob_1>0.55, 'color'] = 'b'

# plot 2D GMM
delta = 0.001
x = np.arange(0.0, 4.0, delta)
y = np.arange(0.0, 3.0, delta)
X, Y = np.meshgrid(x, y)
Z1 = mlab.bivariate_normal(X, Y, 
                           sigmax=np.sqrt(g_2.covars_[0][0][0]), 
                           sigmay=np.sqrt(g_2.covars_[0][1][1]),
                           mux=g_2.means_[0][0],
                           muy=g_2.means_[0][1],
                           sigmaxy=g_2.covars_[0][0][1])
Z2 = mlab.bivariate_normal(X, Y, 
                           sigmax=np.sqrt(g_2.covars_[1][0][0]), 
                           sigmay=np.sqrt(g_2.covars_[1][1][1]),
                           mux=g_2.means_[1][0],
                           muy=g_2.means_[1][1],
                           sigmaxy=g_2.covars_[1][0][1])
# difference of Gaussians
Z = (Z2 - Z1)

wm_vs_summary.plot(kind='scatter',x='wm_wcereb',y='summary_wcereb',c=probs.color)
CS1 = plt.contour(X, Y, Z1, [0.99])
CS2 = plt.contour(X, Y, Z2, [0.99])
#CS = plt.contour(X, Y, Z, [0.0])
#zc = CS.collections[0]
#plt.setp(zc, linewidth=4)

# diffs = pd.DataFrame(abs(probs['prob_0']-probs['prob_1']))
# diffs.plot(kind='hist')

plt.show()

