import pandas as pd
import numpy as np
import random
from collections import Counter, defaultdict
from scipy.stats import f_oneway
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

# make up some randomly distributed data
seed(1234)
npts = 1000
x = uniform(-2,2,npts)
y = uniform(-2,2,npts)
z = gauss(x,y,Sigma=np.asarray([[1.,.5],[0.5,1.]]),mu=np.asarray([0.,0.]))
plot_countour(x,y,z)


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
        covariance_type='diag',
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
Z1 = mlab.normpdf(x,g_1.means_[0][0],np.sqrt(g_1.covars_[0][0]))
Z2 = mlab.normpdf(x,g_1.means_[1][0],np.sqrt(g_1.covars_[1][0]))
Z = (Z2-Z1)
diffpts = zip(x,Z)
diffpts = [(a,b) for a,b in diffpts if a > 0.6 and a < 1.0]
zeropt = sorted(diffpts, key=lambda x: abs(x[1]))[0][0]
summary_bigref.plot(kind='density')
plt.plot(x,Z1,label='gaussian one')
plt.plot(x,Z2,label='gaussian two')
plt.plot(x,Z,label='gaussian diff')
plt.axvline(x=zeropt, label='Threshold @ %s' % zeropt)
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
probs.loc[probs.prob_0>=0.90, 'color'] = 'r'
probs.loc[probs.prob_1>=0.90, 'color'] = 'b'

# plot 2D GMM
delta = 0.01
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

summary_vs_wm.plot(kind='scatter',x='wm_wcereb',y='summary_wcereb',c=probs.color)
CS1 = plt.contour(X, Y, Z1)
CS2 = plt.contour(X, Y, Z2)
CS = plt.contour(X, Y, Z, [0.0])
zc = CS.collections[0]
plt.setp(zc, linewidth=4)

diffs = pd.DataFrame(abs(probs['prob_0']-probs['prob_1']))
diffs.plot(kind='hist')

plt.show()

