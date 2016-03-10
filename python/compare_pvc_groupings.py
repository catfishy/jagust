import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

all_df = pd.read_csv('../datasets/pvc_adni_av1451/allregions_output.csv')
tau_df = pd.read_csv('../datasets/pvc_adni_av1451/tauregions_output.csv')

BRAAK_KEYS = ['BRAAK1','BRAAK2','BRAAK3','BRAAK4','BRAAK5']

rmse_residuals = pd.DataFrame()
braak_stages = pd.DataFrame()

for sid, rows in all_df.groupby('subject'):
    braak_rows = rows[rows.name.isin(BRAAK_KEYS)][['subject','timepoint','pvcval','nonpvcval','name']]
    braak_rows['grouping'] = 'allregions'
    braak_stages = pd.concat((braak_stages,braak_rows))
    rmse_rows = rows[['rmse_residual','timepoint','subject']].drop_duplicates().dropna()
    rmse_rows['grouping'] = 'allregions'
    rmse_residuals = pd.concat((rmse_residuals,rmse_rows))

for sid, rows in tau_df.groupby('subject'):
    braak_rows = rows[rows.name.isin(BRAAK_KEYS)][['subject','timepoint','pvcval','nonpvcval','name']]
    braak_rows['grouping'] = 'tauregions'
    braak_stages = pd.concat((braak_stages,braak_rows))
    rmse_rows = rows[['rmse_residual','timepoint','subject']].drop_duplicates().dropna()
    rmse_rows['grouping'] = 'tauregions'
    rmse_residuals = pd.concat((rmse_residuals,rmse_rows))


tau_rmse = rmse_residuals[rmse_residuals.grouping=='tauregions'].set_index('subject')[['rmse_residual']]
tau_rmse.columns = ['tau_rmse']
all_rmse = rmse_residuals[rmse_residuals.grouping=='allregions'].set_index('subject')[['rmse_residual']]
all_rmse.columns = ['all_rmse']
rmse = all_rmse.merge(tau_rmse, left_index=True, right_index=True)

tau_braak = braak_stages[braak_stages.grouping=='tauregions'].set_index(['subject','name'])[['nonpvcval','pvcval']]
tau_braak.columns = ['tau_nonpvcval','tau_pvcval']
all_braak = braak_stages[braak_stages.grouping=='allregions'].set_index(['subject','name'])[['nonpvcval','pvcval']]
all_braak.columns = ['all_nonpvcval','all_pvcval']
braak = all_braak.merge(tau_braak, left_index=True, right_index=True).dropna().drop_duplicates()

braak.plot(kind='scatter',x='all_pvcval',y='tau_pvcval')
plt.show()

# look at residuals
residuals = {}
for name, rows in all_df.groupby('name'):
    res = list(rows['residual'])
    sizes = list(rows['groupsize'])
    mean = np.mean(res)
    sizes = np.mean(sizes)
    if np.isnan(mean):
        continue
    residuals[name] = (mean,sizes)

ranked = sorted(residuals.items(), key=lambda x: x[1][0], reverse=True)
for r in ranked[:80]:
    print r
