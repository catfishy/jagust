import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

BRAAK_KEYS = ['BRAAK1','BRAAK2','BRAAK3','BRAAK4','BRAAK5']
av45_groupings = {'allregions': '../datasets/pvc_adni_av45/allregions_output.csv',
                  'mostregions': '../datasets/pvc_adni_av45/mostregions_output.csv',
                  'aggregions': '../datasets/pvc_adni_av45/aggregions_output.csv',
                  'tauregions': '../datasets/pvc_adni_av45/tauregions_output.csv'}
av1451_groupings = {'allregions': '../datasets/pvc_adni_av1451/allregions_output.csv',
                  'mostregions': '../datasets/pvc_adni_av1451/mostregions_output.csv',
                  'tauregions': '../datasets/pvc_adni_av1451/tauregions_output.csv'}



rmse_residuals = pd.DataFrame()
braak_stages = pd.DataFrame()
for grouping_name, output_file in av1451_groupings.iteritems():
    df = pd.read_csv(output_file)
    for sid, rows in df.groupby('subject'):
        rmse_rows = rows[['rmse_residual','timepoint','subject']].drop_duplicates().dropna()
        rmse_rows['grouping'] = grouping_name
        rmse_residuals = pd.concat((rmse_residuals,rmse_rows))


for grouping_name in groupings.keys():
    sns.distplot(rmse_residuals[rmse_residuals.grouping==grouping_name].rmse_residual,
                 hist=True,kde=False,bins=100, label=grouping_name)

plt.legend()
plt.show()

ax = sns.violinplot(x='grouping', y='rmse_residual', data=rmse_residuals, scale='width')
plt.show()


