import pandas as pd
import matplotlib.pyplot as plt
from pandas.stats.api import ols
import numpy as np
#pd.options.display.mpl_style = 'default'

old_file = '06_25_15.csv'
new_file = '11_25_15.csv'

key_columns = ['RID', 'EXAMDATE']
columns = ['SUMMARYSUVR_WHOLECEREBNORM','CORTICAL_SUMMARY','BRAIN_STEM', 'WHOLECEREBELLUM',
		   'ERODED_SUBCORTICALWM', 'SUMMARYSUVR_COMPOSITE_REFNORM', 'SUMMARYSUVR_ERODEDWMNORM']

old_df = pd.read_csv(old_file)
new_df = pd.read_csv(new_file)

# add new fields
old_df['CORTICAL_SUMMARY'] = old_df['SUMMARYSUVR_WHOLECEREBNORM'] * old_df['WHOLECEREBELLUM']
old_df['SUMMARYSUVR_ERODEDWMNORM'] = old_df['CORTICAL_SUMMARY'] / old_df['ERODED_SUBCORTICALWM']
new_df['CORTICAL_SUMMARY'] = new_df['SUMMARYSUVR_WHOLECEREBNORM'] * new_df['WHOLECEREBELLUM']
new_df['SUMMARYSUVR_ERODEDWMNORM'] = new_df['CORTICAL_SUMMARY'] / new_df['ERODED_SUBCORTICALWM']

old_df = old_df[key_columns+columns]
old_df.columns = key_columns + ['%s_old' % _ for _ in columns]
new_df = new_df[key_columns+columns]
new_df.columns = key_columns + ['%s_new' % _ for _ in columns]
merged_df = pd.merge(old_df, new_df, on=['RID','EXAMDATE'], how='inner')
bl_df = pd.DataFrame()
v2_df = pd.DataFrame()
v3_df = pd.DataFrame()
for rid, rows in merged_df.groupby('RID'):
	sorteddf = rows.set_index('EXAMDATE').sort().reset_index()
	for i, row in sorteddf.iterrows():
		row_df = pd.DataFrame(row).T
		if i == 0:
			bl_df = pd.concat((bl_df,row_df))
		elif i == 1:
			v2_df = pd.concat((v2_df,row_df))
		elif i == 2:
			v3_df = pd.concat((v3_df,row_df))

bl_df.reset_index(drop=True,inplace=True)
v2_df.reset_index(drop=True,inplace=True)
v3_df.reset_index(drop=True,inplace=True)

props={'facecolor':'wheat','alpha':0.5,'pad':10}
for c in columns:
	fig, axes = plt.subplots(nrows=1, ncols=3)
	oldkey = '%s_old' % c
	newkey = '%s_new' % c
	# bl
	model = ols(x=bl_df[oldkey],y=bl_df[newkey])
	x = model.beta.x
	intercept = model.beta.intercept
	r2 = '{0:.3f}'.format(model.r2)
	equation = 'Y ~ %sx + %s' % ('{0:.3f}'.format(x),'{0:.3f}'.format(intercept))
	ax = axes[0]
	bl_df.plot(ax=ax, kind='scatter',x='%s_old' % c, y = '%s_new' % c)
	xvals = np.array(sorted(list(bl_df[oldkey])))
	ax.plot(xvals, xvals*x + intercept, 'r')
	ax.set_title('BL')
	ax.annotate('%s\nR2=%s' % (equation,r2), xy=(0.05,0.95), xycoords='axes fraction', bbox=props, verticalalignment='top', fontsize=15, color='k')
	# v2
	model = ols(x=v2_df[oldkey],y=v2_df[newkey])
	x = model.beta.x
	intercept = model.beta.intercept
	r2 = '{0:.3f}'.format(model.r2)
	equation = 'Y ~ %sx + %s' % ('{0:.3f}'.format(x),'{0:.3f}'.format(intercept))
	ax = axes[1]
	v2_df.plot(ax=ax, kind='scatter',x='%s_old' % c, y = '%s_new' % c)
	xvals = np.array(sorted(list(v2_df[oldkey])))
	ax.plot(xvals, xvals*x + intercept, 'r')
	ax.set_title('V2')
	ax.annotate('%s\nR2=%s' % (equation,r2), xy=(0.05,0.95), xycoords='axes fraction', bbox=props, verticalalignment='top', fontsize=15, color='k')
	# v3
	model = ols(x=v3_df[oldkey],y=v3_df[newkey])
	x = model.beta.x
	intercept = model.beta.intercept
	r2 = '{0:.3f}'.format(model.r2)
	equation = 'Y ~ %sx + %s' % ('{0:.3f}'.format(x),'{0:.3f}'.format(intercept))
	ax = axes[2]
	v3_df.plot(ax=ax, kind='scatter',x='%s_old' % c, y = '%s_new' % c)
	xvals = np.array(sorted(list(v3_df[oldkey])))
	ax.plot(xvals, xvals*x + intercept, 'r')
	ax.set_title('V3')
	ax.annotate('%s\nR2=%s' % (equation,r2), xy=(0.05,0.95), xycoords='axes fraction', bbox=props, verticalalignment='top', fontsize=15, color='k')

plt.show()