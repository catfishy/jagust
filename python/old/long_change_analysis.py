'''
Look at the components of SUVR change independently (prior to SUVR calculation)
'''
import pandas as pd
import seaborn as sns


from utils import slope, parseDate

INPUT_FILE = '../output/02_19_16/UCBERKELEYAV45_02_19_16_regular_tp.csv'


df = pd.read_csv(INPUT_FILE)
df = df[['RID','EXAMDATE','WHOLECEREBELLUM','FRONTAL','CINGULATE','PARIETAL','TEMPORAL','COMPOSITE','SUMMARYSUVR_WHOLECEREBNORM','BRAIN_STEM']]
df.loc[:,'EXAMDATE'] = df.loc[:,'EXAMDATE'].apply(parseDate)
slopes = []
for rid, rows in df.groupby('RID'):
    if len(rows.index) == 1:
        continue
    rows = rows.sort_values(by='EXAMDATE')
    times = list(rows['EXAMDATE'])
    years = [((t-times[0]).days)/365.0 for t in times]

    brainstem = rows['BRAIN_STEM']
    suvr_slope = slope(zip(years,list(rows['SUMMARYSUVR_WHOLECEREBNORM'])))
    bl_suvr = rows.iloc[0]['SUMMARYSUVR_WHOLECEREBNORM']

    slope_data = {'WHOLECEREBELLUM':slope(zip(years,list(rows['WHOLECEREBELLUM'].divide(brainstem)))),
                  'FRONTAL':slope(zip(years,list(rows['FRONTAL'].divide(brainstem)))),
                  'CINGULATE':slope(zip(years,list(rows['CINGULATE'].divide(brainstem)))),
                  'PARIETAL':slope(zip(years,list(rows['PARIETAL'].divide(brainstem)))),
                  'TEMPORAL':slope(zip(years,list(rows['TEMPORAL'].divide(brainstem)))),
                  'COMPOSITE':slope(zip(years,list(rows['COMPOSITE'].divide(brainstem)))),
                  'SUMMARYSUVR_WHOLECEREBNORM':slope(zip(years,list(rows['SUMMARYSUVR_WHOLECEREBNORM']))),
                  'INCREASER': 1 if suvr_slope > 0 else 0,
                  'BL_SUVR': bl_suvr,
                  'BL_POSITIVE': 1 if bl_suvr > 1.11 else 0,
                  'RID': rid}
    slopes.append(slope_data)

slopes_df = pd.DataFrame(slopes)
slopes_df = slopes_df[slopes_df.INCREASER==0]
sns.lmplot(x='WHOLECEREBELLUM', y='COMPOSITE', data=slopes_df, hue='INCREASER', fit_reg=False)

sns.plt.show()
