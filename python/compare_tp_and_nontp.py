import os
import sys

from scipy.stats import linregress
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from utils import parseDate


NONTP_FILE = '../output/02_19_16/UCBERKELEYAV45_02_19_16_merged_nontp.csv'
TP_FILE = '../output/02_19_16/UCBERKELEYAV45_02_19_16_merged_tp.csv'
MASTER_FILE = '../FDG_AV45_COGdata/FDG_AV45_COGdata_02_22_16.csv'

def convertToLongitudinalChange(av45_df, group_key, time_key, long_keys):
    '''
    Converts dataset into longitudinal AV45 and volumetric change, by subject
    Assumes rows indexed by subject ID
    '''
    slopes = []
    av45_df.loc[:,time_key] = av45_df.loc[:,time_key].apply(parseDate)
    for rid, rows in av45_df.groupby(group_key):
        if len(rows.index) == 1:
            continue
        sorted_rows = rows.sort_values(time_key)
        bl_time = sorted_rows.iloc[0][time_key]
        yrs_after_bl = [(_-bl_time).days/365.0 for _ in sorted_rows[time_key]]
        slope_values = {group_key: rid}
        for lk in long_keys:
            slope, intercept, r, p, stderr = linregress(yrs_after_bl, list(sorted_rows[lk]))
            slope_values[lk] = slope
        slopes.append(slope_values)
    slopes_df = pd.DataFrame(slopes)
    slopes_df.set_index(group_key,inplace=True)
    return slopes_df

def createMeanVarLabels(df, group_key, value_key):
    '''
    Calculates the mean/var of each group, and create labels to be used for plotting
    '''
    dfg = df.groupby(group_key)
    mean_vars = [(k, np.mean(v[value_key]),np.var(v[value_key])) for k,v in dfg]
    labels = ['%s\nMean: %.3f\nVar: %.3f' % (k,m,v) for k,m,v in mean_vars]
    return labels


if __name__ == "__main__":
    nontp_df = pd.read_csv(NONTP_FILE)
    tp_df = pd.read_csv(TP_FILE)
    master_df = pd.read_csv(MASTER_FILE, low_memory=False, header=[0,1])
    master_df.columns = master_df.columns.get_level_values(1)
    master_df.set_index('RID', inplace=True)

    # add some fields
    nontp_df['HIPPOCAMPUS_SIZE'] = nontp_df['LEFT_HIPPOCAMPUS_SIZE'] + nontp_df['RIGHT_HIPPOCAMPUS_SIZE']
    tp_df['HIPPOCAMPUS_SIZE'] = tp_df['LEFT_HIPPOCAMPUS_SIZE'] + tp_df['RIGHT_HIPPOCAMPUS_SIZE']

    # convert to annualized slopes
    long_keys = ['SUMMARYSUVR_WHOLECEREBNORM','HIPPOCAMPUS_SIZE']
    nontp_long_df = convertToLongitudinalChange(nontp_df,'RID','EXAMDATE',long_keys)
    tp_long_df = convertToLongitudinalChange(tp_df,'RID','EXAMDATE',long_keys)

    # stratify by BL Dx and positivity
    to_add_df = master_df[['Init_Diagnosis','AV45_wcereb_BIN1.11']]
    nontp_long_df = nontp_long_df.merge(to_add_df, how='left', left_index=True, right_index=True)
    tp_long_df = tp_long_df.merge(to_add_df, how='left', left_index=True, right_index=True)

    print nontp_long_df
    print tp_long_df

    # violin plots
    plot1 = sns.violinplot(x='Init_Diagnosis',y='SUMMARYSUVR_WHOLECEREBNORM',data=nontp_long_df)
    plot1.set_xticklabels(createMeanVarLabels(nontp_long_df, 'Init_Diagnosis', 'SUMMARYSUVR_WHOLECEREBNORM'))
    plt.title("Summary SUVR (wcereb ref) annualized slope, NON TIMEPOINT SPEC., by Dx")
    plt.show()

    plot2 = sns.violinplot(x='Init_Diagnosis',y='SUMMARYSUVR_WHOLECEREBNORM',data=tp_long_df)
    plot2.set_xticklabels(createMeanVarLabels(tp_long_df, 'Init_Diagnosis', 'SUMMARYSUVR_WHOLECEREBNORM'))
    plt.title("Summary SUVR (wcereb ref) annualized slope, TIMEPOINT SPEC., by Dx")
    plt.show()

    plot3 = sns.violinplot(x='AV45_wcereb_BIN1.11',y='SUMMARYSUVR_WHOLECEREBNORM',data=nontp_long_df)
    plot3.set_xticklabels(createMeanVarLabels(nontp_long_df, 'AV45_wcereb_BIN1.11', 'SUMMARYSUVR_WHOLECEREBNORM'))
    plt.title("Summary SUVR (wcereb ref) annualized slope, NON TIMEPOINT SPEC., by BL Pos")
    plt.show()

    plot4 = sns.violinplot(x='AV45_wcereb_BIN1.11',y='SUMMARYSUVR_WHOLECEREBNORM',data=tp_long_df)
    plot4.set_xticklabels(createMeanVarLabels(tp_long_df, 'AV45_wcereb_BIN1.11', 'SUMMARYSUVR_WHOLECEREBNORM'))
    plt.title("Summary SUVR (wcereb ref) annualized slope, TIMEPOINT SPEC., by BL Pos")
    plt.show()

    # direct comparison scatter plot
    suvr_merged_df = nontp_long_df[['SUMMARYSUVR_WHOLECEREBNORM']].merge(tp_long_df[['SUMMARYSUVR_WHOLECEREBNORM','HIPPOCAMPUS_SIZE']], how='inner', left_index=True, right_index=True, suffixes=('_nontp','_tp'))
    suvr_merged_df.plot(kind='scatter',x='SUMMARYSUVR_WHOLECEREBNORM_nontp',y='SUMMARYSUVR_WHOLECEREBNORM_tp')
    plt.show()

    # regression and add residuals
    slope, intercept, r, p, stderr = linregress(suvr_merged_df['SUMMARYSUVR_WHOLECEREBNORM_nontp'], suvr_merged_df['SUMMARYSUVR_WHOLECEREBNORM_tp'])
    p1 = np.poly1d([slope, intercept])
    suvr_merged_df['reg_fit'] = suvr_merged_df.loc[:,'SUMMARYSUVR_WHOLECEREBNORM_nontp'].apply(p1)
    suvr_merged_df['R'] = suvr_merged_df['SUMMARYSUVR_WHOLECEREBNORM_tp'] - suvr_merged_df['reg_fit']
    suvr_merged_df['R2'] = suvr_merged_df.loc[:,'R'].apply(np.square)

    # plot residuals
    suvr_merged_df.plot(kind='scatter',x='SUMMARYSUVR_WHOLECEREBNORM_nontp',y='R')
    plt.show()
    suvr_merged_df.plot(kind='scatter',x='SUMMARYSUVR_WHOLECEREBNORM_nontp',y='R2')
    plt.show()
    suvr_merged_df.plot(kind='scatter',x='HIPPOCAMPUS_SIZE',y='R2')
    plt.show()

