import os
import sys

from scipy.stats import linregress
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from utils import parseDate, isnan, FRONTAL, PARIETAL, TEMPORAL, CINGULATE, WHOLECEREBELLUM, importFreesurferLookup

LUT_FILE='../FreeSurferColorLUT.txt'

NONTP_FILE = '../output/02_19_16/UCBERKELEYAV45_02_19_16_merged_nontp.csv'
TP_FILE = '../output/02_19_16/UCBERKELEYAV45_02_19_16_merged_tp.csv'
MASTER_FILE = '../FDG_AV45_COGdata/FDG_AV45_COGdata_02_22_16.csv'


def extractLobeVolumes(av45_df):
    '''
    Sums up the voxel counts for regions that fall into
    frontal, parietal, temporal, cingulate
    '''
    lut = importFreesurferLookup(LUT_FILE, flip=False)
    frontal_columns = ['%s_SIZE' % lut[_].replace('-','_').upper() for _ in FRONTAL]
    parietal_columns = ['%s_SIZE' % lut[_].replace('-','_').upper() for _ in PARIETAL]
    temporal_columns = ['%s_SIZE' % lut[_].replace('-','_').upper() for _ in TEMPORAL]
    cingulate_columns = ['%s_SIZE' % lut[_].replace('-','_').upper() for _ in CINGULATE]
    cereb_columns = ['%s_SIZE' % lut[_].replace('-','_').upper() for _ in WHOLECEREBELLUM]
    av45_df['FRONTAL_SIZE'] = av45_df[frontal_columns].sum(axis=1)
    av45_df['PARIETAL_SIZE'] = av45_df[parietal_columns].sum(axis=1)
    av45_df['TEMPORAL_SIZE'] = av45_df[temporal_columns].sum(axis=1)
    av45_df['CINGULATE_SIZE'] = av45_df[cingulate_columns].sum(axis=1)
    av45_df['WHOLECEREBELLUM_SIZE'] = av45_df[cereb_columns].sum(axis=1)
    return av45_df

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
    labels = ['%s\nMean: %.5f\nVar: %.5f' % (k,m,v) for k,m,v in mean_vars]
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
    nontp_df = extractLobeVolumes(nontp_df)
    tp_df = extractLobeVolumes(tp_df)

    # convert to annualized slopes
    long_keys = ['SUMMARYSUVR_WHOLECEREBNORM',
                 'FRONTAL',
                 'PARIETAL',
                 'TEMPORAL',
                 'CINGULATE',
                 'WHOLECEREBELLUM',
                 'FRONTAL_SIZE',
                 'PARIETAL_SIZE',
                 'TEMPORAL_SIZE',
                 'CINGULATE_SIZE',
                 'WHOLECEREBELLUM_SIZE']
    nontp_long_df = convertToLongitudinalChange(nontp_df,'RID','EXAMDATE',long_keys)
    tp_long_df = convertToLongitudinalChange(tp_df,'RID','EXAMDATE',long_keys)

    # stratify by BL Dx and positivity
    def sum_dx_and_bl_pos(row):
        dx = row['Init_Diagnosis']
        pos = row['AV45_wcereb_BIN1.11']
        if isnan(dx) or isnan(pos):
            return np.nan
        if float(pos) == 0.0:
            new_val = '%s_BLNeg' % (dx,)
        elif float(pos) == 1.0:
            new_val = '%s_BLPos' % (dx,)
        else:
            raise Exception("Unknown positivity: %s" % pos)
        return new_val

    master_df['Dx_and_BLPos'] = master_df.apply(sum_dx_and_bl_pos, axis=1)
    to_add_df = master_df[['Init_Diagnosis','AV45_wcereb_BIN1.11','Dx_and_BLPos']]
    nontp_long_df = nontp_long_df.merge(to_add_df, how='left', left_index=True, right_index=True)
    tp_long_df = tp_long_df.merge(to_add_df, how='left', left_index=True, right_index=True)

    # direct comparison scatter plot
    uptake_keys = ['SUMMARYSUVR_WHOLECEREBNORM', 'FRONTAL', 'PARIETAL', 'TEMPORAL', 'CINGULATE', 'WHOLECEREBELLUM']
    suvr_merged_df = nontp_long_df[uptake_keys].merge(tp_long_df, how='inner', left_index=True, right_index=True, suffixes=('_nontp','_tp'))

    for uk in uptake_keys:
        # regression and add residuals
        print uk
        uk_nontp = '%s_nontp' % uk
        uk_tp = '%s_tp' % uk
        slope, intercept, r, p, stderr = linregress(suvr_merged_df[uk_nontp], suvr_merged_df[uk_tp])
        print '\tR2: %s' % (r**2,)
        print '\ty ~ %sx + %s' % (slope, intercept)
        p1 = np.poly1d([slope, intercept])
        suvr_merged_df['%s_regfit' % uk] = suvr_merged_df.loc[:,uk_nontp].apply(p1)
        suvr_merged_df['%s_R' % uk] = suvr_merged_df[uk_tp] - suvr_merged_df['%s_regfit' % uk]

    # weird = suvr_merged_df[suvr_merged_df.R2>0.08]
    # print weird

    ylim = (-0.2,0.2)


    # violinplots
    plt.figure(1)
    plot1 = sns.violinplot(x='Dx_and_BLPos',y='SUMMARYSUVR_WHOLECEREBNORM',data=nontp_long_df)
    plot1.set_xticklabels(createMeanVarLabels(nontp_long_df, 'Dx_and_BLPos', 'SUMMARYSUVR_WHOLECEREBNORM'))
    plt.title("Summary SUVR (wcereb ref) annualized slope, NON TIMEPOINT SPEC., by Dx")

    plt.figure(2)
    plot2 = sns.violinplot(x='Dx_and_BLPos',y='SUMMARYSUVR_WHOLECEREBNORM',data=tp_long_df)
    plot2.set_xticklabels(createMeanVarLabels(tp_long_df, 'Dx_and_BLPos', 'SUMMARYSUVR_WHOLECEREBNORM'))
    plt.title("Summary SUVR (wcereb ref) annualized slope, TIMEPOINT SPEC., by Dx")
    plt.show()

    # scatter plots
    scatter_pairs = [('FRONTAL_SIZE', 'SUMMARYSUVR_WHOLECEREBNORM_R'),
                     ('PARIETAL_SIZE', 'SUMMARYSUVR_WHOLECEREBNORM_R'),
                     ('TEMPORAL_SIZE', 'SUMMARYSUVR_WHOLECEREBNORM_R'),
                     ('CINGULATE_SIZE', 'SUMMARYSUVR_WHOLECEREBNORM_R'),
                     ('WHOLECEREBELLUM_SIZE', 'SUMMARYSUVR_WHOLECEREBNORM_R'),
                     ('SUMMARYSUVR_WHOLECEREBNORM_nontp', 'SUMMARYSUVR_WHOLECEREBNORM_tp'),
                     ('SUMMARYSUVR_WHOLECEREBNORM_nontp', 'SUMMARYSUVR_WHOLECEREBNORM_R')]
    scatter_pairs = [('FRONTAL_SIZE', 'FRONTAL_R'),
                     ('PARIETAL_SIZE', 'PARIETAL_R'),
                     ('TEMPORAL_SIZE', 'TEMPORAL_R'),
                     ('CINGULATE_SIZE', 'CINGULATE_R'),
                     ('WHOLECEREBELLUM_SIZE', 'WHOLECEREBELLUM_R')]
    for fignum, (xkey,ykey) in enumerate(scatter_pairs):
        plt.figure(fignum)
        suvr_merged_df[[xkey,ykey]].plot(kind='scatter', x=xkey, y=ykey, ylim=ylim)
        plt.title('%s vs %s, Annualized Change' % (xkey,ykey))

    plt.show()

