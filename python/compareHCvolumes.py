'''
Compare ICV and HC volumes from our freesurfer results, versus what's reported by UCSF's crosssectional results
'''

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')

mprage_file_0 = '../mr_docs/UCSF/cross_section/UCSFFSX_11_02_15.csv'
mprage_file_1 = '../mr_docs/UCSF/cross_section/UCSFFSX51_11_02_15.csv'
mprage_file_2 = '../mr_docs/UCSF/cross_section/UCSFFSX51_ADNI1_3T_11_02_15.csv'
mprage_df0 = pd.read_csv(mprage_file_0)
mprage_df1 = pd.read_csv(mprage_file_1)
mprage_df2 = pd.read_csv(mprage_file_2)
mprage_df = pd.concat((mprage_df0, mprage_df1, mprage_df2))
mprage_df.set_index('IMAGEUID',inplace=True)
raw_numrows = float(len(mprage_df.index))
mprage_df = mprage_df[mprage_df.OVERALLQC!='Pass']
filtered_numrows = float(len(mprage_df.index))

print 'Passed: %s' % (filtered_numrows/raw_numrows,)

df = pd.read_csv('../adni_av45_freesurfer_nontp_aseg_stats.csv')

comparisons = []
for i, row in df.iterrows():
    mris = eval(row['MRI_inputs'])
    mris = [int(_.replace('I','')) for _ in mris]
    icv = row['EstimatedTotalIntraCranialVol']
    left_hcv = row['Left-Hippocampus']
    right_hcv = row['Right-Hippocampus']
    hcv = left_hcv + right_hcv
    for mri in mris:
        if mri in mprage_df.index:
            ucsf_left_hcv = mprage_df.loc[mri,'ST29SV']
            ucsf_right_hcv = mprage_df.loc[mri,'ST88SV']
            ucsf_hcv = ucsf_left_hcv + ucsf_right_hcv
            ucsf_icv = mprage_df.loc[mri,'ST10CV']
            new_row = {'IMAGEUID': mri, 'icv': icv, 'lhcv': left_hcv, 'rhcv': right_hcv, 'hcv': hcv, 'ucsf_icv': ucsf_icv, 'ucsf_hcv': ucsf_hcv, 'ucsf_lhcv': ucsf_left_hcv, 'ucsf_rhcv': ucsf_right_hcv}
            comparisons.append(new_row)

comp_df = pd.DataFrame(comparisons)
comp_df.dropna(inplace=True)

comp_df.plot(kind='scatter',x='hcv',y='ucsf_hcv')
comp_df.plot(kind='scatter',x='lhcv',y='ucsf_lhcv')
comp_df.plot(kind='scatter',x='rhcv',y='ucsf_rhcv')
comp_df.plot(kind='scatter',x='icv',y='ucsf_icv')
plt.show()

hcv_ols = pd.stats.api.ols(y=comp_df.ucsf_hcv, x=comp_df.hcv)
lhcv_ols = pd.stats.api.ols(y=comp_df.ucsf_lhcv, x=comp_df.lhcv)
rhcv_ols = pd.stats.api.ols(y=comp_df.ucsf_rhcv, x=comp_df.rhcv)
icv_ols = pd.stats.api.ols(y=comp_df.ucsf_icv, x=comp_df.icv)


print 'HCV:'
print '\tR2: %s' % hcv_ols.r2
print '\ty ~ %sx + %s' % tuple(hcv_ols.beta)

print 'LHCV'
print '\tR2: %s' % lhcv_ols.r2
print '\ty ~ %sx + %s' % tuple(lhcv_ols.beta)

print 'RHCV:'
print '\tR2: %s' % rhcv_ols.r2
print '\ty ~ %sx + %s' % tuple(rhcv_ols.beta)

print 'ICV:'
print '\tR2: %s' % icv_ols.r2
print '\ty ~ %sx + %s' % tuple(icv_ols.beta)