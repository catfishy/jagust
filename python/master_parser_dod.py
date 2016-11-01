import re
from collections import defaultdict, Counter
from datetime import datetime, timedelta
from scipy import stats, optimize
import numpy as np

from utils import *

def syncAPOEData(master_df, apoe_file):
    apoe_df = importAPOE(apoe_file, as_df=True)

    headers = ['APOE4BIN']
    def extraction_fn(scrno, subj_rows):
        row = subj_rows.iloc[0]
        data = {'SCRNO': scrno,
                'APOE4BIN': 1 if (row['APGEN1'] == 4 or row['APGEN2'] == 4) else 0}
        return pd.DataFrame([data]).set_index('SCRNO')

    parsed_df = parseSubjectGroups(apoe_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncDemogData(master_df, demog_file):
    demog_df = importDemog(demog_file, as_df=True)

    headers = ['Sex', 'Age', 'Edu']
    def extraction_fn(scrno, subj_rows):
        row = subj_rows.iloc[0]
        if row['PTGENDER'] == 1:
            sex = 'M'
        elif row['PTGENDER'] == 2:
            sex = 'F'
        else:
            raise Exception("Bad Gender: %s" % row['PTGENDER'])
        data = {'SCRNO': scrno,
                'Sex': sex,
                'Age': row['PTAGE'],
                'Edu': row['PTEDUCAT']}
        return pd.DataFrame([data]).set_index('SCRNO')

    parsed_df = parseSubjectGroups(demog_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after='Notes', restrict=True)
    return master_df

def syncCAPSData(master_df, caps_curr_file, caps_lifetime_file):
    caps_df = importCAPS(caps_curr_file, caps_lifetime_file, as_df=True)

    headers = ['CAPS_CURRSCORE', 'CAPS_LIFETIME_SCORE']
    def extraction_fn(scrno, subj_rows):
        row = subj_rows.iloc[0]
        data = {'SCRNO': scrno,
                'CAPS_CURRSCORE': row['curr'],
                'CAPS_LIFETIME_SCORE': row['life']}
        return pd.DataFrame([data]).set_index('SCRNO')

    parsed_df = parseSubjectGroups(caps_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncWMHData(master_df, wmh_file):
    wmh_df = importWMH(wmh_file, as_df=True)
    timepoints = max(Counter(wmh_df.index).values())

    # create header order
    row_indices = range(1,timepoints+1)
    headers = []
    headers += ['WMH_WHITMATHYP_%s' % i for i in row_indices]
    headers += ['WMH_percentOfICV_%s' % i for i in row_indices]
    headers += ['WMH_DATE_%s' % i for i in row_indices]
    headers += ['WMH_postAV45_%s' % i for i in row_indices]
    headers += ['WMH_WHITMATHYP_postAV45_SLOPE','WMH_WHITMATHYP_closest_AV45_BL','WMH_WHITMATHYP_closest_AV45_2']
    headers += ['WMH_percentOfICV_postAV45_SLOPE','WMH_percentOfICV_closest_AV45_BL','WMH_percentOfICV_closest_AV45_2']

    def extraction_fn(scrno, subj_rows):
        av45_date1, av45_date2 = getAV45Dates(scrno, master_df)
        subj_rows['yrDiff'] = subj_rows['EXAMDATE'].apply(lambda x: yrDiff(x,av45_date1))
        # get longitudinal measurements
        subj_rows['SCRNO'] = scrno
        wmh_percent_long = groupLongPivot(subj_rows, 'SCRNO','wmh_percent','WMH_percentOfICV_')
        wmh_hyp_long = groupLongPivot(subj_rows, 'SCRNO','wmh','WMH_WHITMATHYP_')
        wmh_date_long = groupLongPivot(subj_rows, 'SCRNO','EXAMDATE','WMH_DATE_')
        wmh_post_av45 = groupLongPivot(subj_rows, 'SCRNO','yrDiff','WMH_postAV45_')
        data_df = wmh_percent_long.merge(wmh_hyp_long,left_index=True,right_index=True)
        data_df = data_df.merge(wmh_date_long,left_index=True,right_index=True)
        data_df = data_df.merge(wmh_post_av45,left_index=True,right_index=True)
        # get slope
        data_df.loc[scrno,'WMH_WHITMATHYP_postAV45_SLOPE'] = groupSlope(subj_rows,'EXAMDATE','wmh',cutoff_date=av45_date1-timedelta(days=90)) if not isnan(av45_date1) else np.nan
        data_df.loc[scrno,'WMH_percentOfICV_postAV45_SLOPE'] = groupSlope(subj_rows,'EXAMDATE','wmh_percent',cutoff_date=av45_date1-timedelta(days=90)) if not isnan(av45_date1) else np.nan
        # get closest
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'wmh', [av45_date1, av45_date2], day_limit=365/2)
        data_df.loc[scrno,'WMH_WHITMATHYP_closest_AV45_BL'] = closest_vals[0]
        data_df.loc[scrno,'WMH_WHITMATHYP_closest_AV45_2'] = closest_vals[1]
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'wmh_percent', [av45_date1, av45_date2], day_limit=365/2)
        data_df.loc[scrno,'WMH_percentOfICV_closest_AV45_BL'] = closest_vals[0]
        data_df.loc[scrno,'WMH_percentOfICV_closest_AV45_2'] = closest_vals[1]
        return data_df

    parsed_df = parseSubjectGroups(wmh_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncMRIData(master_df, mri_file):
    mri_df = importDODMRI(mri_file, as_df=True)

    headers = ["MRIDATE"]
    def extraction_fn(scrno, subj_rows):
        subj_rows.sort_values(by='EXAMDATE', inplace=True)
        row = subj_rows.iloc[0]
        data = {'SCRNO': scrno,
                'MRIDATE': row['EXAMDATE']}
        return pd.DataFrame([data]).set_index('SCRNO')

    parsed_df = parseSubjectGroups(mri_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df


def syncGDData(master_df, gd_file, registry):
    gd_df = importGD(gd_file, registry=registry, as_df=True)
    timepoints = max(Counter(gd_df.index).values())

    # create header order
    row_indices = range(1,timepoints+1)
    headers = []
    headers += ['GD_%s' % i for i in row_indices]
    headers += ['GD_DATE_%s' % i for i in row_indices]
    headers += ['GD_postAV45_%s' % i for i in row_indices]
    headers += ['GD_postAV45_SLOPE','GD_closest_AV45_BL','GD_closest_AV45_2']

    def extraction_fn(scrno, subj_rows):
        av45_date1, av45_date2 = getAV45Dates(scrno, master_df)
        subj_rows['yrDiff'] = subj_rows['EXAMDATE'].apply(lambda x: yrDiff(x,av45_date1))
        # get longitudinal measurements
        subj_rows['SCRNO'] = scrno
        gd_long = groupLongPivot(subj_rows, 'SCRNO','GDTOTAL','GD_')
        date_long = groupLongPivot(subj_rows, 'SCRNO','EXAMDATE','GD_DATE_')
        yrdiff_long = groupLongPivot(subj_rows, 'SCRNO','yrDiff','GD_postAV45_')
        data_df = gd_long.merge(date_long,left_index=True,right_index=True)
        data_df = data_df.merge(yrdiff_long,left_index=True,right_index=True)
        # get slope
        av45_date1, av45_date2 = getAV45Dates(scrno, master_df)
        data_df.loc[scrno,'GD_postAV45_SLOPE'] = groupSlope(subj_rows,'EXAMDATE','GDTOTAL',cutoff_date=av45_date1-timedelta(days=90)) if not isnan(av45_date1) else np.nan
        # get closest
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'GDTOTAL', [av45_date1, av45_date2], day_limit=365/2)
        data_df.loc[scrno,'GD_closest_AV45_BL'] = closest_vals[0]
        data_df.loc[scrno,'GD_closest_AV45_2'] = closest_vals[1]
        return data_df

    parsed_df = parseSubjectGroups(gd_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncAntidepData(master_df, backmeds_file, registry):
    antidep_df = importDODAntidep(backmeds_file, registry, as_df=True)
    headers = ['ANTIDEP_USE','WHICH_ANTIDEP','SSRI']

    def extraction_fn(scrno, subj_rows):
        subj_rows.sort_values(by='EXAMDATE', inplace=True)
        bl = subj_rows.iloc[0]
        new_data = {'SCRNO': scrno,
                    'ANTIDEP_USE': 1 if bl['ANTIDEP_USE'] else 0,
                    'WHICH_ANTIDEP': bl['WHICH_ANTIDEP'],
                    'SSRI': 1 if bl['SSRI'] else 0}
        return pd.DataFrame([new_data]).set_index('SCRNO')

    parsed_df = parseSubjectGroups(antidep_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncCSFData(master_df, csf_file, registry):
    csf_df = importCSF([csf_file], registry, as_df=True)

    headers = ['CSF_abeta', 'CSF_tau', 'CSF_ptau']
    def extraction_fn(scrno, subj_rows):
        subj_rows.sort_values(by='EXAMDATE', inplace=True)
        row = subj_rows.iloc[0]
        data = {'SCRNO': scrno,
                'CSF_abeta': row['ABETA'],
                'CSF_tau': row['TAU'],
                'CSF_ptau': row['PTAU']}
        return pd.DataFrame([data]).set_index('SCRNO')

    parsed_df = parseSubjectGroups(csf_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncADASData(master_df, adas_file, registry):
    adas_df = importADASCog(None, adas_file, registry, as_df=True)
    timepoints = max(Counter(adas_df.index).values())

    # create header order
    row_indices = range(1,timepoints+1)
    headers = []
    headers += ['ADAS_%s' % i for i in row_indices]
    headers += ['ADAS_DATE_%s' % i for i in row_indices]
    headers += ['ADAS_postAV45_%s' % i for i in row_indices]
    headers += ['ADAS_postAV45_SLOPE','ADAS_closest_AV45_BL','ADAS_closest_AV45_2']

    def extraction_fn(scrno, subj_rows):
        av45_date1, av45_date2 = getAV45Dates(scrno, master_df)
        subj_rows['yrDiff'] = subj_rows['EXAMDATE'].apply(lambda x: yrDiff(x,av45_date1))
        # get longitudinal measurements
        subj_rows['SCRNO'] = scrno
        adas_long = groupLongPivot(subj_rows, 'SCRNO','TOTSCORE','ADAS_')
        adas_date = groupLongPivot(subj_rows, 'SCRNO','EXAMDATE','ADAS_DATE_')
        yrdiff_long = groupLongPivot(subj_rows,'SCRNO','yrDiff','ADAS_postAV45_')
        data_df = adas_long.merge(adas_date,left_index=True,right_index=True)
        data_df = data_df.merge(yrdiff_long,left_index=True,right_index=True)
        # get slope
        data_df.loc[scrno,'ADAS_postAV45_SLOPE'] = groupSlope(subj_rows,'EXAMDATE','TOTSCORE',cutoff_date=av45_date1-timedelta(days=90)) if not isnan(av45_date1) else np.nan
        # get closest
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'TOTSCORE', [av45_date1, av45_date2], day_limit=365/2)
        data_df.loc[scrno,'ADAS_closest_AV45_BL'] = closest_vals[0]
        data_df.loc[scrno,'ADAS_closest_AV45_2'] = closest_vals[1]
        return data_df

    parsed_df = parseSubjectGroups(adas_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df


def syncFAQData(master_df, faq_file, registry):
    faq_df = importFAQ(faq_file, registry=registry, as_df=True)
    timepoints = max(Counter(faq_df.index).values())

    # create header order
    row_indices = range(1,timepoints+1)
    headers = []
    headers += ['FAQ_%s' % i for i in row_indices]
    headers += ['FAQ_DATE_%s' % i for i in row_indices]
    headers += ['FAQ_postAV45_%s' % i for i in row_indices]
    headers += ['FAQ_postAV45_SLOPE', 'FAQ_closest_AV45_BL', 'FAQ_closest_AV45_2']

    def extraction_fn(scrno, subj_rows):
        av45_date1, av45_date2 = getAV45Dates(scrno, master_df)
        subj_rows['yrDiff'] = subj_rows['EXAMDATE'].apply(lambda x: yrDiff(x,av45_date1))
        # get longitudinal measurements
        subj_rows['SCRNO'] = scrno
        faq_long = groupLongPivot(subj_rows, 'SCRNO','FAQTOTAL','FAQ_')
        faq_date = groupLongPivot(subj_rows, 'SCRNO','EXAMDATE','FAQ_DATE_')
        yrdiff_long = groupLongPivot(subj_rows, 'SCRNO','yrDiff','FAQ_postAV45_')
        data_df = faq_long.merge(faq_date,left_index=True,right_index=True)
        data_df = data_df.merge(yrdiff_long,left_index=True,right_index=True)
        # get slope
        data_df.loc[scrno,'FAQ_postAV45_SLOPE'] = groupSlope(subj_rows,'EXAMDATE','FAQTOTAL',cutoff_date=av45_date1-timedelta(days=90)) if not isnan(av45_date1) else np.nan
        # get closest
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'FAQTOTAL', [av45_date1, av45_date2], day_limit=365/2)
        data_df.loc[scrno,'FAQ_closest_AV45_BL'] = closest_vals[0]
        data_df.loc[scrno,'FAQ_closest_AV45_2'] = closest_vals[1]
        return data_df

    parsed_df = parseSubjectGroups(faq_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncDTI(master_df, dti_file, registry):
    df = importDTI(dti_file, registry=registry, as_df=True)
    timepoints = max(Counter(df.index).values())
    fa_cols = [u'FA_CST_L', u'FA_CST_R', u'FA_ICP_L',
               u'FA_ICP_R', u'FA_ML_L', u'FA_ML_R', u'FA_SCP_L', u'FA_SCP_R',
               u'FA_CP_L', u'FA_CP_R', u'FA_ALIC_L', u'FA_ALIC_R', u'FA_PLIC_L',
               u'FA_PLIC_R', u'FA_PTR_L', u'FA_PTR_R', u'FA_ACR_L', u'FA_ACR_R',
               u'FA_SCR_L', u'FA_SCR_R', u'FA_PCR_L', u'FA_PCR_R', u'FA_CGC_L',
               u'FA_CGC_R', u'FA_CGH_L', u'FA_CGH_R', u'FA_FX_ST_L', u'FA_FX_ST_R',
               u'FA_SLF_L', u'FA_SLF_R', u'FA_SFO_L', u'FA_SFO_R', u'FA_SS_L',
               u'FA_SS_R', u'FA_EC_L', u'FA_EC_R', u'FA_UNC_L', u'FA_UNC_R', u'FA_FX',
               u'FA_GCC', u'FA_BCC', u'FA_SCC', u'FA_RLIC_L', u'FA_RLIC_R']
    headers = ['FA_DATE_%s' % (i+1,) for i in range(timepoints)]
    for fa_col in fa_cols:
        headers += ['%s_%s' % (fa_col,i+1) for i in range(timepoints)]

    def extraction_fn(scrno, subj_rows):
        # average rows together if same viscode/EXAMDATE
        subj_rows = subj_rows.groupby(['VISCODE','EXAMDATE']).aggregate(np.mean)
        subj_rows = subj_rows.reset_index().sort_values('EXAMDATE')
        subj_rows['SCRNO'] = scrno
        data_df = groupLongPivot(subj_rows, 'SCRNO','EXAMDATE','FA_DATE_')
        for fa_col in fa_cols:
            cur_long_df = groupLongPivot(subj_rows, 'SCRNO',fa_col,'%s_' % (fa_col,))
            data_df = data_df.merge(cur_long_df,left_index=True,right_index=True)
        return data_df

    parsed_df = parseSubjectGroups(df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncTBIAGE(master_df, tbiage_file, registry):
    df = importTBIAGE(tbiage_file, registry=registry, as_df=True)
    headers = ['TBINJAGE']

    def extraction_fn(scrno, subj_rows):
        subj_rows['SCRNO'] = scrno
        cur_df = subj_rows.head(1)
        cur_df.set_index('SCRNO',inplace=True)
        cur_df = cur_df[['TBINJAGE']]
        return cur_df

    parsed_df = parseSubjectGroups(df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncSleep(master_df, sleep_file, registry):
    df = importSleepQuality(sleep_file, registry=registry, as_df=True)
    timepoints = max(Counter(df.index).values())
    headers = ['PSQI.GLOBAL_%s' % (i+1,) for i in range(timepoints)]
    headers += ['PSQI.DATE_%s' % (i+1,) for i in range(timepoints)]
    headers += ['PSQI.GLOBAL_closest_AV45_BL','PSQI.GLOBAL_closest_AV45_Scan2',
                'PSQI.GLOBAL_closest_AV1451_BL','PSQI.GLOBAL_closest_AV1451_Scan2']

    def extraction_fn(scrno, subj_rows):
        av45_date1, av45_date2 = getAV45Dates(scrno, master_df)
        av1451_date1, av1451_date2 = getAV1451Dates(scrno, master_df)
        # get longitudinal measurements
        subj_rows['SCRNO'] = scrno
        psqi_long = groupLongPivot(subj_rows, 'SCRNO','GLOBAL','PSQI.GLOBAL_')
        psqi_date = groupLongPivot(subj_rows, 'SCRNO','EXAMDATE','PSQI.DATE_')
        data_df = psqi_long.merge(psqi_date,left_index=True,right_index=True)
        # get closest av45
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'GLOBAL', [av45_date1, av45_date2], day_limit=365)
        print "%s Dates: %s, %s" % (scrno,(av45_date1,av45_date2),(av1451_date1,av1451_date2))
        # print closest_vals
        data_df.loc[scrno,'PSQI.GLOBAL_closest_AV45_BL'] = closest_vals[0]
        data_df.loc[scrno,'PSQI.GLOBAL_closest_AV45_Scan2'] = closest_vals[1]
        # get closest av1451
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'GLOBAL', [av1451_date1, av1451_date2], day_limit=365)
        data_df.loc[scrno,'PSQI.GLOBAL_closest_AV1451_BL'] = closest_vals[0]
        data_df.loc[scrno,'PSQI.GLOBAL_closest_AV1451_Scan2'] = closest_vals[1]
        return data_df

    parsed_df = parseSubjectGroups(df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncASIData(master_df, asi_file, registry):
    asi_df = importASI(asi_file, registry, as_df=True)

    headers = ['ASI_1A','ASI_2B']
    def extraction_fn(scrno, subj_rows):
        subj_rows.sort_values(by='EXAMDATE',inplace=True)
        row = subj_rows.iloc[0]
        data = {'SCRNO': scrno,
                'ASI_1A': row['ASI1A'],
                'ASI_2B': row['ASI2B']}
        return pd.DataFrame([data]).set_index('SCRNO')

    parsed_df = parseSubjectGroups(asi_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncCDRData(master_df, cdr_file, registry):
    cdr_df = importCDR(cdr_file, registry=registry, as_df=True)
    timepoints = max(Counter(cdr_df.index).values())

    # create header order
    row_indices = range(1,timepoints+1)
    headers = []
    headers += ['CDR_GLOBAL_%s' % i for i in row_indices]
    headers += ['CDR_SOB_%s' % i for i in row_indices]
    headers += ['CDR_DATE_%s' % i for i in row_indices]
    headers += ['CDR_postAV45_%s' % i for i in row_indices]
    headers += ['CDR_GLOBAL_postAV45_SLOPE', 'CDR_GLOBAL_closest_AV45_BL', 'CDR_GLOBAL_closest_AV45_2']
    headers += ['CDR_SOB_postAV45_SLOPE', 'CDR_SOB_closest_AV45_BL', 'CDR_SOB_closest_AV45_2']

    def extraction_fn(scrno, subj_rows):
        av45_date1, av45_date2 = getAV45Dates(scrno, master_df)
        subj_rows['yrDiff'] = subj_rows['EXAMDATE'].apply(lambda x: yrDiff(x,av45_date1))
        # get longitudinal measurements
        subj_rows['SCRNO'] = scrno
        cdr_long = groupLongPivot(subj_rows, 'SCRNO','CDGLOBAL','CDR_GLOBAL_')
        cdr_sum_long = groupLongPivot(subj_rows, 'SCRNO','CDSOB','CDR_SOB_')
        cdr_date = groupLongPivot(subj_rows, 'SCRNO','EXAMDATE','CDR_DATE_')
        yrdiff_long = groupLongPivot(subj_rows, 'SCRNO','yrDiff','CDR_postAV45_')
        data_df = cdr_long.merge(cdr_sum_long,left_index=True,right_index=True)
        data_df = data_df.merge(cdr_date,left_index=True,right_index=True)
        data_df = data_df.merge(yrdiff_long,left_index=True,right_index=True)
        # get slope
        data_df.loc[scrno,'CDR_GLOBAL_postAV45_SLOPE'] = groupSlope(subj_rows,'EXAMDATE','CDGLOBAL',cutoff_date=av45_date1-timedelta(days=90)) if not isnan(av45_date1) else np.nan
        data_df.loc[scrno,'CDR_SOB_postAV45_SLOPE'] = groupSlope(subj_rows,'EXAMDATE','CDSOB',cutoff_date=av45_date1-timedelta(days=90)) if not isnan(av45_date1) else np.nan
        # get closest
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'CDGLOBAL', [av45_date1, av45_date2], day_limit=365/2)
        data_df.loc[scrno,'CDR_GLOBAL_closest_AV45_BL'] = closest_vals[0]
        data_df.loc[scrno,'CDR_GLOBAL_closest_AV45_2'] = closest_vals[1]
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'CDSOB', [av45_date1, av45_date2], day_limit=365/2)
        data_df.loc[scrno,'CDR_SOB_closest_AV45_BL'] = closest_vals[0]
        data_df.loc[scrno,'CDR_SOB_closest_AV45_2'] = closest_vals[1]
        return data_df

    parsed_df = parseSubjectGroups(cdr_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncMMSEData(master_df, mmse_file, registry):
    mmse_df = importMMSE(mmse_file, registry=registry, as_df=True)
    timepoints = max(Counter(mmse_df.index).values())

    # create header order
    row_indices = range(1,timepoints+1)
    headers = []
    headers += ['MMSE_%s' % i for i in row_indices]
    headers += ['MMSE_DATE_%s' % i for i in row_indices]
    headers += ['MMSE_postAV45_%s' % i for i in row_indices]
    headers += ['MMSE_postAV45_SLOPE', 'MMSE_closest_AV45_BL', 'MMSE_closest_AV45_2']

    def extraction_fn(scrno, subj_rows):
        av45_date1, av45_date2 = getAV45Dates(scrno, master_df)
        subj_rows['yrDiff'] = subj_rows['EXAMDATE'].apply(lambda x: yrDiff(x,av45_date1))
        # get longitudinal measurements
        subj_rows['SCRNO'] = scrno
        mmse_long = groupLongPivot(subj_rows, 'SCRNO','MMSCORE','MMSE_')
        mmse_date = groupLongPivot(subj_rows, 'SCRNO','EXAMDATE','MMSE_DATE_')
        yrdiff_long = groupLongPivot(subj_rows, 'SCRNO','yrDiff','MMSE_postAV45_')
        data_df = mmse_long.merge(mmse_date,left_index=True,right_index=True)
        data_df = data_df.merge(yrdiff_long,left_index=True,right_index=True)
        # get slope
        data_df.loc[scrno,'MMSE_postAV45_SLOPE'] = groupSlope(subj_rows,'EXAMDATE','MMSCORE',cutoff_date=av45_date1-timedelta(days=90)) if not isnan(av45_date1) else np.nan
        # get closest
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'MMSCORE', [av45_date1, av45_date2], day_limit=365/2)
        data_df.loc[scrno,'MMSE_closest_AV45_BL'] = closest_vals[0]
        data_df.loc[scrno,'MMSE_closest_AV45_2'] = closest_vals[1]
        return data_df

    parsed_df = parseSubjectGroups(mmse_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncAVLTData(master_df, avlt_file, registry):
    avlt_df = importAVLT(avlt_file, registry=registry, as_df=True)
    timepoints = max(Counter(avlt_df.index).values())

    # create header order
    row_indices = range(1,timepoints+1)
    headers = []
    headers += ['AVLT_%s' % i for i in row_indices]
    headers += ['AVLT_DATE_%s' % i for i in row_indices]
    headers += ['AVLT_postAV45_%s' % i for i in row_indices]
    headers += ['AVLT_postAV45_SLOPE', 'AVLT_closest_AV45_BL', 'AVLT_closest_AV45_2']

    def extraction_fn(scrno, subj_rows):
        av45_date1, av45_date2 = getAV45Dates(scrno, master_df)
        subj_rows['yrDiff'] = subj_rows['EXAMDATE'].apply(lambda x: yrDiff(x,av45_date1))
        # get longitudinal measurements
        subj_rows['SCRNO'] = scrno
        avlt_long = groupLongPivot(subj_rows, 'SCRNO','TOTS','AVLT_')
        avlt_date = groupLongPivot(subj_rows, 'SCRNO','EXAMDATE','AVLT_DATE_')
        yrdiff_long = groupLongPivot(subj_rows, 'SCRNO','yrDiff','AVLT_postAV45_')
        data_df = avlt_long.merge(avlt_date,left_index=True,right_index=True)
        data_df = data_df.merge(yrdiff_long,left_index=True,right_index=True)
        # get slope
        data_df.loc[scrno,'AVLT_postAV45_SLOPE'] = groupSlope(subj_rows,'EXAMDATE','TOTS',cutoff_date=av45_date1-timedelta(days=90)) if not isnan(av45_date1) else np.nan
        # get closest
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'TOTS', [av45_date1, av45_date2], day_limit=365/2)
        data_df.loc[scrno,'AVLT_closest_AV45_BL'] = closest_vals[0]
        data_df.loc[scrno,'AVLT_closest_AV45_2'] = closest_vals[1]
        return data_df

    parsed_df = parseSubjectGroups(avlt_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=True)
    return master_df

def syncADNIMasterData(master_df, av45_master_file):
    adni_df = importMaster(av45_master_file, as_df=True)
    adni_subj = list(set(master_df[master_df.index < 6000].index) & set(adni_df.index))
    adni_df = adni_df.loc[adni_subj,:]
    if len(adni_df.index) == 0:
        return master_df

    def extraction_fn(scrno, subj_rows):
        row = subj_rows.iloc[0]
        gender = row['Gender']
        if gender == 1:
            gender = 'M'
        elif gender == 2:
            gender = 'F'
        else:
            raise Exception("Unknown Gender: %s" % gender)
        new_data = {'Sex': gender,
                    'Age': row['Age@AV45'],
                    'Edu': row['Edu.(Yrs)'],
                    'APOE4BIN': row['APOE4_BIN'],
                    'ADAS_TOTSCORE': row['ADAS_AV45_1'],
                    'AVLT_total_6mths_Examdate': row['AVLT_AV45_1'],
                    'CAPS_CURRSCORE': '',
                    'CAPS_LIFETIME_SCORE': '',
                    'WMH_WHITMATHYP': row['WMH_WHITMATHYP.1'],
                    'WMH_percentOfICV': row['WMH_percentOfICV.1'],
                    'GDtotal': row['GD_AV45_1'],
                    'MRIDATE': '',
                    'CSF_abeta': row['CSF_ABETA_closest_AV45_1'],
                    'CSF_tau': row['CSF_TAU_closest_AV45_1'],
                    'CSF_ptau': row['CSF_PTAU_closest_AV45_1'],
                    'SCRNO': scrno
                    }
        return pd.DataFrame([new_data]).set_index('SCRNO')

    parsed_df = parseSubjectGroups(adni_df, extraction_fn)
    master_df = master_df.update(parsed_df, overwrite=True)
    return master_df


def syncDiagData(master_df, diag_file):
    diag_df = importADNIDiagnosis(diag_file, as_df=True)
    timepoints = max(Counter(diag_df.index).values())

    # Add comparison date
    comp_date = datetime(year=2016, month=4, day=1)
    comp_key = 'Diag_%s' % (comp_date.strftime('%b%Y'),)
    diag_df['TimeFromCompDate'] = diag_df['EXAMDATE'].apply(lambda x: abs(x-comp_date).days)

    # header order
    headers = ['Diag_%s' % (_+1,) for _ in range(timepoints)]
    headers += ["Diag_Date_%s" % (_+1,) for _ in range(timepoints)]
    headers += [comp_key, 'Diag_closest_AV45_BL', 'Diag_closest_AV45_2']

    # extraction fn
    def extraction_fn(scrno, subj_rows):
        subj_rows.sort_values(by='EXAMDATE', inplace=True)
        subj_rows['SCRNO'] = scrno
        diag_long = groupLongPivot(subj_rows, 'SCRNO','diag','Diag_')
        date_long = groupLongPivot(subj_rows, 'SCRNO','EXAMDATE','Diag_Date_')
        merge_df = diag_long.merge(date_long,left_index=True,right_index=True)
        comp_diag = subj_rows.sort_values(by='TimeFromCompDate').iloc[0]['diag']
        merge_df.loc[scrno,comp_key] = comp_diag

        av45_date1, av45_date2 = getAV45Dates(scrno, master_df)
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'diag', [av45_date1, av45_date2], day_limit=365/2)
        merge_df.loc[scrno,'Diag_closest_AV45_BL'] = closest_vals[0]
        merge_df.loc[scrno,'Diag_closest_AV45_2'] = closest_vals[1]
        return merge_df

    parsed_df = parseSubjectGroups(diag_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after='Notes')
    return master_df

def syncStudyData(master_df, elig_file):
    '''
    GroupNum: 1=C, 2=PTSD, 3=TBI, 4=PTSD+TBI
    GroupNum_PTSD:  1=C, 2=PTSD or PTSD+TBI, 3=TBI
    GroupNum_TBI: 1=C, 2=PTSD, 3=TBI or PTSD+TBI
    '''

    elig_df = importDODEligibility(elig_file, as_df=True)

    headers = ['PTGroup', 'Study', 'GroupNum', 'GroupNum_TBI', 'GroupNum_PTSD']
    def extraction_fn(scrno, subj_rows):
        cohort = subj_rows.iloc[0]['COHORT']
        new_data = {'SCRNO': scrno,
                    'Study': 1}
        if cohort == 1:
            new_data['PTGroup'] = 'PTSD'
            new_data['GroupNum'] = 2
            new_data['GroupNum_PTSD']= 2
            new_data['GroupNum_TBI'] = 2
        elif cohort == 2:
            new_data['PTGroup'] = 'TBI'
            new_data['GroupNum'] = 3
            new_data['GroupNum_PTSD']= 3
            new_data['GroupNum_TBI'] = 3
        elif cohort == 3:
            new_data['PTGroup'] = 'C'
            new_data['GroupNum'] = 1
            new_data['GroupNum_PTSD']= 1
            new_data['GroupNum_TBI'] = 1
        elif cohort == 4:
            new_data['PTGroup'] = 'TBI+PTSD'
            new_data['GroupNum'] = 4
            new_data['GroupNum_PTSD']= 2
            new_data['GroupNum_TBI'] = 3
        else:
            raise Exception("Unknown cohort: %s" % cohort)
        return pd.DataFrame([new_data]).set_index('SCRNO')

    parsed_df = parseSubjectGroups(elig_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after='Notes')
    adni_subj = master_df[master_df.index < 6000].index
    master_df.loc[adni_subj,'PTGroup'] = 'ADNI C'
    master_df.loc[adni_subj,'Study'] = 0
    master_df.loc[adni_subj,'GroupNum'] = 5
    master_df.loc[adni_subj,'GroupNum_PTSD'] = 5
    master_df.loc[adni_subj,'GroupNum_TBI'] = 5
    return master_df

def syncAV1451RoussetData(master_df, av1451_pvc_file):
    timepoints = ['BL','Scan2']
    av1451_df, _ = importRoussetCSV(av1451_pvc_file, as_df=True)

    av1451_df['BRAAKALL'] = df_mean(av1451_df, ['BRAAK1_SIZE','BRAAK2_SIZE','BRAAK3_SIZE','BRAAK4_SIZE','BRAAK5_SIZE','BRAAK6_SIZE'], ['BRAAK1','BRAAK2','BRAAK3','BRAAK4','BRAAK5','BRAAK6'])
    column_translations = {'BRAAK1': 'AV1451_PVC_Braak1_%s_%s',
                           'BRAAK2': 'AV1451_PVC_Braak2_%s_%s',
                           'BRAAK3': 'AV1451_PVC_Braak3_%s_%s',
                           'BRAAK4': 'AV1451_PVC_Braak4_%s_%s',
                           'BRAAK5': 'AV1451_PVC_Braak5_%s_%s',
                           'BRAAK6': 'AV1451_PVC_Braak6_%s_%s',
                           'BRAAK12': 'AV1451_PVC_Braak12_%s_%s',
                           'BRAAK34': 'AV1451_PVC_Braak34_%s_%s',
                           'BRAAK56': 'AV1451_PVC_Braak56_%s_%s',
                           'BRAAKALL': 'AV1451_PVC_BraakAll_%s_%s'}
    column_order = ['BRAAK1','BRAAK2','BRAAK3','BRAAK4','BRAAK5','BRAAK6',
                    'BRAAK12','BRAAK34','BRAAK56','BRAAKALL']
    ref_region_translations = {'CEREBGM': 'CerebGray',
                               'WHOLECEREB': 'WholeCereb'}
    to_merge = []
    for tp in timepoints:
        cur_df = av1451_df[av1451_df['TP'] == tp]
        roi_df = cur_df[column_order]
        for ref in ref_region_translations:
            suvr_df = roi_df.divide(cur_df[ref],axis=0)
            cur_translations = {k: v % (ref_region_translations[ref],tp) for k,v in column_translations.iteritems()}
            suvr_df.rename(columns=cur_translations,inplace=True)
            to_merge.append(suvr_df)

    merged = to_merge[0]
    for df in to_merge[1:]:
        merged = merged.merge(df,left_index=True,right_index=True,how='outer')

    after = [_ for _ in master_df.columns if _.startswith('AV45') and 'PVC' not in _][-1]
    headers = list(merged.columns)
    master_df = updateDataFrame(master_df, merged, headers=headers, after=after, restrict=True)
    return master_df

def syncAV1451MaxData(master_df, av1451_max_file):
    av1451_df = importAV1451(av1451_file, as_df=True)
    timepoints = max(Counter(av1451_df.index).values())

    # create header order
    row_indices = range(1,timepoints+1)
    headers = []
    headers += ['AV1451_%s_Braak12_Max' % i for i in row_indices]
    headers += ['AV1451_%s_Braak34_Max' % i for i in row_indices]
    headers += ['AV1451_%s_Braak56_Max' % i for i in row_indices]

    # create extraction function
    def extraction_fn(scrno, subj_rows):
        new_data = {'SCRNO': scrno}
        subj_rows.sort_values(by='EXAMDATE', inplace=True)

        for i, (idx, row) in enumerate(subj_rows.iterrows()):
            row_i = i + 1
            braak12 = max(row['BRAAK1'],row['BRAAK2'])
            braak34 = max(row['BRAAK3'],row['BRAAK4'])
            braak56 = max(row['BRAAK5'],row['BRAAK6'])
            new_data['AV1451_%s_Braak12_Max' % row_i] = braak12
            new_data['AV1451_%s_Braak34_Max' % row_i] = braak34
            new_data['AV1451_%s_Braak56_Max' % row_i] = braak56

        return pd.DataFrame([new_data]).set_index('SCRNO')

    parsed_df = parseSubjectGroups(av1451_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=False)
    return master_df


def syncAV1451Data(master_df, av1451_file):
    av1451_df = importAV1451(av1451_file, as_df=True)
    timepoints = max(Counter(av1451_df.index).values())

    # create header order
    row_indices = range(1,timepoints+1)
    headers = []
    headers += ['AV1451_%s_DATE' % i for i in row_indices]
    headers += ['AV1451_%s_Braak12_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_Braak34_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_Braak56_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_LEFT_CAUDATE_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_RIGHT_CAUDATE_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_LEFT_CHOROID_PLEXUS_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_RIGHT_CHOROID_PLEXUS_CerebGray' % i for i in row_indices]
    headers += ['AV1451_%s_CEREBELLUMGREYMATTER' % i for i in row_indices]
    headers += ['AV1451_%s_BRAIN_STEM' % i for i in row_indices]
    headers += ['AV1451_%s_WHOLECEREBELLUM' % i for i in row_indices]

    # create extraction function
    def extraction_fn(scrno, subj_rows):
        new_data = {'SCRNO': scrno}
        subj_rows.sort_values(by='EXAMDATE', inplace=True)

        for i, (idx, row) in enumerate(subj_rows.iterrows()):
            row_i = i + 1

            new_data['AV1451_%s_DATE' % row_i] = row['EXAMDATE']
            new_data['AV1451_%s_BRAIN_STEM' % row_i] = float(row['BRAIN_STEM'])

            # get whole cerebellum
            cerebellum_parts = ['RIGHT_CEREBELLUM_CORTEX', 'RIGHT_CEREBELLUM_WHITE_MATTER', 'LEFT_CEREBELLUM_CORTEX', 'LEFT_CEREBELLUM_WHITE_MATTER']
            weights_values = [(float(row["%s_SIZE" % _]), float(row[_])) for _ in cerebellum_parts]
            wcereb = weightedMean(weights_values)
            new_data['AV1451_%s_WHOLECEREBELLUM' % row_i] = wcereb

            # get cerebellar gray
            cerebg_parts = ['RIGHT_CEREBELLUM_CORTEX', 'LEFT_CEREBELLUM_CORTEX']
            weights_values = [(float(row["%s_SIZE" % _]), float(row[_])) for _ in cerebg_parts]
            cerebg = weightedMean(weights_values)
            new_data['AV1451_%s_CEREBELLUMGREYMATTER' % row_i] = cerebg

            braak1 = (row['BRAAK1_SIZE'],row['BRAAK1'])
            braak2 = (row['BRAAK2_SIZE'],row['BRAAK2'])
            braak3 = (row['BRAAK3_SIZE'],row['BRAAK3'])
            braak4 = (row['BRAAK4_SIZE'],row['BRAAK4'])
            braak5 = (row['BRAAK5_SIZE'],row['BRAAK5'])
            braak6 = (row['BRAAK6_SIZE'],row['BRAAK6'])
            braak12 = weightedMean([braak1,braak2]) / cerebg
            braak34 = weightedMean([braak3,braak4]) / cerebg
            braak56 = weightedMean([braak5,braak6]) / cerebg

            # region suvrs
            new_data['AV1451_%s_Braak12_CerebGray' % row_i] = braak12
            new_data['AV1451_%s_Braak34_CerebGray' % row_i] = braak34
            new_data['AV1451_%s_Braak56_CerebGray' % row_i] = braak56
            new_data['AV1451_%s_LEFT_CAUDATE_CerebGray' % row_i] = float(row['LEFT_CAUDATE'])/cerebg
            new_data['AV1451_%s_RIGHT_CAUDATE_CerebGray' % row_i] = float(row['RIGHT_CAUDATE'])/cerebg
            new_data['AV1451_%s_LEFT_PUTAMEN_CerebGray' % row_i] = float(row['LEFT_PUTAMEN'])/cerebg
            new_data['AV1451_%s_RIGHT_PUTAMEN_CerebGray' % row_i] = float(row['RIGHT_PUTAMEN'])/cerebg
            new_data['AV1451_%s_LEFT_CHOROID_PLEXUS_CerebGray' % row_i] = float(row['LEFT_CHOROID_PLEXUS'])/cerebg
            new_data['AV1451_%s_RIGHT_CHOROID_PLEXUS_CerebGray' % row_i] = float(row['RIGHT_CHOROID_PLEXUS'])/cerebg

        return pd.DataFrame([new_data]).set_index('SCRNO')

    parsed_df = parseSubjectGroups(av1451_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=False)
    return master_df


def syncAV45Data(master_df, av45_file, diags_df):
    av45_df = importAV45(av45_file, as_df=True)
    timepoints = max(Counter(av45_df.index).values())
    adni_controls_df = set(diags_df[diags_df.diag.str.match('N')].index)
    # create header order
    row_indices = range(1,timepoints+1)
    headers = []
    for r in row_indices:
        # Dates
        headers.append('AV45_%s_EXAMDATE' % r)
        # SUVRs
        headers.append('AV45_%s_Left-Putamen' % r)
        headers.append('AV45_%s_Right-Putamen' % r)
        headers.append('AV45_%s_Left-Caudate' % r)
        headers.append('AV45_%s_Right-Caudate' % r)
        headers.append('AV45_%s_Left-Pallidum' % r)
        headers.append('AV45_%s_Right-Pallidum' % r)
        headers.append('AV45_%s_Left_BG_avg' % r)
        headers.append('AV45_%s_Right_BG_avg' % r)
        headers.append('AV45_%s_BG_avg' % r)
        headers.append('AV45_%s_comp/wcerb' % r)
        headers.append('AV45_%s_wcerb_bin' % r)
        headers.append('AV45_%s_comp/brainstem' % r)
        headers.append('AV45_%s_comp/bigref' % r)
        headers.append('AV45_%s_comp/cerbg' % r)
        headers.append('AV45_%s_comp/wm70' % r)
        # Assymetry values
        headers.append('AV45_%s_frontal_asymmetry_negvalue_means_R<L' % r)
        headers.append('AV45_%s_cingulate_asymmetry_negvalue_means_R<L' % r)
        headers.append('AV45_%s_parietal_asymmetry_negvalue_means_R<L' % r)
        headers.append('AV45_%s_temporal_asymmetry_negvalue_means_R<L' % r)
        #headers.append('AV45_%s_asym_summary_absvals_negvalue_means_R<L' % r)
        headers.append('AV45_%s_frontal_MR_asymmetry' % r)
        headers.append('AV45_%s_cingulate_MR_asymmetry' % r)
        headers.append('AV45_%s_parietal_MR_asymmetry' % r)
        headers.append('AV45_%s_temporal__MR_asymmetry' % r)
        headers.append('AV45_%s_ventrical_MR_asymmetry' % r)
    # parse AV45 results
    parsed_df = parseSubjectGroups(av45_df, parseAV45SubjectRows)
    # filter old subjects
    valid_adni_subj = set(master_df.index) &  adni_controls_df
    old_subj = valid_adni_subj | set(master_df[master_df.index >= 6000].index)
    master_df = master_df.loc[list(old_subj),:]
    # update
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=False)
    return master_df

def parseAV45SubjectRows(scrno, subj_rows):
    # reference dicts
    left_frontal_keys = ['CTX_LH_CAUDALMIDDLEFRONTAL', 'CTX_LH_LATERALORBITOFRONTAL', 'CTX_LH_MEDIALORBITOFRONTAL', 'CTX_LH_PARSOPERCULARIS',
                    'CTX_LH_PARSORBITALIS', 'CTX_LH_PARSTRIANGULARIS', 'CTX_LH_ROSTRALMIDDLEFRONTAL', 'CTX_LH_SUPERIORFRONTAL', 'CTX_LH_FRONTALPOLE']
    right_frontal_keys = ['CTX_RH_CAUDALMIDDLEFRONTAL', 'CTX_RH_LATERALORBITOFRONTAL', 'CTX_RH_MEDIALORBITOFRONTAL', 'CTX_RH_PARSOPERCULARIS',
                     'CTX_RH_PARSORBITALIS', 'CTX_RH_PARSTRIANGULARIS', 'CTX_RH_ROSTRALMIDDLEFRONTAL', 'CTX_RH_SUPERIORFRONTAL', 'CTX_RH_FRONTALPOLE']
    left_cingulate_keys = ['CTX_LH_CAUDALANTERIORCINGULATE', 'CTX_LH_ISTHMUSCINGULATE', 'CTX_LH_POSTERIORCINGULATE', 'CTX_LH_ROSTRALANTERIORCINGULATE']
    right_cingulate_keys = ['CTX_RH_CAUDALANTERIORCINGULATE', 'CTX_RH_ISTHMUSCINGULATE', 'CTX_RH_POSTERIORCINGULATE', 'CTX_RH_ROSTRALANTERIORCINGULATE']
    left_parietal_keys = ['CTX_LH_INFERIORPARIETAL', 'CTX_LH_PRECUNEUS', 'CTX_LH_SUPERIORPARIETAL', 'CTX_LH_SUPRAMARGINAL']
    right_parietal_keys = ['CTX_RH_INFERIORPARIETAL', 'CTX_RH_PRECUNEUS', 'CTX_RH_SUPERIORPARIETAL', 'CTX_RH_SUPRAMARGINAL']
    left_temporal_keys = ['CTX_LH_MIDDLETEMPORAL', 'CTX_LH_SUPERIORTEMPORAL']
    right_temporal_keys = ['CTX_RH_MIDDLETEMPORAL', 'CTX_RH_SUPERIORTEMPORAL']
    left_bg_keys = ['LEFT_PUTAMEN', 'LEFT_CAUDATE', 'LEFT_PALLIDUM']
    right_bg_keys = ['RIGHT_PUTAMEN', 'RIGHT_CAUDATE', 'RIGHT_PALLIDUM']
    left_ventrical_keys = []
    right_ventrical_keys = []

    data = {'SCRNO': int(scrno)}
    subj_rows.sort_values(by='EXAMDATE', inplace=True)

    for i, (idx, point) in enumerate(subj_rows.iterrows()):
        examdate = point['EXAMDATE']
        # extract necessary values
        wcereb = float(point['WHOLECEREBELLUM'])
        cerebg = float(point['CEREBELLUMGREYMATTER'])
        compositeroi = float(point['COMPOSITE'])
        bigref = float(point['COMPOSITE_REF'])
        wm70 = float(point['ERODED_SUBCORTICALWM'])
        brainstem = float(point['BRAIN_STEM'])
        cingulate = float(point['CINGULATE'])
        frontal = float(point['FRONTAL'])
        parietal = float(point['PARIETAL'])
        temporal = float(point['TEMPORAL'])
        leftputamen = float(point['LEFT_PUTAMEN'])
        rightputamen = float(point['RIGHT_PUTAMEN'])
        leftcaudate = float(point['LEFT_CAUDATE'])
        rightcaudate = float(point['RIGHT_CAUDATE'])
        leftpallidum = float(point['LEFT_PALLIDUM'])
        rightpallidum = float(point['RIGHT_PALLIDUM'])
        left_frontal = np.average([float(point[k]) for k in left_frontal_keys], weights=[float(point["%s_SIZE" % k]) for k in left_frontal_keys])
        right_frontal = np.average([float(point[k]) for k in right_frontal_keys], weights=[float(point["%s_SIZE" % k]) for k in right_frontal_keys])
        left_cingulate = np.average([float(point[k]) for k in left_cingulate_keys], weights=[float(point["%s_SIZE" % k]) for k in left_cingulate_keys])
        right_cingulate = np.average([float(point[k]) for k in right_cingulate_keys], weights=[float(point["%s_SIZE" % k]) for k in right_cingulate_keys])
        left_parietal = np.average([float(point[k]) for k in left_parietal_keys], weights=[float(point["%s_SIZE" % k]) for k in left_parietal_keys])
        right_parietal = np.average([float(point[k]) for k in right_parietal_keys], weights=[float(point["%s_SIZE" % k]) for k in right_parietal_keys])
        left_temporal = np.average([float(point[k]) for k in left_temporal_keys], weights=[float(point["%s_SIZE" % k]) for k in left_temporal_keys])
        right_temporal = np.average([float(point[k]) for k in right_temporal_keys], weights=[float(point["%s_SIZE" % k]) for k in right_temporal_keys])
        left_bg = np.average([float(point[k]) for k in left_bg_keys])
        right_bg = np.average([float(point[k]) for k in right_bg_keys])
        left_frontal_size = np.sum([float(point["%s_SIZE" % k]) for k in left_frontal_keys])
        right_frontal_size = np.sum([float(point["%s_SIZE" % k]) for k in right_frontal_keys])
        left_cingulate_size = np.sum([float(point["%s_SIZE" % k]) for k in left_cingulate_keys])
        right_cingulate_size = np.sum([float(point["%s_SIZE" % k]) for k in right_cingulate_keys])
        left_parietal_size = np.sum([float(point["%s_SIZE" % k]) for k in left_parietal_keys])
        right_parietal_size = np.sum([float(point["%s_SIZE" % k]) for k in right_parietal_keys])
        left_temporal_size = np.sum([float(point["%s_SIZE" % k]) for k in left_temporal_keys])
        right_temporal_size = np.sum([float(point["%s_SIZE" % k]) for k in right_temporal_keys])
        left_ventrical_size = np.sum([float(point["%s_SIZE" % k]) for k in left_ventrical_keys])
        right_ventrical_size = np.sum([float(point["%s_SIZE" % k]) for k in right_ventrical_keys])

        # Dates
        data['AV45_%s_EXAMDATE' % (i+1)] = examdate
        # SUVR values
        data['AV45_%s_Left-Putamen' % (i+1)] = leftputamen/wm70
        data['AV45_%s_Right-Putamen' % (i+1)] = rightputamen/wm70
        data['AV45_%s_Left-Caudate' % (i+1)] = leftcaudate/wm70
        data['AV45_%s_Right-Caudate' % (i+1)] = rightcaudate/wm70
        data['AV45_%s_Left-Pallidum' % (i+1)] = leftpallidum/wm70
        data['AV45_%s_Right-Pallidum' % (i+1)] = rightpallidum/wm70
        data['AV45_%s_Left_BG_avg' % (i+1)] = left_bg/wm70
        data['AV45_%s_Right_BG_avg' % (i+1)] = right_bg/wm70
        data['AV45_%s_BG_avg' % (i+1)] = np.mean([left_bg,right_bg])/wm70
        data['AV45_%s_comp/wcerb' % (i+1)] = compositeroi / wcereb
        data['AV45_%s_wcerb_bin' % (i+1)] = 1 if (compositeroi/wcereb) >= 1.11 else 0 # WHATS THE REAL THRESHOLD
        data['AV45_%s_comp/brainstem' % (i+1)] = compositeroi / brainstem
        data['AV45_%s_comp/bigref' % (i+1)] = compositeroi / bigref
        data['AV45_%s_comp/cerbg' % (i+1)] = compositeroi / cerebg
        data['AV45_%s_comp/wm70' % (i+1)] = compositeroi / wm70
        # Asymmetry values
        data['AV45_%s_frontal_asymmetry_negvalue_means_R<L' % (i+1)] = asymIndex(left_frontal, right_frontal)
        data['AV45_%s_cingulate_asymmetry_negvalue_means_R<L' % (i+1)] = asymIndex(left_cingulate, right_cingulate)
        data['AV45_%s_parietal_asymmetry_negvalue_means_R<L' % (i+1)] = asymIndex(left_parietal, right_parietal)
        data['AV45_%s_temporal_asymmetry_negvalue_means_R<L' % (i+1)] = asymIndex(left_temporal, right_temporal)
        left_summary = np.mean([left_frontal, left_cingulate, left_parietal, left_temporal])
        right_summary = np.mean([right_frontal, right_cingulate, right_parietal, right_temporal])
        data['AV45_%s_frontal_MR_asymmetry' % (i+1)] = asymIndex(left_frontal_size, right_frontal_size)
        data['AV45_%s_cingulate_MR_asymmetry' % (i+1)] = asymIndex(left_cingulate_size, right_cingulate_size)
        data['AV45_%s_parietal_MR_asymmetry' % (i+1)] = asymIndex(left_parietal_size, right_parietal_size)
        data['AV45_%s_temporal__MR_asymmetry' % (i+1)] = asymIndex(left_temporal_size, right_temporal_size)
        data['AV45_%s_ventrical_MR_asymmetry' % (i+1)] = asymIndex(left_ventrical_size, right_ventrical_size)

    return pd.DataFrame([data]).set_index('SCRNO')

def getAV1451Dates(scrno, master_df):
    bl_av1451 = av1451_2 = np.nan
    try:
        row = master_df.loc[scrno]
        bl_av1451 = row.get('AV1451_1_DATE')
        av1451_2 = row.get('AV1451_2_DATE')
    except Exception as e:
        pass
    return (bl_av1451, av1451_2)


def getAV45Dates(scrno, master_df):
    bl_av45 = av45_2 = np.nan
    # print master_df.index
    try:
        row = master_df.loc[scrno]
        bl_av45 = row.get('AV45_1_EXAMDATE')
        av45_2 = row.get('AV45_2_EXAMDATE')
    except Exception as e:
        pass
    return (bl_av45, av45_2)


def runPipeline():
    # syncing pipeline
    try:
        master_df = pd.read_csv(master_file)
        master_df.set_index('SCRNO',inplace=True)
        print len(master_df.index)
    except:
        master_df = pd.DataFrame(columns=['Notes'])
        new_lines = []

    # get diagnoses
    diags_df = extractDiagnosesFromMasterData(importMaster(av45_master_file, as_df=True))
    registry = importDODRegistry(registry_file, as_df=True)

    print "\nSYNCING AV45\n"
    master_df = syncAV45Data(master_df, av45_file, diags_df) # adds new patients
    print "\nSYNCING AV1451\n"
    master_df = syncAV1451Data(master_df, av1451_file) # adds new patients
    print "\nSYNCING AV1451 MAX\n"
    master_df = syncAV1451MaxData(master_df, av1451_max_file)
    print "\nSYNCING AV1451 Rousset\n"
    master_df = syncAV1451RoussetData(master_df, av1451_pvc_file)
    print "\nSYNCING DIAG\n"
    master_df = syncDiagData(master_df, diag_file)
    print "\nSYNCING DTI\n"
    master_df = syncDTI(master_df, dti_file, registry)
    print "\nSYNCING TBIAGE\n"
    master_df = syncTBIAGE(master_df, tbiage_file, registry)
    print "\nSYNCING SLEEP QUALITY\n"
    master_df = syncSleep(master_df, sleep_file, registry)
    print "\nSYNCING ANTIDEP\n"
    master_df = syncAntidepData(master_df, backmeds_file, registry)
    print "\nSYNCING APOE\n"
    master_df = syncAPOEData(master_df, apoe_file)
    print "\nSYNCING ASI\n"
    master_df = syncASIData(master_df, asi_file, registry)
    print "\nSYNCING CDR\n"
    master_df = syncCDRData(master_df, cdr_file, registry)
    print "\nSYNCING FAQ\n"
    master_df = syncFAQData(master_df, faq_file, registry)
    print "\nSYNCING MMSE\n"
    master_df = syncMMSEData(master_df, mmse_file, registry)
    print "\nSYNCING ADAS\n"
    master_df = syncADASData(master_df, adas_file, registry)
    print "\nSYNCING AVLT\n"
    master_df = syncAVLTData(master_df, avlt_file, registry)
    print "\nSYNCING DEMOG\n"
    master_df = syncDemogData(master_df, demog_file)
    print "\nSYNCING CAPS\n"
    master_df = syncCAPSData(master_df, caps_curr_file, caps_lifetime_file)
    print "\nSYNCING WMH\n"
    master_df = syncWMHData(master_df, wmh_file)
    print "\nSYNCING GD\n"
    master_df = syncGDData(master_df, gd_file, registry)
    print "\nSYNCING MRI\n"
    master_df = syncMRIData(master_df, mri_file)
    print "\nSYNCING STUDY\n"
    master_df = syncStudyData(master_df, elig_file)
    print "\nSYNCING CSF\n"
    master_df = syncCSFData(master_df, csf_file, registry)
    print "\nSYNCING ADNI MASTER\n"
    master_df = syncADNIMasterData(master_df, av45_master_file)
    print "\nDUMPING CSV\n"
    dumpDFtoCSV(master_df,output_file,decimal_places=3)

if __name__ == '__main__':
    now = datetime.now()
    pd.options.mode.chained_assignment = None

    # IO files
    master_file = "" # not used
    output_file = "../DOD_DATA_%s.csv" % (now.strftime("%m_%d_%y"))

    # ADNI master file
    av45_master_file = '../FDG_AV45_COGdata/FDG_AV45_COGdata_09_19_16.csv'
    # AV45 file
    av45_file = '../output/10-14-2016/UCBERKELEYAV45_DOD_10-14-2016_regular_nontp.csv'
    # AV1451 file
    av1451_file = '../output/10-14-2016/UCBERKELEYAV1451_DOD_10-14-2016_regular_tp.csv'
    # AV1451 Max file
    av1451_max_file = '../output/10-14-2016/UCBERKELEYAV1451_DOD_MAX_10-14-2016_regular_tp.csv'
    # AV1451 PVC file
    av1451_pvc_file = '../pvc/pvc_dod_av1451/tauskullregions_output.csv'

    # Registry file
    registry_file = "../docs/DOD/REGISTRY.csv"
    # ASI file
    asi_file = "../docs/DOD/ASI.csv"
    # AVLT file
    avlt_file = "../docs/DOD/NEUROBAT.csv"
    # ADAS file
    adas_file = "../docs/DOD/ADAS.csv"
    # FAQ file
    faq_file = "../docs/DOD/FAQ.csv"
    # MMSE file
    mmse_file = "../docs/DOD/MMSE.csv"
    # CDR file
    cdr_file = "../docs/DOD/CDR.csv"
    # Demog file
    demog_file = "../docs/DOD/PTDEMOG.csv"
    # CSF file
    csf_file = "../docs/DOD/UPENNBIOMK.csv"
    # APOE file
    apoe_file = "../docs/DOD/APOERES.csv"
    # CAPS files
    caps_curr_file = "../docs/DOD/CAPSCURR.csv"
    caps_lifetime_file = "../docs/DOD/CAPSLIFE.csv"
    # WMH file
    wmh_file = "../docs/DOD/WMH_VOLUMES.csv"
    # GD file
    gd_file = "../docs/DOD/GDSCALE.csv"
    # MRi file
    mri_file = "../docs/DOD/MRIMETA.csv"
    # Elig file
    elig_file = "../docs/DOD/VAELG.csv"
    # diag file
    diag_file = "../docs/DOD/DXSUM.csv"
    # antidep file
    backmeds_file = '../docs/DOD/BACKMEDS.csv'
    # DTI file
    dti_file = '../docs/DOD/THOMPSONDTI.csv'
    # TBI age file
    tbiage_file = '../docs/DOD/RECTBIINJ.csv'
    # Sleep quality file
    sleep_file = '../docs/DOD/PSQI.csv'

    runPipeline()
