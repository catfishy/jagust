'''

Pull in the master FDG_AV45_COG DATA

remove null char:
tr -d '\000' < file1 > file2

'''
import re
from collections import defaultdict
from datetime import datetime, timedelta
from scipy import optimize
from scipy.stats import linregress
import numpy as np
import itertools

from utils import *



RID_ADDONS = [78,93,108,458,514,691,724,739,821,853,916,1080,1271]

######### SYNCING INFRASTRUCTURE (MOST FUNCTIONS DEPEND ON THESE) ##################

def getAV45Dates(rid, master_df):
    bl_av45 = av45_2 = av45_3 = np.nan
    try:
        row = master_df.loc[rid]
        bl_av45 = row.get('AV45_Date')
        av45_2 = row.get('AV45_2_Date')
        av45_3 = row.get('AV45_3_Date')
    except Exception as e:
        pass
    return (bl_av45, av45_2, av45_3)

def getAV1451Dates(rid, master_df):
    bl_av1451 = av1451_2 = av1451_3 = np.nan
    try:
        row = master_df.loc[rid]
        bl_av1451 = row.get('AV1451_Date')
        av1451_2 = row.get('AV1451_2_Date')
        av1451_3 = row.get('AV1451_3_Date')
    except Exception as e:
        pass
    return (bl_av1451, av1451_2, av1451_3)

def addCategories(output_file):
    colcats_df = pd.read_csv('column_categories.csv')
    colcats_df.set_index('COLUMN',inplace=True)
    df = pd.read_csv(output_file, low_memory=False)
    cats = {k:'Unknown' for k in df.columns}
    for column, category in colcats_df['CATEGORY'].iteritems():
        if '*' in column:
            new_data = {k: category for k in cats.keys() if re.search(column,k.strip())}
            cats.update(new_data)
        elif column in cats:
            cats[column] = category

    # print unknowns
    for k,v in cats.iteritems():
        if v == 'Unknown':
            print k

    df.columns = [[cats[_] for _ in df.columns],df.columns]
    df.to_csv(output_file, index=False)

def mergeSlopes(output_file):
    # merge 2/3 point slope columns into one column (replacing 2pt slopes where 3pt slopes are available)
    # add merged column in after 3 point column
    df = pd.read_csv(output_file, low_memory=False)
    twopoint_cols = [i for i in df.columns if '2pts' in i or '2points' in i]
    threepoint_cols = [i for i in df.columns if '3pts' in i or '3points' in i]
    for two_col in twopoint_cols:
        # find three point col
        new_col = two_col.replace('_2pts','').replace('_2points','')
        if new_col in df.columns:
            continue
        three_col = two_col.replace('2pts','3pts').replace('2points', '3points')
        if three_col not in threepoint_cols:
            print "Could not find 3 point slope column for %s" % (three_col,)
            continue
        three_col_ind = df.columns.get_loc(three_col)
        df.insert(three_col_ind+1, new_col, df.loc[:,two_col])
        df[new_col].update(df[three_col])
    df.to_csv(output_file, index=False)

def manualAddOns(master_df, rid_list):
    to_add = [int(_) for _ in rid_list if int(_) not in master_df.index]
    for rid in to_add:
        master_df.loc[rid] = None
    master_df.sort_index(inplace=True)
    return master_df


######### END SYNCING INFRASTRUCTURE ###############################################


def syncAV1451RoussetResults(master_df, rousset_csv):
    av1451_df, _ = importRoussetCSV(rousset_csv, translate_threshold=None, as_df=True)
    av1451_df['BRAAK12'] = df_mean(av1451_df, ['BRAAK1_SIZE','BRAAK2_SIZE'], ['BRAAK1','BRAAK2']) / av1451_df['CEREBGM']
    av1451_df['BRAAK34'] = df_mean(av1451_df, ['BRAAK3_SIZE','BRAAK4_SIZE'], ['BRAAK3','BRAAK4']) / av1451_df['CEREBGM']
    av1451_df['BRAAK56'] = df_mean(av1451_df, ['BRAAK5_SIZE','BRAAK6_SIZE'], ['BRAAK5','BRAAK6']) / av1451_df['CEREBGM']

    # make pivots
    index_name = av1451_df.index.name
    av1451_df.reset_index(inplace=True)
    braak12_df = av1451_df.pivot(index_name,'TP','BRAAK12')
    braak34_df = av1451_df.pivot(index_name,'TP','BRAAK34')
    braak56_df = av1451_df.pivot(index_name,'TP','BRAAK56')

    # rename columns
    braak12_df.columns = ['AV1451_PVC_Braak12_CerebGray_%s' % _ for _ in braak12_df.columns]
    braak34_df.columns = ['AV1451_PVC_Braak34_CerebGray_%s' % _ for _ in braak34_df.columns]
    braak56_df.columns = ['AV1451_PVC_Braak56_CerebGray_%s' % _ for _ in braak56_df.columns]

    # merge together
    parsed_df = braak12_df.merge(braak34_df,left_index=True,right_index=True)
    parsed_df = parsed_df.merge(braak56_df,left_index=True,right_index=True)

    valid_timepoints = ['BL', 'Scan2', 'Scan3']
    headers = ['AV1451_PVC_Braak12_CerebGray_%s' % tp for tp in valid_timepoints]
    headers += ['AV1451_PVC_Braak34_CerebGray_%s' % tp for tp in valid_timepoints]
    headers += ['AV1451_PVC_Braak56_CerebGray_%s' % tp for tp in valid_timepoints]
    after = [_ for _ in master_df.columns if _.startswith('AV45') and 'PVC' not in _][-1]

    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=after, restrict=True)
    return master_df

def syncAV45RoussetResults(master_df, av45_rousset_csv):
    av45_df, threshold = importRoussetCSV(av45_rousset_csv, translate_threshold=1.11, as_df=True)
    av45_df = av45_df[['TP','COMPOSITE','WHOLECEREB','BIGREF']]
    av45_df['WCEREB_SUVR'] = av45_df['COMPOSITE']/av45_df['WHOLECEREB']
    av45_df['BIGREF_SUVR'] = av45_df['COMPOSITE']/av45_df['BIGREF']
    index_name = av45_df.index.name
    av45_df.reset_index(inplace=True)

    # get pivoted dataframes
    wcereb_suvr_df = av45_df.pivot(index_name,'TP','WCEREB_SUVR')
    wcereb_suvr_pos_df = (wcereb_suvr_df>=threshold).astype(int)
    wcereb_suvr_pos_df[pd.isnull(wcereb_suvr_df)] = np.nan
    bigref_suvr_df = av45_df.pivot(index_name,'TP','BIGREF_SUVR')

    # get scan dates
    dates_df = wcereb_suvr_df.apply(lambda x: pd.Series(dict(enumerate(getAV45Dates(x.name, master_df)))), axis=1)

    # get summary/wcereb slopes
    wcereb_suvr_df = wcereb_suvr_df.merge(dates_df,left_index=True,right_index=True)
    wcereb_two_point_slope = df_slope(wcereb_suvr_df, [0,1], ['BL','Scan2'],take_diff=True, exact=True)
    wcereb_three_point_slope = df_slope(wcereb_suvr_df, [0,1,2], ['BL','Scan2','Scan3'],take_diff=True, exact=True)
    wcereb_suvr_df.drop(dates_df.columns,axis=1,inplace=True)

    # get summary/bigref slopes
    bigref_suvr_df = bigref_suvr_df.merge(dates_df,left_index=True,right_index=True)
    bigref_two_point_slope = df_slope(bigref_suvr_df, [0,1], ['BL','Scan2'],take_diff=True, exact=True)
    bigref_three_point_slope = df_slope(bigref_suvr_df, [0,1,2], ['BL','Scan2','Scan3'],take_diff=True, exact=True)
    bigref_suvr_df.drop(dates_df.columns,axis=1,inplace=True)

    # rename columns
    wcereb_suvr_df.columns = ['AV45_PVC_CorticalSummary_WholeCereb_%s' % _ for _ in wcereb_suvr_df.columns]
    wcereb_suvr_pos_df.columns = ['AV45_PVC_CorticalSummary_WholeCereb_%s_%s' % (threshold,_) for _ in wcereb_suvr_pos_df.columns]
    bigref_suvr_df.columns = ['AV45_PVC_CorticalSummary_BigRef_%s' % _ for _ in bigref_suvr_df.columns]

    # add on slopes
    wcereb_suvr_df['AV45_PVC_CorticalSummary_WholeCereb_slope_2points'] = wcereb_two_point_slope
    wcereb_suvr_df['AV45_PVC_CorticalSummary_WholeCereb_slope_3points'] = wcereb_three_point_slope
    bigref_suvr_df['AV45_PVC_CorticalSummary_BigRef_slope_2points'] = bigref_two_point_slope
    bigref_suvr_df['AV45_PVC_CorticalSummary_BigRef_slope_3points'] = bigref_three_point_slope

    # merge everything together
    parsed_df = wcereb_suvr_df.merge(wcereb_suvr_pos_df,left_index=True,right_index=True)
    parsed_df = parsed_df.merge(bigref_suvr_df,left_index=True,right_index=True)

    # create headers
    valid_timepoints = ['BL', 'Scan2', 'Scan3']
    headers = ['AV45_PVC_CorticalSummary_WholeCereb_%s' % tp for tp in valid_timepoints]
    headers += ['AV45_PVC_CorticalSummary_WholeCereb_%s_%s' % (threshold,tp) for tp in valid_timepoints]
    headers += ['AV45_PVC_CorticalSummary_WholeCereb_slope_2points', 'AV45_PVC_CorticalSummary_WholeCereb_slope_3points']
    headers += ['AV45_PVC_CorticalSummary_BigRef_%s' % tp for tp in valid_timepoints]
    headers += ['AV45_PVC_CorticalSummary_BigRef_slope_2points', 'AV45_PVC_CorticalSummary_BigRef_slope_3points']
    after = [_ for _ in master_df.columns if _.startswith('AV45') and 'PVC' not in _][-1]

    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=after, restrict=True)
    return master_df

def syncFAQData(master_df, faq_file, registry):
    faq_df = importFAQ(faq_file, registry, as_df=True)
    tmpts = max(Counter(faq_df.index).values())

    headers = ['FAQTOTAL.%s' % (i+1) for i in range(tmpts)]
    headers += ['FAQTOTAL_timePostAV45.%s' % (i+1) for i in range(tmpts)]
    headers += ['FAQTOTAL_AV45_6MTHS', 'FAQTOTAL_AV45_DATE',
                'FAQTOTAL_AV45_2_6MTHS', 'FAQTOTAL_AV45_2_DATE',
                'FAQTOTAL_AV45_3_6MTHS', 'FAQTOTAL_AV45_3_DATE',
                'FAQTOTAL_slope']

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)

        # get longitudinal measurements
        subj_rows['RID'] = rid
        faq_long = groupLongPivot(subj_rows, 'RID','FAQTOTAL','FAQTOTAL.')
        date_long = groupLongPivot(subj_rows, 'RID','EXAMDATE','FAQTOTAL_timePostAV45.')
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)

        all_df = pd.concat((faq_long,date_long),axis=1)
        slope = df_slope(all_df,
                         date_long.columns,
                         faq_long.columns,
                         take_diff=False,
                         exact=False)
        all_df['FAQTOTAL_slope'] = slope[rid]

        # closest to AV45 scans
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'FAQTOTAL', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        all_df['FAQTOTAL_AV45_6MTHS'] = closest_vals[0]
        all_df['FAQTOTAL_AV45_DATE'] = closest_dates[0]
        all_df['FAQTOTAL_AV45_2_6MTHS'] = closest_vals[1]
        all_df['FAQTOTAL_AV45_2_DATE'] = closest_dates[1]
        all_df['FAQTOTAL_AV45_3_6MTHS'] = closest_vals[2]
        all_df['FAQTOTAL_AV45_3_DATE'] = closest_dates[2]

        all_df['FAQTOTAL_AV45_DATE'] = pd.to_datetime(all_df['FAQTOTAL_AV45_DATE'])
        all_df['FAQTOTAL_AV45_2_DATE'] = pd.to_datetime(all_df['FAQTOTAL_AV45_2_DATE'])
        all_df['FAQTOTAL_AV45_3_DATE'] = pd.to_datetime(all_df['FAQTOTAL_AV45_3_DATE'])

        return all_df


    parsed_df = parseSubjectGroups(faq_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df


def syncNPIData(master_df, npi_file):
    npi_df = importNPI(npi_file,as_df=True)
    tmpts = max(Counter(npi_df.index).values())

    headers = ['NPITOTAL.%s' % (i+1) for i in range(tmpts)]
    headers += ['NPITOTAL_timePostAV45.%s' % (i+1) for i in range(tmpts)]
    headers += ['NPITOTAL_AV45_6MTHS', 'NPITOTAL_AV45_DATE',
                'NPITOTAL_AV45_2_6MTHS', 'NPITOTAL_AV45_2_DATE',
                'NPITOTAL_AV45_3_6MTHS', 'NPITOTAL_AV45_3_DATE',
                'NPITOTAL_slope']

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)

        # get longitudinal measurements
        subj_rows['RID'] = rid
        npi_long = groupLongPivot(subj_rows, 'RID','NPITOTAL','NPITOTAL.')
        date_long = groupLongPivot(subj_rows, 'RID','EXAMDATE','NPITOTAL_timePostAV45.')
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)

        all_df = pd.concat((npi_long,date_long),axis=1)
        slope = df_slope(all_df,
                         date_long.columns,
                         npi_long.columns,
                         take_diff=False,
                         exact=False)
        all_df['NPITOTAL_slope'] = slope[rid]

        # closest to AV45 scans
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'NPITOTAL', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        all_df['NPITOTAL_AV45_6MTHS'] = closest_vals[0]
        all_df['NPITOTAL_AV45_DATE'] = closest_dates[0]
        all_df['NPITOTAL_AV45_2_6MTHS'] = closest_vals[1]
        all_df['NPITOTAL_AV45_2_DATE'] = closest_dates[1]
        all_df['NPITOTAL_AV45_3_6MTHS'] = closest_vals[2]
        all_df['NPITOTAL_AV45_3_DATE'] = closest_dates[2]

        all_df['NPITOTAL_AV45_DATE'] = pd.to_datetime(all_df['NPITOTAL_AV45_DATE'])
        all_df['NPITOTAL_AV45_2_DATE'] = pd.to_datetime(all_df['NPITOTAL_AV45_2_DATE'])
        all_df['NPITOTAL_AV45_3_DATE'] = pd.to_datetime(all_df['NPITOTAL_AV45_3_DATE'])

        return all_df

    parsed_df = parseSubjectGroups(npi_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df


def syncMHISTData(master_df, mhist_file):
    mhist_df = importMedicalHistory(mhist_file,as_df=True)

    mhist_df = mhist_df[['smoking','diabetes','hyperlipidemia']].astype(int)
    mhist_df.columns = [_.upper() for _ in mhist_df.columns]

    master_df = updateDataFrame(master_df, mhist_df, headers=mhist_df.columns, after='Handedness')
    return master_df

def syncCDRData(master_df, cdr_file, registry):
    cdr_df = importCDR(cdr_file, registry, as_df=True)
    tmpts = max(Counter(cdr_df.index).values())

    headers = ['CDR_GLOBAL_AV45_1','CDR_GLOBAL_AV45_2','CDR_GLOBAL_AV45_3',
               'CDR_MEMORY_AV45_1','CDR_MEMORY_AV45_2','CDR_MEMORY_AV45_3',
               'CDR_ORIENT_AV45_1','CDR_ORIENT_AV45_2','CDR_ORIENT_AV45_3',
               'CDR_JUDGE_AV45_1','CDR_JUDGE_AV45_2','CDR_JUDGE_AV45_3',
               'CDR_COMMUN_AV45_1','CDR_COMMUN_AV45_2','CDR_COMMUN_AV45_3',
               'CDR_HOME_AV45_1','CDR_HOME_AV45_2','CDR_HOME_AV45_3',
               'CDR_CARE_AV45_1','CDR_CARE_AV45_2','CDR_CARE_AV45_3',
               'CDR_GLOBAL_AV1451_1','CDR_GLOBAL_AV1451_2','CDR_GLOBAL_AV1451_3',
               'CDR_MEMORY_AV1451_1','CDR_MEMORY_AV1451_2','CDR_MEMORY_AV1451_3',
               'CDR_ORIENT_AV1451_1','CDR_ORIENT_AV1451_2','CDR_ORIENT_AV1451_3',
               'CDR_JUDGE_AV1451_1','CDR_JUDGE_AV1451_2','CDR_JUDGE_AV1451_3',
               'CDR_COMMUN_AV1451_1','CDR_COMMUN_AV1451_2','CDR_COMMUN_AV1451_3',
               'CDR_HOME_AV1451_1','CDR_HOME_AV1451_2','CDR_HOME_AV1451_3',
               'CDR_CARE_AV1451_1','CDR_CARE_AV1451_2','CDR_CARE_AV1451_3']

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)
        data = {'RID': rid}

        # get closest av45
        cdglobal = groupClosest(subj_rows, 'EXAMDATE', 'CDGLOBAL', [av45_date1,av45_date2,av45_date3],day_limit=365)
        cdcare = groupClosest(subj_rows, 'EXAMDATE', 'CDCARE', [av45_date1,av45_date2,av45_date3],day_limit=365)
        cdhome = groupClosest(subj_rows, 'EXAMDATE', 'CDHOME', [av45_date1,av45_date2,av45_date3],day_limit=365)
        cdcommun = groupClosest(subj_rows, 'EXAMDATE', 'CDCOMMUN', [av45_date1,av45_date2,av45_date3],day_limit=365)
        cdjudge = groupClosest(subj_rows, 'EXAMDATE', 'CDJUDGE', [av45_date1,av45_date2,av45_date3],day_limit=365)
        cdorient = groupClosest(subj_rows, 'EXAMDATE', 'CDORIENT', [av45_date1,av45_date2,av45_date3],day_limit=365)
        cdmemory = groupClosest(subj_rows, 'EXAMDATE', 'CDMEMORY', [av45_date1,av45_date2,av45_date3],day_limit=365)
        data['CDR_GLOBAL_AV45_1'],data['CDR_GLOBAL_AV45_2'],data['CDR_GLOBAL_AV45_3'] = tuple(cdglobal)
        data['CDR_CARE_AV45_1'],data['CDR_CARE_AV45_2'],data['CDR_CARE_AV45_3'] = tuple(cdcare)
        data['CDR_HOME_AV45_1'],data['CDR_HOME_AV45_2'],data['CDR_HOME_AV45_3'] = tuple(cdhome)
        data['CDR_COMMUN_AV45_1'],data['CDR_COMMUN_AV45_2'],data['CDR_COMMUN_AV45_3'] = tuple(cdcommun)
        data['CDR_JUDGE_AV45_1'],data['CDR_JUDGE_AV45_2'],data['CDR_JUDGE_AV45_3'] = tuple(cdjudge)
        data['CDR_ORIENT_AV45_1'],data['CDR_ORIENT_AV45_2'],data['CDR_ORIENT_AV45_3'] = tuple(cdorient)
        data['CDR_MEMORY_AV45_1'],data['CDR_MEMORY_AV45_2'],data['CDR_MEMORY_AV45_3'] = tuple(cdmemory)

        # get closest av1451
        cdglobal = groupClosest(subj_rows, 'EXAMDATE', 'CDGLOBAL', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        cdcare = groupClosest(subj_rows, 'EXAMDATE', 'CDCARE', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        cdhome = groupClosest(subj_rows, 'EXAMDATE', 'CDHOME', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        cdcommun = groupClosest(subj_rows, 'EXAMDATE', 'CDCOMMUN', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        cdjudge = groupClosest(subj_rows, 'EXAMDATE', 'CDJUDGE', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        cdorient = groupClosest(subj_rows, 'EXAMDATE', 'CDORIENT', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        cdmemory = groupClosest(subj_rows, 'EXAMDATE', 'CDMEMORY', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        data['CDR_GLOBAL_AV45_1'],data['CDR_GLOBAL_AV45_2'],data['CDR_GLOBAL_AV45_3'] = tuple(cdglobal)
        data['CDR_CARE_AV45_1'],data['CDR_CARE_AV45_2'],data['CDR_CARE_AV45_3'] = tuple(cdcare)
        data['CDR_HOME_AV45_1'],data['CDR_HOME_AV45_2'],data['CDR_HOME_AV45_3'] = tuple(cdhome)
        data['CDR_COMMUN_AV45_1'],data['CDR_COMMUN_AV45_2'],data['CDR_COMMUN_AV45_3'] = tuple(cdcommun)
        data['CDR_JUDGE_AV45_1'],data['CDR_JUDGE_AV45_2'],data['CDR_JUDGE_AV45_3'] = tuple(cdjudge)
        data['CDR_ORIENT_AV45_1'],data['CDR_ORIENT_AV45_2'],data['CDR_ORIENT_AV45_3'] = tuple(cdorient)
        data['CDR_MEMORY_AV45_1'],data['CDR_MEMORY_AV45_2'],data['CDR_MEMORY_AV45_3'] = tuple(cdmemory)

        return pd.DataFrame([data]).set_index('RID')

    parsed_df = parseSubjectGroups(cdr_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df

def syncGDData(master_df, gd_file, registry):
    gd_df = importGD(gd_file, registry=registry, as_df=True)
    tmpts = max(Counter(gd_df.index).values())

    headers = ['GD_TOTAL.%s' % (i+1) for i in range(tmpts)]
    headers += ['GD_timePostAV45.%s' % (i+1) for i in range(tmpts)]
    headers += ['GD_AV45_1', 'GD_AV45_1_DATE',
                'GD_AV45_2', 'GD_AV45_2_DATE',
                'GD_AV45_3', 'GD_AV45_3_DATE',
                'GD_AV1451_1', 'GD_AV1451_1_DATE',
                'GD_AV1451_2', 'GD_AV1451_2_DATE',
                'GD_AV1451_3', 'GD_AV1451_3_DATE',
                'GD_slope']

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)

        # Get long measurements
        subj_rows['RID'] = rid
        gd_long = groupLongPivot(subj_rows,'RID','GDTOTAL','GD_TOTAL.')
        date_long = groupLongPivot(subj_rows,'RID','EXAMDATE','GD_timePostAV45.')
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)
        all_df = pd.concat((gd_long,date_long),axis=1)

        gd_slope = df_slope(all_df, date_long.columns, gd_long.columns, take_diff=False, exact=False)
        all_df['GD_slope'] = gd_slope[rid]

        # Get closest av45 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'GDTOTAL', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        all_df['GD_AV45_1'] = closest_vals[0]
        all_df['GD_AV45_1_DATE'] = closest_dates[0]
        all_df['GD_AV45_2'] = closest_vals[1]
        all_df['GD_AV45_2_DATE'] = closest_dates[1]
        all_df['GD_AV45_3'] = closest_vals[2]
        all_df['GD_AV45_3_DATE'] = closest_dates[2]
        all_df['GD_AV45_1_DATE'] = pd.to_datetime(all_df['GD_AV45_1_DATE'])
        all_df['GD_AV45_2_DATE'] = pd.to_datetime(all_df['GD_AV45_2_DATE'])
        all_df['GD_AV45_3_DATE'] = pd.to_datetime(all_df['GD_AV45_3_DATE'])

        # Get closest av1451 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'GDTOTAL', [av1451_date1,av1451_date2,av1451_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av1451_date1,av1451_date2,av1451_date3],day_limit=365/2)
        all_df['GD_AV1451_1'] = closest_vals[0]
        all_df['GD_AV1451_1_DATE'] = closest_dates[0]
        all_df['GD_AV1451_2'] = closest_vals[1]
        all_df['GD_AV1451_2_DATE'] = closest_dates[1]
        all_df['GD_AV1451_3'] = closest_vals[2]
        all_df['GD_AV1451_3_DATE'] = closest_dates[2]
        all_df['GD_AV1451_1_DATE'] = pd.to_datetime(all_df['GD_AV1451_1_DATE'])
        all_df['GD_AV1451_2_DATE'] = pd.to_datetime(all_df['GD_AV1451_2_DATE'])
        all_df['GD_AV1451_3_DATE'] = pd.to_datetime(all_df['GD_AV1451_3_DATE'])


        return all_df

    parsed_df = parseSubjectGroups(gd_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df


def syncADASCogData(master_df, adni1_adas_file, adnigo2_adas_file, registry):
    adas_df = importADASCog(adni1_adas_file, adnigo2_adas_file, registry=registry, as_df=True)
    tmpts = max(Counter(adas_df.index).values())

    headers = ['ADAScog.%s' % (i+1) for i in range(tmpts)]
    headers += ['ADAScog_timePostAV45.%s' % (i+1) for i in range(tmpts)]
    headers += ['ADAS_post_AV45_followuptime','ADASslope_postAV45',
                'ADAS_AV45_1','ADAS_AV45_1_DATE',
                'ADAS_AV45_2','ADAS_AV45_2_DATE',
                'ADAS_AV45_3','ADAS_AV45_3_DATE',
                'ADAS_AV1451_1','ADAS_AV1451_1_DATE',
                'ADAS_AV1451_2','ADAS_AV1451_2_DATE',
                'ADAS_AV1451_3','ADAS_AV1451_3_DATE']

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)

        # Get long measurements
        subj_rows['RID'] = rid
        adas_long = groupLongPivot(subj_rows,'RID','TOTSCORE','ADAScog.')
        date_long = groupLongPivot(subj_rows,'RID','EXAMDATE','ADAScog_timePostAV45.')
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)
        all_df = pd.concat((adas_long,date_long),axis=1)

        adas_slope = df_slope(all_df, date_long.columns, adas_long.columns, take_diff=False, exact=False)
        all_df['ADASslope_postAV45'] = adas_slope[rid]
        last_date = subj_rows.iloc[-1]['EXAMDATE']
        all_df['ADAS_post_AV45_followuptime'] = (last_date - av45_date1).days/365.25 if (not isnan(av45_date1) and last_date > av45_date1) else np.nan

        # Get closest av45 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'TOTSCORE', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        all_df['ADAS_AV45_1'] = closest_vals[0]
        all_df['ADAS_AV45_1_DATE'] = closest_dates[0]
        all_df['ADAS_AV45_2'] = closest_vals[1]
        all_df['ADAS_AV45_2_DATE'] = closest_dates[1]
        all_df['ADAS_AV45_3'] = closest_vals[2]
        all_df['ADAS_AV45_3_DATE'] = closest_dates[2]
        all_df['ADAS_AV45_1_DATE'] = pd.to_datetime(all_df['ADAS_AV45_1_DATE'])
        all_df['ADAS_AV45_2_DATE'] = pd.to_datetime(all_df['ADAS_AV45_2_DATE'])
        all_df['ADAS_AV45_3_DATE'] = pd.to_datetime(all_df['ADAS_AV45_3_DATE'])

        # Get closest av1451 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'TOTSCORE', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        all_df['ADAS_AV1451_1'] = closest_vals[0]
        all_df['ADAS_AV1451_1_DATE'] = closest_dates[0]
        all_df['ADAS_AV1451_2'] = closest_vals[1]
        all_df['ADAS_AV1451_2_DATE'] = closest_dates[1]
        all_df['ADAS_AV1451_3'] = closest_vals[2]
        all_df['ADAS_AV1451_3_DATE'] = closest_dates[2]
        all_df['ADAS_AV1451_1_DATE'] = pd.to_datetime(all_df['ADAS_AV1451_1_DATE'])
        all_df['ADAS_AV1451_2_DATE'] = pd.to_datetime(all_df['ADAS_AV1451_2_DATE'])
        all_df['ADAS_AV1451_3_DATE'] = pd.to_datetime(all_df['ADAS_AV1451_3_DATE'])

        return all_df

    parsed_df = parseSubjectGroups(adas_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df


def syncMMSEData(master_df, mmse_file, registry):
    mmse_df = importMMSE(mmse_file, registry=registry, as_df=True)
    tmpts = max(Counter(mmse_df.index).values())

    headers = ['MMSE_SCORE.%s' % (i+1) for i in range(tmpts)]
    headers += ['MMSE_timePostAV45.%s' % (i+1) for i in range(tmpts)]
    headers += ['MMSE_post_AV45_followuptime',
                'MMSEslope_postAV45',
                'MMSE_AV45_1','MMSE_AV45_1_DATE',
                'MMSE_AV45_2','MMSE_AV45_2_DATE',
                'MMSE_AV45_3','MMSE_AV45_3_DATE']

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)

        # Get long measurements
        subj_rows['RID'] = rid
        mmse_long = groupLongPivot(subj_rows,'RID','MMSCORE','MMSE_SCORE.')
        date_long = groupLongPivot(subj_rows,'RID','EXAMDATE','MMSE_timePostAV45.')
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)
        all_df = pd.concat((mmse_long,date_long),axis=1)

        mmse_slope = df_slope(all_df, date_long.columns, mmse_long.columns, take_diff=False, exact=False)
        all_df['MMSEslope_postAV45'] = mmse_slope[rid]
        last_date = subj_rows.iloc[-1]['EXAMDATE']
        all_df['MMSE_post_AV45_followuptime'] = (last_date - av45_date1).days/365.25 if (not isnan(av45_date1) and last_date > av45_date1) else np.nan

        # Get closest av45 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'MMSCORE', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        all_df['MMSE_AV45_1'] = closest_vals[0]
        all_df['MMSE_AV45_1_DATE'] = closest_dates[0]
        all_df['MMSE_AV45_2'] = closest_vals[1]
        all_df['MMSE_AV45_2_DATE'] = closest_dates[1]
        all_df['MMSE_AV45_3'] = closest_vals[2]
        all_df['MMSE_AV45_3_DATE'] = closest_dates[2]
        all_df['MMSE_AV45_1_DATE'] = pd.to_datetime(all_df['MMSE_AV45_1_DATE'])
        all_df['MMSE_AV45_2_DATE'] = pd.to_datetime(all_df['MMSE_AV45_2_DATE'])
        all_df['MMSE_AV45_3_DATE'] = pd.to_datetime(all_df['MMSE_AV45_3_DATE'])

        return all_df

    parsed_df = parseSubjectGroups(mmse_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df

def syncAVLTData(master_df, neuro_battery_file, registry):
    avlt_df = importAVLT(neuro_battery_file, registry=registry, as_df=True)
    tmpts = max(Counter(avlt_df.index).values())

    headers = ['AVLT.%s' % (i+1) for i in range(tmpts)]
    headers += ['AVLT_timePostAV45.%s' % (i+1) for i in range(tmpts)]
    headers += ['AVLT_AV45_1','AVLT_AV45_1_DATE',
                'AVLT_AV45_2','AVLT_AV45_2_DATE',
                'AVLT_AV45_3','AVLT_AV45_3_DATE',
                'AVLT_AV1451_1','AVLT_AV1451_1_DATE',
                'AVLT_AV1451_2','AVLT_AV1451_2_DATE',
                'AVLT_AV1451_3','AVLT_AV1451_3_DATE',
                'AVLT_post_AV45_followuptime',
                'AVLT_slope_postAV45']

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)

        # Get long measurements
        subj_rows['RID'] = rid
        avlt_long = groupLongPivot(subj_rows,'RID','TOTS','AVLT.')
        date_long = groupLongPivot(subj_rows,'RID','EXAMDATE','AVLT_timePostAV45.')
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)
        all_df = pd.concat((avlt_long,date_long),axis=1)

        avlt_slope = df_slope(all_df, date_long.columns, avlt_long.columns, take_diff=False, exact=False)
        all_df['AVLT_slope_postAV45'] = avlt_slope[rid]
        last_date = subj_rows.iloc[-1]['EXAMDATE']
        all_df['AVLT_post_AV45_followuptime'] = (last_date - av45_date1).days/365.25 if (not isnan(av45_date1) and last_date > av45_date1) else np.nan

        # Get closest av45 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'TOTS', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        all_df['AVLT_AV45_1'], all_df['AVLT_AV45_2'], all_df['AVLT_AV45_3'] = tuple(closest_vals)
        all_df['AVLT_AV45_1_DATE'], all_df['AVLT_AV45_2_DATE'], all_df['AVLT_AV45_3_DATE'] = tuple(closest_dates)
        all_df['AVLT_AV45_1_DATE'] = pd.to_datetime(all_df['AVLT_AV45_1_DATE'])
        all_df['AVLT_AV45_2_DATE'] = pd.to_datetime(all_df['AVLT_AV45_2_DATE'])
        all_df['AVLT_AV45_3_DATE'] = pd.to_datetime(all_df['AVLT_AV45_3_DATE'])

        # Get closest av1451 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'TOTS', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        all_df['AVLT_AV1451_1'], all_df['AVLT_AV1451_2'], all_df['AVLT_AV1451_3'] = tuple(closest_vals)
        all_df['AVLT_AV1451_1_DATE'], all_df['AVLT_AV1451_2_DATE'], all_df['AVLT_AV1451_3_DATE'] = tuple(closest_dates)
        all_df['AVLT_AV1451_1_DATE'] = pd.to_datetime(all_df['AVLT_AV1451_1_DATE'])
        all_df['AVLT_AV1451_2_DATE'] = pd.to_datetime(all_df['AVLT_AV1451_2_DATE'])
        all_df['AVLT_AV1451_3_DATE'] = pd.to_datetime(all_df['AVLT_AV1451_3_DATE'])

        return all_df

    parsed_df = parseSubjectGroups(avlt_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df


def syncDemogData(master_df, demog_file, pet_meta_file):
    # get demographics data
    demog_df = importDemog(demog_file,as_df=True)
    # merge in pet scan dates
    av45_meta_df = importPetMETA(pet_meta_file, 'AV45', as_df=True)
    av45_meta_df.rename(columns={1:'AV45_Date',2:'AV45_2_Date',3:'AV45_3_Date',4:'AV45_4_Date'}, inplace=True)
    av1451_meta_df = importPetMETA(pet_meta_file, 'AV1451', as_df=True)
    av1451_meta_df.rename(columns={1:'AV1451_Date',2:'AV1451_2_Date',3:'AV1451_3_Date',4:'AV1451_4_Date'}, inplace=True)
    demog_df = demog_df.merge(av45_meta_df, left_index=True, right_index=True, how='outer')
    demog_df = demog_df.merge(av1451_meta_df, left_index=True, right_index=True, how='outer')

    # calculate derivative columns
    demog_df['Age@AV45'] = demog_df.apply(lambda x: yrDiff(x['AV45_Date'],x['PTDOB']), axis=1)
    demog_df['Age@AV45_2'] = demog_df.apply(lambda x: yrDiff(x['AV45_2_Date'],x['PTDOB']), axis=1)
    demog_df['Age@AV45_3'] = demog_df.apply(lambda x: yrDiff(x['AV45_3_Date'],x['PTDOB']), axis=1)
    demog_df['Age@AV45_4'] = demog_df.apply(lambda x: yrDiff(x['AV45_4_Date'],x['PTDOB']), axis=1)
    demog_df['Age@AV1451'] = demog_df.apply(lambda x: yrDiff(x['AV1451_Date'],x['PTDOB']), axis=1)
    demog_df['Age@AV1451_2'] = demog_df.apply(lambda x: yrDiff(x['AV1451_2_Date'],x['PTDOB']), axis=1)
    demog_df['Age@AV1451_3'] = demog_df.apply(lambda x: yrDiff(x['AV1451_3_Date'],x['PTDOB']), axis=1)
    demog_df['AV45_1_2_Diff'] = demog_df.apply(lambda x: yrDiff(x['AV45_2_Date'],x['AV45_Date']), axis=1)
    demog_df['AV45_1_3_Diff'] = demog_df.apply(lambda x: yrDiff(x['AV45_3_Date'],x['AV45_Date']), axis=1)
    demog_df['AV45_1_4_Diff'] = demog_df.apply(lambda x: yrDiff(x['AV45_4_Date'],x['AV45_Date']), axis=1)
    demog_df['AV1451_1_2_Diff'] = demog_df.apply(lambda x: yrDiff(x['AV1451_2_Date'],x['AV1451_Date']), axis=1)
    demog_df['AV1451_1_3_Diff'] = demog_df.apply(lambda x: yrDiff(x['AV1451_3_Date'],x['AV1451_Date']), axis=1)

    headers = ['AV45_Date','AV45_2_Date','AV45_3_Date','AV45_4_Date',
               'AV1451_Date','AV1451_2_Date','AV1451_3_Date',
               'PTDOB','Age@AV45','Age@AV45_2','Age@AV45_3','Age@AV45_4',
               'Age@AV1451','Age@AV1451_2','Age@AV1451_3',
               'AV45_1_2_Diff','AV45_1_3_Diff','AV45_1_4_Diff',
               'AV1451_1_2_Diff','AV1451_1_3_Diff']
    master_df = updateDataFrame(master_df, demog_df, headers=headers, after='APOE4_NUM', restrict=True)
    return master_df


def syncDiagnosisData(master_df, diag_file, arm_file, registry):
    diag_df = importADNIDiagnosis(diag_file, registry=registry, as_df=True)
    arm_df = importARM(arm_file, as_df=True)

    pivot_date = datetime(day=1, month=2, year=2016)
    pivot_date_closest_diag_key = 'Closest_DX_Feb16'
    pivot_date_closest_date_key = 'DX_Feb16_closestdate'
    headers = ['MCItoADConv','MCItoADConvDate',
               'AV45_MCItoAD_ConvTime','Baseline_MCItoAD_ConvTime',
               'Diag@AV45','Diag@AV45_2','Diag@AV45_3',
               'Diag@AV1451','Diag@AV1451_2','Diag@AV1451_3',
               'FollowupTimetoDX','Baseline','Init_Diagnosis',
               pivot_date_closest_diag_key,pivot_date_closest_date_key]

    # get initial diagnoses
    index_name = arm_df.index.name
    arm_df.reset_index(inplace=True)
    init_diag_df = arm_df[arm_df.groupby(index_name)['USERDATE'].rank(ascending=False) == 1].set_index(index_name)

    # get initial visits
    try:
        bad_vc2 = registry.loc(axis=0)[:,:,'scmri'].index
    except:
        bad_vc2 = []
    try:
        bad_vc = registry.loc(axis=0)[:,'scmri',:].index
    except:
        bad_vc = []
    bad_vc2_idx = registry.index.isin(bad_vc2)
    bad_vc_idx = registry.index.isin(bad_vc)
    bad_idx = bad_vc2_idx | bad_vc_idx
    bl_registry = registry[~bad_idx]
    bl_registry.reset_index(level=0,inplace=True)
    bl_registry = bl_registry.groupby('RID')['EXAMDATE'].aggregate(min)

    # add comparison time diffs
    diag_df['TimeFromCompDate'] = diag_df['EXAMDATE'].apply(lambda x: abs(x-pivot_date).days)

    def extraction_fn(rid, subj_rows):
        try:
            init_diag = init_diag_df.loc[rid,'STATUS']
        except:
            init_diag = None
        bl_visit = bl_registry[rid]
        mci_type = 'LMCI'
        if init_diag in set(['LMCI', 'EMCI', 'SMC']):
            mci_type = init_diag
        subj_rows.sort_values(by='EXAMDATE', inplace=True)
        subj_rows.loc[subj_rows.diag == 'MCI','diag'] = mci_type
        data = {'RID': rid,
                'Init_Diagnosis': init_diag,
                'Baseline': bl_visit,
                'FollowupTimetoDX': (pivot_date-bl_visit).days/365.25}
        # closest to AV45 scans
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'diag', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        data['Diag@AV45'] = closest_vals[0]
        data['Diag@AV45_2'] = closest_vals[1]
        data['Diag@AV45_3'] = closest_vals[2]
        # closest to AV1451 scans
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'diag', [av1451_date1,av1451_date2,av1451_date3],day_limit=365/2)
        data['Diag@AV1451'] = closest_vals[0]
        data['Diag@AV1451_2'] = closest_vals[1]
        data['Diag@AV1451_3'] = closest_vals[2]
        # closest to pivot date
        closest_row = subj_rows.sort_values('TimeFromCompDate').iloc[0]
        data[pivot_date_closest_diag_key] = closest_row['diag']
        data[pivot_date_closest_date_key] = closest_row['EXAMDATE']
        # mci to ad conversion
        conv_row = subj_rows[subj_rows.change == 5].head(1)
        if len(conv_row.index) == 1:
            conv_date = conv_row.iloc[0]['EXAMDATE']
            data['MCItoADConv'] = 1
            data['MCItoADConvDate'] = conv_date
            data['Baseline_MCItoAD_ConvTime'] = (conv_date-bl_visit).days/365.25
            if not isnan(av45_date1):
                data['AV45_MCItoAD_ConvTime'] = (conv_date-av45_date1).days/365.25
        else:
            data['MCItoADConv'] = 0
        # return
        return pd.DataFrame([data]).set_index('RID')

    parsed_df = parseSubjectGroups(diag_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after='AV45_1_4_Diff')
    return master_df

def syncAV1451Data(master_df, av1451_file):
    av1451_df = importAV1451(av1451_file, as_df=True)

    valid_timepoints = ['BL', 'Scan2', 'Scan3']
    headers = ['AV1451_Braak12_CerebGray_%s' % tp for tp in valid_timepoints]
    headers += ['AV1451_Braak34_CerebGray_%s' % tp for tp in valid_timepoints]
    headers += ['AV1451_Braak56_CerebGray_%s' % tp for tp in valid_timepoints]

    after = [_ for _ in master_df.columns if _.startswith('AV45_') and 'PVC' not in _][-1]

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE', inplace=True)
        data = {'RID': int(rid)}
        for i, (idx, row) in enumerate(subj_rows.iterrows()):
            if i == 0:
                tp = 'BL'
            elif i == 1:
                tp = 'Scan2'
            elif i == 2:
                tp = 'Scan3'
            else:
                raise Exception("Raise # of AV1451 timepoints")

            cerebg_left = (row['LEFT_CEREBELLUM_CORTEX_SIZE'],row['LEFT_CEREBELLUM_CORTEX'])
            cerebg_right = (row['RIGHT_CEREBELLUM_CORTEX_SIZE'],row['RIGHT_CEREBELLUM_CORTEX'])
            braak1 = (row['BRAAK1_SIZE'],row['BRAAK1'])
            braak2 = (row['BRAAK2_SIZE'],row['BRAAK2'])
            braak3 = (row['BRAAK3_SIZE'],row['BRAAK3'])
            braak4 = (row['BRAAK4_SIZE'],row['BRAAK4'])
            braak5 = (row['BRAAK5_SIZE'],row['BRAAK5'])
            braak6 = (row['BRAAK6_SIZE'],row['BRAAK6'])

            cerebg = weightedMean([cerebg_left,cerebg_right])
            braak12 = weightedMean([braak1,braak2]) / cerebg
            braak34 = weightedMean([braak3,braak4]) / cerebg
            braak56 = weightedMean([braak5,braak6]) / cerebg

            data['AV1451_Braak12_CerebGray_%s' % tp] = braak12
            data['AV1451_Braak34_CerebGray_%s' % tp] = braak34
            data['AV1451_Braak56_CerebGray_%s' % tp] = braak56
        return pd.DataFrame([data]).set_index('RID')

    parsed_df = parseSubjectGroups(av1451_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=after, restrict=False)
    return master_df


def syncAV45Data(master_df, av45_file, suffix=None):
    '''
    This function does not follow the pattern of the other sync functions because
    the update is accomplished in a nested function, and it also allows
    for the possibility of introducing new subjects via the dataset being synced in
    '''
    av45_df = importAV45(av45_file, as_df=True)

    # parse
    parsed_df = parseSubjectGroups(av45_df,lambda rid,subj_rows: parseAV45Entries(rid, subj_rows, suffix=suffix))
    headers = parsed_df.columns

    # update
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None, restrict=False)
    return master_df


def parseAV45Entries(rid, subj_rows, suffix=None):
    subj_rows.sort_values(by='EXAMDATE', inplace=True)
    exam_times = list(subj_rows['EXAMDATE'])
    exam_timedeltas = [(_-exam_times[0]).days / 365.0 for _ in exam_times]

    wm70_composite_keys = ['AV45_WM70/composite', 'AV45_2_WM70/composite', 'AV45_3_WM70/composite']
    wm70_cerebg_keys = ['AV45_WM70/cerebg', 'AV45_2_WM70/cerebg', 'AV45_3_WM70/cerebg']
    wm70_wcereb_keys = ['AV45_WM70/wcereb', 'AV45_2_WM70/wcereb', 'AV45_3_WM70/wcereb']
    unilateral_keys = ['AV45_LeftPutamen/WM70','AV45_RightPutamen/WM70',
                       'AV45_LeftCaudate/WM70','AV45_RightCaudate/WM70',
                       'AV45_LeftPallidum/WM70','AV45_RightPallidum/WM70']
    bigref_keys = ['AV45_BigRef','AV45_2_BigRef','AV45_3_BigRef'] # assuming composite ROI
    wm70_keys = ['AV45_WM70','AV45_2_WM70','AV45_3_WM70'] # assuming composite ROI
    cerebg_keys = ['AV45_cerebg','AV45_2_cerebg','AV45_3_cerebg'] # assuming composite ROI
    wcereb_keys = ['AV45_wcereb','AV45_2_wcereb','AV45_3_wcereb'] # assuming composite ROI
    brainstem_keys = ['AV45_brainstem','AV45_2_brainstem','AV45_3_brainstem'] # assuming composite ROI
    wmratio_keys = ['AV45_WMratio','AV45_2_WMratio','AV45_3_WMratio']
    frontal_bigref_keys = ['AV45_Frontal/BigRef','AV45_2_Frontal/BigRef','AV45_3_Frontal/BigRef']
    cingulate_bigref_keys = ['AV45_Cingulate/BigRef','AV45_2_Cingulate/BigRef','AV45_3_Cingulate/BigRef']
    parietal_bigref_keys = ['AV45_Parietal/BigRef','AV45_2_Parietal/BigRef','AV45_3_Parietal/BigRef']
    temporal_bigref_keys = ['AV45_Temporal/BigRef','AV45_2_Temporal/BigRef','AV45_3_Temporal/BigRef']

    # generate additional keys and arrange into header list
    all_wm70_composite_keys = wm70_composite_keys + \
                         ["%s_pchange" % _ for _ in wm70_composite_keys[1:]] + \
                         ["%s_pchange_ABS" % _ for _ in wm70_composite_keys[1:]] + \
                         ["%s_diff" % _ for _ in wm70_composite_keys[1:]] + \
                         ["%s_diff_ABS" % _ for _ in wm70_composite_keys[1:]] + \
                         ['AV45_WM70/composite_Slope_2pts', 'AV45_WM70/composite_Slope_3pts'] + \
                         ['AV45_WM70/composite_Slope_1and3', 'AV45_WM70/composite_Slope_2and3']
    all_wm70_cerebg_keys = wm70_cerebg_keys + \
                           ["%s_pchange" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ["%s_pchange_ABS" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ["%s_diff" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ["%s_diff_ABS" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ['AV45_WM70/cerebg_Slope_2pts', 'AV45_WM70/cerebg_Slope_3pts'] + \
                           ['AV45_WM70/cerebg_Slope_1and3', 'AV45_WM70/cerebg_Slope_2and3']
    all_wm70_wcereb_keys = wm70_wcereb_keys + \
                           ["%s_pchange" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ["%s_pchange_ABS" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ["%s_diff" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ["%s_diff_ABS" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ['AV45_WM70/wcereb_Slope_2pts', 'AV45_WM70/wcereb_Slope_3pts'] + \
                           ['AV45_WM70/wcereb_Slope_1and3', 'AV45_WM70/wcereb_Slope_2and3']
    all_unilateral_keys = unilateral_keys
    all_bigref_keys = bigref_keys + \
                      ["%s_pchange" % _ for _ in bigref_keys[1:]] + \
                      ["%s_pchange_ABS" % _ for _ in bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in bigref_keys[1:]] + \
                      ["%s_diff_ABS" % _ for _ in bigref_keys[1:]] + \
                      ["%s_BIN.79" % _ for _ in bigref_keys] + \
                      ['AV45_BigRef_Slope_2pts', 'AV45_BigRef_Slope_3pts'] + \
                      ['AV45_BigRef_Slope_1and3', 'AV45_BigRef_Slope_2and3']
    all_wm70_keys = wm70_keys + \
                    ["%s_pchange" % _ for _ in wm70_keys[1:]] + \
                    ["%s_pchange_ABS" % _ for _ in wm70_keys[1:]] + \
                    ["%s_diff" % _ for _ in wm70_keys[1:]] + \
                    ["%s_diff_ABS" % _ for _ in wm70_keys[1:]] + \
                    ["%s_BIN.62" % _ for _ in wm70_keys] + \
                    ['AV45_WM70_Slope_2pts', 'AV45_WM70_Slope_3pts'] + \
                    ['AV45_WM70_Slope_1and3', 'AV45_WM70_Slope_2and3']
    all_cerebg_keys = cerebg_keys + \
                    ["%s_pchange" % _ for _ in cerebg_keys[1:]] + \
                    ["%s_pchange_ABS" % _ for _ in cerebg_keys[1:]] + \
                    ["%s_diff" % _ for _ in cerebg_keys[1:]] + \
                    ["%s_diff_ABS" % _ for _ in cerebg_keys[1:]] + \
                    ["%s_BIN1.26" % _ for _ in cerebg_keys] + \
                    ['AV45_cerebg_Slope_2pts', 'AV45_cerebg_Slope_3pts'] + \
                    ['AV45_cerebg_Slope_1and3', 'AV45_cerebg_Slope_2and3']
    all_wcereb_keys = wcereb_keys + \
                    ["%s_pchange" % _ for _ in wcereb_keys[1:]] + \
                    ["%s_pchange_ABS" % _ for _ in wcereb_keys[1:]] + \
                    ["%s_diff" % _ for _ in wcereb_keys[1:]] + \
                    ["%s_diff_ABS" % _ for _ in wcereb_keys[1:]] + \
                    ["%s_BIN1.11" % _ for _ in wcereb_keys] + \
                    ['AV45_wcereb_Slope_2pts', 'AV45_wcereb_Slope_3pts'] + \
                    ['AV45_wcereb_Slope_1and3', 'AV45_wcereb_Slope_2and3']
    all_brainstem_keys = brainstem_keys + \
                    ["%s_pchange" % _ for _ in brainstem_keys[1:]] + \
                    ["%s_pchange_ABS" % _ for _ in brainstem_keys[1:]] + \
                    ["%s_diff" % _ for _ in brainstem_keys[1:]] + \
                    ["%s_diff_ABS" % _ for _ in brainstem_keys[1:]] + \
                    ["%s_BIN.79" % _ for _ in brainstem_keys] + \
                    ['AV45_brainstem_Slope_2pts', 'AV45_brainstem_Slope_3pts'] + \
                    ['AV45_brainstem_Slope_1and3', 'AV45_brainstem_Slope_2and3']
    all_wmratio_keys = wmratio_keys + \
                      ["%s_pchange" % _ for _ in wmratio_keys[1:]] + \
                      ["%s_pchange_ABS" % _ for _ in wmratio_keys[1:]] + \
                      ["%s_diff" % _ for _ in wmratio_keys[1:]] + \
                      ["%s_diff_ABS" % _ for _ in wmratio_keys[1:]] + \
                      ['AV45_WMratio_Slope_2pts', 'AV45_WMratio_Slope_3pts'] + \
                      ['AV45_WMratio_Slope_1and3', 'AV45_WMratio_Slope_2and3']
    all_frontal_bigref_keys = frontal_bigref_keys + \
                      ["%s_pchange" % _ for _ in frontal_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in frontal_bigref_keys[1:]] + \
                      ['AV45_Frontal/BigRef_Slope_2pts', 'AV45_Frontal/BigRef_Slope_3pts'] + \
                      ['AV45_Frontal/BigRef_Slope_1and3', 'AV45_Frontal/BigRef_Slope_2and3']
    all_cingulate_bigref_keys = cingulate_bigref_keys + \
                      ["%s_pchange" % _ for _ in cingulate_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in cingulate_bigref_keys[1:]] + \
                      ['AV45_Cingulate/BigRef_Slope_2pts', 'AV45_Cingulate/BigRef_Slope_3pts'] + \
                      ['AV45_Cingulate/BigRef_Slope_1and3', 'AV45_Cingulate/BigRef_Slope_2and3']
    all_parietal_bigref_keys = parietal_bigref_keys + \
                      ["%s_pchange" % _ for _ in parietal_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in parietal_bigref_keys[1:]] + \
                      ['AV45_Parietal/BigRef_Slope_2pts', 'AV45_Parietal/BigRef_Slope_3pts'] + \
                      ['AV45_Parietal/BigRef_Slope_1and3', 'AV45_Parietal/BigRef_Slope_2and3']
    all_temporal_bigref_keys = temporal_bigref_keys + \
                      ["%s_pchange" % _ for _ in temporal_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in temporal_bigref_keys[1:]] + \
                      ['AV45_Temporal/BigRef_Slope_2pts', 'AV45_Temporal/BigRef_Slope_3pts'] + \
                      ['AV45_Temporal/BigRef_Slope_1and3', 'AV45_Temporal/BigRef_Slope_2and3']
    all_av45_key_lists = [all_wm70_composite_keys,
                          all_wm70_cerebg_keys,
                          all_wm70_wcereb_keys,
                          all_unilateral_keys,
                          all_bigref_keys,
                          all_wm70_keys,
                          all_cerebg_keys,
                          all_wcereb_keys,
                          all_brainstem_keys,
                          all_wmratio_keys,
                          all_frontal_bigref_keys,
                          all_cingulate_bigref_keys,
                          all_parietal_bigref_keys,
                          all_temporal_bigref_keys]
    pchange_diff_lists = [all_wm70_composite_keys,
                          all_wm70_cerebg_keys,
                          all_wm70_wcereb_keys,
                          all_bigref_keys,
                          all_wm70_keys,
                          all_cerebg_keys,
                          all_wcereb_keys,
                          all_brainstem_keys,
                          all_wmratio_keys,
                          all_frontal_bigref_keys,
                          all_cingulate_bigref_keys,
                          all_parietal_bigref_keys,
                          all_temporal_bigref_keys]
    abs_lists = [all_wm70_composite_keys,
                 all_wm70_cerebg_keys,
                 all_wm70_wcereb_keys,
                 all_bigref_keys,
                 all_wm70_keys,
                 all_cerebg_keys,
                 all_wcereb_keys,
                 all_brainstem_keys,
                 all_wmratio_keys]
    all_av45_keys = [_ for l in all_av45_key_lists for _ in l]

    data = {k:None for k in all_av45_keys}
    wm70_0 = None
    # fill in values
    for i, (idx,point) in enumerate(subj_rows.iterrows()):
        # extract necessary values
        cerebw = float(point['CEREBELLUMWHITEMATTER'])
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

        # fill in basic keys
        data[wm70_composite_keys[i]] = wm70/compositeroi
        data[wm70_cerebg_keys[i]] = wm70/cerebg
        data[wm70_wcereb_keys[i]] = wm70/wcereb
        data[bigref_keys[i]] = compositeroi/bigref
        data[wm70_keys[i]] = compositeroi/wm70
        data[cerebg_keys[i]] = compositeroi/cerebg
        data[wcereb_keys[i]] = compositeroi/wcereb
        data[brainstem_keys[i]] = compositeroi/brainstem
        data[frontal_bigref_keys[i]] = frontal/bigref
        data[cingulate_bigref_keys[i]] = cingulate/bigref
        data[parietal_bigref_keys[i]] = parietal/bigref
        data[temporal_bigref_keys[i]] = temporal/bigref
        data["%s_BIN.79" % bigref_keys[i]] = 1 if (compositeroi/bigref) > 0.79 else 0
        data["%s_BIN.62" % wm70_keys[i]] = 1 if (compositeroi/wm70) > 0.62 else 0
        data["%s_BIN1.26" % cerebg_keys[i]] = 1 if (compositeroi/cerebg) > 1.26 else 0
        data["%s_BIN1.11" % wcereb_keys[i]] = 1 if (compositeroi/wcereb) > 1.11 else 0
        data["%s_BIN.79" % brainstem_keys[i]] = 1 if (compositeroi/brainstem) > 0.79 else 0


        # fill in derivative keys
        if i == 0:
            # fill in unilateral
            data['AV45_LeftPutamen/WM70'] = leftputamen/wm70
            data['AV45_RightPutamen/WM70'] = rightputamen/wm70
            data['AV45_LeftCaudate/WM70'] = leftcaudate/wm70
            data['AV45_RightCaudate/WM70'] = rightcaudate/wm70
            data['AV45_LeftPallidum/WM70'] = leftpallidum/wm70
            data['AV45_RightPallidum/WM70'] = rightpallidum/wm70
            data[wmratio_keys[i]] = (compositeroi/wcereb)
            wm70_0 = wm70
        elif i == 1 or i == 2:
            data[wmratio_keys[i]] = (compositeroi/wcereb) / (wm70/wm70_0)
            for pl in pchange_diff_lists:
                diff = (data[pl[i]] - data[pl[0]]) / exam_timedeltas[i] # annualized
                data["%s_pchange" % pl[i]] = diff / data[pl[0]]
                data["%s_diff" % pl[i]] = diff
            for al in abs_lists:
                data["%s_pchange_ABS" % al[i]] = abs(data["%s_pchange" % al[i]])
                data["%s_diff_ABS" % al[i]] = abs(data["%s_diff" % al[i]])
            if i == 1:
                times = exam_timedeltas[:2]
                data['AV45_BigRef_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in bigref_keys[:2]])))
                data['AV45_WM70_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in wm70_keys[:2]])))
                data['AV45_cerebg_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in cerebg_keys[:2]])))
                data['AV45_wcereb_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in wcereb_keys[:2]])))
                data['AV45_brainstem_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in brainstem_keys[:2]])))
                data['AV45_WM70/composite_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in wm70_composite_keys[:2]])))
                data['AV45_WM70/cerebg_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in wm70_cerebg_keys[:2]])))
                data['AV45_WM70/wcereb_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in wm70_wcereb_keys[:2]])))
                data['AV45_WMratio_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in wmratio_keys[:2]])))
                data['AV45_Frontal/BigRef_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in frontal_bigref_keys[:2]])))
                data['AV45_Cingulate/BigRef_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in cingulate_bigref_keys[:2]])))
                data['AV45_Parietal/BigRef_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in parietal_bigref_keys[:2]])))
                data['AV45_Temporal/BigRef_Slope_2pts'] = slope(list(zip(times,[data[_] for _ in temporal_bigref_keys[:2]])))
            if i == 2:
                # using all 3 timepoints
                times = exam_timedeltas[:3]
                data['AV45_BigRef_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in bigref_keys[:3]])))
                data['AV45_WM70_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in wm70_keys[:3]])))
                data['AV45_cerebg_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in cerebg_keys[:3]])))
                data['AV45_wcereb_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in wcereb_keys[:3]])))
                data['AV45_brainstem_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in brainstem_keys[:3]])))
                data['AV45_WM70/composite_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in wm70_composite_keys[:3]])))
                data['AV45_WM70/cerebg_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in wm70_cerebg_keys[:3]])))
                data['AV45_WM70/wcereb_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in wm70_wcereb_keys[:3]])))
                data['AV45_WMratio_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in wmratio_keys[:3]])))
                data['AV45_Frontal/BigRef_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in frontal_bigref_keys[:3]])))
                data['AV45_Cingulate/BigRef_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in cingulate_bigref_keys[:3]])))
                data['AV45_Parietal/BigRef_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in parietal_bigref_keys[:3]])))
                data['AV45_Temporal/BigRef_Slope_3pts'] = slope(list(zip(times,[data[_] for _ in temporal_bigref_keys[:3]])))
                # using second and third timepoints
                times = exam_timedeltas[1:3]
                data['AV45_BigRef_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in bigref_keys[1:3]])))
                data['AV45_WM70_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in wm70_keys[1:3]])))
                data['AV45_cerebg_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in cerebg_keys[1:3]])))
                data['AV45_wcereb_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in wcereb_keys[1:3]])))
                data['AV45_brainstem_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in brainstem_keys[1:3]])))
                data['AV45_WM70/composite_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in wm70_composite_keys[1:3]])))
                data['AV45_WM70/cerebg_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in wm70_cerebg_keys[1:3]])))
                data['AV45_WM70/wcereb_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in wm70_wcereb_keys[1:3]])))
                data['AV45_WMratio_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in wmratio_keys[1:3]])))
                data['AV45_Frontal/BigRef_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in frontal_bigref_keys[1:3]])))
                data['AV45_Cingulate/BigRef_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in cingulate_bigref_keys[1:3]])))
                data['AV45_Parietal/BigRef_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in parietal_bigref_keys[1:3]])))
                data['AV45_Temporal/BigRef_Slope_2and3'] = slope(list(zip(times,[data[_] for _ in temporal_bigref_keys[1:3]])))
                # using first and third timepoints
                times = exam_timedeltas[0::2]
                data['AV45_BigRef_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in bigref_keys[0::2]])))
                data['AV45_WM70_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in wm70_keys[0::2]])))
                data['AV45_cerebg_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in cerebg_keys[0::2]])))
                data['AV45_wcereb_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in wcereb_keys[0::2]])))
                data['AV45_brainstem_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in brainstem_keys[0::2]])))
                data['AV45_WM70/composite_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in wm70_composite_keys[0::2]])))
                data['AV45_WM70/cerebg_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in wm70_cerebg_keys[0::2]])))
                data['AV45_WM70/wcereb_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in wm70_wcereb_keys[0::2]])))
                data['AV45_WMratio_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in wmratio_keys[0::2]])))
                data['AV45_Frontal/BigRef_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in frontal_bigref_keys[0::2]])))
                data['AV45_Cingulate/BigRef_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in cingulate_bigref_keys[0::2]])))
                data['AV45_Parietal/BigRef_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in parietal_bigref_keys[0::2]])))
                data['AV45_Temporal/BigRef_Slope_1and3'] = slope(list(zip(times,[data[_] for _ in temporal_bigref_keys[0::2]])))

    if suffix is not None:
        new_keys = []
        for old_key in all_av45_keys:
            new_key = old_key.replace('AV45_','AV45_%s_' % suffix)
            data[new_key] = data.pop(old_key)
            new_keys.append(new_key)
        all_av45_keys = new_keys


    df = pd.DataFrame([data])[all_av45_keys]
    df['RID'] = int(rid)
    df.set_index('RID',inplace=True)
    return df

def syncFDGData(master_df, fdg_file, registry):
    fdg_df = importFDG(fdg_file, registry, as_df=True)
    tmpts = max(fdg_df.reset_index().groupby('RID')['EXAMDATE'].nunique())
    pons_df = fdg_df.reset_index().groupby(['RID','EXAMDATE'])['MEAN'].aggregate(np.mean).reset_index(level=1)

    headers = ['FDG_PONS.%s' % (i+1) for i in range(tmpts)]
    headers += ['FDG_PONS_ReltoAV45.%s' % (i+1) for i in range(tmpts)]
    headers += ['FDG_postAV45_slope','FDG_postAV45_followuptime']
    headers += ['FDG_PONS_AV45_1', 'FDG_PONS_AV45_1_DATE',
                'FDG_PONS_AV45_2', 'FDG_PONS_AV45_2_DATE',
                'FDG_PONS_AV45_3', 'FDG_PONS_AV45_3_DATE',
                'FDG_PONS_AV1451_1','FDG_PONS_AV1451_1_DATE',
                'FDG_PONS_AV1451_2','FDG_PONS_AV1451_2_DATE',
                'FDG_PONS_AV1451_3','FDG_PONS_AV1451_3_DATE',]

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)

        # Get long measurements
        subj_rows['RID'] = rid
        pons_long = groupLongPivot(subj_rows,'RID','MEAN','FDG_PONS.')
        date_long = groupLongPivot(subj_rows,'RID','EXAMDATE','FDG_PONS_ReltoAV45.')
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)
        all_df = pd.concat((pons_long,date_long),axis=1)

        # Get closest av45 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'MEAN', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        all_df['FDG_PONS_AV45_1'] = closest_vals[0]
        all_df['FDG_PONS_AV45_1_DATE'] = closest_dates[0]
        all_df['FDG_PONS_AV45_2'] = closest_vals[1]
        all_df['FDG_PONS_AV45_2_DATE'] = closest_dates[1]
        all_df['FDG_PONS_AV45_3'] = closest_vals[2]
        all_df['FDG_PONS_AV45_3_DATE'] = closest_dates[2]
        all_df['FDG_PONS_AV45_1_DATE'] = pd.to_datetime(all_df['FDG_PONS_AV45_1_DATE'])
        all_df['FDG_PONS_AV45_2_DATE'] = pd.to_datetime(all_df['FDG_PONS_AV45_2_DATE'])
        all_df['FDG_PONS_AV45_3_DATE'] = pd.to_datetime(all_df['FDG_PONS_AV45_3_DATE'])

        # Get closest AV1451 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'HC/ICV_BL', [av1451_date1,av1451_date2,av1451_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av1451_date1,av1451_date2,av1451_date3],day_limit=365/2)
        all_df['FDG_PONS_AV1451_1'] = closest_vals[0]
        all_df['FDG_PONS_AV1451_1_DATE'] = closest_dates[0]
        all_df['FDG_PONS_AV1451_2'] = closest_vals[1]
        all_df['FDG_PONS_AV1451_2_DATE'] = closest_dates[1]
        all_df['FDG_PONS_AV1451_3'] = closest_vals[2]
        all_df['FDG_PONS_AV1451_3_DATE'] = closest_dates[2]
        all_df['FDG_PONS_AV1451_1_DATE'] = pd.to_datetime(all_df['FDG_PONS_AV1451_1_DATE'])
        all_df['FDG_PONS_AV1451_2_DATE'] = pd.to_datetime(all_df['FDG_PONS_AV1451_2_DATE'])
        all_df['FDG_PONS_AV1451_3_DATE'] = pd.to_datetime(all_df['FDG_PONS_AV1451_3_DATE'])

        # Get slopes + interval
        pons_slope = df_slope(all_df, date_long.columns, pons_long.columns, take_diff=False, exact=False)
        all_df['FDG_postAV45_slope'] = pons_slope[rid]
        last_date = subj_rows.iloc[-1]['EXAMDATE']
        all_df['FDG_postAV45_followuptime'] = (last_date - av45_date1).days/365.25 if (not isnan(av45_date1) and last_date > av45_date1) else np.nan

        return all_df

    parsed_df = parseSubjectGroups(pons_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df

def syncTBMSynData(master_df, tbm_file):
    tbm_df = importTBMSyn(tbm_file, as_df=True)
    tmpts = max(Counter(tbm_df.index).values())

    headers = ['TBMSyn_SCORE.%s' % (i+1) for i in range(tmpts)]
    headers += ['TBMSyn_postAV45.%s' % (i+1) for i in range(tmpts)]
    headers += ['TBMSyn_BL_DATE', 'TBMSyn_count']

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)
        bl_examdate = subj_rows.iloc[0]['EXAMDATEBL']
        if isnan(av45_date1) or abs(av45_date1 - bl_examdate).days > 210:
            return pd.DataFrame()

        # Get long measurements
        subj_rows['RID'] = rid
        tbm_long = groupLongPivot(subj_rows,'RID','TBMSYNSCOR','TBMSyn_SCORE.')
        date_long = groupLongPivot(subj_rows,'RID','EXAMDATE','TBMSyn_postAV45.')
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)
        all_df = pd.concat((tbm_long,date_long),axis=1)

        all_df['TBMSyn_count'] = df_count(all_df,date_long.columns,tbm_long.columns)[rid]
        all_df['TBMSyn_BL_DATE'] = bl_examdate

        return all_df

    parsed_df = parseSubjectGroups(tbm_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df

def syncCSFData(master_df, csf_files, registry):
    csf_df = importCSF(csf_files, registry, as_df=True)
    csf_df = csf_df[['EXAMDATE','ABETA','TAU','PTAU']]
    csf_df['ABETA_BIN'] = (csf_df['ABETA'] >= 192).astype(int)
    csf_df['TAU_BIN'] = (csf_df['TAU'] >= 93).astype(int)
    csf_df['PTAU_BIN'] = (csf_df['PTAU'] >= 23).astype(int)
    tmpts = max(Counter(csf_df.index).values())

    headers = ['CSF_ABETA.%s' % (i+1) for i in range(tmpts)]
    headers += ['CSF_ABETA_slope',
                'CSF_ABETA_closest_AV45_1',
                'CSF_ABETA_closest_AV45_2',
                'CSF_ABETA_closest_AV45_3',
                'CSF_ABETA_closest_AV45_1_BIN_192',
                'CSF_ABETA_closest_AV45_2_BIN_192',
                'CSF_ABETA_closest_AV45_3_BIN_192',
                'CSF_ABETA_closest_AV1451_1',
                'CSF_ABETA_closest_AV1451_2',
                'CSF_ABETA_closest_AV1451_3',
                'CSF_ABETA_closest_AV1451_1_BIN_192',
                'CSF_ABETA_closest_AV1451_2_BIN_192',
                'CSF_ABETA_closest_AV1451_3_BIN_192']
    headers += ['CSF_TAU.%s' % (i+1) for i in range(tmpts)]
    headers += ['CSF_TAU_slope',
                'CSF_TAU_closest_AV45_1',
                'CSF_TAU_closest_AV45_2',
                'CSF_TAU_closest_AV45_3',
                'CSF_TAU_closest_AV45_1_BIN_93',
                'CSF_TAU_closest_AV45_2_BIN_93',
                'CSF_TAU_closest_AV45_3_BIN_93',
                'CSF_TAU_closest_AV1451_1',
                'CSF_TAU_closest_AV1451_2',
                'CSF_TAU_closest_AV1451_3',
                'CSF_TAU_closest_AV1451_1_BIN_93',
                'CSF_TAU_closest_AV1451_2_BIN_93',
                'CSF_TAU_closest_AV1451_3_BIN_93']
    headers += ['CSF_PTAU.%s' % (i+1) for i in range(tmpts)]
    headers += ['CSF_PTAU_slope',
                'CSF_PTAU_closest_AV45_1',
                'CSF_PTAU_closest_AV45_2',
                'CSF_PTAU_closest_AV45_3',
                'CSF_PTAU_closest_AV45_1_BIN_23',
                'CSF_PTAU_closest_AV45_2_BIN_23',
                'CSF_PTAU_closest_AV45_3_BIN_23',
                'CSF_PTAU_closest_AV1451_1',
                'CSF_PTAU_closest_AV1451_2',
                'CSF_PTAU_closest_AV1451_3',
                'CSF_PTAU_closest_AV1451_1_BIN_23',
                'CSF_PTAU_closest_AV1451_2_BIN_23',
                'CSF_PTAU_closest_AV1451_3_BIN_23']
    headers += ['CSF_postAV45.%s' % (i+1) for i in range(tmpts)]

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)

        subj_rows['RID'] = rid
        abeta_long = groupLongPivot(subj_rows,'RID','ABETA','CSF_ABETA.')
        tau_long = groupLongPivot(subj_rows,'RID','TAU','CSF_TAU.')
        ptau_long = groupLongPivot(subj_rows,'RID','PTAU','CSF_PTAU.')
        date_long = groupLongPivot(subj_rows,'RID','EXAMDATE','CSF_postAV45.')
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)
        all_df = pd.concat((abeta_long,tau_long,ptau_long,date_long),axis=1)

        closest_ptau = groupClosest(subj_rows, 'EXAMDATE', 'PTAU', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_ptau_bin = groupClosest(subj_rows, 'EXAMDATE', 'PTAU_BIN', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_abeta = groupClosest(subj_rows, 'EXAMDATE', 'ABETA', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_abeta_bin = groupClosest(subj_rows, 'EXAMDATE', 'ABETA_BIN', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_tau = groupClosest(subj_rows, 'EXAMDATE', 'TAU', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_tau_bin = groupClosest(subj_rows, 'EXAMDATE', 'TAU_BIN', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        all_df['CSF_PTAU_closest_AV45_1'], all_df['CSF_PTAU_closest_AV45_2'], all_df['CSF_PTAU_closest_AV45_3'] = tuple(closest_ptau)
        all_df['CSF_PTAU_closest_AV45_1_BIN_23'], all_df['CSF_PTAU_closest_AV45_2_BIN_23'], all_df['CSF_PTAU_closest_AV45_3_BIN_23'] = tuple(closest_ptau_bin)
        all_df['CSF_ABETA_closest_AV45_1'], all_df['CSF_ABETA_closest_AV45_2'], all_df['CSF_ABETA_closest_AV45_3'] = tuple(closest_abeta)
        all_df['CSF_ABETA_closest_AV45_1_BIN_192'], all_df['CSF_ABETA_closest_AV45_2_BIN_192'], all_df['CSF_ABETA_closest_AV45_3_BIN_192'] = tuple(closest_abeta_bin)
        all_df['CSF_TAU_closest_AV45_1'], all_df['CSF_TAU_closest_AV45_2'], all_df['CSF_TAU_closest_AV45_3'] = tuple(closest_tau)
        all_df['CSF_TAU_closest_AV45_1_BIN_93'], all_df['CSF_TAU_closest_AV45_2_BIN_93'], all_df['CSF_TAU_closest_AV45_3_BIN_93'] = tuple(closest_tau_bin)

        closest_ptau = groupClosest(subj_rows, 'EXAMDATE', 'PTAU', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        closest_ptau_bin = groupClosest(subj_rows, 'EXAMDATE', 'PTAU_BIN', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        closest_abeta = groupClosest(subj_rows, 'EXAMDATE', 'ABETA', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        closest_abeta_bin = groupClosest(subj_rows, 'EXAMDATE', 'ABETA_BIN', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        closest_tau = groupClosest(subj_rows, 'EXAMDATE', 'TAU', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        closest_tau_bin = groupClosest(subj_rows, 'EXAMDATE', 'TAU_BIN', [av1451_date1,av1451_date2,av1451_date3],day_limit=365)
        all_df['CSF_PTAU_closest_AV1451_1'], all_df['CSF_PTAU_closest_AV1451_2'], all_df['CSF_PTAU_closest_AV1451_3'] = tuple(closest_ptau)
        all_df['CSF_PTAU_closest_AV1451_1_BIN_23'], all_df['CSF_PTAU_closest_AV1451_2_BIN_23'], all_df['CSF_PTAU_closest_AV1451_3_BIN_23'] = tuple(closest_ptau_bin)
        all_df['CSF_ABETA_closest_AV1451_1'], all_df['CSF_ABETA_closest_AV1451_2'], all_df['CSF_ABETA_closest_AV1451_3'] = tuple(closest_abeta)
        all_df['CSF_ABETA_closest_AV1451_1_BIN_192'], all_df['CSF_ABETA_closest_AV1451_2_BIN_192'], all_df['CSF_ABETA_closest_AV1451_3_BIN_192'] = tuple(closest_abeta_bin)
        all_df['CSF_TAU_closest_AV1451_1'], all_df['CSF_TAU_closest_AV1451_2'], all_df['CSF_TAU_closest_AV1451_3'] = tuple(closest_tau)
        all_df['CSF_TAU_closest_AV1451_1_BIN_93'], all_df['CSF_TAU_closest_AV1451_2_BIN_93'], all_df['CSF_TAU_closest_AV1451_3_BIN_93'] = tuple(closest_tau_bin)

        # get slopes
        ptau_slope = df_slope(all_df, date_long.columns, ptau_long.columns, take_diff=False, exact=False)
        tau_slope = df_slope(all_df, date_long.columns, tau_long.columns, take_diff=False, exact=False)
        abeta_slope = df_slope(all_df, date_long.columns, abeta_long.columns, take_diff=False, exact=False)
        all_df['CSF_PTAU_slope'] = ptau_slope[rid]
        all_df['CSF_ABETA_slope'] = abeta_slope[rid]
        all_df['CSF_TAU_slope'] = tau_slope[rid]

        return all_df

    parsed_df = parseSubjectGroups(csf_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df


def syncUWData(master_df, uw_file, registry):
    uw_df = importUW(uw_file, registry=registry, as_df=True)
    tmpts = max(Counter(uw_df.index).values())

    headers = ['UW_MEM_%s' % (i+1) for i in range(tmpts)]
    headers += ['UW_EF_%s' % (i+1) for i in range(tmpts)]
    headers += ['UW_postAV45_%s' % (i+1) for i in range(tmpts)]
    headers += ['UW_MEM_postAV45_count', 'UW_MEM_slope',
                'UW_MEM_AV45_1','UW_MEM_AV45_2','UW_MEM_AV45_3',
                'UW_MEM_AV1451_1','UW_MEM_AV1451_2','UW_MEM_AV1451_3']
    headers += ['UW_EF_postAV45_count', 'UW_EF_slope',
                'UW_EF_AV45_1','UW_EF_AV45_2','UW_EF_AV45_3',
                'UW_EF_AV1451_1','UW_EF_AV1451_2','UW_EF_AV1451_3']

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)

        # get long measurements
        subj_rows['RID'] = rid
        mem_long = groupLongPivot(subj_rows, 'RID','ADNI_MEM','UW_MEM_')
        ef_long = groupLongPivot(subj_rows, 'RID', 'ADNI_EF','UW_EF_')
        date_long = groupLongPivot(subj_rows, 'RID','EXAMDATE','UW_postAV45_')
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)
        all_df = pd.concat((mem_long,ef_long,date_long),axis=1)

        # get slopes
        mem_slope = df_slope(all_df, date_long.columns, mem_long.columns, take_diff=False, exact=False)
        ef_slope = df_slope(all_df, date_long.columns, ef_long.columns, take_diff=False, exact=False)
        all_df['UW_MEM_slope'] = mem_slope[rid]
        all_df['UW_EF_slope'] = ef_slope[rid]

        # get counts
        all_df['UW_MEM_postAV45_count'] = df_count(all_df, date_long.columns, mem_long.columns)[rid]
        all_df['UW_EF_postAV45_count'] = df_count(all_df, date_long.columns, ef_long.columns)[rid]

        # get closest to av45
        closest_mem = groupClosest(subj_rows, 'EXAMDATE', 'ADNI_MEM', [av45_date1, av45_date2, av45_date3],day_limit=365/2)
        closest_ef = groupClosest(subj_rows, 'EXAMDATE', 'ADNI_EF', [av45_date1, av45_date2, av45_date3],day_limit=365/2)
        all_df['UW_MEM_AV45_1'] = closest_mem[0]
        all_df['UW_MEM_AV45_2'] = closest_mem[1]
        all_df['UW_MEM_AV45_3'] = closest_mem[2]
        all_df['UW_EF_AV45_1'] = closest_ef[0]
        all_df['UW_EF_AV45_2'] = closest_ef[1]
        all_df['UW_EF_AV45_3'] = closest_ef[2]

        # get closest to av1451
        closest_mem = groupClosest(subj_rows, 'EXAMDATE', 'ADNI_MEM', [av1451_date1, av1451_date2, av1451_date3],day_limit=365)
        closest_ef = groupClosest(subj_rows, 'EXAMDATE', 'ADNI_EF', [av1451_date1, av1451_date2, av1451_date3],day_limit=365)
        all_df['UW_MEM_AV1451_1'] = closest_mem[0]
        all_df['UW_MEM_AV1451_2'] = closest_mem[1]
        all_df['UW_MEM_AV1451_3'] = closest_mem[2]
        all_df['UW_EF_AV1451_1'] = closest_ef[0]
        all_df['UW_EF_AV1451_2'] = closest_ef[1]
        all_df['UW_EF_AV1451_3'] = closest_ef[2]

        return all_df

    parsed_df = parseSubjectGroups(uw_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df

def syncUCBFreesurferData(master_df, ucb_fs_volumes):
    fs_df = importFSVolumes(ucb_fs_volumes, as_df=True)
    tmpts = max(Counter(fs_df.index).values())

    headers = ['UCB_FS_HC/ICV_%s' % (i+1) for i in range(tmpts)]
    headers += ['UCB_FS_postAV45_%s' % (i+1) for i in range(tmpts)]
    headers += ['UCB_FS_postAV45_count', 'UCB_FS_postAV45_interval',
                'UCB_FS_HC/ICV_slope', 'UCB_FS_HC/ICVavg_slope',
                'UCB_FS_HC/ICV_AV45_1','UCB_FS_HC/ICV_AV45_1_DATE',
                'UCB_FS_HC/ICV_AV45_2','UCB_FS_HC/ICV_AV45_2_DATE',
                'UCB_FS_HC/ICV_AV45_3','UCB_FS_HC/ICV_AV45_3_DATE',
                'UCB_FS_HC/ICV_AV1451_1','UCB_FS_HC/ICV_AV1451_1_DATE',
                'UCB_FS_HC/ICV_AV1451_2','UCB_FS_HC/ICV_AV1451_2_DATE',
                'UCB_FS_HC/ICV_AV1451_3','UCB_FS_HC/ICV_AV1451_3_DATE']

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)

        # Calculate HC/ICV
        bl_icv = subj_rows.iloc[0]['ICV']
        avg_icv = np.mean(subj_rows['ICV'])
        subj_rows['HC/ICV_BL'] = subj_rows['HCV'] / bl_icv
        subj_rows['HC/ICV_AVG'] = subj_rows['HCV'] / avg_icv

        # Get long measurements
        subj_rows['RID'] = rid
        hcicv_bl_long = groupLongPivot(subj_rows, 'RID','HC/ICV_BL','UCB_FS_HC/ICV_')
        hcicv_avg_long = groupLongPivot(subj_rows, 'RID','HC/ICV_AVG','UCB_FS_HC/ICV_')
        date_long = groupLongPivot(subj_rows, 'RID','EXAMDATE','UCB_FS_postAV45_')
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)
        all_df = pd.concat((hcicv_bl_long,date_long),axis=1)

        # Get closest av45 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'HC/ICV_BL', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        all_df['UCB_FS_HC/ICV_AV45_1'] = closest_vals[0]
        all_df['UCB_FS_HC/ICV_AV45_1_DATE'] = closest_dates[0]
        all_df['UCB_FS_HC/ICV_AV45_2'] = closest_vals[1]
        all_df['UCB_FS_HC/ICV_AV45_2_DATE'] = closest_dates[1]
        all_df['UCB_FS_HC/ICV_AV45_3'] = closest_vals[2]
        all_df['UCB_FS_HC/ICV_AV45_3_DATE'] = closest_dates[2]
        all_df['UCB_FS_HC/ICV_AV45_1_DATE'] = pd.to_datetime(all_df['UCB_FS_HC/ICV_AV45_1_DATE'])
        all_df['UCB_FS_HC/ICV_AV45_2_DATE'] = pd.to_datetime(all_df['UCB_FS_HC/ICV_AV45_2_DATE'])
        all_df['UCB_FS_HC/ICV_AV45_3_DATE'] = pd.to_datetime(all_df['UCB_FS_HC/ICV_AV45_3_DATE'])

        # Get closest AV1451 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'HC/ICV_BL', [av1451_date1,av1451_date2,av1451_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av1451_date1,av1451_date2,av1451_date3],day_limit=365/2)
        all_df['UCB_FS_HC/ICV_AV1451_1'] = closest_vals[0]
        all_df['UCB_FS_HC/ICV_AV1451_1_DATE'] = closest_dates[0]
        all_df['UCB_FS_HC/ICV_AV1451_2'] = closest_vals[1]
        all_df['UCB_FS_HC/ICV_AV1451_2_DATE'] = closest_dates[1]
        all_df['UCB_FS_HC/ICV_AV1451_3'] = closest_vals[2]
        all_df['UCB_FS_HC/ICV_AV1451_3_DATE'] = closest_dates[2]
        all_df['UCB_FS_HC/ICV_AV1451_1_DATE'] = pd.to_datetime(all_df['UCB_FS_HC/ICV_AV1451_1_DATE'])
        all_df['UCB_FS_HC/ICV_AV1451_2_DATE'] = pd.to_datetime(all_df['UCB_FS_HC/ICV_AV1451_2_DATE'])
        all_df['UCB_FS_HC/ICV_AV1451_3_DATE'] = pd.to_datetime(all_df['UCB_FS_HC/ICV_AV1451_3_DATE'])

        # Get slopes
        hcicv_bl_slope = df_slope(all_df, date_long.columns, hcicv_bl_long.columns, take_diff=False, exact=False)
        hcicv_avg_slope = df_slope(all_df, date_long.columns, hcicv_avg_long.columns, take_diff=False, exact=False)
        all_df['UCB_FS_HC/ICV_slope'] = hcicv_bl_slope[rid]
        all_df['UCB_FS_HC/ICVavg_slope'] = hcicv_avg_slope[rid]

        # Get counts + interval
        all_df['UCB_FS_postAV45_count'] = df_count(all_df,date_long.columns,hcicv_bl_long.columns)[rid]
        last_date = subj_rows.iloc[-1]['EXAMDATE']
        all_df['UCB_FS_postAV45_interval'] = (last_date - av45_date1).days/365.25 if (not isnan(av45_date1) and last_date > av45_date1) else np.nan

        return all_df

    parsed_df = parseSubjectGroups(fs_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df


def syncUCSFFreesurferData(master_df, ucsf_files, mprage_file, prefix):
    fs_df_list = [importUCSFFreesurfer(filename, mprage_file, version=fsversion, include_failed=False, as_df=True) for fsversion,filename in ucsf_files]
    fs_df = pd.concat(fs_df_list,axis=0)
    fs_df = fs_df[['EXAMDATE','FLDSTRENG','version','ICV','HCV']]

    # average same day scans if matching in field strength + fs version
    # then, of same day results, take first row
    fs_df = fs_df.reset_index().groupby(['RID','FLDSTRENG','version','EXAMDATE']).aggregate(np.mean).reset_index()
    fs_df = fs_df.groupby(['RID','EXAMDATE']).first().reset_index()
    fs_df.set_index('RID',inplace=True)
    tmpts = max(Counter(fs_df.index).values())

    headers = ['%s_HC/ICV_%s' % (prefix,i+1) for i in range(tmpts)]
    headers += ['%s_postAV45_%s' % (prefix,i+1) for i in range(tmpts)]
    headers += ['%s_postAV45_count' % prefix, '%s_postAV45_interval' % prefix,
                '%s_HC/ICV_slope' % prefix, '%s_HC/ICVavg_slope' % prefix,
                '%s_HC/ICV_AV45_1' % prefix,'%s_HC/ICV_AV45_1_DATE' % prefix,
                '%s_HC/ICV_AV45_2' % prefix,'%s_HC/ICV_AV45_2_DATE' % prefix,
                '%s_HC/ICV_AV45_3' % prefix,'%s_HC/ICV_AV45_3_DATE' % prefix,
                '%s_HC/ICV_AV1451_1' % prefix,'%s_HC/ICV_AV1451_1_DATE' % prefix,
                '%s_HC/ICV_AV1451_2' % prefix,'%s_HC/ICV_AV1451_2_DATE' % prefix,
                '%s_HC/ICV_AV1451_3' % prefix,'%s_HC/ICV_AV1451_3_DATE' % prefix,]
    headers += ['%s_MRI_STRENGTH_%s' % (prefix,i+1) for i in range(tmpts)]
    headers += ['%s_FSVERSION_%s' % (prefix,i+1) for i in range(tmpts)]

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)

        # Filter: only use 1.5T scans if no 3T scans
        scans_15 = subj_rows[subj_rows.FLDSTRENG == 1.5]
        scans_3 = subj_rows[subj_rows.FLDSTRENG == 3.0]
        subj_rows = scans_15 if len(scans_3.index) == 0 else scans_3
        if len(subj_rows) == 0:
            return pd.DataFrame()

        # Calculate HC/ICV
        bl_icv = subj_rows.iloc[0]['ICV']
        avg_icv = np.mean(subj_rows['ICV'])
        subj_rows['HC/ICV_BL'] = subj_rows['HCV'] / bl_icv
        subj_rows['HC/ICV_AVG'] = subj_rows['HCV'] / avg_icv

        # Get long measurements
        subj_rows['RID'] = rid
        hcicv_bl_long = groupLongPivot(subj_rows, 'RID','HC/ICV_BL','%s_HC/ICV_' % prefix)
        hcicv_avg_long = groupLongPivot(subj_rows, 'RID','HC/ICV_AVG','%s_HC/ICV_' % prefix)
        strength_long = groupLongPivot(subj_rows, 'RID','FLDSTRENG','%s_MRI_STRENGTH_' % prefix)
        version_long = groupLongPivot(subj_rows, 'RID','version','%s_FSVERSION_' % prefix)
        date_long = groupLongPivot(subj_rows, 'RID','EXAMDATE','%s_postAV45_' % prefix)
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)
        all_df = pd.concat((hcicv_bl_long,date_long,strength_long,version_long),axis=1)

        # Get closest av45 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'HC/ICV_BL', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        all_df['%s_HC/ICV_AV45_1' % prefix] = closest_vals[0]
        all_df['%s_HC/ICV_AV45_1_DATE' % prefix] = closest_dates[0]
        all_df['%s_HC/ICV_AV45_2' % prefix] = closest_vals[1]
        all_df['%s_HC/ICV_AV45_2_DATE' % prefix] = closest_dates[1]
        all_df['%s_HC/ICV_AV45_3' % prefix] = closest_vals[2]
        all_df['%s_HC/ICV_AV45_3_DATE' % prefix] = closest_dates[2]
        all_df['%s_HC/ICV_AV45_1_DATE' % prefix] = pd.to_datetime(all_df['%s_HC/ICV_AV45_1_DATE' % prefix])
        all_df['%s_HC/ICV_AV45_2_DATE' % prefix] = pd.to_datetime(all_df['%s_HC/ICV_AV45_2_DATE' % prefix])
        all_df['%s_HC/ICV_AV45_3_DATE' % prefix] = pd.to_datetime(all_df['%s_HC/ICV_AV45_3_DATE' % prefix])

        # Get closest AV1451 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'HC/ICV_BL', [av1451_date1,av1451_date2,av1451_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av1451_date1,av1451_date2,av1451_date3],day_limit=365/2)
        all_df['%s_HC/ICV_AV1451_1' % prefix] = closest_vals[0]
        all_df['%s_HC/ICV_AV1451_1_DATE' % prefix] = closest_dates[0]
        all_df['%s_HC/ICV_AV1451_2' % prefix] = closest_vals[1]
        all_df['%s_HC/ICV_AV1451_2_DATE' % prefix] = closest_dates[1]
        all_df['%s_HC/ICV_AV1451_3' % prefix] = closest_vals[2]
        all_df['%s_HC/ICV_AV1451_3_DATE' % prefix] = closest_dates[2]
        all_df['%s_HC/ICV_AV1451_1_DATE' % prefix] = pd.to_datetime(all_df['%s_HC/ICV_AV1451_1_DATE' % prefix])
        all_df['%s_HC/ICV_AV1451_2_DATE' % prefix] = pd.to_datetime(all_df['%s_HC/ICV_AV1451_2_DATE' % prefix])
        all_df['%s_HC/ICV_AV1451_3_DATE' % prefix] = pd.to_datetime(all_df['%s_HC/ICV_AV1451_3_DATE' % prefix])

        # Get slopes
        hcicv_bl_slope = df_slope(all_df, date_long.columns, hcicv_bl_long.columns, take_diff=False, exact=False)
        hcicv_avg_slope = df_slope(all_df, date_long.columns, hcicv_avg_long.columns, take_diff=False, exact=False)
        all_df['%s_HC/ICV_slope' % prefix] = hcicv_bl_slope[rid]
        all_df['%s_HC/ICVavg_slope' % prefix] = hcicv_avg_slope[rid]

        # Get counts + interval
        all_df['%s_postAV45_count' % prefix] = df_count(all_df,date_long.columns,hcicv_bl_long.columns)[rid]
        last_date = subj_rows.iloc[-1]['EXAMDATE']
        all_df['%s_postAV45_interval' % prefix] = (last_date - av45_date1).days/365.25 if (not isnan(av45_date1) and last_date > av45_date1) else np.nan

        return all_df

    parsed_df = parseSubjectGroups(fs_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df

def syncWMHData(master_df, wmh_file):
    wmh_df = importWMH(wmh_file, as_df=True)
    tmpts = max(Counter(wmh_df.index).values())

    headers = ['WMH_percentOfICV.%s' % (i+1) for i in range(tmpts)]
    headers += ['WMH_WHITMATHYP.%s' % (i+1) for i in range(tmpts)]
    headers += ['WMH_postAV45.%s' % (i+1) for i in range(tmpts)]
    headers += ['WMH_percentOfICV_AV45_1', 'WMH_percentOfICV_AV45_1_DATE',
                'WMH_percentOfICV_AV45_2', 'WMH_percentOfICV_AV45_2_DATE',
                'WMH_percentOfICV_AV45_3', 'WMH_percentOfICV_AV45_3_DATE',
                'WMH_percentOfICV_slope']

    def extraction_fn(rid, subj_rows):
        subj_rows.sort_values('EXAMDATE',inplace=True)
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)

        # Get long measurements
        subj_rows['RID'] = rid
        wmh_percent_long = groupLongPivot(subj_rows,'RID','wmh_percent','WMH_percentOfICV.')
        wmh_long = groupLongPivot(subj_rows,'RID','wmh','WMH_WHITMATHYP.')
        date_long = groupLongPivot(subj_rows,'RID','EXAMDATE','WMH_postAV45.')
        date_long = date_long.applymap(lambda x: (x-av45_date1).days/365.25 if not isnan(av45_date1) else np.nan)
        date_long = date_long.applymap(lambda x: x if (not isnan(x) and x >= -90/365.0) else np.nan)
        all_df = pd.concat((wmh_percent_long,wmh_long,date_long),axis=1)

        gd_slope = df_slope(all_df, date_long.columns, wmh_percent_long.columns, take_diff=False, exact=False)
        all_df['WMH_percentOfICV_slope'] = gd_slope[rid]

        # Get closest av45 measurements
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'wmh_percent', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        closest_dates = groupClosest(subj_rows, 'EXAMDATE', 'EXAMDATE', [av45_date1,av45_date2,av45_date3],day_limit=365/2)
        all_df['WMH_percentOfICV_AV45_1'] = closest_vals[0]
        all_df['WMH_percentOfICV_AV45_1_DATE'] = closest_dates[0]
        all_df['WMH_percentOfICV_AV45_2'] = closest_vals[1]
        all_df['WMH_percentOfICV_AV45_2_DATE'] = closest_dates[1]
        all_df['WMH_percentOfICV_AV45_3'] = closest_vals[2]
        all_df['WMH_percentOfICV_AV45_3_DATE'] = closest_dates[2]
        all_df['WMH_percentOfICV_AV45_1_DATE'] = pd.to_datetime(all_df['WMH_percentOfICV_AV45_1_DATE'])
        all_df['WMH_percentOfICV_AV45_2_DATE'] = pd.to_datetime(all_df['WMH_percentOfICV_AV45_2_DATE'])
        all_df['WMH_percentOfICV_AV45_3_DATE'] = pd.to_datetime(all_df['WMH_percentOfICV_AV45_3_DATE'])

        return all_df

    parsed_df = parseSubjectGroups(wmh_df, extraction_fn)
    master_df = updateDataFrame(master_df, parsed_df, headers=headers, after=None)
    return master_df

def eliminateColumns(master_df):
    to_remove = ['LastClinicalVisit',
                 'DIFF_LastClinicalVisit_AV45',
                 'ClinicalVisitClosestto_AV45',
                 'ClinicalVisitClosestto_AV45_2',
                 'datediff_MMSE_AV45',
                 'PostAV45Followup',
                 'abeta_closest_AV45',
                 'abeta_bin192',
                 'abeta_closest_AV45_2',
                 'abeta_2_bin192',
                 'abeta_diff',
                 'abeta_slope',
                 'tau_closest_AV45',
                 'tau_bin93',
                 'tau_closest_AV45_2',
                 'ptau_closest_AV45',
                 'ptau_bin23',
                 'ptau_closest_2',
                 'WMH_CLOSEST_AV45',
                 'WMH_CLOSEST_',
                 'WMH_slope',
                 'ClosestCSF_AV45',
                 'diff_CSF_AV45',
                 'AV45CSF<3months',
                 'ClosestCSF_AV45_2',
                 'diff_CSF_AV45_2',
                 'LastCSF',
                 'LastCSF_AV45_1_TimeDiff',
                 'LastCSFAbeta',
                 'ADAS_slope_all',
                 'Hippocampal_Volume_3months',
                 'Left_Volume',
                 'Right_Volume',
                 'ICV',
                 'Hippocampal_Volume_normalizd',
                 'HippVol_norm_BIN0.0047',
                 'ICV_pchange(yrs)',
                 'Closest_DX_Jun15',
                 'DX_Jun15_closestdate',
                 'MCItoADConv(fromav45)',
                 'MCItoADconv_',
                 'BD MM-YY',
                 '2nd_AV45_Jul2014']
    to_remove_prefixes = ['WHITMATHYP.','ABETA.','TAU.','PTAU.',
                          'EXAM_tau.','ADAScog_DATE','MMSE','AVLT_DATE.',
                          'TBMSyn','AV45_','UW_','FSL_','FSX_','UCB_FS',
                          'FDG','AVLT','ADAS','TIMEpostAV45_ADAS','TIMEpostAV45_MMSE',
                          'MMSCORE','CSF','WMH','TIME_ADAS','TIMEreltoAV45_ADAS',
                          'TIME_AVLT','TIMEreltoAV45_AVLT','TIMEpostAV45_AVLT']
    for prefix in to_remove_prefixes:
        to_remove += [_ for _ in master_df.columns if _.startswith(prefix)]

    # drop columns
    to_remove = [_ for _ in to_remove if _ in master_df.columns]
    master_df.drop(to_remove, axis=1, inplace=True)
    return master_df


def runPipeline():
    master_df = parseCSV(master_file, use_second_line=True, as_df=True)
    master_df.set_index('RID',inplace=True)

    # get references
    registry = importRegistry(registry_file, as_df=True)

    print "\nELIMINATING COLUMNS\n"
    master_df = eliminateColumns(master_df)
    print "\nMANUALLY ADDING SUBJECTS"
    master_df = manualAddOns(master_df, RID_ADDONS)
    print "\nSYNCING AV45 NONTP\n"
    master_df = syncAV45Data(master_df, av45_nontp_file, suffix='NONTP')
    print "\nSYNCING AV45 TP\n"
    master_df = syncAV45Data(master_df, av45_tp_file, suffix='TP')
    print "\nSYNCING DEMOG\n"
    master_df = syncDemogData(master_df, demog_file, pet_meta_file)
    print "\nSYNCING CDR\n"
    master_df = syncCDRData(master_df, cdr_file, registry)
    print "\nSYNCING AV1451 TP\n"
    master_df = syncAV1451Data(master_df, av1451_tp_file)
    print "\nSYNCING AV45 ROUSSET PVC\n"
    master_df = syncAV45RoussetResults(master_df, av45_rousset_csv)
    print "\nSYNCING AV1451 ROUSSET PVC\n"
    master_df = syncAV1451RoussetResults(master_df, av1451_rousset_csv)
    print "\nSYNCING DIAGNOSES\n"
    master_df = syncDiagnosisData(master_df, diagnosis_file, arm_file, registry)
    print "\nSYNCING FAQ\n"
    master_df = syncFAQData(master_df, faq_file, registry)
    print "\nSYNCING NPI\n"
    master_df = syncNPIData(master_df, npi_file)
    print "\nSYNCING MHIST\n"
    master_df = syncMHISTData(master_df, mhist_file)
    print "\nSYNCING UW NEURO\n"
    master_df = syncUWData(master_df, uw_file, registry)
    print "\nSYNCING UCSF LONG FreeSurfer\n"
    master_df = syncUCSFFreesurferData(master_df, ucsf_long_files, mprage_file, 'FSL')
    print "\nSYNCING UCSF CROSS FreeSurfer"
    master_df = syncUCSFFreesurferData(master_df, ucsf_cross_files, mprage_file, 'FSX')
    print "\nSYNCING UCB Freesurfer\n"
    master_df = syncUCBFreesurferData(master_df, ucb_fs_volumes)
    print "\nSYNCING FDG\n"
    master_df = syncFDGData(master_df, fdg_file, registry)
    print "\nSYNCING TBMSYN\n"
    master_df = syncTBMSynData(master_df, tbm_file)
    print "\nSYNCING MMSE\n"
    master_df = syncMMSEData(master_df, mmse_file, registry)
    print "\nSYNCING ADASCOG\n"
    master_df = syncADASCogData(master_df, adni1_adas_file, adnigo2_adas_file, registry)
    print "\nSYNCING AVLT\n"
    master_df = syncAVLTData(master_df, neuro_battery_file, registry)
    print "\nSYNCING CSF\n"
    master_df = syncCSFData(master_df, csf_files, registry)
    print "\nSYNCING GD\n"
    master_df = syncGDData(master_df, gd_file, registry)
    print "\nSYNCING WMH\n"
    master_df = syncWMHData(master_df, wmh_file)

    dumpDFtoCSV(master_df,output_file,decimal_places=3)
    mergeSlopes(output_file)
    addCategories(output_file)

if __name__ == '__main__':
    now = datetime.now()
    pd.options.mode.chained_assignment = None

    # IO files
    master_file = "../FDG_AV45_COGdata/FDG_AV45_COGdata_03_07_16.csv"
    output_file = "../FDG_AV45_COGdata_%s.csv" % (now.strftime("%m_%d_%y"))

    # ADNI docs
    registry_file = "../docs/ADNI/REGISTRY.csv"
    arm_file = "../docs/ADNI/ARM.csv"
    pet_meta_file = "../docs/ADNI/PET_META_LIST.csv"
    demog_file = "../docs/ADNI/PTDEMOG.csv"
    diagnosis_file = '../docs/ADNI/DXSUM_PDXCONV_ADNIALL.csv'
    csf_files = ['../docs/ADNI/UPENNBIOMK5_10_31_13.csv',
                 '../docs/ADNI/UPENNBIOMK6_07_02_13.csv',
                 '../docs/ADNI/UPENNBIOMK7.csv',
                 '../docs/ADNI/UPENNBIOMK8.csv']
    wmh_file = '../docs/ADNI/UCD_ADNI2_WMH_10_26_15.csv'
    fdg_file = '../docs/ADNI/UCBERKELEYFDG_07_30_15.csv'
    mprage_file = '../docs/ADNI/MPRAGEMETA.csv'
    mhist_file = '../docs/ADNI/RECMHIST.csv'
    faq_file = '../docs/ADNI/FAQ.csv'
    npi_file = '../docs/ADNI/NPI.csv'

    # Output files
    av45_tp_file = "../output/06_01_16/UCBERKELEYAV45_06_01_16_regular_tp.csv"
    av45_nontp_file = "../output/06_01_16/UCBERKELEYAV45_06_01_16_regular_nontp.csv"
    av1451_tp_file = "../output/06_01_16/UCBERKELEYAV1451_06_01_16_regular_tp.csv"
    av45_rousset_csv = "../datasets/pvc_adni_av45/aggregions_output.csv"
    av1451_rousset_csv = "../datasets/pvc_adni_av1451/mostregions_output.csv"

    # Cog files
    mmse_file = "../docs/ADNI/MMSE.csv"
    adni1_adas_file = '../docs/ADNI/ADASSCORES.csv'
    adnigo2_adas_file = '../docs/ADNI/ADAS_ADNIGO2.csv'
    neuro_battery_file = '../docs/ADNI/NEUROBAT.csv'
    gd_file = '../docs/ADNI/GDSCALE.csv'
    uw_file = '../docs/ADNI/UWNPSYCHSUM_04_22_16.csv'
    cdr_file = '../docs/ADNI/CDR.csv'

    # MR files
    tbm_file = '../docs/ADNI/MAYOADIRL_MRI_TBMSYN_12_08_15.csv'
    ucsf_long_files = [('4.4','../docs/ADNI/UCSFFSL_02_01_16.csv'),
                       ('5.1','../docs/ADNI/UCSFFSL51Y1_05_04_16.csv')]
    ucsf_cross_files = [('4.3','../docs/ADNI/UCSFFSX_11_02_15.csv'),
                        ('5.1','../docs/ADNI/UCSFFSX51_05_04_16.csv'),
                        ('5.1','../docs/ADNI/UCSFFSX51_ADNI1_3T_02_01_16.csv')]
    ucb_fs_volumes = '../docs/ADNI/adni_av45_fs_volumes_04-05-2016.csv'

    # Run pipeline
    runPipeline()
