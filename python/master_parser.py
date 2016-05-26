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



'''
def getAV45Dates(old_l, patient_pets=None):
    if patient_pets is None:
        patient_pets = []
    # Get AV45 Scan dates
    if old_l.get('AV45_Date','') != '':
        bl_av45 = datetime.strptime(old_l['AV45_Date'], '%m/%d/%y')
    elif len(patient_pets) >= 1:
        bl_av45 = patient_pets[0]
    else:
        bl_av45 = None
    if old_l.get('AV45_2_Date','') != '':
        av45_2 = datetime.strptime(old_l['AV45_2_Date'], '%m/%d/%y')
    elif len(patient_pets) >= 2:
        av45_2 = patient_pets[1]
    else:
        av45_2 = None
    if old_l.get('AV45_3_Date','') != '':
        av45_3 = datetime.strptime(old_l['AV45_3_Date'], '%m/%d/%y')
    elif len(patient_pets) >= 3:
        av45_3 = patient_pets[2]
    else:
        av45_3 = None
    return (bl_av45, av45_2, av45_3)
'''

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

'''
def getAV1451Dates(old_l, patient_pets=None):
    if patient_pets is None:
        patient_pets = []
    # Get AV45 Scan dates
    if old_l.get('AV1451_Date','') != '':
        bl_av45 = datetime.strptime(old_l['AV1451_Date'], '%m/%d/%y')
    elif len(patient_pets) >= 1:
        bl_av45 = patient_pets[0]
    else:
        bl_av45 = None
    if old_l.get('AV1451_2_Date','') != '':
        av45_2 = datetime.strptime(old_l['AV1451_2_Date'], '%m/%d/%y')
    elif len(patient_pets) >= 2:
        av45_2 = patient_pets[1]
    else:
        av45_2 = None
    if old_l.get('AV1451_3_Date','') != '':
        av45_3 = datetime.strptime(old_l['AV1451_3_Date'], '%m/%d/%y')
    elif len(patient_pets) >= 3:
        av45_3 = patient_pets[2]
    else:
        av45_3 = None
    return (bl_av45, av45_2, av45_3)

def wipeKeys(lines, keys):
    empty_dict = {k:'' for k in keys}
    new_lines = []
    for l in lines:
        l.update(empty_dict)
        new_lines.append(l)
    return l
'''

def addCategories(output_file):
    colcats_df = pd.read_csv('column_categories.csv')
    colcats_df.set_index('COLUMN',inplace=True)
    df = pd.read_csv(output_file, low_memory=False)
    cats = []
    for col in df.columns:
        try:
            cat = colcats_df.loc[col.strip(),'CATEGORY']
        except Exception as e:
            cat = 'Unknown'
            print col
        cats.append(cat)
    df.columns = [cats, df.columns]
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

        # get longitudinal measurements
        subj_rows['RID'] = rid
        first_date = subj_rows.iloc[0]['EXAMDATE']
        faq_long = groupLongPivot(subj_rows, 'RID','FAQTOTAL','FAQTOTAL.')
        date_long = groupLongPivot(subj_rows, 'RID','EXAMDATE','FAQTOTAL_timePostAV45.')
        date_long = date_long.applymap(lambda x: (x-first_date).days/365.25)

        all_df = pd.concat((faq_long,date_long),axis=1)
        slope = df_slope(all_df,
                         date_long.columns,
                         faq_long.columns,
                         take_diff=False,
                         exact=False)
        all_df['FAQTOTAL_slope'] = slope[rid]

        # closest to AV45 scans
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
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

        # get longitudinal measurements
        subj_rows['RID'] = rid
        first_date = subj_rows.iloc[0]['EXAMDATE']
        npi_long = groupLongPivot(subj_rows, 'RID','NPITOTAL','NPITOTAL.')
        date_long = groupLongPivot(subj_rows, 'RID','EXAMDATE','NPITOTAL_timePostAV45.')
        date_long = date_long.applymap(lambda x: (x-first_date).days/365.25)

        all_df = pd.concat((npi_long,date_long),axis=1)
        slope = df_slope(all_df,
                         date_long.columns,
                         npi_long.columns,
                         take_diff=False,
                         exact=False)
        all_df['NPITOTAL_slope'] = slope[rid]

        # closest to AV45 scans
        av45_date1, av45_date2, av45_date3 = getAV45Dates(rid, master_df)
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

def syncGDData(old_headers, old_lines, gd_file, registry_file):
    tmpts = 10
    registry = importRegistry(registry_file, include_all=True)
    gd_by_subj = importGD(gd_file, registry=registry)

    to_add_headers = []
    to_add_headers += ['GD_TOTAL.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['GD_timePostAV45.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['GD_AV45_6MTHS', 'GD_AV45_DATE',
                       'GD_AV45_2_6MTHS', 'GD_AV45_2_DATE',
                       'GD_AV45_3_6MTHS', 'GD_AV45_3_DATE',
                       'GD_slope']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        scores = sorted(subj_row, key=lambda x: x['EXAMDATE'])
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        new_subj_data = {}
        
        six_months = 31 * 6
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToScans(scores, bl_av45, av45_2, av45_3, day_limit=six_months)
        if av45_bl_closest:
            new_subj_data['GD_AV45_6MTHS'] = av45_bl_closest['GDTOTAL']
            new_subj_data['GD_AV45_DATE'] = av45_bl_closest['EXAMDATE']
        else:
            new_subj_data['GD_AV45_6MTHS'] = ''
            new_subj_data['GD_AV45_DATE'] = ''
        if av45_2_closest:
            new_subj_data['GD_AV45_2_6MTHS'] = av45_2_closest['GDTOTAL']
            new_subj_data['GD_AV45_2_DATE'] = av45_2_closest['EXAMDATE']
        else:
            new_subj_data['GD_AV45_2_6MTHS'] = ''
            new_subj_data['GD_AV45_2_DATE'] = ''
        if av45_3_closest:
            new_subj_data['GD_AV45_3_6MTHS'] = av45_3_closest['GDTOTAL']
            new_subj_data['GD_AV45_3_DATE'] = av45_3_closest['EXAMDATE']
        else:
            new_subj_data['GD_AV45_3_6MTHS'] = ''
            new_subj_data['GD_AV45_3_DATE'] = ''
        
        post_av45_points = []
        for i in range(tmpts):
            if i < len(scores):
                cur_score = scores[i]
                date = cur_score['EXAMDATE']
                tot_score = cur_score['GDTOTAL']
                try:
                    tot_score = int(tot_score)
                except:
                    tot_score = ''
                if bl_av45 is not None:
                    timediff = (date - bl_av45).days / 365.0
                    if tot_score != '' and timediff > (-90/365.0):
                        post_av45_points.append((timediff, tot_score))
                else:
                    timediff = ''
                new_subj_data['GD_TOTAL.%s' % (i+1)] = tot_score
                new_subj_data['GD_timePostAV45.%s' % (i+1)] = timediff
        # get slope
        gd_slope = slope(post_av45_points)
        new_subj_data['GD_slope'] = gd_slope if gd_slope is not None else ''
        return new_subj_data

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, gd_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)


    return (new_headers, new_lines)


def syncADASCogData(old_headers, old_lines, adni1_adas_file, adnigo2_adas_file, registry_file):
    tmpts=13
    registry = importRegistry(registry_file)
    adas_by_subj = importADASCog(adni1_adas_file, adnigo2_adas_file, registry=registry)

    to_add_headers = []
    to_add_headers += ['ADAScog.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIME_ADAS.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIMEreltoAV45_ADAS.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIMEpostAV45_ADAS.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['ADAS_post_AV45_followuptime','ADASslope_postAV45',
                       'ADAS_3MTH_AV45','ADAS_3MTHS_AV45DATE',
                       'ADAS_AV45_2_3MTHS','ADAS_AV45_2_DATE',
                       'ADAS_AV45_3_3MTHS','ADAS_AV45_3_DATE',
                       'ADAS_AV1451_1','ADAS_AV1451_1_DATE',
                       'ADAS_AV1451_2','ADAS_AV1451_2_DATE',
                       'ADAS_AV1451_3','ADAS_AV1451_3_DATE',]
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        if len(subj_row) > tmpts:
            raise Exception("Increase ADAS timepoints to %s" % len(subj_row))


        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        bl_av1451, av1451_2, av1451_3 = getAV1451Dates(old_l, patient_pets=patient_pets)
        tests = sorted(subj_row, key=lambda x : x['EXAMDATE'])
        
        dates = [t['EXAMDATE'] for t in tests]
        if len(set(dates)) != len(dates):
            print "%s has Dup dates, skipping: %s" % (subj, dates)
            return {}

        first_scan_date = dates[0]
        new_subj_data = {}
        all_slope_points = []
        post_slope_points = []
        max_followup = None

        # Get closest to AV45 scans
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToScans(tests, bl_av45, av45_2, av45_3, day_limit=31*3)
        new_subj_data['ADAS_3MTH_AV45'] = av45_bl_closest['TOTSCORE'] if av45_bl_closest else ''
        new_subj_data['ADAS_3MTHS_AV45DATE'] = av45_bl_closest['EXAMDATE'] if av45_bl_closest else ''
        new_subj_data['ADAS_AV45_2_3MTHS'] = av45_2_closest['TOTSCORE'] if av45_2_closest else ''
        new_subj_data['ADAS_AV45_2_DATE'] = av45_2_closest['EXAMDATE'] if av45_2_closest else ''
        new_subj_data['ADAS_AV45_3_3MTHS'] = av45_3_closest['TOTSCORE'] if av45_3_closest else ''
        new_subj_data['ADAS_AV45_3_DATE'] = av45_3_closest['EXAMDATE'] if av45_3_closest else ''

        # Get closest to AV1451 scans
        av1451_bl_closest, av1451_2_closest, av1451_3_closest = getClosestToScans(tests, bl_av1451, av1451_2, av1451_3)
        new_subj_data['ADAS_AV1451_1'] = av1451_bl_closest['TOTSCORE'] if av1451_bl_closest else ''
        new_subj_data['ADAS_AV1451_1_DATE'] = av1451_bl_closest['EXAMDATE'] if av1451_bl_closest else ''
        new_subj_data['ADAS_AV1451_2'] = av1451_2_closest['TOTSCORE'] if av1451_2_closest else ''
        new_subj_data['ADAS_AV1451_2_DATE'] = av1451_2_closest['EXAMDATE'] if av1451_2_closest else ''
        new_subj_data['ADAS_AV1451_3'] = av1451_3_closest['TOTSCORE'] if av1451_3_closest else ''
        new_subj_data['ADAS_AV1451_3_DATE'] = av1451_3_closest['EXAMDATE'] if av1451_3_closest else ''

        for i in range(len(tests)):
            test_results = tests[i]
            test_date = test_results['EXAMDATE']
            test_score = test_results['TOTSCORE']
            diff_from_first = (test_date-first_scan_date).days/365.0
            all_slope_points.append((diff_from_first, test_score))
            new_subj_data['ADAScog_DATE%s' % (i+1)] = test_date
            new_subj_data['ADAScog.%s' % (i+1)] = test_score
            new_subj_data['TIME_ADAS.%s' % (i+1)] = diff_from_first
            if bl_av45 is not None:
                timediff = (test_date-bl_av45).days / 365.0
                new_subj_data['TIMEreltoAV45_ADAS.%s' % (i+1)] = timediff
                if timediff > -93.0/365.0:
                    max_followup = timediff
                    new_subj_data['TIMEpostAV45_ADAS.%s' % (i+1)] = timediff
                    post_slope_points.append((timediff, test_score))

        new_subj_data['ADAS_post_AV45_followuptime'] = max(max_followup, 0.0) if max_followup is not None else ''
        adas_slope = slope(post_slope_points)
        new_subj_data['ADASslope_postAV45'] = adas_slope if adas_slope is not None else ''
        return new_subj_data

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, adas_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)


    return (new_headers, new_lines)

def syncMMSEData(old_headers, old_lines, mmse_file, registry_file):
    tmpts = 11
    registry = importRegistry(registry_file, include_all=True)
    mmse_by_subj = importMMSE(mmse_file, registry=registry)

    # add these new headers, if not present
    to_add_headers = []
    to_add_headers += ['MMSCORE.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIMEpostAV45_MMSE.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['MMSE_post_AV45_followuptime',
                       'MMSEslope_postAV45',
                       'MMSE_AV45_3MTHS','MMSE_AV45_DATE',
                       'MMSE_AV45_2_3MTHS','MMSE_AV45_2_DATE',
                       'MMSE_AV45_3_3MTHS','MMSE_AV45_3_DATE']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        tests = sorted(subj_row, key=lambda x: x[0])
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)

        new_subj_data = {}
        max_followup = None
        post_av45_points = []
        for i in range(tmpts):
            if i < len(tests):
                test_date, test_results = tests[i]
                test_score = test_results['MMSCORE']
                try:
                    test_score = int(test_score)
                    if test_score < 0:
                        test_score = ''
                except:
                    test_score = ''

                if bl_av45 is not None:
                    timediff = (test_date-bl_av45).days / 365.0
                    if abs(timediff) <= (93.0/365.0) and 'MMSE_AV45_3MTHS' not in new_subj_data:
                        new_subj_data['MMSE_AV45_3MTHS'] = test_score
                        new_subj_data['MMSE_AV45_DATE'] = test_date
                    if timediff > (-93.0/365.0):
                        max_followup = timediff
                        new_subj_data['TIMEpostAV45_MMSE.%s' % (i+1)] = timediff
                        if test_score != '':
                            post_av45_points.append((timediff, test_score))
                if av45_2 is not None:
                    timediff = (test_date-av45_2).days / 365.0
                    if abs(timediff) <= (93.0/365.0) and 'MMSE_AV45_2_3MTHS' not in new_subj_data:
                        new_subj_data['MMSE_AV45_2_3MTHS'] = test_score
                        new_subj_data['MMSE_AV45_2_DATE'] = test_date
                if av45_3 is not None:
                    timediff = (test_date-av45_3).days / 365.0
                    if abs(timediff) <= (93.0/365.0) and 'MMSE_AV45_3_3MTHS' not in new_subj_data:
                        new_subj_data['MMSE_AV45_3_3MTHS'] = test_score
                        new_subj_data['MMSE_AV45_3_DATE'] = test_date

                new_subj_data['MMSE_DATE%s' % (i+1)] = test_date
                new_subj_data['MMSCORE.%s' % (i+1)] = test_score
        new_subj_data['MMSE_post_AV45_followuptime'] = max(max_followup, 0.0) if max_followup is not None else ''

        # get post av45 slope
        mmse_slope = slope(post_av45_points)
        new_subj_data['MMSEslope_postAV45'] = mmse_slope if mmse_slope is not None else ''
        return new_subj_data

    new_lines = []
    for i, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, mmse_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)


    return new_headers, new_lines

def syncAVLTData(old_headers, old_lines, neuro_battery_file, registry_file):
    tmpts = 11
    registry = importRegistry(registry_file, include_all=True)
    avlt_by_subj = importAVLT(neuro_battery_file, registry=registry)

    # no change in headers
    to_add_headers = []
    to_add_headers += ['AVLT.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIME_AVLT.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIMEreltoAV45_AVLT.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['TIMEpostAV45_AVLT.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['AVLT_AV45_1_3MTHS','AVLT_AV45_1_DATE',
                       'AVLT_AV45_2_3MTHS','AVLT_AV45_2_DATE',
                       'AVLT_AV45_3_3MTHS','AVLT_AV45_3_DATE',
                       'AVLT_AV1451_1','AVLT_AV1451_1_DATE',
                       'AVLT_AV1451_2','AVLT_AV1451_2_DATE',
                       'AVLT_AV1451_3','AVLT_AV1451_3_DATE',
                       'AVLT_post_AV45_followuptime',
                       'AVLT_slope_all','AVLTslope_postAV45']

    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        bl_av1451, av1451_2, av1451_3 = getAV1451Dates(old_l, patient_pets=patient_pets)
        tests = sorted(subj_row, key=lambda x : x['EXAMDATE'])
        
        dates = [t['EXAMDATE'] for t in tests]
        if len(set(dates)) != len(dates):
            print "%s has Dup dates, skipping: %s" % (subj, dates)
            return {}
        new_subj_data = {}
        all_slope_points = []
        post_slope_points = []
        max_followup = ''
        first_scan_date = dates[0]

        three_months = 31 * 3
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToScans(tests, bl_av45, av45_2, av45_3, day_limit=three_months)
        new_subj_data['AVLT_AV45_1_3MTHS'] = av45_bl_closest['TOTS'] if av45_bl_closest else ''
        new_subj_data['AVLT_AV45_1_DATE'] = av45_bl_closest['EXAMDATE'] if av45_bl_closest else ''
        new_subj_data['AVLT_AV45_2_3MTHS'] = av45_2_closest['TOTS'] if av45_2_closest else ''
        new_subj_data['AVLT_AV45_2_DATE'] = av45_2_closest['EXAMDATE'] if av45_2_closest else ''
        new_subj_data['AVLT_AV45_3_3MTHS'] = av45_3_closest['TOTS'] if av45_3_closest else ''
        new_subj_data['AVLT_AV45_3_DATE'] = av45_3_closest['EXAMDATE'] if av45_3_closest else ''

        # Get closest to AV1451 scans
        av1451_bl_closest, av1451_2_closest, av1451_3_closest = getClosestToScans(tests, bl_av1451, av1451_2, av1451_3)
        new_subj_data['AVLT_AV1451_1'] = av1451_bl_closest['TOTS'] if av1451_bl_closest else ''
        new_subj_data['AVLT_AV1451_1_DATE'] = av1451_bl_closest['EXAMDATE'] if av1451_bl_closest else ''
        new_subj_data['AVLT_AV1451_2'] = av1451_2_closest['TOTS'] if av1451_2_closest else ''
        new_subj_data['AVLT_AV1451_2_DATE'] = av1451_2_closest['EXAMDATE'] if av1451_2_closest else ''
        new_subj_data['AVLT_AV1451_3'] = av1451_3_closest['TOTS'] if av1451_3_closest else ''
        new_subj_data['AVLT_AV1451_3_DATE'] = av1451_3_closest['EXAMDATE'] if av1451_3_closest else ''

        for i in range(tmpts):
            if i < len(tests):
                test_results = tests[i]
                test_date = test_results['EXAMDATE']
                test_score = test_results['TOTS']
                diff_from_first = (test_date-first_scan_date).days / 365.0
                diff_from_bl = (test_date - bl_av45).days / 365.0 if bl_av45 is not None else ''
                diff_from_scan2 = (test_date - av45_2).days / 365.0 if av45_2 is not None else ''
                if test_score != '':
                    all_slope_points.append((diff_from_first,test_score))
                if diff_from_bl != '':
                    max_followup = diff_from_bl
                    if test_score != '' and diff_from_bl >= (-93.0/365.0):
                        post_slope_points.append((diff_from_first,test_score))
                        new_subj_data['TIMEpostAV45_AVLT.%s' % (i+1)] = diff_from_bl
                new_subj_data['AVLT_DATE.%s' % (i+1)] = test_date
                new_subj_data['AVLT.%s' % (i+1)] = test_score
                new_subj_data['TIME_AVLT.%s' % (i+1)] = diff_from_first
                new_subj_data['TIMEreltoAV45_AVLT.%s' % (i+1)] = diff_from_bl
        new_subj_data['AVLT_post_AV45_followuptime'] = max(max_followup, 0.0) if max_followup is not None else ''

        # get slopes
        all_slope = slope(all_slope_points)
        post_slope = slope(post_slope_points)
        new_subj_data['AVLT_slope_all'] = all_slope if all_slope is not None else ''
        new_subj_data['AVLTslope_postAV45'] = post_slope if post_slope is not None else ''
        return new_subj_data

    new_lines = []
    new_values = 0
    total = 0
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, avlt_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)   
        old_l.update(new_data)
        new_lines.append(old_l)


    return (new_headers, new_lines)

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
               'Diag@AV45_long','Diag@AV45_2_long','Diag@AV45_3_long',
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
        data['Diag@AV45_long'] = closest_vals[0]
        data['Diag@AV45_2_long'] = closest_vals[1]
        data['Diag@AV45_3_long'] = closest_vals[2]
        # closest to AV1451 scans
        av1451_date1, av1451_date2, av1451_date3 = getAV1451Dates(rid, master_df)
        closest_vals = groupClosest(subj_rows, 'EXAMDATE', 'diag', [av1451_date1,av1451_date2,av1451_date3],day_limit=365/2)
        data['Diag@AV1451_long'] = closest_vals[0]
        data['Diag@AV1451_2_long'] = closest_vals[1]
        data['Diag@AV1451_3_long'] = closest_vals[2]
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

'''
def syncDiagnosisData(old_headers, old_lines, diag_file, registry_file, arm_file):
    registry = importRegistry(registry_file)
    diag_by_subj = importADNIDiagnosis(diag_file, registry=registry)
    arm = importARM(arm_file)

    pivot_date = datetime(day=1, month=2, year=2016)
    pivot_date_closest_diag_key = 'Closest_DX_Feb16'
    pivot_date_closest_date_key = 'DX_Feb16_closestdate'

    to_add_headers = ['MCItoADConv','MCItoADConvDate','AV45_MCItoAD_ConvTime','Baseline_MCItoAD_ConvTime',
                      'Diag@AV45_long','Diag@AV45_2_long','Diag@AV45_3_long',
                      'Diag@AV1451','Diag@AV1451_2','Diag@AV1451_3',
                      'FollowupTimetoDX','Baseline','Init_Diagnosis',
                      pivot_date_closest_diag_key,pivot_date_closest_date_key]
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after='AV45_1_3_Diff')
    wipeKeys(old_lines, to_add_headers)

    # find initial diags per subject
    for subj in diag_by_subj.keys():
        init_diag = None
        if subj in arm:
            sorted_arm = sorted(arm[subj], key=lambda x: x['USERDATE'])
            init_diags = list(set([_['STATUS'] for _ in sorted_arm]))
            init_diag = init_diags[-1]
        for point in diag_by_subj[subj]:
            point['init_diag'] = init_diag

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        subj_diags = sorted(subj_row, key=lambda x: x['EXAMDATE'])

        init_diag = subj_diags[0]['init_diag']
        if init_diag is None:
            raise Exception("NO INITIAL DIAGNOSIS FOUND FOR %s" % subj)

        # get mci type
        mci_type = 'LMCI'
        if init_diag in set(['LMCI', 'EMCI', 'SMC']):
            mci_type = init_diag

        # convert codes to diagnoses base on init diag
        for idx in range(len(subj_diags)):
            if idx == 0:
                subj_diags[idx]['diag'] = init_diag
                continue
            change = subj_diags[idx]['change']
            if change in set([1,2,3]): # stable
                subj_diags[idx]['diag'] = subj_diags[idx-1]['diag']
            elif change in set([5,6]): # conversion to AD
                subj_diags[idx]['diag'] = 'AD'
            elif change in set([4]): # conversion to MCI
                subj_diags[idx]['diag'] = mci_type
            elif change in set([8]): # reversion to MCI
                subj_diags[idx]['diag'] = mci_type
            elif change in set([7,9]): # reversion to normal
                subj_diags[idx]['diag'] = 'N'
        
        # find very first visit date in registry (excluding scmri)
        subj_registry = sorted([_ for _ in registry[subj] if _['VISCODE'] != 'scmri' and _['VISCODE2'] != 'scmri'], key=lambda x: x['EXAMDATE'])
        baseline_registry = subj_registry[0]

        # find pivot data point
        sorted_by_pivot = sorted(subj_diags, key=lambda x: abs(x['EXAMDATE'] - pivot_date))
        closest_to_pivot = sorted_by_pivot[0]

        # get closest to av45
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=None)
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToScans(subj_diags, bl_av45, av45_2, av45_3)

        bl_av1451, av1451_2, av1451_3 = getAV1451Dates(old_l, patient_pets=None)
        av1451_bl_closest, av1451_2_closest, av1451_3_closest = getClosestToScans(subj_diags, bl_av1451, av1451_2, av1451_3)

        # find MCI to AD conversion
        mci_to_ad = None
        for diag_row in subj_diags:
            if diag_row['change'] == 5:
                mci_to_ad = diag_row
                break

        new_data = {'MCItoADConv': 1 if mci_to_ad else 0,
                    'MCItoADConvDate': mci_to_ad['EXAMDATE'] if mci_to_ad else '',
                    'AV45_MCItoAD_ConvTime': (mci_to_ad['EXAMDATE'] - bl_av45).days / 365.0 if mci_to_ad is not None and bl_av45 is not None else '',
                    'Baseline_MCItoAD_ConvTime': (mci_to_ad['EXAMDATE'] - baseline_registry['EXAMDATE']).days / 365.0 if mci_to_ad else '',
                    'Diag@AV45_long': av45_bl_closest['diag'] if av45_bl_closest else '',
                    'Diag@AV45_2_long': av45_2_closest['diag'] if av45_2_closest else '',
                    'Diag@AV45_3_long': av45_3_closest['diag'] if av45_3_closest else '',
                    'Diag@AV1451': av1451_bl_closest['diag'] if av1451_bl_closest else '',
                    'Diag@AV1451_2': av1451_2_closest['diag'] if av1451_2_closest else '',
                    'Diag@AV1451_3': av1451_3_closest['diag'] if av1451_3_closest else '',
                    'FollowupTimetoDX': (pivot_date - baseline_registry['EXAMDATE']).days / 365.0,
                    'Baseline': baseline_registry['EXAMDATE'],
                    'Init_Diagnosis': init_diag,
                    pivot_date_closest_diag_key: closest_to_pivot['diag'],
                    pivot_date_closest_date_key: closest_to_pivot['EXAMDATE']}

        return new_data

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, diag_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    return (new_headers, new_lines)
'''

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

'''
def syncAV1451Data(old_headers, old_lines, av1451_file, pet_meta_file):
    av1451_by_subj = importAV1451(av1451_file)
    pet_meta = importPetMETA(pet_meta_file, tracer='AV1451')

    valid_timepoints = ['BL', 'Scan2', 'Scan3']
    to_add_headers = ['AV1451_Date','AV1451_2_Date','AV1451_3_Date']
    to_add_headers += ['AV1451_Braak12_CerebGray_%s' % tp for tp in valid_timepoints]
    to_add_headers += ['AV1451_Braak34_CerebGray_%s' % tp for tp in valid_timepoints]
    to_add_headers += ['AV1451_Braak56_CerebGray_%s' % tp for tp in valid_timepoints]
    after = old_headers[max(i for i,_ in enumerate(old_headers) if _.startswith('AV45_') and 'PVC' not in _)] # last element that contains 'AV45'
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=after)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        subj_row = sorted(subj_row, key=lambda x: x['EXAMDATE'])
        new_data = {}
        for i, av1451_row in enumerate(subj_row):
            if i == 0:
                tp = 'BL'
                date_key = 'AV1451_Date'
            elif i == 1:
                tp = 'Scan2'
                date_key = 'AV1451_2_Date'
            elif i == 2:
                tp = 'Scan3'
                date_key = 'AV1451_3_Date'
            else:
                raise Exception("Raise # of AV1451 timepoints")
            
            cerebg_left = (av1451_row['LEFT_CEREBELLUM_CORTEX_SIZE'],av1451_row['LEFT_CEREBELLUM_CORTEX'])
            cerebg_right = (av1451_row['RIGHT_CEREBELLUM_CORTEX_SIZE'],av1451_row['RIGHT_CEREBELLUM_CORTEX'])
            braak1 = (av1451_row['BRAAK1_SIZE'],av1451_row['BRAAK1'])
            braak2 = (av1451_row['BRAAK2_SIZE'],av1451_row['BRAAK2'])
            braak3 = (av1451_row['BRAAK3_SIZE'],av1451_row['BRAAK3'])
            braak4 = (av1451_row['BRAAK4_SIZE'],av1451_row['BRAAK4'])
            braak5 = (av1451_row['BRAAK5_SIZE'],av1451_row['BRAAK5'])
            braak6 = (av1451_row['BRAAK6_SIZE'],av1451_row['BRAAK6'])

            cerebg = weightedMean([cerebg_left,cerebg_right])
            braak12 = weightedMean([braak1,braak2]) / cerebg
            braak34 = weightedMean([braak3,braak4]) / cerebg
            braak56 = weightedMean([braak5,braak6]) / cerebg

            new_data[date_key] = patient_pets[i]
            new_data['AV1451_Braak12_CerebGray_%s' % tp] = braak12
            new_data['AV1451_Braak34_CerebGray_%s' % tp] = braak34
            new_data['AV1451_Braak56_CerebGray_%s' % tp] = braak56
        return new_data

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, av1451_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=pet_meta)
        old_l.update(new_data)
        new_lines.append(old_l)

    return (new_headers, new_lines)
'''
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

def syncFDGData(old_headers, old_lines, fdg_file, registry_file):
    fdg_by_subj = importFDG(fdg_file)

    to_add_headers = []
    to_add_headers += ['FDG_pons_Vis%s' % (i+1) for i in range(10)]
    to_add_headers += ['FDGVis%s_ReltoAV45' % (i+1) for i in range(10)]
    to_add_headers += ['FDGVis%s_Date' % (i+1) for i in range(10)]
    to_add_headers += ['FDG_postAV45_slope','FDG_postAV45_followuptime']
    to_add_headers += ['FDG_PONS_AV45_6MTHS', 'FDG_PONS_AV45_DATE',
                       'FDG_PONS_AV45_2_6MTHS', 'FDG_PONS_AV45_2_DATE',
                       'FDG_PONS_AV45_3_6MTHS', 'FDG_PONS_AV45_3_DATE',
                       'FDG_PONS_AV1451_1','FDG_PONS_AV1451_1_DATE',
                       'FDG_PONS_AV1451_2','FDG_PONS_AV1451_2_DATE',
                       'FDG_PONS_AV1451_3','FDG_PONS_AV1451_3_DATE',]
    new_headers= rearrangeHeaders(old_headers, to_add_headers, after='FDG_Bin_Baseline')
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        # collapse subj rows per visit
        av45_bl, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        bl_av1451, av1451_2, av1451_3 = getAV1451Dates(old_l, patient_pets=patient_pets)
        new_data = {}

        examdates = list(set([_['EXAMDATE'] for _ in subj_row]))
        examdates = sorted(examdates)
        subj_fdg = []
        for ed in examdates:
            rel_rows = [_ for _ in subj_row if _['EXAMDATE'] == ed]
            means = [_['MEAN'] for _ in rel_rows]
            pons_avg = np.mean(means)
            subj_fdg.append({'EXAMDATE': ed, 'PONS': pons_avg})

        # get closest to av45s
        six_months = 31 * 6
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToScans(subj_fdg, av45_bl, av45_2, av45_3, day_limit=six_months)
        new_data['FDG_PONS_AV45_6MTHS'] = av45_bl_closest['PONS'] if av45_bl_closest else ''
        new_data['FDG_PONS_AV45_DATE'] = av45_bl_closest['EXAMDATE'] if av45_bl_closest else ''
        new_data['FDG_PONS_AV45_2_6MTHS'] = av45_2_closest['PONS'] if av45_2_closest else ''
        new_data['FDG_PONS_AV45_2_DATE'] = av45_2_closest['EXAMDATE'] if av45_2_closest else ''
        new_data['FDG_PONS_AV45_3_6MTHS'] = av45_3_closest['PONS'] if av45_3_closest else ''
        new_data['FDG_PONS_AV45_3_DATE'] = av45_3_closest['EXAMDATE'] if av45_3_closest else ''

        # Get closest to AV1451 scans
        av1451_bl_closest, av1451_2_closest, av1451_3_closest = getClosestToScans(subj_fdg, bl_av1451, av1451_2, av1451_3)
        new_data['FDG_PONS_AV1451_1'] = av1451_bl_closest['PONS'] if av1451_bl_closest else ''
        new_data['FDG_PONS_AV1451_1_DATE'] = av1451_bl_closest['EXAMDATE'] if av1451_bl_closest else ''
        new_data['FDG_PONS_AV1451_2'] = av1451_2_closest['PONS'] if av1451_2_closest else ''
        new_data['FDG_PONS_AV1451_2_DATE'] = av1451_2_closest['EXAMDATE'] if av1451_2_closest else ''
        new_data['FDG_PONS_AV1451_3'] = av1451_3_closest['PONS'] if av1451_3_closest else ''
        new_data['FDG_PONS_AV1451_3_DATE'] = av1451_3_closest['EXAMDATE'] if av1451_3_closest else ''

        slope_points = []
        for i in range(10):
            if i < len(subj_fdg):
                datapoint = subj_fdg[i]
                new_data['FDG_pons_Vis%s' % (i+1)] = datapoint['PONS']
                new_data['FDGVis%s_Date' % (i+1)] = datapoint['EXAMDATE']
                if av45_bl:
                    rel_time = (datapoint['EXAMDATE'] - av45_bl).days / 365.0
                    new_data['FDGVis%s_ReltoAV45' % (i+1)] = rel_time
                    if rel_time >= -(30.0/365.0):
                        slope_points.append((rel_time, datapoint['PONS']))
                else:
                    new_data['FDGVis%s_ReltoAV45' % (i+1)] = ''
            else:
                new_data['FDG_pons_Vis%s' % (i+1)] = ''
                new_data['FDGVis%s_ReltoAV45' % (i+1)] = ''
        if len(slope_points) >= 2:
            raw_dates = [_[0] for _ in slope_points]
            raw_scores = [_[1] for _ in slope_points]
            slope, intercept, r, p, stderr = linregress(raw_dates, raw_scores)
            new_data['FDG_postAV45_slope'] = slope
            new_data['FDG_postAV45_followuptime'] = raw_dates[-1]
        else:
            new_data['FDG_postAV45_slope'] = ''
            new_data['FDG_postAV45_followuptime'] = ''

        return new_data
    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, fdg_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)


    return (new_headers, new_lines)


def syncTBMSynData(old_headers, old_lines, tbm_file, registry_file):
    tbm_by_subj = importTBMSyn(tbm_file)

    # add new headers as needed
    to_add_headers = ['TBMSyn_DATE.%s' % (i+1) for i in range(10)]
    to_add_headers.extend(['TBMSyn_SCORE.%s' % (i+1) for i in range(10)])
    to_add_headers.extend(['TBMSyn_postAV45.%s' % (i+1) for i in range(10)])
    to_add_headers.append('TBMSyn_BL_DATE')
    to_add_headers.append('TBMSyn_count')
    #to_add_headers.append('TBMSyn_SLOPE')
    #to_add_headers.append('TBMSyn_slope_followuptime')
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        subj_tbm = sorted(subj_row, key=lambda x: x['EXAMDATE'])
        av45_bl, av45_2, av45_3 = getAV45Dates(old_l)
        bl_examdate = subj_tbm[0]['EXAMDATEBL']
        if av45_bl is None:
            print "No Baseline AV45 date for %s" % subj
            return {}
        # check that TBMSyn baseline is within range of AV45 baseline
        if abs(av45_bl - bl_examdate).days > 210:
            print "TBM BL and AV45 BL %s days apart for subj %s" % (abs(av45_bl - bl_examdate).days, subj)
            return {}

        new_data = {}
        slope_points = []
        for i in range(10):
            if i < len(subj_tbm):
                datapoint = subj_tbm[i]
                new_data['TBMSyn_DATE.%s' % (i+1)] = datapoint['EXAMDATE']
                new_data['TBMSyn_SCORE.%s' % (i+1)] = datapoint['TBMSYNSCOR']
                timediff = (datapoint['EXAMDATE']-datapoint['EXAMDATEBL']).days / 365.0
                slope_points.append((timediff,datapoint['TBMSYNSCOR']))
                new_data['TBMSyn_postAV45.%s' % (i+1)] = timediff
            else:
                new_data['TBMSyn_DATE.%s' % (i+1)] = ''
                new_data['TBMSyn_SCORE.%s' % (i+1)] = ''
                new_data['TBMSyn_postAV45.%s' % (i+1)] = ''
        '''
        # calculate slope
        if len(slope_points) >= 2:
            raw_dates = [_[0] for _ in slope_points]
            raw_scores = [_[1] for _ in slope_points]
            diff_dates = [0.0] + raw_dates
            diff_scores = [0.0] + raw_scores
            #slope, intercept, r, p, stderr = linregress(diff_dates, diff_scores)
            output = optimize.fmin(lambda b, x, y: ((b*x-y)**2).sum(), x0=0.1, args=(diff_dates, diff_scores))
            new_data['TBMSyn_SLOPE'] = output[0]
            new_data['TBMSyn_slope_followuptime'] = raw_dates[-1]
        else:
            new_data['TBMSyn_SLOPE'] = ''
            new_data['TBMSyn_slope_followuptime'] = ''
        '''
        new_data['TBMSyn_count'] = len(slope_points)
        new_data['TBMSyn_BL_DATE'] = bl_examdate
        return new_data

    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, tbm_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)


    return (new_headers, new_lines)

def syncCSFData(old_headers, old_lines, csf_files, registry_file):
    tmpts = 7
    registry = importRegistry(registry_file)
    csf_by_subj = importCSF(csf_files, registry)
    
    # add new headers as needed
    to_add_headers = []
    to_add_headers += ['CSF_ABETA.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['CSF_ABETApostAV45.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['CSF_ABETA_slope', 
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
    to_add_headers += ['CSF_TAU.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['CSF_TAUpostAV45.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['CSF_TAU_slope',
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
    to_add_headers += ['CSF_PTAU.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['CSF_PTAUpostAV45.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['CSF_PTAU_slope',
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
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    # get dates for each subject csf
    for subj in csf_by_subj.keys():
        subj_csf = csf_by_subj[subj]
        subj_csf_filtered = []
        # find the examdate for each point, then sort
        subj_listings = registry[subj]

        for i in range(len(subj_csf)):
            csf_row = subj_csf[i]
            viscode = csf_row['vc']
            viscode2 = csf_row['vc2']
            # check if exam date already in there
            if csf_row['EXAMDATE'] is not None:
                subj_csf_filtered.append(csf_row)
                continue
            # if not, get date from registry listings
            examdate = None
            for listing in subj_listings:
                if listing['VISCODE'] == viscode or listing['VISCODE2'] == viscode2:
                    examdate = listing['EXAMDATE']
                    break
            if examdate is None:
                print "No exam date for %s (%s, %s)" % (subj, viscode, viscode2)
                continue
            else:
                csf_row['EXAMDATE'] = examdate
                subj_csf_filtered.append(csf_row)
        subj_csf = sorted(subj_csf_filtered, key=lambda x: x['EXAMDATE'])
        csf_by_subj[subj] = subj_csf

    def extraction_fn(subj, subj_csf, old_l, patient_pets):
        if len(subj_csf) > tmpts:
            raise Exception("Raise CSF timepoints to %s" % len(subj_csf))

        # find av45 dates
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        bl_av1451, av1451_2, av1451_3 = getAV1451Dates(old_l, patient_pets=patient_pets)

        new_data = {}
        ptau_slope_points = []
        abeta_slope_points = []
        tau_slope_points = []
        for i in range(len(subj_csf)):
            datapoint = subj_csf[i]
            examdate = datapoint['EXAMDATE']
            timediff = ((examdate-bl_av45).days / 365.0) if bl_av45 else ''
            if timediff <= -(90.0/365.0):
                timediff = ''
            abeta_val = datapoint['abeta']
            ptau_val = datapoint['ptau']
            tau_val = datapoint['tau']
            new_data['CSF_ABETA.%s' % (i+1)] = abeta_val
            new_data['CSF_PTAU.%s' % (i+1)] = ptau_val
            new_data['CSF_TAU.%s' % (i+1)] = tau_val
            new_data['CSF_ABETApostAV45.%s' % (i+1)] = timediff
            new_data['CSF_PTAUpostAV45.%s' % (i+1)] = timediff
            new_data['CSF_TAUpostAV45.%s' % (i+1)] = timediff
            if timediff != '':
                if abeta_val is not None:
                    abeta_slope_points.append((timediff, abeta_val))
                if tau_val is not None:
                    tau_slope_points.append((timediff, tau_val))
                if ptau_val is not None:
                    ptau_slope_points.append((timediff, ptau_val))

        # match up with av45 scans
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToScans(subj_csf, bl_av45, av45_2, av45_3, day_limit=93)
        for i, closest_av45 in enumerate([av45_bl_closest, av45_2_closest, av45_3_closest]):
            if closest_av45 is None:
                continue
            idx = i+1
            tau = closest_av45['tau']
            ptau = closest_av45['ptau']
            abeta = closest_av45['abeta']
            new_data['CSF_TAU_closest_AV45_%s' % idx] = tau
            new_data['CSF_PTAU_closest_AV45_%s' % idx] = ptau 
            new_data['CSF_ABETA_closest_AV45_%s' % idx] = abeta
            if tau is not None:
                new_data['CSF_TAU_closest_AV45_%s_BIN_93' % idx] = 1 if tau >= 93 else 0
            if ptau is not None:
                new_data['CSF_PTAU_closest_AV45_%s_BIN_23' % idx] = 1 if ptau >= 23 else 0
            if abeta is not None:
                new_data['CSF_ABETA_closest_AV45_%s_BIN_192' % idx] = 1 if abeta <= 192 else 0

        # match up with av1451 scans
        av1451_bl_closest, av1451_2_closest, av1451_3_closest = getClosestToScans(subj_csf, bl_av1451, av1451_2, av1451_3)
        for i, closest_av1451 in enumerate([av1451_bl_closest, av1451_2_closest, av1451_3_closest]):
            if closest_av1451 is None:
                continue
            idx = i+1
            tau = closest_av1451['tau']
            ptau = closest_av1451['ptau']
            abeta = closest_av1451['abeta']
            new_data['CSF_TAU_closest_AV1451_%s' % idx] = tau
            new_data['CSF_PTAU_closest_AV1451_%s' % idx] = ptau 
            new_data['CSF_ABETA_closest_AV1451_%s' % idx] = abeta
            if tau is not None:
                new_data['CSF_TAU_closest_AV1451_%s_BIN_93' % idx] = 1 if tau >= 93 else 0
            if ptau is not None:
                new_data['CSF_PTAU_closest_AV1451_%s_BIN_23' % idx] = 1 if ptau >= 23 else 0
            if abeta is not None:
                new_data['CSF_ABETA_closest_AV1451_%s_BIN_192' % idx] = 1 if abeta <= 192 else 0

        # calculate slope
        if len(ptau_slope_points) >= 2:
            ptau_slope = slope(ptau_slope_points)
            new_data['CSF_PTAU_slope'] = ptau_slope if ptau_slope is not None else ''
        if len(tau_slope_points) >= 2:
            tau_slope = slope(tau_slope_points)
            new_data['CSF_TAU_slope'] = tau_slope if tau_slope is not None else ''
        if len(abeta_slope_points) >= 2:
            abeta_slope = slope(abeta_slope_points)
            new_data['CSF_ABETA_slope'] = abeta_slope if abeta_slope is not None else ''
        return new_data


    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, csf_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)


    return (new_headers, new_lines)

def syncUWData(old_headers, old_lines, uw_file, registry_file):
    tmpts=12
    registry = importRegistry(registry_file)
    uw_by_subj = importUW(uw_file, registry=registry)

    to_add_headers = []
    to_add_headers += ['UW_MEM_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['UW_MEM_postAV45_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['UW_MEM_postAV45_count', 'UW_MEM_slope', 'UW_MEM_BL_3months']
    to_add_headers += ['UW_EF_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['UW_EF_postAV45_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['UW_EF_postAV45_count', 'UW_EF_slope', 'UW_EF_BL_3months']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        subj_uw = sorted(subj_row, key=lambda x: x['EXAMDATE'])
        if len(subj_uw) > tmpts:
            raise Exception("INCREASE NUMBER OF UW TIMEPOINTS TO > %s" % (len(subj_uw)))
        
        new_data = {}
        mem_slope_points = []
        ef_slope_points = []
        for i in range(tmpts):
            if i < len(subj_uw):
                datapoint = subj_uw[i]
                examdate = datapoint['EXAMDATE']
                cur_mem = datapoint['ADNI_MEM']
                cur_ef = datapoint['ADNI_EF']
                timediff = ((examdate-bl_av45).days / 365.0) if bl_av45 else ''
                if timediff <= -(90.0/365.0):
                    timediff = ''
                if not np.isnan(cur_mem):
                    new_data['UW_MEM_%s' % (i+1)] = cur_mem
                    new_data['UW_MEM_postAV45_%s' % (i+1)] = timediff
                    if timediff != '':
                        mem_slope_points.append((timediff, cur_mem))
                        if abs(timediff) <= (90.0/365.0) and 'UW_MEM_BL_3months' not in new_data:
                            new_data['UW_MEM_BL_3months'] = cur_mem
                if not np.isnan(cur_ef):
                    new_data['UW_EF_%s' % (i+1)] = cur_ef
                    new_data['UW_EF_postAV45_%s' % (i+1)] = timediff
                    if timediff != '':
                        ef_slope_points.append((timediff, cur_ef))
                        if abs(timediff) <= (90.0/365.0) and 'UW_EF_BL_3months' not in new_data:
                            new_data['UW_EF_BL_3months'] = cur_ef
        # calc slope
        new_data['UW_MEM_postAV45_count'] = len(mem_slope_points)
        new_data['UW_MEM_slope'] = slope(mem_slope_points)
        new_data['UW_EF_postAV45_count'] = len(ef_slope_points)
        new_data['UW_EF_slope'] = slope(ef_slope_points)
        return new_data

    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, uw_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)


    return (new_headers, new_lines)

def syncUCBFreesurferData(old_headers, old_lines, ucb_fs_volumes):
    tmpts = 3
    fs_by_subj = importFSVolumes(ucb_fs_volumes)

    to_add_headers = []
    to_add_headers += ['UCB_FS_HC/ICV_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['UCB_FS_postAV45_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['UCB_FS_postAV45_count', 'UCB_FS_postAV45_interval', 
                       'UCB_FS_HC/ICV_slope', 'UCB_FS_HC/ICVavg_slope', 
                       'UCB_FS_HC/ICV_AV45_1_3MTHS','UCB_FS_HC/ICV_AV45_1_DATE',
                       'UCB_FS_HC/ICV_AV45_2_3MTHS','UCB_FS_HC/ICV_AV45_2_DATE',
                       'UCB_FS_HC/ICV_AV45_3_3MTHS','UCB_FS_HC/ICV_AV45_3_DATE',
                       'UCB_FS_HC/ICV_AV1451_1','UCB_FS_HC/ICV_AV1451_1_DATE',
                       'UCB_FS_HC/ICV_AV1451_2','UCB_FS_HC/ICV_AV1451_2_DATE',
                       'UCB_FS_HC/ICV_AV1451_3','UCB_FS_HC/ICV_AV1451_3_DATE',]
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        bl_av1451, av1451_2, av1451_3 = getAV1451Dates(old_l, patient_pets=patient_pets)
        subj_fs = sorted(subj_row, key=lambda x: x['Date'])
        if len(subj_fs) > tmpts:
            raise Exception("INCREASE NUMBER OF UCB FS TIMEPOINTS TO > %s" % (len(subj_fs)))

        new_data = {}
        slope_points = []
        avg_slope_points = []
        bl_icv = subj_fs[0]['EstimatedTotalIntraCranialVol']
        avg_icv = np.mean([float(_['EstimatedTotalIntraCranialVol']) for _ in subj_fs])

        # get closest to av45s
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToScans(subj_fs, bl_av45, av45_2, av45_3, day_limit=31*3, date_key='Date')
        new_data['UCB_FS_HC/ICV_AV45_1_3MTHS'] = av45_bl_closest['HCV']/bl_icv if av45_bl_closest else ''
        new_data['UCB_FS_HC/ICV_AV45_1_DATE'] = av45_bl_closest['Date'] if av45_bl_closest else ''
        new_data['UCB_FS_HC/ICV_AV45_2_3MTHS'] = av45_2_closest['HCV']/bl_icv if av45_2_closest else ''
        new_data['UCB_FS_HC/ICV_AV45_2_DATE'] = av45_2_closest['Date'] if av45_2_closest else ''
        new_data['UCB_FS_HC/ICV_AV45_3_3MTHS'] = av45_3_closest['HCV']/bl_icv if av45_3_closest else ''
        new_data['UCB_FS_HC/ICV_AV45_3_DATE'] = av45_3_closest['Date'] if av45_3_closest else ''

        # Get closest to AV1451 scans
        av1451_bl_closest, av1451_2_closest, av1451_3_closest = getClosestToScans(subj_fs, bl_av1451, av1451_2, av1451_3, date_key='Date')
        new_data['UCB_FS_HC/ICV_AV1451_1'] = av1451_bl_closest['HCV']/bl_icv if av1451_bl_closest else ''
        new_data['UCB_FS_HC/ICV_AV1451_1_DATE'] = av1451_bl_closest['Date'] if av1451_bl_closest else ''
        new_data['UCB_FS_HC/ICV_AV1451_2'] = av1451_2_closest['HCV']/bl_icv if av1451_2_closest else ''
        new_data['UCB_FS_HC/ICV_AV1451_2_DATE'] = av1451_2_closest['Date'] if av1451_2_closest else ''
        new_data['UCB_FS_HC/ICV_AV1451_3'] = av1451_3_closest['HCV']/bl_icv if av1451_3_closest else ''
        new_data['UCB_FS_HC/ICV_AV1451_3_DATE'] = av1451_3_closest['Date'] if av1451_3_closest else ''

        for i in range(len(subj_fs)):
            datapoint = subj_fs[i]
            examdate = datapoint['Date']
            hc_icv = datapoint['HCV']/bl_icv
            hc_icv_avg = datapoint['HCV']/avg_icv
            new_data['UCB_FS_HC/ICV_%s' % (i+1)] = hc_icv
            timediff = ((examdate-bl_av45).days / 365.0) if bl_av45 else ''
            if timediff != '' and timediff > -(90.0/365.0):
                slope_points.append((timediff, hc_icv))
                avg_slope_points.append((timediff,hc_icv_avg))
                new_data['UCB_FS_postAV45_%s' % (i+1)] = timediff
            
        # calc slope
        new_data['UCB_FS_postAV45_count'] = len(slope_points)
        new_data['UCB_FS_HC/ICV_slope'] = slope(slope_points)
        new_data['UCB_FS_HC/ICVavg_slope'] = slope(avg_slope_points)
        if len(slope_points) > 0:
            new_data['UCB_FS_postAV45_interval'] = max([_[0] for _ in slope_points])
        else:
            new_data['UCB_FS_postAV45_interval'] = 0
        return new_data

    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, fs_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)


    return (new_headers, new_lines)

def syncUCSFFreesurferCrossData(old_headers, old_lines, ucsf_files, mprage_file):
    tmpts = 12

    fs_by_subj_list = [importUCSFFreesurfer(filename, mprage_file, version=fsversion, include_failed=False) for fsversion,filename in ucsf_files]
    fs_by_subj = {}
    for k,v in itertools.chain(*[_.iteritems() for _ in fs_by_subj_list]):
        if k not in fs_by_subj:
            fs_by_subj[k] = v
        else:
            fs_by_subj[k] = fs_by_subj[k] + v

    to_add_headers = []
    to_add_headers += ['FSX_HC/ICV_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FSX_postAV45_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FSX_postAV45_count', 'FSX_postAV45_interval', 
                       'FSX_HC/ICV_slope', 'FSX_HC/ICVavg_slope',
                       'FSX_HC/ICV_AV45_1_3MTHS','FSX_HC/ICV_AV45_1_DATE',
                       'FSX_HC/ICV_AV45_2_3MTHS','FSX_HC/ICV_AV45_2_DATE',
                       'FSX_HC/ICV_AV45_3_3MTHS','FSX_HC/ICV_AV45_3_DATE',
                       'FSX_HC/ICV_AV1451_1','FSX_HC/ICV_AV1451_1_DATE',
                       'FSX_HC/ICV_AV1451_2','FSX_HC/ICV_AV1451_2_DATE',
                       'FSX_HC/ICV_AV1451_3','FSX_HC/ICV_AV1451_3_DATE',]
    to_add_headers += ['FSX_MRI_STRENGTH_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FSX_FSVERSION_%s' % (i+1) for i in range(tmpts)]
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        '''
        If RID<2000, use only 1.5T scans
        If RID>2000, use only 3T scans
        '''
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        bl_av1451, av1451_2, av1451_3 = getAV1451Dates(old_l, patient_pets=patient_pets)

        # FILTER: Only use 1.5T scans if 3T scans aren't present
        onefive_scans = [_ for _ in subj_row if float(_['FLDSTRENG']) == 1.5]
        three_scans = [_ for _ in subj_row if float(_['FLDSTRENG']) == 3.0]
        if len(three_scans) == 0:
            subj_row = onefive_scans
        else:
            subj_row = three_scans
        if len(subj_row) == 0:
            print "%s had all potential scans eliminated by field strength: %s" % (subj,[_['FLDSTRENG'] for _ in subj_row])
            return {}

        # take care of scans on the same date
        by_date = defaultdict(list)
        for scan in subj_row:
            by_date[scan['EXAMDATE']].append(scan)
        subj_fs = []
        for scandate in sorted(by_date.keys()):
            scans_on_date = by_date[scandate]
            if len(scans_on_date) == 1:
                subj_fs.append(scans_on_date[0])
            else:
                # average if possible, else take first one
                if len(set([_['FLDSTRENG'] for _ in scans_on_date])) == 1 and len(set([_['version'] for _ in scans_on_date])) == 1:
                    avg_row = {}
                    avg_row['EXAMDATE'] = scandate
                    avg_row['ICV'] = np.mean([float(_['ICV']) for _ in scans_on_date])
                    avg_row['HCV'] = np.mean([float(_['HCV']) for _ in scans_on_date])
                    avg_row['FLDSTRENG'] = scans_on_date[0]['FLDSTRENG']
                    avg_row['version'] = scans_on_date[0]['version']
                else:
                    avg_row = scans_on_date[0]
                subj_fs.append(avg_row)

        if len(subj_fs) > tmpts:
            raise Exception("INCREASE NUMBER OF FSX TIMEPOINTS TO %s" % (len(subj_fs)))

        new_data = {}
        slope_points = []
        avg_slope_points = []
        bl_icv = subj_fs[0]['ICV']
        avg_icv = np.mean([float(_['ICV']) for _ in subj_fs])

        # get closest to av45s
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToScans(subj_fs, bl_av45, av45_2, av45_3,day_limit=31*3)
        new_data['FSX_HC/ICV_AV45_1_3MTHS'] = av45_bl_closest['HCV']/bl_icv if av45_bl_closest else ''
        new_data['FSX_HC/ICV_AV45_1_DATE'] = av45_bl_closest['EXAMDATE'] if av45_bl_closest else ''
        new_data['FSX_HC/ICV_AV45_2_3MTHS'] = av45_2_closest['HCV']/bl_icv if av45_2_closest else ''
        new_data['FSX_HC/ICV_AV45_2_DATE'] = av45_2_closest['EXAMDATE'] if av45_2_closest else ''
        new_data['FSX_HC/ICV_AV45_3_3MTHS'] = av45_3_closest['HCV']/bl_icv if av45_3_closest else ''
        new_data['FSX_HC/ICV_AV45_3_DATE'] = av45_3_closest['EXAMDATE'] if av45_3_closest else ''

        # Get closest to AV1451 scans
        av1451_bl_closest, av1451_2_closest, av1451_3_closest = getClosestToScans(subj_fs, bl_av1451, av1451_2, av1451_3)
        new_data['FSX_HC/ICV_AV1451_1'] = av1451_bl_closest['HCV']/bl_icv if av1451_bl_closest else ''
        new_data['FSX_HC/ICV_AV1451_1_DATE'] = av1451_bl_closest['EXAMDATE'] if av1451_bl_closest else ''
        new_data['FSX_HC/ICV_AV1451_2'] = av1451_2_closest['HCV']/bl_icv if av1451_2_closest else ''
        new_data['FSX_HC/ICV_AV1451_2_DATE'] = av1451_2_closest['EXAMDATE'] if av1451_2_closest else ''
        new_data['FSX_HC/ICV_AV1451_3'] = av1451_3_closest['HCV']/bl_icv if av1451_3_closest else ''
        new_data['FSX_HC/ICV_AV1451_3_DATE'] = av1451_3_closest['EXAMDATE'] if av1451_3_closest else ''

        for i in range(len(subj_fs)):
            datapoint = subj_fs[i]
            examdate = datapoint['EXAMDATE']
            hc_icv = datapoint['HCV']/bl_icv
            hc_icv_avg = datapoint['HCV']/avg_icv
            timediff = ((examdate-bl_av45).days / 365.0) if bl_av45 else ''
            if timediff <= -90.0/365.0:
                timediff = ''
            if timediff != '':
                slope_points.append((timediff, hc_icv))
                avg_slope_points.append((timediff,hc_icv_avg))
            new_data['FSX_HC/ICV_%s' % (i+1)] = hc_icv
            new_data['FSX_postAV45_%s' % (i+1)] = timediff
            new_data['FSX_FSVERSION_%s' % (i+1)] = datapoint['version']
            new_data['FSX_MRI_STRENGTH_%s' % (i+1)] = datapoint['FLDSTRENG']
    
        # calc slope
        new_data['FSX_postAV45_count'] = len(slope_points)
        new_data['FSX_HC/ICV_slope'] = slope(slope_points)
        new_data['FSX_HC/ICVavg_slope'] = slope(avg_slope_points)
        if len(slope_points) > 0:
            new_data['FSX_postAV45_interval'] = max([_[0] for _ in slope_points])
        else:
            new_data['FSX_postAV45_interval'] = 0
        
        return new_data

    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, fs_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)


    return (new_headers, new_lines)


def syncUCSFFreesurferLongData(old_headers, old_lines, ucsf_files, mprage_file):
    tmpts = 12

    fs_by_subj_list = [importUCSFFreesurfer(filename, mprage_file, version=fsversion, include_failed=False) for fsversion,filename in ucsf_files]
    fs_by_subj = {}
    for k,v in itertools.chain(*[_.iteritems() for _ in fs_by_subj_list]):
        if k not in fs_by_subj:
            fs_by_subj[k] = v
        else:
            fs_by_subj[k] = fs_by_subj[k] + v

    to_add_headers = []
    to_add_headers += ['FSL_HC/ICV_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FSL_postAV45_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FSL_postAV45_count', 'FSL_postAV45_interval', 
                       'FSL_HC/ICV_slope', 'FSL_HC/ICVavg_slope',
                       'FSL_HC/ICV_AV45_1_3MTHS','FSL_HC/ICV_AV45_1_DATE',
                       'FSL_HC/ICV_AV45_2_3MTHS','FSL_HC/ICV_AV45_2_DATE',
                       'FSL_HC/ICV_AV45_3_3MTHS','FSL_HC/ICV_AV45_3_DATE',
                       'FSL_HC/ICV_AV1451_1','FSL_HC/ICV_AV1451_1_DATE',
                       'FSL_HC/ICV_AV1451_2','FSL_HC/ICV_AV1451_2_DATE',
                       'FSL_HC/ICV_AV1451_3','FSL_HC/ICV_AV1451_3_DATE',]
    to_add_headers += ['FSL_MRI_STRENGTH_%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['FSL_FSVERSION_%s' % (i+1) for i in range(tmpts)]
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        bl_av1451, av1451_2, av1451_3 = getAV1451Dates(old_l, patient_pets=patient_pets)
        mristrengths = [_['FLDSTRENG'] for _ in subj_row]
        
        # FILTER: Only use 1.5T scans if 3T scans aren't present
        onefive_scans = [_ for _ in subj_row if float(_['FLDSTRENG']) == 1.5]
        three_scans = [_ for _ in subj_row if float(_['FLDSTRENG']) == 3.0]
        if len(three_scans) == 0:
            subj_row = onefive_scans
        else:
            subj_row = three_scans
        if len(subj_row) == 0:
            print "%s had all potential scans eliminated by field strength: %s" % (subj,[_['FLDSTRENG'] for _ in subj_row])
            return {}

        # take care of scans on the same date
        by_date = defaultdict(list)
        for scan in subj_row:
            by_date[scan['EXAMDATE']].append(scan)
        subj_fs = []
        for scandate in sorted(by_date.keys()):
            scans_on_date = by_date[scandate]
            if len(scans_on_date) == 1:
                subj_fs.append(scans_on_date[0])
            else:
                if len(set([_['FLDSTRENG'] for _ in scans_on_date])) == 1 and len(set([_['version'] for _ in scans_on_date])) == 1:
                    # average
                    avg_row = {}
                    avg_row['EXAMDATE'] = scandate
                    avg_row['ICV'] = np.mean([float(_['ICV']) for _ in scans_on_date])
                    avg_row['HCV'] = np.mean([float(_['HCV']) for _ in scans_on_date])
                    avg_row['FLDSTRENG'] = scans_on_date[0]['FLDSTRENG']
                    avg_row['version'] = scans_on_date[0]['version']
                else:
                    avg_row = scans_on_date[0]
                subj_fs.append(avg_row)

        if len(subj_fs) > tmpts:
            raise Exception("INCREASE NUMBER OF FSL TIMEPOINTS TO > %s" % (len(subj_fs)))

        new_data = {}
        slope_points = []
        avg_slope_points = []
        bl_icv = subj_fs[0]['ICV']
        avg_icv = np.mean([float(_['ICV']) for _ in subj_fs])

        # get closest to av45s
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToScans(subj_fs, bl_av45, av45_2, av45_3,day_limit=31*3)
        new_data['FSL_HC/ICV_AV45_1_3MTHS'] = av45_bl_closest['HCV']/bl_icv if av45_bl_closest else ''
        new_data['FSL_HC/ICV_AV45_1_DATE'] = av45_bl_closest['EXAMDATE'] if av45_bl_closest else ''
        new_data['FSL_HC/ICV_AV45_2_3MTHS'] = av45_2_closest['HCV']/bl_icv if av45_2_closest else ''
        new_data['FSL_HC/ICV_AV45_2_DATE'] = av45_2_closest['EXAMDATE'] if av45_2_closest else ''
        new_data['FSL_HC/ICV_AV45_3_3MTHS'] = av45_3_closest['HCV']/bl_icv if av45_3_closest else ''
        new_data['FSL_HC/ICV_AV45_3_DATE'] = av45_3_closest['EXAMDATE'] if av45_3_closest else ''

        # Get closest to AV1451 scans
        av1451_bl_closest, av1451_2_closest, av1451_3_closest = getClosestToScans(subj_fs, bl_av1451, av1451_2, av1451_3)
        new_data['FSL_HC/ICV_AV1451_1'] = av1451_bl_closest['HCV']/bl_icv if av1451_bl_closest else ''
        new_data['FSL_HC/ICV_AV1451_1_DATE'] = av1451_bl_closest['EXAMDATE'] if av1451_bl_closest else ''
        new_data['FSL_HC/ICV_AV1451_2'] = av1451_2_closest['HCV']/bl_icv if av1451_2_closest else ''
        new_data['FSL_HC/ICV_AV1451_2_DATE'] = av1451_2_closest['EXAMDATE'] if av1451_2_closest else ''
        new_data['FSL_HC/ICV_AV1451_3'] = av1451_3_closest['HCV']/bl_icv if av1451_3_closest else ''
        new_data['FSL_HC/ICV_AV1451_3_DATE'] = av1451_3_closest['EXAMDATE'] if av1451_3_closest else ''

        for i in range(len(subj_fs)):
            datapoint = subj_fs[i]
            examdate = datapoint['EXAMDATE']
            hc_icv = datapoint['HCV']/bl_icv
            hc_icv_avg = datapoint['HCV']/avg_icv
            timediff = ((examdate-bl_av45).days / 365.0) if bl_av45 else ''
            if timediff <= -(90.0/365.0):
                timediff = ''
            if timediff != '':
                avg_slope_points.append((timediff, hc_icv_avg))
                slope_points.append((timediff, hc_icv))
            new_data['FSL_HC/ICV_%s' % (i+1)] = datapoint['HCV']/bl_icv
            new_data['FSL_postAV45_%s' % (i+1)] = timediff
            new_data['FSL_FSVERSION_%s' % (i+1)] = datapoint['version']
            new_data['FSL_MRI_STRENGTH_%s' % (i+1)] = datapoint['FLDSTRENG']
        # Calculate slope
        new_data['FSL_postAV45_count'] = len(slope_points)
        new_data['FSL_HC/ICV_slope'] = slope(slope_points)
        new_data['FSL_HC/ICVavg_slope'] = slope(avg_slope_points)
        if len(slope_points) > 0:
            new_data['FSL_postAV45_interval'] = max([_[0] for _ in slope_points])
        else:
            new_data['FSL_postAV45_interval'] = 0

        return new_data

    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, fs_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)


    return (new_headers, new_lines)


def syncWMHData(old_headers, old_lines, wmh_file, registry_file):
    tmpts = 7
    wmh_by_subj = importWMH(wmh_file)
    registry = importRegistry(registry_file)

    to_add_headers = ['WMH_percentOfICV.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['WMH_postAV45.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['WMH_WHITMATHYP.%s' % (i+1) for i in range(tmpts)]
    to_add_headers += ['WMH_percentOfICV_AV45_6MTHS', 'WMH_percentOfICV_AV45_DATE',
                       'WMH_percentOfICV_AV45_2_6MTHS', 'WMH_percentOfICV_AV45_2_DATE',
                       'WMH_percentOfICV_AV45_3_6MTHS', 'WMH_percentOfICV_AV45_3_DATE',
                       'WMH_percentOfICV_slope']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    wipeKeys(old_lines, to_add_headers)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        subj_wmh = sorted(subj_row, key=lambda x: x['EXAMDATE'])
        new_data = {}

        # get closest to av45s
        six_months = 31 * 6
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToScans(subj_wmh, bl_av45, av45_2, av45_3, day_limit=six_months)
        new_data['WMH_percentOfICV_AV45_6MTHS'] = av45_bl_closest['wmh_percent'] if av45_bl_closest else ''
        new_data['WMH_percentOfICV_AV45_DATE'] = av45_bl_closest['EXAMDATE'] if av45_bl_closest else ''
        new_data['WMH_percentOfICV_AV45_2_6MTHS'] = av45_2_closest['wmh_percent'] if av45_2_closest else ''
        new_data['WMH_percentOfICV_AV45_2_DATE'] = av45_2_closest['EXAMDATE'] if av45_2_closest else ''
        new_data['WMH_percentOfICV_AV45_3_6MTHS'] = av45_3_closest['wmh_percent'] if av45_3_closest else ''
        new_data['WMH_percentOfICV_AV45_3_DATE'] = av45_3_closest['EXAMDATE'] if av45_3_closest else ''
        
        slope_points = []
        for i in range(tmpts):
            if i < len(subj_wmh):
                datapoint = subj_wmh[i]
                examdate = datapoint['EXAMDATE']
                wmh_raw = datapoint['wmh']
                wmh_val = datapoint['wmh_percent']
                timediff = ((examdate-bl_av45).days / 365.0) if bl_av45 else ''
                if timediff <= -(90.0/365.0):
                    timediff = ''
                new_data['WMH_percentOfICV.%s' % (i+1)] = wmh_val
                new_data['WMH_postAV45.%s' % (i+1)] = timediff
                new_data['WMH_WHITMATHYP.%s' % (i+1)] = wmh_raw
                if timediff and wmh_val:
                    slope_points.append((timediff, wmh_val))
            else:
                new_data['WMH_percentOfICV.%s' % (i+1)] = ''
                new_data['WMH_postAV45.%s' % (i+1)] = ''
                new_data['WMH_WHITMATHYP.%s' % (i+1)] = ''
        # calculate slope
        new_data['WMH_percentOfICV_slope'] = slope(slope_points)
        return new_data

    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, wmh_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)


    return (new_headers, new_lines)

def eliminateColumns(master_df):
    to_remove = ['AVLT_3MTHS_AV45',
                 'AVLT_3MTHSAV45_Date',
                 'LastClinicalVisit',
                 'DIFF_LastClinicalVisit_AV45',
                 'ClinicalVisitClosestto_AV45',
                 'ClinicalVisitClosestto_AV45_2',
                 'MMSEclosest_1', 
                 'MMSEclosest_AV45',
                 'datediff_MMSE_AV45', 
                 'MMSE_3mths_AV45',
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
                 'CSF<3mths_AV45_2',
                 'LastCSF',
                 'LastCSF_AV45_1_TimeDiff',
                 'LastCSFAbeta',
                 'ADAS_slope_all',
                 'FDG_slope',
                 'TBMSyn_SLOPE',
                 'TBMSyn_slope_followuptime',
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
                 'CSF_ABETA_closest_AV45',
                 'CSF_ABETA_closest_AV45_BIN_192',
                 'CSF_TAU_closest_AV45_BIN_93',
                 'CSF_PTAU_closest_AV45_BIN_23',
                 'BD MM-YY']
    to_remove += ['WHITMATHYP.%s' % (i+1) for i in range(5)]
    to_remove += ['ABETA.%s' % (i+1) for i in range(7)]
    to_remove += ['TAU.%s' % (i+1) for i in range(7)]
    to_remove += ['PTAU.%s' % (i+1) for i in range(7)]
    to_remove += ['EXAM_tau.%s' % (i+1) for i in range(7)]
    to_remove += ['ADAScog_DATE%s' % (i+1) for i in range(11)]
    to_remove += ['MMSE_DATE%s' % (i+1) for i in range(11)]
    to_remove += ['AVLT_DATE.%s' % (i+1) for i in range(11)]
    to_remove += ['TBMSyn_DATE.%s' % (i+1) for i in range(11)]
    to_remove += [_ for _ in master_df.columns if _.startswith('AV45_')]

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
    # print "\nSYNCING AV45 NONTP\n"
    # master_df = syncAV45Data(master_df, av45_nontp_file, suffix='NONTP')
    # print "\nSYNCING AV45 TP\n"
    # master_df = syncAV45Data(master_df, av45_tp_file, suffix='TP')
    print "\nSYNCING DEMOG\n"
    master_df = syncDemogData(master_df, demog_file, pet_meta_file)
    # print "\nSYNCING AV1451 TP\n"
    # master_df = syncAV1451Data(master_df, av1451_tp_file)
    # print "\nSYNCING AV45 ROUSSET PVC\n"
    # master_df = syncAV45RoussetResults(master_df, av45_rousset_csv)
    # print "\nSYNCING AV1451 ROUSSET PVC\n"
    # master_df = syncAV1451RoussetResults(master_df, av1451_rousset_csv)
    # print "\nSYNCING DIAGNOSES\n"
    # master_df = syncDiagnosisData(master_df, diagnosis_file, arm_file, registry)
    # print "\nSYNCING FAQ\n"
    # master_df = syncFAQData(master_df, faq_file, registry)
    # print "\nSYNCING NPI\n"
    # master_df = syncNPIData(master_df, npi_file)
    # print "\nSYNCING MHIST\n"
    # master_df = syncMHISTData(master_df, mhist_file)
    print "\nSYNCING UW NEURO\n"
    master_df = syncUWData(master_df, uw_file, registry)

    dumpDFtoCSV(master_df,output_file,decimal_places=3)
    sys.exit(1)






    print "\nSYNCING UCSF LONG FreeSurfer\n"
    new_headers, new_lines = syncUCSFFreesurferLongData(new_headers, new_lines, ucsf_long_files, mprage_file)
    print "\nSYNCING UCSF CROSS FreeSurfer"
    new_headers, new_lines = syncUCSFFreesurferCrossData(new_headers, new_lines, ucsf_cross_files, mprage_file)
    print "\nSYNCING UCB Freesurfer\n"
    new_headers, new_lines = syncUCBFreesurferData(new_headers, new_lines, ucb_fs_volumes)
    print "\nSYNCING FDG\n"
    new_headers, new_lines = syncFDGData(new_headers, new_lines, fdg_file, registry_file)
    print "\nSYNCING TBMSYN\n"
    new_headers, new_lines = syncTBMSynData(new_headers, new_lines, tbm_file, registry_file)
    print "\nSYNCING MMSE\n"
    new_headers, new_lines = syncMMSEData(new_headers, new_lines, mmse_file, registry_file)
    print "\nSYNCING ADASCOG\n"
    new_headers, new_lines = syncADASCogData(new_headers, new_lines, adni1_adas_file, adnigo2_adas_file, registry_file)
    print "\nSYNCING AVLT\n"
    new_headers, new_lines = syncAVLTData(new_headers, new_lines, neuro_battery_file, registry_file)
    print "\nSYNCING CSF\n"
    new_headers, new_lines = syncCSFData(new_headers, new_lines, csf_files, registry_file)
    print "\nSYNCING GD\n"
    new_headers, new_lines = syncGDData(new_headers, new_lines, gd_file, registry_file)
    print "\nSYNCING WMH\n"
    new_headers, new_lines = syncWMHData(new_headers, new_lines, wmh_file, registry_file)
    
    # Save + Post Process
    dumpCSV(output_file, new_headers, new_lines)
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
    av45_tp_file = "../output/04_05_16/UCBERKELEYAV45_04_05_16_regular_tp.csv"
    av45_nontp_file = "../output/04_05_16/UCBERKELEYAV45_04_05_16_regular_nontp.csv"
    av1451_tp_file = "../output/04_05_16/UCBERKELEYAV1451_04_05_16_regular_tp.csv"
    av45_rousset_csv = "../datasets/pvc_adni_av45/aggregions_output.csv"
    av1451_rousset_csv = "../datasets/pvc_adni_av1451/tauregions_output.csv"

    # Cog files
    mmse_file = "../docs/ADNI/MMSE.csv"
    adni1_adas_file = '../docs/ADNI/ADASSCORES.csv'
    adnigo2_adas_file = '../docs/ADNI/ADAS_ADNIGO2.csv'
    neuro_battery_file = '../docs/ADNI/NEUROBAT.csv'
    gd_file = '../docs/ADNI/GDSCALE.csv'
    uw_file = '../docs/ADNI/UWNPSYCHSUM_01_12_16.csv'
    
    # MR files
    tbm_file = '../docs/ADNI/MAYOADIRL_MRI_TBMSYN_12_08_15.csv'
    ucsf_long_files = [('4.4','../docs/ADNI/UCSFFSL_02_01_16.csv'),
                       ('5.1','../docs/ADNI/UCSFFSL51Y1_02_24_16.csv')]
    ucsf_cross_files = [('4.3','../docs/ADNI/UCSFFSX_11_02_15.csv'),
                        ('5.1','../docs/ADNI/UCSFFSX51_11_02_15_V2.csv'),
                        ('5.1','../docs/ADNI/UCSFFSX51_ADNI1_3T_02_01_16.csv')]
    ucb_fs_volumes = '../docs/ADNI/adni_av45_fs_volumes_04-05-2016.csv'

    # Run pipeline
    runPipeline()
