'''

Pull in the master FDG_AV45_COG DATA

remove null char:
tr -d '\000' < file1 > file2

'''
import re
from collections import defaultdict
from datetime import datetime, timedelta
from scipy import stats, optimize
import numpy as np

from utils import *

'''
def syncTemplate(old_headers, old_lines, input_file, dump_to=None):
    data = importTemplate(input_file)

    to_add_headers = []
    after = None
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=after)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        pass

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, data, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)
'''


COLUMN_CATEGORIES = {'AV45_INFO': ['Diag@AV45_long',
                                  'Diag@AV45_2_long',
                                  'Diag@AV45_3_long',
                                  'AV45_scan_count',
                                  'Age@AV45',
                                  'Age@AV45_2',
                                  'Age@AV45_3',
                                  'AV45_Date',
                                  'AV45_2_Date',
                                  'AV45_3_Date',
                                  'AV451_time',
                                  'AV45_1_2_Diff (Yrs)',
                                  'AV45_1_3_Diff (yrs)',
                                  'AV45_TP_SPECIFIC_BL',
                                  'AV45_TP_SPECIFIC_SCAN2',
                                  'AV45_TP_SPECIFIC_SCAN3'],
                    'OLD': ['SNAPgroups_BigRef_FDG1.25',
                            'SNAPgroups_BigRef_HippVol0.0047',
                            'PosAgreement_FDG1.25andHippVol0.0047',
                            'NegAgreement_FDG1.25andHippVol0.0047',
                            'FDGpos_1.25orHippVolpos_0.0047',
                            'SNAPgroups_BigRef_FDG1.25_OR_HippVol0.0047',
                            'SNAPgroups_BigRef_FDG1.25_AND_HippVol0.0047',
                            'Notes',
                            'Florbetapir_Steffi_reads',
                            'Florbetapir_Steffi_read_notes',
                            'Florbetapir_Bill_reads',
                            'Florbetapir_Bill_read_notes',
                            'Fam_Dementia_Total',
                            'Fam_AD_Total',
                            'DXDSEV',
                            'DXMDUE (1 = MCI likely from AD, 2 = likely other)',
                            'Other_MCI_Etiology',
                            'MCI_Etiology_Notes',
                            'MCI_nonamnestic_symptoms',
                            'DXDDUE (1 = AD likely, 2 = other)',
                            'DXODES',
                            'DXOOTHSP_Notes',
                            'DXAPP (1 = AD diagnosis probable, 2 = possible)',
                            'DXAPOSS (1=Atypical clinical course or features (specify); 2=Stroke(s); 3=Depression; 4=Delirium; 5=Parkinsonism; 6=Metabolic / Toxic Disorder (specify); 7=Other (specify))',
                            'AD_probable_notes',
                            'Voxelwise_Groups_FDG_AV45_VBM',
                            'csf_concordant_voxelwise_groups'],
                    'COGNITIVE': ['MMSCORE.1',
                                  'MMSCORE.2',
                                  'MMSCORE.3',
                                  'MMSCORE.4',
                                  'MMSCORE.5',
                                  'MMSCORE.6',
                                  'MMSCORE.7',
                                  'MMSCORE.8',
                                  'MMSCORE.9',
                                  'MMSCORE.10',
                                  'MMSCORE.11',
                                  'MMSE_post_AV45_followuptime',
                                  'MMSE_AV45_3MTHS',
                                  'MMSE_AV45_DATE',
                                  'MMSE_AV45_2_3MTHS',
                                  'MMSE_AV45_2_DATE',
                                  'MMSE_AV45_3_3MTHS',
                                  'MMSE_AV45_3_DATE',
                                  'MMSEslope_postAV45',
                                  'TIMEpostAV45_MMSE.1',
                                  'TIMEpostAV45_MMSE.2',
                                  'TIMEpostAV45_MMSE.3',
                                  'TIMEpostAV45_MMSE.4',
                                  'TIMEpostAV45_MMSE.5',
                                  'TIMEpostAV45_MMSE.6',
                                  'TIMEpostAV45_MMSE.7',
                                  'TIMEpostAV45_MMSE.8',
                                  'TIMEpostAV45_MMSE.9',
                                  'TIMEpostAV45_MMSE.10',
                                  'TIMEpostAV45_MMSE.11',
                                  'ADAScog.1',
                                  'ADAScog.2',
                                  'ADAScog.3',
                                  'ADAScog.4',
                                  'ADAScog.5',
                                  'ADAScog.6',
                                  'ADAScog.7',
                                  'ADAScog.8',
                                  'ADAScog.9',
                                  'ADAScog.10',
                                  'ADAScog.11',
                                  'TIME_ADAS.1',
                                  'TIME_ADAS.2',
                                  'TIME_ADAS.3',
                                  'TIME_ADAS.4',
                                  'TIME_ADAS.5',
                                  'TIME_ADAS.6',
                                  'TIME_ADAS.7',
                                  'TIME_ADAS.8',
                                  'TIME_ADAS.9',
                                  'TIME_ADAS.10',
                                  'TIME_ADAS.11',
                                  'TIMEreltoAV45_ADAS.1',
                                  'TIMEreltoAV45_ADAS.2',
                                  'TIMEreltoAV45_ADAS.3',
                                  'TIMEreltoAV45_ADAS.4',
                                  'TIMEreltoAV45_ADAS.5',
                                  'TIMEreltoAV45_ADAS.6',
                                  'TIMEreltoAV45_ADAS.7',
                                  'TIMEreltoAV45_ADAS.8',
                                  'TIMEreltoAV45_ADAS.9',
                                  'TIMEreltoAV45_ADAS.10',
                                  'TIMEreltoAV45_ADAS.11',
                                  'TIMEpostAV45_ADAS.1',
                                  'TIMEpostAV45_ADAS.2',
                                  'TIMEpostAV45_ADAS.3',
                                  'TIMEpostAV45_ADAS.4',
                                  'TIMEpostAV45_ADAS.5',
                                  'TIMEpostAV45_ADAS.6',
                                  'TIMEpostAV45_ADAS.7',
                                  'TIMEpostAV45_ADAS.8',
                                  'TIMEpostAV45_ADAS.9',
                                  'TIMEpostAV45_ADAS.10',
                                  'TIMEpostAV45_ADAS.11',
                                  'ADAS_3MTH_AV45',
                                  'ADAS_post_AV45_followuptime',
                                  'ADASslope_postAV45',
                                  'ADAS_3MTHS_AV45DATE',
                                  'ADAS_AV45_2_3MTHS',
                                  'ADAS_AV45_2_DATE',
                                  'ADAS_AV45_3_3MTHS',
                                  'ADAS_AV45_3_DATE',
                                  'AVLT.1',
                                  'AVLT.2',
                                  'AVLT.3',
                                  'AVLT.4',
                                  'AVLT.5',
                                  'AVLT.6',
                                  'AVLT.7',
                                  'AVLT.8',
                                  'AVLT.9',
                                  'AVLT.10',
                                  'AVLT.11',
                                  'TIME_AVLT.1',
                                  'TIME_AVLT.2',
                                  'TIME_AVLT.3',
                                  'TIME_AVLT.4',
                                  'TIME_AVLT.5',
                                  'TIME_AVLT.6',
                                  'TIME_AVLT.7',
                                  'TIME_AVLT.8',
                                  'TIME_AVLT.9',
                                  'TIME_AVLT.10',
                                  'TIME_AVLT.11',
                                  'TIMEreltoAV45_AVLT.1',
                                  'TIMEreltoAV45_AVLT.2',
                                  'TIMEreltoAV45_AVLT.3',
                                  'TIMEreltoAV45_AVLT.4',
                                  'TIMEreltoAV45_AVLT.5',
                                  'TIMEreltoAV45_AVLT.6',
                                  'TIMEreltoAV45_AVLT.7',
                                  'TIMEreltoAV45_AVLT.8',
                                  'TIMEreltoAV45_AVLT.9',
                                  'TIMEreltoAV45_AVLT.10',
                                  'TIMEreltoAV45_AVLT.11',
                                  'TIMEpostAV45_AVLT.1',
                                  'TIMEpostAV45_AVLT.2',
                                  'TIMEpostAV45_AVLT.3',
                                  'TIMEpostAV45_AVLT.4',
                                  'TIMEpostAV45_AVLT.5',
                                  'TIMEpostAV45_AVLT.6',
                                  'TIMEpostAV45_AVLT.7',
                                  'TIMEpostAV45_AVLT.8',
                                  'TIMEpostAV45_AVLT.9',
                                  'TIMEpostAV45_AVLT.10',
                                  'TIMEpostAV45_AVLT.11',
                                  'AVLT_3MTHS_AV45',
                                  'AVLT_post_AV45_followuptime',
                                  'AVLT_3MTHSAV45_Date',
                                  'AVLT_AV45_2_3MTHS',
                                  'AVLT_AV45_2_DATE',
                                  'AVLT_slope_all',
                                  'AVLTslope_postAV45',
                                  'GD_TOTAL.1',
                                  'GD_TOTAL.2',
                                  'GD_TOTAL.3',
                                  'GD_TOTAL.4',
                                  'GD_TOTAL.5',
                                  'GD_TOTAL.6',
                                  'GD_TOTAL.7',
                                  'GD_TOTAL.8',
                                  'GD_TOTAL.9',
                                  'GD_TOTAL.10',
                                  'GD_timePostAV45.1',
                                  'GD_timePostAV45.2',
                                  'GD_timePostAV45.3',
                                  'GD_timePostAV45.4',
                                  'GD_timePostAV45.5',
                                  'GD_timePostAV45.6',
                                  'GD_timePostAV45.7',
                                  'GD_timePostAV45.8',
                                  'GD_timePostAV45.9',
                                  'GD_timePostAV45.10',
                                  'GD_AV45_6MTHS',
                                  'GD_AV45_DATE',
                                  'GD_AV45_2_6MTHS',
                                  'GD_AV45_2_DATE',
                                  'GD_AV45_3_6MTHS',
                                  'GD_AV45_3_DATE',
                                  'GD_slope'],
                    'FDG': ['FDGvisit_count',
                            'FDG_Bin_Baseline',
                            'FDG_pons_Vis1',
                            'FDG_pons_Vis2',
                            'FDG_pons_Vis3',
                            'FDG_pons_Vis4',
                            'FDG_pons_Vis5',
                            'FDG_pons_Vis6',
                            'FDG_pons_Vis7',
                            'FDG_pons_Vis8',
                            'FDG_pons_Vis9',
                            'FDG_pons_Vis10',
                            'FDGVis1_ReltoAV45',
                            'FDGVis2_ReltoAV45',
                            'FDGVis3_ReltoAV45',
                            'FDGVis4_ReltoAV45',
                            'FDGVis5_ReltoAV45',
                            'FDGVis6_ReltoAV45',
                            'FDGVis7_ReltoAV45',
                            'FDGVis8_ReltoAV45',
                            'FDGVis9_ReltoAV45',
                            'FDGVis10_ReltoAV45',
                            'FDGVis1_Date',
                            'FDGVis2_Date',
                            'FDGVis3_Date',
                            'FDGVis4_Date',
                            'FDGVis5_Date',
                            'FDGVis6_Date',
                            'FDGVis7_Date',
                            'FDGVis8_Date',
                            'FDGVis9_Date',
                            'FDGVis10_Date',
                            'FDG_postAV45_slope',
                            'FDG_postAV45_followuptime',
                            'FDG_wholebrain_Vis1',
                            'FDG_wholebrain_Vis2',
                            'FDG_wholebrain_Vis3',
                            'FDG_wholebrain_Vis4',
                            'FDG_wholebrain_Vis5',
                            'FDG_wholebrain_Vis6',
                            'FDG_wholebrain_Vis7',
                            'FDG_wholebrain_Vis8',
                            'FDG_wholebrain_Vis9',
                            'FDG_wholebrain_Vis10',
                            'FDG_Vis1_Vis2_Time',
                            'AGE_FDG1',
                            'FDG_3mths_AV45_WBRAIN_date',
                            'FDG_3mths_AV45_WBRAIN',
                            'FDG_3months_AV45',
                            'MEAN_FDG_3MONTHS_pons',
                            'FDG_AV45_BIN1.21_pons',
                            'FDG_AV45_BIN1.25_pons',
                            'FDG_3mths_AV45_2',
                            'FDG_number',
                            'MEAN_FDG_AV45_2',
                            'MEAN_FDG_AV45_2_BIN1.21'],
                    'CSF': ['CSF_diff1',
                            'CSF_diff2',
                            'CSF_diff3',
                            'CSF_diff4',
                            'CSF_diff5',
                            'CSF_diff6',
                            'CSF_diff7',
                            'ClosestCSF_AV45.1',
                            'diff_CSF_AV45.1',
                            'AV45CSF<3months.1',
                            'abeta_closest_AV45.1',
                            'abeta_bin192.1',
                            'Closest_AV45_sec',
                            'CSF_diff1.1',
                            'CSF_diff2.1',
                            'CSF_diff3.1',
                            'CSF_diff4.1',
                            'CSF_diff5.1',
                            'CSF_diff6.1',
                            'CSF_diff7.1',
                            'CSF_ABETA.1',
                            'CSF_ABETA.2',
                            'CSF_ABETA.3',
                            'CSF_ABETA.4',
                            'CSF_ABETA.5',
                            'CSF_ABETA.6',
                            'CSF_ABETA.7',
                            'CSF_ABETApostAV45.1',
                            'CSF_ABETApostAV45.2',
                            'CSF_ABETApostAV45.3',
                            'CSF_ABETApostAV45.4',
                            'CSF_ABETApostAV45.5',
                            'CSF_ABETApostAV45.6',
                            'CSF_ABETApostAV45.7',
                            'CSF_ABETA_slope',
                            'CSF_ABETA_closest_AV45',
                            'CSF_ABETA_closest_AV45_2',
                            'CSF_ABETA_closest_AV45_3',
                            'CSF_ABETA_closest_AV45_BIN_192',
                            'CSF_ABETA_closest_AV45_2_BIN_192',
                            'CSF_ABETA_closest_AV45_3_BIN_192',
                            'CSF_TAU.1',
                            'CSF_TAU.2',
                            'CSF_TAU.3',
                            'CSF_TAU.4',
                            'CSF_TAU.5',
                            'CSF_TAU.6',
                            'CSF_TAU.7',
                            'CSF_TAUpostAV45.1',
                            'CSF_TAUpostAV45.2',
                            'CSF_TAUpostAV45.3',
                            'CSF_TAUpostAV45.4',
                            'CSF_TAUpostAV45.5',
                            'CSF_TAUpostAV45.6',
                            'CSF_TAUpostAV45.7',
                            'CSF_TAU_slope',
                            'CSF_TAU_closest_AV45',
                            'CSF_TAU_closest_AV45_2',
                            'CSF_TAU_closest_AV45_3',
                            'CSF_TAU_closest_AV45_BIN_93',
                            'CSF_TAU_closest_AV45_2_BIN_93',
                            'CSF_TAU_closest_AV45_3_BIN_93',
                            'CSF_PTAU.1',
                            'CSF_PTAU.2',
                            'CSF_PTAU.3',
                            'CSF_PTAU.4',
                            'CSF_PTAU.5',
                            'CSF_PTAU.6',
                            'CSF_PTAU.7',
                            'CSF_PTAUpostAV45.1',
                            'CSF_PTAUpostAV45.2',
                            'CSF_PTAUpostAV45.3',
                            'CSF_PTAUpostAV45.4',
                            'CSF_PTAUpostAV45.5',
                            'CSF_PTAUpostAV45.6',
                            'CSF_PTAUpostAV45.7',
                            'CSF_PTAU_slope',
                            'CSF_PTAU_closest_AV45',
                            'CSF_PTAU_closest_AV45_2',
                            'CSF_PTAU_closest_AV45_3',
                            'CSF_PTAU_closest_AV45_BIN_23',
                            'CSF_PTAU_closest_AV45_2_BIN_23',
                            'CSF_PTAU_closest_AV45_3_BIN_23'],
                'CLINICAL': ['MOD_HACHINSKI',
                             'Withdrawal (1=full,2=partial,3=stayedthroughADNI1)',
                             'Withdrawal_date',
                             'AD_Medication',
                             'Depression_Medication',
                             'Hypertension',
                             'Cholesterol',
                             'CDR_sum_boxes',
                             'Gender',
                             'Edu.(Yrs)',
                             'Handedness',
                             'APOE2_BIN',
                             'APOE4_BIN',
                             'APOE4_NUM',
                             'Init_Diagnosis',
                             'Closest_DX_Jun15',
                             'DX_Jun15_closestdate',
                             'MCItoADConv(fromav45)',
                             'MCItoADConvDate',
                             'FollowupTimetoDX',
                             'Baseline_date',
                             'AV45_MCItoAD_ConvTime',
                             'Baseline_MCItoAD_ConvTime',
                             'BD MM-YY'],
                'MRI': ['Hippocampal_Volume_3months',
                        'Left_Volume',
                        'Right_Volume',
                        'ICV',
                        'Hippocampal_Volume_normalizd',
                        'HippVol_norm_BIN0.0047',
                        'ICV_pchange(yrs)',
                        'LatVent_Diff',
                        'FrontalVol_diff',
                        'CingulateVol_diff',
                        'ParietalVol_Diff',
                        'TemporalVol_diff',
                        'LatVent_pchange',
                        'FrontalVol_pchange',
                        'CingulateVol_pchange',
                        'ParietalVol_pchange',
                        'TemporalVol_pchange',
                        'MRI_3ths_AV45',
                        'MRI_3mths_AV45_2',
                        'DTI_AV45_Closest_CC',
                        'ASL_CC_closest',
                        'TBMSyn_DATE.1',
                        'TBMSyn_DATE.2',
                        'TBMSyn_DATE.3',
                        'TBMSyn_DATE.4',
                        'TBMSyn_DATE.5',
                        'TBMSyn_DATE.6',
                        'TBMSyn_DATE.7',
                        'TBMSyn_DATE.8',
                        'TBMSyn_DATE.9',
                        'TBMSyn_DATE.10',
                        'TBMSyn_SCORE.1',
                        'TBMSyn_SCORE.2',
                        'TBMSyn_SCORE.3',
                        'TBMSyn_SCORE.4',
                        'TBMSyn_SCORE.5',
                        'TBMSyn_SCORE.6',
                        'TBMSyn_SCORE.7',
                        'TBMSyn_SCORE.8',
                        'TBMSyn_SCORE.9',
                        'TBMSyn_SCORE.10',
                        'TBMSyn_postAV45.1',
                        'TBMSyn_postAV45.2',
                        'TBMSyn_postAV45.3',
                        'TBMSyn_postAV45.4',
                        'TBMSyn_postAV45.5',
                        'TBMSyn_postAV45.6',
                        'TBMSyn_postAV45.7',
                        'TBMSyn_postAV45.8',
                        'TBMSyn_postAV45.9',
                        'TBMSyn_postAV45.10',
                        'TBMSyn_BL_DATE',
                        'TBMSyn_count',
                        'WMH_percentOfICV.1',
                        'WMH_percentOfICV.2',
                        'WMH_percentOfICV.3',
                        'WMH_percentOfICV.4',
                        'WMH_percentOfICV.5',
                        'WMH_percentOfICV.6',
                        'WMH_percentOfICV.7',
                        'WMH_postAV45.1',
                        'WMH_postAV45.2',
                        'WMH_postAV45.3',
                        'WMH_postAV45.4',
                        'WMH_postAV45.5',
                        'WMH_postAV45.6',
                        'WMH_postAV45.7',
                        'WMH_WHITMATHYP.1',
                        'WMH_WHITMATHYP.2',
                        'WMH_WHITMATHYP.3',
                        'WMH_WHITMATHYP.4',
                        'WMH_WHITMATHYP.5',
                        'WMH_WHITMATHYP.6',
                        'WMH_WHITMATHYP.7',
                        'WMH_slope'],
                'AV45_NONPVC': ['AV45_WM70/composite',
                                'AV45_2_WM70/composite',
                                'AV45_3_WM70/composite',
                                'AV45_2_WM70/composite_pchange',
                                'AV45_3_WM70/composite_pchange',
                                'AV45_2_WM70/composite_pchange_ABS',
                                'AV45_3_WM70/composite_pchange_ABS',
                                'AV45_2_WM70/composite_diff',
                                'AV45_3_WM70/composite_diff',
                                'AV45_2_WM70/composite_diff_ABS',
                                'AV45_3_WM70/composite_diff_ABS',
                                'AV45_WM70/composite_Slope',
                                'AV45_WM70/cerebg',
                                'AV45_2_WM70/cerebg',
                                'AV45_3_WM70/cerebg',
                                'AV45_2_WM70/cerebg_pchange',
                                'AV45_3_WM70/cerebg_pchange',
                                'AV45_2_WM70/cerebg_pchange_ABS',
                                'AV45_3_WM70/cerebg_pchange_ABS',
                                'AV45_2_WM70/cerebg_diff',
                                'AV45_3_WM70/cerebg_diff',
                                'AV45_2_WM70/cerebg_diff_ABS',
                                'AV45_3_WM70/cerebg_diff_ABS',
                                'AV45_WM70/cerebg_Slope',
                                'AV45_WM70/wcereb',
                                'AV45_2_WM70/wcereb',
                                'AV45_3_WM70/wcereb',
                                'AV45_2_WM70/wcereb_pchange',
                                'AV45_3_WM70/wcereb_pchange',
                                'AV45_2_WM70/wcereb_pchange_ABS',
                                'AV45_3_WM70/wcereb_pchange_ABS',
                                'AV45_2_WM70/wcereb_diff',
                                'AV45_3_WM70/wcereb_diff',
                                'AV45_2_WM70/wcereb_diff_ABS',
                                'AV45_3_WM70/wcereb_diff_ABS',
                                'AV45_WM70/wcereb_Slope',
                                'AV45_LeftPutamen/WM70',
                                'AV45_RightPutamen/WM70',
                                'AV45_LeftCaudate/WM70',
                                'AV45_RightCaudate/WM70',
                                'AV45_LeftPallidum/WM70',
                                'AV45_RightPallidum/WM70',
                                'AV45_BigRef',
                                'AV45_2_BigRef',
                                'AV45_3_BigRef',
                                'AV45_2_BigRef_pchange',
                                'AV45_3_BigRef_pchange',
                                'AV45_2_BigRef_pchange_ABS',
                                'AV45_3_BigRef_pchange_ABS',
                                'AV45_2_BigRef_diff',
                                'AV45_3_BigRef_diff',
                                'AV45_2_BigRef_diff_ABS',
                                'AV45_3_BigRef_diff_ABS',
                                'AV45_BigRef_BIN.79',
                                'AV45_2_BigRef_BIN.79',
                                'AV45_3_BigRef_BIN.79',
                                'AV45_BigRef_Slope',
                                'AV45_WM70',
                                'AV45_2_WM70',
                                'AV45_3_WM70',
                                'AV45_2_WM70_pchange',
                                'AV45_3_WM70_pchange',
                                'AV45_2_WM70_pchange_ABS',
                                'AV45_3_WM70_pchange_ABS',
                                'AV45_2_WM70_diff',
                                'AV45_3_WM70_diff',
                                'AV45_2_WM70_diff_ABS',
                                'AV45_3_WM70_diff_ABS',
                                'AV45_WM70_BIN.62',
                                'AV45_2_WM70_BIN.62',
                                'AV45_3_WM70_BIN.62',
                                'AV45_WM70_Slope',
                                'AV45_cerebg',
                                'AV45_2_cerebg',
                                'AV45_3_cerebg',
                                'AV45_2_cerebg_pchange',
                                'AV45_3_cerebg_pchange',
                                'AV45_2_cerebg_pchange_ABS',
                                'AV45_3_cerebg_pchange_ABS',
                                'AV45_2_cerebg_diff',
                                'AV45_3_cerebg_diff',
                                'AV45_2_cerebg_diff_ABS',
                                'AV45_3_cerebg_diff_ABS',
                                'AV45_cerebg_BIN1.26',
                                'AV45_2_cerebg_BIN1.26',
                                'AV45_3_cerebg_BIN1.26',
                                'AV45_cerebg_Slope',
                                'AV45_wcereb',
                                'AV45_2_wcereb',
                                'AV45_3_wcereb',
                                'AV45_2_wcereb_pchange',
                                'AV45_3_wcereb_pchange',
                                'AV45_2_wcereb_pchange_ABS',
                                'AV45_3_wcereb_pchange_ABS',
                                'AV45_2_wcereb_diff',
                                'AV45_3_wcereb_diff',
                                'AV45_2_wcereb_diff_ABS',
                                'AV45_3_wcereb_diff_ABS',
                                'AV45_wcereb_BIN1.11',
                                'AV45_2_wcereb_BIN1.11',
                                'AV45_3_wcereb_BIN1.11',
                                'AV45_wcereb_Slope',
                                'AV45_brainstem',
                                'AV45_2_brainstem',
                                'AV45_3_brainstem',
                                'AV45_2_brainstem_pchange',
                                'AV45_3_brainstem_pchange',
                                'AV45_2_brainstem_pchange_ABS',
                                'AV45_3_brainstem_pchange_ABS',
                                'AV45_2_brainstem_diff',
                                'AV45_3_brainstem_diff',
                                'AV45_2_brainstem_diff_ABS',
                                'AV45_3_brainstem_diff_ABS',
                                'AV45_brainstem_BIN.79',
                                'AV45_2_brainstem_BIN.79',
                                'AV45_3_brainstem_BIN.79',
                                'AV45_brainstem_Slope',
                                'AV45_WMratio',
                                'AV45_2_WMratio',
                                'AV45_3_WMratio',
                                'AV45_2_WMratio_pchange',
                                'AV45_3_WMratio_pchange',
                                'AV45_2_WMratio_pchange_ABS',
                                'AV45_3_WMratio_pchange_ABS',
                                'AV45_2_WMratio_diff',
                                'AV45_3_WMratio_diff',
                                'AV45_2_WMratio_diff_ABS',
                                'AV45_3_WMratio_diff_ABS',
                                'AV45_Frontal/BigRef',
                                'AV45_2_Frontal/BigRef',
                                'AV45_3_Frontal/BigRef',
                                'AV45_2_Frontal/BigRef_pchange',
                                'AV45_3_Frontal/BigRef_pchange',
                                'AV45_2_Frontal/BigRef_diff',
                                'AV45_3_Frontal/BigRef_diff',
                                'AV45_Cingulate/BigRef',
                                'AV45_2_Cingulate/BigRef',
                                'AV45_3_Cingulate/BigRef',
                                'AV45_2_Cingulate/BigRef_pchange',
                                'AV45_3_Cingulate/BigRef_pchange',
                                'AV45_2_Cingulate/BigRef_diff',
                                'AV45_3_Cingulate/BigRef_diff',
                                'AV45_Parietal/BigRef',
                                'AV45_2_Parietal/BigRef',
                                'AV45_3_Parietal/BigRef',
                                'AV45_2_Parietal/BigRef_pchange',
                                'AV45_3_Parietal/BigRef_pchange',
                                'AV45_2_Parietal/BigRef_diff',
                                'AV45_3_Parietal/BigRef_diff',
                                'AV45_Temporal/BigRef',
                                'AV45_2_Temporal/BigRef',
                                'AV45_3_Temporal/BigRef',
                                'AV45_2_Temporal/BigRef_pchange',
                                'AV45_3_Temporal/BigRef_pchange',
                                'AV45_2_Temporal/BigRef_diff',
                                'AV45_3_Temporal/BigRef_diff',
                                'frontal_asymmetry',
                                'cingulate_asymmetry',
                                'parietal_asymmetry',
                                'temporal_asymmetry',
                                'AV45_WM70/composite_Slope_2pts',
                                'AV45_WM70/composite_Slope_3pts',
                                'AV45_WM70/cerebg_Slope_2pts',
                                'AV45_WM70/cerebg_Slope_3pts',
                                'AV45_WM70/wcereb_Slope_2pts',
                                'AV45_WM70/wcereb_Slope_3pts',
                                'AV45_BigRef_Slope_2pts',
                                'AV45_BigRef_Slope_3pts',
                                'AV45_WM70_Slope_2pts',
                                'AV45_WM70_Slope_3pts',
                                'AV45_cerebg_Slope_2pts',
                                'AV45_cerebg_Slope_3pts',
                                'AV45_wcereb_Slope_2pts',
                                'AV45_wcereb_Slope_3pts',
                                'AV45_brainstem_Slope_2pts',
                                'AV45_brainstem_Slope_3pts'],
                'AV45_PVC': ['AV45_PVC_agghigh_CorticalSummary/Wholecereb_BL',
                            'AV45_PVC_agghigh_CorticalSummary/Wholecereb_Scan2',
                            'AV45_PVC_agghigh_CorticalSummary/Wholecereb_Scan3',
                            'AV45_PVC_agghigh_CorticalSummary/Bigref_BL',
                            'AV45_PVC_agghigh_CorticalSummary/Bigref_Scan2',
                            'AV45_PVC_agghigh_CorticalSummary/Bigref_Scan3',
                            'AV45_PVC_agghigh_CorticalSummary_slope',
                            'AV45_PVC_agghigh_CerebGM_BL',
                            'AV45_PVC_agghigh_CerebGM_Scan2',
                            'AV45_PVC_agghigh_CerebGM_Scan3',
                            'AV45_PVC_agghigh_HemiWM/Wholecereb_BL',
                            'AV45_PVC_agghigh_HemiWM/Wholecereb_Scan2',
                            'AV45_PVC_agghigh_HemiWM/Wholecereb_Scan3',
                            'AV45_PVC_grp2_slope',
                            'AV45_PVC_agghigh_Wholecereb_BL',
                            'AV45_PVC_agghigh_Wholecereb_Scan2',
                            'AV45_PVC_agghigh_Wholecereb_Scan3',
                            'AV45_PVC_agghigh_CerebWM_BL',
                            'AV45_PVC_agghigh_CerebWM_Scan2',
                            'AV45_PVC_agghigh_CerebWM_Scan3',
                            'AV45_PVC_group2_CorticalSummary/Wholecereb_BL',
                            'AV45_PVC_group2_CorticalSummary/Wholecereb_Scan2',
                            'AV45_PVC_group2_CorticalSummary/Wholecereb_Scan3',
                            'AV45_PVC_group2_CorticalSummary/Bigref_BL',
                            'AV45_PVC_group2_CorticalSummary/Bigref_Scan2',
                            'AV45_PVC_group2_CorticalSummary/Bigref_Scan3',
                            'AV45_PVC_group2_CorticalSummary_slope',
                            'AV45_PVC_group2_CorticalSummary_slope_allposs',
                            'AV45_PVC_group2_CerebGM_BL',
                            'AV45_PVC_group2_CerebGM_Scan2',
                            'AV45_PVC_group2_CerebGM_Scan3',
                            'AV45_PVC_group2_HemiWM/Wholecereb_BL',
                            'AV45_PVC_group2_HemiWM/Wholecereb_Scan2',
                            'AV45_PVC_group2_HemiWM/Wholecereb_Scan3',
                            'AV45_PVC_group2_Wholecereb_BL',
                            'AV45_PVC_group2_Wholecereb_Scan2',
                            'AV45_PVC_group2_Wholecereb_Scan3',
                            'AV45_PVC_group2_CerebWM_BL',
                            'AV45_PVC_group2_CerebWM_Scan2',
                            'AV45_PVC_group2_CerebWM_Scan3',
                            'AV45_PVC_group4_CorticalSummary/Wholecereb_BL',
                            'AV45_PVC_group4_CorticalSummary/Wholecereb_Scan2',
                            'AV45_PVC_group4_CorticalSummary/Wholecereb_Scan3',
                            'AV45_PVC_group4_CorticalSummary/Bigref_BL',
                            'AV45_PVC_group4_CorticalSummary/Bigref_Scan2',
                            'AV45_PVC_group4_CorticalSummary/Bigref_Scan3',
                            'AV45_PVC_group4_CorticalSummary_slope',
                            'AV45_PVC_group4_CorticalSummary_slope_allposs',
                            'AV45_PVC_group4_CerebGM_BL',
                            'AV45_PVC_group4_CerebGM_Scan2',
                            'AV45_PVC_group4_CerebGM_Scan3',
                            'AV45_PVC_group4_HemiWM/Wholecereb_BL',
                            'AV45_PVC_group4_HemiWM/Wholecereb_Scan2',
                            'AV45_PVC_group4_HemiWM/Wholecereb_Scan3',
                            'AV45_PVC_group4_Wholecereb_BL',
                            'AV45_PVC_group4_Wholecereb_Scan2',
                            'AV45_PVC_group4_Wholecereb_Scan3',
                            'AV45_PVC_group4_CerebWM_BL',
                            'AV45_PVC_group4_CerebWM_Scan2',
                            'AV45_PVC_group4_CerebWM_Scan3',
                            'AV45_PVC_agghigh_CorticalSummary_slope_2points',
                            'AV45_PVC_agghigh_CorticalSummary_slope_3points',
                            'AV45_PVC_agghigh_CorticalSummary/Bigref_slope_2points',
                            'AV45_PVC_agghigh_CorticalSummary/Bigref_slope_3points',
                            'AV45_PVC_agghigh_HemiWM_slope_2points',
                            'AV45_PVC_agghigh_HemiWM_slope_3points',
                            'AV45_PVC_group2_CorticalSummary_slope_2points',
                            'AV45_PVC_group2_CorticalSummary_slope_3points',
                            'AV45_PVC_group2_CorticalSummary/Bigref_slope_2points',
                            'AV45_PVC_group2_CorticalSummary/Bigref_slope_3points',
                            'AV45_PVC_group2_HemiWM_slope_2points',
                            'AV45_PVC_group2_HemiWM_slope_3points',
                            'AV45_PVC_group4_CorticalSummary_slope_2points',
                            'AV45_PVC_group4_CorticalSummary_slope_3points',
                            'AV45_PVC_group4_CorticalSummary/Bigref_slope_2points',
                            'AV45_PVC_group4_CorticalSummary/Bigref_slope_3points',
                            'AV45_PVC_group4_HemiWM_slope_2points',
                            'AV45_PVC_group4_HemiWM_slope_3points',
                            'AV45_PVC_agghigh_accumulator.00790',
                            'AV45_PVC_agghigh_wcereb1.27_BL'], 
                'RID': ['RID']}




def syncRoussetResults(old_headers, old_lines, rousset_matfile, timepoint, dump_to=None):
    assert timepoint in set(['BL', 'Scan2', 'Scan3'])

    group_names, sorted_data = importRoussetResults(rousset_matfile)

    valid_timepoints = ['BL', 'Scan2', 'Scan3']
    valid_groupings = ['agghigh', 'group2', 'group4']
    to_add_headers = []
    for vg in valid_groupings:
        to_add_headers += ['AV45_PVC_%s_CorticalSummary/Wholecereb_%s' % (vg, tp) for tp in valid_timepoints]
        to_add_headers += ['AV45_PVC_%s_CorticalSummary_slope_2points' % (vg), 'AV45_PVC_%s_CorticalSummary_slope_3points' % (vg)]
        to_add_headers += ['AV45_PVC_%s_CorticalSummary/Bigref_%s' % (vg, tp) for tp in valid_timepoints]
        to_add_headers += ['AV45_PVC_%s_CorticalSummary/Bigref_slope_2points' % (vg), 'AV45_PVC_%s_CorticalSummary/Bigref_slope_3points' % (vg)]
        to_add_headers += ['AV45_PVC_%s_HemiWM/Wholecereb_%s' % (vg, tp) for tp in valid_timepoints]
        to_add_headers += ['AV45_PVC_%s_HemiWM_slope_2points' % (vg), 'AV45_PVC_%s_HemiWM_slope_3points' % (vg)]
        to_add_headers += ['AV45_PVC_%s_Wholecereb_%s' % (vg, tp) for tp in valid_timepoints]
        to_add_headers += ['AV45_PVC_%s_CerebGM_%s' % (vg, tp) for tp in valid_timepoints]
        to_add_headers += ['AV45_PVC_%s_CerebWM_%s' % (vg, tp) for tp in valid_timepoints]
        if vg == 'agghigh':
            # SOME HARDCODED THRESHOLDING
            to_add_headers += ['AV45_PVC_agghigh_accumulator.00790', 'AV45_PVC_agghigh_wcereb1.27_BL']
    after = old_headers[max(i for i,_ in enumerate(old_headers) if _.startswith('AV45_') and 'PVC' not in _)] # last element that contains 'AV45'
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=after)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        valid_groupings = ['agghigh', 'group2', 'group4']
        new_subj_data = {}
        for k in subj_row.keys():
            if k in valid_groupings:
                v = subj_row[k]
                wholecereb = float(v['wholecereb']['pvcval'])
                bigref = float(v['bigref']['pvcval'])
                for region, values in v.iteritems():
                    pvcval = float(values.get('pvcval',0.0))
                    if region == 'composite':
                        new_subj_data['AV45_PVC_%s_CorticalSummary/Wholecereb' % k] = pvcval/wholecereb
                        new_subj_data['AV45_PVC_%s_CorticalSummary/Bigref' % k] = pvcval/bigref
                    elif region == 'cerebGM':
                        new_subj_data['AV45_PVC_%s_CerebGM' % k] = pvcval
                    elif region == 'hemiWM':
                        new_subj_data['AV45_PVC_%s_HemiWM/Wholecereb' % k] = pvcval/wholecereb
                    elif region == 'wholecereb':
                        new_subj_data['AV45_PVC_%s_Wholecereb' % k] = pvcval
                    elif region == 'cerebWM':
                        new_subj_data['AV45_PVC_%s_CerebWM' % k] = pvcval
        return new_subj_data

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, sorted_data, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        new_data = {'%s_%s' % (k, timepoint): v for k,v in new_data.iteritems()}
        old_l.update(new_data)
        if timepoint == 'BL':
            bl_val = old_l.get('AV45_PVC_agghigh_CorticalSummary/Wholecereb_BL',0.0)
            old_l['AV45_PVC_agghigh_wcereb1.27_BL'] = 1 if bl_val >= 1.27 else 0
        elif timepoint == 'Scan2':
            # try to calculate slope with 2 points
            for vg in valid_groupings:
                try:
                    diff1 = float(old_l.get('AV45_1_2_Diff (Yrs)'))
                    x = [0.0, diff1]

                    val1 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Wholecereb_BL' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Wholecereb_Scan2' % vg))
                    slope, intercept, r, p, stderr = stats.linregress(x, [val1, val2])
                    old_l['AV45_PVC_%s_CorticalSummary_slope_2points' % (vg)] = slope
                    if vg == 'agghigh':
                        old_l['AV45_PVC_agghigh_accumulator.00790'] = 1 if slope >= 0.0079 else 0

                    val1 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Bigref_BL' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Bigref_Scan2' % vg))
                    slope, intercept, r, p, stderr = stats.linregress(x, [val1, val2])
                    old_l['AV45_PVC_%s_CorticalSummary/Bigref_slope_2points' % (vg)] = slope

                    val1 = float(old_l.get('AV45_PVC_%s_HemiWM/Wholecereb_BL' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_HemiWM/Wholecereb_Scan2' % vg))
                    slope, intercept, r, p, stderr = stats.linregress(x, [val1, val2])
                    old_l['AV45_PVC_%s_HemiWM_slope_2points' % (vg)] = slope
                except Exception as e:
                    old_l['AV45_PVC_%s_CorticalSummary_slope_2points' % (vg)] = ''
                    old_l['AV45_PVC_%s_CorticalSummary/Bigref_slope_2points' % (vg)] = ''
                    old_l['AV45_PVC_%s_HemiWM_slope_2points' % (vg)] = ''
        elif timepoint == 'Scan3':
            # try to calculate slope with 3 points
            for vg in valid_groupings:
                try:
                    diff1 = float(old_l.get('AV45_1_2_Diff (Yrs)'))
                    diff2 = float(old_l.get('AV45_1_3_Diff (yrs)'))
                    x = [0.0, diff1, diff2]

                    val1 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Wholecereb_BL' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Wholecereb_Scan2' % vg))
                    val3 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Wholecereb_Scan3' % vg))
                    slope, intercept, r, p, stderr = stats.linregress(x, [val1, val2, val3])
                    old_l['AV45_PVC_%s_CorticalSummary_slope_3points' % (vg)] = slope
                    if vg == 'agghigh':
                        old_l['AV45_PVC_agghigh_accumulator.00790'] = 1 if slope >= 0.0079 else 0

                    val1 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Bigref_BL' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Bigref_Scan2' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_CorticalSummary/Bigref_Scan3' % vg))
                    slope, intercept, r, p, stderr = stats.linregress(x, [val1, val2, val3])
                    old_l['AV45_PVC_%s_CorticalSummary/Bigref_slope_3points' % (vg)] = slope


                    val1 = float(old_l.get('AV45_PVC_%s_HemiWM/Wholecereb_BL' % vg))
                    val2 = float(old_l.get('AV45_PVC_%s_HemiWM/Wholecereb_Scan2' % vg))
                    val3 = float(old_l.get('AV45_PVC_%s_HemiWM/Wholecereb_Scan3' % vg))
                    slope, intercept, r, p, stderr = stats.linregress(x, [val1, val2, val3])
                    old_l['AV45_PVC_%s_HemiWM_slope_3points' % (vg)] = slope
                except Exception as e:
                    old_l['AV45_PVC_%s_CorticalSummary_slope_3points' % (vg)] = ''
                    old_l['AV45_PVC_%s_CorticalSummary/Bigref_slope_3points' % (vg)] = ''
                    old_l['AV45_PVC_%s_HemiWM_slope_3points' % (vg)] = ''
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)



def syncGDData(old_headers, old_lines, gd_file, registry_file, dump_to=None):
    registry = importRegistry(registry_file)
    gd_by_subj = importGD(gd_file, registry=registry)

    to_add_headers = ['GD_TOTAL.%s' % (i+1) for i in range(10)]
    to_add_headers += ['GD_timePostAV45.%s' % (i+1) for i in range(10)]
    to_add_headers += ['GD_AV45_6MTHS', 'GD_AV45_DATE',
                       'GD_AV45_2_6MTHS', 'GD_AV45_2_DATE',
                       'GD_AV45_3_6MTHS', 'GD_AV45_3_DATE']
    to_add_headers.append('GD_slope')
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        scores = sorted(subj_row, key=lambda x: x['EXAMDATE'])
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        new_subj_data = {}
        six_months = 31 * 6

        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToAV45(scores, bl_av45, av45_2, av45_3, day_limit=six_months)
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
        for i in range(10):
            if i >= len(scores):
                post_diff = ''
                tot_score = ''
            else:
                cur_score = scores[i]
                date = cur_score['EXAMDATE']
                tot_score = cur_score['GDTOTAL']
                try:
                    tot_score = int(tot_score)
                except:
                    pass
                if bl_av45 is not None:
                    post_diff_days = (date - bl_av45).days
                    post_diff_years = post_diff_days / 365.0
                    if isinstance(tot_score, int) and post_diff_days > -30:
                        post_av45_points.append((post_diff_years, tot_score))
                else:
                    post_diff_years = ''
            new_subj_data['GD_TOTAL.%s' % (i+1)] = tot_score
            new_subj_data['GD_timePostAV45.%s' % (i+1)] = post_diff_years
        # get slope
        if len(post_av45_points) >= 2:
            post_av45_points = np.array(post_av45_points)
            post_times = post_av45_points[:,0]
            post_values = post_av45_points[:,1]
            slope, intercept, r, p, stderr = stats.linregress(post_times, post_values)
            new_subj_data['GD_slope'] = slope
        else:
            new_subj_data['GD_slope'] = ''

        return new_subj_data

    new_lines = []
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, gd_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def syncADASCogData(old_headers, old_lines, adni1_adas_file, adnigo2_adas_file, registry_file, dump_to=None):
    registry = importRegistry(registry_file)
    adas_by_subj = importADASCog(adni1_adas_file, adnigo2_adas_file, registry=registry)

    to_add_headers = ['ADAS_AV45_3_3MTHS', 'ADAS_AV45_3_DATE']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after='ADAS_AV45_2_DATE')
    to_add_headers = ['ADAS_post_AV45_followuptime', 'ADASslope_postAV45']
    new_headers = rearrangeHeaders(new_headers, to_add_headers, after='ADAS_3MTH_AV45')

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        tests = sorted(subj_row, key=lambda x : x['EXAMDATE'])
        dates = [t['EXAMDATE'] for t in tests]
        if len(set(dates)) != len(dates):
            print "%s has Dup dates, skipping: %s" % (subj, dates)
            return {}
        first_scan_date = dates[0]

        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)

        new_subj_data = {}
        all_values = []
        post_values = []
        all_times = []
        post_times = []
        max_followup_counter = None
        for i in range(11):
            if i >= len(tests):
                test_date = ''
                test_score = ''
                diff_from_first = ''
            else:
                test_results = tests[i]
                test_date = test_results['EXAMDATE']
                test_score = test_results['TOTSCORE']
                if test_score != '':
                    test_score = float(test_score)
                diff_from_first = (test_date-first_scan_date).days / 365.0
                if test_score != '':
                    all_values.append(test_score)
                    all_times.append(diff_from_first)
            count = i+1
            new_subj_data['ADAScog_DATE%s' % count] = test_date
            new_subj_data['ADAScog.%s' % count] = test_score
            new_subj_data['TIME_ADAS.%s' % count] = diff_from_first
            if bl_av45 is not None and test_date != '':
                rel_time_days = (test_date - bl_av45).days
                rel_time = round(rel_time_days / 365.0, 2)
                new_subj_data['TIMEreltoAV45_ADAS.%s' % count] = rel_time
                if abs(rel_time_days) <= 93 and 'ADAS_3MTH_AV45' not in new_subj_data:
                    new_subj_data['ADAS_3MTH_AV45'] = test_score
                    new_subj_data['ADAS_3MTHS_AV45DATE'] = test_date
                if rel_time >= (-93.0/365.0):
                    if test_score != '':
                        post_values.append(test_score)
                        post_times.append(diff_from_first)
                        max_followup_counter = rel_time
                    new_subj_data['TIMEpostAV45_ADAS.%s' % count] = rel_time
                else:
                    new_subj_data['TIMEpostAV45_ADAS.%s' % count] = ''
            if av45_2 is not None and test_date != '':
                rel_time_days = (test_date - av45_2).days 
                if abs(rel_time_days) <= 93 and 'ADAS_AV45_2_3MTHS' not in new_subj_data:
                    new_subj_data['ADAS_AV45_2_3MTHS'] = test_score
                    new_subj_data['ADAS_AV45_2_DATE'] = test_date
            if av45_3 is not None and test_date != '':
                rel_time_days = (test_date - av45_3).days 
                if abs(rel_time_days) <= 93 and 'ADAS_AV45_3_3MTHS' not in new_subj_data:
                    new_subj_data['ADAS_AV45_3_3MTHS'] = test_score
                    new_subj_data['ADAS_AV45_3_DATE'] = test_date
        if max_followup_counter is not None:
            new_subj_data['ADAS_post_AV45_followuptime'] = max(max_followup_counter, 0.0)
        # fill in the blanks
        fill_in = ['ADAS_3MTH_AV45', 'ADAS_3MTHS_AV45DATE', 
                   'ADAS_AV45_2_3MTHS','ADAS_AV45_2_DATE',
                   'ADAS_AV45_3_3MTHS','ADAS_AV45_3_DATE',
                   'ADAS_post_AV45_followuptime']
        for f_key in fill_in:
            if f_key not in new_subj_data:
                new_subj_data[f_key] = ''
        
        # get slopes
        '''
        if len(all_values) >= 2:
            slope, intercept, r, p, stderr = stats.linregress(all_times, all_values)
            new_subj_data['ADAS_slope_all'] = slope
        else:
            new_subj_data['ADAS_slope_all'] = ''
            '''
        if len(post_values) >= 2:
            slope, intercept, r, p, stderr = stats.linregress(post_times, post_values)
            new_subj_data['ADASslope_postAV45'] = slope
        else:
            new_subj_data['ADASslope_postAV45'] = ''
        return new_subj_data

    new_lines = []
    new_values = 0
    total = 0
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, adas_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)

        # do comparison:
        #print "SUBJ: %s" % subj
        #print "first: %s, second: %s" % (bl_av45, av45_2)
        changed = False
        for k in sorted(new_data.keys()):
            old_value = old_l.get(k, None) # new columns potentially added
            new_value = new_data[k]
            if k.startswith('ADAScog.') and old_value == '' and new_value != '':
                new_values += 1
                changed = True
            #print "\t%s: %s -> %s" % (k, old_value, new_value)
        if changed:
            total += 1

        old_l.update(new_data)
        new_lines.append(old_l)

    print "Total subj changed: %s" % total
    print "Total new tests: %s" % new_values

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def syncMMSEData(old_headers, old_lines, mmse_file, registry_file, dump_to=None):
    registry = importRegistry(registry_file)
    mmse_by_subj = importMMSE(mmse_file, registry=registry)

    # add these new headers, if not present
    to_add_headers = ['MMSE_post_AV45_followuptime',
                      'MMSE_AV45_3MTHS',
                      'MMSE_AV45_DATE',
                      'MMSE_AV45_2_3MTHS',
                      'MMSE_AV45_2_DATE',
                      'MMSE_AV45_3_3MTHS',
                      'MMSE_AV45_3_DATE',
                      'MMSEslope_postAV45']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after='MMSCORE.11')

    def extraction_fn(subj, subj_row, old_line, patient_pets):
        tests = subj_row

        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        new_subj_data = {}
        max_followup_counter = None
        post_av45_points = []
        for i in range(11):
            if i >= len(tests):
                test_date = ''
                test_score = ''
            else:
                test_date, test_results = tests[i]
                if test_results['MMSCORE'] == '':
                    test_score = ''
                else:
                    test_score = int(test_results['MMSCORE'])
                    if test_score < 0:
                        test_score = ''

            count = i+1
            new_subj_data['MMSE_DATE%s' % count] = test_date
            new_subj_data['MMSCORE.%s' % count] = test_score

            # pair up with av45 scan date if possible
            if bl_av45 is not None and test_date != '':
                rel_time_days = (test_date - bl_av45).days 
                if abs(rel_time_days) <= 93 and 'MMSE_AV45_3MTHS' not in new_subj_data:
                    new_subj_data['MMSE_AV45_3MTHS'] = test_score
                    new_subj_data['MMSE_AV45_DATE'] = test_date
                # also try to put value into time post av45 columns
                if rel_time_days >= -93.0:
                    annualized_time = rel_time_days / 365.0
                    new_subj_data['TIMEpostAV45_MMSE.%s' % count] = annualized_time
                    if test_score != '':
                        post_av45_points.append((annualized_time, test_score))
                        max_followup_counter = annualized_time
            if ('TIMEpostAV45_MMSE.%s' % count) not in new_subj_data:
                new_subj_data['TIMEpostAV45_MMSE.%s' % count] = ''
            # pair up with subsequent av45 scans
            if av45_2 is not None and test_date != '':
                rel_time_days = (test_date - av45_2).days 
                if abs(rel_time_days) <= 93 and 'MMSE_AV45_2_3MTHS' not in new_subj_data:
                    new_subj_data['MMSE_AV45_2_3MTHS'] = test_score
                    new_subj_data['MMSE_AV45_2_DATE'] = test_date
            if av45_3 is not None and test_date != '':
                rel_time_days = (test_date - av45_3).days 
                if abs(rel_time_days) <= 93 and 'MMSE_AV45_3_3MTHS' not in new_subj_data:
                    new_subj_data['MMSE_AV45_3_3MTHS'] = test_score
                    new_subj_data['MMSE_AV45_3_DATE'] = test_date
        if max_followup_counter is not None:
            new_subj_data['MMSE_post_AV45_followuptime'] = max(max_followup_counter, 0.0)

        # get post av45 slope
        if len(post_av45_points) >= 2:
            post_times = [_[0] for _ in post_av45_points]
            post_values = [_[1] for _ in post_av45_points]
            slope, intercept, r, p, stderr = stats.linregress(post_times, post_values)
            new_subj_data['MMSEslope_postAV45'] = slope

        # fill in the blanks
        for f_key in to_add_headers:
            if f_key not in new_subj_data:
                new_subj_data[f_key] = ''
        return new_subj_data

    new_lines = []
    new_values = 0
    total = 0
    for i, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, mmse_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)

        # do comparison:
        #print "SUBJ: %s" % subj
        changed = False
        for k in sorted(new_data.keys()):
            old_value = old_l.get(k,'')
            new_value = new_data[k]
            if k.startswith('MMSCORE') and old_value == '' and new_value != '':
                new_values += 1
                changed = True
            #print "\t%s: %s -> %s" % (k, old_value, new_value)
        if changed:
            total += 1

        old_l.update(new_data)
        new_lines.append(old_l)

    print "%s subj with new tests (%s new tests)" % (total, new_values)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return new_headers, new_lines

def syncAVLTData(old_headers, old_lines, neuro_battery_file, registry_file, dump_to=None):
    registry = importRegistry(registry_file)
    avlt_by_subj = importAVLT(neuro_battery_file, registry=registry)

    # no change in headers
    to_add_headers = ['AVLT_post_AV45_followuptime']
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after='AVLT_3MTHS_AV45')

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        tests = sorted(subj_row, key=lambda x : x['EXAMDATE'])
        dates = [t['EXAMDATE'] for t in tests]
        if len(set(dates)) != len(dates):
            print "%s has Dup dates, skipping: %s" % (subj, dates)
            return {}

        first_scan_date = dates[0]
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)

        new_subj_data = {}
        all_values = []
        post_values = []
        all_times = []
        post_times = []
        max_followup_counter = None
        for i in range(11):
            if i >= len(tests):
                test_date = ''
                test_score = ''
                diff_from_first = ''
            else:
                test_results = tests[i]
                test_date = test_results['EXAMDATE']
                test_score = test_results['TOTS']
                diff_from_first = (test_date-first_scan_date).days / 365.0
                if test_score != '':
                    all_values.append(test_score) 
                    all_times.append(diff_from_first)
            count = i+1
            new_subj_data['AVLT_DATE.%s' % count] = test_date
            new_subj_data['AVLT.%s' % count] = test_score
            new_subj_data['TIME_AVLT.%s' % count] = diff_from_first
            new_subj_data['TIMEreltoAV45_AVLT.%s' % count] = ''
            if bl_av45 is not None and test_date != '':
                rel_time_days = (test_date - bl_av45).days
                rel_time = rel_time_days / 365.0
                new_subj_data['TIMEreltoAV45_AVLT.%s' % count] = rel_time
                if abs(rel_time_days) <= 93 and 'AVLT_3MTHS_AV45' not in new_subj_data:
                    new_subj_data['AVLT_3MTHS_AV45'] = test_score
                    new_subj_data['AVLT_3MTHSAV45_Date'] = test_date
                if rel_time >= (-93.0/365.0):
                    if test_score != '':
                        post_values.append(test_score)
                        post_times.append(diff_from_first)
                        max_followup_counter = rel_time
                    new_subj_data['TIMEpostAV45_AVLT.%s' % count] = rel_time
                else:
                    new_subj_data['TIMEpostAV45_AVLT.%s' % count] = ''
            if av45_2 is not None and test_date != '':
                rel_time_days = (test_date - av45_2).days 
                if abs(rel_time_days) <= 93 and 'AVLT_AV45_2_3MTHS' not in new_subj_data:
                    new_subj_data['AVLT_AV45_2_3MTHS'] = test_score
                    new_subj_data['AVLT_AV45_2_DATE'] = test_date
        if max_followup_counter is not None:
            new_subj_data['AVLT_post_AV45_followuptime'] = max(max_followup_counter, 0.0)

        # fill in the blanks
        fill_in = ['AVLT_3MTHS_AV45', 'AVLT_3MTHSAV45_Date', 
                   'AVLT_AV45_2_3MTHS','AVLT_AV45_2_DATE',
                   'AVLT_post_AV45_followuptime']
        for f_key in fill_in:
            if f_key not in new_subj_data:
                new_subj_data[f_key] = ''

        # get slopes
        if len(all_values) >= 2:
            try:
                slope, intercept, r, p, stderr = stats.linregress(all_times, all_values)
            except Exception as e:
                print "%s vs %s" % (all_times, all_values)
                raise e
            new_subj_data['AVLT_slope_all'] = slope
        else:
            new_subj_data['AVLT_slope_all'] = ''
        if len(post_values) >= 2:
            slope, intercept, r, p, stderr = stats.linregress(post_times, post_values)
            new_subj_data['AVLTslope_postAV45'] = slope
        else:
            new_subj_data['AVLTslope_postAV45'] = ''
        return new_subj_data

    new_lines = []
    new_values = 0
    total = 0
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, avlt_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None)

        # do comparison:
        #print "SUBJ: %s" % subj
        #print "first: %s, second: %s" % (bl_av45, av45_2)
        #print "first date: %s" % (first_scan_date)
        changed = False
        new_for_subj = 0
        for k in sorted(new_data.keys()):
            old_value = old_l.get(k,None)
            new_value = new_data[k]
            if k.startswith('AVLT.') and not bool(old_value) and bool(new_value):
                new_values += 1
                new_for_subj += 1
                changed = True
            #print "\t%s: %s -> %s" % (k, old_value, new_value)
        if changed:
            total += 1
        
        old_l.update(new_data)
        new_lines.append(old_l)

    print "Total subj changed: %s" % total
    print "Total new tests: %s" % new_values

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def syncDiagnosisData(old_headers, old_lines, diag_file, registry_file, demog_file, arm_file, pet_meta_file, dump_to=None):
    registry = importRegistry(registry_file)
    diag_by_subj = importADNIDiagnosis(diag_file, registry=registry)
    demogs = importDemog(demog_file)
    arm = importARM(arm_file)
    pet_meta = importPetMETA(pet_meta_file)

    pivot_date = datetime(day=14, month=10, year=2014)
    pivot_date_closest_diag_key = 'Closest_DX_Jun15'
    pivot_date_closest_date_key = 'DX_Jun15_closestdate'
    pivot_diag_regex = re.compile('^Closest_DX_')
    pivot_date_regex = re.compile('^DX_.*_closestdate$')

    # add 'Diag@AV45_3_long'
    new_headers = rearrangeHeaders(old_headers, ['Diag@AV45_3_long'], after='Diag@AV45_2_long')
    # replace old pivot date keys
    for i in range(len(new_headers)):
        if re.search(pivot_diag_regex, new_headers[i]):
            new_headers[i] = pivot_date_closest_diag_key
        elif re.search(pivot_date_regex, new_headers[i]):
            new_headers[i] = pivot_date_closest_date_key

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

        # Get AV45 Scan dates
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        patient_dob = demogs.get(subj, {}).get('dob',None)
        av45_age = ''
        av45_age2 = ''
        av45_age3 = ''
        if patient_dob is not None:
            if bl_av45:
                av45_age = (bl_av45-patient_dob).days / 365.25
            if av45_2:
                av45_age2 = (av45_2-patient_dob).days / 365.25
            if av45_3:
                av45_age3 = (av45_3-patient_dob).days / 365.25
        
        # Date differences between scans
        av45_1_2_diff = ((av45_2 - bl_av45).days/365.0) if (bl_av45 is not None and av45_2 is not None) else ''
        av45_1_3_diff = ((av45_3 - bl_av45).days/365.0) if (bl_av45 is not None and av45_3 is not None) else ''


        # find very first visit date in registry (excluding scmri)
        subj_registry = [_ for _ in registry[subj] if _['VISCODE'] != 'scmri' and _['VISCODE2'] != 'scmri']
        sorted_registry = sorted(subj_registry, key=lambda x: x['EXAMDATE'])
        baseline_registry = sorted_registry[0]

        # find pivot data point
        sorted_by_pivot = sorted(subj_diags, key=lambda x: abs(x['EXAMDATE'] - pivot_date))
        closest_to_pivot = sorted_by_pivot[0]

        new_data = {'AV45_Date': bl_av45,
                    'AV45_2_Date': av45_2,
                    'AV45_3_Date': av45_3,
                    'BD MM-YY': patient_dob,
                    'Age@AV45': av45_age,
                    'Age@AV45_2' : av45_age2,
                    'Age@AV45_3' : av45_age3,
                    'AV45_1_2_Diff (Yrs)': av45_1_2_diff,
                    'AV45_1_3_Diff (yrs)': av45_1_3_diff,
                    'MCItoADConv(fromav45)': '0',
                    'MCItoADConvDate': '',
                    'MCItoADconv_': '',
                    'AV45_MCItoAD_ConvTime': '',
                    'Baseline_MCItoAD_ConvTime': '',
                    'Diag@AV45_long': '',
                    'Diag@AV45_2_long': '',
                    'Diag@AV45_3_long': '',
                    'FollowupTimetoDX': (pivot_date - baseline_registry['EXAMDATE']).days / 365.0,
                    'Baseline': baseline_registry['EXAMDATE'],
                    'Init_Diagnosis': init_diag,
                    pivot_date_closest_diag_key: closest_to_pivot['diag'],
                    pivot_date_closest_date_key: closest_to_pivot['EXAMDATE']}
        
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToAV45(subj_diags, bl_av45, av45_2, av45_3)
        if av45_bl_closest:
            new_data['Diag@AV45_long'] = av45_bl_closest['diag']
        if av45_2_closest:
            new_data['Diag@AV45_2_long'] = av45_2_closest['diag']
        if av45_3_closest:
            new_data['Diag@AV45_3_long'] = av45_3_closest['diag'] 

        for i, diag_row in enumerate(subj_diags):
            # parse change based on init diag
            if diag_row['change'] == 5:
                new_data['MCItoADConv(fromav45)'] = '1'
                new_data['MCItoADConvDate'] = new_data['MCItoADconv_'] = diag_row['EXAMDATE']
                new_data['Baseline_MCItoAD_ConvTime'] = (diag_row['EXAMDATE'] - new_data['Baseline']).days / 365.0
                if bl_av45 is not None:
                    new_data['AV45_MCItoAD_ConvTime'] = (diag_row['EXAMDATE'] - bl_av45).days / 365.0
        return new_data

    new_lines = []
    new_values = 0
    new_conversions = 0
    total = 0
    for linenum, old_l in enumerate(old_lines):
        new_data = updateLine(old_l, diag_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=pet_meta)

        # compare
        '''
        print "SUBJ: %s" % subj
        print "ARM: %s" % arm.get(subj,[])
        print "First AV: %s" % bl_av45
        print "Second AV: %s" % av45_2
        print "Diags: %s" % ([(str(_['EXAMDATE']),_['change'],_['diag']) for _ in subj_diags])
        print 'Reg Dates: %s' % sorted_registry
        '''
        changed = False
        ignore_change = set([pivot_date_closest_diag_key, 
                             pivot_date_closest_date_key, 
                             'FollowupTimetoDX', 
                             'MCItoADConvDate', 
                             'MCItoADconv_',
                             'AV45_MCItoAD_ConvTime',
                             'MCItoADConv(fromav45)',
                             'Baseline_MCItoAD_ConvTime',
                             'AV45_Date',
                             'AV45_2_Date'])
        
        for k in sorted(new_data.keys()):
            old_value = old_l.get(k)
            new_value = new_data.get(k)
            if k not in ignore_change and old_value != new_value:
                changed = True
            if k == 'MCItoADConv(fromav45)' and old_value == '0' and new_value == '1':
                new_conversions += 1

        if changed:
            for k in sorted(new_data.keys()):
                if k not in ignore_change:
                    old_value = old_l.get(k)
                    new_value = new_data.get(k)
                    #print "\t%s: %s -> %s" % (k, old_value, new_value)
            new_values += 1

        old_l.update(new_data)
        new_lines.append(old_l)

    print "Changed: %s" % new_values
    print "Total: %s" % len(new_lines)
    print "New conversions: %s" % new_conversions

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def syncAV45Data(old_headers, old_lines, av45_file, av45_nontp_file, registry_file, dump_to=None):
    '''
    This function does not follow the pattern of the other sync functions because
    the header/data update is accomplished in a nested function, and it also allows
    for the possibility of introducing new subjects via the dataset being synced in
    '''
    registry = importRegistry(registry_file)
    av45_by_subj = importAV45(av45_file, av45_nontp_file=av45_nontp_file, registry=registry)
    new_headers = None
    new_lines = []
    old_subjects = set([])
    for old_l in old_lines:
        # get subject ID
        try:
            subj = int(old_l['RID'])
        except Exception as e:
            continue
        old_subjects.add(subj)
        if subj not in av45_by_subj:
            print "No AV45 data for %s" % subj
            new_lines.append(old_l)
            continue

        updated_headers, new_data = parseAV45Entries(old_headers, av45_by_subj[subj])
        if new_headers is None:
            new_headers = updated_headers
        new_data = convertToCSVDataType(new_data, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # find new subjects - in av45 file but not in old master csv file
    print len(set(av45_by_subj.keys()))
    print len(old_subjects)
    new_subjects = list(set(av45_by_subj.keys()) - old_subjects)
    print new_subjects
    for ns in new_subjects:
        new_subj_row = {k: '' for k in new_headers}
        new_subj_row['RID'] = str(ns)
        updated_headers, new_columns = parseAV45Entries(old_headers, av45_by_subj[ns])
        new_subj_row.update(new_columns)
        new_lines.append(new_subj_row)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def parseAV45Entries(old_headers, subj_av45):
    subj_av45 = sorted(subj_av45, key=lambda x: x['EXAMDATE'])
    exam_times = [_['EXAMDATE'] for _ in subj_av45]
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
    all_tp_keys = ['AV45_TP_SPECIFIC_BL', 'AV45_TP_SPECIFIC_SCAN2', 'AV45_TP_SPECIFIC_SCAN3']
    all_wm70_composite_keys = wm70_composite_keys + \
                         ["%s_pchange" % _ for _ in wm70_composite_keys[1:]] + \
                         ["%s_pchange_ABS" % _ for _ in wm70_composite_keys[1:]] + \
                         ["%s_diff" % _ for _ in wm70_composite_keys[1:]] + \
                         ["%s_diff_ABS" % _ for _ in wm70_composite_keys[1:]] + \
                         ['AV45_WM70/composite_Slope_2pts', 'AV45_WM70/composite_Slope_3pts']
    all_wm70_cerebg_keys = wm70_cerebg_keys + \
                           ["%s_pchange" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ["%s_pchange_ABS" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ["%s_diff" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ["%s_diff_ABS" % _ for _ in wm70_cerebg_keys[1:]] + \
                           ['AV45_WM70/cerebg_Slope_2pts', 'AV45_WM70/cerebg_Slope_3pts']
    all_wm70_wcereb_keys = wm70_wcereb_keys + \
                           ["%s_pchange" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ["%s_pchange_ABS" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ["%s_diff" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ["%s_diff_ABS" % _ for _ in wm70_wcereb_keys[1:]] + \
                           ['AV45_WM70/wcereb_Slope_2pts', 'AV45_WM70/wcereb_Slope_3pts']
    all_unilateral_keys = unilateral_keys
    all_bigref_keys = bigref_keys + \
                      ["%s_pchange" % _ for _ in bigref_keys[1:]] + \
                      ["%s_pchange_ABS" % _ for _ in bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in bigref_keys[1:]] + \
                      ["%s_diff_ABS" % _ for _ in bigref_keys[1:]] + \
                      ["%s_BIN.79" % _ for _ in bigref_keys] + \
                      ['AV45_BigRef_Slope_2pts', 'AV45_BigRef_Slope_3pts']
    all_wm70_keys = wm70_keys + \
                    ["%s_pchange" % _ for _ in wm70_keys[1:]] + \
                    ["%s_pchange_ABS" % _ for _ in wm70_keys[1:]] + \
                    ["%s_diff" % _ for _ in wm70_keys[1:]] + \
                    ["%s_diff_ABS" % _ for _ in wm70_keys[1:]] + \
                    ["%s_BIN.62" % _ for _ in wm70_keys] + \
                    ['AV45_WM70_Slope_2pts', 'AV45_WM70_Slope_3pts']
    all_cerebg_keys = cerebg_keys + \
                    ["%s_pchange" % _ for _ in cerebg_keys[1:]] + \
                    ["%s_pchange_ABS" % _ for _ in cerebg_keys[1:]] + \
                    ["%s_diff" % _ for _ in cerebg_keys[1:]] + \
                    ["%s_diff_ABS" % _ for _ in cerebg_keys[1:]] + \
                    ["%s_BIN1.26" % _ for _ in cerebg_keys] + \
                    ['AV45_cerebg_Slope_2pts', 'AV45_cerebg_Slope_3pts']
    all_wcereb_keys = wcereb_keys + \
                    ["%s_pchange" % _ for _ in wcereb_keys[1:]] + \
                    ["%s_pchange_ABS" % _ for _ in wcereb_keys[1:]] + \
                    ["%s_diff" % _ for _ in wcereb_keys[1:]] + \
                    ["%s_diff_ABS" % _ for _ in wcereb_keys[1:]] + \
                    ["%s_BIN1.11" % _ for _ in wcereb_keys] + \
                    ['AV45_wcereb_Slope_2pts', 'AV45_wcereb_Slope_3pts']
    all_brainstem_keys = brainstem_keys + \
                    ["%s_pchange" % _ for _ in brainstem_keys[1:]] + \
                    ["%s_pchange_ABS" % _ for _ in brainstem_keys[1:]] + \
                    ["%s_diff" % _ for _ in brainstem_keys[1:]] + \
                    ["%s_diff_ABS" % _ for _ in brainstem_keys[1:]] + \
                    ["%s_BIN.79" % _ for _ in brainstem_keys] + \
                    ['AV45_brainstem_Slope_2pts', 'AV45_brainstem_Slope_3pts']
    all_wmratio_keys = wmratio_keys + \
                      ["%s_pchange" % _ for _ in wmratio_keys[1:]] + \
                      ["%s_pchange_ABS" % _ for _ in wmratio_keys[1:]] + \
                      ["%s_diff" % _ for _ in wmratio_keys[1:]] + \
                      ["%s_diff_ABS" % _ for _ in wmratio_keys[1:]]
    all_frontal_bigref_keys = frontal_bigref_keys + \
                      ["%s_pchange" % _ for _ in frontal_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in frontal_bigref_keys[1:]]
    all_cingulate_bigref_keys = cingulate_bigref_keys + \
                      ["%s_pchange" % _ for _ in cingulate_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in cingulate_bigref_keys[1:]]
    all_parietal_bigref_keys = parietal_bigref_keys + \
                      ["%s_pchange" % _ for _ in parietal_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in parietal_bigref_keys[1:]]
    all_temporal_bigref_keys = temporal_bigref_keys + \
                      ["%s_pchange" % _ for _ in temporal_bigref_keys[1:]] + \
                      ["%s_diff" % _ for _ in temporal_bigref_keys[1:]]
    all_av45_key_lists = [all_tp_keys,
                          all_wm70_composite_keys,
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
    new_headers = rearrangeHeaders(old_headers, all_av45_keys, after='AV45_1_3_Diff (yrs)')

    data = {k:'' for k in all_av45_keys}
    wm70_0 = None
    # fill in values
    for i, point in enumerate(subj_av45):
        # extract necessary values
        cerebw = float(point['CEREBELLUMWHITEMATTER'])
        wcereb = float(point['WHOLECEREBELLUM'])
        cerebg = float(point['CEREBELLUMGREYMATTER'])
        compositeroi = float(point['COMPOSITE'])
        bigref = float(point['COMPOSITE_REF'])
        wm70 = float(point['ERODED_SUBCORTICALWM'])
        brainstem = float(point['BRAINSTEM'])
        cingulate = float(point['CINGULATE'])
        frontal = float(point['FRONTAL'])
        parietal = float(point['PARIETAL'])
        temporal = float(point['TEMPORAL'])
        leftputamen = float(point['LEFT-PUTAMEN'])
        rightputamen = float(point['RIGHT-PUTAMEN'])
        leftcaudate = float(point['LEFT-CAUDATE'])
        rightcaudate = float(point['RIGHT-CAUDATE'])
        leftpallidum = float(point['LEFT-PALLIDUM'])
        rightpallidum = float(point['RIGHT-PALLIDUM'])
        tp_spec = 1 if point['TP_SPECIFIC'] else 0

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
            data['AV45_TP_SPECIFIC_BL'] = tp_spec
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
                data['AV45_TP_SPECIFIC_SCAN2'] = tp_spec
                times = exam_timedeltas[:2]
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in bigref_keys[:2]])
                data['AV45_BigRef_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wm70_keys[:2]])
                data['AV45_WM70_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in cerebg_keys[:2]])
                data['AV45_cerebg_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wcereb_keys[:2]])
                data['AV45_wcereb_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in brainstem_keys[:2]])
                data['AV45_brainstem_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wm70_composite_keys[:2]])
                data['AV45_WM70/composite_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wm70_cerebg_keys[:2]])
                data['AV45_WM70/cerebg_Slope_2pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wm70_wcereb_keys[:2]])
                data['AV45_WM70/wcereb_Slope_2pts'] = slope
            if i == 2:
                data['AV45_TP_SPECIFIC_SCAN3'] = tp_spec
                times = exam_timedeltas[:3]
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in bigref_keys[:3]])
                data['AV45_BigRef_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wm70_keys[:3]])
                data['AV45_WM70_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in cerebg_keys[:3]])
                data['AV45_cerebg_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wcereb_keys[:3]])
                data['AV45_wcereb_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in brainstem_keys[:3]])
                data['AV45_brainstem_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wm70_composite_keys[:3]])
                data['AV45_WM70/composite_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wm70_cerebg_keys[:3]])
                data['AV45_WM70/cerebg_Slope_3pts'] = slope
                slope, intercept, r, p, stderr = stats.linregress(times, [data[_] for _ in wm70_wcereb_keys[:3]])
                data['AV45_WM70/wcereb_Slope_3pts'] = slope
    return (new_headers, data)

def syncFDGData(old_headers, old_lines, fdg_file, registry_file, dump_to=None):
    fdg_by_subj = importFDG(fdg_file)

    fdg_columns = []
    fdg_columns += ['FDG_pons_Vis%s' % (i+1) for i in range(10)]
    fdg_columns += ['FDGVis%s_ReltoAV45' % (i+1) for i in range(10)]
    fdg_columns += ['FDGVis%s_Date' % (i+1) for i in range(10)]
    fdg_columns += ['FDG_postAV45_slope','FDG_postAV45_followuptime']
    new_headers= rearrangeHeaders(old_headers, fdg_columns, after='FDG_Bin_Baseline')

    def extraction_fn(subj, subj_row, old_l, patient_pets):
        # collapse subj rows per visit
        av45_bl, av45_2, av45_3 = getAV45Dates(old_l)
        examdates = list(set([_['EXAMDATE'] for _ in subj_row]))
        examdates = sorted(examdates)
        subj_fdg = []
        for ed in examdates:
            rel_rows = [_ for _ in subj_row if _['EXAMDATE'] == ed]
            means = [_['MEAN'] for _ in rel_rows]
            pons_avg = np.mean(means)
            subj_fdg.append({'EXAMDATE': ed, 'PONS': pons_avg})

        new_data = {}
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
            slope, intercept, r, p, stderr = stats.linregress(raw_dates, raw_scores)
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

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def syncTBMSynData(old_headers, old_lines, tbm_file, registry_file, dump_to=None):
    tbm_by_subj = importTBMSyn(tbm_file)

    # add new headers as needed
    tbm_columns = ['TBMSyn_DATE.%s' % (i+1) for i in range(10)]
    tbm_columns.extend(['TBMSyn_SCORE.%s' % (i+1) for i in range(10)])
    tbm_columns.extend(['TBMSyn_postAV45.%s' % (i+1) for i in range(10)])
    tbm_columns.append('TBMSyn_BL_DATE')
    tbm_columns.append('TBMSyn_count')
    #tbm_columns.append('TBMSyn_SLOPE')
    #tbm_columns.append('TBMSyn_slope_followuptime')
    new_headers = rearrangeHeaders(old_headers, tbm_columns, after=None)
    
    def extraction_fn(subj, subj_row, old_l, patient_pets):
        subj_tbm = sorted(subj_row, key=lambda x: x['EXAMDATE'])
        av45_bl, av45_2, av45_3 = getAV45Dates(old_l)
        bl_examdate = subj_tbm[0]['BL_EXAMDATE']
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
                new_data['TBMSyn_SCORE.%s' % (i+1)] = datapoint['SCORE']
                timediff = (datapoint['EXAMDATE']-datapoint['BL_EXAMDATE']).days / 365.0
                slope_points.append((timediff,datapoint['SCORE']))
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
            #slope, intercept, r, p, stderr = stats.linregress(diff_dates, diff_scores)
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

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def syncCSFData(old_headers, old_lines, csf_files, registry_file, dump_to=None):
    registry = importRegistry(registry_file)
    csf_by_subj = importCSF(csf_files, registry)
    
    # add new headers as needed
    csf_headers = []
    csf_headers += ['CSF_ABETA.%s' % (i+1) for i in range(7)]
    csf_headers += ['CSF_ABETApostAV45.%s' % (i+1) for i in range(7)]
    csf_headers += ['CSF_ABETA_slope', 'CSF_ABETA_closest_AV45', 
                    'CSF_ABETA_closest_AV45_2', 'CSF_ABETA_closest_AV45_3',
                    'CSF_ABETA_closest_AV45_BIN_192',
                    'CSF_ABETA_closest_AV45_2_BIN_192',
                    'CSF_ABETA_closest_AV45_3_BIN_192']
    csf_headers += ['CSF_TAU.%s' % (i+1) for i in range(7)]
    csf_headers += ['CSF_TAUpostAV45.%s' % (i+1) for i in range(7)]
    csf_headers += ['CSF_TAU_slope', 'CSF_TAU_closest_AV45', 
                    'CSF_TAU_closest_AV45_2', 'CSF_TAU_closest_AV45_3',
                    'CSF_TAU_closest_AV45_BIN_93',
                    'CSF_TAU_closest_AV45_2_BIN_93',
                    'CSF_TAU_closest_AV45_3_BIN_93']
    csf_headers += ['CSF_PTAU.%s' % (i+1) for i in range(7)]
    csf_headers += ['CSF_PTAUpostAV45.%s' % (i+1) for i in range(7)]
    csf_headers += ['CSF_PTAU_slope', 'CSF_PTAU_closest_AV45', 
                    'CSF_PTAU_closest_AV45_2', 'CSF_PTAU_closest_AV45_3',
                    'CSF_PTAU_closest_AV45_BIN_23',
                    'CSF_PTAU_closest_AV45_2_BIN_23',
                    'CSF_PTAU_closest_AV45_3_BIN_23']
    new_headers = rearrangeHeaders(old_headers, csf_headers, after=None)

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
        # find av45 dates
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)

        new_data = {}
        ptau_slope_points = []
        abeta_slope_points = []
        tau_slope_points = []
        for i in range(7):
            if i < len(subj_csf):
                datapoint = subj_csf[i]
                examdate = datapoint['EXAMDATE']
                timediff = ((examdate-bl_av45).days / 365.0) if bl_av45 else ''
                if timediff <= -(90.0/365.0):
                    timediff = ''
                abeta_val = datapoint['abeta']
                ptau_val = datapoint['ptau']
                tau_val = datapoint['tau']
                abeta_val = float(abeta_val) if abeta_val is not None else ''
                ptau_val = float(ptau_val) if ptau_val is not None else ''
                tau_val = float(tau_val) if tau_val is not None else ''
                new_data['CSF_ABETA.%s' % (i+1)] = abeta_val
                new_data['CSF_PTAU.%s' % (i+1)] = ptau_val
                new_data['CSF_TAU.%s' % (i+1)] = tau_val
                new_data['CSF_ABETApostAV45.%s' % (i+1)] = timediff
                new_data['CSF_PTAUpostAV45.%s' % (i+1)] = timediff
                new_data['CSF_TAUpostAV45.%s' % (i+1)] = timediff
                if timediff:
                    if abeta_val != '':
                        abeta_slope_points.append((timediff, abeta_val))
                    if tau_val != '':
                        tau_slope_points.append((timediff, tau_val))
                    if ptau_val != '':
                        ptau_slope_points.append((timediff, ptau_val))
            else:
                new_data['CSF_ABETA.%s' % (i+1)] = ''
                new_data['CSF_PTAU.%s' % (i+1)] = ''
                new_data['CSF_TAU.%s' % (i+1)] = ''
                new_data['CSF_ABETApostAV45.%s' % (i+1)] = ''
                new_data['CSF_PTAUpostAV45.%s' % (i+1)] = ''
                new_data['CSF_TAUpostAV45.%s' % (i+1)] = ''

        # match up with av45 scans
        av45_bl_closest, av45_2_closest, av45_3_closest = getClosestToAV45(subj_csf, bl_av45, av45_2, av45_3, day_limit=93)
        if av45_bl_closest:
            tau_val = av45_bl_closest['tau']
            ptau_val = av45_bl_closest['ptau']
            abeta_val = av45_bl_closest['abeta']
            abeta_val = float(abeta_val) if abeta_val is not None else ''
            ptau_val = float(ptau_val) if ptau_val is not None else ''
            tau_val = float(tau_val) if tau_val is not None else ''
            new_data['CSF_TAU_closest_AV45'] = tau_val
            new_data['CSF_PTAU_closest_AV45'] = ptau_val
            new_data['CSF_ABETA_closest_AV45'] = abeta_val
            new_data['CSF_TAU_closest_AV45_BIN_93'] = 1 if tau_val >= 93 else 0
            new_data['CSF_PTAU_closest_AV45_BIN_23'] = 1 if ptau_val >= 23 else 0
            new_data['CSF_ABETA_closest_AV45_BIN_192'] = 1 if abeta_val <= 192 else 0
        if av45_2_closest:
            tau_val = av45_2_closest['tau']
            ptau_val = av45_2_closest['ptau']
            abeta_val = av45_2_closest['abeta']
            abeta_val = float(abeta_val) if abeta_val is not None else ''
            ptau_val = float(ptau_val) if ptau_val is not None else ''
            tau_val = float(tau_val) if tau_val is not None else ''
            new_data['CSF_TAU_closest_AV45_2'] = tau_val
            new_data['CSF_PTAU_closest_AV45_2'] = ptau_val
            new_data['CSF_ABETA_closest_AV45_2'] = abeta_val
            new_data['CSF_TAU_closest_AV45_2_BIN_93'] = 1 if tau_val >= 93 else 0
            new_data['CSF_PTAU_closest_AV45_2_BIN_23'] = 1 if ptau_val >= 23 else 0
            new_data['CSF_ABETA_closest_AV45_2_BIN_192'] = 1 if abeta_val <= 192 else 0
        if av45_3_closest:
            tau_val = av45_3_closest['tau']
            ptau_val = av45_3_closest['ptau']
            abeta_val = av45_3_closest['abeta']
            abeta_val = float(abeta_val) if abeta_val is not None else ''
            ptau_val = float(ptau_val) if ptau_val is not None else ''
            tau_val = float(tau_val) if tau_val is not None else ''
            new_data['CSF_TAU_closest_AV45_3'] = tau_val
            new_data['CSF_PTAU_closest_AV45_3'] = ptau_val
            new_data['CSF_ABETA_closest_AV45_3'] = abeta_val
            new_data['CSF_TAU_closest_AV45_3_BIN_93'] = 1 if tau_val >= 93 else 0
            new_data['CSF_PTAU_closest_AV45_3_BIN_23'] = 1 if ptau_val >= 23 else 0
            new_data['CSF_ABETA_closest_AV45_3_BIN_192'] = 1 if abeta_val <= 192 else 0

        # calculate slope
        if len(ptau_slope_points) >= 2:
            raw_dates = [_[0] for _ in ptau_slope_points]
            raw_scores = [_[1] for _ in ptau_slope_points]
            slope, intercept, r, p, stderr = stats.linregress(raw_dates, raw_scores)
            new_data['CSF_PTAU_slope'] = slope
        if len(tau_slope_points) >= 2:
            raw_dates = [_[0] for _ in tau_slope_points]
            raw_scores = [_[1] for _ in tau_slope_points]
            slope, intercept, r, p, stderr = stats.linregress(raw_dates, raw_scores)
            new_data['CSF_TAU_slope'] = slope
        if len(abeta_slope_points) >= 2:
            raw_dates = [_[0] for _ in abeta_slope_points]
            raw_scores = [_[1] for _ in abeta_slope_points]
            slope, intercept, r, p, stderr = stats.linregress(raw_dates, raw_scores)
            new_data['CSF_ABETA_slope'] = slope
        return new_data


    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, csf_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)

def syncWMHData(old_headers, old_lines, wmh_file, registry_file, dump_to=None):
    wmh_by_subj = importWMH(wmh_file)
    registry = importRegistry(registry_file)

    to_add_headers = ['WMH_percentOfICV.%s' % (i+1) for i in range(7)]
    to_add_headers += ['WMH_postAV45.%s' % (i+1) for i in range(7)]
    to_add_headers += ['WMH_WHITMATHYP.%s' % (i+1) for i in range(7)]
    to_add_headers.append('WMH_slope')
    new_headers = rearrangeHeaders(old_headers, to_add_headers, after=None)
    
    def extraction_fn(subj, subj_row, old_l, patient_pets):
        bl_av45, av45_2, av45_3 = getAV45Dates(old_l, patient_pets=patient_pets)
        subj_wmh = sorted(subj_row, key=lambda x: x['EXAMDATE'])
        new_data = {}
        slope_points = []
        for i in range(7):
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
        if len(slope_points) >= 2:
            raw_dates = [_[0] for _ in slope_points]
            raw_scores = [_[1] for _ in slope_points]
            slope, intercept, r, p, stderr = stats.linregress(raw_dates, raw_scores)
            new_data['WMH_slope'] = slope
        return new_data

    new_lines = []
    for old_l in old_lines:
        new_data = updateLine(old_l, wmh_by_subj, extraction_fn, 
                              pid_key='RID', pet_meta=None, decimal_places=5)
        old_l.update(new_data)
        new_lines.append(old_l)

    # dump out
    if dump_to is not None:
        dumpCSV(dump_to, new_headers, new_lines)

    return (new_headers, new_lines)


def getClosestToAV45(points, bl_av45, av45_2, av45_3, day_limit=550):
    used = set()
    # match up with av45 scans
    av45_bl_closest = av45_2_closest = av45_3_closest = None
    if bl_av45 is not None:
        eligible = [p for p in points if p['EXAMDATE'] not in used and (abs(bl_av45 - p['EXAMDATE']).days <= day_limit)]
        cand = sorted(eligible, key=lambda x: abs(bl_av45-x['EXAMDATE']))
        if len(cand) > 0:
            av45_bl_closest = cand[0]
            used.add(av45_bl_closest['EXAMDATE'])
    if av45_2 is not None:
        eligible = [p for p in points if p['EXAMDATE'] not in used and (abs(av45_2 - p['EXAMDATE']).days <= day_limit)]
        cand = sorted(eligible, key=lambda x: abs(av45_2-x['EXAMDATE']))
        if len(cand) > 0:
            av45_2_closest = cand[0]
            used.add(av45_2_closest['EXAMDATE'])
    if av45_3 is not None:
        eligible = [p for p in points if p['EXAMDATE'] not in used and (abs(av45_3 - p['EXAMDATE']).days <= day_limit)]
        cand = sorted(eligible, key=lambda x: abs(av45_3-x['EXAMDATE']))
        if len(cand) > 0:
            av45_3_closest = cand[0]
            used.add(av45_3_closest['EXAMDATE'])
    return av45_bl_closest, av45_2_closest, av45_3_closest


def eliminateColumns(headers, lines):
    to_remove = ['LastClinicalVisit',
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
                 'TBMSyn_slope_followuptime']
    to_remove += ['WHITMATHYP.%s' % (i+1) for i in range(5)]
    to_remove += ['ABETA.%s' % (i+1) for i in range(7)]
    to_remove += ['TAU.%s' % (i+1) for i in range(7)]
    to_remove += ['PTAU.%s' % (i+1) for i in range(7)]
    to_remove += ['EXAM_tau.%s' % (i+1) for i in range(7)]
    to_remove += ['ADAScog_DATE%s' % (i+1) for i in range(11)]
    to_remove += ['MMSE_DATE%s' % (i+1) for i in range(11)]
    to_remove += ['AVLT_DATE.%s' % (i+1) for i in range(11)]
    to_remove += ['TBMSyn_DATE.%s' % (i+1) for i in range(11)]


    for tm in to_remove:
        while tm in headers:
            headers.remove(tm)
    new_lines = []
    for l in lines:
        for tm in to_remove:
            l.pop(tm, None)
        new_lines.append(l)
    return headers, new_lines

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


def addCategories(output_file):
    df = pd.read_csv(output_file, low_memory=False)
    col_dict = dict([(b.strip(),a) for a,b_list in COLUMN_CATEGORIES.iteritems() for b in b_list])
    cats = [col_dict.get(str(_).strip(),'Unknown') for _ in df.columns]
    for a,b in zip(cats, df.columns):
        if a == 'Unknown':
            print b
    df.columns = [cats, df.columns]
    df.to_csv(output_file, index=False)

def runPipeline():
    # syncing pipeline
    new_headers, new_lines = parseCSV(master_file)
    print "\nELIMINATING COLUMNS\n"
    new_headers, new_lines = eliminateColumns(new_headers, new_lines)
    print "\nSYNCING AV45\n"
    new_headers, new_lines = syncAV45Data(new_headers, new_lines, av45_file, av45_nontp_file, registry_file, dump_to=None) # adds new patients
    print "\nSYNCING DIAGNOSES\n"
    new_headers, new_lines = syncDiagnosisData(new_headers, new_lines, diagnosis_file, registry_file, demog_file, arm_file, pet_meta_file, dump_to=None) # refreshes av45 dates
    print "\nSYNCING ROUSSET BL\n"
    new_headers, new_lines = syncRoussetResults(new_headers, new_lines, rousset_matfile_bl_manual, 'BL', dump_to=None)
    new_headers, new_lines = syncRoussetResults(new_headers, new_lines, rousset_matfile_bl_agg, 'BL', dump_to=None)
    print "\nSYNCING ROUSSET SCAN2\n"
    new_headers, new_lines = syncRoussetResults(new_headers, new_lines, rousset_matfile_scan2_manual, 'Scan2', dump_to=None)
    new_headers, new_lines = syncRoussetResults(new_headers, new_lines, rousset_matfile_scan2_agg, 'Scan2', dump_to=None)
    print "\nSYNCING ROUSSET SCAN3\n"
    new_headers, new_lines = syncRoussetResults(new_headers, new_lines, rousset_matfile_scan3_manual, 'Scan3', dump_to=None)
    new_headers, new_lines = syncRoussetResults(new_headers, new_lines, rousset_matfile_scan3_agg, 'Scan3', dump_to=None)
    print "\nSYNCING FDG\n"
    new_headers, new_lines = syncFDGData(new_headers, new_lines, fdg_file, registry_file, dump_to=None)
    print "\nSYNCING TBMSYN\n"
    new_headers, new_lines = syncTBMSynData(new_headers, new_lines, tbm_file, registry_file, dump_to=None)
    print "\nSYNCING MMSE\n"
    new_headers, new_lines = syncMMSEData(new_headers, new_lines, mmse_file, registry_file, dump_to=None)
    print "\nSYNCING ADASCOG\n"
    new_headers, new_lines = syncADASCogData(new_headers, new_lines, adni1_adas_file, adnigo2_adas_file, registry_file, dump_to=None)
    print "\nSYNCING AVLT\n"
    new_headers, new_lines = syncAVLTData(new_headers, new_lines, neuro_battery_file, registry_file, dump_to=None)
    print "\nSYNCING CSF\n"
    new_headers, new_lines = syncCSFData(new_headers, new_lines, csf_files, registry_file, dump_to=None)
    print "\nSYNCING GD\n"
    new_headers, new_lines = syncGDData(new_headers, new_lines, gd_file, registry_file, dump_to=None)
    print "\nSYNCING WMH\n"
    new_headers, new_lines = syncWMHData(new_headers, new_lines, wmh_file, registry_file, dump_to=output_file)

    # add categories
    addCategories(output_file)

if __name__ == '__main__':
    now = datetime.now()
    # IO files
    master_file = "../FDG_AV45_COGdata.csv"
    output_file = "../FDG_AV45_COGdata_%s.csv" % (now.strftime("%m_%d_%y"))
    # LOOKUP files
    registry_file = "../docs/registry_clean.csv"
    arm_file = "../docs/ARM.csv"
    pet_meta_file = "../docs/PET_META_LIST_edited.csv"
    demog_file = "../docs/PTDEMOG.csv"
    # AV45 File
    av45_file = "../output/UCBERKELEYAV45_09_25_15_extra.csv"
    av45_nontp_file = "../output/UCBERKELEYAV45_09_25_15_extra_nontp.csv"
    # MMSE files
    mmse_file = "../cog_tests/MMSE.csv"
    # ADAS-COG files
    adni1_adas_file = '../cog_tests/ADASSCORES.csv'
    adnigo2_adas_file = '../cog_tests/ADAS_ADNIGO2.csv'
    # Neuropsychological Battery file
    neuro_battery_file = '../cog_tests/NEUROBAT.csv'
    # Diagnosis file
    diagnosis_file = '../docs/DXSUM_PDXCONV_ADNIALL.csv'
    # TBMsyn file
    tbm_file = '../mr_docs/Mayo/MAYOADIRL_MRI_TBMSYN_05_07_15.csv'
    # CSF files
    csf_files = ['../docs/UPENNBIOMK5_firstset.csv','../docs/UPENNBIOMK6_secondset.csv','../docs/UPENNBIOMK7_thirdset.csv','../docs/UPENNBIOMK8_fourthset.csv']
    # WMH file
    wmh_file = '../docs/UCD_ADNI2_WMH_03_12_15.csv'
    # GD file
    gd_file = '../docs/GDSCALE.csv'
    # FDG file
    fdg_file = '../docs/UCBERKELEY_FDG_07_29_15.csv'
    # Rousset output files
    rousset_matfile_bl_manual = '../output/Rousset_BL/rousset_outputs_manual.mat'
    rousset_matfile_scan2_manual = '../output/Rousset_Scan2/rousset_outputs_manual.mat'
    rousset_matfile_scan3_manual = '../output/Rousset_Scan3/rousset_outputs_manual.mat'
    rousset_matfile_bl_agg = '../output/Rousset_BL/rousset_outputs_agg.mat'
    rousset_matfile_scan2_agg = '../output/Rousset_Scan2/rousset_outputs_agg.mat'
    rousset_matfile_scan3_agg = '../output/Rousset_Scan3/rousset_outputs_agg.mat'

    runPipeline()
