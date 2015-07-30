from collections import defaultdict
from datetime import datetime, timedelta
from tabulate import tabulate

from utils import *


def checkAvailablePointsPerSubject(pet_data, bsi_data, longfree_data, longfree_data_adni1, crossfree_data, tbm_data, mri_data, master_data, numerical_output):
    '''
    Aggregate into one dictionary
    remove datapoints from before the first av45 scan
    '''

    '''
    interesting = [21, 23, 31, 56, 58, 59, 61, 69, 74, 89, 96, 106, 113, 120, 
                   123, 130, 159, 166, 171, 173, 230, 257, 259, 260, 272, 301, 
                   311, 315, 413, 416, 419, 498, 545, 602, 610, 618, 680, 684, 
                   685, 717, 741, 751, 767, 778, 842, 896, 907, 920, 923, 926, 
                   969, 972, 981, 985, 1016, 1098, 1123, 1169, 1190, 1206, 1232, 
                   1261, 1280, 4003, 4014, 4018, 4020, 4021, 4026, 4028, 4032, 
                   4037, 4041, 4043, 4050, 4060, 4075, 4076, 4081, 4082, 4084, 
                   4086, 4090, 4100, 4103, 4105, 4119, 4120, 4121, 4148, 4150, 
                   4151, 4158, 4164, 4173, 4176, 4177, 4179, 4196, 4198, 4200, 
                   4208, 4213, 4218, 4222, 4224, 4225, 4234, 4254, 4255, 4262, 
                   4266, 4269, 4270, 4275, 4278, 4290, 4291, 4292, 4313, 4320, 
                   4339, 4343, 4345, 4349, 4350, 4352, 4357, 4365, 4367, 4369, 
                   4371, 4372, 4376, 4382, 4384, 4385, 4386, 4387, 4388, 4389, 
                   4391, 4393, 4396, 4399, 4400, 4410, 4422, 4427, 4428, 4429, 
                   4441, 4446, 4448, 4453, 4466, 4469, 4482, 4483, 4485, 4491, 
                   4496, 4499, 4503, 4505, 4506, 4516, 4520, 4545, 4552, 4555, 
                   4559, 4560, 4566, 4576, 4579, 4580, 4585, 4586, 4587, 4598, 
                   4599, 4604, 4607, 4612, 4616, 4620, 4632, 4637, 4643, 4645, 
                   4649, 4739, 4762, 4795, 4832, 4835, 4843, 4855, 4872, 4878, 
                   5023]
    print len(interesting)
    interesting_counts = defaultdict(list)
    '''
    headers = ['Source', 'Diagnosis', 'Phase', 'BL', 'year1', 'year2', 'year3', 'year4', 'year5', 'm3', 'm6', 'm12', 'm24', 'm36', 'm48']
    rows = []

    lookupfile = '../docs/ADNI2_VISITID.csv'
    vc_lookup = createVisitIDLookup(lookupfile)

    adni1_mri_visit_counts = {'N': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": []},
                              'SMC': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": []},
                              'EMCI': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": []},
                              'LMCI': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": []},
                              'AD': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": []}}
    adni2_mri_visit_counts = {'N': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'SMC': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'EMCI': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'LMCI': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'AD': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []}}
    adni1_tbm_visit_counts = {'N': {"year1": [], "year2": [], "year3": [], "year4": [], "year5": []},
                              'SMC': {"year1": [], "year2": [], "year3": [], "year4": [], "year5": []},
                              'EMCI': {"year1": [], "year2": [], "year3": [], "year4": [], "year5": []},
                              'LMCI': {"year1": [], "year2": [], "year3": [], "year4": [], "year5": []},
                              'AD': {"year1": [], "year2": [], "year3": [], "year4": [], "year5": []}}
    adni2_tbm_visit_counts = {'N': {"m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'SMC': {"m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'EMCI': {"m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'LMCI': {"m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'AD': {"m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []}}
    adni1_long_visit_counts = {'N': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": []},
                              'SMC': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": []},
                              'EMCI': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": []},
                              'LMCI': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": []},
                              'AD': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": []}}
    adni2_long_visit_counts = {'N': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'SMC': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'EMCI': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'LMCI': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'AD': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []}}
    adni1_cross_visit_counts = {'N': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": [], "year5": []},
                              'SMC': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": [], "year5": []},
                              'EMCI': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": [], "year5": []},
                              'LMCI': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": [], "year5": []},
                              'AD': {"BL": [], "year1": [], "year2": [], "year3": [], "year4": [], "year5": []}}
    adni2_cross_visit_counts = {'N': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'SMC': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'EMCI': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'LMCI': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'AD': {"BL": [], "m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []}}
    adni1_bsi_visit_counts = {'N': {"year1": [], "year2": [], "year3": [], "year4": []},
                              'SMC': {"year1": [], "year2": [], "year3": [], "year4": []},
                              'EMCI': {"year1": [], "year2": [], "year3": [], "year4": []},
                              'LMCI': {"year1": [], "year2": [], "year3": [], "year4": []},
                              'AD': {"year1": [], "year2": [], "year3": [], "year4": []}}
    adni2_bsi_visit_counts = {'N': {"m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'SMC': {"m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'EMCI': {"m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'LMCI': {"m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []},
                              'AD': {"m3": [], "m6": [], "m12": [], "m24": [], "m36": [], "m48": []}}

    alls = {'bsi': set(),
            'long': set(),
            'cross': set(),
            'tbm': set(),
            'mri': set(),
    }

    points_by_subj = {}
    for subj,v in pet_data.iteritems():
        print "SUBJECT: %s" % subj
        patient_pets = sorted(v)
        bl_av45 = patient_pets[0] - timedelta(days=150)
        subj_points = {}

        if subj in bsi_data:
            subj_bsi = sorted(bsi_data[subj], key=lambda x: x['EXAMDATE'])
        else:
            subj_bsi = []
        subj_points['bsi'] = [_['VISCODE2'].replace('m0','m') for _ in subj_bsi if _['BL_EXAMDATE'] >= bl_av45]

        if subj in longfree_data or subj in longfree_data_adni1:
            all_longfree_data = longfree_data.get(subj,[]) + longfree_data_adni1.get(subj,[])
            subj_long = sorted(all_longfree_data, key=lambda x: x['EXAMDATE'])
        else:
            subj_long = []
        subj_points['long'] = [_['VISCODE2'].replace('m0','m').replace('scmri','bl') for _ in subj_long if _['EXAMDATE'] >= bl_av45]

        subj_cross = sorted(crossfree_data.get(subj,[]), key=lambda x: x['EXAMDATE'])
        subj_points['cross'] = [_['VISCODE2'].replace('m0','m').replace('scmri','bl') for _ in subj_cross if _['EXAMDATE'] >= bl_av45]

        if subj in tbm_data:
            subj_tbm = sorted(tbm_data[subj], key=lambda x: x['EXAMDATE'])
            subj_tbm_bl = tbm_data[subj][0]['VISCODEBL']
        else:
            subj_tbm = []
            subj_tbm_bl = None
        subj_points['tbm'] = [_['VISCODE2'].replace('m0','m') for _ in subj_tbm if _['BL_EXAMDATE'] >= bl_av45]

        if subj in mri_data:
            subj_mri = sorted(list(set([(_['EXAMDATE'],_['vc']) for _ in mri_data[subj] if _['EXAMDATE'] >= bl_av45])))
        else:
            subj_mri = []
        subj_points['mri'] = [_ for _ in subj_mri if _[0] >= bl_av45]
        if len(subj_points['mri']) > 0:
            first_mri_time = subj_points['mri'][0][0]
            mri_diffs = [(a-first_mri_time).days/ 365.0 for a,b in subj_points['mri']]
            mri_vc = [convertVisitName(b) if convertVisitName(b) else b for a,b in subj_points['mri']]
            mri_vc = [vc_lookup.get(subj,{}).get(_,_).replace('m0','m').replace('scmri','bl') for _ in mri_vc]
            if subj >= 2000:
                print subj_mri
                print mri_vc
        else:
            mri_diffs = []
            mri_vc = []
        #print "%s: %s -> %s" % (subj, len(subj_mri), len(subj_points['mri']))

        # remove duplicates from subject points
        for k in subj_points.keys():
            values = subj_points[k]
            seen = set()
            filtered = []
            for v in subj_points[k]:
                if v not in seen:
                    filtered.append(v)
                    seen.add(v)
            subj_points[k] = filtered

        points_by_subj[subj] = subj_points
        diag = master_data.get(subj,{}).get('Init_Diagnosis','Unknown')

        print subj
        print first_mri_time
        print subj_tbm_bl
        for k,v in subj_points.items():
            if subj >= 2000:
                if k != 'mri':
                    alls[k] |= set(v)
                    print "%s: %s" % (k,v)
                if k == 'mri':
                    alls[k] |= set(mri_vc)
            
        # add to adni1 or adnigo/2 counts
        chosen_timepoints = []
        if subj < 2000:
            print subj
            years = [0,1,2,3,4]

            # Match MRI by time
            used = set([])
            for y in years:
                year_found = None
                for mridiff in sorted(mri_diffs, key=lambda x: abs(y-x)):
                    if abs(y - mridiff) < 0.6 and mridiff not in used:
                        year_found = y
                        chosen_timepoints.append((y,mridiff))
                        break
                if year_found is not None and year_found not in used:
                    used.add(mridiff)
                    if year_found == 0 and "BL" not in used:
                        adni1_mri_visit_counts[diag]["BL"].append(subj)
                        used.add("BL")
                    elif "year%s" % year_found not in used:
                        adni1_mri_visit_counts[diag]["year%s" % year_found].append(subj)
                        used.add("year%s" % year_found)

            if len(used)/2 != len(mri_diffs):
                print "FAILED MRI MATCH"
            
            # Match TBM by vc
            used = set()
            for vc in subj_points['tbm']:
                if not vc.startswith('m'):
                    continue
                # get month diff 
                month_diff = int(vc.replace('m','')) - int(subj_tbm_bl.replace('m',''))
                year_key = "year%s" % int(round(month_diff / 12))
                if year_key not in used:
                    adni1_tbm_visit_counts[diag][year_key].append(subj)
                    used.add(year_key)

            # Match BSI by vc
            if len(subj_points['bsi']) > 0:
                raise Exception("Getting BSI data, write code for it")
                used = set()
                for vc in subj_points['bsi']:
                    if not vc.startswith('m'):
                        continue
                    if vc not in used:
                        adni1_bsi_visit_counts[diag][vc].append(subj)
                        used.add(vc)

            # Match Crossfree by vc
            if len(subj_points['cross']) > 0:
                print subj_points['cross']
                cross_bl = subj_points['cross'][0].lower().strip()
                if cross_bl == 'bl':
                    cross_bl = 0
                    print "CHECK: %s" % subj
                else:
                    cross_bl = int(cross_bl.replace('m',''))
                used = set()
                for vc in subj_points['cross']:
                    if not (vc.startswith('m') or vc == 'bl'):
                        continue
                    # get month diff 
                    if vc == 'bl':
                        month_diff = 0
                    else:
                        month_diff = int(vc.replace('m','')) - cross_bl
                    if month_diff == 0:
                        year_key = "BL"
                    else:
                        year_key = "year%s" % int(round(month_diff / 12))
                    if year_key not in used:
                        adni1_cross_visit_counts[diag][year_key].append(subj)
                        used.add(year_key)

            # Match Longfree by vc
            if len(subj_points['long']) > 0:
                print subj_points['long']
                long_bl = subj_points['long'][0].lower().strip()
                if long_bl == 'bl':
                    long_bl = 0
                    print "CHECK: %s" % subj
                else:
                    long_bl = int(long_bl.replace('m',''))
                used = set()
                for vc in subj_points['long']:
                    if not (vc.startswith('m') or vc == 'bl'):
                        continue
                    # get month diff 
                    if vc == 'bl':
                        month_diff = 0
                    else:
                        month_diff = int(vc.replace('m','')) - long_bl
                    if month_diff == 0:
                        year_key = "BL"
                    else:
                        year_key = "year%s" % int(round(month_diff / 12))
                    if year_key not in used:
                        adni1_long_visit_counts[diag][year_key].append(subj)
                        used.add(year_key)

        elif subj >= 2000:
            months = [(0.0,0),(0.25,3),(0.5,6),(1,12),(2,24),(3,36),(4,48)]

            # Match MRI by name, then by time
            '''
            used_months = set([])
            used_cols = set([])
            print mri_diffs
            print mri_vc
            for m_annual, m in months:
                month_found = None
                for mridiff, mrivc in sorted(zip(mri_diffs, mri_vc), key=lambda x: abs(m_annual-x[0])):
                    if mrivc not in used_cols:
                        adni2_mir_visit_counts[diag][vc].append(subj)
                        used_cols.append(vc)
                    if mridiff in used_months:
                        continue
                
                if vc.lower() == 'bl':
                    if 'bl' not in used:
                        adni2_cross_visit_counts[diag]['BL'].append(subj)
                        used.add('bl')
                    continue
                if not vc.startswith('m'):
                    continue
                if vc not in used:
                    adni2_cross_visit_counts[diag][vc].append(subj)
                    used.add(vc)


                    if mrivc.lower() == 'scmri':
                        if 'bl' not in used:
                            chosen_timepoints.append((m,mridiff))
                            adni2_long_visit_counts[diag]['BL'].append(subj)
                            used_months.add(mridiff)
                            break
                    if vc.startswith('m'):
                        adni2_long_visit_counts[diag][vc].append(subj)
                        used.add(vc)

                if month_found is not None:
                    if month_found == 0 and "BL" not in used_cols:
                        adni2_mri_visit_counts[diag]["BL"].append(subj)
                        used_cols.add("BL")
                    elif "m%s" % month_found not in used_cols:
                        adni2_mri_visit_counts[diag]["m%s" % month_found].append(subj)
                        used_cols.add("m%s" % month_found)
            for m_annual, m in months:
                month_found = None
                if m in used_cols:
                    continue
                for mridiff, mrivc in sorted(zip(mri_diffs, mri_vc), key=lambda x: abs(m_annual-x[0])):
                    if mridiff in used_months:
                        continue
                    if abs(m_annual - mridiff) < 0.8:
                        print "MIRACULOUS"
                        month_found = m
                        used_months.add(mridiff)
                        chosen_timepoints.append((m,mridiff))
                        break
                if month_found is not None and month_found not in used_cols:
                    if month_found == 0:
                        adni2_mri_visit_counts[diag]["BL"].append(subj)
                    else:
                        adni2_mri_visit_counts[diag]["m%s" % month_found].append(subj)
            if len(used_months) != len(mri_diffs):
                print 'FAILED MRI MATCH'
            '''

            # Match MRI by vc
            used = set()
            for vc in mri_vc:
                if vc.lower() == 'bl':
                    if 'bl' not in used:
                        adni2_mri_visit_counts[diag]['BL'].append(subj)
                        used.add('bl')
                    continue
                if not vc.startswith('m'):
                    continue
                if vc not in used:
                    adni2_mri_visit_counts[diag][vc].append(subj)
                    used.add(vc)

            # Match TBM by vc
            used = set()
            for vc in subj_points['tbm']:
                if vc not in used:
                    adni2_tbm_visit_counts[diag][vc].append(subj)
                    used.add(vc)

            # Match BSI by vc
            used = set()
            for vc in subj_points['bsi']:
                if vc not in used:
                    adni2_bsi_visit_counts[diag][vc].append(subj)
                    used.add(vc)

            # Match Crossfree by vc
            used = set()
            for vc in subj_points['cross']:
                if vc.lower() == 'bl':
                    if 'bl' not in used:
                        adni2_cross_visit_counts[diag]['BL'].append(subj)
                        used.add('bl')
                    continue
                if not vc.startswith('m'):
                    continue
                if vc not in used:
                    adni2_cross_visit_counts[diag][vc].append(subj)
                    used.add(vc)

            # Match Longfree by vc
            used = set()
            for vc in subj_points['long']:
                if vc.lower() == 'bl':
                    if 'bl' not in used:
                        adni2_long_visit_counts[diag]['BL'].append(subj)
                        used.add('bl')
                    continue
                if not vc.startswith('m'):
                    continue
                if vc not in used:
                    adni2_long_visit_counts[diag][vc].append(subj)
                    used.add(vc)


        '''
        if subj in interesting:
            mritimes = subj_points['mri']
            
            if len(mritimes) >= 5:
                continue
            
            print '\n%s' % subj
            for k,m in subj_points.items():
                if k == 'tbm':
                    print "%s: %s" % (k,m)
                elif k == 'mri':
                    diffs = [(_-m[0]).days/ 365.0 for _ in m]
                    print "%s: %s" % (k,diffs) 
                    interesting_counts[len(diffs)].append(subj)
        '''

    '''
    n_m3_mri = adni2_mri_visit_counts['SMC']['m6']
    n_m3_tbm = adni2_tbm_visit_counts['SMC']['m6']
    print len(n_m3_tbm)
    print len(n_m3_mri)
    print set(n_m3_tbm) - set(n_m3_mri)
    print set(n_m3_mri) - set(n_m3_tbm)

    print len(set(n_m3_mri) & set(n_m3_tbm))
    '''

    print "ALLS"
    print alls


    print '\n\n'
    print "ADNI1 MRI's"
    adni1table = [['Init Diagnosis', 'BL', 'Year1', 'Year2', 'Year3', 'Year4']]
    key_order = ['BL', 'year1', 'year2', 'year3', 'year4']
    for diag_row in ['N', 'SMC', 'EMCI', 'LMCI', 'AD']:
        adni1counts = adni1_mri_visit_counts[diag_row]
        cols = [diag_row] + [len(adni1counts[k]) for k in key_order]
        adni1table.append(cols)
        output_data = {'Source': 'MRI',
                       'Diagnosis': diag_row,
                       'Phase': 'ADNI1',
                       'BL': cols[1],
                       'year1': cols[2],
                       'year2': cols[3],
                       'year3': cols[4],
                       'year4': cols[5],
                       'year5': '',
                       'm3': '',
                       'm6': '',
                       'm12': '',
                       'm24': '',
                       'm36': '',
                       'm48': ''}
        rows.append(output_data)
    print tabulate(adni1table)

    print "ADNIGO/2 MRI's"
    adni2table = [['Init Diagnosis', 'BL', 'm3', 'm6', 'm12', 'm24', 'm36', 'm48']]
    key_order = ['BL', 'm3', 'm6', 'm12', 'm24', 'm36', 'm48']
    for diag_row in ['N', 'SMC', 'EMCI', 'LMCI', 'AD']:
        adni2counts = adni2_mri_visit_counts[diag_row]
        cols = [diag_row] + [len(adni2counts[k]) for k in key_order]
        adni2table.append(cols)
        output_data = {'Source': 'MRI',
                       'Diagnosis': diag_row,
                       'Phase': 'ADNIGO/2',
                       'BL': cols[1],
                       'year1': '',
                       'year2': '',
                       'year3': '',
                       'year4': '',
                       'year5': '',
                       'm3': cols[2],
                       'm6': cols[3],
                       'm12': cols[4],
                       'm24': cols[5],
                       'm36': cols[6],
                       'm48': cols[7]}
        rows.append(output_data)
    print tabulate(adni2table)

    print '\n\n'
    print "ADNI1 TBMSyn's"
    adni1table = [['Init Diagnosis', 'Year1', 'Year2', 'Year3', 'Year4']]
    key_order = ['year1', 'year2', 'year3', 'year4']
    for diag_row in ['N', 'SMC', 'EMCI', 'LMCI', 'AD']:
        adni1counts = adni1_tbm_visit_counts[diag_row]
        cols = [diag_row] + [len(adni1counts[k]) for k in key_order]
        adni1table.append(cols)
        output_data = {'Source': 'TBMSyn',
                       'Diagnosis': diag_row,
                       'Phase': 'ADNI1',
                       'BL': '',
                       'year1': cols[1],
                       'year2': cols[2],
                       'year3': cols[3],
                       'year4': cols[4],
                       'year5': '',
                       'm3': '',
                       'm6': '',
                       'm12': '',
                       'm24': '',
                       'm36': '',
                       'm48': ''}
        rows.append(output_data)
    print tabulate(adni1table)

    print "ADNIGO/2 TBMSyn's"
    adni2table = [['Init Diagnosis', 'm3', 'm6', 'm12', 'm24', 'm36', 'm48']]
    key_order = ['m3', 'm6', 'm12', 'm24', 'm36', 'm48']
    for diag_row in ['N', 'SMC', 'EMCI', 'LMCI', 'AD']:
        adni2counts = adni2_tbm_visit_counts[diag_row]
        cols = [diag_row] + [len(adni2counts[k]) for k in key_order]
        adni2table.append(cols)
        output_data = {'Source': 'TBMSyn',
                       'Diagnosis': diag_row,
                       'Phase': 'ADNIGO/2',
                       'BL': '',
                       'year1': '',
                       'year2': '',
                       'year3': '',
                       'year4': '',
                       'year5': '',
                       'm3': cols[1],
                       'm6': cols[2],
                       'm12': cols[3],
                       'm24': cols[4],
                       'm36': cols[5],
                       'm48': cols[6]}
        rows.append(output_data)
    print tabulate(adni2table)

    print '\n\n'
    print "ADNI1 BSI's"
    adni1table = [['Init Diagnosis', 'Year1', 'Year2', 'Year3', 'Year4']]
    key_order = ['year1', 'year2', 'year3', 'year4']
    for diag_row in ['N', 'SMC', 'EMCI', 'LMCI', 'AD']:
        adni1counts = adni1_bsi_visit_counts[diag_row]
        cols = [diag_row] + [len(adni1counts[k]) for k in key_order]
        adni1table.append(cols)
        output_data = {'Source': 'BSI',
                       'Diagnosis': diag_row,
                       'Phase': 'ADNI1',
                       'BL': '',
                       'year1': cols[1],
                       'year2': cols[2],
                       'year3': cols[3],
                       'year4': cols[4],
                       'year5': '',
                       'm3': '',
                       'm6': '',
                       'm12': '',
                       'm24': '',
                       'm36': '',
                       'm48': ''}
        rows.append(output_data)
    print tabulate(adni1table)

    print "ADNIGO/2 BSI's"
    adni2table = [['Init Diagnosis', 'm3', 'm6', 'm12', 'm24', 'm36', 'm48']]
    key_order = ['m3', 'm6', 'm12', 'm24', 'm36', 'm48']
    for diag_row in ['N', 'SMC', 'EMCI', 'LMCI', 'AD']:
        adni2counts = adni2_bsi_visit_counts[diag_row]
        cols = [diag_row] + [len(adni2counts[k]) for k in key_order]
        adni2table.append(cols)
        output_data = {'Source': 'BSI',
                       'Diagnosis': diag_row,
                       'Phase': 'ADNIGO/2',
                       'BL': '',
                       'year1': '',
                       'year2': '',
                       'year3': '',
                       'year4': '',
                       'year5': '',
                       'm3': cols[1],
                       'm6': cols[2],
                       'm12': cols[3],
                       'm24': cols[4],
                       'm36': cols[5],
                       'm48': cols[6]}
        rows.append(output_data)
    print tabulate(adni2table)

    print '\n\n'
    print "ADNI1 Longitudinal Freesurfer's"
    adni1table = [['Init Diagnosis', 'BL', 'Year1', 'Year2', 'Year3', 'Year4']]
    key_order = ['BL', 'year1', 'year2', 'year3', 'year4']
    for diag_row in ['N', 'SMC', 'EMCI', 'LMCI', 'AD']:
        adni1counts = adni1_long_visit_counts[diag_row]
        cols = [diag_row] + [len(adni1counts[k]) for k in key_order]
        adni1table.append(cols)
        output_data = {'Source': 'Longitudinal Freesurfer',
                       'Diagnosis': diag_row,
                       'Phase': 'ADNI1',
                       'BL': cols[1],
                       'year1': cols[2],
                       'year2': cols[3],
                       'year3': cols[4],
                       'year4': cols[5],
                       'year5': '',
                       'm3': '',
                       'm6': '',
                       'm12': '',
                       'm24': '',
                       'm36': '',
                       'm48': ''}
        rows.append(output_data)
    print tabulate(adni1table)

    print "ADNIGO/2 Longitudinal Freesurfer's"
    adni2table = [['Init Diagnosis', 'BL', 'm3', 'm6', 'm12', 'm24', 'm36', 'm48']]
    key_order = ['BL', 'm3', 'm6', 'm12', 'm24', 'm36', 'm48']
    for diag_row in ['N', 'SMC', 'EMCI', 'LMCI', 'AD']:
        adni2counts = adni2_long_visit_counts[diag_row]
        cols = [diag_row] + [len(adni2counts[k]) for k in key_order]
        adni2table.append(cols)
        output_data = {'Source': 'Longitudinal Freesurfer',
                       'Diagnosis': diag_row,
                       'Phase': 'ADNIGO/2',
                       'BL': cols[1],
                       'year1': '',
                       'year2': '',
                       'year3': '',
                       'year4': '',
                       'year5': '',
                       'm3': cols[2],
                       'm6': cols[3],
                       'm12': cols[4],
                       'm24': cols[5],
                       'm36': cols[6],
                       'm48': cols[7]}
        rows.append(output_data)
    print tabulate(adni2table)

    print '\n\n'
    print "ADNI1 Cross-sectional Freesurfer's"
    adni1table = [['Init Diagnosis', 'BL', 'Year1', 'Year2', 'Year3', 'Year4']]
    key_order = ['BL', 'year1', 'year2', 'year3', 'year4']
    for diag_row in ['N', 'SMC', 'EMCI', 'LMCI', 'AD']:
        adni1counts = adni1_cross_visit_counts[diag_row]
        cols = [diag_row] + [len(adni1counts[k]) for k in key_order]
        adni1table.append(cols)
        output_data = {'Source': 'Cross-sectional Freesurfer',
                       'Diagnosis': diag_row,
                       'Phase': 'ADNI1',
                       'BL': cols[1],
                       'year1': cols[2],
                       'year2': cols[3],
                       'year3': cols[4],
                       'year4': cols[5],
                       'year5': '',
                       'm3': '',
                       'm6': '',
                       'm12': '',
                       'm24': '',
                       'm36': '',
                       'm48': ''}
        rows.append(output_data)
    print tabulate(adni1table)

    print "ADNIGO/2 Cross-sectional Freesurfer's"
    adni2table = [['Init Diagnosis', 'BL', 'm3', 'm6', 'm12', 'm24', 'm36', 'm48']]
    key_order = ['BL', 'm3', 'm6', 'm12', 'm24', 'm36', 'm48']
    for diag_row in ['N', 'SMC', 'EMCI', 'LMCI', 'AD']:
        adni2counts = adni2_cross_visit_counts[diag_row]
        cols = [diag_row] + [len(adni2counts[k]) for k in key_order]
        adni2table.append(cols)
        output_data = {'Source': 'Cross-sectional Freesurfer',
                       'Diagnosis': diag_row,
                       'Phase': 'ADNIGO/2',
                       'BL': cols[1],
                       'year1': '',
                       'year2': '',
                       'year3': '',
                       'year4': '',
                       'year5': '',
                       'm3': cols[2],
                       'm6': cols[3],
                       'm12': cols[4],
                       'm24': cols[5],
                       'm36': cols[6],
                       'm48': cols[7]}
        rows.append(output_data)
    print tabulate(adni2table)

    dumpCSV(numerical_output, headers, rows)

    '''
    interesting_counts = dict(interesting_counts)
    for a,b in interesting_counts.iteritems():
        print a
        print b
    '''

    # count
    tbm_counts = defaultdict(int)
    long_counts = defaultdict(int)
    bsi_counts = defaultdict(int)
    mri_counts = defaultdict(int)

    for k,v in points_by_subj.iteritems():
        #print k
        diag = master_data.get(k,{}).get('Init_Diagnosis','Unknown')
        for a,b in v.iteritems():
            #print '\t%s: %s' % (a,b)
            if a == 'long':
                if len(b) >= 2:
                    long_counts[(diag, len(b)-1)] += 1
                else:
                    long_counts[(diag, 0)] += 1
            elif a == 'tbm':
                tbm_counts[(diag, len(b))] += 1
            elif a == 'bsi':
                bsi_counts[(diag, len(b))] += 1
            elif a == 'mri':
                if len(b) >= 2:
                    mri_counts[(diag, len(b)-1)] += 1
                else:
                    mri_counts[(diag, 0)] += 1


    print '\n\n'
    print "MRI"
    mri_graph = [['Init Diagnosis', 
                  '=0 followups', 
                  '=1 followups', 
                  '=2 followups', 
                  '=3 followups',
                  '=4 followups',
                  '=5 followups',
                  '=6 followups',
                  '=7 followups']]
    mri_graph.append(['CN',
                      mri_counts.get(('N',0),0),
                      mri_counts.get(('N',1),0),
                      mri_counts.get(('N',2),0),
                      mri_counts.get(('N',3),0),
                      mri_counts.get(('N',4),0),
                      mri_counts.get(('N',5),0),
                      mri_counts.get(('N',6),0),
                      mri_counts.get(('N',7),0)])
    mri_graph.append(['SMC',
                      mri_counts.get(('SMC',0),0),
                      mri_counts.get(('SMC',1),0),
                      mri_counts.get(('SMC',2),0),
                      mri_counts.get(('SMC',3),0),
                      mri_counts.get(('SMC',4),0),
                      mri_counts.get(('SMC',5),0),
                      mri_counts.get(('SMC',6),0),
                      mri_counts.get(('SMC',7),0)])
    mri_graph.append(['EMCI',
                      mri_counts.get(('EMCI',0),0),
                      mri_counts.get(('EMCI',1),0),
                      mri_counts.get(('EMCI',2),0),
                      mri_counts.get(('EMCI',3),0),
                      mri_counts.get(('EMCI',4),0),
                      mri_counts.get(('EMCI',5),0),
                      mri_counts.get(('EMCI',6),0),
                      mri_counts.get(('EMCI',7),0)])
    mri_graph.append(['LMCI',
                      mri_counts.get(('LMCI',0),0),
                      mri_counts.get(('LMCI',1),0),
                      mri_counts.get(('LMCI',2),0),
                      mri_counts.get(('LMCI',3),0),
                      mri_counts.get(('LMCI',4),0),
                      mri_counts.get(('LMCI',5),0),
                      mri_counts.get(('LMCI',6),0),
                      mri_counts.get(('LMCI',7),0)])
    mri_graph.append(['AD',
                      mri_counts.get(('AD',0),0),
                      mri_counts.get(('AD',1),0),
                      mri_counts.get(('AD',2),0),
                      mri_counts.get(('AD',3),0),
                      mri_counts.get(('AD',4),0),
                      mri_counts.get(('AD',5),0),
                      mri_counts.get(('AD',6),0),
                      mri_counts.get(('AD',7),0)])
    print tabulate(mri_graph) + '\n'
    print "TBMSyn"
    tbm_graph = [['Init Diagnosis', 
                  '=0 followups', 
                  '=1 followups', 
                  '=2 followups', 
                  '=3 followups',
                  '=4 followups',
                  '=5 followups',
                  '=6 followups',
                  '=7 followups']]
    tbm_graph.append(['CN',
                      tbm_counts.get(('N',0),0),
                      tbm_counts.get(('N',1),0),
                      tbm_counts.get(('N',2),0),
                      tbm_counts.get(('N',3),0),
                      tbm_counts.get(('N',4),0),
                      tbm_counts.get(('N',5),0),
                      tbm_counts.get(('N',6),0),
                      tbm_counts.get(('N',7),0)])
    tbm_graph.append(['SMC',
                      tbm_counts.get(('SMC',0),0),
                      tbm_counts.get(('SMC',1),0),
                      tbm_counts.get(('SMC',2),0),
                      tbm_counts.get(('SMC',3),0),
                      tbm_counts.get(('SMC',4),0),
                      tbm_counts.get(('SMC',5),0),
                      tbm_counts.get(('SMC',6),0),
                      tbm_counts.get(('SMC',7),0)])
    tbm_graph.append(['EMCI',
                      tbm_counts.get(('EMCI',0),0),
                      tbm_counts.get(('EMCI',1),0),
                      tbm_counts.get(('EMCI',2),0),
                      tbm_counts.get(('EMCI',3),0),
                      tbm_counts.get(('EMCI',4),0),
                      tbm_counts.get(('EMCI',5),0),
                      tbm_counts.get(('EMCI',6),0),
                      tbm_counts.get(('EMCI',7),0)])
    tbm_graph.append(['LMCI',
                      tbm_counts.get(('LMCI',0),0),
                      tbm_counts.get(('LMCI',1),0),
                      tbm_counts.get(('LMCI',2),0),
                      tbm_counts.get(('LMCI',3),0),
                      tbm_counts.get(('LMCI',4),0),
                      tbm_counts.get(('LMCI',5),0),
                      tbm_counts.get(('LMCI',6),0),
                      tbm_counts.get(('LMCI',7),0)])
    tbm_graph.append(['AD',
                      tbm_counts.get(('AD',0),0),
                      tbm_counts.get(('AD',1),0),
                      tbm_counts.get(('AD',2),0),
                      tbm_counts.get(('AD',3),0),
                      tbm_counts.get(('AD',4),0),
                      tbm_counts.get(('AD',5),0),
                      tbm_counts.get(('AD',6),0),
                      tbm_counts.get(('AD',7),0)])
    print tabulate(tbm_graph) + '\n'
    print "Longitudinal Freesurfer"
    long_graph = [['Init Diagnosis', 
                  '=0 followups', 
                  '=1 followups', 
                  '=2 followups', 
                  '=3 followups',
                  '=4 followups',
                  '=5 followups',
                  '=6 followups']]
    long_graph.append(['CN',
                      long_counts.get(('N',0),0),
                      long_counts.get(('N',1),0),
                      long_counts.get(('N',2),0),
                      long_counts.get(('N',3),0),
                      long_counts.get(('N',4),0),
                      long_counts.get(('N',5),0),
                      long_counts.get(('N',6),0)])
    long_graph.append(['SMC',
                      long_counts.get(('SMC',0),0),
                      long_counts.get(('SMC',1),0),
                      long_counts.get(('SMC',2),0),
                      long_counts.get(('SMC',3),0),
                      long_counts.get(('SMC',4),0),
                      long_counts.get(('SMC',5),0),
                      long_counts.get(('SMC',6),0)])
    long_graph.append(['EMCI',
                      long_counts.get(('EMCI',0),0),
                      long_counts.get(('EMCI',1),0),
                      long_counts.get(('EMCI',2),0),
                      long_counts.get(('EMCI',3),0),
                      long_counts.get(('EMCI',4),0),
                      long_counts.get(('EMCI',5),0),
                      long_counts.get(('EMCI',6),0)])
    long_graph.append(['LMCI',
                      long_counts.get(('LMCI',0),0),
                      long_counts.get(('LMCI',1),0),
                      long_counts.get(('LMCI',2),0),
                      long_counts.get(('LMCI',3),0),
                      long_counts.get(('LMCI',4),0),
                      long_counts.get(('LMCI',5),0),
                      long_counts.get(('LMCI',6),0)])
    long_graph.append(['AD',
                      long_counts.get(('AD',0),0),
                      long_counts.get(('AD',1),0),
                      long_counts.get(('AD',2),0),
                      long_counts.get(('AD',3),0),
                      long_counts.get(('AD',4),0),
                      long_counts.get(('AD',5),0),
                      long_counts.get(('AD',6),0)])
    print tabulate(long_graph) + '\n'
    print "BSI"
    bsi_graph = [['Init Diagnosis', 
                  '=0 followups', 
                  '=1 followups', 
                  '=2 followups', 
                  '=3 followups',
                  '=4 followups',
                  '=5 followups',
                  '=6 followups']]
    bsi_graph.append(['CN',
                      bsi_counts.get(('N',0),0),
                      bsi_counts.get(('N',1),0),
                      bsi_counts.get(('N',2),0),
                      bsi_counts.get(('N',3),0),
                      bsi_counts.get(('N',4),0),
                      bsi_counts.get(('N',5),0),
                      bsi_counts.get(('N',6),0)])
    bsi_graph.append(['SMC',
                      bsi_counts.get(('SMC',0),0),
                      bsi_counts.get(('SMC',1),0),
                      bsi_counts.get(('SMC',2),0),
                      bsi_counts.get(('SMC',3),0),
                      bsi_counts.get(('SMC',4),0),
                      bsi_counts.get(('SMC',5),0),
                      bsi_counts.get(('SMC',6),0)])
    bsi_graph.append(['EMCI',
                      bsi_counts.get(('EMCI',0),0),
                      bsi_counts.get(('EMCI',1),0),
                      bsi_counts.get(('EMCI',2),0),
                      bsi_counts.get(('EMCI',3),0),
                      bsi_counts.get(('EMCI',4),0),
                      bsi_counts.get(('EMCI',5),0),
                      bsi_counts.get(('EMCI',6),0)])
    bsi_graph.append(['LMCI',
                      bsi_counts.get(('LMCI',0),0),
                      bsi_counts.get(('LMCI',1),0),
                      bsi_counts.get(('LMCI',2),0),
                      bsi_counts.get(('LMCI',3),0),
                      bsi_counts.get(('LMCI',4),0),
                      bsi_counts.get(('LMCI',5),0),
                      bsi_counts.get(('LMCI',6),0)])
    bsi_graph.append(['AD',
                      bsi_counts.get(('AD',0),0),
                      bsi_counts.get(('AD',1),0),
                      bsi_counts.get(('AD',2),0),
                      bsi_counts.get(('AD',3),0),
                      bsi_counts.get(('AD',4),0),
                      bsi_counts.get(('AD',5),0),
                      bsi_counts.get(('AD',6),0)])
    print tabulate(bsi_graph) + '\n'

    print '\n\n'

    print "MRI"
    mri_graph_morethan = [['Init Diagnosis', 
                          '>=0 followups', 
                          '>=1 followups', 
                          '>=2 followups', 
                          '>=3 followups',
                          '>=4 followups',
                          '>=5 followups',
                          '>=6 followups',
                          '>=7 followups']]
    mri_graph_morethan.append(['CN',
                              sum(mri_graph[1][1:]),
                              sum(mri_graph[1][2:]),
                              sum(mri_graph[1][3:]),
                              sum(mri_graph[1][4:]),
                              sum(mri_graph[1][5:]),
                              sum(mri_graph[1][6:]),
                              sum(mri_graph[1][7:]),
                              sum(mri_graph[1][8:])])
    mri_graph_morethan.append(['SMC',
                              sum(mri_graph[2][1:]),
                              sum(mri_graph[2][2:]),
                              sum(mri_graph[2][3:]),
                              sum(mri_graph[2][4:]),
                              sum(mri_graph[2][5:]),
                              sum(mri_graph[2][6:]),
                              sum(mri_graph[2][7:]),
                              sum(mri_graph[2][8:])])
    mri_graph_morethan.append(['EMCI',
                              sum(mri_graph[3][1:]),
                              sum(mri_graph[3][2:]),
                              sum(mri_graph[3][3:]),
                              sum(mri_graph[3][4:]),
                              sum(mri_graph[3][5:]),
                              sum(mri_graph[3][6:]),
                              sum(mri_graph[3][7:]),
                              sum(mri_graph[3][8:])])
    mri_graph_morethan.append(['LMCI',
                              sum(mri_graph[4][1:]),
                              sum(mri_graph[4][2:]),
                              sum(mri_graph[4][3:]),
                              sum(mri_graph[4][4:]),
                              sum(mri_graph[4][5:]),
                              sum(mri_graph[4][6:]),
                              sum(mri_graph[4][7:]),
                              sum(mri_graph[4][8:])])
    mri_graph_morethan.append(['AD',
                              sum(mri_graph[5][1:]),
                              sum(mri_graph[5][2:]),
                              sum(mri_graph[5][3:]),
                              sum(mri_graph[5][4:]),
                              sum(mri_graph[5][5:]),
                              sum(mri_graph[5][6:]),
                              sum(mri_graph[5][7:]),
                              sum(mri_graph[5][8:])])
    print tabulate(mri_graph_morethan) + '\n'
    print "TBMSyn"
    tbm_graph_morethan = [['Init Diagnosis', 
                  '>=0 followups', 
                  '>=1 followups', 
                  '>=2 followups', 
                  '>=3 followups',
                  '>=4 followups',
                  '>=5 followups',
                  '>=6 followups',
                  '>=7 followups']]
    tbm_graph_morethan.append(['CN',
                              sum(tbm_graph[1][1:]),
                              sum(tbm_graph[1][2:]),
                              sum(tbm_graph[1][3:]),
                              sum(tbm_graph[1][4:]),
                              sum(tbm_graph[1][5:]),
                              sum(tbm_graph[1][6:]),
                              sum(tbm_graph[1][7:]),
                              sum(tbm_graph[1][8:])])
    tbm_graph_morethan.append(['SMC',
                              sum(tbm_graph[2][1:]),
                              sum(tbm_graph[2][2:]),
                              sum(tbm_graph[2][3:]),
                              sum(tbm_graph[2][4:]),
                              sum(tbm_graph[2][5:]),
                              sum(tbm_graph[2][6:]),
                              sum(tbm_graph[2][7:]),
                              sum(tbm_graph[2][8:])])
    tbm_graph_morethan.append(['EMCI',
                              sum(tbm_graph[3][1:]),
                              sum(tbm_graph[3][2:]),
                              sum(tbm_graph[3][3:]),
                              sum(tbm_graph[3][4:]),
                              sum(tbm_graph[3][5:]),
                              sum(tbm_graph[3][6:]),
                              sum(tbm_graph[3][7:]),
                              sum(tbm_graph[3][8:])])
    tbm_graph_morethan.append(['LMCI',
                              sum(tbm_graph[4][1:]),
                              sum(tbm_graph[4][2:]),
                              sum(tbm_graph[4][3:]),
                              sum(tbm_graph[4][4:]),
                              sum(tbm_graph[4][5:]),
                              sum(tbm_graph[4][6:]),
                              sum(tbm_graph[4][7:]),
                              sum(tbm_graph[4][8:])])
    tbm_graph_morethan.append(['AD',
                              sum(tbm_graph[5][1:]),
                              sum(tbm_graph[5][2:]),
                              sum(tbm_graph[5][3:]),
                              sum(tbm_graph[5][4:]),
                              sum(tbm_graph[5][5:]),
                              sum(tbm_graph[5][6:]),
                              sum(tbm_graph[5][7:]),
                              sum(tbm_graph[5][8:])])
    print tabulate(tbm_graph_morethan) + '\n'
    print "Longitudinal Freesurfer"
    long_graph_morethan = [['Init Diagnosis', 
                              '>=0 followups', 
                              '>=1 followups', 
                              '>=2 followups', 
                              '>=3 followups',
                              '>=4 followups',
                              '>=5 followups',
                              '>=6 followups']]
    long_graph_morethan.append(['CN',
                              sum(long_graph[1][1:]),
                              sum(long_graph[1][2:]),
                              sum(long_graph[1][3:]),
                              sum(long_graph[1][4:]),
                              sum(long_graph[1][5:]),
                              sum(long_graph[1][6:]),
                              sum(long_graph[1][7:])])
    long_graph_morethan.append(['SMC',
                              sum(long_graph[2][1:]),
                              sum(long_graph[2][2:]),
                              sum(long_graph[2][3:]),
                              sum(long_graph[2][4:]),
                              sum(long_graph[2][5:]),
                              sum(long_graph[2][6:]),
                              sum(long_graph[2][7:])])
    long_graph_morethan.append(['EMCI',
                              sum(long_graph[3][1:]),
                              sum(long_graph[3][2:]),
                              sum(long_graph[3][3:]),
                              sum(long_graph[3][4:]),
                              sum(long_graph[3][5:]),
                              sum(long_graph[3][6:]),
                              sum(long_graph[3][7:])])
    long_graph_morethan.append(['LMCI',
                              sum(long_graph[4][1:]),
                              sum(long_graph[4][2:]),
                              sum(long_graph[4][3:]),
                              sum(long_graph[4][4:]),
                              sum(long_graph[4][5:]),
                              sum(long_graph[4][6:]),
                              sum(long_graph[4][7:])])
    long_graph_morethan.append(['AD',
                              sum(long_graph[5][1:]),
                              sum(long_graph[5][2:]),
                              sum(long_graph[5][3:]),
                              sum(long_graph[5][4:]),
                              sum(long_graph[5][5:]),
                              sum(long_graph[5][6:]),
                              sum(long_graph[5][7:])])
    print tabulate(long_graph_morethan) + '\n'
    print "BSI"
    bsi_graph_morethan = [['Init Diagnosis', 
                  '>=0 followups', 
                  '>=1 followups', 
                  '>=2 followups', 
                  '>=3 followups',
                  '>=4 followups',
                  '>=5 followups',
                  '>=6 followups']]
    bsi_graph_morethan.append(['CN',
                              sum(bsi_graph[1][1:]),
                              sum(bsi_graph[1][2:]),
                              sum(bsi_graph[1][3:]),
                              sum(bsi_graph[1][4:]),
                              sum(bsi_graph[1][5:]),
                              sum(bsi_graph[1][6:]),
                              sum(bsi_graph[1][7:])])
    bsi_graph_morethan.append(['SMC',
                              sum(bsi_graph[2][1:]),
                              sum(bsi_graph[2][2:]),
                              sum(bsi_graph[2][3:]),
                              sum(bsi_graph[2][4:]),
                              sum(bsi_graph[2][5:]),
                              sum(bsi_graph[2][6:]),
                              sum(bsi_graph[2][7:])])
    bsi_graph_morethan.append(['EMCI',
                              sum(bsi_graph[3][1:]),
                              sum(bsi_graph[3][2:]),
                              sum(bsi_graph[3][3:]),
                              sum(bsi_graph[3][4:]),
                              sum(bsi_graph[3][5:]),
                              sum(bsi_graph[3][6:]),
                              sum(bsi_graph[3][7:])])
    bsi_graph_morethan.append(['LMCI',
                              sum(bsi_graph[4][1:]),
                              sum(bsi_graph[4][2:]),
                              sum(bsi_graph[4][3:]),
                              sum(bsi_graph[4][4:]),
                              sum(bsi_graph[4][5:]),
                              sum(bsi_graph[4][6:]),
                              sum(bsi_graph[4][7:])])
    bsi_graph_morethan.append(['AD',
                              sum(bsi_graph[5][1:]),
                              sum(bsi_graph[5][2:]),
                              sum(bsi_graph[5][3:]),
                              sum(bsi_graph[5][4:]),
                              sum(bsi_graph[5][5:]),
                              sum(bsi_graph[5][6:]),
                              sum(bsi_graph[5][7:])])
    print tabulate(bsi_graph_morethan) + '\n'



if __name__ == "__main__":

    include_failed = True

    # Input/output/lookup files
    master_file = "../FDG_AV45_COGdata_07_15_15.csv"
    registry_file = "../docs/registry_clean.csv"
    pet_meta_file = "../docs/PET_META_LIST_edited.csv"
    numerical_output = "../longitudinal_mri_counts.csv"
    
    mri_meta_file = "../docs/MPRAGEMETA.csv"
    #mri_meta_file = "../docs/idaSearch_7_09_2015.csv"
    
    # BSI file
    bsi_file = "../mr_docs/Fox/FOXLABBSI_04_30_15.csv"
    # long freesurfer file
    longfree_file = '../mr_docs/UCSF/longitudinal/UCSFFSL51Y1_08_01_14.csv'
    longfree_adni1_file = "../mr_docs/UCSF/longitudinal/UCSFFSL_05_20_15_ADNI1.csv"
    # cross freesurfer file
    crossfree_file = "../mr_docs/UCSF/cross_section/UCSFFSX51_05_20_15.csv"
    # TBMsyn file
    tbm_file = '../mr_docs/Mayo/MAYOADIRL_MRI_TBMSYN_05_07_15.csv'

    pet_data = importPetMETA(pet_meta_file)
    print 'PET patients: %s' % len(pet_data)
    bsi_data = importBSI(bsi_file, include_failed=include_failed)
    print 'BSI patients: %s' % len(bsi_data)
    longfree_data = importLongitudinalFreesurfer(longfree_file, include_failed=include_failed)
    print 'Longfree patients: %s' % len(longfree_data)
    longfree_data_adni1 = importLongitudinalFreesurfer(longfree_adni1_file, include_failed=include_failed)
    print 'Longfree ADNI1 patients: %s' % len(longfree_data_adni1)
    crossfree_data = importCrossSectionFreesurfer(crossfree_file, include_failed=include_failed)
    print "Crossfree patients: %s" % len(crossfree_data)
    tbm_data = importTBMSyn(tbm_file)
    print 'TBM patients: %s' % len(tbm_data)
    mri_data = importMRI(mri_meta_file)
    print 'MRI patients: %s' % len(mri_data)
    master_data = importMaster(master_file)

    avai_points = checkAvailablePointsPerSubject(pet_data, bsi_data, longfree_data, longfree_data_adni1, crossfree_data, tbm_data, mri_data, master_data, numerical_output)





