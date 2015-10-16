from utils import *
import numpy as np
from collections import defaultdict
from scipy.stats import spearmanr
import itertools
import matplotlib.pyplot as plt
import pylab
from scipy.stats import norm, mannwhitneyu, linregress
from ggplot import *

GROUPS = {'increasing_high': {'AD': [4657, 4660, 4153, 4672, 4696, 4195, 4252, 4258, 4477, 4494, 4500, 4568, 4089], 
                              'EMCI': [4614, 2063, 4624, 2073, 4635, 2077, 2087, 4146, 4149, 4674, 4679, 4680, 2121, 2133, 
                                       4623, 2155, 4764, 2167, 4799, 4235, 2195, 2196, 4259, 2213, 2216, 4272, 2225, 4805, 
                                       4809, 4816, 4309, 2264, 2274, 4844, 4858, 4351, 4868, 4893, 2336, 4404, 4405, 4415, 
                                       4419, 2373, 2381, 4944, 4947, 2190, 2390, 2391, 4447, 2403, 2022, 4974, 2109, 4467, 
                                       4473, 5000, 5007, 4005, 4007, 2037, 4547, 5066, 2045, 2047], 
                              'LMCI': [4611, 4782, 4631, 4636, 4654, 4675, 4170, 4689, 4189, 4197, 4712, 4713, 112, 4714, 
                                       126, 4736, 4715, 142, 4240, 4757, 4796, 4042, 702, 4287, 800, 4294, 4815, 729, 1030, 
                                       1246, 227, 4857, 4363, 4877, 4203, 1300, 4057, 4902, 1074, 4406, 4925, 4414, 4929, 
                                       1346, 4945, 4444, 4035, 378, 906, 4502, 925, 4510, 4521, 4015, 5047, 4538, 4542, 
                                       4243, 4562, 994, 4582], 
                              'N': [4612, 4616, 4120, 31, 545, 4645, 4151, 61, 4176, 4179, 4198, 106, 4225, 130, 166, 173, 
                                    4270, 4278, 4290, 4291, 717, 1232, 230, 4339, 4343, 778, 4385, 4386, 4400, 4422, 842, 
                                    4482, 920, 4003, 4014, 1098, 972, 981, 4566, 985, 4081]}, 
          'stable_high': {'AD': [4501, 4641, 4526, 4192, 4591, 4211, 4215], 
                          'EMCI': [2307, 2055, 2068, 4891, 4765, 2079, 4128, 4392, 2100, 4661, 2376, 2380, 4557, 2130, 4188, 2142], 
                          'LMCI': [4595, 4888, 4250, 1186, 4263, 1066, 4167, 1326, 4936, 4531, 4030, 4928, 4162, 709, 4807, 331, 
                                   4430, 4311, 1117, 887, 4584, 4458, 1004, 1106, 4596, 4729], 
                          'N': [4100, 4371, 4262, 4320, 1190]}, 
          'stable_low': {'AD': [4009], 
                         'EMCI': [2304, 2052, 2074, 4127, 2083, 4160, 4678, 2119, 2146, 4199, 2153, 2164, 4220, 2191, 
                                  2072, 2200, 2201, 2208, 4380, 2220, 4271, 2245, 2247, 4299, 4813, 4301, 4310, 4312, 
                                  2301, 2347, 4356, 4876, 4883, 2332, 4917, 2374, 2378, 2395, 4455, 4476, 4268, 4513, 
                                  4036, 4004, 2010, 4063, 4072, 2027, 2031, 2036], 
                         'LMCI': [4626, 1052, 4138, 4169, 4187, 4722, 1155, 135, 4750, 4244, 4293, 4769, 4806, 200, 
                                  722, 746, 4989, 4381, 4395, 4960, 4842, 4543], 
                         'N': [4103, 4105, 4075, 21, 4632, 4121, 4148, 56, 59, 4158, 4164, 74, 4173, 89, 4200, 4208, 
                               680, 4269, 4275, 1206, 4292, 4832, 4345, 4350, 4352, 4365, 272, 4369, 4372, 4382, 4384, 
                               1169, 315, 4491, 4427, 4429, 4446, 4453, 4598, 4469, 4483, 907, 4499, 4505, 416, 419, 
                               4516, 4021, 4084, 4255, 4545, 4037, 4041, 4043, 4559, 4560, 4579, 4599, 4090, 4604, 4607]}, 
          'increasing_low': {'AD': [4172, 4676], 
                             'EMCI': [2061, 4133, 2093, 4143, 4159, 2116, 4168, 2123, 4621, 4281, 4184, 4185, 2148, 2151, 
                                      4212, 2180, 4742, 2183, 2187, 4237, 4216, 2219, 4780, 2233, 2234, 4285, 2238, 2239, 
                                      2249, 2263, 4328, 4331, 4360, 4874, 4919, 4898, 2357, 2360, 2363, 4926, 2367, 4417, 
                                      2379, 4434, 2394, 4443, 2398, 2405, 2407, 4468, 4986, 4512, 4012, 4536, 2002, 4571, 
                                      2018, 4594, 2042], 
                             'LMCI': [4114, 4115, 1045, 4061, 4122, 4869, 42, 1122, 1118, 4210, 4723, 4214, 4229, 668, 1187, 
                                      679, 4777, 4784, 4300, 4817, 225, 4354, 4873, 4741, 1318, 1352, 867, 1414, 4489, 1418, 
                                      919, 908, 1427, 408, 4767, 4653, 4094], 
                             'N': [4620, 23, 4637, 4643, 4649, 4150, 58, 69, 4177, 602, 610, 618, 113, 4213, 120, 4218, 123, 
                                   4222, 4739, 4234, 4119, 4762, 4254, 159, 685, 4795, 4313, 741, 1261, 751, 767, 1280, 260, 
                                   4357, 4878, 4367, 4387, 4388, 4389, 4393, 4396, 301, 4399, 311, 4410, 4441, 4448, 4466, 
                                   4843, 4485, 4496, 4503, 923, 413, 926, 5023, 4520, 4082, 4018, 4020, 4028, 4032, 4552, 
                                   4060, 4576, 4580, 4585, 4586, 4076, 4086, 1016]}}

# don't consider regions that include white matter, cerebellum, brainstem
REGION_BLACKLIST = [30,62,2,80,81,82,41,77,251,252,253,254,255,1004,2004,85,16,7,8,46,47]
REFERENCE_REGIONS = [7,8,46,47]

def getRegionAverageRankings(sorted_uptakes):
    sorted_regions = [[k for k,v in sorted_values] for sorted_values in sorted_uptakes]
    rankings = defaultdict(list)
    for sorted_list in sorted_uptakes:
        region_ranks = [k for k,v in sorted_list]
        for rank, region in enumerate(region_ranks):
            rankings[region].append(rank)
    avg_rankings = {}
    for k in rankings:
        mean_rank = np.mean(rankings[k])
        std_rank = np.std(rankings[k])
        avg_rankings[k] = (mean_rank, std_rank)
    sorted_avg_rankings = sorted(avg_rankings.iteritems(), key=lambda x:x[1][0])
    return sorted_avg_rankings

def getRankingSimilarity(avg_rank, sorted_uptakes, n_cutoff=None):
    sorted_regions = [[k for k,v in sorted_values] for sorted_values in sorted_uptakes]
    sims = []
    ranking_truth = [_[0] for _ in avg_rank]
    translation = {rt:i for i,rt in enumerate(ranking_truth)}
    translated_truth = [translation[_] for _ in ranking_truth]
    if n_cutoff is None:
        n_cutoff = len(translated_truth)
    for ranking_tocompare in sorted_regions:
        # get score
        translated_tocompare = [translation[_] for _ in ranking_tocompare]
        rho, p = spearmanr(translated_truth[:n_cutoff], translated_tocompare[:n_cutoff])
        sims.append(rho)
    avg_sim = np.mean(sims)
    return avg_sim

def getMeanAveragePrecision(avg_rank, sorted_uptakes, k):
    sorted_regions = [[name for name,v in sorted_values] for sorted_values in sorted_uptakes]

    actual = [_[0] for _ in avg_rank][:k]
    actuals = [actual]*len(sorted_regions)
    sim = mapk(actuals, sorted_regions, k=k)
    return sim

    '''
    sims = []
    for a,b in itertools.permutations(sorted_regions, 2):
        sim = apk(a[:k],b[:k],k=k)
        sims.append(sim)
    return np.mean(sims)
    '''

def removeBlacklistedGroups(data_bl, data_scan2, data_scan3, index_lookup, suvr=True):
    ref_regions = [k for k,v in index_lookup.iteritems() if len(list(set(REFERENCE_REGIONS) & set(list(v)))) > 0]
    regions_to_remove = [k for k,v in index_lookup.iteritems() if len(list(set(REGION_BLACKLIST) & set(list(v)))) > 0]
    # add ref regions
    regions_to_remove.append('whole_cerebellum')
    index_lookup = {k:v for k,v in index_lookup.iteritems() if k not in regions_to_remove}
    for rid,val in data_bl.iteritems():
        ref_value = 1.0
        if suvr:
            ref_value = val.get('whole_cerebellum',1.0)
        new_val = {k:v/ref_value for k,v in val.iteritems() if k not in regions_to_remove}
        data_bl[rid] = new_val
    for rid,val in data_scan2.iteritems():
        ref_value = 1.0
        if suvr:
            ref_value = val.get('whole_cerebellum',1.0)
        new_val = {k:v/ref_value for k,v in val.iteritems() if k not in regions_to_remove}
        data_scan2[rid] = new_val
    for rid,val in data_scan3.iteritems():
        ref_value = 1.0
        if suvr:
            ref_value = val.get('whole_cerebellum',1.0)
        new_val = {k:v/ref_value for k,v in val.iteritems() if k not in regions_to_remove}
        data_scan3[rid] = new_val
    return data_bl, data_scan2, data_scan3, index_lookup


def valuesByFreesurferRegion(values_by_group, index_lookup):
    by_freesurfer_region = {}
    for group, val in values_by_group.iteritems():
        fs_indices = index_lookup[group]
        for fsi in fs_indices:
            by_freesurfer_region[fsi] = val
    return by_freesurfer_region

def regionEffectSizesBetweenGroups(group_one, group_two, data, index_lookup):
    # calculate group one suvr distr
    group_uptakes_one = {k:[] for k in index_lookup}
    for rid in group_one:
        if rid not in data:
            continue
        for k in data[rid]:
            group_uptakes_one[k].append(float(data[rid][k]))
    # calculate group two suvr distr
    group_uptakes_two = {k:[] for k in index_lookup}
    for rid in group_two:
        if rid not in data:
            continue
        for k in data[rid]:
            group_uptakes_two[k].append(float(data[rid][k]))
    # calculate effect sizes
    group_effects = {}
    for k in index_lookup:
        one_uptake = group_uptakes_one[k]
        two_uptake = group_uptakes_two[k]

        u_uptake, pvalue_uptake = mannwhitneyu(one_uptake, two_uptake, use_continuity=True)
        u_max_uptake = len(one_uptake) * len(two_uptake)
        rank_biserial_uptake = 1.0 - (2*u_uptake/u_max_uptake)

        to_save = {'pvalue': pvalue_uptake,'rank_biserial': rank_biserial_uptake}
        group_effects[k] = to_save
    return group_effects


def regionalEffectSizes(subj_group, data_prior, data_post, index_lookup):
    # calculate prior suvr/rank distr
    sorted_uptakes_prior = []
    group_uptakes_prior = {k:[] for k in index_lookup}
    group_ranks_prior = {k:[] for k in index_lookup}
    for rid in subj_group:
        if rid in data_prior:
            sorted_uptakes_prior.append(sorted(data_prior[rid].iteritems(), key=lambda x: x[1], reverse=True))
            for k in data_prior[rid]:
                group_uptakes_prior[k].append(float(data_prior[rid][k]))
    uptakes_prior = {k: (np.mean(v),np.std(v)) for k,v in group_uptakes_prior.iteritems()}
    # calculate prior rank distr
    for sorted_list in sorted_uptakes_prior:
        region_ranks = [k for k,v in sorted_list]
        for rank, region in enumerate(region_ranks):
            group_ranks_prior[region].append(rank)
    ranks_prior = {k: (np.mean(v),np.std(v)) for k,v in group_ranks_prior.iteritems()}

    # calculate post suvr/rank distr
    sorted_uptakes_post = []
    group_uptakes_post = {k:[] for k in index_lookup}
    group_ranks_post = {k:[] for k in index_lookup}
    for rid in subj_group:
        if rid in data_post:
            sorted_uptakes_post.append(sorted(data_post[rid].iteritems(), key=lambda x: x[1], reverse=True))
            for k in data_post[rid]:
                group_uptakes_post[k].append(float(data_post[rid][k]))
    uptakes_post = {k: (np.mean(v),np.std(v)) for k,v in group_uptakes_post.iteritems()}
    # calculate post rank distr
    for sorted_list in sorted_uptakes_post:
        region_ranks = [k for k,v in sorted_list]
        for rank, region in enumerate(region_ranks):
            group_ranks_post[region].append(rank)
    ranks_post = {k: (np.mean(v),np.std(v)) for k,v in group_ranks_post.iteritems()}

    # calculate effect sizes
    group_effects = {}
    for k in index_lookup:
        prior_uptake = group_uptakes_prior[k]
        post_uptake = group_uptakes_post[k]
        prior_rank = group_ranks_prior[k]
        post_rank = group_ranks_post[k]

        u_uptake, pvalue_uptake = mannwhitneyu(prior_uptake, post_uptake, use_continuity=True)
        u_max_uptake = len(prior_uptake) * len(post_uptake)
        rank_biserial_uptake = 1.0 - (2*u_uptake/u_max_uptake)

        u_rank, pvalue_rank = mannwhitneyu(prior_rank, post_rank, use_continuity=True)
        u_max_rank = len(prior_rank) * len(post_rank)
        rank_biserial_rank = 1.0 - (2*u_rank/u_max_rank)

        to_save = {'uptake_effect': {'pvalue': pvalue_uptake,
                                     'rank_biserial': rank_biserial_uptake},
                   'rank_effect': {'pvalue': pvalue_rank,
                                   'rank_biserial': rank_biserial_rank}}
        print to_save
        group_effects[k] = to_save

    return (uptakes_prior, ranks_prior, uptakes_post, ranks_post, group_effects)

def rankSubjects(key_groups, subj, data):
    avg_ranks = {}
    for rid in subj:
        if rid in data:
            sorted_list = sorted(data[rid].iteritems(), key=lambda x: x[1], reverse=True)
            key_group_ranks = [i for i,(k,v) in enumerate(sorted_list) if k in key_groups]
            avg_rank = np.mean(key_group_ranks)
            avg_ranks[rid] = avg_rank
    return avg_ranks

if __name__ == "__main__":
    lut_file = "../FreeSurferColorLUT.txt"
    lut_table = importFreesurferLookup(lut_file)
    master_file = '../FDG_AV45_COGdata_10_15_15.csv'
    master_data = importMaster(master_file)
    diags = {}
    for rid, row in master_data.iteritems():
        diag = row['Init_Diagnosis'].strip()
        diags[rid] = diag

    # # agghigh raw results
    # bl_file = '../raw_agghigh_output_BL.mat'
    # scan2_file = '../raw_agghigh_output_Scan2.mat'
    # scan3_file = '../raw_agghigh_output_Scan3.mat'
    # data_bl, data_scan2, data_scan3, index_lookup = parseRawRousset(bl_file, scan2_file, scan3_file)
    # data_bl, data_scan2, data_scan3, index_lookup = removeBlacklistedGroups(data_bl, data_scan2, data_scan3, index_lookup, suvr=True)

    # # group 4 raw results
    # bl_file = '../raw_group4_output_BL.mat'
    # scan2_file = '../raw_group4_output_Scan2.mat'
    # scan3_file = '../raw_group4_output_Scan3.mat'
    # data_bl, data_scan2, data_scan3, index_lookup = parseRawRousset(bl_file, scan2_file, scan3_file)
    # data_bl, data_scan2, data_scan3, index_lookup = removeBlacklistedGroups(data_bl, data_scan2, data_scan3, index_lookup, suvr=True)

    # tp-specific results
    av45_file = "../output/UCBERKELEYAV45_09_25_15_extra.csv"
    registry_file = "../docs/registry_clean.csv"
    data_bl, data_scan2, data_scan3, index_lookup = parseRawAV45Output(av45_file, registry_file, lut_file)
    data_bl, data_scan2, data_scan3, index_lookup = removeBlacklistedGroups(data_bl, data_scan2, data_scan3, index_lookup, suvr=True)

    data_prior = data_bl
    data_post = data_scan2
    data_post.update(data_scan3)

    all_data = []
    allrids = list(set(data_scan2.keys() + data_scan3.keys()))
    for rid in allrids:
        datarow = {'rid': rid,
                   'diag': diags.get(rid,'')}
        bl_data = data_bl[rid]
        scan2_data = data_scan2.get(rid,{})
        scan3_data = data_scan3.get(rid,{})
        long_data = {}
        if scan2_data:
            yrs = yrs = float(master_data[rid]['AV45_1_2_Diff (Yrs)'])
            for k in index_lookup:
                long_data[k] = scan2_data[k] - bl_data[k] / yrs
        if scan3_data:
            yrs = float(master_data[rid]['AV45_1_3_Diff (yrs)'])
            for k in index_lookup:
                long_data[k] = scan2_data[k] - bl_data[k] / yrs
        for k,v in bl_data.iteritems():
            datarow['%s_bl' % k ] = v
        for k,v in long_data.iteritems():
            datarow['%s_change' % k] = v
        all_data.append(datarow)
    df = pd.DataFrame(all_data)
    print df
    outfile = open('../regional_bl_vs_change.json','w')
    outfile.write(df.to_json())
    sys.exit(1)


    line_data = {k:{'regions': [lut_table[_] for _ in v]} for k,v in index_lookup.iteritems()}

    # # identify subjects who rank first in key_groups out of subj_group
    # subj = GROUPS['increasing_low']['N'] + GROUPS['stable_low']['N']
    # key_groups = ['Cluster17', 'Cluster36', 'Cluster45']
    # avg_ranks = rankSubjects(key_groups, subj, data_scan2)
    # for k,v in sorted(avg_ranks.iteritems(), key=lambda x: x[1]):
    #     print "%s: %s" % (k,v)
    # sys.exit(1)

    # low normals
    subj_group = GROUPS['increasing_low']['N'] + GROUPS['stable_low']['N']
    (uptakes_prior, ranks_prior, uptakes_post, ranks_post, group_effects) = regionalEffectSizes(subj_group, data_prior, data_post, index_lookup)
    for k,v in uptakes_prior.iteritems():
        line_data[k]['LOW_N_%s' % 'uptakes_prior_mean'] = v[0]
        line_data[k]['LOW_N_%s' % 'uptakes_prior_std'] = v[1]
    for k,v in ranks_prior.iteritems():
        line_data[k]['LOW_N_%s' % 'ranks_prior_mean'] = v[0]
        line_data[k]['LOW_N_%s' % 'ranks_prior_std'] = v[1]
    for k,v in uptakes_post.iteritems():
        line_data[k]['LOW_N_%s' % 'uptakes_post_mean'] = v[0]
        line_data[k]['LOW_N_%s' % 'uptakes_post_std'] = v[1]
    for k,v in ranks_post.iteritems():
        line_data[k]['LOW_N_%s' % 'ranks_post_mean'] = v[0]
        line_data[k]['LOW_N_%s' % 'ranks_post_std'] = v[1]
    for k,v in group_effects.iteritems():
        uptake_effect = v['uptake_effect']
        for eff_k, eff_v in uptake_effect.iteritems():
            line_data[k]['LOW_N_uptake_effect_%s' % eff_k] = eff_v
        rank_effect =v['rank_effect']
        for eff_k, eff_v in rank_effect.iteritems():
            line_data[k]['LOW_N_rank_effect_%s' % eff_k] = eff_v

    # low, stable normals
    subj_group = GROUPS['stable_low']['N']
    (uptakes_prior, ranks_prior, uptakes_post, ranks_post, group_effects) = regionalEffectSizes(subj_group, data_prior, data_post, index_lookup)
    for k,v in uptakes_prior.iteritems():
        line_data[k]['STABLE_LOW_N_%s' % 'uptakes_prior_mean'] = v[0]
        line_data[k]['STABLE_LOW_N_%s' % 'uptakes_prior_std'] = v[1]
    for k,v in ranks_prior.iteritems():
        line_data[k]['STABLE_LOW_N_%s' % 'ranks_prior_mean'] = v[0]
        line_data[k]['STABLE_LOW_N_%s' % 'ranks_prior_std'] = v[1]
    for k,v in uptakes_post.iteritems():
        line_data[k]['STABLE_LOW_N_%s' % 'uptakes_post_mean'] = v[0]
        line_data[k]['STABLE_LOW_N_%s' % 'uptakes_post_std'] = v[1]
    for k,v in ranks_post.iteritems():
        line_data[k]['STABLE_LOW_N_%s' % 'ranks_post_mean'] = v[0]
        line_data[k]['STABLE_LOW_N_%s' % 'ranks_post_std'] = v[1]
    for k,v in group_effects.iteritems():
        uptake_effect = v['uptake_effect']
        for eff_k, eff_v in uptake_effect.iteritems():
            line_data[k]['STABLE_LOW_N_uptake_effect_%s' % eff_k] = eff_v
        rank_effect =v['rank_effect']
        for eff_k, eff_v in rank_effect.iteritems():
            line_data[k]['STABLE_LOW_N_rank_effect_%s' % eff_k] = eff_v

    #low, increasing normals
    subj_group = GROUPS['increasing_low']['N']
    (uptakes_prior, ranks_prior, uptakes_post, ranks_post, group_effects) = regionalEffectSizes(subj_group, data_prior, data_post, index_lookup)
    for k,v in uptakes_prior.iteritems():
        line_data[k]['INC_LOW_N_%s' % 'uptakes_prior_mean'] = v[0]
        line_data[k]['INC_LOW_N_%s' % 'uptakes_prior_std'] = v[1]
    for k,v in ranks_prior.iteritems():
        line_data[k]['INC_LOW_N_%s' % 'ranks_prior_mean'] = v[0]
        line_data[k]['INC_LOW_N_%s' % 'ranks_prior_std'] = v[1]
    for k,v in uptakes_post.iteritems():
        line_data[k]['INC_LOW_N_%s' % 'uptakes_post_mean'] = v[0]
        line_data[k]['INC_LOW_N_%s' % 'uptakes_post_std'] = v[1]
    for k,v in ranks_post.iteritems():
        line_data[k]['INC_LOW_N_%s' % 'ranks_post_mean'] = v[0]
        line_data[k]['INC_LOW_N_%s' % 'ranks_post_std'] = v[1]
    for k,v in group_effects.iteritems():
        uptake_effect = v['uptake_effect']
        for eff_k, eff_v in uptake_effect.iteritems():
            line_data[k]['INC_LOW_N_uptake_effect_%s' % eff_k] = eff_v
        rank_effect =v['rank_effect']
        for eff_k, eff_v in rank_effect.iteritems():
            line_data[k]['INC_LOW_N_rank_effect_%s' % eff_k] = eff_v

    #effects between INC_LOW and STABLE_LOW
    inc_vs_low_effects = regionEffectSizesBetweenGroups(GROUPS['increasing_low']['N'], 
                                                        GROUPS['stable_low']['N'], 
                                                        data_prior, 
                                                        index_lookup)
    for k,v in inc_vs_low_effects.iteritems():
        for eff_k, eff_v in v.iteritems():
            line_data[k]['INC_VS_STABLE_uptake_effect_%s' % eff_k] = eff_v

    columns = ['Name'] + sorted(line_data.values()[0].keys())
    lines = []
    for k,v in line_data.iteritems():
        v['Name'] = k
        lines.append(v)
    dumpCSV('../regional_effect_sizes_tpspecific.csv', columns, lines)


    # # write out results to matfile
    out_file = '../output/fake_aparc_inputs/REGION_INDICES.mat'
    values = {k: i+1 for i, (k,v) in enumerate(line_data.iteritems())}
    sio.savemat(out_file, {'regions': [{'inds': index_lookup[k], 'val': v} for k,v in values.iteritems()]})
    out_file = '../output/fake_aparc_inputs/ALL_LOW_N_mean_suvr.mat'
    values = {k: v['LOW_N_uptakes_prior_mean'] for k,v in line_data.iteritems()}
    sio.savemat(out_file, {'regions': [{'inds': index_lookup[k], 'val': v} for k,v in values.iteritems()]})
    out_file = '../output/fake_aparc_inputs/INC_LOW_N_mean_suvr.mat'
    values = {k: v['INC_LOW_N_uptakes_prior_mean'] for k,v in line_data.iteritems()}
    sio.savemat(out_file, {'regions': [{'inds': index_lookup[k], 'val': v} for k,v in values.iteritems()]})
    out_file = '../output/fake_aparc_inputs/STABLE_LOW_N_mean_suvr.mat'
    values = {k: v['STABLE_LOW_N_uptakes_prior_mean'] for k,v in line_data.iteritems()}
    sio.savemat(out_file, {'regions': [{'inds': index_lookup[k], 'val': v} for k,v in values.iteritems()]})
    out_file = '../output/fake_aparc_inputs/ALL_LOW_N_long_effect.mat'
    values = {k: v['LOW_N_uptake_effect_rank_biserial'] for k,v in line_data.iteritems()}
    sio.savemat(out_file, {'regions': [{'inds': index_lookup[k], 'val': v} for k,v in values.iteritems()]})
    out_file = '../output/fake_aparc_inputs/INC_LOW_N_long_effect.mat'
    values = {k: v['INC_LOW_N_uptake_effect_rank_biserial'] for k,v in line_data.iteritems()}
    sio.savemat(out_file, {'regions': [{'inds': index_lookup[k], 'val': v} for k,v in values.iteritems()]})
    out_file = '../output/fake_aparc_inputs/STABLE_LOW_N_long_effect.mat'
    values = {k: v['STABLE_LOW_N_uptake_effect_rank_biserial'] for k,v in line_data.iteritems()}
    sio.savemat(out_file, {'regions': [{'inds': index_lookup[k], 'val': v} for k,v in values.iteritems()]})
    out_file = '../output/fake_aparc_inputs/INC_VS_STABLE_bl_effect.mat'
    values = {k: v['INC_VS_STABLE_uptake_effect_rank_biserial'] for k,v in line_data.iteritems()}
    sio.savemat(out_file, {'regions': [{'inds': index_lookup[k], 'val': v} for k,v in values.iteritems()]})
