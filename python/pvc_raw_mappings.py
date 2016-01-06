from utils import *
import numpy as np
from collections import defaultdict
from scipy.stats import spearmanr
import itertools
import matplotlib.pyplot as plt
import pylab
from scipy.stats import norm, mannwhitneyu, linregress
from ggplot import *
import pandas as pd

REGION_BLACKLIST = [0,30,62,80,81,82,77,251,252,253,254,255,1000,2000,1004,2004,85,24,14,15,72,4,43,75,76]

def generateBathroomTilePlot(bl_vs_change_json):
    df = pd.read_json(bl_vs_change_json)
    summary_regions = ['ctx-lh-parsorbitalis','ctx-rh-parsorbitalis','ctx-rh-lateralorbitofrontal',
                       'ctx-lh-lateralorbitofrontal','ctx-rh-frontalpole','ctx-rh-parstriangularis',
                       'ctx-lh-frontalpole','ctx-lh-parstriangularis','ctx-lh-caudalanteriorcingulate',
                       'ctx-rh-rostralmiddlefrontal','ctx-lh-caudalmiddlefrontal',
                       'ctx-rh-caudalanteriorcingulate','ctx-rh-rostralanteriorcingulate',
                       'ctx-lh-rostralmiddlefrontal','ctx-rh-caudalmiddlefrontal',
                       'ctx-lh-superiorparietal','ctx-rh-isthmuscingulate',
                       'ctx-lh-rostralanteriorcingulate','ctx-rh-parsopercularis',
                       'ctx-rh-superiorparietal','ctx-lh-parsopercularis',
                       'ctx-rh-medialorbitofrontal','ctx-lh-isthmuscingulate',
                       'ctx-lh-supramarginal','ctx-lh-inferiorparietal','ctx-rh-supramarginal',
                       'ctx-lh-superiorfrontal','ctx-rh-superiorfrontal','ctx-rh-middletemporal',
                       'ctx-lh-middletemporal','ctx-rh-inferiorparietal','ctx-rh-superiortemporal',
                       'ctx-lh-posteriorcingulate','ctx-lh-precuneus','ctx-lh-medialorbitofrontal',
                       'ctx-lh-superiortemporal','ctx-rh-posteriorcingulate','ctx-rh-precuneus']
    ordering = {x:i for i,x in enumerate(summary_regions)}
    rank_by = summary_regions # could take subset of cortical summary regions
    subjects = GROUPS['increasing_low']['N']
    df = df[df['rid'].isin(subjects)]

    baseline_keys = ["%s_bl" % _ for _ in rank_by]
    change_keys = ["%s_change" % _ for _ in summary_regions]
    df['rank'] = df[baseline_keys].mean(axis=1)

    keep_keys = ['rid', 'rank'] + change_keys
    df = df[keep_keys]
    df_long = pd.melt(df,id_vars=['rank'],value_vars=change_keys)

    # sort change
    df_long['variable'] = [_.replace('_change','') for _ in df_long['variable']]
    df_long['variable'] = ['%s_%s' % (str(ordering[_]).zfill(2),_) for _ in df_long['variable']]

    print ggplot(aes(x='variable',y='rank'),data=df_long)+geom_tile(aes(fill='value'))+theme(axis_text_x=element_text(angle=270,size=8), axis_text_y=element_text(size=6))

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


def valuesByFreesurferRegion(values_by_group, index_lookup):
    by_freesurfer_region = {}
    for group, val in values_by_group.iteritems():
        fs_indices = index_lookup[group]
        for fsi in fs_indices:
            by_freesurfer_region[fsi] = val
    return by_freesurfer_region

def regionEffectSizesBetweenGroups(group_prefix, group_one, group_two, data, index_lookup):
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
    line_data = defaultdict(dict)
    for k,v in group_effects.iteritems():
        for eff_k, eff_v in v.iteritems():
            line_data[k]['%s_uptake_effect_%s' % (group_prefix,eff_k)] = eff_v
    df = pd.DataFrame(dict(line_data)).T
    df.index.name = 'Region'
    return df

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
            sorted_uptakes_post.append(sorted(data_post[rid].iteritems(), key=lambda x: x[1][0], reverse=True))
            for k,(v,yrs) in data_post[rid].iteritems():
                group_uptakes_post[k].append(float(v))
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

def regionalRanksToLines(subj_group, group_prefix, data_prior, data_post, index_lookup):
    line_data = defaultdict(dict)
    (uptakes_prior, ranks_prior, uptakes_post, ranks_post, group_effects) = regionalEffectSizes(subj_group, data_prior, data_post, index_lookup)
    for k,v in uptakes_prior.iteritems():
        line_data[k]['%s_%s' % (group_prefix,'uptakes_prior_mean')] = v[0]
        line_data[k]['%s_%s' % (group_prefix,'uptakes_prior_std')] = v[1]
    for k,v in ranks_prior.iteritems():
        line_data[k]['%s_%s' % (group_prefix,'ranks_prior_mean')] = v[0]
        line_data[k]['%s_%s' % (group_prefix,'ranks_prior_std')] = v[1]
    for k,v in uptakes_post.iteritems():
        line_data[k]['%s_%s' % (group_prefix,'uptakes_post_mean')] = v[0]
        line_data[k]['%s_%s' % (group_prefix,'uptakes_post_std')] = v[1]
    for k,v in ranks_post.iteritems():
        line_data[k]['%s_%s' % (group_prefix,'ranks_post_mean')] = v[0]
        line_data[k]['%s_%s' % (group_prefix,'ranks_post_std')] = v[1]
    for k,v in group_effects.iteritems():
        uptake_effect = v['uptake_effect']
        for eff_k, eff_v in uptake_effect.iteritems():
            line_data[k]['%s_uptake_effect_%s' % (group_prefix,eff_k)] = eff_v
        rank_effect =v['rank_effect']
        for eff_k, eff_v in rank_effect.iteritems():
            line_data[k]['%s_rank_effect_%s' % (group_prefix,eff_k)] = eff_v
    df = pd.DataFrame(dict(line_data)).T
    df.index.name = 'Region'
    return df

if __name__ == "__main__":
    lut_file = "../FreeSurferColorLUT.txt"
    bilat_translate = bilateralTranslations(lut_file)
    lut_table = importFreesurferLookup(lut_file)
    master_file = '../FDG_AV45_COGdata_11_10_15.csv'
    master_data = importMaster(master_file)
    diags = extractDiagnosesFromMasterData(master_data)

    # # make bathroom tile plot
    # bl_vs_change_json = '../pvcsummary_bl_vs_change.json'
    # generateBathroomTilePlot(bl_vs_change_json)
    # sys.exit(1)

    # # summary pvc results
    # bl_file = '../output/Rousset_BL/raw_summary_output_BL.mat'
    # scan2_file = '../output/Rousset_Scan2/raw_summary_output_Scan2.mat'
    # scan3_file = '../output/Rousset_Scan3/raw_summary_output_Scan3.mat'
    # data_bl, data_scan2, data_scan3, index_lookup = parseRawRousset(bl_file, scan2_file, scan3_file)
    # data_bl, data_scan2, data_scan3, index_lookup = removeBlacklistedGroups(data_bl, data_scan2, data_scan3, index_lookup, suvr=True)

    # # summary tp-specific results
    # av45_file = "../output/UCBERKELEYAV45_09_25_15_extra.csv"
    # registry_file = "../docs/ADNI/REGISTRY.csv"
    # data_bl, data_scan2, data_scan3, index_lookup = parseRawAV45Output(av45_file, registry_file, lut_file)
    # data_bl, data_scan2, data_scan3, index_lookup = removeBlacklistedGroups(data_bl, data_scan2, data_scan3, index_lookup, suvr=True)

    # allregion pvc results
    bl_file = '../output/Rousset_BL/raw_allregions_output_BL.mat'
    scan2_file = '../output/Rousset_Scan2/raw_allregions_output_Scan2.mat'
    scan3_file = '../output/Rousset_Scan3/raw_allregions_output_Scan3.mat'
    data_bl, data_scan2, data_scan3, index_lookup = parseRawRousset(bl_file, scan2_file, scan3_file, translations=bilat_translate)
    data_bl, data_scan2, data_scan3, index_lookup = removeBlacklistedGroups(data_bl, data_scan2, data_scan3, index_lookup, suvr=False)

    # # allregion nontp results
    # av45_file = "../output/UCBERKELEYAV45_11_09_15_allregions_nontp.csv"
    # registry_file = "../docs/ADNI/REGISTRY.csv"
    # data_bl, data_scan2, data_scan3, index_lookup = parseAllRegionOutput(av45_file, lut_file)
    # data_bl, data_scan2, data_scan3, index_lookup = removeBlacklistedGroups(data_bl, data_scan2, data_scan3, index_lookup, suvr=False)

    all_data = []
    for rid, bl_row in data_bl.iteritems():
        master_row = master_data[rid]
        diag_1 = master_row['Diag@AV45_long']
        diag_2 = master_row['Diag@AV45_2_long']
        diag_3 = master_row['Diag@AV45_3_long']
        # bl row
        bl_datarow = {'rid': rid,
                      'diag': diag_1,
                      'timepoint': 'BL'}
        bl_datarow.update({k:bl_row[k] for k in index_lookup})
        all_data.append(bl_datarow)
        if rid in data_scan2:
            scan2_row = data_scan2[rid]
            scan2_datarow = {'rid': rid,
                             'diag': diag_2,
                             'timepoint': 'SCAN2'}
            scan2_datarow.update({k:scan2_row[k] for k in index_lookup})
            all_data.append(scan2_datarow)
        if rid in data_scan3:
            scan3_row = data_scan3[rid]
            scan3_datarow = {'rid': rid,
                             'diag': diag_3,
                             'timepoint': 'SCAN3'}
            scan3_datarow.update({k:scan3_row[k] for k in index_lookup})
            all_data.append(scan3_datarow)
    df = pd.DataFrame(all_data)
    df.to_csv('../datasets/pvc_allregions_uptake_bilateral.csv', index=False)

    data_prior = data_bl
    data_post = {}
    for rid, datarow in data_scan2.iteritems():
        yrs = float(master_data[rid]['AV45_1_2_Diff'])
        new_diag = master_data[rid]['Diag@AV45_2_long']
        withyrs = {k:(v,yrs) for k,v in datarow.iteritems()}
        data_post[rid] = (withyrs, new_diag)
    for rid, datarow in data_scan3.iteritems():
        yrs = float(master_data[rid]['AV45_1_3_Diff'])
        new_diag = master_data[rid]['Diag@AV45_3_long']
        withyrs = {k:(v,yrs) for k,v in datarow.iteritems()}
        data_post[rid] = (withyrs, new_diag)

    all_data = []
    for rid, prior_data in data_prior.iteritems():
        datarow = {'rid': rid,
                   'diag_prior': diags.get(rid,'')}
        prior_regions = {'%s_prior' % k.upper().replace('-','_') : prior_data[k] for k in index_lookup}
        datarow.update(prior_regions)
        if rid in data_post:
            post_data, diag_post = data_post[rid]
            datarow['diag_post'] = diag_post
            for k in index_lookup:
                k_f = k.upper().replace('-','_')
                post_val, yrs = post_data[k]
                datarow['%s_post' % k_f] = post_val
                datarow['yrs'] = yrs
        else:
            datarow['diag_post'] = ''
            datarow['yrs'] = ''
            post_regions = {'%s_post' % k.upper().replace('-','_') : '' for k in index_lookup}
            datarow.update(post_regions)
        all_data.append(datarow)
    df = pd.DataFrame(all_data)
    df.to_csv('../datasets/pvc_allregions_uptake_change_bilateral.csv', index=False)
    sys.exit(1)

