import csv
from collections import defaultdict
import random

from master_parser import parseCSV


def groupByDiagnosis(inputfile):
	headers, lines = parseCSV(inputfile)
	diagnoses = defaultdict(list)
	for l in lines:
		florbetapir = l['DIAG_AV45']
		ad = l['WCERB_BIN']

		if florbetapir == '' or ad == '':
			continue

		ad = int(ad)
		rid = int(l['RID'])
		diagnoses[(florbetapir,ad)].append(rid)
	diagnoses = dict(diagnoses)
	for k in diagnoses:
		diagnoses[k] = list(set(diagnoses[k]))
	return diagnoses


def groupByKey(master_file, key):
	'''
	Blacklist by column Q (reads already done)
	'''
	headers, lines = parseCSV(master_file)
	by_value = defaultdict(list)
	for l in lines:
		rid_string = l['RID']
		if rid_string == '':
			continue
		rid = int(rid_string)
		col_val = l[key]
		by_value[col_val].append(rid)
	by_value = dict(by_value)
	for k in by_value:
		by_value[k] = list(set(by_value[k]))
	return by_value

if __name__ == '__main__':
	diagnosis_file = '../docs/FDGnormalized_voxelwise_alltimepoints.csv'
	master_file = '../FDG_AV45_COGdata_synced.csv'

	diagnoses = groupByDiagnosis(diagnosis_file)
	master_groups = groupByKey(master_file, 'Florbetapir_vis_read')
	blacklist = set(master_groups['1']) | set(master_groups['0'])

	# get AB-, AD+
	ab_neg_AD = diagnoses[('AD',0)]
	ab_neg_AD_filtered = list(set(ab_neg_AD) - blacklist)
	# get AB+, AD+
	ab_pos_AD = diagnoses[('AD', 1)]
	ab_pos_AD_filtered = list(set(ab_pos_AD) - blacklist)
	# get AB-, LMCI
	ab_neg_LMCI = diagnoses[('LMCI',0)]
	ab_neg_LMCI_filtered = list(set(ab_neg_LMCI) - blacklist)

	print "Blacklist: %s\n" % blacklist
	print "AB- AD's: %s\n" % len(ab_neg_AD_filtered)
	print "AB+ AD's: %s\n" % len(ab_pos_AD_filtered)
	print "AB- LMCI's: %s\n" % len(ab_neg_LMCI_filtered)

	# choose every AB- AD, every AB- LMCI, and randomly pick 25 AB+ AD
	tocopy = []
	tocopy += ab_neg_AD_filtered
	tocopy += ab_neg_LMCI_filtered
	tocopy += random.sample(ab_pos_AD_filtered,25)

	print tocopy
	print len(tocopy)