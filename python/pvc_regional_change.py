from utils import *
import matplotlib.pyplot as plt

def extractRegionalValues(data):
	wholecereb = float(data['wholecereb']['pvcval'])
	bigref = float(data['bigref']['pvcval'])
	cingulate_uptake = data['cingulate']['pvcval']/wholecereb
	cingulate_vol = data['cingulate']['size']
	parietal_uptake = data['parietal']['pvcval']/wholecereb
	parietal_vol = data['parietal']['size']
	temporal_uptake = data['temporal']['pvcval']/wholecereb
	temporal_vol = data['temporal']['size']
	frontal_uptake = data['frontal']['pvcval']/wholecereb
	frontal_vol = data['frontal']['size']
	uptakes = {'cingulate': cingulate_uptake,
			   'parietal': parietal_uptake,
			   'temporal': temporal_uptake,
			   'frontal': frontal_uptake}
	sizes = {'cingulate': cingulate_vol,
			   'parietal': parietal_vol,
			   'temporal': temporal_vol,
			   'frontal': frontal_vol}
	return uptakes, sizes

master_file = '../FDG_AV45_COGdata_09_25_15.csv'
master_data = importMaster(master_file)

diags = {}
for rid, row in master_data.iteritems():
    diag = row['Init_Diagnosis'].strip()
    diags[rid] = diag

# Rousset output files
rousset_matfile_bl_manual = '../output/Rousset_BL/rousset_outputs_manual.mat'
rousset_matfile_scan2_manual = '../output/Rousset_Scan2/rousset_outputs_manual.mat'
rousset_matfile_scan3_manual = '../output/Rousset_Scan3/rousset_outputs_manual.mat'
rousset_matfile_bl_agg = '../output/Rousset_BL/rousset_outputs_agg.mat'
rousset_matfile_scan2_agg = '../output/Rousset_Scan2/rousset_outputs_agg.mat'
rousset_matfile_scan3_agg = '../output/Rousset_Scan3/rousset_outputs_agg.mat'

group_bl, data_bl = importRoussetResults(rousset_matfile_bl_manual)
group_scan2, data_scan2 = importRoussetResults(rousset_matfile_scan2_manual)
group_scan3, data_scan3 = importRoussetResults(rousset_matfile_scan3_manual)

# find annualized
yr_diff = {}
vol_frontal_diff = {}
vol_parietal_diff = {}
vol_temporal_diff = {}
vol_cingulate_diff = {}
uptake_frontal_diff = {}
uptake_parietal_diff = {}
uptake_temporal_diff = {}
uptake_cingulate_diff = {}

for rid, data in data_scan2.iteritems():
	yrs = float(master_data[rid]['AV45_1_2_Diff (Yrs)'])
	data = data['group4']
	residuals = data['residuals']['raw']
	scan2_uptakes, scan2_sizes = extractRegionalValues(data)
	bl_uptakes, bl_sizes = extractRegionalValues(data_bl[rid]['group4'])
	vol_frontal_diff[rid] = (scan2_sizes['frontal'] - bl_sizes['frontal']) / yrs
	vol_parietal_diff[rid] = (scan2_sizes['parietal'] - bl_sizes['parietal']) / yrs
	vol_temporal_diff[rid] = (scan2_sizes['temporal'] - bl_sizes['temporal']) / yrs
	vol_cingulate_diff[rid] = (scan2_sizes['cingulate'] - bl_sizes['cingulate']) / yrs
	uptake_frontal_diff[rid] = (scan2_uptakes['frontal'] - bl_uptakes['frontal']) / yrs
	uptake_parietal_diff[rid] = (scan2_uptakes['parietal'] - bl_uptakes['parietal']) / yrs
	uptake_temporal_diff[rid] = (scan2_uptakes['temporal'] - bl_uptakes['temporal']) / yrs
	uptake_cingulate_diff[rid] = (scan2_uptakes['cingulate'] - bl_uptakes['cingulate']) / yrs

for rid, data in data_scan3.iteritems():
	yrs = float(master_data[rid]['AV45_1_3_Diff (yrs)'])
	data = data['group4']
	residuals = data['residuals']['raw']
	scan3_uptakes, scan3_sizes = extractRegionalValues(data)
	bl_uptakes, bl_sizes = extractRegionalValues(data_bl[rid]['group4'])
	vol_frontal_diff[rid] = (scan3_sizes['frontal'] - bl_sizes['frontal']) / yrs
	vol_parietal_diff[rid] = (scan3_sizes['parietal'] - bl_sizes['parietal']) / yrs
	vol_temporal_diff[rid] = (scan3_sizes['temporal'] - bl_sizes['temporal']) / yrs
	vol_cingulate_diff[rid] = (scan3_sizes['cingulate'] - bl_sizes['cingulate']) / yrs
	uptake_frontal_diff[rid] = (scan3_uptakes['frontal'] - bl_uptakes['frontal']) / yrs
	uptake_parietal_diff[rid] = (scan3_uptakes['parietal'] - bl_uptakes['parietal']) / yrs
	uptake_temporal_diff[rid] = (scan3_uptakes['temporal'] - bl_uptakes['temporal']) / yrs
	uptake_cingulate_diff[rid] = (scan3_uptakes['cingulate'] - bl_uptakes['cingulate']) / yrs

all_rids = list(set(data_scan2.keys() + data_scan3.keys()))
colors = []
for rid in all_rids:
	if diags[rid] == 'AD':
		c = 'r'
	elif diags[rid] == 'N':
		c = 'g'
	else:
		c = 'b'
	colors.append(c)

frontal_x = []
frontal_y = []
for rid in all_rids:
	frontal_x.append(vol_frontal_diff[rid])
	frontal_y.append(uptake_frontal_diff[rid])
parietal_x = []
parietal_y = []
for rid in all_rids:
	parietal_x.append(vol_parietal_diff[rid])
	parietal_y.append(uptake_parietal_diff[rid])
temporal_x = []
temporal_y = []
for rid in all_rids:
	temporal_x.append(vol_temporal_diff[rid])
	temporal_y.append(uptake_temporal_diff[rid])
cingulate_x = []
cingulate_y = []
for rid in all_rids:
	cingulate_x.append(vol_cingulate_diff[rid])
	cingulate_y.append(uptake_cingulate_diff[rid])

plt.figure(1)
plt.scatter(frontal_x, frontal_y, c=colors)
plt.xlabel('Annualized Frontal Size Change')
plt.ylabel('Annualized Frontal Uptake Change (post PVC)')
plt.figure(2)
plt.scatter(parietal_x, parietal_y, c=colors)
plt.xlabel('Annualized Parietal Size Change')
plt.ylabel('Annualized Parietal Uptake Change (post PVC)')
plt.figure(3)
plt.scatter(temporal_x, temporal_y, c=colors)
plt.xlabel('Annualized Temporal Size Change')
plt.ylabel('Annualized Temporal Uptake Change (post PVC)')
plt.figure(4)
plt.scatter(cingulate_x, cingulate_y, c=colors)
plt.xlabel('Annualized Cingulate Size Change')
plt.ylabel('Annualized Cingulate Uptake Change (post PVC)')

plt.show()
