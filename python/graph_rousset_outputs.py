import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

from utils import *



def parseResult(arg):
	while type(arg) != np.ndarray or len(arg) <= 1:
		arg = arg[0]
	to_return = []
	for i, a in enumerate(arg):
		fields = a.dtype.names
		data = a[0][0]
		extracted = {}
		for k,v in zip(fields,data):
			while isinstance(v, np.ndarray) and len(v) <= 1:
				v = v[0]
			extracted[k] = v
		to_return.append(extracted)
	return to_return



if __name__ == "__main__":
	normal = [4222,4580,4555,4441,4151,4277]
	ad = [4172,4924,4755,4280,5187,5146]

	rousset_mat = '../rousset_outputs.mat'
	data = loadMATFile(rousset_mat)
	result_list = [parseResult(_) for _ in data['result_list'][0]]
	subj_list = data['subj_list'][0]
	groups = list(data['groups'][0])

	# sort by temporal_suvr_nonpvc
	sorted_data = sorted(zip(subj_list,result_list), key=lambda x: np.mean([_['temporal_suvr_nonpvc'] for _ in x[1]]))

	print [_[0] for _ in sorted_data]

	plt.figure(1)
	plt.subplot(251)
	plt.title('Temporal SUVR')
	for subj, subjdata in sorted_data:
		startval = np.mean([_['temporal_suvr_nonpvc'] for _ in subjdata])
		tograph = [_['temporal_suvr_pvc'] for _ in subjdata]
		color = 'b'
		if subj in normal:
			color = 'g'
		elif subj in ad:
			color = 'r'
		x = [0] + groups
		y = [startval] + tograph
		plt.plot(x,y, color)

	#plt.figure(2)
	plt.subplot(252)
	plt.title('Residuals')
	by_group = defaultdict(list)
	for subj, subjdata in sorted_data:
		tograph = [_['rmse_residual'] for _ in subjdata]
		for i,v in enumerate(tograph):
			by_group[i].append(v)
		color = 'b'
		if subj in normal:
			color = 'g'
		elif subj in ad:
			color = 'r'
		plt.plot(groups, tograph, color)
	avgs = [np.mean(by_group[i]) for i,g in enumerate(groups)]
	stds = [np.std(by_group[i]) for i,g in enumerate(groups)]
	print stds
	plt.errorbar(groups, avgs, stds)
	x1,x2,y1,y2 = plt.axis()
	plt.axis((0,5,y1,y2))

	#plt.figure(3)
	plt.subplot(253)
	plt.title('Whole Cereb Ref')
	for subj, subjdata in sorted_data:
		startval = np.mean([_['ref_nonpvc'] for _ in subjdata])
		tograph = [_['ref_pvc'] for _ in subjdata]
		color = 'b'
		if subj in normal:
			color = 'g'
		elif subj in ad:
			color = 'r'
		x = [0] + groups
		y = [startval] + tograph
		plt.plot(x,y, color)

	#plt.figure(4)
	plt.subplot(254)
	plt.title('Composite SUVR')
	for subj, subjdata in sorted_data:
		startval = np.mean([np.mean([_['parietal_suvr_nonpvc'],_['frontal_suvr_nonpvc'],_['temporal_suvr_nonpvc'],_['cingulate_suvr_nonpvc']])*_['ref_nonpvc'] for _ in subjdata])
		
		print subj
		print "parietal raw avg: %s" % np.mean([_['parietal_suvr_nonpvc']*_['ref_nonpvc'] for _ in subjdata])
		print "frontal raw avg: %s" % np.mean([_['frontal_suvr_nonpvc']*_['ref_nonpvc'] for _ in subjdata])
		print "temporal raw avg: %s" % np.mean([_['temporal_suvr_nonpvc']*_['ref_nonpvc'] for _ in subjdata])
		print "wcereb ref: %s" % np.mean([_['ref_nonpvc'] for _ in subjdata])
		print "cerebgm ref: %s" % np.mean([_['cerebGM_nonpvc'] for _ in subjdata])
		print "cerebwm ref: %s" % np.mean([_['cerebWM_nonpvc'] for _ in subjdata])
		print "brainstem ref: %s" % np.mean([_['brainstem_nonpvc'] for _ in subjdata])
		print "composite suvr avg: %s" % (startval)

		tograph = [np.mean([_['parietal_suvr_pvc'],_['frontal_suvr_pvc'],_['temporal_suvr_pvc'],_['cingulate_suvr_pvc']]) for _ in subjdata]
		color = 'b'
		if subj in normal:
			color = 'g'
		elif subj in ad:
			color = 'r'
		x = [0] + groups
		y = [startval] + tograph
		plt.plot(x,y, color)
	#plt.figure(5)
	plt.subplot(255)
	plt.title('Occipital SUVR')
	for subj, subjdata in sorted_data:
		startval = np.mean([_['occipital_suvr_nonpvc'] for _ in subjdata])
		tograph = [_['occipital_suvr_pvc'] for _ in subjdata]
		color = 'b'
		if subj in normal:
			color = 'g'
		elif subj in ad:
			color = 'r'
		x = [0] + groups
		y = [startval] + tograph
		plt.plot(x,y, color)

	#plt.figure(6)
	plt.subplot(256)
	plt.title('Cingulate SUVR')
	for subj, subjdata in sorted_data:
		startval = np.mean([_['cingulate_suvr_nonpvc'] for _ in subjdata])
		tograph = [_['cingulate_suvr_pvc'] for _ in subjdata]
		color = 'b'
		if subj in normal:
			color = 'g'
		elif subj in ad:
			color = 'r'
		x = [0] + groups
		y = [startval] + tograph
		plt.plot(x,y, color)

	#plt.figure(7)
	plt.subplot(257)
	plt.title('Frontal SUVR')
	for subj, subjdata in sorted_data:
		startval = np.mean([_['frontal_suvr_nonpvc'] for _ in subjdata])
		tograph = [_['frontal_suvr_pvc'] for _ in subjdata]
		color = 'b'
		if subj in normal:
			color = 'g'
		elif subj in ad:
			color = 'r'
		x = [0] + groups
		y = [startval] + tograph
		plt.plot(x,y, color)

	#plt.figure(8)
	plt.subplot(258)
	plt.title('Parietal SUVR')
	for subj, subjdata in sorted_data:
		startval = np.mean([_['parietal_suvr_nonpvc'] for _ in subjdata])
		tograph = [_['parietal_suvr_pvc'] for _ in subjdata]
		color = 'b'
		if subj in normal:
			color = 'g'
		elif subj in ad:
			color = 'r'
		x = [0] + groups
		y = [startval] + tograph
		plt.plot(x,y, color)

	#plt.figure(9)
	plt.subplot(259)
	plt.title('Brainstem')
	for subj, subjdata in sorted_data:
		startval = np.mean([_['brainstem_nonpvc'] for _ in subjdata])
		tograph = [_['brainstem_pvc'] for _ in subjdata]
		color = 'b'
		if subj in normal:
			color = 'g'
		elif subj in ad:
			color = 'r'
		x = [0] + groups
		y = [startval] + tograph
		plt.plot(x,y, color)

	
	#plt.figure(10)
	plt.subplot(2,5,10)
	plt.title('Basal Ganglia SUVR')
	for subj, subjdata in sorted_data:
		startval = np.mean([_['basal_ganglia_suvr_nonpvc'] for _ in subjdata])
		tograph = [_['basal_ganglia_suvr_pvc'] for _ in subjdata]
		color = 'b'
		if subj in normal:
			color = 'g'
		elif subj in ad:
			color = 'r'
		x = [0] + groups
		y = [startval] + tograph
		plt.plot(x,y, color)
	
	plt.show()

	'''
	print result_list
	print subj_list
	print groups
	'''

