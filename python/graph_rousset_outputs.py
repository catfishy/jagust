import numpy as np
import matplotlib.pyplot as plt
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
	plt.title('Temporal SUVR Difference')
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

	plt.figure(2)
	plt.title('Residuals')
	for subj, subjdata in sorted_data:
		tograph = [_['rmse_residual'] for _ in subjdata]
		color = 'b'
		if subj in normal:
			color = 'g'
		elif subj in ad:
			color = 'r'
		plt.plot(groups, tograph, color)

	plt.figure(3)
	plt.title('Whole Cereb Ref SUVR Difference')
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

	plt.figure(4)
	plt.title('Composite SUVR Difference')
	for subj, subjdata in sorted_data:
		startval = np.mean([np.mean([_['parietal_suvr_nonpvc'],_['frontal_suvr_nonpvc'],_['temporal_suvr_nonpvc'],_['occipital_suvr_nonpvc']]) for _ in subjdata])
		tograph = [np.mean([_['parietal_suvr_pvc'],_['frontal_suvr_pvc'],_['temporal_suvr_pvc'],_['occipital_suvr_pvc']]) for _ in subjdata]
		color = 'b'
		if subj in normal:
			color = 'g'
		elif subj in ad:
			color = 'r'
		x = [0] + groups
		y = [startval] + tograph
		plt.plot(x,y, color)

	plt.figure(5)
	plt.title('Occipital SUVR Difference')
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


	plt.show()

	'''
	print result_list
	print subj_list
	print groups
	'''

