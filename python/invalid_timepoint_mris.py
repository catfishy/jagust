import csv
from datetime import datetime
from utils import *

BL=[('003-S-4524', 145), ('005-S-0448', 1029), ('005-S-0572', 549), ('031-S-0294', 729), ('033-S-0734', 364), ('033-S-1116', 879), ('051-S-1072', 211), ('051-S-1123', 133), ('067-S-0257', 345), ('098-S-0269', 360), ('100-S-0035', 629), ('127-S-0684', 1075), ('128-S-0138', 1750), ('130-S-0969', 358), ('137-S-0668', 352), ('021-S-0337', 343), ('029-S-0914', 455), ('114-S-2392', 356), ('016-S-5032', 138), ('073-S-5016', 161), ('073-S-5090', 160), ('027-S-5079', 182), ('027-S-5083', 393), ('128-S-5123', 144), ('003-S-5150', 201), ('082-S-5278', 147), ('057-S-5295', 138)]
Scan2=[('002-S-1268', 364), ('002-S-4213', 358), ('002-S-4746', 386), ('003-S-0907', 429), ('003-S-1122', 356), ('003-S-4119', 337), ('003-S-4136', 607), ('003-S-4555', 499), ('005-S-0572', 1184), ('006-S-0498', 399), ('006-S-1130', 352), ('006-S-4153', 367), ('006-S-4346', 353), ('006-S-4515', 377), ('007-S-2394', 369), ('009-S-1030', 359), ('011-S-4105', 339), ('011-S-4222', 337), ('012-S-4026', 370), ('012-S-4128', 341), ('013-S-2389', 288), ('016-S-2007', 501), ('016-S-4121', 365), ('018-S-2133', 153), ('018-S-4349', 217), ('018-S-4597', 407), ('019-S-4548', 381), ('019-S-4835', 569), ('022-S-0096', 364), ('022-S-4196', 385), ('022-S-4266', 273), ('024-S-4084', 404), ('024-S-4158', 375), ('029-S-2395', 358), ('029-S-4327', 364), ('031-S-4149', 351), ('031-S-4590', 372), ('032-S-0214', 426), ('033-S-4176', 368), ('033-S-4177', 371), ('035-S-2074', 371), ('035-S-4082', 320), ('035-S-4114', 336), ('036-S-2378', 387), ('041-S-4138', 348), ('041-S-4143', 384), ('052-S-0671', 387), ('057-S-1269', 373), ('067-S-0056', 329), ('067-S-0059', 329), ('067-S-4184', 362), ('068-S-4061', 429), ('072-S-2026', 330), ('072-S-4103', 373), ('072-S-4206', 351), ('072-S-4226', 377), ('072-S-4383', 396), ('072-S-4390', 408), ('072-S-4391', 389), ('072-S-4394', 394), ('072-S-4445', 385), ('072-S-4462', 370), ('072-S-4465', 393), ('072-S-4522', 396), ('072-S-4613', 427), ('072-S-4769', 395), ('073-S-4762', 504), ('073-S-4777', 588), ('082-S-2307', 364), ('082-S-4090', 500), ('082-S-4208', 399), ('082-S-4224', 372), ('082-S-4428', 477), ('098-S-0171', 338), ('098-S-0896', 409), ('098-S-2052', 346), ('098-S-4050', 654), ('098-S-4506', 642), ('099-S-4076', 409), ('099-S-4086', 395), ('100-S-4556', 334), ('109-S-4531', 410), ('116-S-4167', 356), ('116-S-4195', 395), ('116-S-4338', 370), ('116-S-4732', 327), ('123-S-2363', 359), ('123-S-4780', 343), ('123-S-4806', 594), ('126-S-4686', 371), ('126-S-4712', 457), ('127-S-0684', 1802), ('127-S-4148', 348), ('127-S-4624', 403), ('127-S-4765', 365), ('127-S-4843', 358), ('128-S-0138', 2541), ('128-S-0225', 139), ('128-S-2036', 393), ('128-S-4586', 547), ('128-S-4832', 699), ('130-S-0285', 384), ('130-S-0969', 329), ('130-S-4817', 355), ('137-S-4587', 374), ('137-S-4631', 358), ('137-S-4632', 418), ('141-S-2333', 712), ('141-S-4426', 464), ('141-S-4711', 378), ('153-S-4133', 373), ('153-S-4159', 387), ('941-S-4376', 505), ('012-S-4849', 596), ('067-S-4782', 565), ('072-S-4871', 577), ('114-S-0166', 397), ('116-S-4855', 553), ('126-S-4891', 364), ('127-S-0259', 385), ('127-S-4844', 350), ('128-S-4842', 376), ('003-S-4872', 531), ('013-S-4917', 140), ('073-S-4795', 370), ('116-S-4898', 526), ('127-S-4928', 295), ('137-S-4815', 547), ('137-S-4816', 563), ('130-S-4925', 360), ('012-S-4987', 386), ('027-S-4966', 472), ('072-S-4941', 583), ('073-S-4986', 595), ('073-S-5023', 370), ('094-S-4858', 386), ('114-S-5047', 360), ('128-S-5066', 378), ('127-S-5132', 364)]
Scan3=[('002-S-1280', 411), ('023-S-2068', 377), ('033-S-1016', 365), ('033-S-1098', 366), ('041-S-4041', 730), ('052-S-1352', 409), ('067-S-0059', 390), ('067-S-2301', 364), ('067-S-2304', 357), ('072-S-2037', 404), ('072-S-2072', 751), ('072-S-2093', 783), ('072-S-2116', 379), ('072-S-2164', 827), ('073-S-2191', 345), ('073-S-2264', 356), ('082-S-2121', 427), ('094-S-2201', 379), ('094-S-2216', 363), ('094-S-2238', 411), ('098-S-2052', 384), ('098-S-2079', 365), ('099-S-2146', 333), ('109-S-4531', 750), ('123-S-0106', 392), ('123-S-1300', 750), ('123-S-2363', 399), ('127-S-0260', 355), ('127-S-0684', 2551), ('127-S-0925', 358), ('127-S-1427', 367), ('127-S-2213', 366), ('127-S-2234', 351), ('128-S-0135', 763), ('128-S-0545', 398), ('128-S-2045', 366), ('128-S-2123', 349), ('128-S-2130', 413), ('128-S-2220', 385), ('129-S-0778', 385), ('129-S-2332', 399), ('131-S-0123', 384), ('137-S-0668', 757), ('137-S-0722', 724), ('137-S-0800', 370), ('137-S-0972', 370), ('137-S-0994', 372), ('137-S-1414', 357)]

file_path = 'invalid_timepoint_mris.csv'
headers = ['RID','Timepoint','PET Date','PET/MRI Interval','Init Diagnosis', 'Timepoint Diagnosis', 
           'MRIMETA_closest_VC', 'MRIMETA_closest_STRENGTH', 'MRIMETA_closest_Interval', 'MRIMETA_closest_Date',
           'MRISEARCH_closest_VC', 'MRISEARCH_closest_DATE']

master_file = "../FDG_AV45_COGdata_08_06_15.csv"
master_dict = importMaster(master_file)

mri_meta_file = "../docs/MPRAGEMETA.csv"
mri_data = importMRI(mri_meta_file, magstrength_filter=None)

mri_search_file = "../docs/idaSearch_8_19_2015.csv"
search_data = importMRI(mri_search_file, magstrength_filter=None)

rows = []
for subj,days in BL:
    rid = int(subj.split('-')[-1])
    master_row = master_dict[rid]
    pet_date = datetime.strptime(master_row['AV45_Date'], "%m/%d/%y")
    meta_rows = mri_data[rid]
    search_rows = search_data.get(rid)
    if not search_rows:
        search_closest = {}
    else:
        search_closest = sorted(search_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]
    
    meta_closest = sorted(meta_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]
    data = {'RID': rid,
            'PET/MRI Interval': days,
            'PET Date': pet_date,
            'Timepoint': 'BL',
            'Init Diagnosis': master_row['Init_Diagnosis'],
            'Timepoint Diagnosis': master_row['Diag@AV45_long'],
            'MRIMETA_closest_VC': meta_closest['vc'],
            'MRIMETA_closest_STRENGTH': meta_closest['strength'], 
            'MRIMETA_closest_Interval': abs(meta_closest['EXAMDATE']-pet_date).days,
            'MRIMETA_closest_Date': meta_closest['EXAMDATE'],
            'MRISEARCH_closest_VC': search_closest.get('vc'), 
            'MRISEARCH_closest_DATE': search_closest.get('EXAMDATE')}
    rows.append(data)
for subj,days in Scan2:
    rid = int(subj.split('-')[-1])
    master_row = master_dict[rid]
    pet_date = datetime.strptime(master_row['AV45_2_Date'], "%m/%d/%y")
    meta_rows = mri_data[rid]
    search_rows = search_data.get(rid)
    if not search_rows:
        search_closest = {}
    else:
        search_closest = sorted(search_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]
    meta_closest = sorted(meta_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]
    data = {'RID': rid,
            'PET/MRI Interval': days,
            'PET Date': pet_date,
            'Timepoint': 'Scan2',
            'Init Diagnosis': master_row['Init_Diagnosis'],
            'Timepoint Diagnosis': master_row['Diag@AV45_2_long'],
            'MRIMETA_closest_VC': meta_closest['vc'],
            'MRIMETA_closest_STRENGTH': meta_closest['strength'], 
            'MRIMETA_closest_Interval': abs(meta_closest['EXAMDATE']-pet_date).days,
            'MRIMETA_closest_Date': meta_closest['EXAMDATE'],
            'MRISEARCH_closest_VC': search_closest.get('vc'), 
            'MRISEARCH_closest_DATE': search_closest.get('EXAMDATE')}
    rows.append(data)
for subj,days in Scan3:
    rid = int(subj.split('-')[-1])
    master_row = master_dict[rid]
    pet_date = datetime.strptime(master_row['AV45_3_Date'], "%m/%d/%y")
    meta_rows = mri_data[rid]
    search_rows = search_data.get(rid)
    if not search_rows:
        search_closest = {}
    else:
        search_closest = sorted(search_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]
    meta_closest = sorted(meta_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]
    data = {'RID': rid,
            'PET/MRI Interval': days,
            'PET Date': pet_date,
            'Timepoint': 'Scan3',
            'Init Diagnosis': master_row['Init_Diagnosis'],
            'Timepoint Diagnosis': master_row['Diag@AV45_3_long'],
            'MRIMETA_closest_VC': meta_closest['vc'],
            'MRIMETA_closest_STRENGTH': meta_closest['strength'], 
            'MRIMETA_closest_Interval': abs(meta_closest['EXAMDATE']-pet_date).days,
            'MRIMETA_closest_Date': meta_closest['EXAMDATE'],
            'MRISEARCH_closest_VC': search_closest.get('vc'), 
            'MRISEARCH_closest_DATE': search_closest.get('EXAMDATE')}
    rows.append(data)

rows = [convertToCSVDataType(_) for _ in rows]

# find rows that are our fault
to_dl = []
for r in rows:
    if r['PET/MRI Interval'] != r['MRIMETA_closest_Interval'] and r['MRISEARCH_closest_DATE'] == r['MRIMETA_closest_Date']:
        to_dl.append(str(r['RID']))

print ','.join(to_dl)

dumpCSV(file_path, headers, rows)


