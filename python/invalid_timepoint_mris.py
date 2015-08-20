import csv
from datetime import datetime
from utils import *

'''
file_path = 'invalid_timepoint_mris.csv'
headers = ['RID','Timepoint','PET Date','PET/MRI Interval','Init Diagnosis', 'Timepoint Diagnosis', 
           'MRIMETA_closest_VC', 'MRIMETA_closest_STRENGTH', 'MRIMETA_closest_Interval', 'MRIMETA_closest_Date',
           'MRISEARCH_closest_VC', 'MRISEARCH_closest_DATE']
'''
file_path = 'invalid_timepoint_mris_meta.csv'
headers = ['RID','Timepoint','PET Date','PET/MRI Interval','Init Diagnosis', 'Timepoint Diagnosis', 
           'MRIMETA_closest_VC', 'MRIMETA_closest_STRENGTH', 'MRIMETA_closest_Interval', 'MRIMETA_closest_Date']


master_file = "../FDG_AV45_COGdata_08_06_15.csv"
master_dict = importMaster(master_file)

mri_meta_file = "../docs/MPRAGEMETA.csv"
mri_data = importMRI(mri_meta_file, magstrength_filter=None)

mri_search_file = "../docs/idaSearch_8_19_2015.csv"
search_data = importMRI(mri_search_file, magstrength_filter=None)

pet_meta_file = "../docs/PET_META_LIST_edited.csv"
pet_data = importPetMETA(pet_meta_file)


# create list of discrepancies from meta files
BL=[]
Scan2=[]
Scan3=[]
for subj, pet_dates in pet_data.iteritems():
    for i, pet_date in enumerate(pet_dates):
        search_rows = search_data.get(subj)
        if not search_rows:
            search_closest_interval = ''
        else:
            search_closest = sorted(search_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]
            search_closest_interval = abs(search_closest['EXAMDATE']-pet_date).days
            if search_closest_interval <= 130:
                continue
        if i == 0:
            BL.append((subj, search_closest_interval))
        elif i == 1:
            Scan2.append((subj, search_closest_interval))
        elif i == 2:
            Scan3.append((subj, search_closest_interval))

print BL
print Scan2
print Scan3

# hardcoded list of discrepancies
'''
BL=[('003-S-4524', 145), ('005-S-0448', 1029), ('005-S-0572', 549), ('031-S-0294', 729), ('033-S-0734', 364), ('033-S-1116', 879), ('051-S-1072', 211), ('051-S-1123', 133), ('067-S-0257', 345), ('098-S-0269', 360), ('100-S-0035', 629), ('127-S-0684', 1075), ('128-S-0138', 1750), ('130-S-0969', 358), ('029-S-0914', 455), ('016-S-5032', 135), ('073-S-5016', 161), ('073-S-5090', 160), ('003-S-5150', 201), ('082-S-5278', 147), ('057-S-5295', 138)]
Scan2=[('002-S-1268', 364), ('002-S-4746', 386), ('003-S-4136', 350), ('003-S-4555', 499), ('005-S-0572', 1184), ('006-S-0498', 399), ('006-S-1130', 352), ('006-S-4346', 353), ('006-S-4515', 377), ('009-S-1030', 359), ('012-S-4026', 370), ('013-S-2389', 288), ('016-S-2007', 501), ('018-S-2133', 153), ('018-S-4349', 217), ('018-S-4597', 407), ('019-S-4548', 381), ('019-S-4835', 569), ('022-S-0096', 364), ('022-S-4196', 385), ('022-S-4266', 273), ('029-S-4327', 364), ('031-S-4590', 372), ('032-S-0214', 426), ('052-S-0671', 387), ('057-S-1269', 373), ('067-S-0056', 329), ('067-S-0059', 329), ('072-S-2026', 330), ('072-S-4206', 351), ('072-S-4226', 377), ('072-S-4383', 396), ('072-S-4390', 408), ('072-S-4391', 389), ('072-S-4394', 394), ('072-S-4445', 385), ('072-S-4462', 370), ('072-S-4465', 393), ('072-S-4522', 396), ('072-S-4613', 427), ('082-S-4224', 372), ('082-S-4428', 477), ('098-S-0171', 338), ('098-S-0896', 409), ('098-S-2052', 346), ('098-S-4050', 654), ('098-S-4506', 236), ('100-S-4556', 334), ('109-S-4531', 329), ('116-S-4338', 370), ('116-S-4732', 327), ('126-S-4686', 371), ('127-S-0684', 1802), ('128-S-0138', 2541), ('128-S-0225', 139), ('130-S-0285', 384), ('137-S-4587', 374), ('141-S-2333', 712), ('141-S-4426', 464), ('141-S-4711', 378), ('941-S-4376', 505), ('012-S-4849', 414), ('072-S-4871', 577), ('116-S-4855', 386), ('127-S-0259', 385), ('003-S-4872', 531), ('012-S-4987', 386), ('027-S-4966', 472), ('072-S-4941', 583), ('127-S-5132', 364)]
Scan3=[('002-S-1280', 411), ('023-S-2068', 377), ('033-S-1016', 365), ('033-S-1098', 366), ('041-S-4041', 730), ('052-S-1352', 409), ('072-S-2072', 751), ('072-S-2093', 783), ('072-S-2116', 379), ('072-S-2164', 827), ('073-S-2191', 345), ('073-S-2264', 356), ('123-S-1300', 377), ('123-S-2363', 399), ('127-S-0260', 355), ('127-S-0684', 2551), ('128-S-0135', 763), ('128-S-0545', 398)]
'''


rows = []
for subj,days in BL:
    try:
        rid = int(subj.split('-')[-1])
    except Exception as e:
        rid = int(subj)
    master_row = master_dict[rid]
    pet_date = datetime.strptime(master_row['AV45_Date'], "%m/%d/%y")
    
    search_rows = search_data.get(rid)
    if not search_rows:
        search_closest = {}
    else:
        search_closest = sorted(search_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]
    
    meta_rows = mri_data.get(rid)
    if not meta_rows:
        meta_closest = {}
        meta_interval = ''
    else:
        meta_closest = sorted(meta_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]
        meta_interval = abs(meta_closest['EXAMDATE']-pet_date).days

    data = {'RID': rid,
            'PET/MRI Interval': days,
            'PET Date': pet_date,
            'Timepoint': 'BL',
            'Init Diagnosis': master_row['Init_Diagnosis'],
            'Timepoint Diagnosis': master_row['Diag@AV45_long'],
            'MRIMETA_closest_VC': meta_closest.get('vc'),
            'MRIMETA_closest_STRENGTH': meta_closest.get('strength'), 
            'MRIMETA_closest_Interval': meta_interval,
            'MRIMETA_closest_Date': meta_closest.get('EXAMDATE'),
            'MRISEARCH_closest_VC': search_closest.get('vc'), 
            'MRISEARCH_closest_DATE': search_closest.get('EXAMDATE')}
    rows.append(data)
for subj,days in Scan2:
    try:
        rid = int(subj.split('-')[-1])
    except Exception as e:
        rid = int(subj)
    master_row = master_dict[rid]
    pet_date = datetime.strptime(master_row['AV45_2_Date'], "%m/%d/%y")
    search_rows = search_data.get(rid)
    if not search_rows:
        search_closest = {}
    else:
        search_closest = sorted(search_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]

    meta_rows = mri_data.get(rid)
    if not meta_rows:
        meta_closest = {}
        meta_interval = ''
    else:
        meta_closest = sorted(meta_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]
        meta_interval = abs(meta_closest['EXAMDATE']-pet_date).days

    data = {'RID': rid,
            'PET/MRI Interval': days,
            'PET Date': pet_date,
            'Timepoint': 'Scan2',
            'Init Diagnosis': master_row['Init_Diagnosis'],
            'Timepoint Diagnosis': master_row['Diag@AV45_2_long'],
            'MRIMETA_closest_VC': meta_closest.get('vc'),
            'MRIMETA_closest_STRENGTH': meta_closest.get('strength'), 
            'MRIMETA_closest_Interval': meta_interval,
            'MRIMETA_closest_Date': meta_closest.get('EXAMDATE'),
            'MRISEARCH_closest_VC': search_closest.get('vc'), 
            'MRISEARCH_closest_DATE': search_closest.get('EXAMDATE')}
    rows.append(data)
for subj,days in Scan3:
    try:
        rid = int(subj.split('-')[-1])
    except Exception as e:
        rid = int(subj)
    master_row = master_dict[rid]
    pet_date = datetime.strptime(master_row['AV45_3_Date'], "%m/%d/%y")

    search_rows = search_data.get(rid)
    if not search_rows:
        search_closest = {}
    else:
        search_closest = sorted(search_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]

    meta_rows = mri_data.get(rid)
    if not meta_rows:
        meta_closest = {}
        meta_interval = ''
    else:
        meta_closest = sorted(meta_rows, key=lambda x: abs(x['EXAMDATE']-pet_date).days)[0]
        meta_interval = abs(meta_closest['EXAMDATE']-pet_date).days

    data = {'RID': rid,
            'PET/MRI Interval': days,
            'PET Date': pet_date,
            'Timepoint': 'Scan3',
            'Init Diagnosis': master_row['Init_Diagnosis'],
            'Timepoint Diagnosis': master_row['Diag@AV45_3_long'],
            'MRIMETA_closest_VC': meta_closest.get('vc'),
            'MRIMETA_closest_STRENGTH': meta_closest.get('strength'), 
            'MRIMETA_closest_Interval': meta_interval,
            'MRIMETA_closest_Date': meta_closest.get('EXAMDATE'),
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


