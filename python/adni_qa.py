'''
Compares ADNI MRIS and PETS
- MPRAGE META file versus idaSearch (preprocessed scans that are downloadable)
- PETMETA file versus idaSearch
'''
import pandas as pd


def compareMRILists(search_list, meta_list, conducted_list):
    conducted_df = pd.read_csv(conducted_list)
    conducted_df = conducted_df[conducted_df.MMCONDCT==1]

    meta_df = pd.read_csv(meta_list)
    adni1_meta_df = meta_df[meta_df.Sequence.str.contains('N3;')]
    adni2_meta_df = meta_df[meta_df.Sequence.str.contains('N3m')]

    search_df = pd.read_csv(search_list)
    adni1_search_df = search_df[search_df.Description.str.contains('N3;')]
    adni2_search_df = search_df[search_df.Description.str.contains('N3m')]

    meta_images = set(adni1_meta_df['ImageUID']) | set(adni2_meta_df['ImageUID'])
    search_images = set(adni1_search_df['Image ID']) | set(adni2_search_df['Image ID'])

    print "Processed"
    print "META - SEARCH: %s" % (len(meta_images-search_images))
    print "SEARCH - META: %s" % (len(search_images-meta_images))

    # find originals with no processed
    series_ids = meta_df.groupby('SeriesID')
    orig_only = set()
    for series_id, rows in meta_df.groupby('SeriesID'):
        processed = rows[rows['Orig/Proc'] == 'Processed']
        original = rows[rows['Orig/Proc'] == 'Original']
        if len(processed) == 0:
            orig_only.add(series_id)
    orig_only_df = meta_df[meta_df.SeriesID.isin(orig_only)]
    adni2_orig_only_df = orig_only_df[orig_only_df.Visit.str.contains('ADNI2')]
    adni1_orig_only_df = orig_only_df[orig_only_df.Visit.str.contains('ADNI1')]

    valid_sequences = pd.DataFrame()
    for seq in ['MPRAGE','MP RAGE','MP-RAGE']:
        valid = adni2_orig_only_df[adni2_orig_only_df.Sequence.str.contains(seq)]
        valid_sequences = pd.concat([valid_sequences,valid])
    return valid_sequences

def comparePETLists(search_list, meta_list):
    meta_df = pd.read_csv(meta_list)
    meta_df = meta_df[meta_df.PMCONDCT==1]
    meta_df.loc[:,'EXAMDATE'] = pd.to_datetime(meta_df.loc[:,'EXAMDATE'])
    meta_by_rid = {int(k):v for k,v in meta_df.groupby('RID')}

    search_df = pd.read_csv(search_list)
    search_df.loc[:,'Study Date'] = pd.to_datetime(search_df.loc[:,'Study Date'])
    search_by_rid = {int(k.split('_')[-1]):v for k,v in search_df.groupby('Subject ID')}

    missing_in_search = pd.DataFrame()
    for rid, rows in meta_by_rid.iteritems():
        if rid not in search_by_rid:
            print "RID MISSING IN SEARCH %s" % rid
            missing_in_search = pd.concat((missing_in_search, rows))
            continue
        search_rows = search_by_rid[rid]
        if len(search_rows) != len(rows):
            print "DIFFERENT NUM SCANS %s" % rid
            print "META: %s, SEARCH: %s" % (len(rows), len(search_rows))
            search_dates = search_rows['Study Date']
            missing_rows = rows[~rows.EXAMDATE.isin(search_dates)]
            missing_in_search = pd.concat((missing_in_search, missing_rows))

    return missing_in_search

def scansConducted(MRIMETA, MRI3META, AV45META):
    av45_meta_df = pd.read_csv(AV45META)
    av45_meta_df = av45_meta_df[av45_meta_df.PMCONDCT==1]
    av45_meta_df.loc[:,'EXAMDATE'] = pd.to_datetime(av45_meta_df.loc[:,'EXAMDATE'])
    av45_meta = {k: list(set(v['EXAMDATE'])) for k,v in av45_meta_df.groupby('RID')}

    mri_meta_df = pd.read_csv(MRIMETA)
    mri_meta_df = mri_meta_df[mri_meta_df.MMCONDCT==1]
    mri_meta_df.loc[:,'EXAMDATE'] = pd.to_datetime(mri_meta_df.loc[:,'EXAMDATE'])
    mri_meta = {k: list(v['EXAMDATE']) for k,v in mri_meta_df.groupby('RID')}

    mri3_meta_df = pd.read_csv(MRI3META)
    mri3_meta_df = mri3_meta_df[mri3_meta_df.MMCONDCT==1]
    mri3_meta_df.loc[:,'EXAMDATE'] = pd.to_datetime(mri3_meta_df.loc[:,'EXAMDATE'])
    mri3_meta = {k: list(v['EXAMDATE']) for k,v in mri3_meta_df.groupby('RID')}

    return (mri_meta, mri3_meta, av45_meta)

def scansQCED(MRIMETA, PETMETA):
    # for mprage
    mri_df = pd.read_csv(MRIMETA)
    mri_df = mri_df[mri_df['Orig/Proc']=='Original']
    mri_df.loc[:,'ScanDate'] = pd.to_datetime(mri_df.loc[:,'ScanDate'])
    mri_df['Sequence'] = mri_df.Sequence.str.lower()
    # mri_df = mri_df[~mri_df.Sequence.str.contains('fspgr')]
    # mri_df = mri_df[mri_df.Sequence.str.contains('rage')]
    '''
    valid_sequences = pd.DataFrame()
    for seq in ['fspgr', 'rage']:
        valids = mri_df[mri_df.Sequence.str.contains(seq)]
        valid_sequences = pd.concat((valid_sequences, valids))
    mri_df = valid_sequences
    '''
    mri15_df = mri_df[mri_df.MagStrength==1.5]
    mri3_df = pd.concat((mri_df[mri_df.MagStrength==3],mri_df[mri_df.MagStrength==2.9]))
    mri_by_rid = {int(k.split('_')[-1]):list(v['ScanDate']) for k,v in mri15_df.groupby('SubjectID')}
    mri3_by_rid = {int(k.split('_')[-1]):list(v['ScanDate']) for k,v in mri3_df.groupby('SubjectID')}

    av45_df = pd.read_csv(PETMETA)
    av45_df = av45_df[av45_df['Orig/Proc']=='Original']
    av45_df.loc[:,'ScanDate'] = pd.to_datetime(av45_df.loc[:,'Scan Date'])
    av45_df['Sequence'] = av45_df.Sequence.str.lower()
    av45_df = av45_df[av45_df.Sequence.str.contains('av45')]
    av45_by_rid = {int(k.split('_')[-1]):list(set(v['Scan Date'])) for k,v in av45_df.groupby('Subject')}

    return mri_by_rid, mri3_by_rid, av45_by_rid

def scansSearchable(MRIMETA, PETMETA):
    mri_df = pd.read_csv(MRIMETA)
    adni1_mri_df = mri_df[mri_df.Description.str.contains('N3;')]
    adni2_mri_df = mri_df[mri_df.Description.str.contains('N3m')]
    mri_df = pd.concat((adni1_mri_df, adni2_mri_df))
    mri_df.loc[:,'Study Date'] = pd.to_datetime(mri_df.loc[:,'Study Date'])
    mri_by_rid = {int(k.split('_')[-1]):list(v['Study Date']) for k,v in mri_df.groupby('Subject ID')}

    av45_df = pd.read_csv(PETMETA)
    av45_df.loc[:,'Study Date'] = pd.to_datetime(av45_df.loc[:,'Study Date'])
    av45_by_rid = {int(k.split('_')[-1]):list(set(v['Study Date'])) for k,v in av45_df.groupby('Subject ID')}

    return mri_by_rid, av45_by_rid

def valLength(inputdict):
    val_lengths = [len(v) for v in inputdict.values()]
    return sum(val_lengths)

if __name__ == "__main__":
    MRI_CONDUCTED = "../docs/ADNI/MRIMETA.csv"
    MRI3_CONDUCTED = "../docs/ADNI/MRI3META.csv"
    PET_CONDUCTED = "../docs/ADNI/AV45META.csv"

    MRI_QCED = "../docs/ADNI/MPRAGEMETA.csv"
    PET_QCED = "../docs/ADNI/PET_META_LIST.csv"

    MRI_SEARCH = "../docs/ADNI/idaSearch_1_04_2016.csv"
    PET_SEARCH = "../docs/ADNI/idaSearch_1_07_2016.csv"

    mri_cond, mri3_cond, av45_cond = scansConducted(MRI_CONDUCTED, MRI3_CONDUCTED, PET_CONDUCTED)
    mri_qced, mri3_qced, av45_qced = scansQCED(MRI_QCED, PET_QCED)
    mri_search, av45_search = scansSearchable(MRI_SEARCH, PET_SEARCH)

    qc_keys = list(set(mri_qced.keys()) | set(mri3_qced.keys()))
    for k in qc_keys:
        mri_qc = list(mri_qced.get(k,[]))
        mri3_qc = list(mri3_qced.get(k,[]))
        searched = list(mri_search.get(k,[]))
        if len(mri_qc) + len(mri3_qc) != len(searched):
            print k
            print "mri: %s" % mri_qc
            print "mri3: %s" % mri3_qc
            print "search: %s" % searched

    # print "AV45 COND: %s" % valLength(av45_cond)
    # print "AV45 QCED: %s" % valLength(av45_qced)
    # print "AV45 SEAR: %s" % valLength(av45_search)
    print "MRI COND: %s" % valLength(mri_cond)
    print "MRI3 COND: %s" % valLength(mri3_cond)
    print "MRI(all) COND: %s" % (valLength(mri_cond) + valLength(mri3_cond))
    print "MRI QCED: %s" % valLength(mri_qced)
    print "MRI3 QCED: %s" % valLength(mri3_qced)
    print "MRI(all) QCED: %s" % (valLength(mri_qced) + valLength(mri3_qced))
    print "MRI SEAR: %s" % valLength(mri_search)
