'''
Compares ADNI MRIS
- MPRAGE META file versus idaSearch (preprocessed scans that are downloadable)
'''

import pandas as pd


def compareMRILists(search_list, meta_list):
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


if __name__ == "__main__":
    IDA_LIST = "../docs/ADNI/idaSearch_1_04_2016.csv"
    META_LIST = "../docs/ADNI/MPRAGEMETA.csv"
    compareMRILists(IDA_LIST, META_LIST) 