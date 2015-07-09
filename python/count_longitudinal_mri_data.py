from collections import defaultdict
from datetime import datetime, timedelta
from tabulate import tabulate

from utils import *


def checkAvailablePointsPerSubject(pet_data, bsi_data, longfree_data, tbm_data, mri_data, master_data):
    '''
    Aggregate into one dictionary
    remove datapoints from before the first av45 scan
    '''
    points_by_subj = {}
    for subj,v in pet_data.iteritems():
        patient_pets = sorted(v)
        bl_av45 = patient_pets[0] - timedelta(days=90)
        subj_points = {}

        if subj in bsi_data:
            subj_bsi = sorted(bsi_data[subj], key=lambda x: x['EXAMDATE'])
        else:
            subj_bsi = []
        subj_points['bsi'] = [_['VISCODE2'] for _ in subj_bsi if _['EXAMDATE'] >= bl_av45]

        if subj in longfree_data:
            subj_long = sorted(longfree_data[subj], key=lambda x: x['EXAMDATE'])
        else:
            subj_long = []
        subj_points['long'] = [_['VISCODE2'] for _ in subj_long if _['EXAMDATE'] >= bl_av45]

        if subj in tbm_data:
            subj_tbm = sorted(tbm_data[subj], key=lambda x: x['EXAMDATE'])
        else:
            subj_tbm = []
        subj_points['tbm'] = [_['VISCODE2'] for _ in subj_tbm if _['EXAMDATE'] >= bl_av45]

        if subj in mri_data:
            subj_mri = sorted(mri_data[subj], key=lambda x: x['EXAMDATE'])
        else:
            subj_mri = []
        subj_points['mri'] = list(set([_['EXAMDATE'] for _ in subj_mri if _['EXAMDATE'] >= bl_av45]))
        
        #print "%s: %s -> %s" % (subj, len(subj_mri), len(subj_points['mri']))
        points_by_subj[subj] = subj_points

    # count
    tbm_counts = defaultdict(int)
    long_counts = defaultdict(int)
    bsi_counts = defaultdict(int)
    mri_counts = defaultdict(int)

    for k,v in points_by_subj.iteritems():
        #print k
        diag = master_data.get(k,{}).get('Init_Diagnosis','Unknown')
        for a,b in v.iteritems():
            #print '\t%s: %s' % (a,b)
            if a == 'long':
                if len(b) >= 2:
                    long_counts[(diag, len(b)-1)] += 1
                else:
                    long_counts[(diag, 0)] += 1
            elif a == 'tbm':
                tbm_counts[(diag, len(b))] += 1
            elif a == 'bsi':
                bsi_counts[(diag, len(b))] += 1
            elif a == 'mri':
                if len(b) >= 2:
                    mri_counts[(diag, len(b)-1)] += 1
                else:
                    mri_counts[(diag, 0)] += 1
    print tbm_counts.keys()

    '''
    for k,v in mri_counts.iteritems():
        print "%s, %s" % (k,v)
    '''

    print '\n\n'
    print "MRI"
    mri_graph = [['Init Diagnosis', 
                  '=0 followups', 
                  '=1 followups', 
                  '=2 followups', 
                  '=3 followups',
                  '=4 followups',
                  '=5 followups',
                  '=6 followups',
                  '=7 followups']]
    mri_graph.append(['CN',
                      mri_counts.get(('N',0),0),
                      mri_counts.get(('N',1),0),
                      mri_counts.get(('N',2),0),
                      mri_counts.get(('N',3),0),
                      mri_counts.get(('N',4),0),
                      mri_counts.get(('N',5),0),
                      mri_counts.get(('N',6),0),
                      mri_counts.get(('N',7),0)])
    mri_graph.append(['SMC',
                      mri_counts.get(('SMC',0),0),
                      mri_counts.get(('SMC',1),0),
                      mri_counts.get(('SMC',2),0),
                      mri_counts.get(('SMC',3),0),
                      mri_counts.get(('SMC',4),0),
                      mri_counts.get(('SMC',5),0),
                      mri_counts.get(('SMC',6),0),
                      mri_counts.get(('SMC',7),0)])
    mri_graph.append(['EMCI',
                      mri_counts.get(('EMCI',0),0),
                      mri_counts.get(('EMCI',1),0),
                      mri_counts.get(('EMCI',2),0),
                      mri_counts.get(('EMCI',3),0),
                      mri_counts.get(('EMCI',4),0),
                      mri_counts.get(('EMCI',5),0),
                      mri_counts.get(('EMCI',6),0),
                      mri_counts.get(('EMCI',7),0)])
    mri_graph.append(['LMCI',
                      mri_counts.get(('LMCI',0),0),
                      mri_counts.get(('LMCI',1),0),
                      mri_counts.get(('LMCI',2),0),
                      mri_counts.get(('LMCI',3),0),
                      mri_counts.get(('LMCI',4),0),
                      mri_counts.get(('LMCI',5),0),
                      mri_counts.get(('LMCI',6),0),
                      mri_counts.get(('LMCI',7),0)])
    mri_graph.append(['AD',
                      mri_counts.get(('AD',0),0),
                      mri_counts.get(('AD',1),0),
                      mri_counts.get(('AD',2),0),
                      mri_counts.get(('AD',3),0),
                      mri_counts.get(('AD',4),0),
                      mri_counts.get(('AD',5),0),
                      mri_counts.get(('AD',6),0),
                      mri_counts.get(('AD',7),0)])
    print tabulate(mri_graph) + '\n'
    print "TBMSyn"
    tbm_graph = [['Init Diagnosis', 
                  '=0 followups', 
                  '=1 followups', 
                  '=2 followups', 
                  '=3 followups',
                  '=4 followups',
                  '=5 followups',
                  '=6 followups',
                  '=7 followups']]
    tbm_graph.append(['CN',
                      tbm_counts.get(('N',0),0),
                      tbm_counts.get(('N',1),0),
                      tbm_counts.get(('N',2),0),
                      tbm_counts.get(('N',3),0),
                      tbm_counts.get(('N',4),0),
                      tbm_counts.get(('N',5),0),
                      tbm_counts.get(('N',6),0),
                      tbm_counts.get(('N',7),0)])
    tbm_graph.append(['SMC',
                      tbm_counts.get(('SMC',0),0),
                      tbm_counts.get(('SMC',1),0),
                      tbm_counts.get(('SMC',2),0),
                      tbm_counts.get(('SMC',3),0),
                      tbm_counts.get(('SMC',4),0),
                      tbm_counts.get(('SMC',5),0),
                      tbm_counts.get(('SMC',6),0),
                      tbm_counts.get(('SMC',7),0)])
    tbm_graph.append(['EMCI',
                      tbm_counts.get(('EMCI',0),0),
                      tbm_counts.get(('EMCI',1),0),
                      tbm_counts.get(('EMCI',2),0),
                      tbm_counts.get(('EMCI',3),0),
                      tbm_counts.get(('EMCI',4),0),
                      tbm_counts.get(('EMCI',5),0),
                      tbm_counts.get(('EMCI',6),0),
                      tbm_counts.get(('EMCI',7),0)])
    tbm_graph.append(['LMCI',
                      tbm_counts.get(('LMCI',0),0),
                      tbm_counts.get(('LMCI',1),0),
                      tbm_counts.get(('LMCI',2),0),
                      tbm_counts.get(('LMCI',3),0),
                      tbm_counts.get(('LMCI',4),0),
                      tbm_counts.get(('LMCI',5),0),
                      tbm_counts.get(('LMCI',6),0),
                      tbm_counts.get(('LMCI',7),0)])
    tbm_graph.append(['AD',
                      tbm_counts.get(('AD',0),0),
                      tbm_counts.get(('AD',1),0),
                      tbm_counts.get(('AD',2),0),
                      tbm_counts.get(('AD',3),0),
                      tbm_counts.get(('AD',4),0),
                      tbm_counts.get(('AD',5),0),
                      tbm_counts.get(('AD',6),0),
                      tbm_counts.get(('AD',7),0)])
    print tabulate(tbm_graph) + '\n'
    print "Longitudinal Freesurfer"
    long_graph = [['Init Diagnosis', 
                  '=0 followups', 
                  '=1 followups', 
                  '=2 followups', 
                  '=3 followups',
                  '=4 followups',
                  '=5 followups',
                  '=6 followups']]
    long_graph.append(['CN',
                      long_counts.get(('N',0),0),
                      long_counts.get(('N',1),0),
                      long_counts.get(('N',2),0),
                      long_counts.get(('N',3),0),
                      long_counts.get(('N',4),0),
                      long_counts.get(('N',5),0),
                      long_counts.get(('N',6),0)])
    long_graph.append(['SMC',
                      long_counts.get(('SMC',0),0),
                      long_counts.get(('SMC',1),0),
                      long_counts.get(('SMC',2),0),
                      long_counts.get(('SMC',3),0),
                      long_counts.get(('SMC',4),0),
                      long_counts.get(('SMC',5),0),
                      long_counts.get(('SMC',6),0)])
    long_graph.append(['EMCI',
                      long_counts.get(('EMCI',0),0),
                      long_counts.get(('EMCI',1),0),
                      long_counts.get(('EMCI',2),0),
                      long_counts.get(('EMCI',3),0),
                      long_counts.get(('EMCI',4),0),
                      long_counts.get(('EMCI',5),0),
                      long_counts.get(('EMCI',6),0)])
    long_graph.append(['LMCI',
                      long_counts.get(('LMCI',0),0),
                      long_counts.get(('LMCI',1),0),
                      long_counts.get(('LMCI',2),0),
                      long_counts.get(('LMCI',3),0),
                      long_counts.get(('LMCI',4),0),
                      long_counts.get(('LMCI',5),0),
                      long_counts.get(('LMCI',6),0)])
    long_graph.append(['AD',
                      long_counts.get(('AD',0),0),
                      long_counts.get(('AD',1),0),
                      long_counts.get(('AD',2),0),
                      long_counts.get(('AD',3),0),
                      long_counts.get(('AD',4),0),
                      long_counts.get(('AD',5),0),
                      long_counts.get(('AD',6),0)])
    print tabulate(long_graph) + '\n'
    print "BSI"
    bsi_graph = [['Init Diagnosis', 
                  '=0 followups', 
                  '=1 followups', 
                  '=2 followups', 
                  '=3 followups',
                  '=4 followups',
                  '=5 followups',
                  '=6 followups']]
    bsi_graph.append(['CN',
                      bsi_counts.get(('N',0),0),
                      bsi_counts.get(('N',1),0),
                      bsi_counts.get(('N',2),0),
                      bsi_counts.get(('N',3),0),
                      bsi_counts.get(('N',4),0),
                      bsi_counts.get(('N',5),0),
                      bsi_counts.get(('N',6),0)])
    bsi_graph.append(['SMC',
                      bsi_counts.get(('SMC',0),0),
                      bsi_counts.get(('SMC',1),0),
                      bsi_counts.get(('SMC',2),0),
                      bsi_counts.get(('SMC',3),0),
                      bsi_counts.get(('SMC',4),0),
                      bsi_counts.get(('SMC',5),0),
                      bsi_counts.get(('SMC',6),0)])
    bsi_graph.append(['EMCI',
                      bsi_counts.get(('EMCI',0),0),
                      bsi_counts.get(('EMCI',1),0),
                      bsi_counts.get(('EMCI',2),0),
                      bsi_counts.get(('EMCI',3),0),
                      bsi_counts.get(('EMCI',4),0),
                      bsi_counts.get(('EMCI',5),0),
                      bsi_counts.get(('EMCI',6),0)])
    bsi_graph.append(['LMCI',
                      bsi_counts.get(('LMCI',0),0),
                      bsi_counts.get(('LMCI',1),0),
                      bsi_counts.get(('LMCI',2),0),
                      bsi_counts.get(('LMCI',3),0),
                      bsi_counts.get(('LMCI',4),0),
                      bsi_counts.get(('LMCI',5),0),
                      bsi_counts.get(('LMCI',6),0)])
    bsi_graph.append(['AD',
                      bsi_counts.get(('AD',0),0),
                      bsi_counts.get(('AD',1),0),
                      bsi_counts.get(('AD',2),0),
                      bsi_counts.get(('AD',3),0),
                      bsi_counts.get(('AD',4),0),
                      bsi_counts.get(('AD',5),0),
                      bsi_counts.get(('AD',6),0)])
    print tabulate(bsi_graph) + '\n'

    print '\n\n'

    print "MRI"
    mri_graph_morethan = [['Init Diagnosis', 
                          '>=0 followups', 
                          '>=1 followups', 
                          '>=2 followups', 
                          '>=3 followups',
                          '>=4 followups',
                          '>=5 followups',
                          '>=6 followups',
                          '>=7 followups']]
    mri_graph_morethan.append(['CN',
                              sum(mri_graph[1][1:]),
                              sum(mri_graph[1][2:]),
                              sum(mri_graph[1][3:]),
                              sum(mri_graph[1][4:]),
                              sum(mri_graph[1][5:]),
                              sum(mri_graph[1][6:]),
                              sum(mri_graph[1][7:]),
                              sum(mri_graph[1][8:])])
    mri_graph_morethan.append(['SMC',
                              sum(mri_graph[2][1:]),
                              sum(mri_graph[2][2:]),
                              sum(mri_graph[2][3:]),
                              sum(mri_graph[2][4:]),
                              sum(mri_graph[2][5:]),
                              sum(mri_graph[2][6:]),
                              sum(mri_graph[2][7:]),
                              sum(mri_graph[2][8:])])
    mri_graph_morethan.append(['EMCI',
                              sum(mri_graph[3][1:]),
                              sum(mri_graph[3][2:]),
                              sum(mri_graph[3][3:]),
                              sum(mri_graph[3][4:]),
                              sum(mri_graph[3][5:]),
                              sum(mri_graph[3][6:]),
                              sum(mri_graph[3][7:]),
                              sum(mri_graph[3][8:])])
    mri_graph_morethan.append(['LMCI',
                              sum(mri_graph[4][1:]),
                              sum(mri_graph[4][2:]),
                              sum(mri_graph[4][3:]),
                              sum(mri_graph[4][4:]),
                              sum(mri_graph[4][5:]),
                              sum(mri_graph[4][6:]),
                              sum(mri_graph[4][7:]),
                              sum(mri_graph[4][8:])])
    mri_graph_morethan.append(['AD',
                              sum(mri_graph[5][1:]),
                              sum(mri_graph[5][2:]),
                              sum(mri_graph[5][3:]),
                              sum(mri_graph[5][4:]),
                              sum(mri_graph[5][5:]),
                              sum(mri_graph[5][6:]),
                              sum(mri_graph[5][7:]),
                              sum(mri_graph[5][8:])])
    print tabulate(mri_graph_morethan) + '\n'
    print "TBMSyn"
    tbm_graph_morethan = [['Init Diagnosis', 
                  '>=0 followups', 
                  '>=1 followups', 
                  '>=2 followups', 
                  '>=3 followups',
                  '>=4 followups',
                  '>=5 followups',
                  '>=6 followups',
                  '>=7 followups']]
    tbm_graph_morethan.append(['CN',
                              sum(tbm_graph[1][1:]),
                              sum(tbm_graph[1][2:]),
                              sum(tbm_graph[1][3:]),
                              sum(tbm_graph[1][4:]),
                              sum(tbm_graph[1][5:]),
                              sum(tbm_graph[1][6:]),
                              sum(tbm_graph[1][7:]),
                              sum(tbm_graph[1][8:])])
    tbm_graph_morethan.append(['SMC',
                              sum(tbm_graph[2][1:]),
                              sum(tbm_graph[2][2:]),
                              sum(tbm_graph[2][3:]),
                              sum(tbm_graph[2][4:]),
                              sum(tbm_graph[2][5:]),
                              sum(tbm_graph[2][6:]),
                              sum(tbm_graph[2][7:]),
                              sum(tbm_graph[2][8:])])
    tbm_graph_morethan.append(['EMCI',
                              sum(tbm_graph[3][1:]),
                              sum(tbm_graph[3][2:]),
                              sum(tbm_graph[3][3:]),
                              sum(tbm_graph[3][4:]),
                              sum(tbm_graph[3][5:]),
                              sum(tbm_graph[3][6:]),
                              sum(tbm_graph[3][7:]),
                              sum(tbm_graph[3][8:])])
    tbm_graph_morethan.append(['LMCI',
                              sum(tbm_graph[4][1:]),
                              sum(tbm_graph[4][2:]),
                              sum(tbm_graph[4][3:]),
                              sum(tbm_graph[4][4:]),
                              sum(tbm_graph[4][5:]),
                              sum(tbm_graph[4][6:]),
                              sum(tbm_graph[4][7:]),
                              sum(tbm_graph[4][8:])])
    tbm_graph_morethan.append(['AD',
                              sum(tbm_graph[5][1:]),
                              sum(tbm_graph[5][2:]),
                              sum(tbm_graph[5][3:]),
                              sum(tbm_graph[5][4:]),
                              sum(tbm_graph[5][5:]),
                              sum(tbm_graph[5][6:]),
                              sum(tbm_graph[5][7:]),
                              sum(tbm_graph[5][8:])])
    print tabulate(tbm_graph_morethan) + '\n'
    print "Longitudinal Freesurfer"
    long_graph_morethan = [['Init Diagnosis', 
                              '>=0 followups', 
                              '>=1 followups', 
                              '>=2 followups', 
                              '>=3 followups',
                              '>=4 followups',
                              '>=5 followups',
                              '>=6 followups']]
    long_graph_morethan.append(['CN',
                              sum(long_graph[1][1:]),
                              sum(long_graph[1][2:]),
                              sum(long_graph[1][3:]),
                              sum(long_graph[1][4:]),
                              sum(long_graph[1][5:]),
                              sum(long_graph[1][6:]),
                              sum(long_graph[1][7:])])
    long_graph_morethan.append(['SMC',
                              sum(long_graph[2][1:]),
                              sum(long_graph[2][2:]),
                              sum(long_graph[2][3:]),
                              sum(long_graph[2][4:]),
                              sum(long_graph[2][5:]),
                              sum(long_graph[2][6:]),
                              sum(long_graph[2][7:])])
    long_graph_morethan.append(['EMCI',
                              sum(long_graph[3][1:]),
                              sum(long_graph[3][2:]),
                              sum(long_graph[3][3:]),
                              sum(long_graph[3][4:]),
                              sum(long_graph[3][5:]),
                              sum(long_graph[3][6:]),
                              sum(long_graph[3][7:])])
    long_graph_morethan.append(['LMCI',
                              sum(long_graph[4][1:]),
                              sum(long_graph[4][2:]),
                              sum(long_graph[4][3:]),
                              sum(long_graph[4][4:]),
                              sum(long_graph[4][5:]),
                              sum(long_graph[4][6:]),
                              sum(long_graph[4][7:])])
    long_graph_morethan.append(['AD',
                              sum(long_graph[5][1:]),
                              sum(long_graph[5][2:]),
                              sum(long_graph[5][3:]),
                              sum(long_graph[5][4:]),
                              sum(long_graph[5][5:]),
                              sum(long_graph[5][6:]),
                              sum(long_graph[5][7:])])
    print tabulate(long_graph_morethan) + '\n'
    print "BSI"
    bsi_graph_morethan = [['Init Diagnosis', 
                  '>=0 followups', 
                  '>=1 followups', 
                  '>=2 followups', 
                  '>=3 followups',
                  '>=4 followups',
                  '>=5 followups',
                  '>=6 followups']]
    bsi_graph_morethan.append(['CN',
                              sum(bsi_graph[1][1:]),
                              sum(bsi_graph[1][2:]),
                              sum(bsi_graph[1][3:]),
                              sum(bsi_graph[1][4:]),
                              sum(bsi_graph[1][5:]),
                              sum(bsi_graph[1][6:]),
                              sum(bsi_graph[1][7:])])
    bsi_graph_morethan.append(['SMC',
                              sum(bsi_graph[2][1:]),
                              sum(bsi_graph[2][2:]),
                              sum(bsi_graph[2][3:]),
                              sum(bsi_graph[2][4:]),
                              sum(bsi_graph[2][5:]),
                              sum(bsi_graph[2][6:]),
                              sum(bsi_graph[2][7:])])
    bsi_graph_morethan.append(['EMCI',
                              sum(bsi_graph[3][1:]),
                              sum(bsi_graph[3][2:]),
                              sum(bsi_graph[3][3:]),
                              sum(bsi_graph[3][4:]),
                              sum(bsi_graph[3][5:]),
                              sum(bsi_graph[3][6:]),
                              sum(bsi_graph[3][7:])])
    bsi_graph_morethan.append(['LMCI',
                              sum(bsi_graph[4][1:]),
                              sum(bsi_graph[4][2:]),
                              sum(bsi_graph[4][3:]),
                              sum(bsi_graph[4][4:]),
                              sum(bsi_graph[4][5:]),
                              sum(bsi_graph[4][6:]),
                              sum(bsi_graph[4][7:])])
    bsi_graph_morethan.append(['AD',
                              sum(bsi_graph[5][1:]),
                              sum(bsi_graph[5][2:]),
                              sum(bsi_graph[5][3:]),
                              sum(bsi_graph[5][4:]),
                              sum(bsi_graph[5][5:]),
                              sum(bsi_graph[5][6:]),
                              sum(bsi_graph[5][7:])])
    print tabulate(bsi_graph_morethan) + '\n'



if __name__ == "__main__":

    include_failed = False

    # Input/output/lookup files
    master_file = "../FDG_AV45_COGdata_synced.csv"
    registry_file = "../docs/registry_clean.csv"
    pet_meta_file = "../docs/PET_META_LIST_edited.csv"
    mri_meta_file = "../docs/MPRAGEMETA.csv"
    # BSI file
    bsi_file = "../mr_docs/Fox/FOXLABBSI_04_30_15.csv"
    # long freesurfer file
    longfree_file = '../mr_docs/UCSF/longitudinal/UCSFFSL51Y1_08_01_14.csv'
    # TBMsyn file
    tbm_file = '../mr_docs/Mayo/MAYOADIRL_MRI_TBMSYN_05_07_15.csv'

    pet_data = importPetMETA(pet_meta_file)
    print 'PET patients: %s' % len(pet_data)
    bsi_data = importBSI(bsi_file, include_failed=include_failed)
    print 'BSI patients: %s' % len(bsi_data)
    longfree_data = importLongitudinalFreesurfer(longfree_file, include_failed=include_failed)
    print 'Longfree patients: %s' % len(longfree_data)
    tbm_data = importTBMSyn(tbm_file)
    print 'TBM patients: %s' % len(tbm_data)
    mri_data = importMRI(mri_meta_file)
    print 'MRI patients: %s' % len(mri_data)
    master_data = importMaster(master_file)

    avai_points = checkAvailablePointsPerSubject(pet_data, bsi_data, longfree_data, tbm_data, mri_data, master_data)





