from collections import defaultdict
import pylab as P
import numpy as np
from utils import *




def main(master_file):
    data = importMaster(master_file)
    # split by diagnosis at first av45 scan: 'Diag@AV45_long'
    by_diag = defaultdict(dict)
    for subj, d in data.iteritems():
        if int(subj) < 2000:
            continue
        diag = d['Diag@AV45_long'].strip()
        wcereb = d['AV45_wcereb'].strip()
        if diag == '':
            continue
        if wcereb == '':
            continue
        by_diag[diag][subj] = float(wcereb)

    AD = dict(by_diag['AD'])
    N = dict(by_diag['N'])

    # sort by whole cereb value: AV45_wcereb
    sorted_AD = sorted(AD.items(), key=lambda x: x[1])
    sorted_AD_vals = [_[1] for _ in sorted_AD]
    sorted_N = sorted(N.items(), key=lambda x: x[1])
    sorted_N_vals = [_[1] for _ in sorted_N]


    # print histograms
    P.figure()
    n, bins, patches = P.hist([sorted_AD_vals], 
                              70, 
                              histtype='bar',
                              label=['AD'])
    print "AD"
    for _ in zip(bins,n):
        print "\t%s" % (_,)
    P.legend()
    P.figure()
    n, bins, patches = P.hist([sorted_N_vals], 
                              70, 
                              histtype='bar',
                              label=['N'])
    print "N"
    for _ in zip(bins,n):
        print "\t%s" % (_,)

    P.legend()
    #P.show()

    # grab high, low, and around threshold (1.11)
    N_low = [(a,_) for a,_ in sorted_N if _ >= 1.009 and _ <= 1.0257]
    N_mid = [(a,_) for a,_ in sorted_N if _ >= 1.1087 and _ <= 1.1253]
    N_high = [(a,_) for a,_ in sorted_N if _ >= 1.3577 and _ <= 1.391]
    AD_low = [(a,_) for a,_ in sorted_AD if _ >= 0.9921 and _ <= 1.0259]
    AD_mid = [(a,_) for a,_ in sorted_AD if _ >= 1.09 and _ <= 1.125]
    AD_high = [(a,_) for a,_ in sorted_AD if _ >= 1.38 and _ <= 1.42]

    print "N_LOW"
    for _ in N_low: print _
    print "N_MID"
    for _ in N_mid: print _
    print "N_HIGH"
    for _ in N_high: print _
    print "AD_LOW"
    for _ in AD_low: print _
    print "AD_MID"
    for _ in AD_mid: print _
    print "AD_HIGH"
    for _ in AD_high: print _

if __name__ == "__main__":
    master = "../FDG_AV45_COGdata_07_28_15.csv"
    main(master)