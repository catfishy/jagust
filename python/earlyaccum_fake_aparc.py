from utils import *

lut_file = "../FreeSurferColorLUT.txt"

lut_table = importFreesurferLookup(lut_file)
index_lookup = importFreesurferLookup(lut_file,flip=True)

pattern = {}
pattern_left = {}
pattern_right = {}

for roi in SUMMARY:
    pattern[lut_table[roi]] = 1
    if roi < 2000:
        pattern_left[lut_table[roi]] = 1
    else:
        pattern_right[lut_table[roi]] = 1
for roi in AV45_EARLYACCUM_MORE:
    pattern[lut_table[roi]] = 2
    if roi < 2000:
        pattern_left[lut_table[roi]] = 2
    else:
        pattern_right[lut_table[roi]] = 2
for roi in AV45_EARLYACCUM:
    pattern[lut_table[roi]] = 3
    if roi < 2000:
        pattern_left[lut_table[roi]] = 3
    else:
        pattern_right[lut_table[roi]] = 3


saveFakeAparcInput("/usr/local/jagust/earlyaccum_cortsummary_aparc",pattern,index_lookup)
saveFakeAparcInput("/usr/local/jagust/earlyaccum_cortsummary_aparc_left",pattern_left,index_lookup)
saveFakeAparcInput("/usr/local/jagust/earlyaccum_cortsummary_aparc_right",pattern_right,index_lookup)
