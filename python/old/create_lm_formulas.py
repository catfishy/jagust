import pandas as pd
import sys

av45_pos_input_files = {'full': '/usr/local/jagust/glmnet_av45pos_coef_full.csv',
                        'regions': '/usr/local/jagust/glmnet_av45pos_coef_regions.csv',
                        'onlyregions': '/usr/local/jagust/glmnet_av45pos_coef_onlyregions.csv',
                        'braak': '/usr/local/jagust/glmnet_av45pos_coef_braak.csv',
                        'patterns': '/usr/local/jagust/glmnet_av45pos_coef_patterns.csv'}
av45_neg_input_files = {'full': '/usr/local/jagust/glmnet_av45neg_coef_full.csv',
                        'regions': '/usr/local/jagust/glmnet_av45neg_coef_regions.csv',
                        'onlyregions': '/usr/local/jagust/glmnet_av45neg_coef_onlyregions.csv',
                        'braak': '/usr/local/jagust/glmnet_av45neg_coef_braak.csv',
                        'patterns': '/usr/local/jagust/glmnet_av45neg_coef_patterns.csv'}
all_subj_input_files = {'full': '/usr/local/jagust/glmnet_coef_full.csv',
                        'regions': '/usr/local/jagust/glmnet_coef_regions.csv',
                        'onlyregions': '/usr/local/jagust/glmnet_coef_onlyregions.csv',
                        'braak': '/usr/local/jagust/glmnet_coef_braak.csv',
                        'patterns': '/usr/local/jagust/glmnet_coef_patterns.csv'}



def cleanup(varname):
    return varname.replace('Gender1','Gender').replace('APOE4_BIN1','APOE4_BIN')

def createLMFormula(target, coeff):
    included = list(coeff[coeff!=0].index)
    included = map(cleanup,included)
    if len(included) == 0:
        return
    formula = "%s ~ %s" % (target,' + '.join(included))
    return formula

formulas = []
for name, filepath in av45_pos_input_files.iteritems():
    df = pd.read_csv(filepath)
    df.set_index('Unnamed: 0', inplace=True)
    for target, row in df.iterrows():
        lmform = createLMFormula(target,row)
        full_name = "%s.%s" % (name,target)
        formulas.append('"%s.%s"="%s"' % (name,target,lmform))
print "AV45 Pos"
# print "av45_pos_lm_list = c()"
# for form in formulas:
#     print "av45_pos_lm_list = c(av45_pos_lm_list, %s)" % form
print "c(%s)" % (', '.join(formulas))



formulas = []
for name, filepath in all_subj_input_files.iteritems():
    df = pd.read_csv(filepath)
    df.set_index('Unnamed: 0', inplace=True)
    for target, row in df.iterrows():
        lmform = createLMFormula(target,row)
        if lmform is None:
            formulas.append('"%s"=""' % (target,))
        else:
            formulas.append('"%s.%s"="%s"' % (name,target,lmform))
print '\n\n'
print "All Subj"
# print "all_subj_lm_list = c()"
# for form in formulas:
#     print "all_subj_lm_list = c(all_subj_lm_list, %s)" % form
print "c(%s)" % (', '.join(formulas))

formulas = []
for name, filepath in av45_neg_input_files.iteritems():
    df = pd.read_csv(filepath)
    df.set_index('Unnamed: 0', inplace=True)
    for target, row in df.iterrows():
        lmform = createLMFormula(target,row)
        if lmform is None:
            formulas.append('"%s.%s"=""' % (name,target,))
        else:
            formulas.append('"%s.%s"="%s"' % (name,target,lmform))
print '\n\n'
print "AV45 Neg Subj"
# print "all_subj_lm_list = c()"
# for form in formulas:
#     print "all_subj_lm_list = c(all_subj_lm_list, %s)" % form
print "c(%s)" % (', '.join(formulas))
