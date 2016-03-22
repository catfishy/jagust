library(lme4)
library(coefplot2)
library(ggplot2)
library(lmerTest)
library(pbkrtest)
library(multcomp)
library(contrast)
library(xtable)
library(sjPlot)
library(splines)
library(car)
library(stats)
library(gdata)

# Import data
df_av45 = read.csv('dpgmm_alpha14.36_bilateral_AV45_ALL_longdata_continuous_slope.csv')

# demean pattern variables + binary variables
bin_variables = c('Gender','CORTICAL_SUMMARY_POSITIVE','APOE2_BIN','APOE4_BIN')
for (i in names(df_av45)){
  if (startsWith(i,'X') | i %in% bin_variables){
    df_av45[,eval(i)] = (df_av45[,eval(i)] * 2) - 1.0
  }
}

# Convert to factors
to_factor = c('CORTICAL_SUMMARY_POSITIVE','RID','group','diag_prior','APOE4_BIN','APOE2_BIN','Gender')
for (i in to_factor){
  df_av45[,eval(i)] = factor(df_av45[,eval(i)])
}

# add indicator variables
df_av45$AD = factor(as.numeric(df_av45$diag_prior == 'AD'))
df_av45$EMCI = factor(as.numeric(df_av45$diag_prior == 'EMCI'))
df_av45$LMCI = factor(as.numeric(df_av45$diag_prior == 'LMCI'))
df_av45$N = factor(as.numeric(df_av45$diag_prior %in% c('N','SMC')))


fm_av45_cs = glm(AD ~ CORTICAL_SUMMARY_prior, family='binomial', data=df_av45)
fm_av45_patterns = glm(AD ~ 
                            X23 +
                            X18 +
                            X42 +
                            X12 +
                            X16 +
                            X6 +
                            X19 +
                            X7 +
                            X1 +
                            X0, family='binomial', data=df_av45)

cs_summary = summary(fm_av45_cs)
patterns_summary = summary(fm_av45_patterns)