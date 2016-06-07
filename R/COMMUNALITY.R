df_av45 = read.csv('nsfa/av45_factor_loadings.csv')
df_av1451 = read.csv('nsfa/av1451_factor_loadings.csv')
num_vars_av45 = dim(df_av45)[2]
num_vars_av1451 = dim(df_av1451)[2]

# square
df_av45.squared = df_av45^2
df_av1451.squared = df_av1451^2

# variance explained per factor
df_av45.varpercent = rowSums(df_av45.squared)/num_vars_av45
df_av1451.varpercent = rowSums(df_av1451.squared)/num_vars_av1451

# communality
df_av45.communality = colSums(df_av45.squared)
df_av1451.communality = colSums(df_av1451.squared)
df_av45.comm_total = sum(df_av45.communality)/num_vars_av45
df_av1451.comm_total = sum(df_av1451.communality)/num_vars_av1451

# plot
df_av45.comm_total
df_av1451.comm_total
plot(df_av45.communality)
plot(df_av1451.communality)
plot(df_av45.varpercent)
plot(df_av1451.varpercent)


