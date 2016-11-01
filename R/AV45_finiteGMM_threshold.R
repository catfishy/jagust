library(mixtools)

df = read.csv('~/Google Drive/ADNI_shared/Andy Documentation/FDG_AV45_COG_data/FDG_AV45_COGdata_10_28_16.csv',skip=1)
df = df[complete.cases(df$AV45_NONTP_1_BigRef),]
df$AV45_NONTP_1_BigRef_BIN.82 = (df$AV45_NONTP_1_BigRef >= 0.82)

# only take negative normals at baseline
# df = df[(df$Diag.AV45=='N' & df$AV45_NONTP_1_BigRef_BIN.82==0),]

# only take negatives at baseline
df = df[(df$AV45_NONTP_1_BigRef_BIN.82==0),]

# extract values
bigref.vals = df$AV45_NONTP_1_BigRef
twopt.slope.vals = df$AV45_NONTP_BigRef_Slope_2pts
twopt.slope.vals = twopt.slope.vals[!is.na(twopt.slope.vals)]
twopt.slope.vals = twopt.slope.vals[twopt.slope.vals < 0.05]
threept.slope.vals = df$AV45_NONTP_BigRef_Slope_3pts
threept.slope.vals = threept.slope.vals[!is.na(threept.slope.vals)]
threept.slope.vals = threept.slope.vals[threept.slope.vals < 0.05]
allpt.slope.vals = df$AV45_NONTP_BigRef_Slope
allpt.slope.vals = allpt.slope.vals[!is.na(allpt.slope.vals)]
allpt.slope.vals = allpt.slope.vals[allpt.slope.vals < 0.05]
hist(bigref.vals,breaks=100)
hist(twopt.slope.vals,breaks=100)
hist(threept.slope.vals,breaks=100)
hist(allpt.slope.vals,breaks=100)

# gaussian mixture model for cross sectional threshold
mix.results = normalmixEM(bigref.vals)
plot(mix.results, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8)
mus = mix.results$mu
sigmas = mix.results$sigma
lambdas = mix.results$lambda
f <- function(x) dnorm(x, m=mus[1], sd=sigmas[1]) * lambdas[1] - dnorm(x, m=mus[2], sd=sigmas[2]) * lambdas[2]
root.results = uniroot(f, interval=mus)
root.results$root

# gaussian mixture model for 2 pts slope
mix.results = normalmixEM(twopt.slope.vals)
plot(mix.results, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8)
mus = mix.results$mu
sigmas = mix.results$sigma
lambdas = mix.results$lambda
f <- function(x) dnorm(x, m=mus[1], sd=sigmas[1]) * lambdas[1] - dnorm(x, m=mus[2], sd=sigmas[2]) * lambdas[2]
root.results = uniroot(f, interval=c(0,0.06))
root.results$root

# gaussian mixture model for three points slope
mix.results = normalmixEM(threept.slope.vals)
plot(mix.results, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8)
mus = mix.results$mu
sigmas = mix.results$sigma
lambdas = mix.results$lambda
f <- function(x) dnorm(x, m=mus[1], sd=sigmas[1]) * lambdas[1] - dnorm(x, m=mus[2], sd=sigmas[2]) * lambdas[2]
root.results = uniroot(f, interval=c(-0.001,0.06))
root.results$root

# gaussian mixture model for all points slope
mix.results = normalmixEM(allpt.slope.vals)
plot(mix.results, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8)
mus = mix.results$mu
sigmas = mix.results$sigma
lambdas = mix.results$lambda
f <- function(x) dnorm(x, m=mus[1], sd=sigmas[1]) * lambdas[1] - dnorm(x, m=mus[2], sd=sigmas[2]) * lambdas[2]
root.results = uniroot(f, interval=c(0,0.06))
root.results$root
