library(lme4)
library(ggplot2)
library(pbkrtest)
library(contrast)
library(xtable)
library(car)
library(stats)
library(gdata)
library(psych)
library(reshape2)
library(LambertW)
library(nnet)
library(DAAG)
library(caret)
library(cvTools)
library(VGAM)
library(lmtest)
library(languageR)
library(stringr)

# FUNCTIONS

isPatternColumn = function(i){
  if (startsWith(i,pattern_prefix)) return(TRUE) else return(FALSE)
}
isPatternColumn = Vectorize(isPatternColumn)

lm.addvar = function(var.name) {
  paste('+',paste(var.name,'*','APOE4_BIN',sep=''))
}

lme.addvar = function(var.name) {
  paste('+',paste(var.name,'*','time',sep=''))
}

save.printout = function(output_file, obj) {
  sink(output_file); print(obj, correlation=TRUE); sink(file=NULL)
}

save.plot = function(output_file, plot_fn) {
  pdf(file=output_file);plot_fn();dev.off();
}


plot.model = function(model) {
  plot(model)
}

fm.relimp = function(model) {
  calc.relimp(model, type=c('lmg'), rela=TRUE)
}

rmse = function(m, o) {
  sqrt(mean((m-o)^2))
}

to.long = function(df, time_col_prefix, value_col_prefix) {
  # Keep relevant columns
  time_columns = Filter(function(i){startsWith(i,time_col_prefix)}, names(df))
  value_columns = Filter(function(i){startsWith(i,value_col_prefix)}, names(df))
  df = df[c(demog_columns,av45_columns,pattern_columns,time_columns,value_columns)]
  # Convert to long format
  df_time_wide = df[c(demog_columns,av45_columns,pattern_columns,time_columns)]
  colnames(df_time_wide) = gsub(time_col_prefix,'TP',names(df_time_wide))
  df_value_wide = df[c(demog_columns,av45_columns,pattern_columns,value_columns)]
  colnames(df_value_wide) = gsub(value_col_prefix,'TP',names(df_value_wide))
  df_time_long = melt(df_time_wide, 
                      id.vars=c(demog_columns,av45_columns,pattern_columns),
                      measure.vars=Filter(function(x){startsWith(x,'TP')},names(df_time_wide)),
                      variable.name='timepoint',
                      value.name='time')
  df_value_long = melt(df_value_wide,
                       id.vars=c(demog_columns,av45_columns,pattern_columns),
                       measure.vars=Filter(function(x){startsWith(x,'TP')},names(df_value_wide)),
                       variable.name='timepoint',
                       value.name='value')
  merge_on = c(demog_columns,av45_columns,pattern_columns,'timepoint')
  df_long = merge(df_time_long,df_value_long,merge_on)
  df_long[complete.cases(df_long[,names(df_long)]),]
}

# Training models

run.rfe = function(form, var.response, dataset, min_size) {
  x = as.data.frame(model.matrix(as.formula(form),dataset))[,-1]
  nzv_cols = nearZeroVar(x)
  if (length(nzv_cols) > 0) {
    x = x[, -nzv_cols]
  }
  corr_cols = findCorrelation(cor(x),.9)
  if (length(corr_cols) > 0) {
    x = x[, -corr_cols]
  }
  
  colnames(x) = lapply(colnames(x), make.names)
  rownames(x) = NULL
  y = as.numeric(dataset[,var.response])
  
  ctrl = rfeControl(functions = lmFuncs, 
                    method = "repeatedcv", 
                    number = 10,
                    repeats = 5,
                    rerank = TRUE,
                    verbose = FALSE)
  set.seed(1337)
  rfe.output = rfe(x, 
                   y, 
                   sizes = c(min_size:ncol(x)),
                   rfeControl = ctrl,
                   metric = 'Rsquared')
  rfe.output
}

