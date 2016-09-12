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
library(bootstrap)
library(zoo)
library(scales)

# FUNCTIONS

coef.heatplot = function(lars_coef) {
  base_size = 30
  lars_coef[,'(Intercept)'] = NULL
  lars_coef.m = melt(as.matrix(lars_coef))
  p = ggplot(lars_coef.m, aes(Var2, Var1)) + 
    geom_tile(aes(fill = value), 
              colour='grey50') + 
    scale_fill_gradient2(low="red", 
                         high="blue", 
                         limits=c(-1,1), 
                         oob=squish,
                         breaks=c(-1,-0.5,0,0.5,1),
                         labels=c("<= -1",'-0.5','0','0.5','>= 1')) + 
    labs(x="",y="") + 
    theme_bw() +
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed(ratio=1) + 
    theme(panel.border=element_blank(),
          axis.ticks = element_blank(),
          legend.key.size = unit(2, "cm"),
          legend.title=element_blank(),
          legend.text = element_text(size=base_size),
          axis.text.y = element_text(size=base_size),
          axis.text.x = element_text(size=base_size, 
                                     angle = 300, 
                                     hjust = 0))
  p
}


create.lm.from.coefcsv = function(coef_csv) {
  df = read.csv(coef_csv,row.names='X')
  df$X.Intercept. = NULL
  formulas = c()
  for(i in 1:nrow(df)) {
    row = df[i,]
    target = rownames(row)
    chosen = colnames(row)[row != 0]
    if (length(chosen) == 0L) {
      formulas[[target]] = ""
    } else {
      args = paste(sapply(chosen,clean.varname),collapse=' + ')
      lmform = paste(target,'~',args)
      formulas[[target]] = lmform
    }
  }
  formulas
}

clean.varname = function(name) {
  name = gsub('APOE4_BIN1','APOE4_BIN',name)
  name = gsub('Gender1','Gender',name)
  name
}


get.lars.coeff = function(x,y,target){
  full.lars.model = lars(x,y,type='lar')
  full.lars.coef = coef(full.lars.model, s=which.min(summary(full.lars.model)$Cp), mode='step')
  coef.df = as.data.frame(as.matrix(full.lars.coef))
  colnames(coef.df) = c(target)
  t(coef.df)
}

get.glmnet.coeff = function(x,y,target,use_min=FALSE,family='gaussian'){
  num_pts = dim(x)[1]
  cvfit = glmnet::cv.glmnet(x, y, 
                            nfolds=num_pts,
                            family=family)
#   localminima = rollapply(as.zoo(cvfit$cvm), 3, function(x) which.min(x)==2)
#   minima.lambda = cvfit$lambda[localminima]
  min.lambda = cvfit$lambda.min
  se.lambda = cvfit$lambda.1se
  mid.lambda = (min.lambda + se.lambda)/2
  if (use_min) {
    lambda = min.lambda
  } else {
    lambda = mid.lambda
  }
  coef.df = as.data.frame(as.matrix(coef(cvfit,s=lambda)))
  colnames(coef.df) = c(target)
  t(coef.df)
}

theta.fit <- function(x,y){lsfit(x,y)}

theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef} 

r2.shrinkage = function(form, target, data){
  fit = lm(form,data)
  X = getxy(form,data,FALSE)
  y = as.numeric(data[,target])
  results <- crossval(X,y,theta.fit,theta.predict,ngroup=10)
  raw.r2 = cor(y, fit$fitted.values)**2 # raw R2 
  cv.r2 = cor(y,results$cv.fit)**2 # cross-validated R2
  c(raw.r2,cv.r2)
}

r2.shrinkage.mean = function(form, target, data){
  results = c()
  for (i in 1:100) {
    test_r2 = r2.shrinkage(form,target,data)[2]
    results = c(results,test_r2)
  }
  mean(results)
}



isFSColumn = function(i){
  if (startsWith(i,'LEFT') || startsWith(i,'RIGHT') || startsWith(i,'CTX')) return(TRUE) else return(FALSE)
}
isFSColumn = Vectorize((isFSColumn))

isSizeColumn = function(i){
  if (length(grep("SIZE",i))>0) return(TRUE) else return(FALSE)
}
isSizeColumn = Vectorize((isSizeColumn))

isPatternColumn = function(i){
  if (startsWith(i,'NSFA')) return(TRUE) else return(FALSE)
}
isPatternColumn = Vectorize(isPatternColumn)

isNaiveColumn = function(i){
  if (startsWith(i,'SCORE_NSFA')) return(TRUE) else return(FALSE)
}
isNaiveColumn = Vectorize(isNaiveColumn)

isScan2Column = function(i){
  if (startsWith(i,'SCAN2_NSFA')) return(TRUE) else return(FALSE)
}
isScan2Column = Vectorize(isScan2Column)

isScan3Column = function(i){
  if (startsWith(i,'SCAN3_NSFA')) return(TRUE) else return(FALSE)
}
isScan3Column = Vectorize(isScan3Column)

r2adj = function(r2, n, p) {
  1 - (1 - r2) * ((n - 1)/(n - p - 1))
}

likelihood.test = function(null_form, test_form, data) {
  liketest = anova(lm(null_form,data),lm(test_form,data))
  liketest$`Pr(>F)`[2]
}


lm.addvar = function(var.name) {
  paste('+',paste(var.name,'*','APOE4_BIN',sep=''))
}

lm.createformula = function(target, form_str) {
  paste(target,'~',form_str)
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

mlm.testvar = function(var.name) {
  #CORTICAL_SUMMARY_prior*APOE4_BIN + I(CORTICAL_SUMMARY_prior^2)*APOE4_BIN + 
  #CORTICAL_SUMMARY_prior*APOE4_BIN + positive_prior*APOE4_BIN
  base_str = paste(target,"~","APOE4_BIN + Age.AV45 + Gender + Edu..Yrs.")
  add_str = lm.addvar(var.name)
  form_base = as.formula(base_str)
  form = as.formula(paste(base_str,add_str))
  fm = vglm(as.formula(form), family=multinomial(refLevel=1), data=df_av45)
  fm_base = vglm(as.formula(form_base), family=multinomial(refLevel=1), data=df_av45)
  like = VGAM::lrtest(fm, fm_base)
  like.p = like@Body$`Pr(>Chisq)`[2]
}

mlm.cv = function(dataset, form, target) {
  k = 10
  folds = cvFolds(nrow(dataset), K=k)
  holdoutpred = rep(0,nrow(dataset))
  for (i in 1:k) {
    train = dataset[folds$subsets[folds$which != i],]
    validation = dataset[folds$subsets[folds$which == i],]
    newlm = vglm(form, family=multinomial(refLevel=1), data=train)
    newprobs = VGAM::predict(newlm, validation, type='response')
    if (is.null(colnames(newprobs))) {
      newpred = apply(newprobs,1,which.max)
    } else {
      newpred = colnames(newprobs)[apply(newprobs,1,which.max)]
    }
    holdoutpred[folds$subsets[folds$which ==i]] = newpred
  }
  responses = as.numeric(dataset[,eval(target)])
  # find accuracy
  ctab = xtabs(as.formula(paste('~',target,'+ holdoutpred')), data=dataset)
  sum(diag(ctab))/sum(ctab)
}

lme.testvar = function(var.name) {
  #diag_prior*time + CORTICAL_SUMMARY_prior*time + 
  base_str = paste('value',"~","APOE4_BIN*time + Age.AV45 + Gender + Edu..Yrs.")
  add_str = lme.addvar(var.name)
  random_str = '+ (1 + time | RID)'
  form_base = as.formula(paste(base_str,random_str))
  form = as.formula(paste(base_str,add_str,random_str))
  fm = lmer(form,data=df_long)
  fm_base = lmer(form_base,df_long)
  like = anova(fm_base,fm)
  like.p = like$`Pr(>Chisq)`[2]
}

lme.cv = function(dataset, form) {
  k = 20
  subjects = levels(as.factor(dataset$RID))
  folds = cvFolds(length(subjects), K=k)
  holdoutpred = rep(0,nrow(dataset))
  for (i in 1:k) {
    train_subjects = subjects[folds$subsets[folds$which != i]]
    validation_subjects = subjects[folds$subsets[folds$which == i]]
    train = dataset[dataset$RID %in% train_subjects,]
    validation_indices = dataset$RID %in% validation_subjects
    validation = dataset[validation_indices,]
    newlm = lmer(as.formula(form),df_long)
    newpred = predict(newlm, newdata=validation)
    holdoutpred[validation_indices] = newpred
  }
  responses = dataset[,'value']
  rmse(holdoutpred,round(responses))
}


# Training models

run.glmnet = function(form, dataset, metric) {
  set.seed(825)
  fitControl = trainControl(method = "boot632",
                            number = 25,
                            repeats = 2)
  grid = expand.grid(.alpha = seq(0,1,length=50),
                     .lambda = 10^seq(1.5,-5,length=100))
  model = train(as.formula(form),
                data = dataset,
                method='glmnet',
                trControl = fitControl,
                tuneGrid = grid,
                preProcess=c('center', 'scale'),
                metric=metric)
  coefs = predict.glmnet(model$finalModel,type='coefficients',s=model$bestTune$lambda)
  results = model$results
  best_alpha = subset(results, alpha == model$bestTune$alpha)
  best_metric = subset(results, alpha == model$bestTune$alpha & lambda == model$bestTune$lambda)
  model 
}

run.lasso = function(form, dataset, metric) {
  set.seed(825)
  fitControl = trainControl(method = "boot632",
                            number = 25,
                            repeats = 5)
  grid = expand.grid(.fraction = seq(0,1,length=1000))
  model = train(as.formula(form),
                data = dataset,
                method='lasso',
                trControl = fitControl,
                tuneGrid = grid,
                preProcess=c('center', 'scale'),
                metric=metric)
  #preProcess=c('center', 'scale'))
  coefs = predict.enet(model$finalModel, type='coefficients',s=model$bestTune$fraction, mode='fraction')
  results = model$results
  best_alpha = subset(results, alpha == model$bestTune$alpha)
  best_metric = subset(results, fraction == model$bestTune$fraction)
  model 
}

getxy = function(form, dataset, remove) {
  x = as.matrix(as.data.frame(model.matrix(as.formula(form),dataset))[,-1])
  if (remove) {
    nzv_cols = nearZeroVar(x)
    if (length(nzv_cols) > 0) {
      x = x[, -nzv_cols]
    }
    if (dim(x)[2] >1) {
      corr_cols = findCorrelation(cor(x),.95)
      if (length(corr_cols) > 0) {
        x = x[, -corr_cols]
      }
    }
  }
  colnames(x) = lapply(colnames(x), make.names)
  rownames(x) = NULL
  as.matrix(x)
}


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

