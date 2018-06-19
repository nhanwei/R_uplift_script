eval_lift = function(pred, gt, bin = 10)
{
  re = data.frame(pred = pred, gt = gt)
  re = re[order(-re$pred),]
  stat = NULL
  for (i in 1:bin)
  {
    if (i == bin)
    {
      tt =  dim(re)[1]
    }
    else
    {
      tt = round((i * dim(re)[1]) / bin)
    }
    #print(sum(re[1:tt,]$gt) / round((i * dim(re)[1]) / 10) / mean(re$gt))
    l_ind = round(((i-1) * dim(re)[1]) / bin)
    obs =   tt - l_ind + 1
    event = sum(re[l_ind:tt,]$gt)
    non_event = obs - event
    cum_obs = tt
    cum_event = sum(re[1:tt,]$gt)
    cum_non_event = cum_obs - sum(re[1:tt,]$gt)
    cum_event_perc = cum_event/sum(re$gt)
    cum_non_event_perc =cum_non_event/(dim(re)[1]-sum(re$gt))
    ks = cum_event_perc - cum_non_event_perc
    lift = sum(re[1:tt,]$gt) / round((i * dim(re)[1]) / bin) / mean(re$gt)
    temp = data.frame(
      index = i, obs = obs, events = event, Non_event = non_event, cum_obs = cum_obs,
      cum_event = cum_event, cum_non_event = cum_non_event, cum_event_perc  = cum_event_perc,
      cum_non_event_perc = cum_non_event_perc, ks = ks, lift = lift
    )
    if (is.null(stat)) {
      stat = temp
    }
    else{
      stat = rbind(stat, temp)
    }
  }
  stat
}


eval_gain_ratio = function(pred, gt, bin = 10)
{
  re = data.frame(pred = pred, gt = gt)
  re = re[order(-re$pred),]

  for (i in 1:10)
  {
    print(sum(re[1:round((i * dim(re)[1]) / 10),]$gt) / sum(re$gt))

  }
}

eval_gain_ratio_num = function(pred, gt, bin = 10)
{
  re = data.frame(pred = pred, gt = gt)
  re = re[order(-re$pred),]

  for (i in 1:10)
  {
    le = round(i * 500)
    nu = sum(re[1:round(i * 500),]$gt) / sum(re$gt)
    ca = sum(re[1:round(i * 500),]$gt)
    temp = data.frame(le = le, nu = nu, ca = ca)
    print(temp)
  }
}


get_all_missing_col = function(sampleData)
{
  dput(colnames(sampleData)[sapply(sampleData, function(x)all(is.na(x)))])
  colnames(sampleData)[sapply(sampleData, function(x)all(is.na(x)))]

}

get_all_equal_col = function(sampleData)
{
  dput(colnames(sampleData)[sapply(sampleData, function(x) length(unique(x))<2)])
  colnames(sampleData)[sapply(sampleData, function(x)length(unique(x))<2)]
}



gbm_plot = function(gbmModel, top_n)
{
  imp_fea_names = names(relative.influence(
    gbmModel, n.trees = gbmModel$n.trees, scale. = T, sort. = T
  )[1:top_n])
  for (cl in imp_fea_names)
  {
    i = which(gbmModel$var.names == cl)
    print(cl)
    temp = plot(
      gbmModel,
      i.var = i,
      n.trees = gbmModel$n.trees,
      continuous.resolution = 200,
      return.grid = TRUE,
      type = "response", xlim = c(0,2)
    )
    jpeg(paste0('graph/',cl,'.jpg'))
    if (!is.factor(temp[,names(temp)[1]]))
    {
      if (max(temp[,names(temp)[1]]) - min(temp[,names(temp)[1]]) > 1000)
      {
        plot(
          x = log(1 + temp[,names(temp)[1]], base = 10),y = temp$y,type = "l"
          ,xlab = paste0(names(temp)[1],'(log_scale)'), ylab = 'probability'
        )
      }else{
        plot(
          x = temp[,names(temp)[1]], base = 10,y = temp$y,type = "l"
          ,xlab = paste0(names(temp)[1]), ylab = 'probability'
        )
      }
    }else{
      plot(
        gbmModel,
        i.var = i,
        n.trees = gbmModel$n.trees,
        continuous.resolution = 200,
        return.grid = FALSE,
        type = "response", xlim = c(0,2)
      )

    }
    dev.off()
  }
}



replace_missing_all = function(data)
{
  for (i in 1:dim(data)[2])
  {
    print(paste0('replace the columns ',colnames(data)[i]))
    data = replace_missing(data, colnames(data)[i])
  }
  data
}

replace_missing = function(data, clnames)
{
  i = which(colnames(data) == clnames)
  if (is.factor(data[,i]))
  {
    data[, i] = as.character(data[, i])
    data[is.na(data[, i]), i] = 'NA'
    data[, i] = as.factor(data[, i])
  }
  else{
    data[is.na(data[, i]), i] = 0
  }
  data
}


numeric_to_factor = function(data,excl_cl)
{
  for (i in 1:dim(data)[2])
  {
    if (is.numeric(data[,i]))
    {
      if(!colnames(data)[i] %in% excl_cl)
      {
        print(colnames(data)[i])
        print(i)
        temp = data[,i]
        temp = rep('M', dim(data)[1])
        temp[data[,i]<=quantile(data[,i],0.25)] = paste0( "L(<=",quantile(data[,i],0.25),")")
        temp[data[,i]> quantile(data[,i],0.75)] = paste0( "H(>",quantile(data[,i],0.75),")")
        data[,i] = temp
      }
    }
  }
  data
}

feature_report2 = function(sampleData,var,target)
{
  library(reshape2)
  fm = as.formula(paste(var,'~',target))
  bivar_result = dcast(sampleData, fm)
  bivar_result$Total = bivar_result[,'0'] + bivar_result[,'1']
  #bivar_result$Event_rate = percent(bivar_result[,'1']/bivar_result$Total)
  bivar_result$Event_rate = round(bivar_result[,'1']/bivar_result$Total,2)
  bivar_result$group = var
  bivar_result
}

bin_value = function(sampleData, fea, bin = 5)
{

  library(Hmisc)
  for(f in fea)
  {
    if(is.numeric(sampleData[,f]))
    {
      print(paste0('processing column ', f))
      sampleData$temp = cut2(sampleData[,f], g = bin)
      sampleData$temp = factor(sampleData$temp, levels=c(levels(sampleData$temp), 'NA'))
      sampleData$temp[is.na(sampleData$temp)] = 'NA'
      colnames(sampleData)[which(colnames(sampleData) == 'temp')] = paste0(f,'_bin',bin)
    }
  }
  sampleData
}

fea_clean = function(sampleData, fea)
{
  t = unlist(sapply(sampleData[,fea], function(x) if(is.factor(x)) {length(levels(x))>32} ))
  fea = names(t)[t== F]
}


gbm_imp = function(sampleData, fea, iter = 1)
{
  g1 = gbm(event ~.,         # formula
           data = sampleData[, fea],                   # dataset
           distribution = "bernoulli",     # see the help for other choices
           n.trees = 100,                # number of trees
           shrinkage = 0.1,              # shrinkage or learning rate,
           interaction.depth = 2,         # 1: additive model, 2: two-way interactions, etc.
           n.minobsinnode = 200,         # minimum total weight needed in each node
           verbose = TRUE               #  print out progress
  )
  g1_fea = names(relative.influence(g1,   g1$n.trees, sort. = T))
  g1_weight = relative.influence(g1,   g1$n.trees, sort. = T)
  gbm_fea_report = data.frame(fea = g1_fea, importance = g1_weight, rank = rank(-g1_weight))
  #colnames(gbm_fea_report)[1] = paste0('fea_', iter)
  colnames(gbm_fea_report)[2] = paste0('importance_', iter)
  colnames(gbm_fea_report)[3] = paste0('rank_', iter)
  gbm_fea_report
}


lmequ_extract = function(best.model)
{
  eq = paste0(best.model$coefficients[1],' + ',
              paste(best.model$coefficients[2:length(best.model$coefficients)],names(best.model$coefficients[2:length(best.model$coefficients)]),collapse = ' + ',sep = '*'))
  print(eq)
  eq
}
