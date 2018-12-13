library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(cetcolor)

set.seed(2018)

results_to_df <- function(x){
  colnames_ = names(x)
  colnames = character(0)
  df = data.frame(placeholder_='a')
  for (name in colnames_){
    if ((is.character(x[[name]]) | is.numeric(x[[name]]) | is.logical(x[[name]]))){
      if (length(x[[name]])==1)
        df[,name] = x[[name]]
      else{
        for (i in 1:length(x[[name]])){
          df[,paste(name,  i, sep='__')] = x[[name]][i]
        }
      }
    }
    else if (!is.null(names(x[[name]]))){
      subnames = names(x[[name]])
      for (subname in subnames){
        if ((is.character(x[[name]][[subname]]) | is.numeric(x[[name]][[subname]]) | is.logical(x[[name]][[subname]]))){
          if (length(x[[name]][[subname]]) == 1){
            df[,paste(name, subname, sep='__')] = x[[name]][[subname]]
          }
          else{
            for (i in 1:length(x[[name]][[subname]])){
              df[,paste(name, subname, i, sep='__')] = x[[name]][[subname]][i]
            }
          }
        }
      }
    }
  }
  df[,'placeholder_'] = NULL
  return(df)
}


## EXPERIMENT 1 - TWO IID NORMAL DISTRIBUTIONS, SAME PARAMETERS
experiment_1_trial <- function(min_n=5, max_n = 500, step_size=1, target_p=0.05, ...){
  x1 = rnorm(min_n)
  x2 = rnorm(min_n)
  test_results = t.test(x1,x2,...)
  sample_size = min_n
  while ((test_results$p.value > target_p) & sample_size < max_n){
    sample_size = sample_size + step_size
    x1 = c(x1, rnorm(step_size))
    x2 = c(x2, rnorm(step_size))
    test_results = t.test(x1,x2)
  }
  return(list(sample_size=sample_size, test_results=test_results, hacked=test_results$p.value <= target_p, min_n=min_n, max_n=max_n, step_size=step_size))
}

# creates 3x10,000 matrix, row names are order of list
experiment_1_raw = replicate(1000, experiment_1_trial(), simplify=FALSE)
experiment_1_df = ldply(experiment_1_raw, results_to_df)

experiment_1b_raw = replicate(1000, experiment_1_trial(step_size=2), simplify=FALSE)
experiment_1b_df = ldply(experiment_1b_raw, results_to_df)

experiment_1c_raw = replicate(1000, experiment_1_trial(step_size=5), simplify=FALSE)
experiment_1c_df = ldply(experiment_1c_raw, results_to_df)

experiment_1d_raw = replicate(1000, experiment_1_trial(step_size=10), simplify=FALSE)
experiment_1d_df = ldply(experiment_1d_raw, results_to_df)

experiment_1e_raw = replicate(1000, experiment_1_trial(step_size=100), simplify=FALSE)
experiment_1e_df = ldply(experiment_1e_raw, results_to_df)

# results - in general, about 43% for pure p-hacking, slowly lowers to a minimum of 5% per additional attempt, roughly

## EXPERIMENT 2 - TWO IID NORMAL DISTRIBUTIONS, DIFFERENT MEANS

experiment_2_trial <- function(mu1=10, mu2=11, sd1=1, sd2=1, min_diff=2, min_n=5, max_n = 500, step_size=1, target_p=0.05, ...){
  x1 = rnorm(min_n)
  x2 = rnorm(min_n)
  test_results = t.test(x1,x2,...)
  sample_size = min_n
  while ((test_results$p.value > target_p) & sample_size < max_n){
    sample_size = sample_size + step_size
    x1 = c(x1, rnorm(step_size))
    x2 = c(x2, rnorm(step_size))
    test_results = t.test(x1,x2)
  }
  return(list(sample_size=sample_size, test_results=test_results, hacked=test_results$p.value <= target_p, min_n=min_n, max_n=max_n, step_size=step_size))
}

## EXPERIMENT 3 - normal distribution test

experiment_3_trial <- function(n_sample=5, mu_min = 0, mu_max = 20, sd_min = 0.5, sd_max=5, 
                               mu_func = runif(1, mu_min, mu_max), sd_func=runif(1, sd_min, sd_max),
                               sd_grid_d=0.01, mu_grid_d=0.02, data=NULL){
  if (is.null(data)){
    mu = eval(mu_func)
    sigma = eval(sd_func)
    data = rnorm(n_sample, mu, sigma)
    N = length(data)
  }
  else{
    N = length(data)
    mu = mean(data)
    sigma=sd(data)
    mu_min=mu-3*sigma/sqrt(N)
    mu_max=mu+3*sigma/sqrt(N)
    sigma_min=sigma*0.5^(1/sqrt(N))
    sigma_max=sigma*2^(1/sqrt(N))
  }
    
  # evaluate grid
  param_grid = expand.grid(
    mu=seq(mu_min - (mu_max-mu_min)/5, mu_max + (mu_max-mu_min)/5, mu_grid_d),
    sigma=seq(sd_min*0.8, sd_max*1.2, sd_grid_d)
  )
  param_grid$probs = apply(param_grid, 1, function(x){sum(dnorm(data, x[1], x[2], log=TRUE))})
  
  param_grid$density = exp(param_grid$probs)/sum(exp(param_grid$probs))/(mu_grid_d*sd_grid_d)
  summary_data = data.frame(sigma_s=sd(data), mu_s=mean(data), mu=mu, sigma=sigma)
  return(list(grid=param_grid, mu=mu, sigma=sigma, data=data, sigma_s = sd(data), mu_s = mean(data), summary_data=summary_data))
  
}

times <- data.frame(
  time=c(
    12.167,11.617,13.133,15.250,14.900, # 65%
    17.133,14.400,12.083,16.767,17.550, # 50%
    13.867,14.567,15.467,14.300,16.667 # 0%
  ),
  protection = factor(rep(c('sandbag','chunk','none'), each=5), levels=c('none','sandbag','chunk')),
  coverage = rep(c(0.65, 0.5, 0), each=5)
)

e3_trial_a = experiment_3_trial(data=times$time[1:5])

ggplot(e3_trial_a$grid) + geom_raster(aes(x=mu, y=sigma, fill=1/(0.01*0.02) * exp(probs)/sum(exp(probs)))) + 
  geom_point(data=e3_trial_a$summary_data, aes(x=mu_s, y=sigma_s, color='black')) + scale_color_identity() + 
  scale_fill_gradientn('prob density', colors=cet_pal(7, 'inferno'))

e3_trial_b = experiment_3_trial(data=times$time[6:10])

ggplot(e3_trial_b$grid) + geom_raster(aes(x=mu, y=sigma, fill=1/(0.01*0.02) * exp(probs)/sum(exp(probs)))) + 
  geom_point(data=e3_trial_b$summary_data, aes(x=mu_s, y=sigma_s, color='black')) + scale_color_identity() + 
  scale_fill_gradientn('prob density', colors=cet_pal(7, 'inferno'))

e3_trial_c = experiment_3_trial(data=times$time[11:15])

ggplot(e3_trial_c$grid) + geom_raster(aes(x=mu, y=sigma, fill=1/(0.01*0.02) * exp(probs)/sum(exp(probs)))) + 
  geom_point(data=e3_trial_c$summary_data, aes(x=mu_s, y=sigma_s, color='black')) + scale_color_identity() + 
  scale_fill_gradientn('prob density', colors=cet_pal(7, 'inferno'))

# same mean + var, different data; use anscombe quartet

p1 = anscombe$y1
p2 = anscombe$y2

e3_trial_d = experiment_3_trial(data=p1)

ggplot(e3_trial_d$grid) + geom_raster(aes(x=mu, y=sigma, fill=1/(0.01*0.02) * exp(probs)/sum(exp(probs)))) + 
  geom_point(data=e3_trial_d$summary_data, aes(x=mu_s, y=sigma_s, color='black')) + scale_color_identity() + 
  scale_fill_gradientn('prob density', colors=cet_pal(7, 'inferno'))

e3_trial_e = experiment_3_trial(data=p2)

ggplot(e3_trial_e$grid) + geom_raster(aes(x=mu, y=sigma, fill=1/(0.01*0.02) * exp(probs)/sum(exp(probs)))) + 
  geom_point(data=e3_trial_e$summary_data, aes(x=mu_s, y=sigma_s, color='black')) + scale_color_identity() + 
  scale_fill_gradientn('prob density', colors=cet_pal(7, 'inferno'))

## EXPERIMENT 4 - USING VARIANCE AND RESAMPLING

produce_param_grid <- function(n_sample=5, mu_min = 0, mu_max = 20, var_min = 0.5, var_max=5, 
                               mu_func = runif(1, mu_min, mu_max), var_func=runif(1, var_min, var_max),
                               var_grid_d=0.01, mu_grid_d=0.02, data=NULL){
  if (is.null(data)){
    mu = eval(mu_func)
    sigma2 = eval(var_func)
    data = rnorm(n_sample, mu, sigma2)
    N = length(data)
  }
  else{
    N = length(data)
    mu = mean(data)
    sigma2=var(data)
    mu_min=mu-2*sqrt(sigma2)/sqrt(N)
    mu_max=mu+2*sqrt(sigma2)/sqrt(N)
    var_min=sigma2*0.5^(4/sqrt(N))
    var_max=sigma2*2^(4/sqrt(N))
  }
  
  # evaluate grid
  param_grid = expand.grid(
    mu=seq(mu_min - (mu_max-mu_min)/5, mu_max + (mu_max-mu_min)/5, mu_grid_d),
    sigma2=seq(var_min*0.75, var_max*1.3, var_grid_d)
  )
  param_grid$probs = apply(param_grid, 1, function(x){sum(dnorm(data, x[1], sqrt(x[2]), log=TRUE))})
  
  param_grid$density = exp(param_grid$probs)/sum(exp(param_grid$probs))/(mu_grid_d*var_grid_d)
  summary_data = data.frame(sigma2_s=var(data), mu_s=mean(data), mu=mu, sigma2=sigma2)
  return(list(grid=param_grid, mu=mu, sigma2=sigma2, data=data, sigma2_s = var(data), mu_s = mean(data), summary_data=summary_data))
  
}

e5_test0_grid = produce_param_grid(data=c(1,1,5,5))

ggplot(e5_test0_grid$grid) + geom_raster(aes(x=mu, y=sigma2, fill=density)) + 
  geom_point(data=e5_test0_grid$summary_data, aes(x=mu_s, y=sigma2_s, color='black')) + scale_color_identity() + 
  scale_fill_gradientn('prob density', colors=cet_pal(7, 'inferno'))


sample_from_grid <- function(grid, n=1, replace=FALSE){
  if ('grid' %in% names(grid) & class(grid)=='list')
    grid = grid$grid
  rows = sample_n(grid, size=n, replace=replace, weight=grid$density)
  return(rows)
}

test_params = sample_from_grid(e5_test0_grid, n=100)

experiment_4_trial(data1=NULL, data2=NULL, mu1=NULL, mu2=NULL, var1=NULL, var2=NULL, n=5, target_significance=0.95, target_difference=1){
  if (!is.null(data1)){
    # use data as source distribution and run sampling trials
  }
  else{
    # generate data and then run hypothesis
    if (mu1==mu2){
      # false-positive testing mode
      
    }
    else{
      # regular testing
    }
  }
}
