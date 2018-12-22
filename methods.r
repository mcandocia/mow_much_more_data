library(plyr)
library(dplyr)
library(reshape2)
library(mvtnorm)
library(MASS)


## UTILITY
decorrect <- function(n) return((n-1)/n)

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

generate_from_mixed_distribution <- function(n, mu, sigma, class_variables){
  if (!is.null(mu) & length(mu) > 0)
    data = as.data.frame(rmvnorm(n, mean=mu, sigma=sigma))
  else
    data = data.frame(null_value=rep(NA, n))
  
  if (length(class_variables)> 0){
    for (v in names(class_variables)){
      data[,v] = rbinom(n, 1, class_variables[v])
    }
  }
  return(data)
}



sample_from_grid <- function(grid, n=1, replace=FALSE){
  if ('grid' %in% names(grid) & class(grid)=='list')
    grid = grid$grid
  rows = sample_n(grid, size=n, replace=replace, weight=grid$density)
  return(rows)
}

## EXACT LIKELIHOOD


produce_normal_param_grid <- function(n_sample=5, mu_min = 0, mu_max = 20, var_min = 0.5, var_max=5, 
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

## 2-SAMPLE TEST (SIMPLE)

# (1) create distribution 

# (2) for each possible parameter combination, calculate a few quantile values

# (3) find cutoff value that proportionally achieves desired p-value based on sum of normally-distributed variables

## REJECTION SAMPLING - regression 

# (0) Determine which parameter you are looking for significance

# (1) Establish general form/bounds of multivariate/mixed distribution for input data

# (2) Establish covariance matrix of low-sample-size regression model as well as coefficient estimates

# (3) determine maximum likelihood value for distribution, then randomly sample values of parameters using 
#     rejection sampling, using maximum likelihood as normalization constant

# (4) for each sample of independent variable distribution parameters, run dozens of trials of (5)-(6)

# (5) Significance roughly increases with 

# REGERSSION 

cormat = matrix(c(
  1,0.3,-0.3,
  0.3,1,0,
  -0.3,0,1
),
  nrow=3,
  ncol=3
)

regression_data = as.data.frame(rmvnorm(15, sigma=cormat))
names(regression_data) = c('x1','x2','x3')
regression_data = regression_data %>% mutate(y = x1+0.5*x2 -0.3*x3 + rnorm(15)*1.2)

model = lm(y~x1+x2+x3, data=regression_data)
(regression_result=summary(lm(y~x1+x2+x3, data=regression_data)))

