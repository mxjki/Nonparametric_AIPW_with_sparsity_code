library(listdtr)
library(grplasso)
library(foreach)
library(parallel)
library(doParallel)
library(kernlab)
source("Auxiliary.R")


i = 1

np_index = 2
model_index = 2
np_set = data.frame(n = c(800, 1000,  800, 1000),
                    p = c(400, 400, 2000, 2000))

n.run = 800
p.run = 400

run.x <- matrix(runif(n.run*p.run,-0.5,0.5),n.run,p.run)


f1<-run.x[,1] 
f2<-run.x[,2]
f3<-run.x[,3]
f4<-run.x[,4]
run.y<-5*f1 +6*f2+ 4*f3 + 4*f4+ rnorm(n.run,0,1)
mis.prob = min.fun((1+ exp(0.1-2*run.x[,1] - 2*run.x[,3]))^(-1),1) # Missing probability

delta = rbinom(n.run,1,mis.prob) # Missing indicator 
index = c(NA,1:p.run) # Index for the covariate

## First, consider the group lasso for logistic model.
x_inc_inter = cbind(1,run.x) # Construct a design matrix
lambda <- lambdamax(x_inc_inter, y = delta, index = index, penscale = sqrt,
                    model = LogReg()) * 0.5^(0:10) # Tuning parameter sets 
neg_log = rep(0,11)
for (j in 1:5){ # 5-fold CV to select the tuning parameter for group lasso.
  print(j)
  test_index = ((j-1)*n.run/5+1):(j*n.run/5)
  train_index = setdiff(1:n.run,test_index)
  
  x_train_j = x_inc_inter[train_index,]
  delta_train_j = delta[train_index]
  
  x_test_j = x_inc_inter[test_index,]
  delta_test_j = delta[test_index]
  
  
  fit <- grplasso(x_train_j, y = delta_train_j, index = index, lambda = lambda, model = LogReg(),
                  penscale = sqrt,
                  control = grpl.control(update.hess = "lambda", trace = 0))
  
  predict_j = predict(fit,newdata = x_test_j,type = "response")
  
  neg_log = neg_log + apply(log(predict_j) * delta_test_j,2,mean) + apply(log(1-predict_j) * (1-delta_test_j),2,mean)
}
lambda_resp = lambda[which.max(neg_log)] # Selected tuning parameter 

fit_resp <- grplasso(x_inc_inter, y = delta, index = index, lambda = lambda_resp, model = LogReg(),
                     penscale = sqrt,
                     control = grpl.control(update.hess = "lambda", trace = 0))

### Obtain the fitted response probability.
fitted_resp = fit_resp$fitted # The fitted response rate
resp_selected = (0:p.run)[!(fit_resp$coefficients == 0)] # The selected features
### Next, we fit KRR using observed data
resp_index = delta == 1

x_obs = run.x[resp_index,]
y_obs = run.y[resp_index]
n_obs = sum(resp_index)
weight_obs = (fitted_resp[resp_index])^(-1)

x_miss = run.x[!resp_index,]
y_miss = run.y[!resp_index]

KRR_fit = main.fun(x_obs,y_obs,x_miss) #We do not need the weights here.
KRR_sel = (1:p.run)[KRR_fit[1][[1]]==1] #The selected features.
lamb_sel = KRR_fit$lam_sel

x_obs_sel = x_obs[,KRR_sel]
x_mis_sel = x_miss[,KRR_sel]
obj <- krr(x_obs_sel, y_obs) # Fit a kernel model only using the selected covariates. 

fit_obs_simple_krr = predict(obj,x_obs_sel)
fit_mis_simple_krr = predict(obj,x_mis_sel)
m_hat_simple_krr = c(fit_obs_simple_krr,fit_mis_simple_krr)

m_hat_simple_krr[1:n_obs] = fit_obs_simple_krr  + (y_obs - fit_obs_simple_krr)*weight_obs
est_simple_krr = mean(m_hat_simple_krr) # The proposed method
var_simple_krr = var(m_hat_simple_krr)/n.run # Variance estimator