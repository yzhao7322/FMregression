### A statistic that tests the null hypothesis that the mean of a data of interest statistically deviates from the nominal value.
#   i.e., testing the significance of risk premium.
# input: obj_data - data of interest
#        nullvalue - the nominal value under the null hypothesis
# output: p-value
functionl_t <- function(obj_data,nullvalue){
  # the test T statistics
  H0_Tn<-function(z,nullvalue){
    point_grid=nrow(z)
    n=ncol(z)
    int_approx=function(x){
      temp_n=NROW(x)
      return((1/temp_n)*sum(x))}
    Tn=n*int_approx(colMeans(z)-nullvalue)^2
    return(Tn)
  }

  # get the critical values
  Tn_cv<-function(z,cv_N){
    N=ncol(z)
    point_grid=nrow(z)

    # autocovariance operator from -N to N
    cov_est<-function(x,h){
      N=ncol(x)
      c_est=x[,1:(N-h)]%*%t(x[,(h+1):N])/(N-h)
      return(c_est)
    }

    Gam=0
    Gam=Gam+cov_est(z,0)#long_run_cov2(z, C0 = 3, H = 3)#c

    eig=eigen(Gam)
    eigvals=eig$values
    eigvals=as.vector(eigvals)
    vect=eigvals/point_grid

    int_approx=function(x){
      temp_n=NROW(x)
      return((1/temp_n)*sum(x))}

    lim_sum=matrix(0,cv_N,1)

    lim_sum=0
    for (j in 1:length(vect)){
      lim_sum=lim_sum+vect[j]*rnorm(cv_N,mean=0,sd=1)}

    cv=quantile(lim_sum,probs=c(0.90,0.95,0.99))
    return(list(cv,lim_sum))
  }

  cv_N=5000 # sample size for approximating critical values
  Tn_stat=H0_Tn(obj_data,nullvalue)
  limit=Tn_cv(obj_data,cv_N)
  emplim=limit[[2]]

  return(1-ecdf(emplim)(Tn_stat))
}

