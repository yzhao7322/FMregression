# This file documents the R code for the functional FAMA-MACBETH regression
# Reference: Zhao, Y. (2021). Validating intra-day risk premium in cross-sectional return curves. Finance Research Letters, 102020.

# install.packages(freqdom)
# install.packages(foreign)
# install.packages(reshape2)
# install.packages(data.table)
# install.packages(lattice)
# install.packages(MASS)
# install.packages(memisc)
# install.packages(fda)
# install.packages(rainbow)

#########################################
###### concurrent functional linear model
#########################################
####### step 1
coef_cons = coef_mktbeta = coef_smb = coef_hml = coef_mom = list()

# assume there are two samples of 2015 and 2016
nyear=c(2015,2016)

for (v in 1:2){
  target_year = nyear[v]
  const <- rep(1, dim(eval(parse(text = paste0("idr_ast_",target_year,"_", 1)))$coef)[2] )
  target_mktbeta <- eval(parse(text = paste0("fa_",target_year,"_", 1)))
  target_smb <- eval(parse(text = paste0("fa_",target_year,"_", 2)))
  target_hml <- eval(parse(text = paste0("fa_",target_year,"_", 3)))
  target_mom <- eval(parse(text = paste0("fa_",target_year,"_", 4)))
  target_rfr <- eval(parse(text = paste0("fa_",target_year,"_", 5)))

  c_cons = c_mktbeta = c_smb = c_hml = c_mom = matrix(NA,nbasis,N)
  for (j in 1:N){
     target_asset <- eval(parse(text = paste0("idr_ast_",target_year,"_", j)))
     xfdlist <- list(const=const, mktbeta=target_mktbeta, smb=target_smb, hml=target_hml, mom=target_mom)
     beta0 <- with(target_asset, fd(basisobj=basis, fdnames=fdnames))
     beta1 <- with(target_mktbeta, fd(basisobj=basis, fdnames=fdnames))
     beta2 <- with(target_smb, fd(basisobj=basis, fdnames=fdnames))
     beta3 <- with(target_hml, fd(basisobj=basis, fdnames=fdnames))
     beta4 <- with(target_mom, fd(basisobj=basis, fdnames=fdnames))
     betalist <- list(const=fdPar(beta0), mktbeta=fdPar(beta1), smb=fdPar(beta2),hml=fdPar(beta3),mom=fdPar(beta4))
     yfd = target_asset - target_rfr
     fRegressout <- fRegress(yfd, xfdlist, betalist)
     
     c_cons[,j] = fRegressout$betaestlist$const$fd$coefs
     c_mktbeta[,j] = fRegressout$betaestlist$mktbeta$fd$coefs
     c_smb[,j] = fRegressout$betaestlist$smb$fd$coefs
     c_hml[,j] = fRegressout$betaestlist$hml$fd$coefs
     c_mom[,j] = fRegressout$betaestlist$mom$fd$coefs
  }
     coef_cons[[v]] = c_cons
     coef_mktbeta[[v]] = c_mktbeta
     coef_smb[[v]] = c_smb
     coef_hml[[v]] = c_hml
     coef_mom[[v]] = c_mom
}



####### step 2
gamma_cons = gamma_mktbeta = gamma_smb = gamma_hml = gamma_mom = list()
nbasis2 = 48
for (v in 1:2){
  target_year = nyear[v]
  T = dim(eval(parse(text = paste0("idr_",target_year)))[[1]])[2]

  const <- rep(1, dim(eval(parse(text = paste0("idr_tran_",target_year,"_", 1)))$coef)[2] )
  target_mktbeta <- Data2fd(y=as.matrix(coef_mktbeta[[v]]),  argvals=seq(0, 1, len = nbasis),
                      create.bspline.basis(nbasis = nbasis2))
  target_smb <- Data2fd(y=as.matrix(coef_smb[[v]]),  argvals=seq(0, 1, len = nbasis),
                      create.bspline.basis(nbasis = nbasis2))
  target_hml <- Data2fd(y=as.matrix(coef_hml[[v]]),  argvals=seq(0, 1, len = nbasis),
                      create.bspline.basis(nbasis = nbasis2))
  target_mom <- Data2fd(y=as.matrix(coef_mom[[v]]),  argvals=seq(0, 1, len = nbasis),
                      create.bspline.basis(nbasis = nbasis2))
  target_rfr <- eval(parse(text = paste0("fa_",target_year,"_", 5)))
  
  g_cons = g_mktbeta = g_smb = g_hml = g_mom = matrix(NA,nbasis2,T)
  for (j in 1:T){
     target_asset <- eval(parse(text = paste0("idr_tran_",target_year,"_", j)))
     yfd = apply(target_asset$coefs,2,function(x){x-target_rfr$coefs[,j]})
     yfd = Data2fd(y=as.matrix(yfd),  argvals=seq(0, 1, len = nbasis),
                      create.bspline.basis(nbasis = nbasis2))
     xfdlist <- list(const=const, mktbeta=target_mktbeta, smb=target_smb, hml=target_hml, mom=target_mom)
     beta0 <- with(yfd, fd(basisobj=basis, fdnames=fdnames))
     beta1 <- with(target_mktbeta, fd(basisobj=basis, fdnames=fdnames))
     beta2 <- with(target_smb, fd(basisobj=basis, fdnames=fdnames))
     beta3 <- with(target_hml, fd(basisobj=basis, fdnames=fdnames))
     beta4 <- with(target_mom, fd(basisobj=basis, fdnames=fdnames))
     betalist <- list(const=fdPar(beta0), mktbeta=fdPar(beta1), smb=fdPar(beta2),hml=fdPar(beta3),mom=fdPar(beta4))
     fRegressout1 <- fRegress(yfd, xfdlist, betalist)

     g_cons[,j] = fRegressout1$betaestlist$const$fd$coefs
     g_mktbeta[,j] = fRegressout1$betaestlist$mktbeta$fd$coefs
     g_smb[,j] = fRegressout1$betaestlist$smb$fd$coefs
     g_hml[,j] = fRegressout1$betaestlist$hml$fd$coefs
     g_mom[,j] = fRegressout1$betaestlist$mom$fd$coefs
  }
  gamma_cons[[v]] = g_cons
  gamma_mktbeta[[v]] = g_mktbeta
  gamma_smb[[v]] = g_smb
  gamma_hml[[v]] = g_hml
  gamma_mom[[v]] = g_mom
}
