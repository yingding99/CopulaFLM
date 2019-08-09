flm_surv_beta_smooth_only = function(time, delta, mode = "Additive", geno, pos, order, basis, base = "bspline")  
{
  geno[is.na(geno)]=0
  ### define genotyping matrix for Dom and Rec modes ###
  ### For Dom mode, redefine geno[i,j] = 1 if geno[i,j] = 1 or 2
  ### For Rec mode, redefine geno[i,j] = 1 if geno[i,j] = 2 
  geno_X = geno * 0
  for (i in 1:nrow(geno))
    for (j in 1:ncol(geno))
    {
      if (mode == "Dom")
      {
        if (geno[i, j] == 1 || geno[i, j] == 2)
          geno_X[i,j] = 1
      }
      else if ( mode == "Rec")
        if (geno[i, j] == 2)
          geno_X[i, j] = 1
    }
  
  if  ( mode == "Rec" || mode == "Dom")
    geno = geno_X     
  
  idx     = is.na(time)
  time    = time[!idx]
  delta   = delta[!idx]
  geno    = geno[!idx,]
  
  dqr     = qr(geno)
  index   = dqr$pivot[1:dqr$rank]
  geno    = geno[, index]
  pos     = pos[index]
  nsample = nrow(geno)
  nsnp    = ncol(geno)
  
  if(max(pos) > 1) 
  {
    pos = (pos - min(pos)) / (max(pos) - min(pos))
  }
  
  if (base ==  "bspline"){
    betabasis  = create.bspline.basis(norder = order, nbasis = basis)
  } else if (base == "fspline"){
    betabasis  = create.fourier.basis(c(0,1), nbasis = basis)
  }else { }
  
  B = eval.basis(pos, betabasis)
  
  UJ = geno %*% B
  
  ### Make sure UJ has full rank of bbasis or fbasis ###
  UJdqr   = qr(UJ)
  UJindex = UJdqr$pivot[1:UJdqr$rank]
  UJ      = UJ[, UJindex]
  ###
  
  return(UJ)
}


FR_beta_smooth_only_dat<-function(data, geno, pos, order, basis, base){
  UJ = flm_surv_beta_smooth_only(time=data$obs_time, delta=data$status, mode = "Additive", geno, pos, order, basis, base = "bspline")
  
  dat<-cbind(data,UJ)
  SNP_name<-names(dat)[(ncol(dat)-basis+1):ncol(dat)]
  
  indata1<-dat[dat$enum==1,]
  indata2<-dat[dat$enum==2,]
  
  return(list(data=dat,basis=basis,SNP_name=SNP_name))
}


###############################################################
#############   The Following Code Are for Comparison  ########
###############################################################
flm_surv_beta_smooth_only_cov = function(time, delta, mode = "Additive", geno, pos,order,basis,covariate,base = "bspline",type="lrt",cluster=NULL)  
{
  geno[is.na(geno)]=0
  ### define genotyping matrix for Dom and Rec modes ###
  ### For Dom mode, redefine geno[i,j] = 1 if geno[i,j] = 1 or 2
  ### For Rec mode, redefine geno[i,j] = 1 if geno[i,j] = 2 
  geno_X = geno * 0
  for (i in 1:nrow(geno))
    for (j in 1:ncol(geno))
    {
      if (mode == "Dom")
      {
        if (geno[i, j] == 1 || geno[i, j] == 2)
          geno_X[i,j] = 1
      }
      else if ( mode == "Rec")
        if (geno[i, j] == 2)
          geno_X[i, j] = 1
    }
  
  if  ( mode == "Rec" || mode == "Dom")
    geno = geno_X     
  
  idx     = is.na(time)
  time    = time[!idx]
  delta   = delta[!idx]
  geno    = geno[!idx,]
  
  dqr     = qr(geno)
  index   = dqr$pivot[1:dqr$rank]
  geno    = geno[, index]
  pos     = pos[index]
  nsample = nrow(geno)
  nsnp    = ncol(geno)
  
  if(max(pos) > 1) 
  {
    pos = (pos - min(pos)) / (max(pos) - min(pos))
  }
  
  if (base ==  "bspline"){
    betabasis  = create.bspline.basis(norder = order, nbasis = basis)
  } else if (base == "fspline"){
    betabasis  = create.fourier.basis(c(0,1), nbasis = basis)
  }else { }
  
  B = eval.basis(pos, betabasis)
  
  UJ = geno %*% B 
  
  ### Make sure UJ has full rank of bbasis or fbasis ###
  UJdqr   = qr(UJ)
  UJindex = UJdqr$pivot[1:UJdqr$rank]
  UJ      = UJ[, UJindex]
  ###
  
  pval = list()
  if(length(cluster)>0){
    fit        = coxph(Surv(time, delta) ~ as.matrix( covariate ) + as.matrix(UJ) + cluster(cluster))
    wald<-coxph.wtest(fit$var[(ncol(as.matrix(covariate))+1):(nrow(fit$var)),(ncol(as.matrix(covariate))+1):(nrow(fit$var))],
                      fit$coefficients[(ncol(as.matrix(covariate))+1):(nrow(fit$var))])
    pval = signif(pchisq(wald$test,df=wald$df,lower.tail = F),4)
    pval
  } else if (type == "lrt") {
    fit        = coxph(Surv(time, delta) ~ as.matrix( covariate ) + as.matrix(UJ))
    pval = signif(pchisq(anova(fit, test = "Chisq")[3, 2],df=basis,lower.tail = F),4)
    pval
  } else if (type == "wald") {
    fit        = coxph(Surv(time, delta) ~ as.matrix( covariate ) + as.matrix(UJ))
    wald<-coxph.wtest(fit$var[(ncol(as.matrix(covariate))+1):(nrow(fit$var)),(ncol(as.matrix(covariate))+1):(nrow(fit$var))],
                      fit$coefficients[(ncol(as.matrix(covariate))+1):(nrow(fit$var))])
    pval = signif(pchisq(wald$test,df=wald$df,lower.tail = F),4)
    pval
  }
}


