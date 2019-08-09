indata1<-read.csv("lefteye.csv")
indata2<-read.csv("righteye.csv")
dat<-rbind(indata1,indata2)
dat$ind<-ifelse(dat$EYE=="L",1,2)
dat$id<-dat$SampleID
dat$obs_time<-dat$time

##################################################################################
###############                   Table 3: 4 regions                  ############
##################################################################################
library(fda)
library(CopulaCenR)
library(seqMeta)
library(MASS)
library(survival)
source("fun_copFLM.R")
order=4
basisnum=c(5,6,7)
base="bspline"
var_list=c("ENROLLAGE_std","SevScaleBL_std")
gene1<-read.csv("CFH_5000.csv")
gene2<-read.csv("C2_5000.csv")
gene3<-read.csv("ARMS2region_5000.csv")
gene4<-read.csv("C3_5000.csv")

genels<-list(gene1,gene2,gene3,gene4)

p.cop<-rep(NA,3*4)
p.FLMRst<-rep(NA,3*4)
p.uniFLM<-rep(NA,3*4)
skat.all<-rep(NA,4)
burden.all<-rep(NA,4)

for ( i in 1:4){
  gene<-genels[[i]]
  geno<-as.data.frame(t(gene[,-c(1:5)]))
  geno[is.na(geno)]=0
  geno_all<-geno[,(colMeans(geno)/2)>=0.01]
  snp_all<-as.character(gene$SNP)[(colMeans(geno)/2)>=0.01]
  pos_all<-gene$Pos[(colMeans(geno)/2)>=0.01]
  geno = as.matrix(rbind(geno_all,geno_all))
  pos=c(pos_all,pos_all)
  for (j in 1:3){
    copula_flm_mod<-FR_beta_smooth_only_dat(data=dat,geno=geno,pos=pos,order=order,basis=basisnum[j],base="bspline")
    copula2_null <- rc_par_copula(data = copula_flm_mod$data, var_list = var_list, copula = "Clayton", m.dist = "Weibull")
    p.cop[3*i-3+j]<-score_copula(object = copula2_null, var_score = copula_flm_mod$SNP_name)[2]
    p.FLMRst[3*i-3+j]<-flm_surv_beta_smooth_only_cov(time=dat$obs_time, delta=dat$status, mode = "Additive", geno = geno, pos=pos, order=order, basis=basisnum[j], covariate = dat[,var_list], base = "bspline",cluster=dat$id)
    p.uniFLM[3*i-3+j]<-flm_surv_beta_smooth_only_cov(indata1$time, indata1$status, mode = "Additive", geno=as.matrix(geno_all), pos=pos_all, order, basisnum[j], covariate = as.matrix(indata1[,var_list]), base = "bspline",cluster=NULL)
  }
  VariantInfo <- data.frame()
  VariantInfo[1:dim(geno_all)[2], "Chr"] <- rep(1, dim(geno_all)[2])
  VariantInfo[1:dim(geno_all)[2], "Name"] <- names(geno_all)
  VariantInfo[1:dim(geno_all)[2], "gene"] <- rep("gene1", dim(geno_all)[2])
  
  skatobj<-prepCox(Z = as.matrix(geno_all), Surv(time, status)~ SevScaleBL_std + ENROLLAGE_std, SNPInfo = VariantInfo, data = indata1)
  skat.all[i] <- skatMeta(skatobj, SNPInfo = VariantInfo)$p
  burden.all[i] <- burdenMeta(skatobj, SNPInfo = VariantInfo)$p
}

tab.FLM<-as.data.frame(cbind(c(rep("CFH",3),rep("C2",3),rep("ARMS2region",3),rep("C3",3)),rep(c(5,6,7),4),p.cop,p.FLMRst,p.uniFLM))
tab.skat.burden<-as.data.frame(cbind(skat.all,burden.all))
names(tab.skat.burden)<-c("skat.all","burden.all")
row.names(tab.skat.burden)<-c("CFH","C2","ARMS2region","C3")



####################################################################################
#################             Web Figure 1: Gamma Plot                   ###############
####################################################################################
library(fda)
library(CopulaCenR)
order=4
basis=6
base="bspline"
var_list=c("ENROLLAGE_std","SevScaleBL_std")

require(splines)
betabasis  = create.bspline.basis(norder = order, nbasis = basis)
x<-seq(0,1,length.out=100)
B = eval.basis(x, betabasis)

gammadat<-matrix(rep(NA,51*5),ncol=5)

for (i in 1:4){
  gene<-genels[[i]]
  geno<-as.data.frame(t(gene[,-c(1:5)]))
  geno[is.na(geno)]=0
  geno_all<-geno[,(colMeans(geno)/2)>=0.01]
  snp_all<-as.character(gene$SNP)[(colMeans(geno)/2)>=0.01]
  pos_all<-gene$Pos[(colMeans(geno)/2)>=0.01]
  geno = as.matrix(rbind(geno_all,geno_all))
  pos=c(pos_all,pos_all)
  copula_flm_mod<-FR_beta_smooth_only_dat(data=dat,geno=geno,pos=pos,order=order,basis=basis,base="bspline")
  cop_model<-rc_par_copula(data = copula_flm_mod$data, var_list = c(var_list,copula_flm_mod$SNP_name), copula = "Clayton", m.dist = "Weibull",method="BFGS")
  coef<-cop_model$summary[5:(5+basis-1),1]
  gammadat[,1]<-predict(interpSpline(x,B%*%coef))$x
  gammadat[,(1+i)]<-predict(interpSpline(x,B%*%coef))$y
}


matplot(gammadat[,1], gammadat[,-1],type="l",pch=1,col = c(2:4,"burlywood4"),lty=2:5,xlab="standardized position",ylab="")
legend("topleft", legend = c("CFH","C2","ARMS2","C3"), col=c(2:4,"burlywood4"), pch=1,lty=2:5)
abline(h=0,lty=1,col="black")

