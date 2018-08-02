#################################################################################
#  sims.r                                                                       # 
#                                                                               #
#  Simulation code for paper                                                    #
# "Canonical causal diagrams to guide the treatment of missing data in          #
#  epidemiological studies"                                                     #
#                                                                               #
#  Margarita Moreno-Betancur, 18 May 2018                                       #
#                                                                               #
#################################################################################

## Load required packages, simulation model parameters and auxiliary functions

rm(list = ls())
library(boot)
library(mice)
parms<-read.csv("parms.csv")
source("sim_fncs.R")

## Set meta-parameters

set.seed(29387)
n<-200
M<-50
Iter<-5
sim<-1000
sce<-"Scenario1"   # "Scenario1" "Scenario2" "Scenario3"

## Modify simulation model parameters according to scenario

if(sce!="Scenario1") 
{
  nam<-(substr(names(parms),1,1)%in%c("T",LETTERS[1:10])&
               substr(names(parms),4,4)!="I" & substr(names(parms),5,5)!="I"&
          substr(names(parms),4,4)!="W" & substr(names(parms),5,5)!="W")|
        names(parms)%in%c("deltaZ2", "deltaX" , "deltaY" )
  
  parms[nam]<-parms[nam]*ifelse(sce=="Scenario2",2,3)

}
  
attach(parms)
  

## Create datasets to collect results

resCD<-data.frame()
resT<-data.frame()
resA<-data.frame()
resB<-data.frame()
resC<-data.frame()
resD<-data.frame()
resE<-data.frame()
resF<-data.frame()
resG<-data.frame()
resH<-data.frame()
resI<-data.frame()
resJ<-data.frame()

rpmT<-data.frame()
rpmA<-data.frame()
rpmB<-data.frame()
rpmC<-data.frame()
rpmD<-data.frame()
rpmE<-data.frame()
rpmF<-data.frame()
rpmG<-data.frame()
rpmH<-data.frame()
rpmI<-data.frame()
rpmJ<-data.frame()


## Start simulation loop

for(i in 1:sim)
{

  ### GENERATE DATA FROM c-DAG

  U<-rnorm(n,0,1)
  Z1<-rbinom(n,1,inv.logit(Z1_Int+Z1_U*U))
  Z2<-rbinom(n,1,inv.logit(Z2_Int+Z2_U*U))
  X<-rbinom(n,1,inv.logit(X_Int+X_Z1*Z1+X_Z2*Z2)) 
  Y<-rnorm(n,mean=Y_Int+Y_Z1*Z1+Y_Z2*Z2+Y_X*X,sd=Y_SD) 

  ### COMPLETE DATA ANALYSIS

  resCD<-rbind(resCD,c(mean(X),mean(Y),coef(lm(Y~X+Z1+Z2))["X"],NA,NA,NA))

  ### GENERATE MISSING DATA AND ANALYSE FOR EACH m-DAG
  
  for(let in c("T",LETTERS[1:10]))
  {
    W<-rnorm(n,0,1)

    if(let=="T")
    {
      MZ2<-rbinom(n,1,inv.logit(TZ2_Int+TZ2_W*W)) 
      MX<-rbinom(n,1,inv.logit(TX_Int+TX_W*W)) 
      MY<-rbinom(n,1,inv.logit(TY_Int+TY_W*W)) 
      
       
    }else if(let=="A")
    {
      MZ2<-rbinom(n,1,inv.logit(AZ2_Int+AZ2_Z1*Z1+AZ2_W*W))
      MX<-rbinom(n,1,inv.logit(AX_Int+AX_Z1*Z1+AX_W*W))
      MY<-rbinom(n,1,inv.logit(AY_Int+AY_Z1*Z1+AY_W*W))

    }else if(let=="B")
    {
      MZ2<-rbinom(n,1,inv.logit(BZ2_Int+BZ2_Z1*Z1+BZ2_X*X+BZ2_W*W))
      MX<-rbinom(n,1,inv.logit(BX_Int+BX_Z1*Z1+BX_Z2*Z2+BX_W*W))
      MY<-rbinom(n,1,inv.logit(BY_Int+BY_Z1*Z1+BY_Z2*Z2+BY_X*X+BY_W*W))

    }else if(let=="C")
    {
      MZ2<-rbinom(n,1,inv.logit(CZ2_Int+CZ2_Z1*Z1+CZ2_X*X+CZ2_Y*Y+CZ2_W*W))
      MX<-rbinom(n,1,inv.logit(CX_Int+CX_Z1*Z1+CX_Z2*Z2+CX_Y*Y+CX_W*W))
      MY<-rbinom(n,1,inv.logit(CY_Int+CY_Z1*Z1+CY_Z2*Z2+CY_X*X+CY_W*W))
      

    }else if(let=="D")
    {

      MZ2<-rbinom(n,1,inv.logit(DZ2_Int+DZ2_Z1*Z1+deltaZ2*Z2+DZ2_W*W))
      MX<-rbinom(n,1,inv.logit(DX_Int+DX_Z1*Z1+deltaX*X+DX_W*W))
      MY<-rbinom(n,1,inv.logit(DY_Int+DY_Z1*Z1+DY_W*W))

    }else if(let=="E")
    {
      MZ2<-rbinom(n,1,inv.logit(EZ2_Int+EZ2_Z1*Z1+EZ2_X*X+deltaZ2*Z2+EZ2_W*W))
      MX<-rbinom(n,1,inv.logit(EX_Int+EX_Z1*Z1+EX_Z2*Z2+deltaX*X+EX_W*W))
      MY<-rbinom(n,1,inv.logit(EY_Int+EY_Z1*Z1+EY_Z2*Z2+EY_X*X+EY_W*W))
      

    }else if(let=="F")
    {

      MZ2<-rbinom(n,1,inv.logit(FZ2_Int+FZ2_Z1*Z1+deltaZ2*Z2+FZ2_Y*Y+FZ2_W*W))
      MX<-rbinom(n,1,inv.logit(FX_Int+FX_Z1*Z1++deltaX*X+FX_Y*Y+FX_W*W))
      MY<-rbinom(n,1,inv.logit(FY_Int+FY_Z1*Z1+FY_W*W))
      
    }else if(let=="G")
    {
      MZ2<-rbinom(n,1,inv.logit(GZ2_Int+GZ2_Z1*Z1+GZ2_X*X+GZ2_W*W))
      MX<-rbinom(n,1,inv.logit(GX_Int+GX_Z1*Z1+GX_Z2*Z2+GX_W*W))
      MY<-rbinom(n,1,inv.logit(GY_Int+GY_Z1*Z1+GY_Z2*Z2+GY_X*X+deltaY*Y+GY_W*W))

    }else if(let=="H")
    {

      MZ2<-rbinom(n,1,inv.logit(HZ2_Int+HZ2_Z1*Z1+HZ2_X*X+HZ2_Y*Y+HZ2_W*W))
      MX<-rbinom(n,1,inv.logit(HX_Int+HX_Z1*Z1+HX_Z2*Z2+HX_Y*Y+HX_W*W))
      MY<-rbinom(n,1,inv.logit(HY_Int+HY_Z1*Z1+HY_Z2*Z2+HY_X*X+deltaY*Y+HY_W*W))
      

    }else if(let=="I")
    {
      MZ2<-rbinom(n,1,inv.logit(IZ2_Int+IZ2_Z1*Z1+IZ2_X*X+IZ2_Y*Y+deltaZ2*Z2+IZ2_W*W))
      MX<-rbinom(n,1,inv.logit(IX_Int+IX_Z1*Z1+IX_Z2*Z2+IX_Y*Y+deltaX*X+IX_W*W))
      MY<-rbinom(n,1,inv.logit(IY_Int+IY_Z1*Z1+IY_Z2*Z2+IY_X*X+IY_W*W))
      
    }else if(let=="J")
    {
      MZ2<-rbinom(n,1,inv.logit(JZ2_Int+JZ2_Z1*Z1+deltaZ2*Z2+JZ2_W*W))
      MX<-rbinom(n,1,inv.logit(JX_Int+JX_Z1*Z1+deltaX*X+JX_W*W))
      MY<-rbinom(n,1,inv.logit(JY_Int+JY_Z1*Z1+deltaY*Y+JY_W*W))
      
      
    }

    ## Analyse data
    nr<-paste("res",let,sep="")
    assign(nr,analyse(get(nr)))
    
    ## Calculate proportions of missingness
    nm<-paste("rpm",let,sep="")
    assign(nm,missprop(get(nm)))

  }
  print(i)
  flush.console()
}

## Calculate performance indicators

RES<-data.frame()
RES<-rbind(RES,res_table0(resCD))
bench<-unlist(RES[1,c(1,5,9)])
for(ll in c("T",LETTERS[1:10]))
{
  rr<-get(paste("res",ll,sep=""))
  RES<-rbind(RES,res_table(rr[seq(1,2*sim-1,2),]))
  RES<-rbind(RES,res_table(rr[seq(2,2*sim,2),]))
}
RES[,c(1,5,9)]<-format(round(RES[,c(1,5,9)],2),trim=T)
RES[,c(2,6,10)]<-format(round(RES[,c(2,6,10)],1),trim=T)
RES[,c(3,7,11)]<-format(round(RES[,c(3,7,11)],2),trim=T)
RES[,c(4,8,12)]<-format(round(RES[,c(4,8,12)],1),trim=T)

names(RES)<-c("MeanX","%Bias","SE","StdBias","MeanY","%Bias","SE","StdBias","RegCoeff","%Bias","SE","StdBias")
RES<-cbind(DAG=c("Complete data","Trivial",rep("",1),"(a)",rep("",1),"(b)",rep("",1),"(c)",rep("",1),"(d)",rep("",1),"(e)",rep("",1),
                 "(f)",rep("",1),"(g)",rep("",1),"(h)",rep("",1),"(i)",rep("",1),"(j)",rep("",1)),
           Strategy=c("",rep(c("Complete case","MICE"),11)),RES)

write.csv(RES, paste("RES_sim_",sce,".csv",sep=""))

## Calculate mean fraction of missing information per parameter

FMI<-data.frame()
for(ll in c("T",LETTERS[1:10]))
{
  rr<-get(paste("res",ll,sep=""))
  FMI<-rbind(FMI,colMeans(rr[seq(2,2*sim,2),4:6]))
}
names(FMI)<-c("MeanX","MeanY","RegCoeff")
write.csv(FMI, paste("FMI_sim_",sce,".csv",sep=""))

## Calculate summary stats about missing data

RPM<-data.frame()
for(ll in c("T",LETTERS[1:10]))
{
  rr<-get(paste("rpm",ll,sep=""))
  RPM<-rbind(RPM,colMeans(rr))
}
names(RPM)<-c("Z2","X","Y","Any")
write.csv(RPM, paste("RPM_sim_",sce,".csv",sep=""))
