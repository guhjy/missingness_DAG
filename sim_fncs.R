#################################################################################
#  sim_fncs.R                                                                   # 
#                                                                               #
#  Auxiliary functions for simulation code for paper                            #
# "Canonical causal diagrams to guide the treatment of missing data in          #
#  epidemiological studies"                                                     #
#                                                                               #
#  Margarita Moreno-Betancur, 18 May 2018                                       #
#                                                                               #
#################################################################################


# Apply MICE
apply_mice<-function(res,dat)
{
dat$Z1<-as.factor(dat$Z1)
dat$Z2<-as.factor(dat$Z2)
dat$X<-as.factor(dat$X)
imp<-mice(dat,m=M, method=c("","logreg","logreg","norm"), maxit=Iter,printFlag=F)
res<-rbind(res,c(pool(with(imp,lm(I(X==1)~1)))$qbar,pool(with(imp,lm(Y~1)))$qbar,
                   pool(with(imp,lm(Y~X+Z1+Z2)))$qbar[2],
                   pool(with(imp,lm(I(X==1)~1)))$fmi,pool(with(imp,lm(Y~1)))$fmi,
                   pool(with(imp,lm(Y~X+Z1+Z2)))$fmi[2]))
return(res)
}

# Apply avalaible case analysis
apply_CV<-function(res,dat)
{
  res<-rbind(res,c(mean(dat$X,na.rm=T),mean(dat$Y,na.rm=T),
                   coef(lm(Y~X+Z1+Z2,data=dat,na.action=na.omit))["X"],NA,NA,NA))
  return(res)
}

# Wrapper function to analyse
analyse<-function(res)
{
dat<-data.frame(Z1,Z2,X,Y)
dat$Z2[MZ2==1]<-NA
dat$X[MX==1]<-NA
dat$Y[MY==1]<-NA
if(all(is.na(dat$Z2))|all(is.na(dat$X))|all(is.na(dat$Y))|
   all(apply(dat[,c(2,3,4)],1,function(x)return(any(is.na(x))))))
{
res<-rbind(res,rep(NA,6))
res<-rbind(res,rep(NA,6))} else
{
res<-apply_CV(res,dat)
res<-apply_mice(res,dat)}
return(res)
} 

# Get benchmark for bias calculations from complete data analysis  
res_table0<-function(res)
{
  resRow<-vector()
  for(k in 1:3)
    resRow<-c(resRow,mean(res[,k],na.rm=T),0,sd(res[,k],na.rm=T),0)
  return(resRow)}

# Calculate performance indicators for each approach
res_table<-function(res)
{
  resRow<-vector()
  for(k in 1:3)
    resRow<-c(resRow,mean(res[,k],na.rm=T),100*(mean(res[,k],na.rm=T)-bench[k])/bench[k],sd(res[,k],na.rm=T),
              100*(mean(res[,k],na.rm=T)-bench[k])/sd(res[,k],na.rm=T))
  return(resRow)}

# Calculate proportions of missingness
missprop<-function(rpm)
{
  dat<-data.frame(Z1,Z2,X,Y)
  dat$Z2[MZ2==1]<-NA
  dat$X[MX==1]<-NA
  dat$Y[MY==1]<-NA
  rpm<-rbind(rpm,c(sum(is.na(dat$Z2)),sum(is.na(dat$X)),sum(is.na(dat$Y)),
                   sum(is.na(dat$Z2)|is.na(dat$X)|is.na(dat$Y)))/nrow(dat))
  return(rpm)
} 