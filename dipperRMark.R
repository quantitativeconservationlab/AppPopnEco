##### dipper example for CJS in MARK ##########################################
##### load packages
library(RMark)
#### load data ----
data(dipper)
##### data process -------
dipper.process <-process.data(dipper,model="CJS",begin.time=1980,groups="sex")
dipper.ddl <-make.design.data(dipper.process)

dipper.ddl$Phi

# add a flood parameter on the phi ddl:
dipper.ddl$Phi$Flood=0
dipper.ddl$Phi$Flood[dipper.ddl$Phi$time==1981 | dipper.ddl$Phi$time==1982]=1
dipper.ddl$p$Flood=0
dipper.ddl$p$Flood[dipper.ddl$p$time==1982]=1

#now modify the age class in ddl for phi:
dipper.ddl=add.design.data(dipper.process, dipper.ddl,
                           parameter="Phi", type="age", bins=c(0,1,3,6),name="ageclass")
#check it worked:
summary(dipper.ddl$Phi)
# a "(" means the interval is open on the left which means that value is not
# included in the interval. Whereas a square bracket ("[" or "]") is for a closed interval which means the
# interval end point is included

#shift intervals to the left:
dipper.ddl <- add.design.data(dipper.process, dipper.ddl,
                           parameter="Phi", type="age",
                           bins=c(0,1,3,6),name="ageclass",
                           right=FALSE,replace=TRUE)

# In many situations the additional design data are simply covariates to be used in place of
# occasion/time effects. Examples are effort, weather, or observers which vary for occasions and may
# be useful to simplify modeling of capture probability rather than time-varying parameters. For this
# situation, the function merge.occasion.data was created

#dataframe of covariates to be added to ddl:
df <- data.frame( time = c(1980:1986), effort = c(10,5,2,8,1,2,3) )
#combining:
#dipper.ddl=merge.occasion.data(dipper.process,dipper.ddl,"p",df)
dipper.ddl$p <- merge_design.covariates(ddl = dipper.ddl$p, df = df, 
                    bygroup = FALSE, bytime = TRUE )
#check
summary(dipper.ddl$p)


# Females only
dipper.fem.ch.gof <- dipper$ch[dipper$sex == "Female"] %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(dipper[dipper$sex == "Female",]))

dipper.fem.ch.gof[1:10,]

rep(1,nrow(dipper[dipper$sex == "Female",]))
nrow( dipper.fem.ch.gof)
test2ct_fem <- test2ct(dipper.fem.ch.gof, 
                       rep(1,nrow(dipper[dipper$sex == "Female",]))) 

##########################################################################
# GOF
data.dir<-"C:/Users/mike/Dropbox/teaching/software course/wk6"

data.dir<-"C:/document and settings/conroy/my documents/Dropbox/teaching/software course/wk6"
setwd(data.dir)
require(RMark)
data(dipper)
#formulas for general ("global") CJS model (sex*time)
Phi.sex.T=list(formula=~sex*time)
p.sex.T=list(formula=~sex*time)
Phi.t=list(formula=~time)
p.t=list(formula=~time)

dipper.processed=process.data(dipper,model="CJS",begin.time=1980,groups="sex")
dipper.ddl=make.design.data(dipper.processed)
global<-mark(dipper.processed,dipper.ddl,model.parameters=list(Phi=Phi.sex.T,p=p.sex.T))
results<-collect.models()

#deviance / df
c.hat<-global$results$deviance/global$results$deviance.df
c.hat

#export RMark data to MARK format
export.chdata(dipper.processed, filename="dipper",replace=T)


#RELEASE GOF

release.gof(dipper.processed)

#################################BOOTSTRAPPING
#LOAD THESE FUNCTIONS FIRST!!

#function to get number of new releases for each group* occasion
Marked<-function(data=dipper,n.occasions=7,groups=Groups)
{
  group<-data[,2]
  marked<-matrix(nrow=length(Groups),ncol=n.occasions)
  for(g in 1:length(groups))
  {
    data_<-subset(data,group==groups[g])
    data_
    ch<-data_$ch
    for(i in 1:n.occasions)
    {
      ch1<-ch[(as.numeric(substr(ch,1,i)))==1]
      marked[g,i]<-length(ch1)
    }
  }
  return(marked)
}
#get # marked




#####################BOOSTRAPPING

#simulate CJS data for 1 group
simul.cjs<-function(phi,p,marked)
{
  n.occasions<-length(p)+1
  Phi<-matrix(phi,n.occasions-1,nrow=sum(marked),byrow=T)
  P<-matrix(p,n.occasions-1,nrow=sum(marked),byrow=T)
  
  #n.occasions<-dim(Phi)[2]+1
  CH<-matrix(0,ncol=n.occasions,nrow=sum(marked))
  #define a vector with marking occasion
  mark.occ<-rep(1:length(marked),marked[1:length(marked)])
  #fill in CH
  for (i in 1:sum(marked))
  {
    CH[i,mark.occ[i]]<-1
    if (mark.occ[i]==n.occasions) next
    for(t in (mark.occ[i]+1):n.occasions)
    {
      #survive?
      sur<-rbinom(1,1,Phi[i,t-1])
      if(sur==0) break #move to next
      #recaptured?
      rp<-rbinom(1,1,P[i,t-1])
      if(rp==1) CH[i,t]<-1
    } #t
  } #i
  return(CH)
}
###function to create capture history character strings
pasty<-function(x) 
{
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n)
  {
    out[i]<-paste(x[i,],collapse="")
  }
  return(out)
}

### get parameter estimates and number markedfrom Dipper global model 
Groups<-c("Female","Male")
n.occasions<-7
params<-global$results$real
Phi.fem<-params$estimate[1:6]
Phi.mal<-params$estimate[7:12]
p.fem<-params$estimate[13:18]
p.mal<-params$estimate[19:24]

marked<-Marked(data=dipper,n.occasions=7,groups=Groups)

#estimates go into simulation within loop
sims<-function(Phi.fem,Phi.mal,p.fem,p.mal,marked,reps)
{
  deviance<-dim(reps)
  for (i in 1:reps)
  {
    cat("iteration = ", iter <- i, "\n")
    sim.fem<-simul.cjs(Phi.fem,p.fem,marked[1,])
    sim.mal<-simul.cjs(Phi.mal,p.mal,marked[2,])
    fem.hist<-data.frame(ch=pasty(sim.fem),sex="Female")
    mal.hist<-data.frame(ch=pasty(sim.mal),sex="Male")
    sim.data<-rbind(fem.hist,mal.hist)
    sim.processed=process.data(sim.data,model="CJS",groups="sex")
    sim.ddl=make.design.data(sim.processed)
    global.sim<-mark(sim.processed,sim.ddl,model.parameters=list(Phi=Phi.sex.T,p=p.sex.T),output=F,silent=T)
    
    
    deviance[i]<-global.sim$results$deviance
  }
  out<-list(deviance.mean=mean(deviance),deviance.025=quantile(deviance,0.025),deviance.975=quantile(deviance,0.975))
}

sim.out<-sims(Phi.fem,Phi.mal,p.fem,p.mal,marked,reps=100)
sim.out
#compare data deviance and simulation

data.deviance<-global$results$deviance
data.deviance
sim.ci<-c(sim.out$deviance.025,sim.out$deviance.975)
cat("data deviance = ",data.deviance, "simulation mean = ", sim.out$deviance.mean, "simulation 95%CI = ", sim.ci,  "\n")

#modified c.hat
c.hat<-data.deviance/sim.out$deviance.mean
cat("modified c.hat=",c.hat,  "\n")


#adjust by c-hat to form QAIC
results.adj<-adjust.chat(chat=c.hat,results)

#clean up MARK temporary files
rm(list=ls())
cleanup(ask=F)




############## end of script ################################################