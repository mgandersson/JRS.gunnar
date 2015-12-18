#Functions used in the surveillance script

################################################################################
###################  define functions ##########################################
################################################################################

make.time.series  <- function(countvector,t.start="1.01.2006",t.end="31.12.2013",level="week"){
  #countvector <- my.gridinfo[[100]]$counts
  #i<-2
  my.startdate<- as.Date(strptime(t.start, "%d.%m.%Y"))
  my.enddate<- as.Date(strptime(t.end, "%d.%m.%Y"))
  counts<- rep(0,(my.enddate-my.startdate+1))
  week.counts <- rep(0,max(timeline$week_abs))
  my.dates<-list()
  #my.diffs<-list()
  if(!is.null(countvector)){
  for(i in 1:length(countvector)){
    #check against timeline
    my.dates[[i]]<- as.Date(strptime(countvector[[i]], "%d.%m.%Y"))
    #my.dates[[i]]<- as.Date(strptime(countvector[[i]], "%Y-%m-%d"))
    my.absweek <- timeline$week_abs[match( my.dates[[i]], as.Date(strptime(timeline$Date, "%Y-%m-%d") ))]
    #my.absweek <- timeline$week_abs[match( my.dates[[i]], as.Date(strptime(timeline$Date, "%d.%m.%Y") ))]
   
    my.diff<- as.Date(strptime(countvector[i], "%d.%m.%Y"))-my.startdate+1
    #my.diff<- as.Date(strptime(countvector[i], "%Y-%m-%d"))-my.startdate+1
    counts[as.numeric(my.diff)]<- counts[as.numeric(my.diff)]+1
    week.counts[my.absweek]<- week.counts[my.absweek]+1
 
    }
  
  }
return(list(my.dates,counts,week.counts))
}



## make.time.series <- function(countvector,t.start="1.01.2006",t.end="31.12.2013",level="week"){
##   #countvector <- my.gridinfo[[100]]$counts
##   #i<-2
##   my.startdate<- as.Date(strptime(t.start, "%d.%m.%Y"))
##   my.enddate<- as.Date(strptime(t.end, "%d.%m.%Y"))
##   counts<- rep(0,(my.enddate-my.startdate+1))
##   week.counts <- rep(0,max(timeline$week_abs))
##   my.dates<-list()
##   #my.diffs<-list()
##   if(!is.null(countvector)){
##   for(i in 1:length(countvector)){
##     #check against timeline
##     #my.dates[[i]]<- as.Date(strptime(countvector[[i]], "%d.%m.%Y"))
##     my.dates[[i]]<- as.Date(strptime(countvector[[i]], "%Y-%m-%d"))
##     my.absweek <- timeline$week_abs[match( my.dates[[i]], as.Date(strptime(timeline$Date, "%Y-%m-%d") ))]
    
##     #my.diff<- as.Date(strptime(countvector[i], "%d.%m.%Y"))-my.startdate+1
##     my.diff<- as.Date(strptime(countvector[i], "%Y-%m-%d"))-my.startdate+1
##     counts[as.numeric(my.diff)]<- counts[as.numeric(my.diff)]+1
##     week.counts[my.absweek]<- week.counts[my.absweek]+1
  
##     }
   
##   }
## return(list(my.dates,counts,week.counts))
## }

### Use one make.time.series with pattern recognition.
   #grepl(".", countvector[[1]]) 
###

make.time.series2 <- function(countvector,t.start="1.01.2006",t.end="31.12.2013",level="week"){
  #countvector <- my.gridinfo[[100]]$counts
  #i<-2
  my.startdate<- as.Date(strptime(t.start, "%d.%m.%Y"))
  my.enddate<- as.Date(strptime(t.end, "%d.%m.%Y"))
  counts<- rep(0,(my.enddate-my.startdate+1))
  week.counts <- rep(0,max(timeline$week_abs))
  my.dates<-list()
  #my.diffs<-list()
  if(!is.null(countvector)){
  for(i in 1:length(countvector)){
    #check against timeline
    my.dates[[i]]<- as.Date(strptime(countvector[[i]], "%d.%m.%Y"))
    #my.dates[[i]]<- as.Date(strptime(countvector[[i]], "%Y-%m-%d"))
    my.absweek <- timeline$week_abs[match( my.dates[[i]], as.Date(strptime(timeline$Date, "%Y-%m-%d") ))]
    
    my.diff<- as.Date(strptime(countvector[i], "%d.%m.%Y"))-my.startdate+1
    #my.diff<- as.Date(strptime(countvector[i], "%Y-%m-%d"))-my.startdate+1
    counts[as.numeric(my.diff)]<- counts[as.numeric(my.diff)]+1
    week.counts[my.absweek]<- week.counts[my.absweek]+1
  
    }
   
  }
return(list(my.dates,counts,week.counts))
}

#############################################################
###    fit model                                           ##
#############################################################
make.fit <- function(counts,loghistmean,t,season.vector,logpop, time.unit = "week", distribution ="poisson"){
# prepare variables
  loghistmean <- loghistMean.X[[i]]  
  if (time.unit =="week"){
    divisor <- 53
  }else{
    divisor <- 365
  }
  
  mysin <- sin(2*pi*season.vector/divisor)
  mycos <- cos(2*pi*season.vector/divisor)

   # define model 
   predictors <-"mysin + mycos+loghistmean"
   model = paste("counts~",predictors) 

   # builds model string for fitting Poisson and Negminom model
   zeromodel = paste("counts~",predictors,"|",predictors) #builds string for fitting zero inflated negative binomial model

   if(distribution=="poisson"){
    fit <- glm(  glm(as.formula(model),distribution)) 
    #fits2[[i]] <- glm( counts ~ sin + periods, "poisson") 
    base.theta <- 100000
    predicts <- predict.glm(fit,type="response", se.fit=TRUE,newdata=as.data.frame(t))

    #calculate upper and lower conficdence intervals
    intervals <- rbind( qpois(0.95,predicts$fit,lower.tail = TRUE, log.p = FALSE),qpois(0.95,predicts$fit,lower.tail = FALSE,log.p = FALSE))
  
    } else if(distribution=="NB"){
       #fits[[i]] <- glm.nb( counts ~ periods) #Negbinom
       fit <- glm.nb( counts ~ sin, trace = TRUE) #Negbinom 
       predicts <- predict(fit,type="response", se.fit=TRUE,newdata=as.data.frame(t))
  
    } else if(distribution=="ZINB"){
       #fits[[i]] <- glm.nb( counts ~ periods) #Negbinom
       #zeroinfl.control(method = "BFGS", maxit = 10000, trace = FALSE,
       #EM = TRUE,  start = list(count = c(0.001,0.001), zero = 0.001))
       fit <- zeroinfl(as.formula(zeromodel),dist="negbin", EM=TRUE,link = "log")    #for the first model for example

    } else if(distribution=="ZIPois"){
       #fits[[i]] <- glm.nb( counts ~ periods) #Negbinom
       #zeroinfl.control(method = "BFGS", maxit = 10000, trace = FALSE,
       #                 EM = FALSE,  start = list(count = c(0.001,0.001), zero = c(0.001,0.001)))
       fit <- zeroinfl(as.formula(zeromodel),dist="poisson", EM=FALSE,link = "log",start = list(count = c(0.001,0.001), zero = c(0.001,0.001)))    #for the first model for example
       base.theta <- 100000
    }

  #return(list(fit,predicts,intervals,base.theta))
  return(list(fit,predicts)) 
}

###################prob of evidence ########################
trunk.geom =function(n,p){
  #p<- 0.90
  #n<- 1
  f0<-  p
  my.pi <- 1/(1-f0)  # f(0) = 1-(1-p)^1  == p
  Fx<-  1-(1-p)^(abs(n)+1)
  #Fxminus1 <-ifelse (n==0,0,Fxminus1 <- 1-(1-p)^(abs(n-1)+1))
  Fxminus1 <- 1-(1-p)^(abs(n-1)+1)
  fx <- Fx - Fxminus1
  fztx <- my.pi * fx
  return(ifelse(n==0,0,fztx))
}

zero.trunk.geom <- function(n,p,p.nonzero){
  p.if.nonzero <- trunk.geom(n,p)
  p <- p.nonzero*p.if.nonzero + (1-p.nonzero)*as.numeric(n==0)
  return(p)
}

################a function for probability of n counts from distribution ########################
p.n <- function(n,lambda,theta=FALSE,p.zero=0,distr="geometrical",print=FALSE){
  if(distr=="poisson"){
    p.nonzero <- ppois(n-1,lambda = lambda,lower=FALSE)-ppois(n,lambda = lambda,lower=FALSE)
    p <- (1-p.zero)*p.nonzero + p.zero*as.numeric(n==0)  
  }else{
    #ZINB model
    p.nonzero <- pnbinom(n-1, mu = lambda, size=theta,lower=FALSE)-pnbinom(n,mu = lambda, size=theta,lower=FALSE) 
    p <- (1-p.zero)*p.nonzero + p.zero*as.numeric(n==0)    
  } 
  if(print==TRUE){
    return(paste(n,lambda,theta,p.zero,"returns",p))
  }else{
    return(p)
  }
}

################a function for joint probability.########################
joint.prob <- function(p1,p2){
  #p1 and p2 are vectors of probabilities
  #calcultae probability of observe n cases during outbreak
  p.joint <- rep(0,length(p1)+length(p2)) #motsvarar pncase
  for(c in 1:length(p.joint)){
    
    for(b in 1:length(p1)){
      basecases<-b-1# number of cases corresponding to probability vector p1
      totcases <-c-1# number of cases corresponding to probability vector p.joint
      outcases <- totcases-basecases
      o <- outcases +1  # number of cases corresponding to probability vector p2
      #ptotalcase[b,c]<-0
      if (o > 0){
        p.joint[c]<- p.joint[c] + p1[b]*p2[o]
      }   
    }
  }
  return(p.joint[1:nmax])
}

############ probability of n given outbreak given Trunk Geom baseline and NB outbreak
prob.out.TG.NB<- function(n,base.p, base.nonzero,out.mean,out.theta,full=FALSE,nmax=40){   
  n.count <- seq(1:nmax)-1
  p.base.temp <- zero.trunk.geom(n.count,base.p,base.nonzero )
  p.out.temp <- p.n(n.count,out.mean,out.theta,p.zero=0,distr="NB")
  p.tot.temp <- joint.prob(p.base.temp,p.out.temp)
  if(full==FALSE){
    return(c(p.base.temp[n+1],p.out.temp[n+1],p.tot.temp[n+1]))
  }else{
    return(list(p.base.temp,p.out.temp,p.tot.temp))
  }
}

############ probability of n given outbreak given Trunk Geom baseline and Geom outbreak
prob.out.TG.Geom<- function(n,base.p, base.nonzero,out.p,full=FALSE,nmax=40){
  
  n.count <- seq(1:nmax)-1
  p.base.temp <- zero.trunk.geom(n.count,base.p,base.nonzero )
  #p.out.temp <- p.n(n.count,out.mean,out.theta,p.zero=0,distr="NB")
  p.out.temp <- out.p*((1-out.p)^n.count)
  p.tot.temp <- joint.prob(p.base.temp,p.out.temp)
  if(full==FALSE){
    return(c(p.base.temp[n+1],p.out.temp[n+1],p.tot.temp[n+1]))
  }else{
    return(list(p.base.temp,p.out.temp,p.tot.temp))
  }
}
##

####
my.out.params.from.pop = function (pop=1256){ #for Geom
  my.clin <- pop*0.8*0.125*0.4
  my.resp <- my.clin*0.9*0.25
  my.neur <- my.clin*0.1*0.5
  presp <-1/(my.resp+1)
  pneur <-1/(my.neur+1)
  return(cbind(presp,pneur))
}

my.calculate.v =function(i,nmax,outbreak.ID="none"){ #nmax=upper limit Ncases
  #get data gor gridpoint
  tempout.X   <- list()  
  tempcount.X <- list()  
  for(j in 1:length(all.grid.counts.X)){
    tempcounts.X[[j]] <- all.grid.counts.X[[j]][which(all.grid.gridindex == i)]
    if(!outbreak.ID=="none"){
    #tempout.neur <- eval(parse(text=paste("my.gridinfo[[i]]$",outbreak.ID,"_neuro_counts",sep="")))
    #tempout.resp <- eval(parse(text=paste("my.gridinfo[[i]]$",outbreak.ID,"_respiratory_counts",sep="")))
    tempout.X[[j]] <- eval(parse(text=paste(paste("my.gridinfo[[i]]",outbreak.ID,"_decl_",projects[[j]],sep=""))))
    }
   }
  
  tempcounts.X[[j]] <-tempcounts.X[[j]]+tempout.X[[j]]

  temp.p.X <- list(); temp.nonzero.X <- list(); p.base.X <-list()
  for(j in 1:length(all.grid.counts.X)){
    temp.p.X[[j]] <- all.grid.geom.param.X[[j]][which(all.grid.gridindex == i)]
    temp.nonzero.X[[j]] <- all.grid.geom.nonzero[[X]][which(all.grid.gridindex == i)]
    p.base.X[[j]] <-  zero.trunk.geom(tempcounts.X[[j]],temp.p.X[[j]],temp.nonzero.X[[j]])
   }

  for(j in 1:length(all.grid.counts.X)){
   eval(parse(text = paste(paste("v.", projects[[j]],".i",sep=""), "<- c()")))
      tmp1[[j]] <- eval(parse(text = paste(paste("v.", projects[[j]],".i",sep=""), "<- c()")))
   v.tot.i <- c()
  #estimate V for each week in grid i
    nmaxr.X <- list(); probs.ij.X <-list()
    for(q in 1:length(tx)){
      nmaxr.X[[j]] <- min(nmax,max(tempcounts.X[[j]]))+2      
      if(outbreaktype == "Negbinom"){
        probs.ij.X[[j]] <-prob.out.TG.NB(tempcounts.X[[j]][q],all.grid.geom.param.X[[j]][which(all.grid.gridindex == i)][q],all.grid.geom.nonzero.X[[j]][which(all.grid.gridindex == i)][q],outbreak.mean.X[[j]],outbreak.theta.X[[j]],full=FALSE,nmax=nmaxr)
      }
    
      if(outbreaktype == "Geom.from.pop"){
        #prob.out.TG.Geom<- function(n,base.p, base.nonzero,out.p,full=FALSE,nmax=40){
        probs.ij.X[[j]] <- prob.out.TG.Geom(tempcounts.X[[j]][q],all.grid.geom.param.X[[j]][which(all.grid.gridindex == i)][q],all.grid.geom.nonzero.X[[j]][which(all.grid.gridindex == i)][q],my.gridinfo[[i]][["out.geom.p"]][1],full=FALSE,nmax=nmaxr)
      }
    
      tmp1[[j]][q]  <- probs.ij.X[[j]][3]/probs.ij.X[[j]][1]
    }
       
    v.tot.i <- sum(tmp1[[j]][q])
  }

  mylist<-list()
  for(j in 1:length(all.grid.counts.X)){
    mylist[[paste("v.", projects[[j]], sep="")]]<- log10(tmp1[[j]])
  }
  
  v.tot.i <- sum(mylist)
  mylist[["v.tot"]]<- log10(v.tot.i)
  
  return(mylist)
}

trunk.mean <-function(p){return(1/p)} #http://en.wikipedia.org/wiki/Geometric_distribution
