#=====================
# User input
#=====================

#capture history
ch <- data.frame(
  ANIMAL_ID=c(1,2,3,4,5,6,7,7,8,9,10,11,12,13,14,5,15,7,16,2,17,18,5,16,14,19,20,21,22,2,14,4,3,23,21,11,24,25,26,27,28,29,16,18,30,24,1,26,31,32,33,34,3,2,10,2,35,33,36,37,28,35,23,31,3,4,34,14,38,6,39,11,40,41,42,43,44,45,46,17,25,14,14,47,48,23,49,50,35,51,52),
  LOC_ID=c(12,26,29,29,30,50,54,55,58,59,74,78,78,84,15,39,49,53,56,74,76,77,37,43,44,44,53,56,62,64,25,29,30,45,64,78,79,17,25,29,35,54,56,67,69,70,12,12,13,17,18,27,30,64,64,65,82,19,43,46,54,82,5,12,19,19,26,26,27,48,53,78,87,33,48,54,58,68,69,76,8,15,25,33,44,44,48,48,82,88,99),
  DAY=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10)
)

#trap coordinates
traps <- data.frame(
  x=c(0,4.44,8.89,13.33,17.78,22.22,26.67,31.11,35.56,40,0,4.44,8.89,13.33,17.78,22.22,26.67,31.11,35.56,40,0,4.44,8.89,13.33,17.78,22.22,26.67,31.11,35.56,40,0,4.44,8.89,13.33,17.78,22.22,26.67,31.11,35.56,40,0,4.44,8.89,13.33,17.78,22.22,26.67,31.11,35.56,40,0,4.44,8.89,13.33,17.78,22.22,26.67,31.11,35.56,40,0,4.44,8.89,13.33,17.78,22.22,26.67,31.11,35.56,40,0,4.44,8.89,13.33,17.78,22.22,26.67,31.11,35.56,40,0,4.44,8.89,13.33,17.78,22.22,26.67,31.11,35.56,40,0,4.44,8.89,13.33,17.78,22.22,26.67,31.11,35.56,40),
  y=c(0,0,0,0,0,0,0,0,0,0,4.44,4.44,4.44,4.44,4.44,4.44,4.44,4.44,4.44,4.44,8.89,8.89,8.89,8.89,8.89,8.89,8.89,8.89,8.89,8.89,13.33,13.33,13.33,13.33,13.33,13.33,13.33,13.33,13.33,13.33,17.78,17.78,17.78,17.78,17.78,17.78,17.78,17.78,17.78,17.78,22.22,22.22,22.22,22.22,22.22,22.22,22.22,22.22,22.22,22.22,26.67,26.67,26.67,26.67,26.67,26.67,26.67,26.67,26.67,26.67,31.11,31.11,31.11,31.11,31.11,31.11,31.11,31.11,31.11,31.11,35.56,35.56,35.56,35.56,35.56,35.56,35.56,35.56,35.56,35.56,40,40,40,40,40,40,40,40,40,40) 
)

#number of trapping occasions (days)
nDays <- 10

#buffer zone around trapping area
#(in terms of percentage of trapping area)
bufferZonePercent <- 0.10

#------------------------
# MCMC parameters
#------------------------
chainLength <- 5000
burnin <- 3000

#data augmentation upperbound
M <- 120

#parameter scale values for MCMC
sigma_scale <- .5
lam_scale <- 0.05
sigma_pi_scale <- 1
psi_alpha_scale <- .1
h_scale <- 1

#parameter inits for MCMC
sigma <- runif(1,1,5)
lam <- runif(1,.05,.3)
psi <- runif(1,.2, 0.8)
sigma_pi <- runif(1,sigma,10)
psi_alpha <- runif(1,0.5,.9)


#==================
# define functions
#==================

#categorical sampler
rcat <- function(n,p){
  return(sample(1:length(p),size=n,replace=TRUE,prob=p))
}

#distance function
distance <- function(x1,y1,x2,y2){
  return(sqrt((x1-x2)^2+(y1-y2)^2))
}

#gets index for given tiger id, trap and day
getIndex <- function(id,trap,day){
  return((id-1)*nTraps*nDays+(trap-1)*nDays+day)
}

#========================
# prepare data for MCMC
#========================

#number of individuals captures
n <- length(unique(ch$ANIMAL_ID))

#number of traps
nTraps <- nrow(traps)

#add trap coordinates
ch$X_COORD <- traps$x[ch$LOC_ID]
ch$Y_COORD <- traps$y[ch$LOC_ID]

#change the captured animal's IDs to 1:n
ch$ANIMAL_ID <- match(ch$ANIMAL_ID,unique(ch$ANIMAL_ID))
observed_IDs <- 1:n
unobserved_IDs <- (n+1):M

#minimum rectangle area containing traps
areaWidth <- max(traps$x)-min(traps$x)
areaHeight <- max(traps$y)-min(traps$y)

#assign trap location ID's
traps$LOC_ID <- 1:nrow(traps)

#width and height of buffer zone around trapping area
bufferWidth <- areaWidth*bufferZonePercent
bufferHeight <- areaHeight*bufferZonePercent

#borders of the trapping area (including buffer zone)
leftBorder <- min(traps$x) - bufferWidth
rightBorder <- leftBorder + areaWidth + bufferWidth
bottomBorder <- min(traps$y) - bufferHeight
topBorder <- bottomBorder + areaHeight + bufferHeight

#format data for MCMC for all 'M' tigers
trapCaps <- data.frame(ANIMAL_ID=rep(1:M,each=nDays*nTraps),
                       LOC_ID=rep(rep(1:nTraps,each=nDays),times=M),
                       DAY=rep(1:nDays,times=nTraps*M),
                       cap=rep(0,nTraps*M*nDays))

#fill caps column
for(i in 1:nrow(ch)){
  id.i <- ch$ANIMAL_ID[i]
  trap.i <- ch$LOC_ID[i]
  day.i <- ch$DAY[i]
  trapCaps$cap[getIndex(id.i,trap.i,day.i)] <- 1
}

#observation vectors
y <- trapCaps$cap
y.n <- trapCaps$cap[which(trapCaps$ANIMAL_ID<=n)]




#=============
# MCMC
#=============

#lowerbound and upperbounds for sigma_pi
minCutoff <- 1e-300
sigma_pi_min <- sqrt(min(areaWidth,areaHeight)^2/(2*log(1/minCutoff)))
sigma_pi_max <- sqrt((areaWidth)^2+(areaHeight)^2)

#unknown individual's home range centers inits
unobs_hx_inits <- runif(M-n,leftBorder,rightBorder)
unobs_hy_inits <- runif(M-n,bottomBorder,topBorder)

#observed individual's home range centers inits
obs_hx_est <- rep(NA,n)
obs_hy_est <- rep(NA,n)
for(i in observed_IDs){
  tigerCaps.i <- ch[which(ch$ANIMAL_ID==i),]
  obs_hx_est[i] <- mean(traps$x[tigerCaps.i$LOC_ID])
  obs_hy_est[i] <- mean(traps$y[tigerCaps.i$LOC_ID])
}
obs_hx_inits <- obs_hx_est
obs_hy_inits <- obs_hy_est

#all M individual's home range centers inits and home range center distance init
hx <- c(obs_hx_inits,unobs_hx_inits)
hy <- c(obs_hy_inits,unobs_hy_inits)
H_dist2 <- as.matrix(dist(data.frame(hx,hy),upper=T,diag=T))^2

#w inits
w <- c(rep(1,n),rbinom(M-n, 1, psi*(1-(n/M))))


#pi matrix init
Pi_num <- exp(-H_dist2/(2*sigma_pi^2))
diag(Pi_num) <- 0
Pi_num_w <- t(t(Pi_num)*w)
Pi <- Pi_num_w/rowSums(Pi_num_w)

#attraction inits
Attraction <- matrix(NA,nrow=n,ncol=nDays)
for(i in 1:n){
  Attraction[i,] <- rcat(nDays,Pi[i,])
}

#alpha matrix init
Alpha <- matrix(rbinom(n*nDays,1,psi_alpha),nrow=n,ncol=nDays)

#initialize lla (log-like for attractions)
llaRowIndex <- rep(c(1:nrow(Attraction)),each=nDays)
lla <- 0
for(i in 1:n){
  for(t in 1:nDays){
    if(Alpha[i,t]==0){
      lla <- lla + log(Pi[i,Attraction[i,t]])
    }
  }
}

#occasion level activity center for observed tigers inits
S <- matrix(NA,nrow=n,ncol=nDays)
for(i in 1:n){
  S[i,] <- Alpha[i,]*i + (1-Alpha[i,])*Attraction[i,]
}

#squared distance between traps and home range centers
unobservedS <- matrix(rep((n+1):M,each=nDays),nrow=M-n,byrow=T)
dist2 <- (rep(hx[as.vector(t(rbind(S,unobservedS)))],
              each=nTraps)-traps$x[trapCaps$LOC_ID])^2+
         (rep(hy[as.vector(t(rbind(S,unobservedS)))],
              each=nTraps)-traps$y[trapCaps$LOC_ID])^2

#capture probabilty
mu <- 1-exp(-lam*exp(-dist2/(2*sigma^2)))
mu[mu < 1e-30] <- 1e-30

#indexes which w={1,0}
w1Index <- which(trapCaps$ANIMAL_ID%in%which(w==1))
w0Index <- which(trapCaps$ANIMAL_ID%in%which(w==0))

#log-like for observations
lly.vec <- dbinom(y, 1, mu*rep(w,each=nTraps*nDays), log = TRUE)
lly <- sum(lly.vec)

#indexes for individuals and days
#used to call subsets of lly and y
#quicker than calculating during the MCMC
index.i.mat <- matrix(NA,nrow=M,ncol=nDays*nTraps)
index.it.list <- replicate(M, matrix(NA,nrow=nDays,ncol=nTraps), simplify=F)
y.list <- replicate(M, matrix(NA,nrow=nDays,ncol=nTraps), simplify=F)
for(i in 1:M){
  index.i.mat[i,] <- which(trapCaps$ANIMAL_ID==i)
  for(t in 1:nDays){
    index.it.list[[i]][t,] <- getIndex(i,1:nTraps,t)
    y.list[[i]][t,] <- y[getIndex(i,1:nTraps,t)]
  }
}


#====================
# Run MCMC
#====================
chain <- matrix(nrow = chainLength, ncol = 6)
colnames(chain) <- c("sigma", "lam", "psi", "N",'sigma_pi','psi_alpha')
chain[1,] <- c(sigma,lam,psi,sum(w),sigma_pi,psi_alpha)
startTime <- proc.time()
for (iter in 2:chainLength) {
  if (iter%%100 == 0) {
    print(paste0('iteration=',iter,
                 ' elapsed_time=',(proc.time()-startTime)[3],
                 ' N=',sum(w),
                 ' sigma=',round(sigma,2),
                 ' lam=',round(lam,3),
                 ' sigma_pi=',round(sigma_pi,2),
                 ' psi_alpha=',round(psi_alpha,2)))
      par(mfrow=c(4,2))
      plot(chain[,4],type='l',main='N_est traceplot')
      hist(chain[,4],main='N_est  histogram')
      plot(chain[,1],type='l',main='sigma_est traceplot')
      hist(chain[,1],main='sigma_est  histogram')
      plot(chain[,2],type='l',main='lambda_est traceplot')
      hist(chain[,2],main='lambda_est  histogram')
      plot(chain[,5],type='l',main='sigmaPi_est traceplot')
      hist(chain[,5],main='sigmaPi_est  histogram')
      par(mfrow=c(1,1))
  }
  
  #-----------------
  # sigma
  #-----------------
  sigma_prop <- rnorm(1, sigma, sigma_scale)
  if (sigma_prop > 0) {

    #proposal mu
    mu.w1_prop <- 1-exp(-lam*exp(-dist2[w1Index]/(2*sigma_prop^2)))
    
    #log-like's for data (y)
    lly.vec_prop <- dbinom(y[w1Index],1,mu.w1_prop, log = TRUE)
    lly_prop <- sum(lly.vec_prop)
    
    #M-H step
    if (runif(1) < exp(lly_prop - lly)) {
      mu[w1Index] <- mu.w1_prop
      lly.vec[w1Index] <- lly.vec_prop
      lly <- lly_prop
      sigma <- sigma_prop
    }
  }

  #----------------------
  # lambda
  #----------------------
  lam_prop <- rnorm(1, lam, lam_scale)
  if (lam_prop > 0) {
    
    #proposal mu
    mu.w1_prop <- 1-exp(-lam_prop*exp(-dist2[w1Index]/(2*sigma^2)))
    
    #log-like's for the data (y)
    lly.vec_prop <- dbinom(y[w1Index],1,mu.w1_prop,log=T)
    lly_prop <- sum(lly.vec_prop)
    
    #M-H step
    if (runif(1) < exp(lly_prop - lly)) {
      mu[w1Index] <- mu.w1_prop
      lly.vec[w1Index] <- lly.vec_prop
      lly <- lly_prop
      lam <- lam_prop
    }
  }
  
  #-----------
  # sigma_pi
  #-----------
  sigma_pi_prop <- rnorm(1, sigma_pi, sigma_pi_scale)
  if (sigma_pi_prop>sigma_pi_min & sigma_pi_prop<sigma_pi_max) {
    
    #proposal pi matrix
    Pi_num_prop <- exp(-H_dist2/(2*sigma_pi_prop^2))
    diag(Pi_num_prop) <- 0
    Pi_num_w_prop <- t(t(Pi_num_prop)*w)
    Pi_prop <- Pi_num_w_prop/rowSums(Pi_num_w_prop)
    
    #log-like's for attractions
    lla_prop <- 0
    for(i in 1:n){
      for(t in 1:nDays){
        if(Alpha[i,t]==0){
          lla_prop <- lla_prop + log(Pi_prop[i,Attraction[i,t]])
        }
      }
    }

    #M-H step
    if (runif(1) < exp(lla_prop-lla)) {
      sigma_pi <- sigma_pi_prop
      Pi_num <- Pi_num_prop
      Pi_num_w <- Pi_num_w_prop 
      Pi <- Pi_prop
      lla <- lla_prop
    }
  }
  
  #---------------------------
  # unobserved w's
  #---------------------------
  for (i in (n+1):M) {
    w_prop.i <- ifelse(w[i] == 0, 1, 0)
    
    #log-like's for w_i's prior
    llPrior <- dbinom(w[i], 1, psi, log = TRUE)
    llPrior_prop <- dbinom(w_prop.i, 1, psi, log = TRUE)
    
    #proposed pi matrix
    w_prop <- w #vector of w's including proposed w_i 
    w_prop[i] <- w_prop.i
    Pi_num_w_prop <- Pi_num_w
    Pi_num_w_prop[,i] <- Pi_num[,i]*w_prop.i
    Pi_prop <- Pi_num_w_prop/rowSums(Pi_num_w_prop)
    
    #log-like's for attractions
    lla_prop <- 0
    for(j in 1:n){
      for(t in 1:nDays){
        if(Alpha[j,t]==0){
          lla_prop <- lla_prop + log(Pi_prop[j,Attraction[j,t]])
        }
      }
    }
    
    #indexes for i^th individual
    index.i <- index.i.mat[i,]
    
    #log-like's for data (y)
    if(w_prop.i==1){
      mu.i_prop <- 1-exp(-lam*exp(-dist2[index.i]/(2*sigma^2)))
    }else{
      #if w=0 don't calculate mu
      mu.i_prop <- 0
    }
    #observation log-like for i^th individual
    lly.i.vec_prop <- log(1-mu.i_prop)
    lly.i_prop <- sum(lly.i.vec_prop)
    lly.i <- sum(lly.vec[index.i])

    #M-H step
    if (runif(1) < exp((lly.i_prop + llPrior_prop + lla_prop) - 
                       (lly.i + llPrior + lla))) {
      w[i] <- w_prop.i
      Pi_num_w <- Pi_num_w_prop 
      Pi <- Pi_prop
      lla <- lla_prop
      mu[index.i] <- mu.i_prop
      lly.vec[index.i] <- lly.i.vec_prop
    }
  }
  lly <- sum(lly.vec)
  w1Index <- which(trapCaps$ANIMAL_ID%in%which(w==1))
  w0Index <- which(trapCaps$ANIMAL_ID%in%which(w==0))
  
  #------------
  # psi
  #------------
  psi <- rbeta(1, 1 + sum(w), 1 + M - sum(w))
  
  #------------------------------
  # observed home range centers
  #------------------------------
  for (i in 1:n) {
    hx_prop.i <- rnorm(1, hx[i], h_scale)
    hy_prop.i <- rnorm(1, hy[i], h_scale)
    if(hx_prop.i>leftBorder & hx_prop.i<rightBorder & 
       hy_prop.i>bottomBorder & hy_prop.i<topBorder) {
      
      #distances to traps from proposal home range centers
      dist2_prop.i1 <- (hx_prop.i-traps$x)^2+
                       (hy_prop.i-traps$y)^2
      mu_prop.i1 <- 1-exp(-lam*exp(-dist2_prop.i1/(2*sigma^2)))

      #vector of home range centers including the proposals 
      #sx_prop.i and sy_prop.i
      hx_prop<-hx[1:M]; hx_prop[i]<-hx_prop.i
      hy_prop<-hy[1:M]; hy_prop[i]<-hy_prop.i
      
      #calculate proposal pi matrix
      H_dist2_prop.i <- (hx_prop.i-hx_prop)^2+(hy_prop.i-hy_prop)^2
      Pi_num_w_prop <- Pi_num_w
      Pi_num_w_prop[,i] <- exp(-H_dist2_prop.i/(2*sigma_pi^2))
      Pi_num_w_prop[i,] <- Pi_num_w_prop[,i]*w
      Pi_num_w_prop[i,i] <- 0
      Pi_prop <- Pi_num_w_prop/rowSums(Pi_num_w_prop)
      
      #log-like's for attractions
      lla_prop <- 0
      attracted.indexes <- c()
      lly.i.vec_prop <- c()
      dist2_prop.i <- c()
      for(j in 1:n){
        for(t in 1:nDays){
          if(Alpha[j,t]==0){
            lla_prop <- lla_prop + log(Pi_prop[j,Attraction[j,t]])
          }

          #keep track of indexes associated with individuals attracted
          #to the i^th individual
          if(S[j,t]==i){
            attracted.indexes <- c(attracted.indexes,index.it.list[[j]][t,])
            dist2_prop.i <- c(dist2_prop.i,dist2_prop.i1)
          }
        }
      }
      
      #calculate observation log-like for individuals attracted
      #to the i^th individual
      if(length(attracted.indexes)>0){
        mu_prop.i <- 1-exp(-lam*exp(-dist2_prop.i/(2*sigma^2)))
        lly.i.vec_prop <- dbinom(y[attracted.indexes],1,mu_prop.i,log=TRUE)
        lly.i_prop <- sum(lly.i.vec_prop)
        lly.i <- sum(lly.vec[attracted.indexes])
      }else{
        
        #if there are no attractions to the i^th individual
        #this will remove observation log-like from M-H ratio
        lly.i_prop <- 0
        lly.i <- 0
      }
      
      #M-H step
      if (runif(1) < exp(lla_prop+lly.i_prop-
                         (lla+lly.i))) {
        hx[i] <- hx_prop.i
        hy[i] <- hy_prop.i
        
        #if there were attractions to tiger i
        #update mu, dist2 and lly, otherwise 
        #they shouldn't be updated now, they will
        #be updated when the home range centers they're 
        #attracted to are updated
        if(length(attracted.indexes)>0){
          mu[attracted.indexes] <- mu_prop.i
          dist2[attracted.indexes] <- dist2_prop.i
          lly.vec[attracted.indexes] <- lly.i.vec_prop
        }
        
        H_dist2[i,] <- H_dist2_prop.i
        H_dist2[,i] <- H_dist2_prop.i
        Pi_num[i,] <- Pi_num_w_prop[,i]
        Pi_num[,i] <- Pi_num_w_prop[,i]
        Pi_num_w <- Pi_num_w_prop
        Pi <- Pi_prop
        lla <- lla_prop
      }
    }
  }
  lly <- sum(lly.vec)
  
  
  #----------------------------------
  # Alpha matrix: attraction boolean
  #----------------------------------
  for(i in 1:n){
    for(t in 1:nDays){
      
      #prop for alpha for individual during a single day
      alpha_prop.it <- ifelse(Alpha[i,t]==1,0,1)
      
      #log-like's for alpha.it's prior
      llAlpha_prior <- dbinom(Alpha[i,t],1,psi_alpha,
                              log=TRUE)
      llAlpha_prior_prop <- dbinom(alpha_prop.it,1,psi_alpha,
                                   log=TRUE)
      
      #propose attraction if i^th individual
      #is away from home range center on t^th day
      #otherwise set activity center to home range center
      if(alpha_prop.it==0){
        attraction_prop.it <- rcat(1,Pi[i,])
        S_prop.it <- attraction_prop.it
      }else{
        S_prop.it <- i
      }
        
      #squared distances for tiger i on day t
      #for newly proposed home range center
      dist2_prop.it <- (hx[S_prop.it]-traps$x)^2 + 
                       (hy[S_prop.it]-traps$y)^2
      
      #proposed mu for tiger i on day t
      mu.it_prop <- 1 - exp(-lam*exp(-dist2_prop.it/(2*sigma^2)))
      
      #likelihood for observation y for tiger i on day t
      index.it <- index.it.list[[i]][t,]
      lly.it.vec_prop <- dbinom(y.list[[i]][t,],1,mu.it_prop,log=TRUE)
      lly.it_prop <- sum(lly.it.vec_prop)
      lly.it <- sum(lly.vec[index.it])

      #M-H step
      if (runif(1) < exp((lly.it_prop + llAlpha_prior_prop)- 
                         (lly.it + llAlpha_prior))) {
        Alpha[i,t] <- alpha_prop.it
        S[i,t] <- S_prop.it
        dist2[index.it] <- dist2_prop.it
        mu[index.it] <- mu.it_prop
        lly.vec[index.it] <- lly.it.vec_prop
        
        #update the attraction if alpha.it==0,
        #otherwise u.it=i and this attraction
        #won't be referenced
        if(Alpha[i,t]==0){
          Attraction[i,t] <- attraction_prop.it
        }
      }
    }
  }
  #update lly and lla now that all attractions have
  #been proposed
  lly <- sum(lly.vec)
  lla <- 0
  for(i in 1:n){
    for(t in 1:nDays){
      if(Alpha[i,t]==0){
        lla <- lla + log(Pi[i,Attraction[i,t]])
      }
    }
  }
  
  #------------
  # psi_alpha
  #------------
  sumAlpha <- sum(Alpha)
  psi_alpha <- rbeta(1, 1 + sumAlpha, 1 + n*nDays - sumAlpha)
  
  #update chain
  chain[iter, ] <- c(sigma, lam, psi, sum(w),sigma_pi,psi_alpha)
}


#==========
# results
#==========
quant <- apply(chain,2,quantile,probs=c(0.025,0.5,0.975))
estimates <- cbind(as.matrix(apply(X=chain[burnin:chainLength,],FUN=mean,MARGIN=2)),as.matrix(apply(X=chain[burnin:chainLength,],FUN=sd,MARGIN=2)),round(t(quant),3))
colnames(estimates) <- c('mean','sd',"2.5%","median","97.5%")
estimates
