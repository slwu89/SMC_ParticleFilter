######fitting multiple MCMC chains in parallel######
######Sean Wu 4/26/2016#############################
####################################################


#load libraries
library(adaptivetau)
library(plyr)
library(reshape2)
library(ggplot2)
library(parallel)
library(foreach)
library(fitR)

#load data
setwd("~/Desktop/ParticleFilter")
flu_dat <- read.csv("SwineFlu2009DataVer2.csv")
names(flu_dat) <- c("time","I")
source("mcmc_visualization.R")

###Load stochastic model, particle filter###
############################################

#stochastic SIR simulation using ssa.adaptivetau
SEIR_stoch <- function(theta,init_state,times){
  
  #create transition matrix for Markov process
  trans <- list(
    c(S=-1,E=1), #susceptible to exposed
    c(E=-1,I=1), #exposed to infectious
    c(I=-1,R=1) #infectious to recovered
  )
  
  #create transition rates for Markov process at time t
  rates <- function(current_state,theta,t){
    
    #define parameters
    beta <- theta[["R0"]] / theta[["inf_dur"]]
    tau <- 1/theta[["lat_dur"]]
    rec <- 1/theta[["inf_dur"]]
    
    #define states
    S <- current_state[["S"]]
    E <- current_state[["E"]]
    I <- current_state[["I"]]
    R <- current_state[["R"]]
    N <- S + E + I + R
    
    #return rates
    rates <- c(
      beta * S * I/N, #susceptible to exposed
      tau * E, #exposed to infectious
      rec * I #infectious to recovered
    )
    
    return(rates)
  }
  
  #run the simulation
  mod_out <- as.data.frame(ssa.adaptivetau(init.values=init_state,transitions=trans,rateFunc=rates,params=theta,tf=diff(range(times))))
  
  #need to interpolate continuous time Markov process output to discrete time points
  mod_out$time <- mod_out$time + min(times)
  
  #interpolate output to integer values
  int_out <- apply(mod_out[,-1],2,function(i){
    approx(x=mod_out[,1],y=i,xout=times,method="constant")$y
  })
  
  #return simulation trajectory
  mod_traj <- cbind(time=times,int_out)
  return(as.data.frame(mod_traj))
}

#particle filter
particle_filter <- function(model,theta,init_state,data,n_particles,progress=TRUE){
  
  #initialize system
  loglike <- 0 #marginal log-likelihood of theta
  particle_current_state <- rep(list(init_state),n_particles) #current state of each particle
  particle_traj <- rep(list(data.frame(t(c(time=0,init_state)))),n_particles) #entire trajectory of each particle
  particle_weight <- rep(1/n_particles,length=n_particles) #weight of each particle
  
  #run particle filter for only times in the observed data
  time_points <- data$time
  
  #progress bar
  if(progress){
    progress_bar <- txtProgressBar(min=1,max=length(time_points))
  }
  
  #run the particle filter
  for(i in seq_along(time_points)){
    
    #run inner loop only over single timestep indexed by times
    times <- c(ifelse(i==1,0,time_points[i-1]),time_points[i]) #timestep to propagate particles over
    data_point <- unlist(data[i,]) #vector of the current observed data
    
    #resample particles
    na_weight <- which(is.na(particle_weight))
    particle_weight[na_weight] <- 0
    if(!all(particle_weight==0)){
      index_resampled <- sample(x=n_particles,size=n_particles,replace=TRUE,prob=particle_weight)
    }else{
      warning("All particles depleted at step ",i," of SMC. Return margLogLike = -Inf for theta: ",printNamedVector(theta), call.=FALSE)
      return(list(margLogLike=-Inf,traj=NA,traj.weight=NA))
      break
    }
    
    #update particles after resampling
    particle_traj <- particle_traj[index_resampled]
    particle_current_state <- particle_current_state[index_resampled]
    
    #propagate particles (function that works over particle_current_state)
    propagate <- foreach(current_state=particle_current_state,.packages=c("adaptivetau")) %dopar% {
      #run model over current time step
      traj <- model(theta=theta,init_state=unlist(current_state),times=times)
      
      #extract trajectories at latest time step
      model_point <- traj[nrow(traj),2:ncol(traj)]
      
      #compute particle weight
      weight <- dpois(x=data_point[["I"]],lambda=model_point[["E"]]) #incidence is modeled as a Poisson process
      
      return(list(state=model_point,weight=weight))
    }
    
    #collect parallel trajectories of each particle
    particle_current_state <- llply(propagate,function(x) {x$state})
    particle_weight <- unlist(llply(propagate,function(x) {x$weight}))
    particle_traj <- llply(seq_along(propagate),function(j) {rbind(particle_traj[[j]],c(data_point["time"],propagate[[j]]$state))})
    
    #update marginal log-likelihood of theta
    loglike <- loglike + log(mean(particle_weight))
    
    #advance progress bar
    if(progress){
      setTxtProgressBar(progress_bar,i)
    }
  }
  
  #end progress bar
  if(progress){
    close(progress_bar)
    print(paste("Particle Filter Complete! Marginal log-likelihood of theta is",round(loglike,6)))
  }
  
  #return marginal log-likelihood, filtered trajectories, weights of each trajectory
  output <- list(loglike=loglike,trajectory=particle_traj,traj_weight=particle_weight/sum(particle_weight))
  return(output)
}


###MCMC fitting of SEIR model to PA H1N1 outbreak data###
#########################################################

init_theta_seir <- c(R0=2,lat_dur=1.5,inf_dur=3)
init_state_seir <- c(S=369,E=0,I=1,R=0)
low_lim <- c(R0=1,lat_dur=0,inf_dur=0)

prior_serir <- function(theta){
  prior_R0 <- dunif(theta[["R0"]],min=1,max=50,log=TRUE)
  prior_lat_dur <- dunif(theta[["lat_dur"]],min=0,max=20,log=TRUE)
  prior_inf_dur <- log(dtruncnorm(theta[["inf_dur"]],a=0,b=Inf,mean=3,sd=1))
  ans <- sum(prior_R0,prior_lat_dur,prior_inf_dur)
  return(ans)
}

posterior_seir <- function(theta){
  log_prior <- prior_serir(theta=theta)
  log_likelihood <- particle_filter(model=SEIR_stoch,theta=theta,init_state=init_state_seir,data=flu_dat,n_particles=50,progress=TRUE)$loglike
  log_posterior <- log_prior + log_likelihood
  return(log_posterior)
}

#create parallel cluster and set different RNG seed for each node
cl <- parallel::makeCluster(spec=detectCores(),type="PSOCK")
parallel::clusterSetRNGStream(cl=cl,iseed=NULL)

#run 6 parallel chains
system.time(seir_mcmc <- mclapply(X=1:6,FUN=function(x){
  mcmcMH(target=posterior_seir,init.theta=init_theta_seir,n.iterations=5e4,limits=list(lower=low_lim))
},mc.preschedule=FALSE,mc.cores=6,mc.silent=FALSE))

#stop cluster
stopCluster(cl)
rm(cl)

#save output
saveRDS(object=seir_mcmc,file="mcmc_out.rds")