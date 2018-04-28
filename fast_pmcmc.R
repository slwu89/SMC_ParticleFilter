######Fast Particle Filter in R######
######Sean Wu 3/14/2016##############
#####################################

library(adaptivetau)
library(plyr)
library(reshape2)
library(ggplot2)
library(parallel)
library(foreach)
library(doSNOW)
library(fitR)

###SIR Model###
###############

#stochastic SIR simulation using ssa.adaptivetau
SIR_stoch <- function(theta,init_state,times){
  
  #create transition matrix for Markov process
  trans <- list(
    c(S=-1,I=1), #infection
    c(I=-1,R=1) #recovery
  )
  
  #create transition rates for Markov process at time t
  rates <- function(current_state,theta,t){
    
    #define parameters
    beta <- theta[["R0"]] / theta[["inf_dur"]]
    gamma <- 1/theta[["inf_dur"]]
    
    #define states
    S <- current_state[["S"]]
    I <- current_state[["I"]]
    R <- current_state[["R"]]
    N <- S + I + R
    
    #return rates
    rates <- c(
      beta * S * I/N, #infection rate
      gamma * I #recovery rate
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
  
  mod_traj <- cbind(time=times,int_out)
  return(as.data.frame(mod_traj))
}

###SEIRL Model###
#################

SEIRL_stoch <- function(theta,init_state,times){
  
  #create transition matrix for Markov process
  trans <- list(
    c(S=-1,E=1), #exposed
    c(E=-1,I=1), #infectious
    c(I=-1,R=1), #recovery
    c(R=-1,L=1), #conferred immunity
    c(R=-1,S=1) #no conferred immunity
  )
  
  #create transition rates for Markov process at time t
  rates <- function(current_state,theta,t){
    
    #define parameters
    beta <- theta[["R0"]] / theta[["inf_dur"]]
    f <- 1 / theta[["lat_dur"]]
    rec <- 1 / theta[["inf_dur"]]
    tau <- 1 / theta[["lat_dur"]]
    alpha <- theta[["alpha"]]
    
    #define states
    S <- current_state[["S"]]
    E <- current_state[["E"]]
    I <- current_state[["I"]]
    R <- current_state[["R"]]
    L <- current_state[["L"]]
    N <- S + E + I + R + L
    
    #return rates
    rates <- c(
      beta * S * I/N, #rate of infection (latent)
      f * E, #rate of infection (infectious)
      rec * I, #rate of recovery
      alpha * tau * R, #rate of conferred immunity
      (1 - alpha) * tau * R #rate of reversion to susceptible
    )
  }
  
  #run the simulation
  mod_out <- as.data.frame(ssa.adaptivetau(init.values=init_state,transitions=trans,rateFunc=rates,params=theta,tf=diff(range(times))))
  
  #need to interpolate continuous time Markov process output to discrete time points
  mod_out$time <- mod_out$time + min(times)
  
  #interpolate output to integer values
  int_out <- apply(mod_out[,-1],2,function(i){
    approx(x=mod_out[,1],y=i,xout=times,method="constant")$y
  })
  
  mod_traj <- cbind(time=times,int_out)
  return(as.data.frame(mod_traj))
}

###SEIRWL Model###
##################
SEIRWL_stoch <- function(theta,init_state,times){
  
  #create transition matrix for Markov process
  trans <- list(
    c(S=-1), #exposure
    c(E=1), #exposure (with FOI calculated from S and W)
    c(E=-1,I=1), #infectious
    c(I=-1,R=1), #recovery
    c(R=-1,W=1), #window of re-infection
    c(W=-1,L=1), #conferred immunity
    c(W=-1,E=1) #re-infection
  )
  
  #create transition rates for Markov process at time t
  rates <- function(current_state,theta,t){
    
    #define parameters
    beta <- theta[["R0"]] / theta["inf_dur"]
    f <- 1 / theta[["lat_dur"]]
    rec <- 1 / theta[["inf_dur"]]
    gamma <- 1 / theta[["imm_dur"]]
    tau <- 1 / theta[["win_dur"]]
    
    #define states
    S <- current_state[["S"]]
    E <- current_state[["E"]]
    I <- current_state[["I"]]
    R <- current_state[["R"]]
    W <- current_state[["W"]]
    L <- current_state[["L"]]
    N <- S + E + I + R + W + L
    
    #return rates
    rates <- c(
      beta * S * I/N, #rate of exposure
      beta * (S+W) * I/N, #rate of exposure (with FOI calculated from S and W)
      f * E, #rate of infectiousness
      rec * I, #rate of recovery
      gamma * R, #rate of entry to window of re-infection
      tau * L, #rate of conferred immunity
      beta * W * I/N #rate of re-infection
    )
  }
  
  #run the simulation
  mod_out <- as.data.frame(ssa.adaptivetau(init.values=init_state,transitions=trans,rateFunc=rates,params=theta,tf=diff(range(times))))
  
  #need to interpolate continuous time Markov process output to discrete time points
  mod_out$time <- mod_out$time + min(times)
  
  #interpolate output to integer values
  int_out <- apply(mod_out[,-1],2,function(i){
    approx(x=mod_out[,1],y=i,xout=times,method="constant")$y
  })
  
  mod_traj <- cbind(time=times,int_out)
  return(as.data.frame(mod_traj))
}


###Particle Filter###
#####################

#parallel implementation of particle filter algorithm
particle_filter <- function(model,theta,init_state,data,n_particles,progress=TRUE,vis=FALSE){
  
  #redefine variables in function environment for parallel backend
  theta <- theta
  model <- model

  #initialize system
  loglike <- 0 #marginal log-likelihood of theta
  particle_current_state <- rep(list(init_state),n_particles) #current state of each particle
  particle_traj <- rep(list(data.frame(t(c(time=0,init_state)))),n_particles) #entire trajectory of each particle
  particle_weight <- rep(1/n_particles,length=n_particles) #weight of each particle
  
  #run particle filter for only times in the observed data
  time_points <- data$time
  
  #create list to store trajectories for visualization of particle filter
  if(vis){
    vis_dat <- list()
  }
  
  #progress bar
  if(progress){
    progress_bar <- txtProgressBar(min=1,max=length(time_points))
  }
  
  #run the particle filter
  for(i in seq_along(time_points)){
    
    #run inner loop only over single timestep indexed by times
    times <- c(ifelse(i==1,0,time_points[i-1]),time_points[i])
    data_point <- unlist(data[i,]) #vector of the current observed data
    
    #resample particles
    na_weight <- which(is.na(particle_weight))
    particle_weight[na_weight] <- 0
    if(!all(particle_weight==0)){
      index_resampled <- sample(x=n_particles,size=n_particles,replace=TRUE,prob=particle_weight)
    }else{
      warning("All particles depleted at step ",i," of SMC. Return margLogLike = -Inf for theta: ",printNamedVector(theta), call.=FALSE)
      return(list(margLogLike=-Inf,traj=NA,traj.weight=NA))
    }
    
    #particle filter visualization output
    if(vis){
      vis_input <- list(current=particle_current_state,weight=particle_weight,resamp=index_resampled)
      vis_dat[[i]] <- vis_input
    }
    
    #update particles after resampling
    particle_traj <- particle_traj[index_resampled]
    particle_current_state <- particle_current_state[index_resampled]
    
    #propagate particles (function that works over particle_current_state)
    propagate <- foreach(current_state=particle_current_state,.packages=c("adaptivetau")) %dopar% {
      #run model over current time step
      traj <- model(theta=theta,init_state=unlist(current_state),times=times)
      
      #extract trajectories at latest time step
      #model_point <- traj[nrow(traj),c("S","I","R")]
      model_point <- traj[nrow(traj),2:ncol(traj)]
      
      #compute particle weight
      weight <- dpois(x=data_point[["I"]],lambda=model_point[["I"]]) #incidence is modeled as a Poisson process
      
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
  if(vis){
    output <- list(loglike=loglike,trajectory=particle_traj,traj_weight=particle_weight/sum(particle_weight),vis_dat=vis_dat)  
  }
  output <- list(loglike=loglike,trajectory=particle_traj,traj_weight=particle_weight/sum(particle_weight))
  return(output)
}

#initialize variables
dir <- "/Users/slwu89/Desktop/git/StochasticInference/"
# prev_dat <- read.csv(paste0("SIRPrevData.csv"))

init_theta_sir <- c(R0=1.5,inf_dur=5)
init_theta_seirl <- c(R0=7,lat_dur=1,inf_dur=4,imm_dur=10,alpha=0.5)
init_theta_seirwl <- c(R0=5,lat_dur=1,inf_dur=4,imm_dur=10,win_dur=5)

init_state_tdc_sir <- c(S=279,I=2,R=3)
init_state_tdc_seirl <- c(S=279,E=0,I=2,R=3,L=0)
init_state_tdc_seirwl <- c(S=279,E=0,I=2,R=3,W=0,L=0)

SEIRL_stoch(init_theta_seirl,init_state_tdc_seirl,times=1:59)
SEIRWL_stoch(init_theta_seirwl,init_state_tdc_seirwl,times=1:59)

tdc_dat <- read.csv(paste0(dir,"FluTdCIncData.csv"))
names(tdc_dat) <- c("time","I")

#register parallel backend
nodes <- parallel::detectCores()
cl <- parallel::makePSOCKcluster(nodes)
doSNOW::registerDoSNOW(cl)

system.time(particle_filter(model=SIR_stoch,theta=init_theta_sir,init_state=init_state_tdc_sir,data=tdc_dat,n_particles=50,progress=TRUE,vis=FALSE))
system.time(particle_filter(model=SEIRL_stoch,theta=init_theta_seirl,init_state=init_state_tdc_seirl,data=tdc_dat,n_particles=50,progress=TRUE,vis=FALSE))
system.time(particle_filter(model=SEIRWL_stoch,theta=init_theta_seirwl,init_state=init_state_tdc_seirwl,data=tdc_dat,n_particles=50,progress=TRUE,vis=FALSE))

stopCluster(cl)
rm(cl)

###Particle MCMC###
###################

#prior distribution functions
prior_seirwl <- function(theta) {
  prior_R0 <- dunif(theta[["R0"]],min=1,max=20,log=TRUE)
  prior_lat_dur <- dunif(theta[["lat_dur"]],min=0,max=7,log=TRUE)
  prior_inf_dur <- dunif(theta[["inf_dur"]],min=0,max=7,log=TRUE)
  prior_imm_dur <- dunif(theta[["imm_dur"]],min=0,max=14,log=TRUE)
  prior_win_dur <- dunif(theta[["win_dur"]],min=0,max=14,log=TRUE)
  ans <- sum(prior_R0,prior_lat_dur,prior_inf_dur,prior_imm_dur,prior_win_dur)
  return(ans)
}

prior_seirl <- function(theta) {
  prior_R0 <- dunif(theta[["R0"]],min=1,max=20,log=TRUE)
  prior_lat_dur <- dunif(theta[["lat_dur"]],min=0,max=7,log=TRUE)
  prior_inf_dur <- dunif(theta[["inf_dur"]],min=0,max=7,log=TRUE)
  prior_imm_dur <- dunif(theta[["imm_dur"]],min=0,max=14,log=TRUE)
  prior_alpha <- dunif(theta[["alpha"]],min=0,max=1,log=TRUE)
  ans <- sum(prior_R0,prior_lat_dur,prior_inf_dur,prior_imm_dur,prior_alpha)
  return(ans)
}

prior_sir <- function(theta){
  prior_R0 <- dunif(theta[["R0"]],min=1,max=40,log=TRUE)
  prior_inf_dur <- dunif(theta[["inf_dur"]],min=0,max=30,log=TRUE)
  ans <- sum(prior_R0,prior_inf_dur)
  return(ans)
}


#posterior target distribution
posterior_smc <- function(theta){
  log_prior <- prior_seirwl(theta=theta)
  log_likelihood <- particle_filter(model=SEIRWL_stoch,theta=theta,init_state=init_state_tdc_seirwl,data=tdc_dat,n_particles=20,progress=TRUE,vis=FALSE)$loglike
  log_posterior <- log_prior + log_likelihood
  return(log_posterior)
}

posterior_sir <- function(theta){
  log_prior <- prior_sir(theta=theta)
  log_likelihood <- particle_filter(model=SIR_stoch,theta=theta,init_state=init_state_sir,data=prev_dat,n_particles=20,progress=TRUE,vis=FALSE)$loglike
  log_posterior <- log_prior + log_likelihood
  return(log_posterior)
}

#Adaptive Metropolis-Hastings P-MCMC Sampling
low_lim <- c(R0=1,lat_dur=0,inf_dur=0,imm_dur=0,win_dur=0)
up_lim <- c(R0=Inf,inf_dur=Inf)

pmcmc <- mcmcMH(target=posterior_smc,init.theta=init_theta_seirwl,n.iterations=5000,limits=list(lower=low_lim),adapt.size.start=500,adapt.shape.start=800,print.info.every=10,verbose=TRUE)

#sir version
low_lim <- c(R0=1,inf_dur=0)
up_lim <- c(R0=Inf,inf_dur=Inf)
pmcmc_sir <- mcmcMH(target=posterior_sir,init.theta=init_theta_sir,n.iterations=1e3,limits=list(lower=low_lim),adapt.size.start=200,adapt.shape.start=500,print.info.every=10,verbose=TRUE)


###Particle Filter Visualization###
###################################

#extract data from particle filter output
pf_vis <- particle_filter(model=SIR_stoch,theta=init_theta_sir,init_state=init_state_sir,data=swine_dat,n_particles=10,progress=TRUE)
pf_vis <- pf_vis$vis_dat

#extract incidence trajectories of each particle from results
vis_traj <- t(sapply(pf_vis,function(x) unlist(x$current)))
vis_traj <- vis_traj[,which(colnames(vis_traj) == "I")]
#extract resampling weights of each particle from results
vis_weight <- t(sapply(pf_vis,function(x) x$weight))
vis_weight <- t(apply(vis_weight,1,function(x) x / sum(x)))
#extract resampling indicies of each particle from results
vis_resamp <- t(sapply(pf_vis,function(x) x$resamp))

#create data frame for resampled trajectories (geom_line)
vis_traj_resamp <- matrix(NA,nrow=37,ncol=10)
for(i in 1:37){
  vis_traj_resamp[i,] <- vis_traj[i,vis_resamp[i,]]
}

#create data frame for weighted particles (geom_point)
vis_traj <- melt(vis_traj)
vis_traj$weight <- melt(vis_weight)$value
names(vis_traj)[4] <- "Weight"

ggplot() +
  geom_point(data=vis_traj,aes(x=Var1,y=value,size=Weight),alpha=0.5,colour="darkblue") +
  geom_line(data=melt(vis_traj_resamp),aes(x=Var1,y=value,group=Var2),colour="darkblue") +
  geom_line(aes(x=time,y=I),colour="tomato",data=SIRPrevData) +
  labs(x="Time",y="Incidence",title="Particle Filter Visualization") +
  #guides(size=FALSE) +
  theme_bw() +
  theme(legend.position="bottom")