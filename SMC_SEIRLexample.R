######Particle Filter Example for Tristan da Cunha Influenza Outbreak######
######Sean Wu 4/22/2016####################################################
###########################################################################


###load packages and data###
############################

tdc_dat <- read.csv(file="FluTdCIncData.csv")
names(tdc_dat) <- c("time","I")

library(adaptivetau)
library(plyr)
library(reshape2)
library(ggplot2)
library(parallel)
library(foreach)
library(fitR)

#if you cannot install fitR, please uncomment and run the following
# install.packages("devtools")
# library(devtools)
# install_github("sbfnk/fitR")


###Stochastic SEIRL Model (Erlang-distributed waiting times)###
###############################################################

SEIRL_stoch <- function(theta,init_state,times){
  
  #create transition matrix for Markov process
  trans <- list(
    c(S=-1,E1=1), #susceptible to exposed 1
    c(E1=-1,E2=1), #exposed 1 to exposed 2
    c(E2=-1,I1=1),  #exposed 2 to infectious 1
    c(I1=-1,I2=1), #infectious 1 to infectious 2
    c(I2=-1,R1=1), #infectious 2 to recovered 1
    c(R1=-1,R2=1), #recovered 1 to recovered 2
    c(R2=-1,R3=1), #recovered 2 to recovered 3
    c(R3=-1,R4=1), #recovered 3 to recovered 4
    c(R4=-1), #recovered 4 to susceptible and immune
    c(S=1), #recovered 4 to susceptible
    c(L=1) #recovered 4 to immune
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
    E1 <- current_state[["E1"]]
    E2 <- current_state[["E2"]]
    I1 <- current_state[["I1"]]
    I2 <- current_state[["I2"]]
    R1 <- current_state[["R1"]]
    R2 <- current_state[["R2"]]
    R3 <- current_state[["R3"]]
    R4 <- current_state[["R4"]]
    L <- current_state[["L"]]
    I <- I1 + I2
    N <- S + E1 + E2 + I1 + I2 + R1 + R2 + R3 + R4 + L
    
    #return rates
    rates <- c(
      beta * S * I/N, #susceptible to exposed 1
      2 * f * E1, #exposed 1 to exposed 2
      2 * f * E2, #exposred 2 to infectious 1
      2 * rec * I1, #infectious 1 to infectious 2
      2 * rec * I2, #infectious 2 to recovered 1
      4 * tau * R1, #recovered 1 to recovered 2
      4 * tau * R2, #recovered 2 to recovered 3
      4 * tau * R3, #recovered 3 to recovered 4
      4 * tau * R4, #out of recovered 4
      4 * (1-alpha) * tau * R4, #recovered 4 to susceptible
      4 * alpha * tau * R4 #recovered 4 to immune
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
  
  #return simulation trajectory
  mod_traj <- cbind(time=times,int_out)
  return(as.data.frame(mod_traj))
}


###Particle Filter (SMC)###
###########################

#if you want to run the algorithm in parallel, please register a parallel backend prior
#to running the particle filter; otherwise it will run in serial
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
      weight <- dpois(x=data_point[["I"]],lambda=model_point[["E1"]]) #incidence is modeled as a Poisson process
      
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


###run SMC on Tristan da Cunha outbreak data with guessed parameter values###
#############################################################################

#initial parameters and state variables
tdc_init_state <- c(S=279,E1=0,E2=0,I1=2,I2=0,R1=3,R2=0,R3=0,R4=0,L=0)
tdc_theta_seirl <- c(R0=7,lat_dur=1,inf_dur=4,imm_dur=10,alpha=0.5)

#test the stochastic SEIRL model
SEIRL_stoch(theta=tdc_theta_seirl,init_state=tdc_init_state,1:59)

#register parallel backend
cl <- parallel::makeCluster(spec=detectCores(),type="PSOCK")
registerDoParallel(cl)

#run the particle filter (adjust the number of particles as needed)
system.time(smc_out <- particle_filter(model=SEIRL_stoch,
                  init_state=tdc_init_state,
                  theta=tdc_theta_seirl,
                  data=tdc_dat,
                  n_particles=250,
                  progress=TRUE))

#close parallel backend
stopCluster(cl)
rm(cl)


###function to visualize trajectories###
########################################

#type is a character indicating which compartment you wish to view a summarized time series for;
#data is the Tristan da Cunha data frame (tdc_dat); it will add observed data points to plot
plot_smc <- function(smc,data=NULL,type=NULL,colour="steelblue"){
  
  if(!is.character(type)) {
    stop("Please enter a character of the compartment you want to plot!")
  }
  
  #prepare data
  traj <- smc$trajectory
  names(traj) <- 1:length(traj)
  
  #extract compartment of interest
  traj_com <- ldply(traj,function(x){
    type_traj <- ddply(x,"time",function(x) x[type])
    names(type_traj)[ncol(type_traj)] <- "V1"
    type_traj <- join(x,type_traj,by="time")
    return(type_traj)
  },.id="replicate",.progress="text")
  
  #compute summary statistics
  traj_CI <- ddply(traj_com,c("time"),function(x){
    tmp_CI <- as.data.frame(t(quantile(x$V1,prob=c(0.025, 0.25, 0.5, 0.75, 0.975))))
    names(tmp_CI) <- c("low_95", "low_50", "median", "up_50", "up_95")
    tmp_CI$mean <- mean(x$V1)
    return(tmp_CI)
  })
  
  #prepare summary statistics for plotting
  traj_CI_line <- melt(traj_CI[c("time","mean","median")],id.vars="time")
  traj_CI_area <- melt(traj_CI[c("time","low_95","low_50","up_50","up_95")],id.vars="time")
  traj_CI_area$type <- sapply(traj_CI_area$variable, function(x) {str_split(x, "_")[[1]][1]})
  traj_CI_area$CI <- sapply(traj_CI_area$variable, function(x) {str_split(x, "_")[[1]][2]})
  traj_CI_area$variable <- NULL
  traj_CI_area <- dcast(traj_CI_area, paste0("time","+CI~type"))
  
  #plot summary statistics
  p <- ggplot(traj_CI_area)
  p <- p + geom_ribbon(data = traj_CI_area, aes_string(x = "time", ymin = "low", ymax = "up", alpha = "CI"), fill = colour)
  p <- p + geom_line(data = traj_CI_line, aes_string(x = "time", y = "value", linetype = "variable"), colour = colour)
  p <- p + scale_alpha_manual("Percentile", values = c("95" = 0.25, "50" = 0.45), labels = c("95" = "95th", "50" = "50th"))
  p <- p + scale_linetype("Stats")
  p <- p + guides(linetype = guide_legend(order = 1))
  
  #add data (if present)
  if(!is.null(data)){
    print("Adding data points")
    p <- p + geom_point(data=data,aes(x=time,y=I),colour="black")
  }
  
  #print to screen
  p <- p + theme_bw() + theme(legend.position = "top", legend.box = "horizontal")
  print(p)
}

#plot the filtered trajectories against data for compartment E1
plot_smc(smc_out,data=tdc_dat,type="E1",colour="steelblue")


###run Particle-MCMC to estimate theta###
#########################################

#warning: extremely time-consuming
prior_serirl <- function(theta){
  prior_R0 <- dunif(theta[["R0"]],min=1,max=50,log=TRUE)
  prior_lat_dur <- dunif(theta[["lat_dur"]],min=0,max=20,log=TRUE)
  prior_inf_dur <- log(dtruncnorm(theta[["inf_dur"]],a=0,b=Inf,mean=3,sd=1))
  prior_imm_dur <- dunif(theta[["imm_dur"]],min=0,max=20,log=TRUE)
  prior_alpha <- dunif(theta[["alpha"]],min=0,max=1,log=TRUE)
  ans <- sum(prior_R0,prior_lat_dur,prior_inf_dur,prior_imm_dur,prior_alpha)
  return(ans)
}

posterior_seirl <- function(theta){
  log_prior <- prior_serirl(theta=theta)
  log_likelihood <- particle_filter(model=SEIRL_stoch,theta=theta,init_state=tdc_init_state,data=tdc_dat,n_particles=20,progress=TRUE)$loglike
  log_posterior <- log_prior + log_likelihood
  return(log_posterior)
}

cl <- parallel::makeCluster(spec=detectCores(),type="PSOCK")
registerDoParallel(cl)

low_lim <- c(R0=1,lat_dur=0,inf_dur=0,imm_dur=0,alpha=0)
pmcmc_seirl <- mcmcMH(target=posterior_seirl,init.theta=tdc_theta_seirl,limits=list(lower=low_lim),n.iterations=1e3,adapt.size.start=100,adapt.shape.start=250,print.info.every=10,verbose=FALSE)

stopCluster(cl)
rm(cl)