######MCMC VISUALIZATION######
######Sean Wu 4/27/2016#######
##############################


#####################################
###Visualize Trace of P-MCMC Chain###
#####################################

#load libraries
library(adaptivetau)
library(plyr)
library(reshape2)
library(ggplot2)
library(parallel)
library(foreach)
library(fitR)

#load data
flu_dat <- read.csv("FluTdCIncData.csv")
names(flu_dat) <- c("time","I")


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

#####################################
###Visualize Trace of P-MCMC Chain###
#####################################


mcmcMH <- function(target, init.theta, proposal.sd = NULL,
                   n.iterations, covmat = NULL,
                   limits=list(lower = NULL, upper = NULL),
                   adapt.size.start = NULL, adapt.size.cooling = 0.99,
                   adapt.shape.start = NULL,
                   print.info.every = n.iterations/100,
                   verbose = FALSE, max.scaling.sd = 50,
                   acceptance.rate.weight = NULL,
                   acceptance.window = NULL) {
  
  # initialise theta
  theta.current <- init.theta
  theta.propose <- init.theta
  
  # extract theta of gaussian proposal
  covmat.proposal <- covmat
  lower.proposal <- limits$lower
  upper.proposal <- limits$upper
  
  # reorder vector and matrix by names, set to default if necessary
  theta.names <- names(init.theta)
  if (!is.null(proposal.sd) && is.null(names(proposal.sd))) {
    names(proposal.sd) <- theta.names
  }
  
  if (is.null(covmat.proposal)) {
    if (is.null(proposal.sd)) {
      proposal.sd <- init.theta/5
    }
    covmat.proposal <-
      matrix(diag(proposal.sd[theta.names]^2, nrow = length(theta.names)),
             nrow = length(theta.names),
             dimnames = list(theta.names, theta.names))
  } else {
    covmat.proposal <- covmat.proposal[theta.names,theta.names]
  }
  
  if (is.null(lower.proposal)) {
    lower.proposal <- init.theta
    lower.proposal[] <- -Inf
  } else {
    lower.proposal <- lower.proposal[theta.names]
  }
  
  if (is.null(upper.proposal)) {
    upper.proposal <- init.theta
    upper.proposal[] <- Inf
  } else {
    upper.proposal <- upper.proposal[theta.names]
  }
  
  # covmat init
  covmat.proposal.init <- covmat.proposal
  
  adapting.size <- FALSE # will be set to TRUE once we start
  # adapting the size
  
  adapting.shape <- FALSE # will be set to TRUE once we start
  # adapting the shape
  
  # find estimated theta
  theta.estimated.names <- names(which(diag(covmat.proposal) > 0))
  
  # evaluate target at theta init
  target.theta.current <- target(theta.current)
  
  # if return value is a vector, set log.density and trace
  if (class(target.theta.current) == "numeric") {
    suppressWarnings(target.theta.current$log.density <-
                       target.theta.current)
    suppressWarnings(target.theta.current$trace <-
                       theta.current)
  }
  
  if (!is.null(print.info.every)) {
    message("Init: ", printNamedVector(theta.current[theta.estimated.names]),
            ", target: ", target.theta.current[["log.density"]])
  }
  
  trace <- data.frame(t(target.theta.current[["trace"]]))
  
  # acceptance rate
  acceptance.rate <- 0
  if (!is.null(acceptance.window)) {
    acceptance.series <- c()
  }
  
  # scaling factor for covmat size
  scaling.sd  <- 1
  
  # scaling multiplier
  scaling.multiplier <- 1
  
  # empirical covariance matrix (0 everywhere initially)
  covmat.empirical <- covmat.proposal
  covmat.empirical[,] <- 0
  
  # empirical mean vector
  theta.mean <- theta.current
  
  ########################
  ###VISUALIZATION CODE###
  ########################
  #set up plot area
  par(mfrow=c(3,1))
  
  #for R0
  plot(0,type="n",xlim=c(0,n.iterations),ylim=c(0,50),ylab="R0",xlab="")
  points(0,theta.current[[1]],pch=16)
  
  #for lat_dur
  plot(0,type="n",xlim=c(0,n.iterations),ylim=c(0,20),ylab="Latent Duration",xlab="")
  points(0,theta.current[[2]],pch=16)
  
  #for inf_dur
  plot(0,type="n",xlim=c(0,n.iterations),ylim=c(0,10),ylab="Infectious Duration",xlab="")
  points(0,theta.current[[3]],pch=16)
  
  # if print.info.every is null never print info
  if (is.null(print.info.every)) {
    print.info.every <- n.iterations + 1
  }
  
  start_iteration_time <- Sys.time()
  
  for (i.iteration in seq_len(n.iterations)) {
    
    # adaptive step
    if (!is.null(adapt.size.start) && i.iteration >= adapt.size.start &&
        (is.null(adapt.shape.start) || acceptance.rate*i.iteration < adapt.shape.start)) {
      if (!adapting.size) {
        message("\n---> Start adapting size of covariance matrix")
        adapting.size <- TRUE
      }
      # adapt size of covmat until we get enough accepted jumps
      scaling.multiplier <- exp(adapt.size.cooling^(i.iteration-adapt.size.start) * (acceptance.rate - 0.234))
      scaling.sd <- scaling.sd * scaling.multiplier
      scaling.sd <- min(c(scaling.sd,max.scaling.sd))
      # only scale if it doesn't reduce the covariance matrix to 0
      covmat.proposal.new <- scaling.sd^2*covmat.proposal.init
      if (!(any(diag(covmat.proposal.new)[theta.estimated.names] <
                .Machine$double.eps))) {
        covmat.proposal <- covmat.proposal.new
      }
      
    } else if (!is.null(adapt.shape.start) &&
               acceptance.rate*i.iteration >= adapt.shape.start) {
      if (!adapting.shape) {
        message("\n---> Start adapting shape of covariance matrix")
        # flush.console()
        adapting.shape <- TRUE
      }
      # adapt shape of covmat using optimal scaling factor for multivariate target distributions
      scaling.sd <- 2.38/sqrt(length(theta.estimated.names))
      
      covmat.proposal <- scaling.sd^2 * covmat.empirical
    }
    
    # print info
    if (i.iteration %% ceiling(print.info.every) == 0) {
      ## end_iteration_time <- Sys.time()
      state.mcmc <- trace[nrow(trace),]
      ## suppressMessages(time.estimation <- round(as.period((end_iteration_time-start_iteration_time)*10000/round(print.info.every))))
      ## message("Iteration: ",i.iteration,"/",n.iterations,", ETA: ",time.estimation,", acceptance rate: ",sprintf("%.3f",acceptance.rate),appendLF=FALSE)
      message("Iteration: ",i.iteration,"/", n.iterations,
              ", acceptance rate: ",
              sprintf("%.3f",acceptance.rate), appendLF=FALSE)
      if (!is.null(adapt.size.start) || !is.null(adapt.shape.start)) {
        message(", scaling.sd: ", sprintf("%.3f", scaling.sd),
                ", scaling.multiplier: ", sprintf("%.3f", scaling.multiplier),
                appendLF=FALSE)
      }
      message(", state: ",printNamedVector(state.mcmc))
      ## start_iteration_time <- end_iteration_time
    }
    
    # propose another parameter set
    if (any(diag(covmat.proposal)[theta.estimated.names] <
            .Machine$double.eps)) {
      print(covmat.proposal[theta.estimated.names,theta.estimated.names])
      stop("non-positive definite covmat",call.=FALSE)
    }
    if (length(theta.estimated.names) > 0) {
      theta.propose[theta.estimated.names] <-
        as.vector(rtmvnorm(1,
                           mean =
                             theta.current[theta.estimated.names],
                           sigma =
                             covmat.proposal[theta.estimated.names,theta.estimated.names],
                           lower =
                             lower.proposal[theta.estimated.names],
                           upper = upper.proposal[theta.estimated.names]))
    }
    
    # evaluate posterior of proposed parameter
    target.theta.propose <- target(theta.propose)
    # if return value is a vector, set log.density and trace
    if (class(target.theta.propose) == "numeric") {
      suppressWarnings(target.theta.propose$log.density <-
                         target.theta.propose)
      suppressWarnings(target.theta.propose$trace <-
                         theta.propose)
    }
    
    if (!is.finite(target.theta.propose$log.density)) {
      # if posterior is 0 then do not compute anything else and don't accept
      log.acceptance <- -Inf
      
    }else{
      
      # compute Metropolis-Hastings ratio (acceptance probability)
      log.acceptance <- target.theta.propose$log.density -
        target.theta.current$log.density
      log.acceptance <- log.acceptance +
        dtmvnorm(x = theta.current[theta.estimated.names],
                 mean =
                   theta.propose[theta.estimated.names],
                 sigma =
                   covmat.proposal[theta.estimated.names,
                                   theta.estimated.names],
                 lower =
                   lower.proposal[theta.estimated.names],
                 upper =
                   upper.proposal[theta.estimated.names],
                 log = TRUE)
      log.acceptance <- log.acceptance -
        dtmvnorm(x = theta.propose[theta.estimated.names],
                 mean = theta.current[theta.estimated.names],
                 sigma =
                   covmat.proposal[theta.estimated.names,
                                   theta.estimated.names],
                 lower =
                   lower.proposal[theta.estimated.names],
                 upper =
                   upper.proposal[theta.estimated.names],
                 log = TRUE)
      
    }
    
    if (verbose) {
      message("Propose: ", theta.propose[theta.estimated.names],
              ", target: ", target.theta.propose[["log.density"]],
              ", acc prob: ", exp(log.acceptance), ", ",
              appendLF = FALSE)
    }
    
    if (is.accepted <- (log(runif (1)) < log.acceptance)) {
      # accept proposed parameter set
      theta.current <- theta.propose
      target.theta.current <- target.theta.propose
      if (verbose) {
        message("accepted")
      }
    } else if (verbose) {
      message("rejected")
    }
    trace <- rbind(trace,c(target.theta.current$trace))
    
    # update acceptance rate
    if (i.iteration == 1) {
      acceptance.rate <- is.accepted
    } else {
      if (is.null(acceptance.rate.weight)) {
        if (is.null(acceptance.window)) {
          acceptance.rate <- acceptance.rate +
            (is.accepted - acceptance.rate) / i.iteration
        } else {
          acceptance.series <- c(is.accepted, acceptance.series)
          if (length(acceptance.series) > acceptance.window) {
            acceptance.series <- acceptance.series[-length(acceptance.series)]
          }
          acceptance.rate <- mean(acceptance.series)
        }
      } else {
        acceptance.rate <- acceptance.rate * (1 - acceptance.rate.weight) +
          is.accepted * acceptance.rate.weight
      }
    }
    
    # update empirical covariance matrix
    tmp <- updateCovmat(covmat.empirical, theta.mean,
                        theta.current, i.iteration)
    covmat.empirical <- tmp$covmat
    theta.mean <- tmp$theta.mean
    
    ########################
    ###VISUALIZATION CODE###
    ########################
    #R0
    par(mfg=c(1,1))
    points(i.iteration,theta.current[[1]],pch=16)
    
    #lat_dur
    par(mfg=c(2,1))
    points(i.iteration,theta.current[[2]],pch=16)
    
    #inf_dur
    par(mfg=c(3,1))
    points(i.iteration,theta.current[[3]],pch=16)
    
  }
  
  par(mfrow=c(1,1))
  
  return(list(trace = trace,
              acceptance.rate = acceptance.rate,
              covmat.empirical = covmat.empirical))
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

mcmc_out <- mcmcMH(target=posterior_seir,init.theta=init_theta_seir,n.iterations=100,limits=list(lower=low_lim))
