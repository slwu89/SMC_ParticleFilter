######Additional Stochastic Compartmental Model Examples######
######Sean Wu 4/23/2016#######################################
##############################################################

library(adaptivetau)

###Stochastic SIR Model###
##########################

#stochastic SIR simulation using ssa.adaptivetau
SIR_stoch <- function(theta,init_state,times){
  
  #create transition matrix for Markov process
  trans <- list(
    c(S=-1,I=1), #susceptible to infectious
    c(I=-1,R=1) #infectious to recovered
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
      beta * S * I/N, #susceptible to infectious
      gamma * I #infectious to recovered
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

#initial parameters and state variables
init_theta_sir <- c(R0=1.5,inf_dur=5)
init_state_sir <- c(S=999,I=1,R=0)

#test the stochastic SIR model
SIR_stoch(theta=init_theta_sir,init_state=init_state_sir,1:59)


###SEIRWL Model (Erlang-distributed waiting times)###
#####################################################
SEIRWL_stoch <- function(theta,init_state,times){
  
  #create transition matrix for Markov process
  trans <- list(
    c(S=-1,E1=1), #susceptible to exposed 1
    c(E1=-1,E2=1), #exposed 1 to exposed 2
    c(E2=-1,I1=1), #exposed 2 to infectious 1
    c(I1=-1,I2=1), #infectious 1 to infectious 2
    c(I2=-1,R1=1), #infectious 2 to recovered 1
    c(R1=-1,R2=1), #recovered 1 to recovered 2
    c(R2=-1,R3=1), #recovered 2 to recovered 3
    c(R3=-1,R4=1), #recovered 3 to recovered 4
    c(R4=-1,W=1), #recovered 4 to window of re-infection
    c(W=-1,L=1), #window of re-infection to immunity
    c(W=-1,E1=1) #window of re-infection to exposed 1
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
    E1 <- current_state[["E1"]]
    E2 <- current_state[["E2"]]
    I1 <- current_state[["I1"]]
    I2 <- current_state[["I2"]]
    R1 <- current_state[["R1"]]
    R2 <- current_state[["R2"]]
    R3 <- current_state[["R3"]]
    R4 <- current_state[["R4"]]
    W <- current_state[["W"]]
    L <- current_state[["L"]]
    I <- I1 + I2
    N <- S + E1 + E2 + I1 + I2 + R1 + R2 + R3 + R4 + W + L
    
    #return rates
    rates <- c(
      beta * (I/N) * S, #susceptible to exposed 1
      2 * f * E1, #exposed 1 to exposed 2
      2 * f * E2, #exposed 2 to infectious 1
      2 * rec * I1, #infectious 1 to infectious 2
      2 * rec * I2, #infectious 2 to recovered 1
      4 * gamma * R1, #recovered 1 to recovered 2
      4 * gamma * R2, #recovered 2 to recovered 3
      4 * gamma * R3, #recovered 3 to recovered 4
      4 * gamma * R4, #recovered 4 to window of re-infection
      tau * W, #window of re-infection to immunity
      beta * (I/N) * W #window of re-infection to exposed 1
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

#initial parameters and state variables
tdc_init_state <- c(S=279,E1=0,E2=0,I1=2,I2=0,R1=3,R2=0,R3=0,R4=0,W=0,L=0)
tdc_theta_seirwl <- c(R0=7,lat_dur=1,inf_dur=4,imm_dur=10,win_dur=0.5)

#test the stochastic SEIRWL model
SEIRWL_stoch(theta=tdc_theta_seirwl,init_state=tdc_init_state,times=1:59)
