Rcpp::sourceCpp('Desktop/git/StochasticInference/ph252b_example/gillespieSSA.cpp')



#create transition matrix for Markov process
transSEIRL <- matrix(
  data = c(
    -1,1,0,0,0,0,0,0,0,0, #susceptible to exposed 1
    0,-1,1,0,0,0,0,0,0,0, #exposed 1 to exposed 2
    0,0,-1,1,0,0,0,0,0,0, #exposed 2 to infectious 1
    0,0,0,-1,1,0,0,0,0,0, #infectious 1 to infectious 2
    0,0,0,0,-1,1,0,0,0,0, #infectious 2 to recovered 1
    0,0,0,0,0,-1,1,0,0,0, #recovered 1 to recovered 2
    0,0,0,0,0,0,-1,1,0,0, #recovered 2 to recovered 3
    0,0,0,0,0,0,0,-1,1,0, #recovered 3 to recovered 4
    0,0,0,0,0,0,0,0,-1,0, #recovered 4 to susceptible and immune
    1,0,0,0,0,0,0,0,0,0, #recovered 4 to susceptible
    0,0,0,0,0,0,0,0,0,1 #recovered 4 to immune
  ),nrow = 11,ncol = 10,byrow = TRUE,dimnames = list(NULL,c("S","E1","E2","I1","I2","R1","R2","R3","R4","L"))
)

tdc_init_state <- c(S=279, E1=0, E2=0, I1=2, I2=0, R1=3, R2=0, R3=0, R4=0, L=0)
theta <- c(R0=9.53,lat_dur=1.96,inf_dur=2.96,imm_dur=10.5,alpha=0.56)

tdcOut = gillespieDirect(theta = theta,initState = tdc_init_state,trans = transSEIRL,tEnd = 100,seed = 42,info = FALSE)

tdcOut = Reduce(f = "rbind",x = tdcOut$trace)
matplot(tdcOut,type="l",lty=1)

# stochasticSEIRL
# 
stochasticSEIRL <- function(theta, initState, transMat, times, seed){
  
  modOut = gillespieDirect(theta = theta,initState = initState,trans = transMat,tEnd = times[length(times)]-1,seed = seed,info = FALSE)
  modTraj = Reduce(f = "rbind",x = modOut$trace)
  
  # interpolate continuous time Markov process output to discrete time points
  modOut$times = modOut$times + min(times)
  intOut = apply(X = modTraj,MARGIN = 2,FUN = function(x){
    approx(x = modOut$times,y = x,xout = times,method = "constant")$y
  })
  
  return(data.frame(time=times,intOut))
}







