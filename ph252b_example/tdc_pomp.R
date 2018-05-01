################################################################################
# 
#   PH252D (Spring 2018): Modelling the Dynamics of Infectious Disease Processes
#   Inference for Stochastic Models with pomp
#   Sean Wu (slwu89@berkeley.edu)
# 
################################################################################

library(tidyverse)
library(pomp)
library(fitR)
# data("FluTdC1971")

dir <- "/Users/slwu89/Desktop/git/StochasticInference/"
FluTdC1971 <- read.csv2(paste0(dir,"FluTdCIncData.csv"),header = T,sep = ",")

################################################################################
#   SEIRWL (window of re-infection)
################################################################################

# ODE version of "window of re-infection" model
SEIRWL_ode_c <- Csnippet("
                         
                         double rates[6];                  /* rates of movement between compartments */
                        
                         /* total population size */  
                         double N = S + E + I + R + W + L;

                         /* derived parameters */
                         double beta = R0 / inf_d;         /* effective contact rate */
                         double f = 1 / lat_d;
                         double r = 1 / inf_d;
                         double gamma = 1 /imm_d;
                         double tau = 1 / win_d;
                         double lambda = beta * (I / N);
                         
                         /* rates of movement between compartments */
                         rates[0] = lambda * S;    /* force of infection on susceptibles */
                         rates[1] = f * E;                 /* progression to infectiousness */
                         rates[2] = r * I;                 /* recovery */
                         rates[3] = gamma * R;             /* progression to re-infectiousness */
                         rates[4] = tau * W;               /* permanant recovery */
                         rates[5] = lambda * W;    /* force of infection on window-of-reinfection */
                         
                         /* dx/dt */
                         DS = -rates[0];
                         DE = rates[0] + rates[5] - rates[1];
                         DI = rates[1] - rates[2];
                         DR = rates[2] - rates[3];
                         DW = rates[3] - rates[4] - rates[5];
                         DL = rates[4];
                         DInc = rates[0] + rates[5];       /* process model case incidence */
                         ")

# stochastic Markovian process model of "window of re-infection" model
SEIRWL_stochastic_c <- Csnippet("
                                
                                double rates[6];                  /* rates that events occur at */
                                double dX[6];                     /* vector to hold sampled event sizes */
                                
                                /* total population size */  
                                double N = S + E + I + R + W + L;
                                
                                /* derived parameters */
                                double beta = R0 / inf_d;         /* effective contact rate */
                                double f = 1 / lat_d;
                                double r = 1 / inf_d;
                                double gamma = 1 /imm_d;
                                double tau = 1 / win_d;
                                double lambda = beta * (I / N);
                                
                                /* rates of each event occuring */
                                rates[0] = lambda * S;    /* force of infection on susceptibles */
                                rates[1] = f * E;                 /* progression to infectiousness */
                                rates[2] = r * I;                 /* recovery */
                                rates[3] = gamma * R;             /* progression to re-infectiousness */
                                rates[4] = tau * W;               /* permanant recovery */
                                rates[5] = lambda * W;    /* force of infection on window-of-reinfection */
                                
                                /* sample events */
                                reulermultinom(1, S, &rates[0], dt, &dX[0]); /* S -> E */
                                reulermultinom(1, E, &rates[1], dt, &dX[1]); /* E -> I */
                                reulermultinom(1, I, &rates[2], dt, &dX[2]); /* I -> R */
                                reulermultinom(1, R, &rates[3], dt, &dX[3]); /* R -> W */
                                reulermultinom(2, W, &rates[4], dt, &dX[4]); /* W -> L */
                                // reulermultinom(1, W, &rates[5], dt, &dX[5]); /* W -> E */
                                
                                /* update state variables */
                                S += -dX[0];
                                E += dX[0] + dX[5] - dX[1];
                                I += dX[1] - dX[2];
                                R += dX[2] - dX[3];
                                W += dX[3] - dX[4] - dX[5];
                                L += dX[4];
                                Inc += dX[0] + dX[5];             /* process model case incidence */
                                ")

# stochastic measurement model
SEIRWL_rmeasurement_c <- Csnippet("
                                  inc = rpois(rho * I > 0 ? rho * I : 0);
                                  ")

# point likelihood 
SEIRWL_dmeasurement_c <- Csnippet("
                                  lik = dpois(inc, rho * I > 0 ? rho * I : 0, give_log);
                                  ")

# transform parameters to unconstrained space
SEIRWL_toEst_c <- Csnippet("
                           TR0     = log(R0);
                           Tinf_d  = log(inf_d);
                           Tlat_d  = log(lat_d);
                           Timm_d  = log(imm_d);
                           Twin_d  = log(win_d);
                           Trho    = logit(rho);
                           ")

# transform parameters to natural space
SEIRWL_fromEst_c <- Csnippet("
                             TR0     = exp(R0);
                             Tinf_d  = exp(inf_d);
                             Tlat_d  = exp(lat_d);
                             Timm_d  = exp(imm_d);
                             Twin_d  = exp(win_d);
                             Trho    = expit(rho);
                             ")

# prior distribution for MCMC
SEIRWL_prior_c <- Csnippet("
                           lik =  dunif(R0, 1, 50, 1);
                           lik += dunif(inf_d, 0, 20, 1);
                           lik += dunif(lat_d, 0, 20, 1);
                           lik += dunif(imm_d, 0, 20, 1);
                           lik += dunif(win_d, 0, 20, 1);
                           lik += dunif(rho, 0, 1, 1);
                           lik = give_log ? lik: exp(lik);
                           ")

# initial state of hidden process
SEIRWL_initializer_c <- Csnippet("
                                 S = 279;
                                 E = 0;
                                 I = 2;
                                 R = 3;
                                 W = 0;
                                 L = 0;
                                 ")

# pomp object
SEIRWL_pomp <- pomp(data = FluTdC1971[, c("time", "inc")],
                   rprocess = euler.sim(step.fun = SEIRWL_stochastic_c,
                                        delta.t = 0.1),
                   rmeasure = SEIRWL_rmeasurement_c,
                   dmeasure = SEIRWL_dmeasurement_c,
                   skeleton = vectorfield(SEIRWL_ode_c),
                   initializer = SEIRWL_initializer_c,
                   times = "time",
                   t0 = 1,
                   zeronames = "Inc", 
                   paramnames = c("R0","inf_d","lat_d","imm_d","win_d","rho"),
                   statenames = c("S","E","I","R","W","L","Inc"),
                   obsnames = c("inc"))

theta <- c(R0=5,inf_d=2.5,lat_d=2,imm_d=4,win_d=2,rho=0.8)
sim <- simulate(SEIRWL_pomp,params=theta,nsim=50,as.data.frame=TRUE,states=TRUE)
traj <- trajectory(SEIRWL_pomp,params=theta)
pf <-  pfilter(SEIRWL_pomp,Np=1000,params=c(R0=5,inf_d=2,lat_d=1,imm_d=5,win_d=2,rho=0.8))
 
SEIRWL_pomp %>% 
  simulate(params=theta,nsim=9,as.data.frame=TRUE,include.data=TRUE) %>%
  ggplot(aes(x=time,y=inc,group=sim,color=(sim=="data")))+
  guides(color=FALSE)+
  geom_line()+facet_wrap(~sim,ncol=2)


prop.sd <- c(R0=0.5,inf_d=1,lat_d=1,imm_d=1,win_d=1,rho=0.1)
pmcmc_out <- pmcmc(SEIRWL_pomp, Nmcmc = 1e3, Np = 5e2,start = theta,verbose=TRUE, proposal = mvn.rw.adaptive(rw.sd = prop.sd, 
                                                                               scale.start = 100, shape.start = 200), max.fail = Inf)

lhs <- sobolDesign(lower=c(R0=1,inf_d=1,lat_d=1,imm_d=1,win_d=1,rho=0.1),
                   upper=c(R0=10,inf_d=20,lat_d=20,imm_d=20,win_d=20,rho=0.99),nseq=100)

mif_out <- mif2(SEIRWL_pomp, Nmif = 100, Np = 500, cooling.fraction.50 = 0.01, cooling.type = "hyperbolic",
                rw.sd = rw.sd(R0=0.5,inf_d=0.5,lat_d=0.5,imm_d=0.5,win_d=0.5,rho=0.1), start = theta, verbose = TRUE)

ggplot(data=melt(conv.rec(mif_out)),
       aes(x=iteration,y=value))+
  geom_line()+
  guides(color=FALSE)+
  facet_wrap(~variable,scales="free_y")+
  theme_bw()

SEITL.mf.sim <- simulate(mif_out, nsim = 10, as.data.frame = TRUE, include.data = TRUE)
plotTraj(SEITL.mf.sim, data = FluTdC1971, state.names = "inc")
