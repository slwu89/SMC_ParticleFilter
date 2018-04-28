################################################################################
# 
#   PH252D (Spring 2018): Modelling the Dynamics of Infectious Disease Processes
#   Inference for Stochastic Models with pomp
#   Sean Wu (slwu89@berkeley.edu)
# 
################################################################################

library(pomp)

################################################################################
#   SEIRWL (window of re-infection)
################################################################################

# ODE version of "window of re-infection" model
SEIRWL_ode_c <- Csnippet("
  
  double rates[6];                  /* rates of movement between compartments */
  
  /* derived parameters */
  double beta = R0 / inf_d;         /* effective contact rate */
  double f = 1 / lat_d;
  double r = 1 / inf_dur;
  double gamma = 1 /imm_d;
  double tau = 1 / win_d;

  /* total population size */  
  double N = S + E + I + R + W + L;

  /* rates of movement between compartments */
  rates[0] = (beta * I / N) * S;    /* force of infection on susceptibles */
  rates[1] = f * E;                 /* progression to infectiousness */
  rates[2] = r * I;                 /* recovery */
  rates[3] = gamma * R;             /* progression to re-infectiousness */
  rates[4] = tau * W;               /* permanant recovery */
  rates[5] = (beta * I / N) * W;    /* force of infection on window-of-reinfection */
  
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

  /* derived parameters */
  double beta = R0 / inf_d;         /* effective contact rate */
  double f = 1 / lat_d;
  double r = 1 / inf_dur;
  double gamma = 1 /imm_d;
  double tau = 1 / win_d;
  
  /* total population size */  
  double N = S + E + I + R + W + L;

  /* rates of each event occuring */
  rates[0] = (beta * I / N) * S;    /* S -> E: force of infection on susceptibles */
  rates[1] = f * E;                 /* E -> I: progression to infectiousness */
  rates[2] = r * I;                 /* I -> R: recovery */
  rates[3] = gamma * R;             /* R -> W: progression to re-infectiousness */
  rates[4] = tau * W;               /* W -> L: permanant recovery */
  rates[5] = (beta * I / N) * W;    /* W -> E: force of infection on window-of-reinfection */

  /* sample events */
  reulermultinom(1, S, &rate[0], dt, &dX[0]); /* S -> E */
  reulermultinom(1, E, &rate[1], dt, &dX[1]); /* E -> I */
  reulermultinom(1, I, &rate[2], dt, &dX[2]); /* I -> R */
  reulermultinom(1, R, &rate[3], dt, &dX[3]); /* R -> W */
  reulermultinom(1, W, &rate[4], dt, &dX[4]); /* W -> L */
  reulermultinom(1, W, &rate[5], dt, &dX[5]); /* W -> E */

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
  obs = rpois(rho * Inc > 0 ? rho * Inc : 0);
")

# point likelihood 
SEIRWL_dmeasurement_c <- Csnippet("
  lik = dpois(obs, rho * Inc > 0 ? rho * Inc : 0, give_log);
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

# pomp object
SEIRWL_pomp <- pomp(data = FluTdC1971[, c("time", "obs")],
                   rprocess = euler.sim(step.fun = SEIRWL_stochastic_c,
                                        delta.t = 0.1),
                   rmeasure = SEIRWL_rmeasurement_c,
                   dmeasure = SEIRWL_dmeasurement_c,
                   skeleton = vectorfield(SEIRWL_ode_c),
                   times = "time",
                   t0 = 1,
                   zeronames = "Inc", 
                   paramnames = c("R0", "D_inf", "D_lat", "D_imm", "alpha", "rho"),
                   statenames = c("S", "E", "I", "T", "L", "Inc"),
                   obsnames = c("obs"))


pmcmc(
  pomp(ou2,dprior=Csnippet("
                           lik = dnorm(alpha_2,-0.5,1,1) + dnorm(alpha_3,0.3,1,1);
                           lik = (give_log) ? lik : exp(lik);"),
       paramnames=c("alpha.2","alpha.3")),
  Nmcmc=2000,Np=500,verbose=TRUE,
  proposal=mvn.rw.adaptive(rw.sd=c(alpha.2=0.01,alpha.3=0.01),
                           scale.start=200,shape.start=100)) -> chain
continue(chain,Nmcmc=2000,proposal=mvn.rw(covmat(chain))) -> chain
plot(chain)
chain <- pmcmc(chain)
plot(chain)
