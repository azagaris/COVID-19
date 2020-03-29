###################################################################################################
# 'RK4.integrator.R' [Antonios Zagaris, 2020]
#
# This routine integrates the simple SIR model investigated in covid19.Rmd over a single timestep.
# The model is
#   d(S,I,R)/dt = (-beta*S*I/N , beta*S*I/N - gamma*I , gamma*I) , with N the total population. [1]
# Here, however, I work with the scaled version in which (s,i,r) = (S,I,R)/N and thus
#   d(s,i,r)/dt = (-beta*s*i , beta*s*i - gamma*i , gamma*i). [2]
# The initial system state (s.o,i.o,r.o) acts as input for the routine, as do the parameters beta
# and gamma and the timestep dt.
#
#   INPUTS
# * tsir.o = d.f. holding the values of the initial time (t) and the initial system state (s,i,r)
#            [(t,s,i,r) x 1 d.f.]
# * params = d.f. holding the values of the parameters entering the model [(beta,gamma,Ro) x 1 d.f.]
# * dt = timestep over which the scaled system is integrated [double]
#
#   OUTPUTS
# * RK4.integrator = d.f. holding the final time and system state [(t,s,i,r) x 1 d.f.]
#
#    DEPENDENCIES
# -> simulate.epidemic.R
# <- sir.vf.R
###################################################################################################
RK4.integrator <- function(tsir.o , params , dt)
{
  tsir.f <- tsir.o # preallocate output as input (to inherit colnames)
  
  sir.o <- tsir.o[,c('s','i','r')] # isolate sir variables (i.e. exclude time)
  k.1 <- sir.vf(sir.o , params) # vf at  original system state
  sir.1 <- sir.o +  0.5 * dt * k.1 # 1st intermediate system state

  k.2 <- sir.vf(sir.1 , params) # vf at 1st intermediate system state
  sir.2 <- sir.o +  0.5 * dt * k.2 # 2nd intermediate system state

  k.3 <- sir.vf(sir.2 , params) # vf at 2nd intermediate system state
  sir.3 <- sir.o +  dt * k.3 # 3rd intermediate system state

  k.4 <- sir.vf(sir.3 , params) # vf at 3rd intermediate system state
  sir.f <- sir.o +  1/6 * dt * (k.1 + 2*k.2 + 2*k.3 + k.4) # final system state
  
  tsir.f$t <- dt + tsir.o$t # update time
  tsir.f[,c('s','i','r')] <- sir.f # update system state
  
  return(tsir.f)
}