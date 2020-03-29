###################################################################################################
# 'simulate.epidemic.R' [Antonios Zagaris, 2020]
#
# This routine simulates the simple SIR model investigated in covid19.Rmd over a time period,
# specifically until I reaches its peak. The model is
#   d(S,I,R)/dt = (-beta*S*I/N , beta*S*I/N - gamma*I , gamma*I) , with N the total population. [1]
# However, I work throughout with the scaled version in which (s,i,r) = (S,I,R)/N and thus
#   d(s,i,r)/dt = (-beta*s*i , beta*s*i - gamma*i , gamma*i). [2]
# The initial system state (s.o,i.o,r.o) acts as input for the routine, as do the parameters beta
# and gamma and the timestep dt. The simulation ends when s <= gamma/beta = 1/Ro. Provided that
# dt << slowest timescale, this is close to the precise time t* where I is maximized.
#
# Note that the while loop in the code comes without the usual bells and whistles surrounding a
# potentially infinite loop, so use wisely. Also not that outputting the final sate would be
# sufficient to track duration and max no. of infecteds, so you may modify accordingly.
#
#   INPUTS
# * tsir.o = d.f. holding the values of the initial time (t) and the initial system state (s,i,r)
#            [(t,s,i,r) x 1 d.f.]
# * params = d.f. holding the values of the parameters entering the model [(beta,gamma,Ro) x 1 d.f.]
# * dt = timestep over which the scaled system is integrated [double]
#
#   OUTPUTS
# * simulate.epidemic = d.f. holding the times and system states for the entire duration of the
#                       epidemic [(t,s,i,r) x N.t d.f.]
#
#    DEPENDENCIES
# -> tab.hosp.R
# <- RK4.integrator.R
###################################################################################################
simulate.epidemic <- function(tsir.o , params , dt)
{
  s.trgt <- 1/params$Ro # target s-value (triggerring simulation end)

  tsir.f <- tsir.o # initialize trajectory

  while(tsir.f$s[nrow(tsir.f)] > s.trgt)
  { # iterate until susceptibles equal/drop below target value
    tsir.new <- RK4.integrator(tsir.f[nrow(tsir.f),] , params , dt) # integrate ODEs over a timestep
    tsir.f <- rbind(tsir.f , tsir.new , make.row.names=FALSE) # append simulated point to trajectory
  }

  return(tsir.f)
}