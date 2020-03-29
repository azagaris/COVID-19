###################################################################################################
# 'tab.hosp.R' [Antonios Zagaris, 2020]
#
# This routine simulates a simple SIR model for disease transmission, see covid19.Rmd for details.
# Simulation is carried out until the specific time instant where S/N = 1/Ro, where N is the total
# population and Ro the basic reproduction number. The simulation is performed for a grid of Ro
# values spanning the interval [1,maxRo] and the expected number of hospitalizations recorded, as
# are he generating epidemiological parameters and the time (in days) at which the simulation ends.
#
#   INPUTS
# None.
#
#   OUTPUTS
# * out.df = d.f. tabulating expected number of hospitalizations [(Ro, beta,slope,t.star,max.Ihosp) d.f.]
#
#    DEPENDENCIES
# <- settings.R
# <- simulate.epidemic.R
###################################################################################################
# Set working directory and load modules
setwd('/Users/herzog/Documents/Machine/EXTRAMURAL/2016-2020/20-COVID19Modeling-LinkedIn')
source('settings.R') # load parametric settings
source('simulate.epidemic.R') # load epidemic simulator
source('RK4.integrator.R') # load epidemic simulator
source('sir.vf.R') # load RK4 integrator

# Preallocate output
out.df <- data.frame(matrix(NA,nrow=N.Ro,ncol=5))
colnames(out.df) <- c('Ro','beta','slope','t.star','max.Ihosp')# assign colnames

## SIMULATION
# Sanity checks
if(any(params$beta < 0) || any(params$gamma < 0) || p < 1)
{ # abort if no recovery (Ro = oo) or no infection (Ro = 0) or proportion under 1
  stop('ABORTED - no point simulating unless infection and recovery rates are positive [tab.maxI.vs.Ro]')
}
if(sir.o$s < 0 || sir.o$i < 0 || sir.o$r < 0)
{ # abort if initial state has a negative component
  stop('ABORTED - initial state has negative components [tab.maxI.vs.Ro]')
}
# Simulate parameter sets
for(idx in 1:N.Ro)
{
  params.c <- params[idx,] # current parameter set
  tsir.c <- simulate.epidemic(tsir.o , params.c , dt)
  out.df$Ro[idx] <- params.c$Ro # store current Ro value
  out.df$beta[idx] <- params.c$gamma * params.c$Ro # store current beta value
  out.df$slope[idx] <- params.c$beta - params.c$gamma # store current slope value
  out.df$t.star[idx] <- tsir.c$t[nrow(tsir.c)] # store time to max hospitalization
  out.df$max.Ihosp[idx] <- round((0.25/p)*N*tsir.c$i[nrow(tsir.c)]) # store max no. of hospitalizations
}