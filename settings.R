###################################################################################################
# 'settings.R' [Antonios Zagaris, 2020]
#
# This routine loads the settings for the simulator tab.maxI.vs.Ro feeding into covid19.Rmd.
#
#   INPUTS
# None.
#
#   OUTPUTS
# Various parameter setting, see inline.
#
#    DEPENDENCIES
# -> tab.hosp.R
###################################################################################################
## SETTINGS
# Basic parameters
eps <- 1E-12 # tolerance in checking (in)equality
N <- 17.18E6 # NL population
p <- 10 # ratio of true infecteds to officially recorded infecteds
q <- 0.25 # ratio of officially ecorded infecteds that are hospitalized
beds <- 37753 # numbers of beds in NL
# Epidemiological parameters
params <- data.frame('beta'=0.27 , 'gamma'=0.07) # infection and recovery rates [1/day]
params$Ro <- params$beta / params$gamma # basic reproduction number
N.Ro <- 100 # no. of Ro-values in Ro-grid
aux <- data.frame(matrix(NA , nrow=1+N.Ro , ncol=3)) # expand matrix
colnames(aux) <- c('beta','gamma','Ro') # fix colnames
aux$Ro <- seq(from = 1 , to = params$Ro , length.out = 1+N.Ro) # Ro-grid
aux$beta <- params$gamma * aux$Ro # beta-grid
aux$gamma <- params$gamma # gamma remains constant
aux <- aux[-1,] # remove parameter setting with Ro=1
row.names(aux) <- 1:nrow(aux) # renumber rows after removing 1st one
params <- aux # parameter-grid
dt <- 0.1 # timestep
# Initial condition
SIR.o <- data.frame('S'=NA , 'I'=4204*p , 'R'=0) # initial condition (time zero for simulation)
SIR.o$S <- N - SIR.o$I - SIR.o$R # no. of susceptibles
sir.o <- SIR.o / N # compartmental proportions E [0,1]
colnames(sir.o) <- c('s','i','r') # fix colnames
tsir.o <- data.frame('t'=0 , sir.o) # append time
