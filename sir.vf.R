###################################################################################################
# 'sir.vf.R' [Antonios Zagaris, 2020]
#
# This routine outputs the vector field for the simple SIR model investigated in covid19.Rmd. The
# model there is
#   d(S,I,R)/dt = (-beta*S*I/N , beta*S*I/N - gamma*I , gamma*I) , with N the total population. [1]
# Here, I have written (s,i,r) = (S,I,R)/N so that
#   d(s,i,r)/dt = (-beta*s*i , beta*s*i - gamma*i , gamma*i). [2]
# The system state (s,i,r) acts as input for the routine and the vector in the RHS is the output.
#
#   INPUTS
# * sir = d.f. holding the values of the scaled system variables [(s,i,r) x 1 d.f.]
# * params = d.f. holding the values of the parameters entering the model [(beta,gamma,Ro) x 1 d.f.]
#
#   OUTPUTS
# * sir.vf = d.f. holding the vector field in the RHS of [2] above [(s,i,r) x 1 d.f.]
#
#    DEPENDENCIES
# -> RK4.ntegrator.R
###################################################################################################
sir.vf <- function(sir,params)
{
  vf <- data.frame('s'=NA , 'i'=NA , 'r'=NA) # preallocate vf
  vf$s <- - params$beta * sir$i * sir$s  # s-component of vf
  vf$i <- + sir$i * (params$beta * sir$s - params$gamma)  # i-component of vf
  vf$r <- + params$gamma * sir$i  # r-component of vf
  return(vf)
}