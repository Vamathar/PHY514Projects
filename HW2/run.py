from Diplotocus import diplotocus as dp
from PoisEqSolver_class import Poisson
import math

# Poisson equation source function
def source(x,y):
    expx=(x-0.5)*(x-0.5)
    expy=(y-0.5)*(y-0.5)
    return -4.*math.pi*math.exp(-16.*(expx+expy))

# y=0 and y=1 boundary values are the same in the given problem
def ybound(x):
    return 1.-x

def main():
    # Declare inputs for Poisson equation solver class  
    xminmax = (0.,1.)
    yminmax = (0.,1.)
    boundx = (1.,0.)
    boundy = (ybound,ybound)
    steps = 10

    # Instantiate Poisson equation solver
    SolveMe = Poisson(source,xminmax,yminmax,steps,boundx,boundy)

    # Relaxation method
    title = 'Poisson Equation Solution using Relaxation Method'
    RelaxIters = SolveMe.Relaxation(0.01)
    x,y,z = SolveMe.Grids()
    Plotter = dp(title)
    Plotter.Contour(x,y,z,30,'plasma')



main()
