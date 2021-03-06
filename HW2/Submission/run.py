from Diplotocus import diplotocus as dp
from PoisEqSolver_class import Poisson
from Heat_class import Heat
import math

# Poisson equation source function
def source(x,y):
    expx=(x-0.5)*(x-0.5)
    expy=(y-0.5)*(y-0.5)
    return -4.*math.pi*math.exp(-16.*(expx+expy))

# y=0 and y=1 boundary values are the same problems 1 & 3
def ybound(x):
    return 1.-x

# Initial temperature in problem 2
def Tinit(x):
    if x > -.1 and x < .1:
        return 1.
    else:
        return 0.
def main():
    # Declare inputs for Poisson equation solver class  
    xminmax = (0.,1.)
    yminmax = (0.,1.)
    boundx = (1.,0.)
    boundy = (ybound,ybound)
    steps = 10

    # Instantiate Poisson equation solver
    SolveMe = Poisson(source,xminmax,yminmax,steps,boundx,boundy)

################################################################################
#   Problem 1.1
################################################################################
    # Relaxation Method
    convergence = 0.01
    title = 'Poisson Equation Solution using Relaxation Method'
    RelaxIters = SolveMe.Relaxation(convergence)
    print('Iterations to solution: ',RelaxIters)
    x,y,z = SolveMe.Grids()
    Plotter = dp(title)
    Plotter.Contour(x,y,z,30,'plasma','Problem1Part1Plot.png')


################################################################################
#   Problem 1.2
################################################################################
    # Gauss-Seidel Method
    convergence = 0.01
    w = 1.7
    title = 'Poisson Equation Solution using G-S Method'
    RelaxIters = SolveMe.GaussSeidel(w,convergence)
    print('Iterations to solution: ',RelaxIters)
    x,y,z = SolveMe.Grids()
    Plotter = dp(title)
    Plotter.Contour(x,y,z,30,'plasma','Problem1Part2Plot.png')


################################################################################
#   Problem 3
################################################################################
    # Multi-Grid Method
    convergence = 0.01
    w = 1.
    FinalSteps = steps*10
    title = 'Poisson Equation Solution using G-S Multi-Grid Method'
    RelaxIters = SolveMe.MultiGrid(w,convergence,FinalSteps,maxiter=10000)
    x,y,z = SolveMe.Grids()
    Plotter = dp(title)
    Plotter.Contour(x,y,z,30,'plasma','Problem3Plot.png')



################################################################################
#   Problem 2
################################################################################
    xmm = (-1.,1.)
    xsteps = 20
    tmax = 10.
    tsteps = 10
    dx = abs(xmm[1]-xmm[0])/(xsteps-1)
    dt = tmax/(tsteps-1)
    # Ensure stability of solution with following diffusion constant value
    const = 0.25*dx*dx/dt
    
    # Solve 1 dimensional heat equation
    toasty = Heat(Tinit,const,xmm,xsteps,tmax,tsteps)
    toasty.LineSolve()
    title = 'Heat Equation Solution Using Method of Lines'
    Plotter = dp(title,'t','x')
    Plotter.Contour(toasty.t,toasty.x,toasty.Temp,30,'plasma','Problem2Plot.png')
    
main()
