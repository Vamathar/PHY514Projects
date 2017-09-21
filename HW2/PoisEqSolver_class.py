################################################################################
# Description: 2D Poisson equation solver 
# Purpose: For problems 1 and 3 on PHY514 HW2
# Author: Noah Green
################################################################################
from numba import jit
import numpy as np
import math

class Poisson:
    def __init__(self,source,xminmax,yminmax,steps,boundsx,boundsy):
        self.__src = source
        self.__xmm = xminmax
        self.__ymm = yminmax
        self.__steps = steps
        self.__bx = boundsx
        self.__by = boundsy

        # Dictionary of possible errors
        self.__errors = {
                'source':'Error: source must be a function.',
                'xminmax':'Error: invalid minimum and maximum x values.',
                'yminmax':'Error: invalid minimum and maximum y values.',
                'boundsx':'Error: invalid boundary values in x.',
                'boundsy':'Error: invalid boundary values in y.'
                }
        # Verify inputs are valid
        self.__verify()

        # Make x and y coordinates
        self.__xc = np.linspace(self.__xmm[0],self.__xmm[1],self.__steps)
        self.__yc = np.linspace(self.__ymm[0],self.__ymm[1],self.__steps)
        self.__dx = self.__xc[1]-self.__xc[0]
        self.__dx2 = self.__dx*self.__dx
        
        # Initialize phigrid with proper boundary values
        self.phigrid = False
        self.__phigridInit()
        self.__solved = False # Indicates if a solver has been used
        

    def __verify(self):
        """
        Checks that class variables are proper types and formats.
        """
        errors = list()
        if type(self.__src)!=type(self.__verify):
            errors.append('source')

    def __phigridInit(self):
        """
        Initializes phigrid with proper boundary values
        """
        self.phigrid = list()
        i = 0
        for x in self.__xc:
            column = list()
            j=0
            for y in self.__yc:
                if i == 0:
                    column.append(self.__bx[0])
                elif i == (self.__steps-1):
                    column.append(self.__bx[1])
                elif j == 0:
                    column.append(self.__by[0](self.__xc[i]))
                else:
                    column.append(self.__by[1](self.__xc[i]))
                j+=1
            self.phigrid.append(column)
            i+=1
        self.phigrid = np.array(self.phigrid)
    
    def Grids(self):
        """
        Returns grids of x,y,z coordinates

        Input: none
        Returns: list(np.array(steps,steps),np.array(steps,steps),np.array(steps,steps)), i.e. [xgrid,ygrid,phigrid]
        """
        xgrid,ygrid = np.meshgrid(self.__xc,self.__yc)
        return [xgrid,ygrid,self.phigrid]

    @jit
    def Relaxation(self,delta,maxiter=1000):
        """
        Solves 2D Poisson equation via relaxation mehtod.

        Input:
            delta: float, max difference after one relaxation step
            maxiter: maximum iterations allowed
        Returns: integer, number of iterations to converge or maxiter

        Note: Access solution through the Poisson.phigrid variable
        """
        # Reinitialize phigrid if another solver has already been used
        if self.__solved:
            self.__phigridinit()
        self.__solved = True

        diff = delta+1.
        iterations = 0
        newphigrid = self.phigrid.copy()
        while (diff > delta and iterations < maxiter):
            i = -1
            # Apply relaxation step
            for x in self.__xc:
                i+=1
                j=-1
                for y in self.__yc:
                    j+=1
                    boundary = i == 0\
                            or i == (self.__steps - 1)\
                            or j == 0\
                            or j == (self.__steps - 1)
                    if not boundary:
                        phixp1 = self.phigrid[i+1][j]
                        phixm1 = self.phigrid[i-1][j]
                        phiyp1 = self.phigrid[i][j+1]
                        phiym1 = self.phigrid[i][j-1]
                        avg = (phixp1+phixm1+phiyp1+phiym1)*0.25
                        val = avg\
                                +math.pi*self.__src(x,y)
                        newphigrid[i][j] = val 

            diffphigrid = abs(newphigrid - self.phigrid)
            #print(diffphigrid)
            diff = diffphigrid.max()
            self.phigrid = newphigrid.copy()
            iterations+=1
        if iterations >= maxiter:
            print("Warning: Relaxation solver did not converge.")
        return iterations

    def GaussSeidel(self):
        """
        Solves 2D Poisson equation via Gauss-Seidel method.
        """
        return 1
    def MultiGrid(self):
        """
        Solves 2D Poisson equation via Multi-Grid method.
        """
        return 1
