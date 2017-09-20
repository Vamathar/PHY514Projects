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
        
        # Initialize zgrid with proper boundary values
        self.zgrid = False
        self.__zgridInit()
        self.__solved = False # Indicates if a solver has been used
        

    def __verify(self):
        """
        Checks that class variables are proper types and formats.
        """
        errors = list()
        if type(self.__src)!=type(self.__verify):
            errors.append('source')

    def __zgridInit(self):
        """
        Initializes zgrid with proper boundary values
        """
        self.zgrid = list()
        i = 0
        j = 0
        for x in self.__xc:
            column = list()
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
            self.zgrid.append(column)
            i+=1
        self.zgrid = np.array(self.zgrid)
    
    def Grids(self):
        """
        Returns grids of x,y,z coordinates

        Input: none
        Returns: list(np.array(steps,steps),np.array(steps,steps),np.array(steps,steps)), i.e. [xgrid,ygrid,zgrid]
        """
        xgrid,ygrid = np.meshgrid(self.__xc,self.__yc)
        return [xgrid,ygrid,self.zgrid]

    @jit
    def Relaxation(self,delta,maxiter=1000):
        """
        Solves 2D Poisson equation via relaxation mehtod.

        Input:
            delta: float, max difference after one relaxation step
            maxiter: maximum iterations allowed
        Returns: integer, number of iterations to converge or maxiter

        Note: Access solution through the Poisson.zgrid variable
        """
        # Reinitialize zgrid if another solver has already been used
        if self.__solved:
            self.__zgridinit()
        self.__solved = True

        diff = delta+1.
        iterations = 0
        newzgrid = self.zgrid
        print('length = ',len(self.__xc),len(self.__yc))
        while diff > delta and iterations < maxiter:
            i = 0
            j = 0
            # Apply relaxation step
            for x in self.__xc:
                for y in self.__yc:
                    boundary = i==0\
                            or i == (self.__steps - 1)\
                            or j == 0\
                            or j == (self.__steps - 1)

                    if not boundary:
                        newzgrid[i][j] = 0.25*(\
                                self.zgrid[i+1][j]+self.zgrid[i-1][j]\
                                +self.zgrid[i][j+1]+self.zgrid[i][j-1])\
                                +math.pi*self.__src(self.__xc[i],self.__yc[j])\
                                *self.__dx2
                    print('i,j = ',i,j)
                    j+=1
                i+=1
            diffzgrid = abs(newzgrid - self.zgrid)
            diff = diffzgrid.max()
            self.zgrid = newzgrid
            iterations+=1
            print('diff,delta = ',diff,delta)

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
