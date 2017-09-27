################################################################################
# Description: 1D heat equation solver 
# Purpose: For problem 2 on PHY514 HW2
# Author: Noah Green
################################################################################
import numpy as np
import math

class Heat:
    def __init__(self,Tinit,const,xmm,xsteps,tmax,tsteps):
        self.__const = const
        self.__xmm = xmm
        self.__xsteps = xsteps
        self.__tmax = tmax
        self.__tsteps = tsteps
        x = np.linspace(self.__xmm[0],self.__xmm[1],self.__xsteps)
        t = np.linspace(0.,self.__tmax,self.__tsteps)
        self.__dx = abs(x[1]-x[0])
        self.__dx2 = self.__dx*self.__dx
        self.__dt = abs(t[1]-t[0])
        self.x,self.t = np.meshgrid(x,t)
        self.Temp = self.x.copy()
        for i in range(self.Temp.shape[0]):
            for j in range(self.Temp.shape[1]):
                self.Temp[i][j] = Tinit(self.Temp[i][j])
    def LineSolve(self):
        """
        Solves the 1D heat equation using method of lines

        Input: none
        Returns: nothing
        
        Note: access solution through public class variables
        """
        for i in range(self.Temp.shape[0]-1):
            for j in range(self.Temp.shape[1]-1):
                if j != 0:
                    self.Temp[i+1][j] = self.Temp[i][j]\
                            +(self.__const*self.__dt/self.__dx2)\
                            *(self.Temp[i][j+1]+\
                            self.Temp[i][j-1]-2.*self.Temp[i][j])
                
