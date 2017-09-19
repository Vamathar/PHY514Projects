################################################################################
# Author: Noah Green
# Class: PHY 514
# Assignment: Homework 1
# Description: Harmonic oscillator numerical integrator class for homework 1.
###############################################################################

import math
import numpy as np

class HO_int:
    def __init__(self,x0,v0,omega,delta_t,tfinal):
        """
        Initialize harmonic oscillator class

        Input:
            x0: float, initial position
            v0: float, initial velocity
            omega: float, angular frequency
            delta_t: float, size of time step
            tfinal: float, final time
        """
        self.__dt = delta_t
        self.__x0 = x0
        self.__v0 = v0
        self.__w = omega
        self.__w2 = omega*omega
        self.__tf = tfinal
        self.time = self.__time__()

    def __time__(self):
        """
        Returns list of equally spaced times up to final time
        
        Input: none

        Returns: list(float), i.e. [times]
        """
        N = int(self.__tf/self.__dt)+1
        time = list()
        tval = 0.
        for i in range(N):
            time.append(tval)
            tval += self.__dt
        return time

    def Analytical(self):
        """
        Returns list of positions and velocities using analytical solution.

        Input: none

        Returns: list(list(float),list(float)), \
                i.e. [[position list],[velocity list]]
        """
        xlist = list()
        vlist = list()
        for i in range(len(self.time)):
            xlist.append((self.__v0/self.__w)*math.sin(self.__w*self.time[i])\
                    +self.__x0*math.cos(self.__w*self.time[i]))
            vlist.append(self.__v0*math.cos(self.__w*self.time[i])-self.__x0\
                    *self.__w*math.sin(self.__w*self.time[i]))
        return (xlist,vlist)

    def ForwardEuler(self):
        """
        Returns list of positions and velocities using forward Euler algorithm.

        Input: none

        Returns: list(list(float),list(float)), \
                i.e. [[position list],[velocity list]]
        """
        xlist = list()
        vlist = list()
        xlist.append(self.__x0)
        vlist.append(self.__v0)

        for i in range(len(self.time)-1):
            xlist.append(xlist[i]+vlist[i]*self.__dt)
            vlist.append(vlist[i]-self.__w2*xlist[i]*self.__dt)
        return (xlist,vlist)

    def BackwardEuler(self):
        """
        Returns list of positions and velocities using backward Euler algorithm.

        Input: none

        Returns: list(list(float),list(float)), \
                i.e. [[position list],[velocity list]]
        """
        xlist = list()
        vlist = list()
        xlist.append(self.__x0)
        vlist.append(self.__v0)
        for i in range(len(self.time)-1):
            xplusone = xlist[i]+vlist[i]*self.__dt-self.__w2*self.__dt**2*xlist[i]
            xlist.append(xplusone)
            vlist.append(vlist[i]-self.__w2*xplusone*self.__dt)
        return (xlist,vlist)

    def RungeKutta(self):
        """
        Returns list of positions and velocities using RK4 algorithm.

        Input: none

        Returns: list(list(float),list(float)), \
                i.e. [[position list],[velocity list]]
        """
        xlist = list()
        vlist = list()
        xlist.append(self.__x0)
        vlist.append(self.__v0)
        for i in range(len(self.time)-1):
            k1v = -self.__w2*xlist[i]*self.__dt
            k1x = vlist[i]*self.__dt
           
            k2v = -self.__w2*(xlist[i]+k1x/2.)*self.__dt
            k2x = (vlist[i]+k1v/2.)*self.__dt

            k3v = -self.__w2*(xlist[i]+k2x/2.)*self.__dt
            k3x = (vlist[i]+k2v/2.)*self.__dt

            k4v = -self.__w2*(xlist[i]+k3x)*self.__dt
            k4x = (vlist[i]+k3v)*self.__dt

            vlist.append(vlist[i]+(k1v+2.*k2v+2.*k3v+k4v)/6.)
            xlist.append(xlist[i]+(k1x+2.*k2x+2.*k3x+k4x)/6.)
        return (xlist,vlist)

    def LeapFrog(self):
        """
        Returns list of positions and velocities using leapfrog algorithm.

        Input: none

        Returns: list(list(float),list(float)), \
                i.e. [[position list],[velocity list]]
        """
        xlist = list()
        vlist = list()
        xlist.append(self.__x0)
        vlist.append(self.__v0)
        for i in range(len(self.time)-1):
            if i == 0: # Jumpstart with a forward Euler step
                vlist.append(vlist[i]-0.5*self.__w2*xlist[i]*self.__dt)
            else:
                vlist.append(vlist[i]-self.__w2*xlist[i]*self.__dt)
            xlist.append(xlist[i]+vlist[i+1]*self.__dt)
        return (xlist,vlist)

