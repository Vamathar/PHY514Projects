################################################################################
# Name: Noah Green
# Class: PHY 514
# Assignment: HW1 problem 2
# Description: Class for shooting problme
################################################################################
import math

class shoot:
    def __init__(self,distance,v0):
        self.d = distance
        self.v0 = v0
        self.g = 9.81

    def analytical_height(self,theta,dt):
        """
        Analytically calculates height of cannon ball until it hits the gound.

        Input:
            dt: time step
            
        Returns:list(list(float),list(float)), [[time],[height]]
        """
        tf = 2.*v0*math.sin(theta)/g
        ylist = list()
        tlist = list()
        y = 0.
        t = 0.
        ylist.append(y)
        tlist.append(t)
        while y >= 0.:
            t += dt
            y = v0*math.sin(theta)*t-0.5*g*t**2
            ylist.append(y)
            tlist.append(t)
        return [tlist,ylist]
