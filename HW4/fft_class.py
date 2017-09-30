################################################################################
# Author: Noah Green
# Assignment: PHY 514 Homework 4
# Description: Fast Fourier Transform Program
################################################################################

import numpy as np
import numpy.fft as ft
import time

class fourier:
    def __init__(self,vec,n=False):
        """
        Constructor for fourier transform class

        Input:
            vec: 1 dimensional array-like object
            n: int, number of terms to evaluate in vec
        """
        if not n:
            n = len(vec)
        vec = np.array(vec)
        if len(vec) < n:
            m = n-len(vec)
            add = np.zeros(m)
            vec.append(add)
        
        self.v = vec
        self.n = n
        self.t = 0.

    def dft(self):
        """
        Performs discrete Fourier transform
        """
        ft = np.zeros(self.n,dtype=complex)
        start  = time.time()
        const = 2.*np.pi/self.n
        for k in range(self.n):
            re = 0.
            im = 0.
            for i in range(self.n):
                val = const*k*i
                re += self.v[i]*np.cos(val)
                im -= self.v[i]*np.sin(val)
            ft[k] = np.complex(re,im)
        end = time.time()
        self.t = end-start
        return ft

    def fft(self):
        """
        Performs fast Fourier transform using Cooley Tukey algorithm.
        """
        

