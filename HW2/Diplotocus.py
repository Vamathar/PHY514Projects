################################################################################
# Class for making various plots
# Author: Noah Green
################################################################################
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator 

class diplotocus:
    def __init__(self,title='Title',xlabel='x',ylabel='y'):
        self.title=title
        self.xlabel=xlabel
        self.ylabel=ylabel


        self.__errors = {
                'xlabel':'Error: invalid x label.',
                'ylabel':'Error: invalid y label.',
                'title':'Error: invalid title.',
                }
        self.__verify() 

    def __verify(self):
        
        if type(self.title) == str:
            errors = list()
            if type(self.xlabel)!=str:
                errors.append('xlabel')
            if type(self.ylabel)!=str:
                errors.append('ylabel')
            
            if len(errors)!=0:
                error_out = ''
                for name in errors:
                    print(self.__errors[name])
                sys.exit(1)
        # todo: add multiple plot functionality
        # if type(self.title)==list or type(self.title)==tuple:
        else:
            sys.exit(self.__errors['title'])
    def Contour(self,x,y,z,res,colormap,filename):
        fig,axes = plt.subplots(1,1)
        levels = MaxNLocator(nbins = res).tick_values(z.max(),z.min())
        cmap = plt.get_cmap(colormap)
        axes.set_title(self.title)
        axes.set_xlabel(self.xlabel)
        axes.set_ylabel(self.ylabel)
        conplot = axes.contourf(x,y,z,levels=levels,cmap=cmap)
        fig.colorbar(conplot,ax=axes)
        fig.savefig(filename,format='png')

