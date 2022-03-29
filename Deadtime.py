# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 15:02:54 2014

@author: joshua
"""
import numpy as np

#
#  Classes
#
class Deadtime:
    """Deadtime Block"""
    #
    # class members
    #
    bins = 0 # sample bins for process data
    #
    # class methods
    #
    def __init__(self, samples, init_out = 0.0):
        """Initialize deadtime, samples are samples per deadtime"""
        self.bins = max(int(samples),1)
        self._memory = []
        for i in range(self.bins):
            self._memory.append(init_out) 
        
    def Iterate(self, input):
        """Iterate the Shift Register"""
        self._memory = np.roll(self._memory, 1)
        self._memory[0] = input
        return self._memory[-1]
    
    def GetSample(self,sample):
        """Get Sample Less than Max Delay"""
        return self._memory[max(0,min(self.bins,int(sample)))]