# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 19:05:09 2014

@author: joshua
"""
import numpy as np

#
#  Classes
#
class IntProcModel:
    """Generalized Integrator Process Model"""
    #
    # class members
    #
    Kp = 0.0 #process gain in per seconds
    
    #
    # class methods
    #
    def __init__(self, proc_gain, init_pv = 0.0,
                 min_pv = 0.0, max_pv = 100.0):
        """Initialize Process Model, Process Gain in per unit time"""
        self.Kp = proc_gain
        self._pv = [init_pv,0.0]
        self._max = max_pv
        self._min = min_pv
        
    #
    # Iterate process model
    #
    def Iterate(self, input, dt, load = 0.0):
        """Delta T in same units as Kp"""
        self._pv = np.roll(self._pv,1) 
        self._pv[0] = sorted([self._max, self._min, \
            self.Kp*(input - load)*dt + self.pv[1]])[1]
        return self._pv[0]
        