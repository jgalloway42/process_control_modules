# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 19:05:09 2014

@author: joshua
"""
import numpy as np

#
#  Classes
#
class SelfRegProcModel:
    """First Order Self Regulator Process Model"""
    #
    # class members
    #
    Kp = 0.0 # process gain
    Tc = 0.0 # process first order time constant
    
    #
    # class methods
    #
    def __init__(self, proc_gain, time_const, init_pv = 0.0,
              min_pv = 0.0, max_pv = 100.0):
        """Initialize Class"""
        self.Kp = proc_gain
        self.Tc = time_const
        self._pv = [init_pv, 0.0]
        self._min = min_pv
        self._max = max_pv
        
    def Iterate(self, input, dt):
        """Iterate the process model by Delta T,
        Delta T (dt) in same time units as time constant"""
        self._pv =np.roll(self._pv, 1)
        self._pv[0] = sorted([self._min, self._max,
            (self._Kp*input + self.Tc/dt*self._pv[1])/(1.0+self.Tc/dt)])[1]
        return self._pv[0]