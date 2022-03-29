# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 15:02:54 2014

@author: joshua
"""
import numpy as np

#
#  Classes
#
class LeadLag:
    """Lead Lag with Gain"""
    #
    # class members
    #
    K = 0.0 # gain
    T1 = 0.0 # lead time constant
    T2 = 0.0 # lag time constant
    #
    # class methods
    #
    def __init__(self,gain, lead_tc, lag_tc, init_in = 0.0,
                 init_out = 0.0, min_out = 0.0, max_out = 100.0):
        """Initialize Class"""
        self.K = gain
        self.T1 = lead_tc
        self.T2 = lag_tc
        self._out = [init_out,0.0]
        self._in = [init_in,0.0]
        self._min = min_out
        self._max = max_out
        
    def Iterate(self, input, dt):
        """Iterate the Lead Lag by Delta T,
        Delta T (dt) in same time units as time constant"""
        self._in = np.roll(self._in, 1)
        self._in[0] = input
        self._out = np.roll(self._out,1)
        uk = self.K*((self._in[0]*(dt + self.T1)) - self.T1*self._in[1])
        uk = (uk + self.T2*self._out[1])/(dt + self.T2)
        self._out[0] = sorted([self._max,self._min,uk])
        return self._out[0]
        
    #
    # set parameters
    #
    def SetParams(self, gain, lead_tc, lag_tc):
        self.K = gain
        self.T1 = lead_tc
        self.T2 = lag_tc
        return [self.K,self.T1,self.T2]