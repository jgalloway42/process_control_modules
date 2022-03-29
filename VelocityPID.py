# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 19:05:09 2014

@author: joshua
"""
import numpy as np

# Constants
MAX_OUT = 100.0 # percent
MIN_OUT = 0.0 # percent

# convert to percent
def toPct(span,var):
    if span > 0.0:
        return var/span*100.0
    else:
        return 0.0

# calc error
def error(action,spp,pvp):
    if action:
        return pvp - spp
    else:
        return spp - pvp

# shift error register and add new error
def shftErr(e,ek):
    e = np.roll(e,1)
    e[0] = ek
    return e
#
#  Classes
#
class VelocityPID:
    """Velocity PID algorithm"""
    #
    # class members
    #
    kc = 0.0 # controller gain
    Ti = 0.0 # sec/rpt
    Td = 0.0 # sec
    err = [0,0,0] # error shift register
    
    
    #
    # class methods
    #
    def __init__(self, span, init_pv, init_sp, gain, reset, 
                 rate, dir_act = False,init_out = 0.0, ):
        """Initialize PID, reset and rate in seconds"""
        
        #initialize pid and setup 
        self.kc = gain
        self.Ti = max(1.0,reset) # reset must be positive
        self.Td = rate
        self._direct = dir_act #pid action
        self._pv_span = span #pv span
        self._out = init_out # output in %
        self._pv = init_pv # pv in eu
        self._pvp = toPct(span,init_pv) #pv in % of span
        self._sp = init_sp # sp in eu
        self._spp = toPct(span,init_sp) # sp in % span
        self._errp = error(dir_act,self._spp,self._pvp) # error in % span
        self.err = shftErr(self.err,self._errp)

    # 
    # Iterate PID Algorithm
    #
    def Iterate(self, pv, sp, dt):
        """Iterate Classical Velocity PID one time sample"""
        if dt <= 0.0:
            return self._out
            
        # calculate pid parameters
        param = [0,0,0]
        param[0] = self.kc*(1.0 + self.Td/dt)
        param[1] = self.kc*(-1.0)*(1.0 + 2.0*self.Td/dt - dt/self.Ti)
        param[2] = self.kc*(self.Td/dt)
        
        # update object data
        self._pv = pv
        self._sp = sp
        self._pvp = toPct(self._pv_span,pv)
        self._spp = toPct(self._pv_span,sp)
        self._errp = error(self._direct,self._spp,self._pvp)
        self.err = shftErr(self.err,self._errp)
        
        #calculate change in output
        dOut = 0.0
        for i in range(self.err.__len__()):
            dOut = dOut + self.err[i]*param[i]
            
        self._out += dOut
        self._out = sorted([self._out,MAX_OUT,MIN_OUT])[1]
        return self._out
        
    #
    # set tuning parameters
    #
    def SetTuning(self,kc,Ti,Td):
        self.kc = kc
        self.Ti = max(1.0,Ti)
        self.Td = Td
        return [self.kc,self.TI,self.Td]