# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 19:05:09 2014

@author: joshua
Models from Modeling and Simulation in Chemical Engineering, Process control
section
"""
import numpy as np

#
#  Classes
#
class LiquidValve:
    """Generalized Liquid Valve Model"""
    #
    # class members
    #
    MaxCv = 0.0  # max valve Cv
    A0 = 0.05 # 1/rangeability of valve
    _pct = 0.0 # Percent valve opening
    types = ['Linear','EqPct'] # type of valve curves
    _type = types[0] # Linear valve 

    
    #
    # class methods
    #
    def __init__(self, max_cv, curve = 'Linear', rangeability = 50.0, init_pct = 0.0,):
        """Initialize Liquid Valve"""
        self.MaxCv = max_cv
        # rangeability is ratio max controllable/min controllable flow
        self.A0 = 1./max(2.,rangeability)
        self._pct = init_pct
        i = 0 # index of valve curve type
        try:
            i = self.types.index(curve)
        except:
            self._type = self.types[0]
        else:
            self._type = self.types[i]
            
        
    #
    # Calculate flow
    #
    def GetFlow(self, Pu, Pd, pct):
        """ Calculate flow from Up and Downstream 
        Pressures and % Valve Opening"""
        self._pct = sorted([0.,pct,100.])[1] # limit valve opening to 0 to 100%
        
        if (self._type == self.types[1]):
            #Equal Percentage
           A  = np.power(self.A0, 1.- self._pct/100.)
        else:
            # Linear
            A = self.A0 + (1. - self.A0)*self._pct/100.
        
        return A*self.MaxCv*np.sqrt(max(0.0,Pu-Pd))
            
    #
    # Calculate Cv at a % open
    #
    def GetCV(self, pct):
        """ Returns CV at a given Valve Opening %"""
        self._pct = sorted([0.,pct,100.])[1] # limit valve opening to 0 to 100%
        
        if (self._type == self.types[1]):
            #Equal Percentage
           A  = np.power(self.A0, 1.- self._pct/100.)
        else:
            # Linear
            A = self.A0 + (1. - self.A0)*self._pct/100.
        
        return A*self.MaxCv
    
    #
    # Set Max CV 
    #
    def SetMaxCV(self,MaxCV):
        """ Set the Max CV """
        self.MaxCV = max(0.,MaxCV)
        return self.MaxCV
    
    #
    # Calculate Delta P
    #
    def GetDP(self, flow, pct):
        """ Returns Pressure Drop at a given flow and Valve Opening %"""
        CV = self.GetCV(pct)
        return flow**2./max(0.001,CV)**2.
                    
import matplotlib.pyplot as plt
#
# if run as a script create test object
#
if __name__ == "__main__":
    r = 20. # 20:1 rangeability
    MaxCV = 200.  # GPM/sqrt(psi)
    pct = np.linspace(0.,100.,dtype='float')
    lv = LiquidValve(MaxCV,curve='EqPct',rangeability=r)
    op = np.zeros_like(pct)
    for i in range(pct.shape[0]):
        op[i] = lv.GetCV(pct[i])
    plt.plot(pct,op,'ro',label='EqPct Valve')
    plt.xlabel('% Open')
    plt.ylabel("Cv")
    plt.grid()
    plt.show()