# -*- coding: utf-8 -*-
"""
Created on Fri Apr 07 11:17:13 2017

generic centrifugal pump model asssumes NPSHR is met

@author: joshua
"""

import numpy as np

#
# Functions
#
def limit(x):
    """ Limit to 0 to 100 % """
    return min(100.,max(0.0,x))


def bhpPCT(FlowPCT):
    """ % Max Brake Horse Power Consumed at % Max Flow"""
    A2, A1, A0 = 93.77, 0.02490, 15.55
    FlowPCT = limit(FlowPCT)
    return limit(A2*(1.-np.exp(-A1*FlowPCT) + A0))

def headPCT(FlowPCT):
    """ % Max Head Developed at % Max Flow """
    A = np.array([-0.0057, 0.122, 100.])
    FlowPCT = limit(FlowPCT)
    return limit(np.polyval(A,FlowPCT))

def flowPCT(HeadPCT):
    """ % Max Flow at Developed Head % """
    A = np.array([-0.0057, 0.122, 100.])
    HeadPCT = limit(HeadPCT)
    h = np.poly1d(A)  # Create Head Polynomial
    # solve for roots at given head and take positive one
    return limit((h-HeadPCT).roots.max())

def effPCT(FlowPCT):
    """ % Efficiecy at % Max Flow """
    A = np.array([-0.0193, 2.51, 0.0])
    FlowPCT = limit(FlowPCT)
    return limit(np.polyval(A,FlowPCT))    

def calcLiqPower(Flow, Head, Eff = 70., SG = 1.0):
    """ Calculate Liquid Power Flow in GPM, Head in Feet H2O
    Pump Efficiency in Percent; return power in HP"""
    return Flow*Head*max(SG,0.001)/3960./max(Eff/100.,0.01)

#
#  Classes
#
class CentPump:
    """Generalized Liquid Valve Model"""
    #
    # class members
    #
    MaxHead = 100.0         # Maximum head developed (psi)
    MaxFlow = 100.0         # Maximum flow developed (gpm)
   
    #
    # class methods
    #
    def __init__(self, MaxFlow, MaxHead):
        """Initialize Centrifugal Pump
           Flow in GPM, Head in feet of water column """
        self.MaxFlow = MaxFlow
        self.MaxHead = MaxHead
    
    def SetMaxFlow(self, Flow):
        """Set Max Flow for Instance in GPM"""
        self.MaxFlow = Flow
        return Flow
    
    def SetMaxHead(self, Head):
        """Set Max Head in Feet of H2O"""
        self.MaxHead = Head
        return Head        
        
    def GetHead(self, Flow, SpeedPCT = 100., SG = 1.0):
        """ Get the Developed Head given flow and
        speed of pump if VFD driven"""
        # Adjust for VFD speed using affinity laws
        adjMF = max(1.,SpeedPCT/100.*self.MaxFlow)
        adjMH = max(1.,(SpeedPCT/100.)**2.*self.MaxHead)
        return headPCT(Flow/adjMF*100.)/100.*adjMH 
    
    def GetHeadPSI(self, Flow, SpeedPCT = 100., SG = 1.0):
        """ Get the Developed Head in PSI given flow and
        speed of pump if VFD driven"""
        return self.feetToPsi(self.GetHead(Flow,SpeedPCT,SG)) 
        
    def GetEFF(self, Flow, SpeedPCT = 100.):
        """ Get the Developed Head given flow and
        speed of pump if VFD driven"""
        # Adjust for VFD speed using affinity laws
        adjMF = SpeedPCT/100.*self.MaxFlow
        return effPCT(Flow/adjMF*100.)
    
    def GetFlow(self, Head, SpeedPCT = 100.):
        """ Get Flow given Head and speed of pump
           if driven by a VFD"""
        # Adjust for VFD speed using affinity laws
        adjMF = max(1.,SpeedPCT/100.*self.MaxFlow)
        adjMH = max(1.,(SpeedPCT/100.)**2.*self.MaxHead)
        return flowPCT(Head/adjMH*100.)/100.*adjMF
    
    def GetAmps(self, Flow, Head, SpeedPCT = 100., \
                LineVoltage = 480., PF = 90., MotorEff = 96.2):
        """ Get Amps given Flow, Head and electrical char """
        # Get power input at shaft in Watts
        adjMF = max(1.,SpeedPCT/100.*self.MaxFlow)
        p = calcLiqPower(Flow,Head,self.GetEFF(Flow/adjMF*100))*747.7
        p = p*100./max(1.,MotorEff) # get motor power input
        # calc amps
        return p/max(LineVoltage,120.)/max(PF,1.)
    
    def GetDeltaT(self, Flow, Head, SpeedPCT = 100., SpecificHeat = 1.):
        """ Get Delta T in Liquid given Flow in gpm, 
        Head in feet H2O, and Specific Heat Capacity
        of the liquid in BTU/(lbm-degF) return Temp Rise in degF """
        adjMF = max(1.,SpeedPCT/100.*self.MaxFlow)
        eff = max(1.,self.GetEFF(Flow/adjMF*100.))
        return Head*(1-eff/100.)/(788.*max(0.001,SpecificHeat)*eff/100.) 
    
    def psiToFeet(self, p, SG = 1.):
        """convert psi to feet of head"""
        return 2.30666*p/max(SG,0.001)

    def feetToPsi(self, p, SG = 1.):
        """convert feet of head to psi"""
        return p/2.30666*max(SG,0.001)


#
# if run as a script create test object
#
import matplotlib.pyplot as plt
if __name__ == "__main__":
    # test stub
    MaxFlow, MaxPSI= 5000., 150. # gpm and psi
    
        
    cp = CentPump(MaxFlow,1.)
    cp.SetMaxHead(cp.psiToFeet(MaxPSI))
    count = 100
    w = 2*np.pi*2./float(count)
    t = np.linspace(0.,count,count)
    f = np.zeros(count) + 0.5*MaxFlow
    speed = np.zeros_like(f)
        
    for i in range(np.shape(f)[0]):
        f[i] = 0.5*MaxFlow*(1.-np.sin(w*t[i]))
        speed[i] = 50.*(1.-np.cos(w*t[i]))
    
    h = np.zeros_like(f)
    eff = np.zeros_like(f)
    amps = np.zeros_like(f)
    dt = np.zeros_like(f)
    finv = np.zeros_like(f)
    
    hs = np.zeros_like(f)
    effs = np.zeros_like(f)
    ampss = np.zeros_like(f)
    dts = np.zeros_like(f)
    finvs = np.zeros_like(f)
    
    # Modulate Flow
    for i in range(np.shape(h)[0]):
        h[i] = cp.GetHead(f[i])
        eff[i] = cp.GetEFF(f[i])
        amps[i] = cp.GetAmps(f[i],h[i])
        dt[i] = cp.GetDeltaT(f[i],h[i])
        finv[i] = cp.GetFlow(h[i])
    
    # Modulate Speed
    sf = 0.5*MaxFlow
    for i in range(np.shape(hs)[0]):
        hs[i] = cp.GetHead(sf, speed[i])
        effs[i] = cp.GetEFF(sf, speed[i])
        ampss[i] = cp.GetAmps(sf,hs[i],speed[i])
        dts[i] = cp.GetDeltaT(sf,hs[i],speed[i])
        finvs[i] = cp.GetFlow(hs[i],speed[i])
    
    plt.figure()
    plt.subplot(6,2,1)
    plt.plot(t,f,'go',label='Flow')
    plt.xlabel('Time')
    plt.ylabel('GPM')
    plt.legend(loc='best')

    plt.subplot(6,2,2)
    plt.plot(t,h,'ro',label='Head')
    plt.xlabel('Time')
    plt.ylabel('FT H2O')
    plt.legend(loc='best')

    plt.subplot(6,2,3)
    plt.plot(t,eff,'ko',label="Efficiency")
    plt.xlabel('Time')
    plt.ylabel('%')
    plt.legend(loc='best')

    plt.subplot(6,2,4)
    plt.plot(t,amps,'mo',label='Amps')
    plt.xlabel('Time')
    plt.ylabel('Amps')
    plt.legend(loc='best')
    
    plt.subplot(6,2,5)
    plt.plot(t,dt,'yo',label='Delta T')
    plt.xlabel('Time')
    plt.ylabel('Deg F')
    plt.legend(loc='best')
        
    plt.subplot(6,2,6)
    plt.plot(t,finv,'bo',label='Flow Inverse')
    plt.xlabel('Time')
    plt.ylabel("GPM")
    plt.legend(loc='best')
    
    plt.subplot(6,2,7)
    plt.plot(t,speed,'gx',label='Speed')
    plt.xlabel('Time')
    plt.ylabel('%')
    plt.legend(loc='best')

    plt.subplot(6,2,8)
    plt.plot(t,hs,'rx',label='Head')
    plt.xlabel('Time')
    plt.ylabel('FT H2O')
    plt.legend(loc='best')

    plt.subplot(6,2,9)
    plt.plot(t,effs,'kx',label="Efficiency")
    plt.xlabel('Time')
    plt.ylabel('%')
    plt.legend(loc='best')

    plt.subplot(6,2,10)
    plt.plot(t,ampss,'mx',label='Amps')
    plt.xlabel('Time')
    plt.ylabel('Amps')
    plt.legend(loc='best')
    
    plt.subplot(6,2,11)
    plt.plot(t,dts,'yx',label='Delta T')
    plt.xlabel('Time')
    plt.ylabel('Deg F')
    plt.legend(loc='best')
        
    plt.subplot(6,2,12)
    plt.plot(t,finvs,'bx',label='Flow Inverse')
    plt.xlabel('Time')
    plt.ylabel("GPM")
    plt.legend(loc='best')

